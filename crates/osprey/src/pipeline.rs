//! Main analysis pipeline for Osprey
//!
//! This module orchestrates the complete analysis workflow.
//!
//! ## Multi-File RT Calibration Strategy
//!
//! When RT calibration is enabled (default), the pipeline uses a multi-file strategy:
//!
//! ### First File: Calibration Discovery
//! 1. Use ALL library peptides (no sampling)
//! 2. Assume library RT range ≈ mzML RT range
//! 3. Wide initial tolerance (25% of gradient range)
//! 4. Detect peaks and record (library_RT, measured_apex_RT) pairs
//! 5. Fit LOESS calibration curve
//! 6. Calculate residual SD for tight tolerance
//!
//! ### Subsequent Files: Calibrated Search
//! 1. Reuse calibration from first file (same experiment → similar LC conditions)
//! 2. Use tight RT tolerance (3× residual SD from first file)
//! 3. Run full regression with calibrated candidate selection
//!
//! **Rationale**: Files within the same experiment have similar LC conditions,
//! so calibrating once with the first file and reusing for subsequent files
//! is both efficient and accurate.
//!
//! ## Candidate Selection
//!
//! Candidate selection uses:
//! - **Isolation window**: From mzML file (defines which precursors are fragmented)
//! - **RT tolerance**: Calibrated or fallback (if calibration disabled/fails)
//!
//! Note: No additional precursor m/z tolerance is applied - the isolation window
//! from the mzML is sufficient.

use indicatif::{ProgressBar, ProgressStyle};
use osprey_chromatography::{
    PeakDetector, RTCalibration, RTCalibrator, RTCalibratorConfig,
    // Full calibration types
    CalibrationParams, CalibrationMetadata, MzCalibration, RTCalibrationParams, RTCalibrationMethod,
    IsolationScheme, MzQCData, calculate_mz_calibration, save_calibration, load_calibration,
    calibration_path_for_input,
};
use osprey_core::{
    BinConfig, DecoyMethod as CoreDecoyMethod, FeatureSet, LibraryEntry,
    MS1Spectrum, OspreyConfig, OspreyError, RegressionResult, ResolutionMode, Result, Spectrum,
};
use osprey_fdr::{FdrController, MokapotRunner, PsmFeatures};
use osprey_io::{load_library, load_all_spectra, BlibWriter, MS1Index};
use osprey_regression::{Binner, BinnedLibrary, BinnedSpectraCache, DesignMatrixBuilder, OptimizedSolver, RidgeSolver};
use osprey_scoring::{
    batch::{
        BatchScorer, MS1SpectrumLookup, PreprocessedLibrary, PreprocessedSpectra,
        run_libcosine_calibration_scoring_with_ms1,
    },
    DecoyGenerator, DecoyMethod, Enzyme, FeatureExtractor,
};

#[cfg(feature = "streaming")]
use osprey_scoring::pipeline::PipelineConfig;

/// Wrapper to implement MS1SpectrumLookup for MS1Index
struct MS1IndexWrapper<'a>(&'a MS1Index);

impl<'a> MS1SpectrumLookup for MS1IndexWrapper<'a> {
    fn find_nearest(&self, retention_time: f64) -> Option<&MS1Spectrum> {
        self.0.find_nearest(retention_time)
    }
}
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;

use osprey_scoring::batch::{CalibrationMatch, pair_calibration_matches};

/// Extract isolation window scheme from the first cycle of MS2 spectra
///
/// Returns an IsolationScheme describing the DIA windows used in the acquisition.
fn extract_isolation_scheme(spectra: &[Spectrum]) -> Option<IsolationScheme> {
    if spectra.is_empty() {
        return None;
    }

    // Find the first complete cycle by looking for repeated window centers
    let mut first_cycle_windows: Vec<(f64, f64)> = Vec::new();
    let mut seen_centers: std::collections::HashSet<i32> = std::collections::HashSet::new();

    for spectrum in spectra {
        let center = spectrum.isolation_window.center;
        let width = spectrum.isolation_window.width();

        // Round center to integer to detect repeats (handles slight variations)
        let center_key = (center * 10.0).round() as i32;

        if seen_centers.contains(&center_key) {
            // We've completed a cycle
            break;
        }

        seen_centers.insert(center_key);
        first_cycle_windows.push((center, width));
    }

    if first_cycle_windows.is_empty() {
        return None;
    }

    // Sort windows by center m/z
    first_cycle_windows.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mz_min = first_cycle_windows.first().map(|(c, _)| *c).unwrap_or(0.0);
    let mz_max = first_cycle_windows.last().map(|(c, _)| *c).unwrap_or(0.0);

    // Calculate typical width and check uniformity
    let widths: Vec<f64> = first_cycle_windows.iter().map(|(_, w)| *w).collect();
    let typical_width = widths.iter().sum::<f64>() / widths.len() as f64;
    let uniform_width = widths.iter().all(|w| (w - typical_width).abs() < 0.5);

    Some(IsolationScheme {
        num_windows: first_cycle_windows.len(),
        mz_min,
        mz_max,
        typical_width,
        uniform_width,
        windows: first_cycle_windows,
    })
}

/// Run the complete Osprey analysis pipeline
pub fn run_analysis(config: OspreyConfig) -> Result<Vec<RegressionResult>> {
    log::info!("Starting Osprey analysis");

    // Validate config
    if config.input_files.is_empty() {
        return Err(OspreyError::ConfigError("No input files specified".into()));
    }

    // Load library
    log::info!(
        "Loading spectral library from {:?}",
        config.library_source.path()
    );
    let mut library = load_library(&config.library_source)?;
    log::info!("Loaded {} library entries", library.len());

    if library.is_empty() {
        return Err(OspreyError::LibraryLoadError("Library is empty".into()));
    }

    // Generate decoys if not already in library
    if !config.decoys_in_library {
        log::info!("Generating decoys using {:?} method", config.decoy_method);

        let decoy_method = match config.decoy_method {
            CoreDecoyMethod::Reverse => DecoyMethod::Reverse,
            CoreDecoyMethod::Shuffle => DecoyMethod::Shuffle,
            CoreDecoyMethod::FromLibrary => {
                // Shouldn't reach here since decoys_in_library is false
                DecoyMethod::Reverse
            }
        };

        // Generate decoys with collision detection (pyXcorrDIA approach)
        let generator = DecoyGenerator::with_enzyme(decoy_method, Enzyme::Trypsin);

        // Only generate decoys for targets
        let targets: Vec<LibraryEntry> = library.iter().filter(|e| !e.is_decoy).cloned().collect();
        let n_original_targets = targets.len();

        // Generate decoys with collision detection
        // This filters out targets that can't have unique decoys
        let (valid_targets, decoys, _stats) = generator.generate_all_with_collision_detection(&targets);

        log::info!(
            "Generated {} decoys from {} targets ({} excluded due to collisions)",
            decoys.len(),
            n_original_targets,
            n_original_targets - valid_targets.len()
        );

        // Replace library with valid targets + decoys
        library = valid_targets;
        library.extend(decoys);
        log::info!("Total library size: {} (targets + decoys)", library.len());
    } else {
        // Count existing decoys
        let n_decoys = library.iter().filter(|e| e.is_decoy).count();
        let n_targets = library.len() - n_decoys;
        log::info!(
            "Library contains {} targets and {} decoys",
            n_targets,
            n_decoys
        );
    }

    // Set up binning based on resolution mode
    let bin_config = match config.resolution_mode {
        ResolutionMode::UnitResolution => BinConfig::unit_resolution(),
        ResolutionMode::HRAM => BinConfig::hram(),
        ResolutionMode::Auto => {
            // Default to unit resolution for now
            // TODO: Auto-detect from data
            BinConfig::unit_resolution()
        }
    };

    let binner = Binner::new(bin_config);

    // Get regularization lambda
    let lambda = match &config.regularization_lambda {
        osprey_core::RegularizationSetting::Fixed(l) => *l,
        _ => 1.0, // Default lambda for now
    };

    // Create optimized f32 solver (replaces old f64 RidgeSolver)
    let optimized_solver = OptimizedSolver::new(lambda as f32);

    // Build lookup structure for candidates
    let library_by_mz = build_mz_index(&library);

    // Pre-bin library once for efficient regression (major optimization)
    log::info!("Pre-binning library for optimized regression...");
    let binned_library = BinnedLibrary::new(&library, &binner);

    // Determine if we should use windowed scoring for large libraries
    // Windowed scoring is more memory-efficient for libraries > 100K entries
    const WINDOWED_SCORING_THRESHOLD: usize = 100_000;
    let use_windowed_scoring = library.len() > WINDOWED_SCORING_THRESHOLD;

    // Preprocess library for batch scoring (only for small libraries)
    let batch_scorer = BatchScorer::new();
    let preprocessed_library: Option<PreprocessedLibrary> = if use_windowed_scoring {
        log::info!(
            "Library has {} entries (> {}), using windowed batch scoring",
            library.len(),
            WINDOWED_SCORING_THRESHOLD
        );
        None
    } else {
        log::info!("Preprocessing library for batch scoring...");
        let lib = batch_scorer.preprocess_library(&library);
        log::info!(
            "Preprocessed {} library entries ({}×{} matrix)",
            lib.len(),
            lib.matrix.nrows(),
            lib.matrix.ncols()
        );
        Some(lib)
    };

    // Process each input file with multi-file calibration strategy:
    // - First file: calibrate with ALL peptides using wide tolerance
    // - Subsequent files: reuse calibration from first file
    let mut all_results = Vec::new();
    let mut all_spectra = Vec::new();
    // Cached preprocessed spectra for potential reuse in future scoring passes
    let mut _all_preprocessed_spectra: Option<PreprocessedSpectra> = None;
    let mut shared_calibration: Option<(RTCalibration, f64, CalibrationParams)> = None; // (rt_calibration, tolerance, full_params)

    for (file_idx, input_file) in config.input_files.iter().enumerate() {
        log::info!(
            "Processing file {}/{}: {}",
            file_idx + 1,
            config.input_files.len(),
            input_file.display()
        );

        let (file_results, file_spectra, file_preprocessed_spectra, calibration_result) =
            process_file_with_calibration(
                input_file,
                &library,
                &library_by_mz,
                &binned_library,
                &binner,
                &optimized_solver,
                &config,
                shared_calibration.as_ref(),
                preprocessed_library.as_ref(),
                &batch_scorer,
                use_windowed_scoring,
            )?;

        // Store calibration from first file for subsequent files
        if file_idx == 0 {
            if let Some((cal, tol, cal_params)) = calibration_result {
                log::info!(
                    "Storing calibration from first file for subsequent files (tolerance: {:.2} min)",
                    tol
                );
                shared_calibration = Some((cal, tol, cal_params));
            }
            // Store preprocessed spectra from first file for potential reuse
            _all_preprocessed_spectra = file_preprocessed_spectra;
        }

        all_results.extend(file_results);
        all_spectra.extend(file_spectra);
    }

    log::info!(
        "Analysis complete. {} total regression results across {} files",
        all_results.len(),
        config.input_files.len()
    );

    // ========== Score each file and run two-level FDR ==========

    // Score each file's results
    let mut per_file_scored: Vec<(String, Vec<ScoredEntry>)> = Vec::new();

    // Group results and spectra by file for per-file scoring
    // For now, treat all results as coming from a combined analysis
    // TODO: Track which results came from which file for proper per-file FDR
    let combined_file_name = config.input_files.first()
        .and_then(|p| p.file_stem())
        .and_then(|s| s.to_str())
        .unwrap_or("osprey")
        .to_string();

    log::info!("Scoring peptides and extracting features...");
    // Extract CalibrationParams for m/z correction during scoring
    let calibration_params = shared_calibration.as_ref().map(|(_, _, params)| params);
    let scored_entries = score_run(&library, &all_results, &all_spectra, &combined_file_name, calibration_params)?;
    per_file_scored.push((combined_file_name, scored_entries));

    // Run two-level FDR control (run-level + experiment-level)
    let output_dir = config.output_blib.parent().unwrap_or(std::path::Path::new("."));
    let experiment_entries = run_two_level_fdr(
        &mut per_file_scored,
        &library,
        config.run_fdr,
        config.run_fdr, // Use same FDR for experiment level for now
        output_dir,
    )?;

    // Filter to passing entries for blib output
    let passing_entries: Vec<ScoredEntry> = experiment_entries
        .iter()
        .filter(|e| !e.is_decoy && e.experiment_qvalue <= config.run_fdr)
        .cloned()
        .collect();

    log::info!(
        "Final results: {} precursors passing {}% experiment-level FDR",
        passing_entries.len(),
        config.run_fdr * 100.0
    );

    // Write blib output with Mokapot q-values
    if !passing_entries.is_empty() {
        log::info!("Writing blib to {}", config.output_blib.display());
        write_blib_output_with_scores(&config, &library, &passing_entries, &config.input_files)?;
    } else {
        log::warn!("No peptides passed FDR threshold, skipping blib output");
    }

    // Write output report if specified
    if let Some(report_path) = &config.output_report {
        log::info!("Writing report to {}", report_path.display());
        write_scored_report(report_path, &library, &passing_entries)?;
    }

    Ok(all_results)
}

/// Process a single mzML file with RT calibration
///
/// ## Multi-File Strategy
///
/// - **First file** (`existing_calibration` is None): Calibrate using ALL peptides
///   with wide tolerance, return calibration for reuse
/// - **Subsequent files** (`existing_calibration` is Some): Skip calibration discovery,
///   use provided calibration directly
///
/// Returns: (results, spectra, preprocessed_spectra, Option<(calibration, tolerance)>)
/// - spectra are returned for spectral scoring in FDR computation
/// - preprocessed_spectra is returned for reuse in subsequent scoring
/// - calibration is returned only for first file
fn process_file_with_calibration(
    path: &std::path::Path,
    library: &[LibraryEntry],
    library_by_mz: &HashMap<i32, Vec<usize>>,
    binned_library: &BinnedLibrary,
    binner: &Binner,
    optimized_solver: &OptimizedSolver,
    config: &OspreyConfig,
    existing_calibration: Option<&(RTCalibration, f64, CalibrationParams)>,
    preprocessed_library: Option<&PreprocessedLibrary>,
    batch_scorer: &BatchScorer,
    use_windowed_scoring: bool,
) -> Result<(
    Vec<RegressionResult>,
    Vec<Spectrum>,
    Option<PreprocessedSpectra>,
    Option<(RTCalibration, f64, CalibrationParams)>,
)> {
    // Load both MS1 and MS2 spectra in a single pass (more efficient than reading twice)
    let (spectra, ms1_index) = load_all_spectra(path)?;

    if spectra.is_empty() {
        return Ok((Vec::new(), Vec::new(), None, None));
    }

    // Check if RT calibration is enabled
    if !config.rt_calibration.enabled {
        log::info!(
            "RT calibration disabled, using fallback tolerance: {:.2} min",
            config.rt_calibration.fallback_rt_tolerance
        );
        let results = process_spectra_optimized(
            &spectra,
            library,
            library_by_mz,
            binned_library,
            binner,
            optimized_solver,
            config.rt_calibration.fallback_rt_tolerance,
            config.max_candidates_per_spectrum,
            None,
        )?;
        return Ok((results, spectra, None, None));
    }

    // Check if we have an existing calibration from a previous file
    if let Some((cal, tol, _cal_params)) = existing_calibration {
        log::info!(
            "Using calibration from first file (tolerance: {:.2} min)",
            tol
        );

        let results = process_spectra_optimized(
            &spectra,
            library,
            library_by_mz,
            binned_library,
            binner,
            optimized_solver,
            *tol,
            config.max_candidates_per_spectrum,
            Some(cal),
        )?;
        return Ok((results, spectra, None, None)); // Don't return calibration/preprocessing for subsequent files
    }

    // First file: Run calibration discovery with ALL peptides
    log::info!("First file: Calibration Discovery (using all peptides)");

    // Check if a valid calibration file already exists for this input file
    if let Some(output_dir) = config.output_blib.parent() {
        let cal_path = calibration_path_for_input(path, output_dir);
        if cal_path.exists() {
            log::info!("Found existing calibration file: {}", cal_path.display());

            match load_calibration(&cal_path) {
                Ok(cal_params) => {
                    // Check if calibration has model data for RT reconstruction
                    if cal_params.rt_calibration.has_model_data() && cal_params.is_calibrated() {
                        log::info!("Reusing existing RT calibration from file");

                        // Reconstruct RTCalibration from model params
                        let model_params = cal_params.rt_calibration.model_params.as_ref().unwrap();
                        match RTCalibration::from_model_params(model_params, cal_params.rt_calibration.residual_sd) {
                            Ok(rt_cal) => {
                                let tolerance = cal_params.rt_calibration.residual_sd * config.rt_calibration.rt_tolerance_factor;
                                let tolerance = tolerance.max(0.5);

                                log::info!(
                                    "Using cached calibration: {} points, R²={:.4}, tolerance={:.2} min",
                                    cal_params.rt_calibration.n_points,
                                    cal_params.rt_calibration.r_squared,
                                    tolerance
                                );

                                // Log calibration summary
                                cal_params.log_summary();

                                // Run full search with cached calibration
                                let results = process_spectra_optimized(
                                    &spectra,
                                    library,
                                    library_by_mz,
                                    binned_library,
                                    binner,
                                    optimized_solver,
                                    tolerance,
                                    config.max_candidates_per_spectrum,
                                    Some(&rt_cal),
                                )?;

                                return Ok((results, spectra, None, Some((rt_cal, tolerance, cal_params))));
                            }
                            Err(e) => {
                                log::warn!("Failed to reconstruct calibration from cached file: {}. Re-running calibration.", e);
                            }
                        }
                    } else {
                        log::info!("Cached calibration missing model data. Re-running calibration.");
                    }
                }
                Err(e) => {
                    log::warn!("Failed to load cached calibration: {}. Re-running calibration.", e);
                }
            }
        }
    }

    // Run calibration using appropriate method based on mode and library size
    #[cfg(feature = "streaming")]
    let calibration_result = if config.streaming {
        // Use streaming pipeline for memory-efficient calibration
        log::info!("Using streaming mode for calibration discovery");
        run_calibration_discovery_streaming(path, library, config)
    } else if use_windowed_scoring {
        // Use windowed scoring for large libraries (memory-efficient)
        run_calibration_discovery_windowed(library, &spectra, &ms1_index, config)
    } else if let Some(preproc_lib) = preprocessed_library {
        // Preprocess spectra for batch scoring (reused later)
        log::info!("Preprocessing {} spectra for batch scoring...", spectra.len());
        let preprocessed_spectra = batch_scorer.preprocess_spectra(&spectra);

        run_calibration_discovery_with_cache(
            library,
            preproc_lib,
            &preprocessed_spectra,
            config,
        )
    } else {
        // Fallback: use windowed scoring
        run_calibration_discovery_windowed(library, &spectra, &ms1_index, config)
    };

    #[cfg(not(feature = "streaming"))]
    let calibration_result = if use_windowed_scoring {
        // Use windowed scoring for large libraries (memory-efficient)
        run_calibration_discovery_windowed(library, &spectra, &ms1_index, config)
    } else if let Some(preproc_lib) = preprocessed_library {
        // Preprocess spectra for batch scoring (reused later)
        log::info!("Preprocessing {} spectra for batch scoring...", spectra.len());
        let preprocessed_spectra = batch_scorer.preprocess_spectra(&spectra);

        run_calibration_discovery_with_cache(
            library,
            preproc_lib,
            &preprocessed_spectra,
            config,
        )
    } else {
        // Fallback: use windowed scoring
        run_calibration_discovery_windowed(library, &spectra, &ms1_index, config)
    };

    // Determine RT tolerance and extract calibration
    let (rt_tolerance, calibration_opt, calibration_params_opt) = match calibration_result {
        Ok((rt_cal, cal_params)) => {
            // Use calibration residual SD × factor as tolerance
            let tolerance = cal_params.rt_calibration.residual_sd * config.rt_calibration.rt_tolerance_factor;
            let tolerance = tolerance.max(0.5); // Minimum 0.5 min tolerance
            log::info!(
                "Using calibrated RT tolerance: {:.2} min ({}× residual SD)",
                tolerance,
                config.rt_calibration.rt_tolerance_factor
            );

            // Save calibration to JSON using input file name (for reuse)
            if let Some(output_dir) = config.output_blib.parent() {
                let cal_path = calibration_path_for_input(path, output_dir);
                if let Err(e) = save_calibration(&cal_params, &cal_path) {
                    log::warn!("Failed to save calibration: {}", e);
                } else {
                    log::info!("Saved calibration to: {}", cal_path.display());
                }
            }

            (tolerance, Some(rt_cal), Some(cal_params))
        }
        Err(e) => {
            log::warn!(
                "Calibration failed: {}. Using fallback tolerance: {:.2} min",
                e,
                config.rt_calibration.fallback_rt_tolerance
            );
            (config.rt_calibration.fallback_rt_tolerance, None, None)
        }
    };

    // Log calibration summary if available
    if let Some(ref params) = calibration_params_opt {
        params.log_summary();
    }

    // Calibrated full search
    log::info!("Calibrated Full Search");

    let results = process_spectra_optimized(
        &spectra,
        library,
        library_by_mz,
        binned_library,
        binner,
        optimized_solver,
        rt_tolerance,
        config.max_candidates_per_spectrum,
        calibration_opt.as_ref(),
    )?;

    // Return calibration for first file to be reused
    // Note: preprocessed_spectra is only available when not using windowed scoring
    let calibration_to_share = match (calibration_opt, calibration_params_opt) {
        (Some(rt_cal), Some(cal_params)) => Some((rt_cal, rt_tolerance, cal_params)),
        _ => None,
    };
    Ok((results, spectra, None, calibration_to_share))
}

/// Run RT calibration discovery phase using BLAS-accelerated LibCosine scoring
///
/// This uses target-decoy competition to identify high-confidence detections,
/// then uses only peptides passing 1% FDR for RT calibration.
///
/// Process:
/// 1. Preprocess all library entries and spectra (once)
/// 2. Score all pairs via BLAS matrix multiplication (10-20× faster)
/// 3. Find best matches with RT filtering
/// 4. Target-decoy competition → q-values
/// 5. Use only peptides at 1% FDR for calibration
/// 6. Fit LOESS on (library_RT, measured_RT) pairs
///
/// Note: This is the standalone version. For cached preprocessing, use
/// `run_calibration_discovery_with_cache` instead.
#[allow(dead_code)]
fn run_calibration_discovery(
    spectra: &[Spectrum],
    library: &[LibraryEntry],
    _library_by_mz: &HashMap<i32, Vec<usize>>,
    _matrix_builder: &DesignMatrixBuilder,
    _solver: &RidgeSolver,
    config: &OspreyConfig,
) -> Result<RTCalibration> {
    let rt_config = &config.rt_calibration;

    // Calculate library RT range
    let library_rts: Vec<f64> = library
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.retention_time)
        .collect();
    let min_rt = library_rts.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_rt = library_rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let rt_range = max_rt - min_rt;

    // Get all entries with fragments (both targets and decoys)
    let entries_with_fragments: Vec<LibraryEntry> = library
        .iter()
        .filter(|e| !e.fragments.is_empty())
        .cloned()
        .collect();

    let n_targets = entries_with_fragments.iter().filter(|e| !e.is_decoy).count();
    let n_decoys = entries_with_fragments.iter().filter(|e| e.is_decoy).count();

    log::info!(
        "RT calibration: scoring {} targets + {} decoys (RT range: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        min_rt,
        max_rt
    );

    // === BLAS-accelerated batch scoring ===
    log::info!("Preprocessing library entries for batch scoring...");
    let batch_scorer = BatchScorer::new();
    let preprocessed_library = batch_scorer.preprocess_library(&entries_with_fragments);

    log::info!("Preprocessing {} spectra for batch scoring...", spectra.len());
    let preprocessed_spectra = batch_scorer.preprocess_spectra(spectra);

    // Get measured RT range from spectra
    let meas_min_rt = preprocessed_spectra
        .retention_times
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let meas_max_rt = preprocessed_spectra
        .retention_times
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let meas_rt_range = meas_max_rt - meas_min_rt;

    log::info!(
        "Measured RT range: {:.2}-{:.2} min ({:.2} min total)",
        meas_min_rt,
        meas_max_rt,
        meas_rt_range
    );

    // Check if RT ranges are similar enough to skip linear mapping
    // If slope would be close to 1.0 (within 20%), use library RTs directly
    // This avoids introducing mapping errors when library and measured RTs are already comparable
    let candidate_slope = if rt_range > 0.0 { meas_rt_range / rt_range } else { 1.0 };
    let use_direct_rts = (0.8..=1.2).contains(&candidate_slope);

    let (rt_slope, rt_intercept) = if use_direct_rts {
        log::info!("RT ranges are similar (slope would be {:.2}), using library RTs directly", candidate_slope);
        (1.0, 0.0)
    } else if rt_range > 0.0 {
        let slope = candidate_slope;
        let intercept = meas_min_rt - slope * min_rt;
        log::info!(
            "RT mapping: measured_RT = {:.4} × library_RT + {:.4}",
            slope,
            intercept
        );
        (slope, intercept)
    } else {
        // Degenerate case: all library entries have same RT
        (1.0, 0.0)
    };

    // RT tolerance depends on whether linear mapping was needed
    // - Similar RT scales (no mapping): use 20% of RT range
    // - Different RT scales (linear mapping): use 50% of RT range
    let tolerance_fraction = if use_direct_rts { 0.2 } else { 0.5 };
    let initial_tolerance = meas_rt_range * tolerance_fraction;
    let mapped_tolerance = initial_tolerance; // Already in measured RT scale

    log::info!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min mzML range)",
        initial_tolerance,
        tolerance_fraction * 100.0,
        meas_rt_range
    );

    // Build library RT lookup (mapped only if ranges differ significantly)
    let library_rts_mapped: Vec<f64> = preprocessed_library
        .entry_ids
        .iter()
        .filter_map(|id| {
            entries_with_fragments
                .iter()
                .find(|e| e.id == *id)
                .map(|e| rt_slope * e.retention_time + rt_intercept)
        })
        .collect();

    log::info!(
        "Computing LibCosine scores ({} × {} = {} pairs via BLAS)...",
        preprocessed_library.len(),
        preprocessed_spectra.len(),
        preprocessed_library.len() * preprocessed_spectra.len()
    );

    // Find best matches with RT filtering using BLAS
    // Now using library RTs mapped to measured RT scale
    let matches = batch_scorer.find_best_matches_with_rt_filter(
        &preprocessed_library,
        &preprocessed_spectra,
        &library_rts_mapped,
        mapped_tolerance,
    );

    log::info!("BLAS scoring complete: {} matches found", matches.len());

    // Build lookup for is_decoy
    let id_to_entry: HashMap<u32, &LibraryEntry> =
        entries_with_fragments.iter().map(|e| (e.id, e)).collect();

    // Convert matches to scoring format: (library_rt, measured_rt, score, is_decoy)
    let all_scores: Vec<(f64, f64, f64, bool)> = matches
        .into_iter()
        .filter_map(|(entry_id, _spec_idx, score, measured_rt)| {
            let entry = id_to_entry.get(&entry_id)?;
            Some((entry.retention_time, measured_rt, score, entry.is_decoy))
        })
        .collect();

    // Separate targets and decoys for FDR computation
    let mut target_results: Vec<(f64, f64, f64)> = Vec::new(); // (lib_rt, meas_rt, score)
    let mut decoy_scores: Vec<f64> = Vec::new();

    for (lib_rt, meas_rt, score, is_decoy) in all_scores {
        if is_decoy {
            decoy_scores.push(score);
        } else {
            target_results.push((lib_rt, meas_rt, score));
        }
    }

    log::info!(
        "LibCosine results: {} targets, {} decoys with scores > 0",
        target_results.len(),
        decoy_scores.len()
    );

    // Compute q-values using target-decoy competition
    let calibration_fdr = 0.01; // 1% FDR for calibration
    let fdr_controller = FdrController::new(calibration_fdr);

    let target_scores: Vec<f64> = target_results.iter().map(|(_, _, s)| *s).collect();
    let target_qvalues = fdr_controller.compute_qvalues(&target_scores, &decoy_scores)?;

    // Filter to 1% FDR and extract calibration points
    let mut library_rts_detected: Vec<f64> = Vec::new();
    let mut measured_rts_detected: Vec<f64> = Vec::new();

    for (i, (lib_rt, meas_rt, _score)) in target_results.iter().enumerate() {
        if i < target_qvalues.len() && target_qvalues[i] <= calibration_fdr {
            library_rts_detected.push(*lib_rt);
            measured_rts_detected.push(*meas_rt);
        }
    }

    log::info!(
        "RT calibration: {} peptides at 1% FDR (from {} targets, {} decoys)",
        library_rts_detected.len(),
        target_results.len(),
        decoy_scores.len()
    );

    if library_rts_detected.len() < rt_config.min_calibration_points {
        return Err(OspreyError::ConfigError(format!(
            "Insufficient calibration points: {} < {} required",
            library_rts_detected.len(),
            rt_config.min_calibration_points
        )));
    }

    // Fit LOESS calibration
    let calibrator_config = RTCalibratorConfig {
        bandwidth: rt_config.loess_bandwidth,
        degree: 1,
        min_points: rt_config.min_calibration_points,
        robustness_iter: 2,
    };
    let calibrator = RTCalibrator::with_config(calibrator_config);

    calibrator.fit(&library_rts_detected, &measured_rts_detected)
}

/// Run full calibration discovery using windowed batch scoring
///
/// This is the memory-efficient version for large libraries. Instead of
/// preprocessing the entire library at once, it groups spectra by isolation
/// window and only preprocesses library entries within each window.
///
/// For a typical 3mz DIA window with a 500 m/z library range and 3M precursors,
/// each window only contains ~18K precursors, reducing memory by ~170×.
///
/// Returns both the RTCalibration (for prediction) and CalibrationParams (for saving/logging).
fn run_calibration_discovery_windowed(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    ms1_index: &MS1Index,
    config: &OspreyConfig,
) -> Result<(RTCalibration, CalibrationParams)> {
    let rt_config = &config.rt_calibration;

    // Calculate initial wide tolerance based on library RT range
    let library_rts: Vec<f64> = library
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.retention_time)
        .collect();
    let lib_min_rt = library_rts.iter().cloned().fold(f64::INFINITY, f64::min);
    let lib_max_rt = library_rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let lib_rt_range = lib_max_rt - lib_min_rt;

    // Calculate mzML RT range from spectra
    let mzml_min_rt = spectra.iter().map(|s| s.retention_time).fold(f64::INFINITY, f64::min);
    let mzml_max_rt = spectra.iter().map(|s| s.retention_time).fold(f64::NEG_INFINITY, f64::max);
    let mzml_rt_range = mzml_max_rt - mzml_min_rt;

    // Create linear RT mapping from library RT to expected mzML RT
    // If the library is in iRT (e.g., 0-100) and mzML is in minutes (e.g., 0-30),
    // we map library RT to expected measured RT using a simple linear transformation.
    // If ranges are similar (within 2x), assume library is already in the same units.
    let (rt_slope, rt_intercept, use_linear_mapping) = {
        let ranges_similar = lib_rt_range > 0.0
            && mzml_rt_range > 0.0
            && (lib_rt_range / mzml_rt_range).max(mzml_rt_range / lib_rt_range) < 2.0
            && (lib_min_rt - mzml_min_rt).abs() < lib_rt_range * 0.5;

        if ranges_similar {
            // Library and mzML are in similar units, use identity mapping
            log::info!(
                "RT mapping: library and mzML ranges are similar, using identity mapping"
            );
            (1.0, 0.0, false)
        } else {
            // Linear mapping: expected_rt = slope * lib_rt + intercept
            // Maps (lib_min_rt, lib_max_rt) -> (mzml_min_rt, mzml_max_rt)
            let slope = if lib_rt_range > 0.0 { mzml_rt_range / lib_rt_range } else { 1.0 };
            let intercept = mzml_min_rt - slope * lib_min_rt;
            log::info!(
                "RT mapping: linear transform from library ({:.2}-{:.2}) to mzML ({:.2}-{:.2} min)",
                lib_min_rt, lib_max_rt, mzml_min_rt, mzml_max_rt
            );
            log::info!(
                "RT mapping: expected_rt = {:.4} * library_rt + {:.4}",
                slope, intercept
            );
            (slope, intercept, true)
        }
    };

    // RT tolerance depends on whether linear mapping was needed
    // - Similar RT scales (no mapping): use 20% of RT range
    // - Different RT scales (linear mapping): use 50% of RT range
    let tolerance_fraction = if use_linear_mapping { 0.5 } else { 0.2 };
    let initial_tolerance = mzml_rt_range * tolerance_fraction;

    let n_targets = library.iter().filter(|e| !e.is_decoy).count();
    let n_decoys = library.iter().filter(|e| e.is_decoy).count();
    let n_with_fragments = library.iter().filter(|e| !e.fragments.is_empty()).count();

    log::info!(
        "Calibration: library has {} targets + {} decoys = {} total (library RT: {:.1}-{:.1}, mzML RT: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        n_targets + n_decoys,
        lib_min_rt,
        lib_max_rt,
        mzml_min_rt,
        mzml_max_rt
    );

    if n_with_fragments < library.len() {
        log::warn!(
            "{} entries have no fragments and will not be scored",
            library.len() - n_with_fragments
        );
    }

    // Use full library for calibration (XCorr via BLAS is fast enough)
    let calibration_library: Vec<LibraryEntry> = library.to_vec();
    log::info!(
        "Using {} entries for calibration",
        calibration_library.len()
    );

    log::info!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min mzML range)",
        initial_tolerance,
        tolerance_fraction * 100.0,
        mzml_rt_range
    );
    log::info!(
        "Fragment tolerance: {} {}",
        config.fragment_tolerance.tolerance,
        match config.fragment_tolerance.unit {
            osprey_core::ToleranceUnit::Ppm => "ppm",
            osprey_core::ToleranceUnit::Mz => "Th",
        }
    );

    // MS1 spectra were loaded in single pass with MS2 - use them for precursor mass calibration
    let has_ms1 = !ms1_index.is_empty();
    if has_ms1 {
        log::info!(
            "Using {} MS1 spectra for precursor calibration",
            ms1_index.len()
        );
    } else {
        log::warn!(
            "No MS1 spectra available. Falling back to isolation window center for precursor m/z."
        );
    }

    // Run LibCosine calibration scoring (ppm-based peak matching, NO binning)
    // Uses MS1 spectra to extract M+0 isotope peak for accurate precursor mass calibration
    let matches = if has_ms1 {
        run_libcosine_calibration_scoring_with_ms1(
            &calibration_library,
            spectra,
            Some(&MS1IndexWrapper(ms1_index)),
            config.fragment_tolerance,
            config.precursor_tolerance.tolerance, // ppm tolerance for isotope peak matching
            initial_tolerance,
        )
    } else {
        // Fallback: no MS1 spectra available
        run_libcosine_calibration_scoring_with_ms1::<MS1IndexWrapper>(
            &calibration_library,
            spectra,
            None,
            config.fragment_tolerance,
            config.precursor_tolerance.tolerance,
            initial_tolerance,
        )
    };

    log::info!("LibCosine scoring complete: {} matches found", matches.len());

    // Log MS2 error statistics
    let total_ms2_errors: usize = matches.iter().map(|m| m.ms2_mass_errors.len()).sum();
    let avg_n_matched: f64 = if !matches.is_empty() {
        matches.iter().map(|m| m.n_matched_fragments as f64).sum::<f64>() / matches.len() as f64
    } else {
        0.0
    };
    log::info!(
        "MS2 statistics: {} total fragment matches ({:.1} avg per entry)",
        total_ms2_errors,
        avg_n_matched
    );

    // Log MS1 (precursor) coverage (how many matches have M+0 peaks)
    // NOTE: Actual MS1 calibration statistics are computed AFTER FDR filtering
    let ms1_count_all = matches.iter().filter(|m| m.ms1_error.is_some()).count();
    log::info!(
        "MS1 extraction: {} of {} matches have M+0 peak (before FDR)",
        ms1_count_all,
        matches.len()
    );

    // Count targets and decoys
    let n_targets = matches.iter().filter(|m| !m.is_decoy).count();
    let n_decoys = matches.iter().filter(|m| m.is_decoy).count();

    log::info!(
        "LibCosine results: {} targets, {} decoys with scores > 0",
        n_targets,
        n_decoys
    );

    // Write debug CSV with all scored peptides for diagnosis (paired target-decoy format)
    // Use output_blib parent directory, or fall back to input file directory, or current directory
    let debug_dir = config.output_blib.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .map(|p| p.to_path_buf())
        .or_else(|| config.input_files.first().and_then(|f| f.parent()).map(|p| p.to_path_buf()))
        .unwrap_or_else(|| std::path::PathBuf::from("."));

    let debug_path = debug_dir.join("calibration_debug.csv");

    // Create linear RT mapping function for initial expected_rt estimation
    // This maps library RT to expected measured RT using the linear transform computed above
    let linear_rt_mapping = |lib_rt: f64| -> f64 {
        rt_slope * lib_rt + rt_intercept
    };

    // Pass the linear RT mapping for initial delta_rt calculation
    // This gives a crude estimate before LOESS calibration is fitted
    let expected_rt_fn: Option<&dyn Fn(f64) -> f64> = if use_linear_mapping {
        Some(&linear_rt_mapping)
    } else {
        // If ranges are similar, still use identity mapping so delta_rt is meaningful
        Some(&linear_rt_mapping)
    };

    if let Err(e) = write_calibration_debug_csv(&matches, &calibration_library, &debug_path, expected_rt_fn) {
        log::warn!("Failed to write calibration debug CSV: {}", e);
    } else {
        log::info!("Wrote paired calibration debug CSV to: {}", debug_path.display());
    }

    // Run target-decoy competition using E-value (Comet-style, pyXcorrDIA approach)
    // Each target competes with its paired decoy. Lower E-value wins. Ties go to decoy.
    // Then sort winners by -E-value (so lower E-value = higher score) and compute FDR
    let calibration_fdr = 0.01; // 1% FDR for calibration
    let fdr_controller = FdrController::new(calibration_fdr);

    // Convert matches to competition format: (item, -evalue, is_decoy, entry_id)
    // We negate E-value because FdrController uses "higher is better"
    // but E-value is "lower is better"
    let competition_input = matches.iter().map(|m| {
        (m.clone(), -m.evalue, m.is_decoy, m.entry_id)
    });

    let competition_result = fdr_controller.compete_and_filter(competition_input);

    log::info!(
        "Competition: {} target wins, {} decoy wins",
        competition_result.n_target_wins,
        competition_result.n_decoy_wins
    );

    // Extract calibration points + mass errors from passing targets
    let mut library_rts_detected: Vec<f64> = Vec::new();
    let mut measured_rts_detected: Vec<f64> = Vec::new();
    let mut mz_qc_data = MzQCData::new(config.fragment_tolerance.unit);

    for m in &competition_result.passing_targets {
        library_rts_detected.push(m.library_rt);
        measured_rts_detected.push(m.measured_rt);

        // Collect MS1 error for mass calibration (ppm or Th depending on resolution)
        if let Some(ms1_error) = m.ms1_error {
            mz_qc_data.add_ms1_error(ms1_error);
        }

        // Collect MS2 mass errors from matched fragments (for calibration)
        for &ms2_error in &m.ms2_mass_errors {
            mz_qc_data.add_ms2_error(ms2_error);
        }
    }

    let num_confident_peptides = library_rts_detected.len();

    log::info!(
        "Calibration: {} peptides at {:.0}% FDR (from {} target wins, {} decoy wins)",
        num_confident_peptides,
        calibration_fdr * 100.0,
        competition_result.n_target_wins,
        competition_result.n_decoy_wins
    );

    if num_confident_peptides < rt_config.min_calibration_points {
        return Err(OspreyError::ConfigError(format!(
            "Insufficient calibration points: {} < {} required",
            num_confident_peptides,
            rt_config.min_calibration_points
        )));
    }

    // Fit LOESS RT calibration
    let calibrator_config = RTCalibratorConfig {
        bandwidth: rt_config.loess_bandwidth,
        degree: 1,
        min_points: rt_config.min_calibration_points,
        robustness_iter: 2,
    };
    let calibrator = RTCalibrator::with_config(calibrator_config);
    let rt_calibration = calibrator.fit(&library_rts_detected, &measured_rts_detected)?;
    let rt_stats = rt_calibration.stats();

    // Calculate mass calibration from collected MS1 errors
    let (ms1_calibration, ms2_calibration) = calculate_mz_calibration(&mz_qc_data);

    log::info!(
        "MS1 calibration: mean={:.2} {}, SD={:.2} {} (from {} observations)",
        ms1_calibration.mean,
        ms1_calibration.unit,
        ms1_calibration.sd,
        ms1_calibration.unit,
        mz_qc_data.n_ms1()
    );

    // Extract isolation window scheme from spectra
    let isolation_scheme = extract_isolation_scheme(spectra);

    // Build full CalibrationParams
    let calibration_params = CalibrationParams {
        metadata: CalibrationMetadata {
            num_confident_peptides,
            num_sampled_precursors: matches.len(),
            calibration_successful: true,
            timestamp: chrono::Utc::now().to_rfc3339(),
            isolation_scheme,
        },
        ms1_calibration,
        ms2_calibration,
        rt_calibration: RTCalibrationParams {
            method: RTCalibrationMethod::LOESS,
            residual_sd: rt_stats.residual_std,
            n_points: rt_stats.n_points,
            r_squared: rt_stats.r_squared,
            model_params: Some(rt_calibration.export_model_params()),
        },
    };

    Ok((rt_calibration, calibration_params))
}

/// Run calibration discovery using pre-cached preprocessed data
///
/// This version reuses preprocessed library and spectra matrices,
/// avoiding redundant preprocessing for each file.
///
/// Note: This version does not collect MS1 PPM errors (mass calibration not available).
/// For full calibration including mass, use `run_calibration_discovery_windowed`.
fn run_calibration_discovery_with_cache(
    library: &[LibraryEntry],
    preprocessed_library: &PreprocessedLibrary,
    preprocessed_spectra: &PreprocessedSpectra,
    config: &OspreyConfig,
) -> Result<(RTCalibration, CalibrationParams)> {
    let rt_config = &config.rt_calibration;

    // Calculate initial wide tolerance based on library RT range
    let library_rts: Vec<f64> = library
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.retention_time)
        .collect();
    let min_rt = library_rts.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_rt = library_rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let rt_range = max_rt - min_rt;

    let n_targets = library.iter().filter(|e| !e.is_decoy && !e.fragments.is_empty()).count();
    let n_decoys = library.iter().filter(|e| e.is_decoy && !e.fragments.is_empty()).count();

    log::info!(
        "RT calibration (cached): scoring {} targets + {} decoys (RT range: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        min_rt,
        max_rt
    );

    // Get measured RT range from spectra
    let meas_min_rt = preprocessed_spectra
        .retention_times
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min);
    let meas_max_rt = preprocessed_spectra
        .retention_times
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let meas_rt_range = meas_max_rt - meas_min_rt;

    log::info!(
        "Measured RT range: {:.2}-{:.2} min ({:.2} min total)",
        meas_min_rt,
        meas_max_rt,
        meas_rt_range
    );

    // Check if RT ranges are similar enough to skip linear mapping
    // If slope would be close to 1.0 (within 20%), use library RTs directly
    // This avoids introducing mapping errors when library and measured RTs are already comparable
    let candidate_slope = if rt_range > 0.0 { meas_rt_range / rt_range } else { 1.0 };
    let use_direct_rts = (0.8..=1.2).contains(&candidate_slope);

    let (rt_slope, rt_intercept) = if use_direct_rts {
        log::info!("RT ranges are similar (slope would be {:.2}), using library RTs directly", candidate_slope);
        (1.0, 0.0)
    } else if rt_range > 0.0 {
        let slope = candidate_slope;
        let intercept = meas_min_rt - slope * min_rt;
        log::info!(
            "RT mapping: measured_RT = {:.4} × library_RT + {:.4}",
            slope,
            intercept
        );
        (slope, intercept)
    } else {
        (1.0, 0.0)
    };

    // RT tolerance depends on whether linear mapping was needed
    // - Similar RT scales (no mapping): use 20% of RT range
    // - Different RT scales (linear mapping): use 50% of RT range
    let tolerance_fraction = if use_direct_rts { 0.2 } else { 0.5 };
    let initial_tolerance = meas_rt_range * tolerance_fraction;
    let mapped_tolerance = initial_tolerance; // Already in measured RT scale

    log::info!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min mzML range)",
        initial_tolerance,
        tolerance_fraction * 100.0,
        meas_rt_range
    );

    // Build library RT lookup (mapped only if ranges differ significantly)
    let library_rts_mapped: Vec<f64> = preprocessed_library
        .entry_ids
        .iter()
        .filter_map(|id| {
            library
                .iter()
                .find(|e| e.id == *id)
                .map(|e| rt_slope * e.retention_time + rt_intercept)
        })
        .collect();

    log::info!(
        "Computing LibCosine scores ({} × {} = {} pairs via BLAS, using cached matrices)...",
        preprocessed_library.len(),
        preprocessed_spectra.len(),
        preprocessed_library.len() * preprocessed_spectra.len()
    );

    // Score all pairs using BLAS matrix multiplication
    // Now using library RTs mapped to measured RT scale
    let batch_scorer = BatchScorer::new();
    let matches = batch_scorer.find_best_matches_with_rt_filter(
        preprocessed_library,
        preprocessed_spectra,
        &library_rts_mapped,
        mapped_tolerance,
    );

    log::info!("BLAS scoring complete: {} matches found", matches.len());

    // Build lookup for is_decoy
    let id_to_entry: HashMap<u32, &LibraryEntry> =
        library.iter().map(|e| (e.id, e)).collect();

    // Convert matches to scoring format: (library_rt, measured_rt, score, is_decoy)
    let all_scores: Vec<(f64, f64, f64, bool)> = matches
        .into_iter()
        .filter_map(|(entry_id, _spec_idx, score, measured_rt)| {
            let entry = id_to_entry.get(&entry_id)?;
            Some((entry.retention_time, measured_rt, score, entry.is_decoy))
        })
        .collect();

    // Store total matches count before iteration moves the vector
    let num_sampled_precursors = all_scores.len();

    // Separate targets and decoys for FDR computation
    let mut target_results: Vec<(f64, f64, f64)> = Vec::new();
    let mut decoy_scores: Vec<f64> = Vec::new();

    for (lib_rt, meas_rt, score, is_decoy) in all_scores {
        if is_decoy {
            decoy_scores.push(score);
        } else {
            target_results.push((lib_rt, meas_rt, score));
        }
    }

    log::info!(
        "LibCosine results: {} targets, {} decoys with scores > 0",
        target_results.len(),
        decoy_scores.len()
    );

    // Compute q-values using target-decoy competition
    let calibration_fdr = 0.01;
    let fdr_controller = FdrController::new(calibration_fdr);

    let target_scores: Vec<f64> = target_results.iter().map(|(_, _, s)| *s).collect();
    let target_qvalues = fdr_controller.compute_qvalues(&target_scores, &decoy_scores)?;

    // Filter to 1% FDR and extract calibration points
    let mut library_rts_detected: Vec<f64> = Vec::new();
    let mut measured_rts_detected: Vec<f64> = Vec::new();

    for (i, (lib_rt, meas_rt, _score)) in target_results.iter().enumerate() {
        if i < target_qvalues.len() && target_qvalues[i] <= calibration_fdr {
            library_rts_detected.push(*lib_rt);
            measured_rts_detected.push(*meas_rt);
        }
    }

    let num_confident_peptides = library_rts_detected.len();

    log::info!(
        "Calibration: {} peptides at 1% FDR (from {} targets, {} decoys)",
        num_confident_peptides,
        target_results.len(),
        decoy_scores.len()
    );

    if num_confident_peptides < rt_config.min_calibration_points {
        return Err(OspreyError::ConfigError(format!(
            "Insufficient calibration points: {} < {} required",
            num_confident_peptides,
            rt_config.min_calibration_points
        )));
    }

    // Fit LOESS calibration
    let calibrator_config = RTCalibratorConfig {
        bandwidth: rt_config.loess_bandwidth,
        degree: 1,
        min_points: rt_config.min_calibration_points,
        robustness_iter: 2,
    };
    let calibrator = RTCalibrator::with_config(calibrator_config);
    let rt_calibration = calibrator.fit(&library_rts_detected, &measured_rts_detected)?;
    let rt_stats = rt_calibration.stats();

    // Build CalibrationParams (mass calibration not available in cached path)
    let calibration_params = CalibrationParams {
        metadata: CalibrationMetadata {
            num_confident_peptides,
            num_sampled_precursors,
            calibration_successful: true,
            timestamp: chrono::Utc::now().to_rfc3339(),
            isolation_scheme: None,
        },
        ms1_calibration: MzCalibration::uncalibrated(),
        ms2_calibration: MzCalibration::uncalibrated(),
        rt_calibration: RTCalibrationParams {
            method: RTCalibrationMethod::LOESS,
            residual_sd: rt_stats.residual_std,
            n_points: rt_stats.n_points,
            r_squared: rt_stats.r_squared,
            model_params: Some(rt_calibration.export_model_params()),
        },
    };

    log::info!(
        "RT calibration successful: {} points, R²={:.4}, residual_SD={:.3} min",
        rt_stats.n_points,
        rt_stats.r_squared,
        rt_stats.residual_std
    );
    log::warn!("Note: Mass calibration not available in cached preprocessing path");

    Ok((rt_calibration, calibration_params))
}

/// Optimized process_spectra using pre-binned library and f32 operations
///
/// Key optimizations:
/// - Uses pre-binned library (avoids redundant binning)
/// - Uses f32 instead of f64 (faster BLAS, less memory)
/// - Caches binned observed spectra
fn process_spectra_optimized(
    spectra: &[Spectrum],
    library: &[LibraryEntry],
    library_by_mz: &HashMap<i32, Vec<usize>>,
    binned_library: &BinnedLibrary,
    binner: &Binner,
    solver: &OptimizedSolver,
    rt_tolerance: f64,
    max_candidates: usize,
    calibration: Option<&RTCalibration>,
) -> Result<Vec<RegressionResult>> {
    // Pre-bin all observed spectra once
    log::info!("Pre-binning {} observed spectra...", spectra.len());
    let mut spectra_cache = BinnedSpectraCache::with_capacity(binner.n_bins(), spectra.len());
    spectra_cache.populate(spectra, binner);

    // Set up progress bar
    let pb = ProgressBar::new(spectra.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} spectra (optimized)",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Process spectra in parallel
    // Note: We can't mutate spectra_cache in parallel, so we use get() which is immutable
    let results: Vec<RegressionResult> = spectra
        .par_iter()
        .filter_map(|spectrum| {
            // Select candidates
            let candidates = select_candidates_with_calibration(
                spectrum,
                library,
                library_by_mz,
                rt_tolerance,
                max_candidates,
                None,
                calibration,
            );

            if candidates.is_empty() {
                pb.inc(1);
                return None;
            }

            // Get pre-binned observed spectrum
            let observed = spectra_cache.get(spectrum.scan_number)?;
            let observed_norm = observed.dot(observed);

            // Extract design matrix from pre-binned library (no rebinning!)
            let (design_matrix, library_ids) = binned_library.extract_design_matrix(&candidates);

            // Solve regression with f32
            let coefficients = match solver.solve_nonnegative(&design_matrix, observed, None) {
                Ok(c) => c,
                Err(_) => {
                    pb.inc(1);
                    return None;
                }
            };

            let coefficient_sum: f32 = coefficients.iter().sum();

            // Compute residual
            let predicted = design_matrix.dot(&coefficients);
            let diff = observed - &predicted;
            let residual = diff.dot(&diff);

            // Filter to non-zero coefficients and convert to RegressionResult
            let mut result = RegressionResult::new(spectrum.scan_number, spectrum.retention_time);
            for (id, coef) in library_ids.iter().zip(coefficients.iter()) {
                if *coef > 1e-6 {
                    result.library_ids.push(*id);
                    result.coefficients.push(*coef as f64);
                }
            }

            // Set contextual metrics
            result.residual = residual as f64;
            result.n_candidates = candidates.len() as u32;
            result.coefficient_sum = coefficient_sum as f64;
            result.observed_norm = observed_norm as f64;

            pb.inc(1);

            if result.coefficients.is_empty() {
                None
            } else {
                Some(result)
            }
        })
        .collect();

    pb.finish_with_message("Done");

    log::info!(
        "Generated {} regression results with non-zero coefficients (optimized)",
        results.len()
    );

    Ok(results)
}

/// Scored library entry after FDR control
#[derive(Debug, Clone)]
struct ScoredEntry {
    /// Library entry index
    lib_idx: usize,
    /// PSM identifier (for matching with Mokapot results)
    psm_id: String,
    /// Source file name
    #[allow(dead_code)]
    file_name: String,
    /// Coefficient time series
    rt_coef_pairs: Vec<(f64, f64)>,
    /// Extracted features
    features: FeatureSet,
    /// Primary score (LibCosine)
    score: f64,
    /// Run-level q-value (from Mokapot on single file)
    run_qvalue: f64,
    /// Experiment-level q-value (from Mokapot across files)
    experiment_qvalue: f64,
    /// Posterior error probability (from Mokapot)
    pep: f64,
    /// Is decoy
    is_decoy: bool,
}

/// Results from processing a single run (mzML file)
#[derive(Debug)]
#[allow(dead_code)]
struct RunLevelResults {
    /// File name
    file_name: String,
    /// All scored entries (targets + decoys)
    scored_entries: Vec<ScoredEntry>,
    /// Path to the PIN file written for this run
    pin_path: std::path::PathBuf,
}

/// Experiment-level aggregated result (best PSM per peptide across runs)
#[derive(Debug, Clone)]
#[allow(dead_code)]
struct ExperimentPsm {
    /// Modified peptide sequence (unique identifier)
    peptide: String,
    /// Best scoring entry across all runs
    best_entry: ScoredEntry,
    /// Number of runs where this peptide was detected
    n_runs_detected: u32,
}

/// Score peptides and extract features for a single run
///
/// Uses spectral similarity scoring (LibCosine) as the primary score.
/// Returns all scored entries (targets + decoys) for subsequent FDR computation.
///
/// If calibration is provided, applies m/z correction to spectra and uses
/// calibrated tolerance (3×SD) for fragment matching.
fn score_run(
    library: &[LibraryEntry],
    results: &[RegressionResult],
    spectra: &[Spectrum],
    file_name: &str,
    calibration: Option<&CalibrationParams>,
) -> Result<Vec<ScoredEntry>> {
    use osprey_scoring::{RegressionContext, SpectralScorer};
    use osprey_chromatography::calibration::{apply_spectrum_calibration, calibrated_tolerance};
    use osprey_core::ToleranceUnit;

    // Build ID-to-index map
    let id_to_index: HashMap<u32, usize> = library
        .iter()
        .enumerate()
        .map(|(idx, entry)| (entry.id, idx))
        .collect();

    // Aggregate results by library entry ID
    // Store both (RT, coef) pairs and indices to RegressionResults for context building
    let mut entry_data: HashMap<u32, Vec<(f64, f64)>> = HashMap::new(); // (RT, coef)
    let mut entry_result_indices: HashMap<u32, Vec<usize>> = HashMap::new(); // indices into results
    for (result_idx, result) in results.iter().enumerate() {
        for (lib_id, coef) in result.library_ids.iter().zip(result.coefficients.iter()) {
            entry_data
                .entry(*lib_id)
                .or_default()
                .push((result.retention_time, *coef));
            entry_result_indices
                .entry(*lib_id)
                .or_default()
                .push(result_idx);
        }
    }

    // Extract features and create scored entries (parallelized)
    let feature_extractor = FeatureExtractor::new();

    // Configure SpectralScorer with calibrated tolerance if available
    // Use 3×SD from MS2 calibration as the tolerance for fragment matching
    // This is unit-aware: uses Th for unit resolution, ppm for HRAM
    let spectral_scorer = if let Some(cal) = calibration {
        // Get calibrated tolerance in the appropriate unit
        let (tol_value, tol_unit) = calibrated_tolerance(
            &cal.ms2_calibration,
            20.0,  // default 20 ppm for HRAM
            ToleranceUnit::Ppm,
        );

        match tol_unit {
            ToleranceUnit::Mz => {
                log::debug!("Using calibrated fragment tolerance: {:.4} Th (3×SD)", tol_value);
                SpectralScorer::new().with_tolerance_da(tol_value)
            }
            ToleranceUnit::Ppm => {
                log::debug!("Using calibrated fragment tolerance: {:.2} ppm (3×SD)", tol_value);
                SpectralScorer::new().with_tolerance_ppm(tol_value)
            }
        }
    } else {
        SpectralScorer::new()
    };

    // Pre-calibrate all spectra once upfront (if calibration available)
    // This avoids cloning/calibrating spectra inside the per-peptide loop
    let calibrated_spectra: Vec<Spectrum> = if let Some(cal) = calibration {
        log::debug!("Pre-calibrating {} spectra with MS2 calibration", spectra.len());
        spectra
            .iter()
            .map(|s| apply_spectrum_calibration(s, &cal.ms2_calibration))
            .collect()
    } else {
        spectra.to_vec()
    };

    let peak_detector = PeakDetector::new().with_min_height(0.05);

    // Collect library IDs for parallel iteration
    let lib_ids: Vec<u32> = entry_data.keys().copied().collect();

    let scored_entries: Vec<ScoredEntry> = lib_ids
        .par_iter()
        .filter_map(|lib_id| {
            let data_points = entry_data.get(lib_id)?;
            if data_points.len() < 3 {
                return None;
            }

            let lib_idx = id_to_index.get(lib_id).copied()?;
            let entry = &library[lib_idx];

            // Sort by RT
            let mut sorted_data = data_points.clone();
            sorted_data.sort_by(|a, b| a.0.total_cmp(&b.0));

            // Use sorted data directly as RT-coef pairs for peak detection
            let rt_coef_pairs: Vec<(f64, f64)> = sorted_data.clone();

            // Detect peaks
            let peaks = peak_detector.detect(&rt_coef_pairs);
            if peaks.is_empty() {
                return None;
            }

            // Find apex RT (RT with maximum coefficient)
            let apex_rt = sorted_data
                .iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .map(|(rt, _)| *rt)
                .unwrap_or(entry.retention_time);

            // Get RT range from coefficient series for peak region filtering
            let rt_min = rt_coef_pairs.iter().map(|(rt, _)| *rt).fold(f64::INFINITY, f64::min);
            let rt_max = rt_coef_pairs.iter().map(|(rt, _)| *rt).fold(f64::NEG_INFINITY, f64::max);

            // Filter pre-calibrated spectra to those in the peak region:
            // 1. Isolation window contains the peptide's precursor m/z
            // 2. RT is within the coefficient series range
            // Note: We clone filtered spectra here, but the expensive calibration
            // was done once upfront, so this is just a cheap memory copy
            let peak_region_spectra: Vec<Spectrum> = calibrated_spectra
                .iter()
                .filter(|s| {
                    s.contains_precursor(entry.precursor_mz)
                        && s.retention_time >= rt_min
                        && s.retention_time <= rt_max
                })
                .cloned()
                .collect();

            // Build regression context from the RegressionResults that contain this peptide
            let regression_context = entry_result_indices
                .get(lib_id)
                .map(|indices| {
                    let result_refs: Vec<&RegressionResult> =
                        indices.iter().map(|&idx| &results[idx]).collect();
                    RegressionContext::from_results(&result_refs, *lib_id)
                });

            // Extract features with both mixed and deconvoluted scoring
            // - Mixed: spectral scores from raw observed spectrum at apex
            // - Deconvoluted: spectral scores from coefficient-weighted aggregated spectrum (apex ± 2)
            let mut features = feature_extractor.extract_with_deconvolution(
                entry,
                &rt_coef_pairs,
                &peak_region_spectra,
                Some(entry.retention_time),
            );

            // Apply contextual features from regression context
            if let Some(ctx) = &regression_context {
                FeatureExtractor::apply_regression_context(&mut features, ctx);
            }

            // Find apex spectrum from peak region (closest to apex RT) for spectral scoring
            let apex_spectrum = peak_region_spectra
                .iter()
                .min_by(|a, b| {
                    (a.retention_time - apex_rt)
                        .abs()
                        .partial_cmp(&(b.retention_time - apex_rt).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

            // Compute spectral score using LibCosine (primary scoring metric)
            let spectral_score = if let Some(spectrum) = apex_spectrum {
                spectral_scorer.lib_cosine(spectrum, entry)
            } else {
                osprey_scoring::SpectralScore::default()
            };

            // Primary score is LibCosine (spectral similarity)
            // Fall back to peak_apex if no spectral score
            let score = if spectral_score.lib_cosine > 0.0 {
                spectral_score.lib_cosine
            } else {
                features.peak_apex * 0.1 // Scale down chromatographic score
            };

            // Check if this is a decoy (high bit set in ID)
            let is_decoy = *lib_id & 0x80000000 != 0;

            // Generate PSM ID for matching with Mokapot results
            let psm_id = format!("{}_{}", file_name, lib_id);

            Some(ScoredEntry {
                lib_idx,
                psm_id,
                file_name: file_name.to_string(),
                rt_coef_pairs,
                features,
                score,
                run_qvalue: 1.0,        // Will be set by Mokapot
                experiment_qvalue: 1.0, // Will be set by experiment-level Mokapot
                pep: 1.0,               // Will be set by Mokapot
                is_decoy,
            })
        })
        .collect();

    log::info!(
        "Scored {} entries ({} targets, {} decoys) for {}",
        scored_entries.len(),
        scored_entries.iter().filter(|e| !e.is_decoy).count(),
        scored_entries.iter().filter(|e| e.is_decoy).count(),
        file_name
    );

    // Return ALL entries (targets + decoys) for subsequent Mokapot FDR
    Ok(scored_entries)
}

/// Compute FDR using Mokapot semi-supervised learning
///
/// This function provides an alternative FDR method using Mokapot,
/// which trains a semi-supervised model on the features to better
/// separate true and false detections.
#[allow(dead_code)]
fn compute_fdr_with_mokapot(
    library: &[LibraryEntry],
    scored_entries: &mut [ScoredEntry],
    file_name: &str,
    output_dir: &std::path::Path,
    train_fdr: f64,
    test_fdr: f64,
) -> Result<()> {
    let mokapot = MokapotRunner::new()
        .with_train_fdr(train_fdr)
        .with_test_fdr(test_fdr);

    // Check if mokapot is available
    if !mokapot.is_available() {
        log::warn!("Mokapot not available, falling back to simple target-decoy competition");
        log::warn!("Install mokapot with: pip install mokapot");
        return Ok(());
    }

    // Convert scored entries to PSM features for mokapot
    let psm_features: Vec<PsmFeatures> = scored_entries
        .iter()
        .enumerate()
        .map(|(idx, entry)| {
            let lib_entry = &library[entry.lib_idx];

            // Use library index as scan_number for mokapot fold splitting
            // In DIA/Osprey, each PSM represents a peptide detection across a chromatographic
            // peak (many spectra), not a single spectrum. The library index uniquely identifies
            // each precursor (peptide + charge), ensuring fold splitting groups by precursor
            // to avoid data leakage during cross-validation.
            PsmFeatures {
                psm_id: format!("{}_{}", file_name, idx),
                peptide: lib_entry.modified_sequence.clone(),
                proteins: lib_entry.protein_ids.clone(),
                scan_number: entry.lib_idx as u32,
                file_name: file_name.to_string(),
                charge: lib_entry.charge,
                is_decoy: entry.is_decoy,
                features: entry.features.clone(),
                initial_score: Some(entry.score),
            }
        })
        .collect();

    // Create output directory for mokapot
    let mokapot_dir = output_dir.join("mokapot");
    std::fs::create_dir_all(&mokapot_dir)?;

    // Write PIN file
    let pin_file = mokapot_dir.join(format!("{}.pin", file_name));
    mokapot.write_pin(&psm_features, &pin_file)?;

    log::info!("Wrote {} precursors to PIN file for mokapot", psm_features.len());

    // Run mokapot
    let results = mokapot.run(&pin_file, &mokapot_dir)?;

    // Build lookup from psm_id to mokapot results
    let result_map: std::collections::HashMap<String, (f64, f64)> = results
        .into_iter()
        .map(|r| (r.psm_id, (r.q_value, r.pep)))
        .collect();

    // Update run-level q-values and PEPs from mokapot results
    let mut updated = 0;
    for entry in scored_entries.iter_mut() {
        if let Some(&(q_value, pep)) = result_map.get(&entry.psm_id) {
            entry.run_qvalue = q_value;
            entry.pep = pep;
            updated += 1;
        }
    }

    log::info!("Updated {} entries with run-level Mokapot results", updated);

    Ok(())
}

/// Run two-level FDR control using Mokapot
///
/// This implements the two-level FDR strategy:
/// 1. Run-level FDR: Run Mokapot on each file independently → run_qvalue, pep
/// 2. Experiment-level FDR: Keep best PSM per peptide, run Mokapot again → experiment_qvalue
///
/// Returns all entries with both run-level and experiment-level q-values populated.
fn run_two_level_fdr(
    per_file_results: &mut Vec<(String, Vec<ScoredEntry>)>,
    library: &[LibraryEntry],
    run_fdr: f64,
    experiment_fdr: f64,
    output_dir: &std::path::Path,
) -> Result<Vec<ScoredEntry>> {
    let mokapot = MokapotRunner::new()
        .with_train_fdr(run_fdr)
        .with_test_fdr(run_fdr);

    let mokapot_available = mokapot.is_available();
    if !mokapot_available {
        log::warn!("Mokapot not available, using simple target-decoy competition");
        log::warn!("Install mokapot with: pip install mokapot");
    }

    // Create output directory for mokapot files
    let mokapot_dir = output_dir.join("mokapot");
    std::fs::create_dir_all(&mokapot_dir)?;

    // ========== STEP 1: Run-level FDR ==========
    log::info!("Step 1: Run-level FDR control at {}%", run_fdr * 100.0);

    for (file_name, entries) in per_file_results.iter_mut() {
        if mokapot_available {
            // Run Mokapot for this file
            compute_fdr_with_mokapot(
                library,
                entries,
                file_name,
                &mokapot_dir,
                run_fdr,
                run_fdr,
            )?;
        } else {
            // Fallback: simple target-decoy competition
            apply_simple_fdr(entries, run_fdr)?;
        }

        // Report per-file statistics
        let passing = entries.iter().filter(|e| !e.is_decoy && e.run_qvalue <= run_fdr).count();
        log::info!("  {}: {} precursors passing {}% run-level FDR", file_name, passing, run_fdr * 100.0);
    }

    // ========== STEP 2: Experiment-level FDR ==========
    // Skip experiment-level FDR if only one replicate - just use run-level results
    if per_file_results.len() == 1 {
        log::info!("Single replicate - skipping experiment-level FDR (using run-level results)");
        let (_, entries) = per_file_results.iter().next().unwrap();
        let mut experiment_entries: Vec<ScoredEntry> = entries.clone();
        // Copy run-level q-values to experiment-level
        for entry in experiment_entries.iter_mut() {
            entry.experiment_qvalue = entry.run_qvalue;
        }
        return Ok(experiment_entries);
    }

    log::info!("Step 2: Experiment-level FDR control at {}%", experiment_fdr * 100.0);

    // Collect all entries across files
    let all_entries: Vec<ScoredEntry> = per_file_results
        .iter()
        .flat_map(|(_, entries)| entries.iter().cloned())
        .collect();

    // Aggregate: keep best PSM per precursor (peptide + charge) by score
    let mut best_per_precursor: HashMap<String, ScoredEntry> = HashMap::new();
    for entry in all_entries.iter() {
        let lib_entry = &library[entry.lib_idx];
        // Use peptide + charge as the precursor key
        let precursor_key = format!("{}_{}", lib_entry.modified_sequence, lib_entry.charge);
        best_per_precursor
            .entry(precursor_key)
            .and_modify(|existing| {
                if entry.score > existing.score {
                    *existing = entry.clone();
                }
            })
            .or_insert_with(|| entry.clone());
    }

    let mut experiment_entries: Vec<ScoredEntry> = best_per_precursor.into_values().collect();
    log::info!(
        "Aggregated to {} unique precursors ({} targets, {} decoys)",
        experiment_entries.len(),
        experiment_entries.iter().filter(|e| !e.is_decoy).count(),
        experiment_entries.iter().filter(|e| e.is_decoy).count()
    );

    if mokapot_available {
        // Run experiment-level Mokapot
        run_experiment_level_mokapot(
            &mut experiment_entries,
            library,
            &mokapot_dir,
            experiment_fdr,
        )?;
    } else {
        // Fallback: simple target-decoy competition for experiment level
        apply_simple_fdr_experiment(&mut experiment_entries, experiment_fdr)?;
    }

    // Report experiment-level statistics
    let passing = experiment_entries.iter()
        .filter(|e| !e.is_decoy && e.experiment_qvalue <= experiment_fdr)
        .count();
    log::info!(
        "Experiment-level: {} precursors passing {}% FDR",
        passing,
        experiment_fdr * 100.0
    );

    Ok(experiment_entries)
}

/// Run Mokapot at experiment level (on aggregated best PSMs)
fn run_experiment_level_mokapot(
    entries: &mut [ScoredEntry],
    library: &[LibraryEntry],
    output_dir: &std::path::Path,
    fdr: f64,
) -> Result<()> {
    let mokapot = MokapotRunner::new()
        .with_train_fdr(fdr)
        .with_test_fdr(fdr);

    // Convert to PSM features with experiment-level IDs
    let psm_features: Vec<PsmFeatures> = entries
        .iter()
        .map(|entry| {
            let lib_entry = &library[entry.lib_idx];
            PsmFeatures {
                psm_id: format!("exp_{}", entry.psm_id),
                peptide: lib_entry.modified_sequence.clone(),
                proteins: lib_entry.protein_ids.clone(),
                scan_number: 0, // Not relevant for experiment level
                file_name: "experiment".to_string(),
                charge: lib_entry.charge,
                is_decoy: entry.is_decoy,
                features: entry.features.clone(),
                initial_score: Some(entry.score),
            }
        })
        .collect();

    // Write experiment-level PIN file
    let pin_file = output_dir.join("experiment.pin");
    mokapot.write_pin(&psm_features, &pin_file)?;
    log::info!("Wrote experiment-level PIN file: {}", pin_file.display());

    // Run Mokapot
    let results = mokapot.run(&pin_file, &output_dir.to_path_buf())?;

    // Build lookup from psm_id to results
    let result_map: HashMap<String, (f64, f64)> = results
        .into_iter()
        .map(|r| (r.psm_id, (r.q_value, r.pep)))
        .collect();

    // Update experiment-level q-values
    let mut updated = 0;
    for entry in entries.iter_mut() {
        let exp_psm_id = format!("exp_{}", entry.psm_id);
        if let Some(&(q_value, _pep)) = result_map.get(&exp_psm_id) {
            entry.experiment_qvalue = q_value;
            updated += 1;
        }
    }

    log::info!("Updated {} entries with experiment-level Mokapot results", updated);
    Ok(())
}

/// Fallback: Apply simple target-decoy FDR to run-level q-values
fn apply_simple_fdr(entries: &mut [ScoredEntry], fdr_threshold: f64) -> Result<()> {
    let fdr_controller = FdrController::new(fdr_threshold);

    let target_scores: Vec<f64> = entries
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.score)
        .collect();

    let decoy_scores: Vec<f64> = entries
        .iter()
        .filter(|e| e.is_decoy)
        .map(|e| e.score)
        .collect();

    let qvalues = fdr_controller.compute_qvalues(&target_scores, &decoy_scores)?;

    let mut target_idx = 0;
    for entry in entries.iter_mut() {
        if !entry.is_decoy {
            if target_idx < qvalues.len() {
                entry.run_qvalue = qvalues[target_idx];
                target_idx += 1;
            }
        }
    }

    Ok(())
}

/// Fallback: Apply simple target-decoy FDR to experiment-level q-values
fn apply_simple_fdr_experiment(entries: &mut [ScoredEntry], fdr_threshold: f64) -> Result<()> {
    let fdr_controller = FdrController::new(fdr_threshold);

    let target_scores: Vec<f64> = entries
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.score)
        .collect();

    let decoy_scores: Vec<f64> = entries
        .iter()
        .filter(|e| e.is_decoy)
        .map(|e| e.score)
        .collect();

    let qvalues = fdr_controller.compute_qvalues(&target_scores, &decoy_scores)?;

    let mut target_idx = 0;
    for entry in entries.iter_mut() {
        if !entry.is_decoy {
            if target_idx < qvalues.len() {
                entry.experiment_qvalue = qvalues[target_idx];
                target_idx += 1;
            }
        }
    }

    Ok(())
}

/// Write scored results to blib format for Skyline
fn write_blib_output_with_scores(
    config: &OspreyConfig,
    library: &[LibraryEntry],
    scored_entries: &[ScoredEntry],
    input_files: &[std::path::PathBuf],
) -> Result<()> {
    let mut writer = BlibWriter::create(&config.output_blib)?;

    // Add metadata
    writer.add_metadata("osprey_version", env!("CARGO_PKG_VERSION"))?;
    writer.add_metadata("rt_calibration_enabled", &config.rt_calibration.enabled.to_string())?;
    writer.add_metadata("run_fdr", &config.run_fdr.to_string())?;
    writer.add_metadata("fdr_method", "target_decoy_competition")?;

    // Add source files - build lookup by filename string for matching with ScoredEntry.file_name
    let mut file_name_to_id: HashMap<String, i64> = HashMap::new();
    for input_file in input_files {
        let file_name = input_file.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        let file_id = writer.add_source_file(file_name)?;
        file_name_to_id.insert(file_name.to_string(), file_id);
    }

    // Begin batch transaction for much faster writes
    writer.begin_batch()?;

    // Write spectra for passing entries
    for scored in scored_entries {
        let entry = &library[scored.lib_idx];

        // Find peak boundaries
        let max_coef = scored.rt_coef_pairs.iter().map(|(_, c)| *c).fold(0.0f64, f64::max);
        let threshold = max_coef * 0.1;

        let above_threshold: Vec<_> = scored.rt_coef_pairs.iter()
            .filter(|(_, c)| *c >= threshold)
            .collect();

        if above_threshold.is_empty() {
            continue;
        }

        let start_rt = above_threshold.first().map(|(rt, _)| *rt).unwrap_or(0.0);
        let end_rt = above_threshold.last().map(|(rt, _)| *rt).unwrap_or(0.0);
        let (apex_rt, _apex_coef) = scored.rt_coef_pairs.iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, c)| (*rt, *c))
            .unwrap_or((0.0, 0.0));

        // Look up file_id from the entry's file_name
        let file_id = *file_name_to_id.get(&scored.file_name).unwrap_or(&1);

        // Get fragment peaks from library entry
        let mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
        let intensities: Vec<f32> = entry.fragments.iter().map(|f| f.relative_intensity).collect();

        // Use 1 - q_value as the score (higher is better, blib convention)
        let blib_score = 1.0 - scored.experiment_qvalue;

        let ref_id = writer.add_spectrum(
            &entry.sequence,
            &entry.modified_sequence,
            entry.precursor_mz,
            entry.charge as i32,
            apex_rt,
            start_rt,
            end_rt,
            &mzs,
            &intensities,
            blib_score,
            file_id,
        )?;

        // Add peak boundaries - use file_name from the scored entry
        let boundaries = osprey_core::PeakBoundaries {
            start_rt,
            end_rt,
            apex_rt,
            apex_coefficient: scored.score,
            integrated_area: scored.features.peak_area,
            peak_quality: osprey_core::PeakQuality::default(),
        };

        writer.add_peak_boundaries(ref_id, &scored.file_name, &boundaries)?;

        // Add run scores - use run_qvalue for per-file q-value, not experiment q-value
        // run_qvalue: per-file FDR control result
        // experiment_qvalue: experiment-level FDR (used for filtering what goes in blib)
        // blib_score: 1 - experiment_qvalue (used as RefSpectra.score)
        writer.add_run_scores(
            ref_id,
            &scored.file_name,
            scored.run_qvalue,      // Run-level q-value
            scored.score,           // Discriminant score
            scored.pep,             // Posterior error probability
        )?;

        // Add experiment-level scores
        // For now, n_runs_detected = 1 since we only have the detection from one file
        // TODO: Track detections across files when multi-file per-file tracking is implemented
        writer.add_experiment_scores(
            ref_id,
            scored.experiment_qvalue,
            1,                      // n_runs_detected (currently always 1)
            input_files.len() as i32, // n_runs_searched
        )?;
    }

    // Commit the batch transaction
    writer.commit()?;
    writer.finalize()?;

    log::info!("Wrote {} spectra to blib", scored_entries.len());

    Ok(())
}

/// Write scored report to TSV
fn write_scored_report(
    path: &std::path::Path,
    library: &[LibraryEntry],
    scored_entries: &[ScoredEntry],
) -> Result<()> {
    use std::io::Write;

    let mut file = std::fs::File::create(path)?;

    // Write header
    writeln!(
        file,
        "modified_sequence\tprecursor_mz\tcharge\tapex_rt\tpeak_apex\tpeak_area\tn_scans\tq_value"
    )?;

    // Write entries
    for scored in scored_entries {
        let entry = &library[scored.lib_idx];
        let apex_rt = scored.rt_coef_pairs.iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, _)| *rt)
            .unwrap_or(0.0);

        writeln!(
            file,
            "{}\t{:.4}\t{}\t{:.2}\t{:.4}\t{:.4}\t{}\t{:.6}",
            entry.modified_sequence,
            entry.precursor_mz,
            entry.charge,
            apex_rt,
            scored.features.peak_apex,
            scored.features.peak_area,
            scored.features.n_contributing_scans,
            scored.experiment_qvalue
        )?;
    }

    Ok(())
}

/// Write scored entries to PIN file for Mokapot semi-supervised FDR control
///
/// Includes BOTH targets AND decoys, which is required for Mokapot's
/// semi-supervised learning approach.
#[allow(dead_code)]
fn write_pin_output(
    pin_path: &std::path::Path,
    library: &[LibraryEntry],
    scored_entries: &[ScoredEntry],
    file_name: &str,
) -> Result<()> {
    let mokapot_runner = MokapotRunner::new();

    // Convert ScoredEntries to PsmFeatures for all entries (targets AND decoys)
    let psm_features: Vec<PsmFeatures> = scored_entries
        .iter()
        .enumerate()
        .map(|(idx, scored)| {
            let entry = &library[scored.lib_idx];

            // Find apex scan from coefficient series (highest coefficient)
            let apex_scan = scored
                .rt_coef_pairs
                .iter()
                .enumerate()
                .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                .map(|(scan_idx, _)| scan_idx as u32)
                .unwrap_or(0);

            PsmFeatures {
                psm_id: format!("{}_{}", file_name, idx),
                peptide: entry.modified_sequence.clone(),
                proteins: entry.protein_ids.clone(),
                scan_number: apex_scan,
                file_name: file_name.to_string(),
                charge: entry.charge,
                is_decoy: scored.is_decoy,
                features: scored.features.clone(),
                initial_score: Some(scored.score),
            }
        })
        .collect();

    // Write PIN file
    mokapot_runner.write_pin(&psm_features, pin_path)?;

    log::info!(
        "Wrote PIN file with {} precursors ({} targets, {} decoys)",
        psm_features.len(),
        psm_features.iter().filter(|p| !p.is_decoy).count(),
        psm_features.iter().filter(|p| p.is_decoy).count()
    );

    Ok(())
}

/// Write results to blib format for Skyline (legacy, without FDR)
#[allow(dead_code)]
fn write_blib_output(
    config: &OspreyConfig,
    library: &[LibraryEntry],
    results: &[RegressionResult],
    input_files: &[std::path::PathBuf],
) -> Result<()> {
    let mut writer = BlibWriter::create(&config.output_blib)?;

    // Add metadata
    writer.add_metadata("osprey_version", env!("CARGO_PKG_VERSION"))?;
    writer.add_metadata("rt_calibration_enabled", &config.rt_calibration.enabled.to_string())?;
    writer.add_metadata("rt_tolerance_factor", &config.rt_calibration.rt_tolerance_factor.to_string())?;
    writer.add_metadata("fallback_rt_tolerance", &config.rt_calibration.fallback_rt_tolerance.to_string())?;
    writer.add_metadata("run_fdr", &config.run_fdr.to_string())?;

    // Add source files
    let mut file_ids = std::collections::HashMap::new();
    for input_file in input_files {
        let file_name = input_file.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        let file_id = writer.add_source_file(file_name)?;
        file_ids.insert(input_file.clone(), file_id);
    }

    // Begin batch transaction for much faster writes
    writer.begin_batch()?;

    // Build ID-to-index map for looking up library entries
    let id_to_index: HashMap<u32, usize> = library
        .iter()
        .enumerate()
        .map(|(idx, entry)| (entry.id, idx))
        .collect();

    // Aggregate results by library entry ID
    // For each library entry that has non-zero coefficients, create a spectrum entry
    let mut entry_coefficients: HashMap<u32, Vec<(f64, f64)>> = HashMap::new();

    for result in results {
        for (lib_id, coef) in result.library_ids.iter().zip(result.coefficients.iter()) {
            entry_coefficients
                .entry(*lib_id)
                .or_default()
                .push((result.retention_time, *coef));
        }
    }

    // Write spectra for entries with detections
    for (lib_id, rt_coef_pairs) in entry_coefficients.iter() {
        if rt_coef_pairs.is_empty() {
            continue;
        }

        // Look up the library entry by ID
        let lib_idx = match id_to_index.get(lib_id) {
            Some(&idx) => idx,
            None => continue,
        };

        let entry = &library[lib_idx];

        // Find peak boundaries from coefficient time series
        // Simple approach: find region where coefficient is above threshold
        let mut sorted_pairs = rt_coef_pairs.clone();
        sorted_pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

        let max_coef = sorted_pairs.iter().map(|(_, c)| *c).fold(0.0f64, f64::max);
        let threshold = max_coef * 0.1; // 10% of max

        // Find start, apex, end
        let above_threshold: Vec<_> = sorted_pairs.iter()
            .filter(|(_, c)| *c >= threshold)
            .collect();

        if above_threshold.is_empty() {
            continue;
        }

        let start_rt = above_threshold.first().map(|(rt, _)| *rt).unwrap_or(0.0);
        let end_rt = above_threshold.last().map(|(rt, _)| *rt).unwrap_or(0.0);
        let (apex_rt, apex_coef) = sorted_pairs.iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, c)| (*rt, *c))
            .unwrap_or((0.0, 0.0));

        // Use the first input file's ID (in a real implementation, track per-file)
        let file_id = *file_ids.values().next().unwrap_or(&1);

        // Get fragment peaks from library entry
        let mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
        let intensities: Vec<f32> = entry.fragments.iter().map(|f| f.relative_intensity).collect();

        // Compute simple score (will be replaced with proper FDR-controlled scoring)
        let score = apex_coef.min(1.0); // Placeholder score

        // Add spectrum
        let ref_id = writer.add_spectrum(
            &entry.sequence,
            &entry.modified_sequence,
            entry.precursor_mz,
            entry.charge as i32,
            apex_rt,
            start_rt,
            end_rt,
            &mzs,
            &intensities,
            score,
            file_id,
        )?;

        // Add peak boundaries
        let boundaries = osprey_core::PeakBoundaries {
            start_rt,
            end_rt,
            apex_rt,
            apex_coefficient: apex_coef,
            integrated_area: sorted_pairs.iter().map(|(_, c)| c).sum(),
            peak_quality: osprey_core::PeakQuality::default(),
        };

        let file_name = input_files.first()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        writer.add_peak_boundaries(ref_id, file_name, &boundaries)?;

        // Add run scores (placeholder - will be replaced with proper FDR)
        writer.add_run_scores(ref_id, file_name, score, apex_coef, 1.0 - score)?;

        // Optionally export coefficient time series
        if config.export_coefficients {
            for (rt, coef) in &sorted_pairs {
                // Find scan number (approximate - would need to track this properly)
                let scan = (rt * 100.0) as u32;
                writer.add_coefficient(ref_id, file_name, scan, *rt, *coef)?;
            }
        }
    }

    // Commit the batch transaction
    writer.commit()?;
    writer.finalize()?;

    log::info!("Wrote {} spectra to blib", entry_coefficients.len());

    Ok(())
}

/// Build an index of library entries by precursor m/z for fast lookup
fn build_mz_index(library: &[LibraryEntry]) -> HashMap<i32, Vec<usize>> {
    let mut index: HashMap<i32, Vec<usize>> = HashMap::new();

    for (idx, entry) in library.iter().enumerate() {
        // Round m/z to nearest integer for binning
        let mz_bin = entry.precursor_mz.round() as i32;
        // Also add to adjacent bins for tolerance
        for offset in -1..=1 {
            index.entry(mz_bin + offset).or_default().push(idx);
        }
    }

    index
}

/// Select candidate library entries for a spectrum
///
/// Candidate selection uses:
/// 1. Isolation window from mzML (defines which precursors are fragmented)
/// 2. RT tolerance (calibrated or fallback)
/// 3. Optional library filter (for calibration discovery phase)
/// 4. Optional RT calibration (converts library RT to expected measured RT)
///
/// Note: The isolation window check is sufficient for precursor filtering.
/// No additional precursor m/z tolerance is applied.
fn select_candidates_with_calibration(
    spectrum: &Spectrum,
    library: &[LibraryEntry],
    library_by_mz: &HashMap<i32, Vec<usize>>,
    rt_tolerance: f64,
    max_candidates: usize,
    library_filter: Option<&std::collections::HashSet<usize>>,
    calibration: Option<&RTCalibration>,
) -> Vec<usize> {
    let mut candidates = Vec::new();
    let mut seen = std::collections::HashSet::new();

    // Get the range of bins to search based on isolation window
    let lower_bin = spectrum.isolation_window.lower_bound().floor() as i32;
    let upper_bin = spectrum.isolation_window.upper_bound().ceil() as i32;

    // Search all bins that could contain candidates within the isolation window
    for bin in lower_bin..=upper_bin {
        if let Some(indices) = library_by_mz.get(&bin) {
            for &idx in indices {
                // Skip if already processed
                if seen.contains(&idx) {
                    continue;
                }
                seen.insert(idx);

                // Skip if not in library filter (for calibration discovery)
                if let Some(filter) = library_filter {
                    if !filter.contains(&idx) {
                        continue;
                    }
                }

                let entry = &library[idx];

                // Check if precursor falls within isolation window
                // This is the primary precursor filter - isolation window from mzML
                // defines which precursors were actually fragmented
                if !spectrum.isolation_window.contains(entry.precursor_mz) {
                    continue;
                }

                // Apply RT tolerance with optional calibration
                let expected_rt = if let Some(cal) = calibration {
                    // Convert library RT to expected measured RT using calibration
                    cal.predict(entry.retention_time)
                } else {
                    // Use library RT directly
                    entry.retention_time
                };

                // Use local tolerance if calibration available, otherwise global
                let effective_tolerance = if let Some(cal) = calibration {
                    // rt_tolerance = global_residual_sd * factor
                    // Recover factor to apply to local residual: factor = rt_tolerance / global_residual_sd
                    let factor = rt_tolerance / cal.residual_std().max(0.001);
                    cal.local_tolerance(entry.retention_time, factor, 0.25)
                } else {
                    rt_tolerance
                };

                if (expected_rt - spectrum.retention_time).abs() > effective_tolerance {
                    continue;
                }

                candidates.push(idx);
            }
        }
    }

    // Limit to max candidates (keep those closest in RT)
    if candidates.len() > max_candidates {
        candidates.sort_by(|&a, &b| {
            let rt_a = if let Some(cal) = calibration {
                (cal.predict(library[a].retention_time) - spectrum.retention_time).abs()
            } else {
                (library[a].retention_time - spectrum.retention_time).abs()
            };
            let rt_b = if let Some(cal) = calibration {
                (cal.predict(library[b].retention_time) - spectrum.retention_time).abs()
            } else {
                (library[b].retention_time - spectrum.retention_time).abs()
            };
            rt_a.total_cmp(&rt_b)
        });
        candidates.truncate(max_candidates);
    }

    candidates
}

/// Select candidate library entries for a spectrum (simple version for tests)
#[cfg(test)]
fn select_candidates(
    spectrum: &Spectrum,
    library: &[LibraryEntry],
    library_by_mz: &HashMap<i32, Vec<usize>>,
    rt_tolerance: f64,
    max_candidates: usize,
) -> Vec<usize> {
    select_candidates_with_calibration(
        spectrum,
        library,
        library_by_mz,
        rt_tolerance,
        max_candidates,
        None,
        None,
    )
}

/// Write calibration debug CSV for diagnosis of scoring results
///
/// This produces a paired CSV with target and decoy results on the same row,
/// making it easy to analyze what separates targets from decoys.
fn write_calibration_debug_csv(
    matches: &[CalibrationMatch],
    _library: &[LibraryEntry],
    output_path: &std::path::Path,
    expected_rt_fn: Option<&dyn Fn(f64) -> f64>,
) -> Result<()> {
    use std::fs::File;

    let file = File::create(output_path).map_err(|e| {
        OspreyError::OutputError(format!("Failed to create debug CSV: {}", e))
    })?;
    let mut writer = std::io::BufWriter::new(file);

    // Pair targets with their decoys
    let paired = pair_calibration_matches(matches, expected_rt_fn);

    log::info!("Writing {} paired target-decoy results to debug CSV", paired.len());

    // Write header with paired columns
    // Note: sorted by winning_evalue ascending (best matches first, lower E-value = better)
    writeln!(
        writer,
        "target_entry_id,charge,target_sequence,decoy_sequence,\
         winning_evalue,target_evalue,decoy_evalue,target_wins_evalue,\
         winning_xcorr,target_xcorr,decoy_xcorr,\
         winning_libcosine,target_libcosine,decoy_libcosine,\
         target_isotope_score,decoy_isotope_score,\
         target_precursor_error_ppm,decoy_precursor_error_ppm,\
         target_rt,decoy_rt,library_rt,expected_rt,target_delta_rt,decoy_delta_rt,\
         target_matched_frags,decoy_matched_frags,target_wins"
    ).map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

    // Write each paired result (sorted by winning_evalue ascending - best first)
    for p in &paired {
        writeln!(
            writer,
            "{},{},{},{},{:.2e},{:.2e},{:.2e},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{},{},{:.2},{:.2},{:.2},{},{:.2},{:.2},{},{},{}",
            p.target_entry_id,
            p.charge,
            p.target_sequence,
            p.decoy_sequence,
            p.winning_evalue,
            p.target_evalue,
            p.decoy_evalue,
            p.target_wins_evalue,
            p.winning_xcorr,
            p.target_xcorr,
            p.decoy_xcorr,
            p.winning_libcosine,
            p.target_libcosine,
            p.decoy_libcosine,
            p.target_isotope_score.map(|v| format!("{:.4}", v)).unwrap_or_default(),
            p.decoy_isotope_score.map(|v| format!("{:.4}", v)).unwrap_or_default(),
            p.target_precursor_error_ppm.map(|v| format!("{:.2}", v)).unwrap_or_default(),
            p.decoy_precursor_error_ppm.map(|v| format!("{:.2}", v)).unwrap_or_default(),
            p.target_rt,
            p.decoy_rt,
            p.library_rt,
            p.expected_rt.map(|v| format!("{:.2}", v)).unwrap_or_default(),
            p.target_delta_rt,
            p.decoy_delta_rt,
            p.target_matched_frags,
            p.decoy_matched_frags,
            p.target_wins
        ).map_err(|e| OspreyError::OutputError(format!("Failed to write row: {}", e)))?;
    }

    Ok(())
}

/// Write unpaired calibration debug CSV (original format for backwards compatibility)
#[allow(dead_code)]
fn write_calibration_debug_csv_unpaired(
    matches: &[CalibrationMatch],
    library: &[LibraryEntry],
    output_path: &std::path::Path,
) -> Result<()> {
    use std::fs::File;

    let file = File::create(output_path).map_err(|e| {
        OspreyError::OutputError(format!("Failed to create debug CSV: {}", e))
    })?;
    let mut writer = std::io::BufWriter::new(file);

    // Write header
    writeln!(writer, "entry_id,peptide,charge,precursor_mz,library_rt,measured_rt,score,xcorr,isotope_score,n_matched,is_decoy")
        .map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

    // Build lookup
    let id_to_entry: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    // Write each match
    for m in matches {
        if let Some(entry) = id_to_entry.get(&m.entry_id) {
            writeln!(
                writer,
                "{},{},{},{:.4},{:.2},{:.2},{:.6},{:.6},{},{},{}",
                m.entry_id,
                entry.sequence,
                entry.charge,
                entry.precursor_mz,
                entry.retention_time,
                m.measured_rt,
                m.score,
                m.xcorr_score,
                m.isotope_cosine_score.map(|v| format!("{:.4}", v)).unwrap_or_default(),
                m.n_matched_fragments,
                m.is_decoy
            ).map_err(|e| OspreyError::OutputError(format!("Failed to write row: {}", e)))?;
        }
    }

    Ok(())
}

/// Run calibration discovery using streaming pipeline
///
/// This is a memory-efficient alternative that streams spectra from the mzML file
/// instead of loading all spectra into memory at once.
///
/// Requires the `streaming` feature.
#[cfg(feature = "streaming")]
fn run_calibration_discovery_streaming(
    path: &std::path::Path,
    library: &[LibraryEntry],
    config: &OspreyConfig,
) -> Result<(RTCalibration, CalibrationParams)> {
    use tokio::runtime::Runtime;

    let rt_config = &config.rt_calibration;

    // Calculate library RT range
    let library_rts: Vec<f64> = library
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.retention_time)
        .collect();
    let min_rt = library_rts.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_rt = library_rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let rt_range = max_rt - min_rt;

    // For streaming, use 20% RT tolerance (conservative default)
    // We don't have mzML RT range until after streaming, so we assume similar scales
    let tolerance_fraction = 0.2;
    let initial_tolerance = rt_range * tolerance_fraction;

    let n_targets = library.iter().filter(|e| !e.is_decoy && !e.fragments.is_empty()).count();
    let n_decoys = library.iter().filter(|e| e.is_decoy && !e.fragments.is_empty()).count();

    log::info!(
        "Streaming calibration: scoring {} targets + {} decoys (RT range: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        min_rt,
        max_rt
    );

    // Use full library for calibration (XCorr via BLAS is fast enough)
    let calibration_library: Vec<LibraryEntry> = library.to_vec();
    log::info!(
        "Using full library for streaming calibration ({} entries)",
        calibration_library.len()
    );

    // Configure streaming pipeline
    let pipeline_config = PipelineConfig {
        num_preprocessing_threads: config.n_threads,
        channel_buffer_size: config.n_threads * 2,
        rt_tolerance: initial_tolerance,
        ..Default::default()
    };

    // Run streaming pipeline in tokio runtime
    let rt = Runtime::new().map_err(|e| {
        OspreyError::InternalError(format!("Failed to create tokio runtime: {}", e))
    })?;

    let matches = rt.block_on(async {
        osprey_scoring::pipeline::run_streaming_pipeline(path, &calibration_library, pipeline_config).await
    }).map_err(|e| {
        OspreyError::InternalError(format!("Streaming pipeline failed: {}", e))
    })?;

    log::info!("Streaming calibration: {} matches found", matches.len());

    // Count targets and decoys
    let n_target_matches = matches.iter().filter(|m| !m.is_decoy).count();
    let n_decoy_matches = matches.iter().filter(|m| m.is_decoy).count();

    log::info!(
        "Streaming results: {} targets, {} decoys with scores > 0",
        n_target_matches,
        n_decoy_matches
    );

    // Run target-decoy competition using E-value (Comet-style)
    // Lower E-value is better, so we negate for FdrController which uses "higher is better"
    let calibration_fdr = 0.01;
    let fdr_controller = FdrController::new(calibration_fdr);

    let competition_input = matches.iter().map(|m| {
        (m.clone(), -m.evalue, m.is_decoy, m.entry_id)
    });

    let competition_result = fdr_controller.compete_and_filter(competition_input);

    log::info!(
        "Competition: {} target wins, {} decoy wins",
        competition_result.n_target_wins,
        competition_result.n_decoy_wins
    );

    // Extract calibration points from passing targets
    let mut library_rts_detected: Vec<f64> = Vec::new();
    let mut measured_rts_detected: Vec<f64> = Vec::new();

    for m in &competition_result.passing_targets {
        library_rts_detected.push(m.library_rt);
        measured_rts_detected.push(m.measured_rt);
    }

    let num_confident_peptides = library_rts_detected.len();

    log::info!(
        "Streaming calibration: {} peptides at {:.0}% FDR",
        num_confident_peptides,
        calibration_fdr * 100.0
    );

    if num_confident_peptides < rt_config.min_calibration_points {
        return Err(OspreyError::ConfigError(format!(
            "Insufficient calibration points: {} < {} required",
            num_confident_peptides,
            rt_config.min_calibration_points
        )));
    }

    // Fit LOESS RT calibration
    let calibrator_config = RTCalibratorConfig {
        bandwidth: rt_config.loess_bandwidth,
        degree: 1,
        min_points: rt_config.min_calibration_points,
        robustness_iter: 2,
    };
    let calibrator = RTCalibrator::with_config(calibrator_config);
    let rt_calibration = calibrator.fit(&library_rts_detected, &measured_rts_detected)?;
    let rt_stats = rt_calibration.stats();

    // Build CalibrationParams (mass calibration not available in streaming mode)
    let calibration_params = CalibrationParams {
        metadata: CalibrationMetadata {
            num_confident_peptides,
            num_sampled_precursors: matches.len(),
            calibration_successful: true,
            timestamp: chrono::Utc::now().to_rfc3339(),
            isolation_scheme: None,
        },
        ms1_calibration: MzCalibration::uncalibrated(),
        ms2_calibration: MzCalibration::uncalibrated(),
        rt_calibration: RTCalibrationParams {
            method: RTCalibrationMethod::LOESS,
            residual_sd: rt_stats.residual_std,
            n_points: rt_stats.n_points,
            r_squared: rt_stats.r_squared,
            model_params: Some(rt_calibration.export_model_params()),
        },
    };

    Ok((rt_calibration, calibration_params))
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::IsolationWindow;

    fn make_test_spectrum() -> Spectrum {
        Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 500.0],
            intensities: vec![100.0, 200.0, 300.0],
        }
    }

    fn make_test_entry(id: u32, mz: f64, rt: f64) -> LibraryEntry {
        LibraryEntry::new(id, "PEPTIDE".into(), "PEPTIDE".into(), 2, mz, rt)
    }

    #[test]
    fn test_select_candidates() {
        let spectrum = make_test_spectrum();
        let library = vec![
            make_test_entry(0, 500.0, 10.0),  // Should match - in window and correct RT
            make_test_entry(1, 505.0, 10.0),  // In window and correct RT - should match
            make_test_entry(2, 520.0, 10.0),  // Outside isolation window
            make_test_entry(3, 500.0, 20.0),  // In window but wrong RT
        ];
        let mz_index = build_mz_index(&library);

        // With rt_tolerance=2.0, entries 0 and 1 should match (in isolation window and RT)
        let candidates = select_candidates(&spectrum, &library, &mz_index, 2.0, 100);

        assert!(candidates.contains(&0), "Entry 0 should match (in window, correct RT)");
        assert!(candidates.contains(&1), "Entry 1 should match (in window, correct RT)");
        assert!(!candidates.contains(&2), "Entry 2 should not match (outside window)");
        assert!(!candidates.contains(&3), "Entry 3 should not match (wrong RT)");
    }
}
