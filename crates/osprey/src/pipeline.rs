//! Main analysis pipeline for Osprey
//!
//! This module orchestrates the complete analysis workflow.
//!
//! ## Per-File RT Calibration Strategy
//!
//! When RT calibration is enabled (default), each file is calibrated independently:
//!
//! ### Per File: Calibration Discovery
//! 1. Check for cached calibration JSON on disk (keyed by input filename)
//! 2. If valid cached calibration exists, reuse it (fast path)
//! 3. Otherwise: use ALL library peptides (no sampling)
//! 4. Assume library RT range ≈ mzML RT range
//! 5. Wide initial tolerance (25% of gradient range)
//! 6. Detect peaks and record (library_RT, measured_apex_RT) pairs
//! 7. Fit LOESS calibration curve
//! 8. Calculate residual SD for tight tolerance
//! 9. Save calibration JSON to disk (for future reuse)
//!
//! **Rationale**: Each file may have slightly different LC conditions,
//! so independent calibration produces better accuracy. The per-file
//! JSON caching avoids redundant re-calibration on subsequent runs.
//!
//! ## Candidate Selection
//!
//! Candidate selection uses:
//! - **Isolation window**: From mzML file (defines which precursors are fragmented)
//! - **RT tolerance**: Calibrated or fallback (if calibration disabled/fails)
//!
//! Note: No additional precursor m/z tolerance is applied - the isolation window
//! from the mzML is sufficient.

use arrow::array::{ArrayRef, BooleanBuilder, Float32Builder, StringBuilder, UInt8Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use indicatif::{ProgressBar, ProgressStyle};
use osprey_chromatography::{
    calculate_mz_calibration,
    calibration_path_for_input,
    load_calibration,
    save_calibration,
    CalibrationMetadata,
    // Full calibration types
    CalibrationParams,
    IsolationScheme,
    MzCalibration,
    MzQCData,
    PeakDetector,
    RTCalibration,
    RTCalibrationMethod,
    RTCalibrationParams,
    RTCalibrator,
    RTCalibratorConfig,
};
use osprey_core::{
    BinConfig, DecoyMethod as CoreDecoyMethod, FeatureSet, FragmentToleranceConfig, LibraryEntry,
    LibraryFragment, MS1Spectrum, OspreyConfig, OspreyError, RegressionResult, Result, Spectrum,
    ToleranceUnit,
};
use osprey_fdr::{FdrController, MokapotRunner, PsmFeatures};
use osprey_io::{load_all_spectra, load_library, BlibWriter, MS1Index};
use osprey_regression::{Binner, DesignMatrixBuilder, OptimizedSolver, RidgeSolver};
use osprey_scoring::{
    batch::{
        run_coelution_calibration_scoring, sample_library_for_calibration, BatchScorer,
        MS1SpectrumLookup, PreprocessedLibrary, PreprocessedSpectra,
    },
    has_top3_fragment_match, DecoyGenerator, DecoyMethod, Enzyme, FeatureExtractor,
};

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

use osprey_scoring::batch::{pair_calibration_matches, CalibrationMatch};

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
        let (valid_targets, decoys, _stats) =
            generator.generate_all_with_collision_detection(&targets);

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

    // Ridge regression always uses unit resolution bins (~1 Th, 2000 bins).
    // HRAM precision (ppm-based matching) is applied separately in:
    //   - Fragment pre-filter (has_top3_fragment_match with ppm tolerance)
    //   - Spectral scoring (SpectralScorer with ppm tolerance)
    //   - Mass accuracy features (ppm-based error computation)
    // Using 0.02 Th bins (100K bins) for regression creates impractically large
    // dense design matrices (~400MB per spectrum per thread).
    let bin_config = BinConfig::unit_resolution();

    let binner = Binner::new(bin_config);

    // Get regularization lambda
    let lambda = match &config.regularization_lambda {
        osprey_core::RegularizationSetting::Fixed(l) => *l,
        _ => 1.0, // Default lambda for now
    };

    // Create optimized f32 solver (replaces old f64 RidgeSolver)
    let optimized_solver = OptimizedSolver::new(lambda as f32);

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

    // Process each input file with per-file calibration:
    // - Each file runs its own calibration discovery (or loads from cache)
    //
    // MEMORY-EFFICIENT DESIGN:
    // - Score immediately after each file's regression
    // - Write PIN file immediately after scoring
    // - Keep only scored entries in memory (small, ~5MB per file)
    // - Free regression results and spectra after each file
    let mut per_file_scored: Vec<(String, Vec<ScoredEntry>)> = Vec::new();
    let mut pin_files: HashMap<String, std::path::PathBuf> = HashMap::new();
    // Create mokapot output directory upfront
    let output_dir = config
        .output_blib
        .parent()
        .unwrap_or(std::path::Path::new("."));
    let mokapot_dir = output_dir.join("mokapot");
    std::fs::create_dir_all(&mokapot_dir)?;

    // Build mokapot runner for PIN file writing
    let mokapot = MokapotRunner::new().with_test_fdr(config.run_fdr);

    for (file_idx, input_file) in config.input_files.iter().enumerate() {
        log::info!(
            "===== Processing file {}/{}: {} =====",
            file_idx + 1,
            config.input_files.len(),
            input_file.display()
        );

        // Extract file name for identification
        let file_name = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Run regression for this file
        let (file_results, file_spectra, _file_preprocessed_spectra, calibration_result) =
            process_file_with_calibration(
                input_file,
                &library,
                &binner,
                &optimized_solver,
                &config,
                preprocessed_library.as_ref(),
                &batch_scorer,
                use_windowed_scoring,
            )?;

        // Extract per-file CalibrationParams for m/z correction during scoring
        let calibration_params = calibration_result.as_ref().map(|(_, _, params)| params);

        // === SCORE IMMEDIATELY after regression ===
        // This avoids keeping regression results in memory across all files
        let scored_entries = score_run(
            &library,
            &file_results,
            &file_spectra,
            &file_name,
            calibration_params,
            *binner.config(),
        )?;

        log::info!(
            "Scored {} entries ({} targets, {} decoys) for {}",
            scored_entries.len(),
            scored_entries.iter().filter(|e| !e.is_decoy).count(),
            scored_entries.iter().filter(|e| e.is_decoy).count(),
            file_name
        );

        // Deduplicate double-counted peptides sharing fragment ions within the same isolation window
        let scored_entries = deduplicate_double_counting(
            scored_entries,
            &library,
            &file_spectra,
            calibration_params,
            &config,
        );

        // Export coefficient matrix to parquet if requested
        if config.export_coefficients {
            write_coefficient_parquet(
                &scored_entries,
                &library,
                &file_results,
                &file_spectra,
                input_file,
            )?;
        }

        // === WRITE PIN FILE IMMEDIATELY ===
        // Create PSM features and write to disk, freeing memory
        let psm_features = create_psm_features_for_file(&file_name, &scored_entries, &library);
        let pin_path = mokapot.write_single_pin_file(&file_name, &psm_features, &mokapot_dir)?;
        pin_files.insert(file_name.clone(), pin_path);

        // Keep only scored entries (small) for FDR update later
        // Regression results and spectra are dropped here, freeing ~150MB per file
        per_file_scored.push((file_name, scored_entries));
    }

    // Calculate total results for logging
    let total_scored: usize = per_file_scored.iter().map(|(_, s)| s.len()).sum();
    log::info!(
        "Analysis complete. {} total scored entries across {} files",
        total_scored,
        config.input_files.len()
    );

    // Filter out empty files (mokapot can't handle empty PIN files)
    let non_empty_files: Vec<(String, Vec<ScoredEntry>)> = per_file_scored
        .into_iter()
        .filter(|(name, entries)| {
            if entries.is_empty() {
                log::warn!(
                    "Skipping file '{}' with 0 scored entries for FDR control",
                    name
                );
                false
            } else {
                true
            }
        })
        .collect();

    let non_empty_pin_files: HashMap<String, std::path::PathBuf> = pin_files
        .into_iter()
        .filter(|(name, _)| non_empty_files.iter().any(|(n, _)| n == name))
        .collect();

    if non_empty_files.is_empty() {
        return Err(OspreyError::config(
            "No scored entries found across all files. Cannot perform FDR control.",
        ));
    }

    log::info!(
        "====== Beginning Mokapot Post-Processing on {} Files ======",
        non_empty_files.len()
    );
    log::info!("Using {} pre-written PIN files", non_empty_pin_files.len());

    // Run two-level FDR control (run-level + experiment-level)
    // PIN files are already written, pass them to avoid regeneration
    let experiment_entries = run_two_level_fdr(
        &mut non_empty_files.clone(),
        &library,
        config.run_fdr,
        config.experiment_fdr,
        output_dir,
        Some(non_empty_pin_files),
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

    // Note: With memory-efficient processing, regression results are freed after scoring each file.
    // The primary output is the blib file; return empty vec for backward compatibility.
    // Callers needing raw regression results should process files individually.
    Ok(Vec::new())
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
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn process_file_with_calibration(
    path: &std::path::Path,
    library: &[LibraryEntry],
    binner: &Binner,
    optimized_solver: &OptimizedSolver,
    config: &OspreyConfig,
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
            binner,
            optimized_solver,
            config.rt_calibration.fallback_rt_tolerance,
            config.max_candidates_per_spectrum,
            config.rt_calibration.min_rt_tolerance,
            None,
            None,
            config.fragment_tolerance,
        )?;
        return Ok((results, spectra, None, None));
    }

    // Run calibration discovery with ALL peptides
    log::info!("Calibration Discovery (using all peptides)");

    // Check if a valid calibration file already exists alongside the input file
    if let Some(input_dir) = path.parent() {
        let cal_path = calibration_path_for_input(path, input_dir);
        if cal_path.exists() {
            log::info!("Found existing calibration file: {}", cal_path.display());

            match load_calibration(&cal_path) {
                Ok(cal_params) => {
                    // Check if calibration has model data for RT reconstruction
                    if cal_params.rt_calibration.has_model_data() && cal_params.is_calibrated() {
                        log::info!("Reusing existing RT calibration from file");

                        // Reconstruct RTCalibration from model params
                        let model_params = cal_params.rt_calibration.model_params.as_ref().unwrap();
                        match RTCalibration::from_model_params(
                            model_params,
                            cal_params.rt_calibration.residual_sd,
                        ) {
                            Ok(rt_cal) => {
                                let tolerance = cal_params.rt_calibration.residual_sd
                                    * config.rt_calibration.rt_tolerance_factor;
                                let tolerance =
                                    tolerance.max(config.rt_calibration.min_rt_tolerance);

                                log::info!(
                                    "Using cached calibration: {} points, R²={:.4}, tolerance={:.2} min",
                                    cal_params.rt_calibration.n_points,
                                    cal_params.rt_calibration.r_squared,
                                    tolerance
                                );

                                // Log calibration summary
                                cal_params.log_summary();

                                // Run full search with cached calibration (using calibrated MS2 m/z)
                                let results = process_spectra_optimized(
                                    &spectra,
                                    library,
                                    binner,
                                    optimized_solver,
                                    tolerance,
                                    config.max_candidates_per_spectrum,
                                    config.rt_calibration.min_rt_tolerance,
                                    Some(&rt_cal),
                                    Some(&cal_params.ms2_calibration),
                                    config.fragment_tolerance,
                                )?;

                                return Ok((
                                    results,
                                    spectra,
                                    None,
                                    Some((rt_cal, tolerance, cal_params)),
                                ));
                            }
                            Err(e) => {
                                log::warn!("Failed to reconstruct calibration from cached file: {}. Re-running calibration.", e);
                            }
                        }
                    } else {
                        log::info!(
                            "Cached calibration missing model data. Re-running calibration."
                        );
                    }
                }
                Err(e) => {
                    log::warn!(
                        "Failed to load cached calibration: {}. Re-running calibration.",
                        e
                    );
                }
            }
        }
    }

    // Run calibration using appropriate method based on mode and library size
    let calibration_result = if use_windowed_scoring {
        // Use windowed scoring for large libraries (memory-efficient)
        run_calibration_discovery_windowed(library, &spectra, &ms1_index, config)
    } else if let Some(preproc_lib) = preprocessed_library {
        // Preprocess spectra for batch scoring (reused later)
        log::info!(
            "Preprocessing {} spectra for batch scoring...",
            spectra.len()
        );
        let preprocessed_spectra = batch_scorer.preprocess_spectra(&spectra);

        run_calibration_discovery_with_cache(library, preproc_lib, &preprocessed_spectra, config)
    } else {
        // Fallback: use windowed scoring
        run_calibration_discovery_windowed(library, &spectra, &ms1_index, config)
    };

    // Determine RT tolerance and extract calibration
    let (rt_tolerance, calibration_opt, calibration_params_opt) = match calibration_result {
        Ok((rt_cal, cal_params)) => {
            // Use calibration residual SD × factor as tolerance
            let tolerance =
                cal_params.rt_calibration.residual_sd * config.rt_calibration.rt_tolerance_factor;
            let tolerance = tolerance.max(config.rt_calibration.min_rt_tolerance);
            log::info!(
                "Using calibrated RT tolerance: {:.2} min ({}× residual SD)",
                tolerance,
                config.rt_calibration.rt_tolerance_factor
            );

            // Save calibration to JSON alongside input file (for reuse)
            if let Some(input_dir) = path.parent() {
                let cal_path = calibration_path_for_input(path, input_dir);
                if let Err(e) = save_calibration(&cal_params, &cal_path) {
                    log::warn!("Failed to save calibration: {}", e);
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

    // Calibrated full search (using calibrated MS2 m/z)
    log::info!("Calibrated Full Search");

    // Get MS2 calibration if available
    let ms2_cal_ref = calibration_params_opt.as_ref().map(|p| &p.ms2_calibration);

    let results = process_spectra_optimized(
        &spectra,
        library,
        binner,
        optimized_solver,
        rt_tolerance,
        config.max_candidates_per_spectrum,
        config.rt_calibration.min_rt_tolerance,
        calibration_opt.as_ref(),
        ms2_cal_ref,
        config.fragment_tolerance,
    )?;

    // Return calibration for this file (used during scoring)
    // Note: preprocessed_spectra is only available when not using windowed scoring
    let calibration_to_share = match (calibration_opt, calibration_params_opt) {
        (Some(rt_cal), Some(cal_params)) => Some((rt_cal, rt_tolerance, cal_params)),
        _ => None,
    };
    Ok((results, spectra, None, calibration_to_share))
}

/// Run RT calibration discovery phase using BLAS-accelerated XCorr scoring
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
    let max_rt = library_rts
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let rt_range = max_rt - min_rt;

    // Get all entries with fragments (both targets and decoys)
    let entries_with_fragments: Vec<LibraryEntry> = library
        .iter()
        .filter(|e| !e.fragments.is_empty())
        .cloned()
        .collect();

    let n_targets = entries_with_fragments
        .iter()
        .filter(|e| !e.is_decoy)
        .count();
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

    log::info!(
        "Preprocessing {} spectra for batch scoring...",
        spectra.len()
    );
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
    let candidate_slope = if rt_range > 0.0 {
        meas_rt_range / rt_range
    } else {
        1.0
    };
    let use_direct_rts = (0.8..=1.2).contains(&candidate_slope);

    let (rt_slope, rt_intercept) = if use_direct_rts {
        log::info!(
            "RT ranges are similar (slope would be {:.2}), using library RTs directly",
            candidate_slope
        );
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
        "Computing XCorr scores ({} × {} = {} pairs via BLAS)...",
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
        "Calibration results: {} targets, {} decoys with scores > 0",
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
    let lib_max_rt = library_rts
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let lib_rt_range = lib_max_rt - lib_min_rt;

    // Calculate mzML RT range from spectra
    let mzml_min_rt = spectra
        .iter()
        .map(|s| s.retention_time)
        .fold(f64::INFINITY, f64::min);
    let mzml_max_rt = spectra
        .iter()
        .map(|s| s.retention_time)
        .fold(f64::NEG_INFINITY, f64::max);
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
            log::info!("RT mapping: library and mzML ranges are similar, using identity mapping");
            (1.0, 0.0, false)
        } else {
            // Linear mapping: expected_rt = slope * lib_rt + intercept
            // Maps (lib_min_rt, lib_max_rt) -> (mzml_min_rt, mzml_max_rt)
            let slope = if lib_rt_range > 0.0 {
                mzml_rt_range / lib_rt_range
            } else {
                1.0
            };
            let intercept = mzml_min_rt - slope * lib_min_rt;
            log::info!(
                "RT mapping: linear transform from library ({:.2}-{:.2}) to mzML ({:.2}-{:.2} min)",
                lib_min_rt,
                lib_max_rt,
                mzml_min_rt,
                mzml_max_rt
            );
            log::info!(
                "RT mapping: expected_rt = {:.4} * library_rt + {:.4}",
                slope,
                intercept
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

    // Calibration sampling with retry loop
    // Attempt 1: sample calibration_sample_size targets (default 5000)
    // Attempt 2: expand by retry_factor (default 3×)
    // Attempt 3: use ALL library entries (guaranteed fallback)
    // Matches accumulate across attempts — FDR runs on the combined set.
    let sample_size = rt_config.calibration_sample_size;
    let retry_factor = rt_config.calibration_retry_factor;
    let max_attempts: usize = if sample_size > 0 && retry_factor > 1.0 {
        3
    } else {
        1
    };
    let mut current_sample_size = sample_size;

    // Accumulate best match per entry across all attempts
    let mut accumulated_matches: HashMap<u32, CalibrationMatch> = HashMap::new();

    for attempt in 1..=max_attempts {
        let calibration_library =
            sample_library_for_calibration(library, current_sample_size, 42 + attempt as u64);

        log::info!(
            "Calibration attempt {}/{}: {} entries ({})",
            attempt,
            max_attempts,
            calibration_library.len(),
            if current_sample_size == 0 {
                "ALL targets".to_string()
            } else {
                format!("{} targets sampled", current_sample_size)
            }
        );

        // Run co-elution calibration scoring with fragment XIC correlation
        let new_matches = if has_ms1 {
            run_coelution_calibration_scoring(
                &calibration_library,
                spectra,
                Some(&MS1IndexWrapper(ms1_index)),
                config.fragment_tolerance,
                config.precursor_tolerance.tolerance,
                initial_tolerance,
            )
        } else {
            run_coelution_calibration_scoring::<MS1IndexWrapper>(
                &calibration_library,
                spectra,
                None,
                config.fragment_tolerance,
                config.precursor_tolerance.tolerance,
                initial_tolerance,
            )
        };

        // Accumulate matches: keep best score per entry_id across all attempts (higher is better)
        let mut n_new = 0usize;
        let mut n_improved = 0usize;
        for m in new_matches {
            match accumulated_matches.entry(m.entry_id) {
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(m);
                    n_new += 1;
                }
                std::collections::hash_map::Entry::Occupied(mut e) => {
                    if m.score > e.get().score {
                        e.insert(m);
                        n_improved += 1;
                    }
                }
            }
        }

        if attempt > 1 {
            log::info!(
                "Accumulated {} new entries, {} improved ({} total unique entries)",
                n_new,
                n_improved,
                accumulated_matches.len()
            );
        }

        // Write debug CSV
        let debug_dir = config
            .output_blib
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .map(|p| p.to_path_buf())
            .or_else(|| {
                config
                    .input_files
                    .first()
                    .and_then(|f| f.parent())
                    .map(|p| p.to_path_buf())
            })
            .unwrap_or_else(|| std::path::PathBuf::from("."));

        let debug_path = debug_dir.join("calibration_debug.csv");

        let linear_rt_mapping = |lib_rt: f64| -> f64 { rt_slope * lib_rt + rt_intercept };
        let expected_rt_fn: Option<&dyn Fn(f64) -> f64> = Some(&linear_rt_mapping);

        let matches_owned: Vec<CalibrationMatch> = accumulated_matches.values().cloned().collect();
        if let Err(e) = write_calibration_debug_csv(
            &matches_owned,
            &calibration_library,
            &debug_path,
            expected_rt_fn,
        ) {
            log::warn!("Failed to write calibration debug CSV: {}", e);
        } else {
            log::debug!(
                "Wrote paired calibration debug CSV to: {}",
                debug_path.display()
            );
        }

        // LDA-based scoring on accumulated matches
        // Only use isotope feature for HRAM mode (isotopes not separated at unit resolution)
        let use_isotope_feature =
            matches!(config.resolution_mode, osprey_core::ResolutionMode::HRAM) && has_ms1;
        let calibration_fdr = 0.01;

        // Convert to Vec for LDA training
        let mut all_matches: Vec<CalibrationMatch> =
            accumulated_matches.values().cloned().collect();

        // Train LDA and score
        let _n_passing = osprey_scoring::calibration_ml::train_and_score_calibration(
            &mut all_matches,
            use_isotope_feature,
        )?;

        // Filter to passing targets (q-value <= 1% FDR, not decoys)
        // CRITICAL: Also require minimum S/N for RT calibration quality
        // LDA may select peptides with good spectral scores but poor peak quality
        const MIN_SNR_FOR_RT_CAL: f64 = 5.0; // Require reasonable peak signal-to-noise

        let passing_targets: Vec<&CalibrationMatch> = all_matches
            .iter()
            .filter(|m| {
                !m.is_decoy
                    && m.q_value <= calibration_fdr
                    && m.signal_to_noise >= MIN_SNR_FOR_RT_CAL
            })
            .collect();

        // Count wins for logging
        let n_target_wins = all_matches
            .iter()
            .filter(|m| !m.is_decoy && m.q_value <= calibration_fdr)
            .count();
        let n_decoy_wins = all_matches
            .iter()
            .filter(|m| m.is_decoy && m.q_value <= calibration_fdr)
            .count();

        log::info!(
            "Competition: {} target wins, {} decoy wins",
            n_target_wins,
            n_decoy_wins
        );

        // Log S/N filter impact
        if passing_targets.len() < n_target_wins {
            log::info!(
                "  RT quality filter: {} → {} peptides (removed {} with S/N < {:.1})",
                n_target_wins,
                passing_targets.len(),
                n_target_wins - passing_targets.len(),
                MIN_SNR_FOR_RT_CAL
            );
        }

        // Extract calibration points + mass errors from passing targets
        let mut library_rts_detected: Vec<f64> = Vec::new();
        let mut measured_rts_detected: Vec<f64> = Vec::new();
        let mut mz_qc_data = MzQCData::new(config.fragment_tolerance.unit);

        for m in &passing_targets {
            library_rts_detected.push(m.library_rt);
            measured_rts_detected.push(m.measured_rt);

            if let Some(ms1_error) = m.ms1_error {
                mz_qc_data.add_ms1_error(ms1_error);
            }
            for &ms2_error in &m.ms2_mass_errors {
                mz_qc_data.add_ms2_error(ms2_error);
            }
        }

        let num_confident_peptides = library_rts_detected.len();

        log::info!(
            "Calibration: {} peptides at {:.0}% FDR (from {} target wins, {} decoy wins)",
            num_confident_peptides,
            calibration_fdr * 100.0,
            n_target_wins,
            n_decoy_wins
        );

        // Check if we have enough calibration points
        // Absolute minimum for LOESS to work (regardless of config setting)
        const ABSOLUTE_MIN_CALIBRATION_POINTS: usize = 50;

        if num_confident_peptides < rt_config.min_calibration_points {
            if attempt < max_attempts {
                // Not the final attempt - retry with more targets
                // Determine next sample size
                if attempt + 1 == max_attempts {
                    // Final attempt always uses ALL library entries
                    current_sample_size = 0;
                } else {
                    let new_size = (current_sample_size as f64 * retry_factor) as usize;
                    let n_total_targets = library.iter().filter(|e| !e.is_decoy).count();
                    if new_size >= n_total_targets {
                        current_sample_size = 0;
                    } else {
                        current_sample_size = new_size;
                    }
                }
                log::warn!(
                    "Calibration attempt {} found only {} confident peptides (need {}). \
                     Retrying with {} targets...",
                    attempt,
                    num_confident_peptides,
                    rt_config.min_calibration_points,
                    if current_sample_size == 0 {
                        "ALL".to_string()
                    } else {
                        current_sample_size.to_string()
                    }
                );
                continue;
            } else if num_confident_peptides >= ABSOLUTE_MIN_CALIBRATION_POINTS {
                // Final attempt: Use what we have if >= absolute minimum
                log::warn!(
                    "Calibration: Using {} peptides (below target of {} but above minimum of {})",
                    num_confident_peptides,
                    rt_config.min_calibration_points,
                    ABSOLUTE_MIN_CALIBRATION_POINTS
                );
                // Continue to calibration with reduced point count
            } else {
                // Final attempt: Not enough points even for absolute minimum
                return Err(OspreyError::ConfigError(format!(
                    "Insufficient calibration points: {} < {} absolute minimum (after {} attempts)",
                    num_confident_peptides, ABSOLUTE_MIN_CALIBRATION_POINTS, attempt
                )));
            }
        }

        // Fit LOESS RT calibration
        // Use actual number of points or config minimum, whichever is smaller
        let effective_min_points = num_confident_peptides.min(rt_config.min_calibration_points);
        let calibrator_config = RTCalibratorConfig {
            bandwidth: rt_config.loess_bandwidth,
            degree: 1,
            min_points: effective_min_points,
            robustness_iter: 2,
        };
        let calibrator = RTCalibrator::with_config(calibrator_config);
        let rt_calibration = calibrator.fit(&library_rts_detected, &measured_rts_detected)?;
        let rt_stats = rt_calibration.stats();

        // Calculate mass calibration
        let (ms1_calibration, ms2_calibration) = calculate_mz_calibration(&mz_qc_data);

        // Extract isolation window scheme from spectra
        let isolation_scheme = extract_isolation_scheme(spectra);

        // Build full CalibrationParams
        let calibration_params = CalibrationParams {
            metadata: CalibrationMetadata {
                num_confident_peptides,
                num_sampled_precursors: accumulated_matches.len(),
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

        return Ok((rt_calibration, calibration_params));
    }

    // Should not reach here (loop always returns or errors)
    unreachable!("Calibration retry loop exited without result")
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
    let max_rt = library_rts
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let rt_range = max_rt - min_rt;

    let n_targets = library
        .iter()
        .filter(|e| !e.is_decoy && !e.fragments.is_empty())
        .count();
    let n_decoys = library
        .iter()
        .filter(|e| e.is_decoy && !e.fragments.is_empty())
        .count();

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
    let candidate_slope = if rt_range > 0.0 {
        meas_rt_range / rt_range
    } else {
        1.0
    };
    let use_direct_rts = (0.8..=1.2).contains(&candidate_slope);

    let (rt_slope, rt_intercept) = if use_direct_rts {
        log::info!(
            "RT ranges are similar (slope would be {:.2}), using library RTs directly",
            candidate_slope
        );
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
        "Computing XCorr scores ({} × {} = {} pairs via BLAS, using cached matrices)...",
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
    let id_to_entry: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

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
        "Calibration results: {} targets, {} decoys with scores > 0",
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
            num_confident_peptides, rt_config.min_calibration_points
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
/// - Applies MS2 calibration to spectra before binning (if available)
#[allow(clippy::too_many_arguments)]
fn process_spectra_optimized(
    spectra: &[Spectrum],
    library: &[LibraryEntry],
    binner: &Binner,
    solver: &OptimizedSolver,
    rt_tolerance: f64,
    max_candidates: usize,
    min_rt_tolerance: f64,
    calibration: Option<&RTCalibration>,
    ms2_calibration: Option<&MzCalibration>,
    base_fragment_tolerance: FragmentToleranceConfig,
) -> Result<Vec<RegressionResult>> {
    use ndarray::{Array1, Array2, ShapeBuilder};
    use osprey_chromatography::calibration::apply_spectrum_calibration;

    // Apply MS2 m/z calibration to spectra if available
    // This shifts each centroid by the mean offset so the error distribution is centered at 0
    let calibrated_spectra: Vec<Spectrum> = if let Some(ms2_cal) = ms2_calibration {
        if ms2_cal.calibrated {
            // Correction is: corrected = observed - mean
            // So if mean = -0.007, we add 0.007 to shift peaks UP
            log::info!(
                "Applying MS2 calibration: mean error = {:.4} {} → correcting by {:+.4} {} (using 3×SD = {:.3} {} tolerance)",
                ms2_cal.mean,
                ms2_cal.unit,
                -ms2_cal.mean,  // Show actual correction applied
                ms2_cal.unit,
                3.0 * ms2_cal.sd,
                ms2_cal.unit
            );
            spectra
                .iter()
                .map(|s| apply_spectrum_calibration(s, ms2_cal))
                .collect()
        } else {
            spectra.to_vec()
        }
    } else {
        spectra.to_vec()
    };

    // Determine fragment tolerance for top-3 pre-filter
    // Use calibrated 3×SD when available, otherwise use base tolerance
    let effective_fragment_tolerance = if let Some(ms2_cal) = ms2_calibration {
        if ms2_cal.calibrated {
            let tol_3sd = 3.0 * ms2_cal.sd;
            let unit = if ms2_cal.unit == "Th" {
                ToleranceUnit::Mz
            } else {
                ToleranceUnit::Ppm
            };
            // Apply minimum tolerance floor
            let min_tol = if unit == ToleranceUnit::Mz { 0.05 } else { 1.0 };
            FragmentToleranceConfig {
                tolerance: tol_3sd.max(min_tol),
                unit,
            }
        } else {
            base_fragment_tolerance
        }
    } else {
        base_fragment_tolerance
    };

    let n_bins = binner.n_bins();

    // Build RT-sorted m/z index (pre-computes expected_rt from calibration)
    let mz_rt_index = MzRTIndex::build(library, calibration);

    // Search tolerance = rt_tolerance (3× residual_SD) applied uniformly across RT range
    let search_tolerance = rt_tolerance;

    log::info!(
        "Regression: {:.4} m/z bins, {} bins, {} spectra (search_tol={:.2} min)",
        binner.config().bin_width,
        n_bins,
        calibrated_spectra.len(),
        search_tolerance
    );

    // Set up progress bar
    let pb = ProgressBar::new(calibrated_spectra.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} spectra",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Diagnostic counters: verify candidate selection respects RT tolerance
    let rt_violations = std::sync::atomic::AtomicU64::new(0);
    let rt_violation_max = std::sync::atomic::AtomicU64::new(0); // stored as f64 bits

    // Process spectra in parallel (using calibrated m/z values)
    // Each spectrum: bin observed → filter candidates → bin candidates on-the-fly → solve
    let results: Vec<RegressionResult> = calibrated_spectra
        .par_iter()
        .filter_map(|spectrum| {
            // Select candidates using RT-sorted binary search
            let candidates = select_candidates_with_calibration(
                spectrum,
                library,
                &mz_rt_index,
                rt_tolerance,
                max_candidates,
                min_rt_tolerance,
                None,
                calibration,
                Some(effective_fragment_tolerance),
                search_tolerance,
            );

            if candidates.is_empty() {
                pb.inc(1);
                return None;
            }

            // Diagnostic: verify all candidates are within RT tolerance
            for &idx in &candidates {
                let entry = &library[idx];
                let expected_rt = if let Some(cal) = calibration {
                    cal.predict(entry.retention_time)
                } else {
                    entry.retention_time
                };
                let rt_diff = (expected_rt - spectrum.retention_time).abs();
                if rt_diff > search_tolerance + 0.01 {
                    let count = rt_violations.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    // Track max violation
                    let bits = rt_diff.to_bits();
                    rt_violation_max.fetch_max(bits, std::sync::atomic::Ordering::Relaxed);
                    // Log first 10 violations with full detail
                    if count < 10 {
                        log::warn!(
                            "RT violation #{}: {} (idx={}, id={}) expected_rt={:.4} lib_rt={:.4} spectrum_rt={:.4} diff={:.4} tol={:.4} scan={}",
                            count + 1,
                            entry.modified_sequence,
                            idx,
                            entry.id,
                            expected_rt,
                            entry.retention_time,
                            spectrum.retention_time,
                            rt_diff,
                            search_tolerance,
                            spectrum.scan_number
                        );
                    }
                }
            }

            // Bin observed spectrum
            let observed_dense = binner.bin_spectrum_dense_f32(spectrum);
            let observed_vec = Array1::from_vec(observed_dense);

            // Build design matrix in column-major (Fortran) order for BLAS-optimal column access.
            // The CD-NNLS solver accesses columns repeatedly (col.dot(&residual), residual updates),
            // so contiguous columns give cache-friendly BLAS sdot/saxpy operations.
            let n_candidates = candidates.len();
            let mut design_matrix = Array2::<f32>::zeros((n_bins, n_candidates).f());
            let mut library_ids = Vec::with_capacity(n_candidates);

            for (col, &idx) in candidates.iter().enumerate() {
                let entry = &library[idx];
                // Write directly into contiguous column slice (column-major = contiguous columns)
                let mut col_view = design_matrix.column_mut(col);
                let col_slice = col_view
                    .as_slice_mut()
                    .expect("column-major column must be contiguous");
                binner.bin_library_entry_into(entry, col_slice);
                library_ids.push(entry.id);
            }

            let observed_norm = observed_vec.dot(&observed_vec);

            // Solve regression with f32
            let coefficients = match solver.solve_nonnegative(&design_matrix, &observed_vec, None) {
                Ok(c) => c,
                Err(_) => {
                    pb.inc(1);
                    return None;
                }
            };

            let coefficient_sum: f32 = coefficients.iter().sum();

            // Compute residual
            let predicted = design_matrix.dot(&coefficients);
            let diff = &observed_vec - &predicted;
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

    // Report RT tolerance violations (diagnostic)
    let n_violations = rt_violations.load(std::sync::atomic::Ordering::Relaxed);
    if n_violations > 0 {
        let max_bits = rt_violation_max.load(std::sync::atomic::Ordering::Relaxed);
        let max_diff = f64::from_bits(max_bits);
        log::warn!(
            "RT TOLERANCE VIOLATION: {} candidates selected outside search_tolerance ({:.2} min). Max RT diff: {:.2} min. Bug in candidate selection!",
            n_violations, search_tolerance, max_diff
        );
    } else {
        log::info!(
            "RT tolerance check: all candidates within {:.2} min tolerance (no violations)",
            search_tolerance
        );
    }

    log::info!(
        "Generated {} regression results with non-zero coefficients",
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
    /// mzML scan number at coefficient apex
    apex_scan_number: u32,
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
    /// Fragment-based peak boundaries from co-eluting XICs (start_rt, end_rt, fwhm)
    fragment_peak_bounds: Option<(f64, f64, f64)>,
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
    bin_config: BinConfig,
) -> Result<Vec<ScoredEntry>> {
    use osprey_chromatography::calibration::{apply_spectrum_calibration, calibrated_tolerance};
    use osprey_core::ToleranceUnit;
    use osprey_scoring::{RegressionContext, SpectralScorer};

    // Build ID-to-index map and check for duplicate IDs
    let id_to_index: HashMap<u32, usize> = library
        .iter()
        .enumerate()
        .map(|(idx, entry)| (entry.id, idx))
        .collect();

    // Diagnostic: check if any IDs are duplicated (which would cause id_to_index to lose entries)
    if id_to_index.len() != library.len() {
        let mut id_counts: HashMap<u32, Vec<usize>> = HashMap::new();
        for (idx, entry) in library.iter().enumerate() {
            id_counts.entry(entry.id).or_default().push(idx);
        }
        let duplicates: Vec<_> = id_counts
            .iter()
            .filter(|(_, indices)| indices.len() > 1)
            .take(10)
            .collect();
        log::error!(
            "DUPLICATE IDS: library has {} entries but only {} unique IDs ({} collisions). First duplicates:",
            library.len(), id_to_index.len(), library.len() - id_to_index.len()
        );
        for (id, indices) in &duplicates {
            let entries_info: Vec<String> = indices
                .iter()
                .map(|&idx| {
                    format!(
                        "idx={} seq={} rt={:.3}",
                        idx, library[idx].modified_sequence, library[idx].retention_time
                    )
                })
                .collect();
            log::error!("  id={}: [{}]", id, entries_info.join(", "));
        }
    } else {
        log::info!("ID uniqueness check: all {} IDs are unique", library.len());
    }

    // Aggregate results by library entry ID
    // Store both (RT, coef) pairs and indices to RegressionResults for context building
    let mut entry_data: HashMap<u32, Vec<(f64, f64)>> = HashMap::new(); // (RT, coef)
    let mut entry_result_indices: HashMap<u32, Vec<usize>> = HashMap::new(); // indices into results
    let mut aggregation_violations = 0u64;
    for (result_idx, result) in results.iter().enumerate() {
        for (lib_id, coef) in result.library_ids.iter().zip(result.coefficients.iter()) {
            // Diagnostic: verify that each lib_id's expected_rt is near this result's RT
            if let Some(&lib_idx) = id_to_index.get(lib_id) {
                let entry = &library[lib_idx];
                let rt_diff = (entry.retention_time - result.retention_time).abs();
                if rt_diff > 5.0 && aggregation_violations < 20 {
                    aggregation_violations += 1;
                    log::warn!(
                        "AGGREGATION MISMATCH #{}: lib_id={} ({}) lib_rt={:.3} result_rt={:.3} diff={:.2} coef={:.4} scan={}",
                        aggregation_violations, lib_id, entry.modified_sequence,
                        entry.retention_time, result.retention_time, rt_diff, coef, result.scan_number
                    );
                }
            }
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
    if aggregation_violations > 0 {
        // Count total
        let _total_violations: u64 = results
            .iter()
            .flat_map(|r| r.library_ids.iter().zip(r.coefficients.iter()))
            .filter(|(lib_id, _coef)| {
                if let Some(&_lib_idx) = id_to_index.get(lib_id) {
                    // Use a wider tolerance to check: can't compare to search_tolerance here but
                    // 5 min should be well beyond any legitimate tolerance
                    false // placeholder: just counted above
                } else {
                    false
                }
            })
            .count() as u64;
        log::warn!(
            "Total aggregation mismatches (|lib_rt - result_rt| > 5 min): {} (showed first 20)",
            aggregation_violations
        );
    }

    // Extract features and create scored entries (parallelized)
    let feature_extractor = FeatureExtractor::new();

    // Configure SpectralScorer with calibrated tolerance if available
    // Use 3×SD from MS2 calibration as the tolerance for fragment matching
    // This is unit-aware: uses Th for unit resolution, ppm for HRAM
    let (spectral_scorer, mass_accuracy_unit) = if let Some(cal) = calibration {
        // Get calibrated tolerance in the appropriate unit
        let (tol_value, tol_unit) = calibrated_tolerance(
            &cal.ms2_calibration,
            20.0, // default 20 ppm for HRAM
            ToleranceUnit::Ppm,
        );

        let scorer = match tol_unit {
            ToleranceUnit::Mz => {
                log::debug!(
                    "Using calibrated fragment tolerance: {:.4} Th (3×SD)",
                    tol_value
                );
                SpectralScorer::with_bin_config(bin_config).with_tolerance_da(tol_value)
            }
            ToleranceUnit::Ppm => {
                log::debug!(
                    "Using calibrated fragment tolerance: {:.2} ppm (3×SD)",
                    tol_value
                );
                SpectralScorer::with_bin_config(bin_config).with_tolerance_ppm(tol_value)
            }
        };
        (scorer, tol_unit)
    } else {
        // Default: unit resolution (Th)
        (
            SpectralScorer::with_bin_config(bin_config),
            ToleranceUnit::Mz,
        )
    };

    // Pre-calibrate all spectra once upfront (if calibration available)
    // This avoids cloning/calibrating spectra inside the per-peptide loop
    let calibrated_spectra: Vec<Spectrum> = if let Some(cal) = calibration {
        log::debug!(
            "Pre-calibrating {} spectra with MS2 calibration",
            spectra.len()
        );
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

    // Set up progress bar for scoring
    let pb = ProgressBar::new(lib_ids.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} peptides",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    let scored_entries: Vec<ScoredEntry> = lib_ids
        .par_iter()
        .filter_map(|lib_id| {
            let result = (|| {
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

                // Find apex RT (RT with maximum coefficient) and corresponding scan number
                let apex_rt = sorted_data
                    .iter()
                    .max_by(|a, b| a.1.total_cmp(&b.1))
                    .map(|(rt, _)| *rt)
                    .unwrap_or(entry.retention_time);

                // Get the mzML scan number at the apex (data_points and result_indices are in lockstep)
                let apex_scan_number = {
                    let (apex_idx, _) = data_points
                        .iter()
                        .enumerate()
                        .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                        .unwrap();
                    let result_indices = entry_result_indices.get(lib_id).unwrap();
                    results[result_indices[apex_idx]].scan_number
                };

                // Get RT range from coefficient series for peak region filtering
                let rt_min = rt_coef_pairs
                    .iter()
                    .map(|(rt, _)| *rt)
                    .fold(f64::INFINITY, f64::min);
                let rt_max = rt_coef_pairs
                    .iter()
                    .map(|(rt, _)| *rt)
                    .fold(f64::NEG_INFINITY, f64::max);

                // Filter pre-calibrated spectra to those in the peak region:
                // 1. Isolation window contains the peptide's precursor m/z
                // 2. RT is within the coefficient series range
                let peak_region_spectra: Vec<&Spectrum> = calibrated_spectra
                    .iter()
                    .filter(|s| {
                        s.contains_precursor(entry.precursor_mz)
                            && s.retention_time >= rt_min
                            && s.retention_time <= rt_max
                    })
                    .collect();

                // Build regression context from the RegressionResults that contain this peptide
                let regression_context = entry_result_indices.get(lib_id).map(|indices| {
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

                // Fragment co-elution correlation: correlate raw fragment XICs against
                // the regression coefficient time series (deconvolved elution profile)
                let (coelution_sum, coelution_min, n_coeluting) =
                    osprey_scoring::compute_fragment_coelution(
                        &entry.fragments,
                        &rt_coef_pairs,
                        &peak_region_spectra,
                        spectral_scorer.tolerance_da(),
                        spectral_scorer.tolerance_ppm(),
                    );
                features.fragment_coelution_sum = coelution_sum;
                features.fragment_coelution_min = coelution_min;
                features.n_coeluting_fragments = n_coeluting;

                // Compute fragment-based FWHM from co-eluting XICs for blib boundaries
                let fragment_peak_bounds = osprey_scoring::compute_fragment_fwhm(
                    &entry.fragments,
                    &rt_coef_pairs,
                    &peak_region_spectra,
                    spectral_scorer.tolerance_da(),
                    spectral_scorer.tolerance_ppm(),
                );

                // Find apex spectrum from peak region (closest to apex RT) for spectral scoring
                let apex_spectrum = peak_region_spectra.iter().min_by(|a, b| {
                    (a.retention_time - apex_rt)
                        .abs()
                        .partial_cmp(&(b.retention_time - apex_rt).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

                // Per-fragment mass accuracy at apex spectrum
                if let Some(spectrum) = apex_spectrum {
                    let matches = spectral_scorer.match_fragments(spectrum, entry);
                    let (acc_mean, acc_std) =
                        osprey_scoring::compute_mass_accuracy(&matches, mass_accuracy_unit);
                    features.mass_accuracy_mean = acc_mean;
                    features.mass_accuracy_std = acc_std;
                }

                // Primary score is LibCosine from extract_with_deconvolution (already computed)
                // Fall back to peak_apex if no spectral score
                let score = if features.dot_product > 0.0 {
                    features.dot_product
                } else {
                    features.peak_apex * 0.1 // Scale down chromatographic score
                };

                // Check if this is a decoy (high bit set in ID)
                let is_decoy = *lib_id & 0x80000000 != 0;

                // Generate PSM ID for matching with Mokapot results
                let psm_id = format!("{}_{}", file_name, lib_id);

                // Log extreme RT deviations for debugging
                if features.rt_deviation.abs() > 5.0 && !is_decoy {
                    let rt_min_data = sorted_data.iter().map(|(rt, _)| *rt).fold(f64::INFINITY, f64::min);
                    let rt_max_data = sorted_data.iter().map(|(rt, _)| *rt).fold(f64::NEG_INFINITY, f64::max);
                    let coef_max = sorted_data.iter().map(|(_, c)| *c).fold(f64::NEG_INFINITY, f64::max);
                    log::warn!(
                        "RT outlier: {} (id={}, lib_idx={}) apex_rt={:.3} lib_rt={:.3} rt_dev={:.2} scan={} n_scans={} coef_series_rt=[{:.3}..{:.3}] max_coef={:.4} mass_acc={:.1}",
                        entry.modified_sequence, lib_id, lib_idx, apex_rt,
                        entry.retention_time, features.rt_deviation,
                        apex_scan_number, sorted_data.len(),
                        rt_min_data, rt_max_data, coef_max,
                        features.mass_accuracy_mean
                    );
                    // For first few outliers, log all data points
                    if sorted_data.len() <= 10 {
                        for (i, (rt, coef)) in sorted_data.iter().enumerate() {
                            log::warn!(
                                "  data[{}]: rt={:.4} coef={:.6}",
                                i, rt, coef
                            );
                        }
                    }
                }

                Some(ScoredEntry {
                    lib_idx,
                    psm_id,
                    file_name: file_name.to_string(),
                    apex_scan_number,
                    rt_coef_pairs,
                    features,
                    score,
                    run_qvalue: 1.0,        // Will be set by Mokapot
                    experiment_qvalue: 1.0, // Will be set by experiment-level Mokapot
                    pep: 1.0,               // Will be set by Mokapot
                    is_decoy,
                    fragment_peak_bounds,
                })
            })();
            pb.inc(1);
            result
        })
        .collect();

    pb.finish_with_message("Done");

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
        .map(|entry| {
            let lib_entry = &library[entry.lib_idx];

            // Use library index as scan_number for mokapot fold splitting
            // In DIA/Osprey, each PSM represents a peptide detection across a chromatographic
            // peak (many spectra), not a single spectrum. The library index uniquely identifies
            // each precursor (peptide + charge), ensuring fold splitting groups by precursor
            // to avoid data leakage during cross-validation.
            PsmFeatures {
                psm_id: entry.psm_id.clone(),
                peptide: lib_entry.modified_sequence.clone(),
                proteins: lib_entry.protein_ids.clone(),
                scan_number: entry.apex_scan_number,
                file_name: file_name.to_string(),
                charge: lib_entry.charge,
                is_decoy: entry.is_decoy,
                features: entry.features.clone(),
                initial_score: Some(entry.score),
            }
        })
        .collect();

    // Create replicate-specific output directory for mokapot
    // Structure: mokapot/replicate_{file_name}/
    let replicate_dir = output_dir.join(format!("replicate_{}", file_name));
    std::fs::create_dir_all(&replicate_dir)?;

    // Write PIN file
    let pin_file = replicate_dir.join(format!("{}.pin", file_name));
    mokapot.write_pin(&psm_features, &pin_file)?;

    log::info!(
        "Wrote {} precursors to PIN file: {}",
        psm_features.len(),
        pin_file.display()
    );

    // Run mokapot (outputs go to replicate subfolder)
    let results = mokapot.run(&pin_file, &replicate_dir)?;

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

/// Run two-level FDR control using Mokapot with multi-file support
///
/// This implements the two-level FDR strategy using mokapot CLI:
/// 1. Run-level FDR: mokapot WITHOUT --aggregate + WITH --save_models
///    → trains ONE model from all files, reports per-file results
/// 2. Experiment-level FDR: mokapot WITH --aggregate + WITH --load_models
///    → reuses model from Step 1, reports combined experiment results
///
/// # Arguments
/// * `per_file_results` - Per-file scored entries (will be updated with q-values)
/// * `library` - Full library for reference
/// * `run_fdr` - Run-level FDR threshold
/// * `experiment_fdr` - Experiment-level FDR threshold
/// * `output_dir` - Output directory for mokapot files
/// * `pin_files` - Optional pre-written PIN file paths. If provided, skips PIN file generation.
///
/// Returns all entries with both run-level and experiment-level q-values populated.
fn run_two_level_fdr(
    per_file_results: &mut [(String, Vec<ScoredEntry>)],
    library: &[LibraryEntry],
    run_fdr: f64,
    experiment_fdr: f64,
    output_dir: &std::path::Path,
    pin_files: Option<HashMap<String, std::path::PathBuf>>,
) -> Result<Vec<ScoredEntry>> {
    // Build mokapot runner
    let mokapot = MokapotRunner::new().with_test_fdr(run_fdr);

    let mokapot_available = mokapot.is_available();
    if !mokapot_available {
        log::warn!("Mokapot not available, using simple target-decoy competition");
        log::warn!("Install mokapot with: pip install mokapot");
    }

    // Create output directory for mokapot files
    let mokapot_dir = output_dir.join("mokapot");
    std::fs::create_dir_all(&mokapot_dir)?;

    // For single file, use simpler approach
    if per_file_results.len() == 1 {
        let (file_name, entries) = per_file_results.iter_mut().next().unwrap();

        if mokapot_available {
            compute_fdr_with_mokapot(library, entries, file_name, &mokapot_dir, run_fdr, run_fdr)?;
        } else {
            apply_simple_fdr(entries, run_fdr)?;
        }

        // Report single file statistics
        report_file_statistics(file_name, entries, library, run_fdr);

        // Copy run-level to experiment-level for single file
        log::info!("Single file - using run-level results as experiment-level");
        let mut experiment_entries: Vec<ScoredEntry> = entries.clone();
        for entry in experiment_entries.iter_mut() {
            entry.experiment_qvalue = entry.run_qvalue;
        }
        return Ok(experiment_entries);
    }

    // ========== Multi-file analysis with two-step mokapot ==========
    if !mokapot_available {
        // Fallback: simple target-decoy competition
        log::info!("Running simple target-decoy competition (mokapot not available)");

        for (file_name, entries) in per_file_results.iter_mut() {
            apply_simple_fdr(entries, run_fdr)?;
            report_file_statistics(file_name, entries, library, run_fdr);
        }

        // Combine for experiment level
        let mut experiment_entries: Vec<ScoredEntry> = Vec::new();
        for (_, entries) in per_file_results.iter() {
            for entry in entries.iter() {
                experiment_entries.push(entry.clone());
            }
        }
        apply_simple_fdr_experiment(&mut experiment_entries, experiment_fdr)?;
        report_experiment_statistics(&experiment_entries, library, experiment_fdr);
        return Ok(experiment_entries);
    }

    // ========== Run two-step mokapot analysis ==========
    // Use pre-written PIN files if provided, otherwise generate them
    let pin_files = if let Some(files) = pin_files {
        files
    } else {
        // Fallback: collect PSM features and write PIN files
        log::info!("Generating PIN files from scored entries...");
        let psms_by_file = collect_psm_features_by_file(per_file_results, library);
        let files = mokapot.write_pin_files(&psms_by_file, &mokapot_dir)?;
        log::info!("Wrote {} PIN files", files.len());
        files
    };

    // Run two-step analysis: trains ONE model, gets both per-file and experiment results
    let (per_file_mokapot_results, experiment_mokapot_results) =
        mokapot.run_two_step_analysis(&pin_files, &mokapot_dir)?;

    // ========== Update per-file entries with Step 1 (run-level) results ==========
    for (file_name, entries) in per_file_results.iter_mut() {
        if let Some(results) = per_file_mokapot_results.get(file_name) {
            let result_map: HashMap<String, (f64, f64)> = results
                .iter()
                .map(|r| (r.psm_id.clone(), (r.q_value, r.pep)))
                .collect();

            for entry in entries.iter_mut() {
                if let Some(&(q_value, pep)) = result_map.get(&entry.psm_id) {
                    entry.run_qvalue = q_value;
                    entry.pep = pep;
                }
            }
        }
    }

    // Per-file results already reported by mokapot.rs, skip redundant output here

    // ========== Build experiment entries with Step 2 results ==========

    // Build PSM ID to q-value map from experiment results
    let experiment_result_map: HashMap<String, (f64, f64)> = experiment_mokapot_results
        .iter()
        .map(|r| (r.psm_id.clone(), (r.q_value, r.pep)))
        .collect();

    // Collect all entries with experiment-level q-values
    let mut experiment_entries: Vec<ScoredEntry> = Vec::new();
    for (_, entries) in per_file_results.iter() {
        for entry in entries.iter() {
            let mut entry_copy = entry.clone();

            if entry.is_decoy {
                entry_copy.experiment_qvalue = 1.0;
            } else if let Some(&(q_value, _pep)) = experiment_result_map.get(&entry.psm_id) {
                entry_copy.experiment_qvalue = q_value;
            } else {
                entry_copy.experiment_qvalue = 1.0;
            }

            experiment_entries.push(entry_copy);
        }
    }

    Ok(experiment_entries)
}

/// Report per-file statistics (precursors and peptides)
fn report_file_statistics(
    file_name: &str,
    entries: &[ScoredEntry],
    library: &[LibraryEntry],
    fdr: f64,
) {
    let passing_entries: Vec<_> = entries
        .iter()
        .filter(|e| !e.is_decoy && e.run_qvalue <= fdr)
        .collect();

    let precursor_count = passing_entries.len();

    let unique_peptides: std::collections::HashSet<_> = passing_entries
        .iter()
        .map(|e| &library[e.lib_idx].modified_sequence)
        .collect();
    let peptide_count = unique_peptides.len();

    log::info!(
        "  {}: {} precursors, {} peptides",
        file_name,
        precursor_count,
        peptide_count
    );
}

/// Report experiment-level statistics (precursors and peptides)
fn report_experiment_statistics(entries: &[ScoredEntry], library: &[LibraryEntry], fdr: f64) {
    let passing_entries: Vec<_> = entries
        .iter()
        .filter(|e| !e.is_decoy && e.experiment_qvalue <= fdr)
        .collect();

    let precursor_count = passing_entries.len();

    let unique_peptides: std::collections::HashSet<_> = passing_entries
        .iter()
        .map(|e| &library[e.lib_idx].modified_sequence)
        .collect();
    let peptide_count = unique_peptides.len();

    log::info!(
        "  Experiment: {} precursors, {} peptides",
        precursor_count,
        peptide_count
    );
}

/// Count how many of entry A's top-N fragments match entry B's top-N fragments
/// within the given m/z tolerance.
fn count_topn_fragment_overlap(
    frags_a: &[LibraryFragment],
    frags_b: &[LibraryFragment],
    n: usize,
    tolerance: f64,
    unit: ToleranceUnit,
) -> usize {
    // Get top N by intensity for each entry
    let top_a = top_n_fragment_mzs(frags_a, n);
    let mut top_b = top_n_fragment_mzs(frags_b, n);
    top_b.sort_by(|a, b| a.total_cmp(b)); // sort for binary search

    let mut matches = 0;
    for &mz_a in &top_a {
        let tol_da = match unit {
            ToleranceUnit::Ppm => mz_a * tolerance / 1e6,
            ToleranceUnit::Mz => tolerance,
        };
        let lower = mz_a - tol_da;
        let upper = mz_a + tol_da;
        let idx = top_b.partition_point(|&mz| mz < lower);
        if idx < top_b.len() && top_b[idx] <= upper {
            matches += 1;
        }
    }
    matches
}

/// Get the m/z values of the top N fragments by intensity.
fn top_n_fragment_mzs(fragments: &[LibraryFragment], n: usize) -> Vec<f64> {
    if fragments.len() <= n {
        return fragments.iter().map(|f| f.mz).collect();
    }
    let mut indexed: Vec<(f64, f32)> = fragments
        .iter()
        .map(|f| (f.mz, f.relative_intensity))
        .collect();
    indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
    indexed.iter().take(n).map(|(mz, _)| *mz).collect()
}

/// Remove double-counted peptides sharing fragment ions within the same isolation window.
///
/// Two entries are considered double-counted if:
/// 1. Their precursor m/z falls within the same DIA isolation window
/// 2. Their apex RTs are within ±5 spectra of each other
/// 3. ≥50% of their top 6 fragment ions match within calibrated m/z tolerance
///
/// The entry with the lower peak_apex coefficient is removed.
fn deduplicate_double_counting(
    scored_entries: Vec<ScoredEntry>,
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    calibration_params: Option<&CalibrationParams>,
    config: &OspreyConfig,
) -> Vec<ScoredEntry> {
    let original_count = scored_entries.len();

    // 1. Extract unique isolation windows
    let scheme = extract_isolation_scheme(spectra);
    let windows = match &scheme {
        Some(s) => &s.windows,
        None => {
            log::warn!("Could not extract isolation scheme, skipping deduplication");
            return scored_entries;
        }
    };

    // 2. Compute effective fragment tolerance
    let frag_tolerance = if let Some(cal) = calibration_params {
        if cal.ms2_calibration.calibrated {
            let tol_3sd = 3.0 * cal.ms2_calibration.sd;
            let unit = if cal.ms2_calibration.unit == "Th" {
                ToleranceUnit::Mz
            } else {
                ToleranceUnit::Ppm
            };
            let min_tol = if unit == ToleranceUnit::Mz { 0.05 } else { 1.0 };
            FragmentToleranceConfig {
                tolerance: tol_3sd.max(min_tol),
                unit,
            }
        } else {
            config.fragment_tolerance
        }
    } else {
        config.fragment_tolerance
    };

    // 3. Compute median scan interval for "±5 spectra" RT neighborhood
    let rt_neighborhood = {
        let mut intervals: Vec<f64> = Vec::new();
        let mut sorted_rts: Vec<f64> = spectra.iter().map(|s| s.retention_time).collect();
        sorted_rts.sort_by(|a, b| a.total_cmp(b));
        sorted_rts.dedup();
        for w in sorted_rts.windows(2) {
            intervals.push(w[1] - w[0]);
        }
        let median_interval = if intervals.is_empty() {
            0.05 // fallback: 3 seconds
        } else {
            intervals.sort_by(|a, b| a.total_cmp(b));
            intervals[intervals.len() / 2]
        };
        5.0 * median_interval
    };

    // 4. Pre-compute apex RT for each entry
    let apex_rts: Vec<f64> = scored_entries
        .iter()
        .map(|e| {
            e.rt_coef_pairs
                .iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .map(|(rt, _)| *rt)
                .unwrap_or(0.0)
        })
        .collect();

    // 5. For each isolation window, find and deduplicate overlapping entries
    let mut removed = vec![false; scored_entries.len()];

    for &(center, width) in windows {
        let win_lower = center - width / 2.0;
        let win_upper = center + width / 2.0;

        // Collect indices of entries in this window
        let mut window_indices: Vec<usize> = scored_entries
            .iter()
            .enumerate()
            .filter(|(i, e)| {
                !removed[*i] && {
                    let pmz = library[e.lib_idx].precursor_mz;
                    pmz >= win_lower && pmz <= win_upper
                }
            })
            .map(|(i, _)| i)
            .collect();

        // Sort by peak_apex descending (greedy: best entries survive)
        window_indices.sort_by(|&a, &b| {
            scored_entries[b]
                .features
                .peak_apex
                .total_cmp(&scored_entries[a].features.peak_apex)
        });

        // Compare pairs: for each entry, check against subsequent entries
        for i_pos in 0..window_indices.len() {
            let idx_a = window_indices[i_pos];
            if removed[idx_a] {
                continue;
            }
            let apex_a = apex_rts[idx_a];

            for &idx_b in &window_indices[(i_pos + 1)..] {
                if removed[idx_b] {
                    continue;
                }

                // Check apex RT proximity (±5 spectra)
                let apex_b = apex_rts[idx_b];
                if (apex_a - apex_b).abs() > rt_neighborhood {
                    continue;
                }

                // Check top-6 fragment overlap
                let entry_a = &library[scored_entries[idx_a].lib_idx];
                let entry_b = &library[scored_entries[idx_b].lib_idx];
                let overlap = count_topn_fragment_overlap(
                    &entry_a.fragments,
                    &entry_b.fragments,
                    6,
                    frag_tolerance.tolerance,
                    frag_tolerance.unit,
                );

                // 50% of 6 = 3 matching fragments
                let min_a = entry_a.fragments.len().min(6);
                let min_b = entry_b.fragments.len().min(6);
                let threshold = (min_a.min(min_b) as f64 * 0.5).ceil() as usize;

                if overlap >= threshold {
                    // Remove the lower-scoring entry (idx_b, since sorted descending)
                    removed[idx_b] = true;
                }
            }
        }
    }

    let removed_count = removed.iter().filter(|&&r| r).count();
    if removed_count > 0 {
        log::info!(
            "Deduplication: removed {} double-counted entries ({} remaining)",
            removed_count,
            original_count - removed_count
        );
    }

    scored_entries
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !removed[*i])
        .map(|(_, e)| e)
        .collect()
}

/// Write coefficient matrix to a Parquet file.
///
/// Each row is a scored peptide, each data column is a spectrum.
/// Column names encode scan number, precursor m/z, and RT:
///   `scan_{scan_number}_mz_{precursor_mz:.1}_rt_{rt:.2}`
///
/// The file is written alongside the input mzML as `{stem}.coefficients.parquet`.
fn write_coefficient_parquet(
    scored_entries: &[ScoredEntry],
    library: &[LibraryEntry],
    file_results: &[RegressionResult],
    spectra: &[Spectrum],
    input_path: &std::path::Path,
) -> Result<()> {
    use parquet::arrow::ArrowWriter;
    use parquet::basic::Compression;
    use parquet::file::properties::WriterProperties;
    use std::collections::HashMap;

    if scored_entries.is_empty() {
        return Ok(());
    }

    let n_entries = scored_entries.len();

    // Build spectrum info sorted by scan number for column ordering
    let mut spec_info: Vec<(u32, f64, f64)> = spectra
        .iter()
        .map(|s| (s.scan_number, s.precursor_mz, s.retention_time))
        .collect();
    spec_info.sort_by_key(|&(scan, _, _)| scan);
    spec_info.dedup_by_key(|s| s.0);

    // Map lib_idx → row index in scored_entries for fast lookup
    let entry_row: HashMap<usize, usize> = scored_entries
        .iter()
        .enumerate()
        .map(|(row, e)| (e.lib_idx, row))
        .collect();

    // Build schema: metadata columns + one Float32 column per spectrum
    let mut fields: Vec<Field> = vec![
        Field::new("peptide_sequence", DataType::Utf8, false),
        Field::new("modified_sequence", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("precursor_mz", DataType::Float32, false),
        Field::new("protein", DataType::Utf8, false),
        Field::new("is_decoy", DataType::Boolean, false),
    ];

    let col_names: Vec<String> = spec_info
        .iter()
        .map(|&(scan, mz, rt)| format!("scan_{}_mz_{:.1}_rt_{:.2}", scan, mz, rt))
        .collect();

    for name in &col_names {
        fields.push(Field::new(name, DataType::Float32, false));
    }

    let schema = std::sync::Arc::new(Schema::new(fields));

    // Build metadata columns
    let mut peptide_builder = StringBuilder::with_capacity(n_entries, n_entries * 20);
    let mut modified_builder = StringBuilder::with_capacity(n_entries, n_entries * 30);
    let mut charge_builder = UInt8Builder::with_capacity(n_entries);
    let mut pmz_builder = Float32Builder::with_capacity(n_entries);
    let mut protein_builder = StringBuilder::with_capacity(n_entries, n_entries * 15);
    let mut decoy_builder = BooleanBuilder::with_capacity(n_entries);

    for entry in scored_entries {
        let lib = &library[entry.lib_idx];
        peptide_builder.append_value(&lib.sequence);
        modified_builder.append_value(&lib.modified_sequence);
        charge_builder.append_value(lib.charge);
        pmz_builder.append_value(lib.precursor_mz as f32);
        protein_builder.append_value(lib.protein_ids.join(";"));
        decoy_builder.append_value(entry.is_decoy);
    }

    // Map scan_number → index in spec_info for building data columns
    let scan_to_col: HashMap<u32, usize> = spec_info
        .iter()
        .enumerate()
        .map(|(i, &(scan, _, _))| (scan, i))
        .collect();

    // Pre-build a matrix: n_entries × n_spectra, row-major
    // Initialize all zeros
    let n_spectra = spec_info.len();
    let mut coef_matrix = vec![0.0f32; n_entries * n_spectra];

    // Fill from regression results (more efficient than scored_entries.rt_coef_pairs
    // because we get exact scan_number mapping)
    for result in file_results {
        if let Some(&col_idx) = scan_to_col.get(&result.scan_number) {
            for (lib_id, coef) in result.library_ids.iter().zip(result.coefficients.iter()) {
                let lib_idx = (*lib_id & 0x7FFFFFFF) as usize; // mask off decoy high bit
                if let Some(&row_idx) = entry_row.get(&lib_idx) {
                    coef_matrix[row_idx * n_spectra + col_idx] = *coef as f32;
                }
            }
        }
    }

    // Build Arrow arrays
    let mut columns: Vec<ArrayRef> = vec![
        std::sync::Arc::new(peptide_builder.finish()),
        std::sync::Arc::new(modified_builder.finish()),
        std::sync::Arc::new(charge_builder.finish()),
        std::sync::Arc::new(pmz_builder.finish()),
        std::sync::Arc::new(protein_builder.finish()),
        std::sync::Arc::new(decoy_builder.finish()),
    ];

    // Add spectrum columns from the matrix
    for col_idx in 0..n_spectra {
        let mut builder = Float32Builder::with_capacity(n_entries);
        for row_idx in 0..n_entries {
            builder.append_value(coef_matrix[row_idx * n_spectra + col_idx]);
        }
        columns.push(std::sync::Arc::new(builder.finish()));
    }

    // Write parquet file
    let output_path = input_path.with_extension("coefficients.parquet");
    let file = std::fs::File::create(&output_path).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to create parquet file {}: {}",
            output_path.display(),
            e
        ))
    })?;

    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();

    let batch = RecordBatch::try_new(schema.clone(), columns).map_err(|e| {
        OspreyError::OutputError(format!("Failed to create Arrow RecordBatch: {}", e))
    })?;

    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| OspreyError::OutputError(format!("Failed to create Parquet writer: {}", e)))?;

    writer
        .write(&batch)
        .map_err(|e| OspreyError::OutputError(format!("Failed to write Parquet batch: {}", e)))?;

    writer
        .close()
        .map_err(|e| OspreyError::OutputError(format!("Failed to close Parquet writer: {}", e)))?;

    log::info!(
        "Wrote coefficient matrix ({} peptides × {} spectra) to {}",
        n_entries,
        n_spectra,
        output_path.display()
    );

    Ok(())
}

/// Create PSM features for a single file
///
/// This is used for memory-efficient processing where each file's PSM features
/// are created immediately after scoring and written to a PIN file.
fn create_psm_features_for_file(
    file_name: &str,
    entries: &[ScoredEntry],
    library: &[LibraryEntry],
) -> Vec<PsmFeatures> {
    entries
        .iter()
        .map(|entry| {
            let lib_entry = &library[entry.lib_idx];
            PsmFeatures {
                psm_id: entry.psm_id.clone(),
                peptide: lib_entry.modified_sequence.clone(),
                proteins: lib_entry.protein_ids.clone(),
                scan_number: entry.apex_scan_number,
                file_name: file_name.to_string(),
                charge: lib_entry.charge,
                is_decoy: entry.is_decoy,
                features: entry.features.clone(),
                initial_score: Some(entry.score),
            }
        })
        .collect()
}

/// Collect PSM features by file for multi-file mokapot (legacy - used when PIN files not pre-written)
#[allow(dead_code)]
fn collect_psm_features_by_file(
    per_file_results: &[(String, Vec<ScoredEntry>)],
    library: &[LibraryEntry],
) -> HashMap<String, Vec<PsmFeatures>> {
    let mut psms_by_file = HashMap::new();

    for (file_name, entries) in per_file_results.iter() {
        let psm_features = create_psm_features_for_file(file_name, entries, library);
        psms_by_file.insert(file_name.clone(), psm_features);
    }

    psms_by_file
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
        if !entry.is_decoy && target_idx < qvalues.len() {
            entry.run_qvalue = qvalues[target_idx];
            target_idx += 1;
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
        if !entry.is_decoy && target_idx < qvalues.len() {
            entry.experiment_qvalue = qvalues[target_idx];
            target_idx += 1;
        }
    }

    Ok(())
}

/// Compute peak boundaries and apex for a scored entry, ensuring apex is within bounds
fn compute_peak_boundaries(scored: &ScoredEntry) -> Option<(f64, f64, f64)> {
    // Find apex from coefficient time series
    let (coef_apex_rt, apex_coef) = scored
        .rt_coef_pairs
        .iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(rt, c)| (*rt, *c))
        .unwrap_or((0.0, 0.0));

    if apex_coef <= 0.0 {
        return None;
    }

    let (start_rt, end_rt) = if let Some((start, end, _fwhm)) = scored.fragment_peak_bounds {
        // Fragment-based boundaries from co-eluting XICs (preferred)
        (start, end)
    } else if let Some((fwhm, _, _)) =
        osprey_scoring::compute_fwhm_interpolated(&scored.rt_coef_pairs)
    {
        // Fallback: coefficient series FWHM with 95% Gaussian boundaries
        let sigma = fwhm / 2.355;
        (coef_apex_rt - 1.96 * sigma, coef_apex_rt + 1.96 * sigma)
    } else {
        // Final fallback: use first/last non-zero coefficient as boundaries
        let first_rt = scored
            .rt_coef_pairs
            .iter()
            .find(|(_, c)| *c > 0.0)
            .map(|(rt, _)| *rt)
            .unwrap_or(coef_apex_rt);
        let last_rt = scored
            .rt_coef_pairs
            .iter()
            .rev()
            .find(|(_, c)| *c > 0.0)
            .map(|(rt, _)| *rt)
            .unwrap_or(coef_apex_rt);
        (first_rt, last_rt)
    };

    // Ensure apex RT is within [startTime, endTime]
    // If coefficient apex falls outside fragment bounds, use the highest-coefficient
    // point within the bounds, or clamp to the midpoint
    let apex_rt = if coef_apex_rt >= start_rt && coef_apex_rt <= end_rt {
        coef_apex_rt
    } else {
        // Find the highest coefficient within the bounds
        scored
            .rt_coef_pairs
            .iter()
            .filter(|(rt, _)| *rt >= start_rt && *rt <= end_rt)
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, _)| *rt)
            .unwrap_or((start_rt + end_rt) / 2.0)
    };

    Some((apex_rt, start_rt, end_rt))
}

/// Write scored results to blib format for Skyline
///
/// Groups entries by precursor (lib_idx) to produce one RefSpectra row per unique
/// precursor, with per-run data in the RetentionTimes table. The best-scoring run's
/// data is used for the RefSpectra entry. Entries from each run where the precursor
/// passed run-level FDR are included in RetentionTimes.
fn write_blib_output_with_scores(
    config: &OspreyConfig,
    library: &[LibraryEntry],
    scored_entries: &[ScoredEntry],
    input_files: &[std::path::PathBuf],
) -> Result<()> {
    let mut writer = BlibWriter::create(&config.output_blib)?;

    // Add metadata
    writer.add_metadata("osprey_version", env!("CARGO_PKG_VERSION"))?;
    writer.add_metadata(
        "rt_calibration_enabled",
        &config.rt_calibration.enabled.to_string(),
    )?;
    writer.add_metadata("run_fdr", &config.run_fdr.to_string())?;
    writer.add_metadata("fdr_method", "target_decoy_competition")?;

    // Determine library name for idFileName in SpectrumSourceFiles
    let library_name = config
        .library_source
        .path()
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("library")
        .to_string();

    // Add source files - build lookup by file stem for matching with ScoredEntry.file_name
    // Use relative paths from blib location for cross-platform compatibility (WSL2 → Windows)
    let blib_dir = config.output_blib.parent();
    let mut file_stem_to_id: HashMap<String, i64> = HashMap::new();
    for input_file in input_files {
        let file_name = input_file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        let file_stem = input_file
            .file_stem()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        // Compute relative path from blib directory to mzML file
        // Falls back to just filename if relative path can't be computed
        let source_path = if let Some(blib_parent) = blib_dir {
            pathdiff::diff_paths(input_file, blib_parent)
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_else(|| file_name.to_string())
        } else {
            file_name.to_string()
        };
        let file_id = writer.add_source_file(&source_path, &library_name, config.run_fdr)?;
        file_stem_to_id.insert(file_stem.to_string(), file_id);
    }

    // Begin batch transaction for much faster writes
    writer.begin_batch()?;

    // Group passing entries by precursor (lib_idx) so we write one RefSpectra per
    // unique precursor, with per-run entries in RetentionTimes
    let mut precursor_entries: HashMap<usize, Vec<&ScoredEntry>> = HashMap::new();
    for scored in scored_entries {
        precursor_entries
            .entry(scored.lib_idx)
            .or_default()
            .push(scored);
    }

    let mut n_written = 0usize;

    for (lib_idx, entries) in &precursor_entries {
        let entry = &library[*lib_idx];

        // Find the best-scoring run for this precursor (lowest run q-value)
        let best = entries
            .iter()
            .min_by(|a, b| a.run_qvalue.total_cmp(&b.run_qvalue))
            .unwrap();

        // Compute boundaries for the best run
        let (apex_rt, start_rt, end_rt) = match compute_peak_boundaries(best) {
            Some(bounds) => bounds,
            None => continue,
        };

        // Count runs where this precursor passed run-level FDR
        let runs_passing: Vec<&&ScoredEntry> = entries
            .iter()
            .filter(|e| e.run_qvalue <= config.run_fdr)
            .collect();
        let n_runs_detected = runs_passing.len() as i32;

        // Look up file_id for the best run
        let file_id = *file_stem_to_id.get(&best.file_name).unwrap_or(&1);

        // Get fragment peaks from library entry
        let mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
        let intensities: Vec<f32> = entry
            .fragments
            .iter()
            .map(|f| f.relative_intensity)
            .collect();

        // Compute total ion current as sum of log intensities (matching DIA-NN convention)
        let total_ion_current: f64 = intensities
            .iter()
            .filter(|&&i| i > 0.0)
            .map(|&i| (i as f64).ln())
            .sum();

        // Score is the raw q-value (lower is better) — Skyline GENERIC Q-VALUE convention
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
            best.experiment_qvalue,
            file_id,
            n_runs_detected,
            total_ion_current,
        )?;

        // Add modifications if present
        if !entry.modifications.is_empty() {
            writer.add_modifications(ref_id, &entry.modifications)?;
        }

        // Add protein mappings if present
        if !entry.protein_ids.is_empty() {
            writer.add_protein_mapping(ref_id, &entry.protein_ids)?;
        }

        // Write RetentionTimes entries for each run where this precursor passed run FDR
        for scored in &runs_passing {
            let run_file_id = *file_stem_to_id.get(&scored.file_name).unwrap_or(&1);
            let is_best: bool = {
                let a: *const ScoredEntry = **scored;
                let b: *const ScoredEntry = *best;
                std::ptr::eq(a, b)
            };

            if let Some((run_apex, run_start, run_end)) = compute_peak_boundaries(scored) {
                writer.add_retention_time(
                    ref_id,
                    run_file_id,
                    run_apex,
                    run_start,
                    run_end,
                    scored.run_qvalue,
                    is_best,
                )?;
            }
        }

        // Add peak boundaries to Osprey extension table
        let boundaries = osprey_core::PeakBoundaries {
            start_rt,
            end_rt,
            apex_rt,
            apex_coefficient: best.score,
            integrated_area: best.features.peak_area,
            peak_quality: osprey_core::PeakQuality::default(),
        };
        writer.add_peak_boundaries(ref_id, &best.file_name, &boundaries)?;

        // Add run scores for best run
        writer.add_run_scores(
            ref_id,
            &best.file_name,
            best.run_qvalue,
            best.score,
            best.pep,
        )?;

        // Add experiment-level scores
        writer.add_experiment_scores(
            ref_id,
            best.experiment_qvalue,
            n_runs_detected,
            input_files.len() as i32,
        )?;

        n_written += 1;
    }

    // Commit the batch transaction
    writer.commit()?;
    writer.finalize()?;

    log::info!(
        "Wrote {} unique precursors to blib (from {} total entries across {} files)",
        n_written,
        scored_entries.len(),
        input_files.len()
    );

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
        let apex_rt = scored
            .rt_coef_pairs
            .iter()
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
        .map(|scored| {
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
                psm_id: scored.psm_id.clone(),
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
    writer.add_metadata(
        "rt_calibration_enabled",
        &config.rt_calibration.enabled.to_string(),
    )?;
    writer.add_metadata(
        "rt_tolerance_factor",
        &config.rt_calibration.rt_tolerance_factor.to_string(),
    )?;
    writer.add_metadata(
        "fallback_rt_tolerance",
        &config.rt_calibration.fallback_rt_tolerance.to_string(),
    )?;
    writer.add_metadata("run_fdr", &config.run_fdr.to_string())?;

    // Add source files - use relative path for cross-platform compatibility (WSL2 → Windows)
    let blib_dir = config.output_blib.parent();
    let mut file_ids = std::collections::HashMap::new();
    for input_file in input_files {
        let file_name = input_file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        // Compute relative path from blib directory to mzML file
        let source_path = if let Some(blib_parent) = blib_dir {
            pathdiff::diff_paths(input_file, blib_parent)
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_else(|| file_name.to_string())
        } else {
            file_name.to_string()
        };
        let file_id = writer.add_source_file(&source_path, &source_path, 0.0)?;
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
        let above_threshold: Vec<_> = sorted_pairs
            .iter()
            .filter(|(_, c)| *c >= threshold)
            .collect();

        if above_threshold.is_empty() {
            continue;
        }

        let start_rt = above_threshold.first().map(|(rt, _)| *rt).unwrap_or(0.0);
        let end_rt = above_threshold.last().map(|(rt, _)| *rt).unwrap_or(0.0);
        let (apex_rt, apex_coef) = sorted_pairs
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, c)| (*rt, *c))
            .unwrap_or((0.0, 0.0));

        // Use the first input file's ID (in a real implementation, track per-file)
        let file_id = *file_ids.values().next().unwrap_or(&1);

        // Get fragment peaks from library entry
        let mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
        let intensities: Vec<f32> = entry
            .fragments
            .iter()
            .map(|f| f.relative_intensity)
            .collect();

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
            1,   // copies
            0.0, // total_ion_current
        )?;

        // Add modifications if present
        if !entry.modifications.is_empty() {
            writer.add_modifications(ref_id, &entry.modifications)?;
        }

        // Add protein mappings if present
        if !entry.protein_ids.is_empty() {
            writer.add_protein_mapping(ref_id, &entry.protein_ids)?;
        }

        // Add peak boundaries
        let boundaries = osprey_core::PeakBoundaries {
            start_rt,
            end_rt,
            apex_rt,
            apex_coefficient: apex_coef,
            integrated_area: sorted_pairs.iter().map(|(_, c)| c).sum(),
            peak_quality: osprey_core::PeakQuality::default(),
        };

        let file_name = input_files
            .first()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        writer.add_peak_boundaries(ref_id, file_name, &boundaries)?;

        // Add run scores (placeholder - will be replaced with proper FDR)
        writer.add_run_scores(ref_id, file_name, score, apex_coef, 1.0 - score)?;
    }

    // Commit the batch transaction
    writer.commit()?;
    writer.finalize()?;

    log::info!("Wrote {} spectra to blib", entry_coefficients.len());

    Ok(())
}

/// RT-sorted index for fast candidate selection via binary search.
///
/// Each m/z bin contains (expected_rt, library_index) pairs sorted by expected_rt.
/// When calibration is provided, expected_rt = cal.predict(library_rt); otherwise
/// expected_rt = library_rt. Sorting enables binary search for the RT window,
/// reducing candidate selection from O(n) to O(log n + k) per m/z bin.
struct MzRTIndex {
    bins: HashMap<i32, Vec<(f64, usize)>>,
}

impl MzRTIndex {
    /// Build index from library, optionally using RT calibration to pre-compute expected RTs.
    fn build(library: &[LibraryEntry], calibration: Option<&RTCalibration>) -> Self {
        let mut bins: HashMap<i32, Vec<(f64, usize)>> = HashMap::new();

        for (idx, entry) in library.iter().enumerate() {
            let expected_rt = if let Some(cal) = calibration {
                cal.predict(entry.retention_time)
            } else {
                entry.retention_time
            };

            let mz_bin = entry.precursor_mz.round() as i32;
            for offset in -1..=1 {
                bins.entry(mz_bin + offset)
                    .or_default()
                    .push((expected_rt, idx));
            }
        }

        // Sort each bin by expected RT for binary search
        for entries in bins.values_mut() {
            entries.sort_by(|a, b| a.0.total_cmp(&b.0));
        }

        MzRTIndex { bins }
    }

    fn get(&self, mz_bin: i32) -> Option<&[(f64, usize)]> {
        self.bins.get(&mz_bin).map(|v| v.as_slice())
    }

    #[cfg(test)]
    fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }
}

/// Build an index of library entries by precursor m/z (unsorted, for tests)
#[cfg(test)]
fn build_mz_index(library: &[LibraryEntry]) -> MzRTIndex {
    MzRTIndex::build(library, None)
}

/// Select candidate library entries for a spectrum using RT-sorted binary search.
///
/// Candidate selection uses:
/// 1. Isolation window from mzML (defines which precursors are fragmented)
/// 2. RT tolerance — binary search on pre-sorted expected_rt within each m/z bin
/// 3. Optional library filter (for calibration discovery phase)
/// 4. Optional RT calibration (local tolerance refinement after binary search)
/// 5. Optional top-3 fragment pre-filter (requires at least 1 of top 3 peaks in spectrum)
///
/// The MzRTIndex stores entries sorted by expected_rt (calibrated if available),
/// enabling O(log n + k) candidate selection instead of O(n) per m/z bin.
#[allow(clippy::too_many_arguments)]
fn select_candidates_with_calibration(
    spectrum: &Spectrum,
    library: &[LibraryEntry],
    mz_rt_index: &MzRTIndex,
    rt_tolerance: f64,
    max_candidates: usize,
    min_rt_tolerance: f64,
    library_filter: Option<&std::collections::HashSet<usize>>,
    calibration: Option<&RTCalibration>,
    fragment_tolerance: Option<FragmentToleranceConfig>,
    search_tolerance: f64,
) -> Vec<usize> {
    let mut candidates = Vec::new();
    let mut seen = std::collections::HashSet::new();

    let spectrum_rt = spectrum.retention_time;

    // Get the range of bins to search based on isolation window
    let lower_bin = spectrum.isolation_window.lower_bound().floor() as i32;
    let upper_bin = spectrum.isolation_window.upper_bound().ceil() as i32;

    // Binary search bounds for RT
    let rt_low = spectrum_rt - search_tolerance;
    let rt_high = spectrum_rt + search_tolerance;

    // Search all bins that could contain candidates within the isolation window
    for bin in lower_bin..=upper_bin {
        if let Some(entries) = mz_rt_index.get(bin) {
            // Binary search to find the start of the RT range
            let start = entries.partition_point(|&(rt, _)| rt < rt_low);
            // Binary search to find the end of the RT range
            let end = entries.partition_point(|&(rt, _)| rt <= rt_high);

            for &(expected_rt, idx) in &entries[start..end] {
                // Skip if already processed (entry can appear in adjacent m/z bins)
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
                if !spectrum.isolation_window.contains(entry.precursor_mz) {
                    continue;
                }

                // Exact local tolerance check (binary search used conservative global bound)
                if let Some(cal) = calibration {
                    let factor = rt_tolerance / cal.residual_std().max(0.001);
                    let effective_tolerance =
                        cal.local_tolerance(entry.retention_time, factor, min_rt_tolerance);
                    if (expected_rt - spectrum_rt).abs() > effective_tolerance {
                        continue;
                    }
                }

                // Apply top-3 fragment pre-filter if enabled
                if let Some(ref frag_tol) = fragment_tolerance {
                    if !has_top3_fragment_match(
                        &entry.fragments,
                        &spectrum.mzs,
                        frag_tol.tolerance,
                        frag_tol.unit,
                    ) {
                        continue;
                    }
                }

                candidates.push(idx);
            }
        }
    }

    // Limit to max candidates (keep those closest in RT)
    if candidates.len() > max_candidates {
        candidates.sort_by(|&a, &b| {
            let rt_a = if let Some(cal) = calibration {
                (cal.predict(library[a].retention_time) - spectrum_rt).abs()
            } else {
                (library[a].retention_time - spectrum_rt).abs()
            };
            let rt_b = if let Some(cal) = calibration {
                (cal.predict(library[b].retention_time) - spectrum_rt).abs()
            } else {
                (library[b].retention_time - spectrum_rt).abs()
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
    library_by_mz: &MzRTIndex,
    rt_tolerance: f64,
    max_candidates: usize,
) -> Vec<usize> {
    select_candidates_with_calibration(
        spectrum,
        library,
        library_by_mz,
        rt_tolerance,
        max_candidates,
        0.1, // Default min_rt_tolerance for tests
        None,
        None,
        None,         // No fragment pre-filter for tests
        rt_tolerance, // search_tolerance = rt_tolerance for uncalibrated
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

    let file = File::create(output_path)
        .map_err(|e| OspreyError::OutputError(format!("Failed to create debug CSV: {}", e)))?;
    let mut writer = std::io::BufWriter::new(file);

    // Pair targets with their decoys
    let paired = pair_calibration_matches(matches, expected_rt_fn);

    log::debug!(
        "Writing {} paired target-decoy results to debug CSV",
        paired.len()
    );

    // Write header with paired columns
    // Note: sorted by winning_evalue ascending (best matches first, lower E-value = better)
    writeln!(
        writer,
        "target_entry_id,charge,target_sequence,decoy_sequence,\
         winning_evalue,target_evalue,decoy_evalue,target_wins_evalue,\
         winning_xcorr,target_xcorr,decoy_xcorr,\
         winning_hyperscore,target_hyperscore,decoy_hyperscore,\
         target_n_b,target_n_y,decoy_n_b,decoy_n_y,\
         target_isotope_score,decoy_isotope_score,\
         target_precursor_error_ppm,decoy_precursor_error_ppm,\
         target_rt,decoy_rt,library_rt,expected_rt,target_delta_rt,decoy_delta_rt,\
         target_matched_frags,decoy_matched_frags,target_wins"
    )
    .map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

    // Write each paired result (sorted by winning_evalue ascending - best first)
    for p in &paired {
        writeln!(
            writer,
            "{},{},{},{},{:.2e},{:.2e},{:.2e},{},{:.6},{:.6},{:.6},{:.4},{:.4},{:.4},{},{},{},{},{},{},{},{},{:.2},{:.2},{:.2},{},{:.2},{:.2},{},{},{}",
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
            p.winning_hyperscore,
            p.target_hyperscore,
            p.decoy_hyperscore,
            p.target_n_b,
            p.target_n_y,
            p.decoy_n_b,
            p.decoy_n_y,
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

    let file = File::create(output_path)
        .map_err(|e| OspreyError::OutputError(format!("Failed to create debug CSV: {}", e)))?;
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
                m.isotope_cosine_score
                    .map(|v| format!("{:.4}", v))
                    .unwrap_or_default(),
                m.n_matched_fragments,
                m.is_decoy
            )
            .map_err(|e| OspreyError::OutputError(format!("Failed to write row: {}", e)))?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{FragmentAnnotation, IsolationWindow};

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

    /// Creates a library entry with explicit fragment ions at given m/z values.
    fn make_entry_with_fragments(id: u32, mz: f64, rt: f64, frag_mzs: &[f64]) -> LibraryEntry {
        let mut entry = LibraryEntry::new(
            id,
            format!("PEPTIDE{}", id),
            format!("PEPTIDE{}", id),
            2,
            mz,
            rt,
        );
        entry.fragments = frag_mzs
            .iter()
            .enumerate()
            .map(|(i, &fmz)| LibraryFragment {
                mz: fmz,
                relative_intensity: 100.0 - (i as f32 * 10.0), // decreasing intensity
                annotation: FragmentAnnotation::default(),
            })
            .collect();
        entry
    }

    /// Creates a set of DIA-like spectra with the given isolation window center and width.
    fn make_dia_spectra(center: f64, width: f64, n_scans: usize) -> Vec<Spectrum> {
        (0..n_scans)
            .map(|i| Spectrum {
                scan_number: i as u32 + 1,
                retention_time: 10.0 + (i as f64 * 0.05), // 3-second cycle time
                precursor_mz: center,
                isolation_window: IsolationWindow::symmetric(center, width / 2.0),
                mzs: vec![300.0, 400.0, 500.0],
                intensities: vec![100.0, 200.0, 300.0],
            })
            .collect()
    }

    /// Verifies candidate selection includes entries within the isolation window and RT tolerance.
    #[test]
    fn test_select_candidates() {
        let spectrum = make_test_spectrum();
        let library = vec![
            make_test_entry(0, 500.0, 10.0), // Should match - in window and correct RT
            make_test_entry(1, 505.0, 10.0), // In window and correct RT - should match
            make_test_entry(2, 520.0, 10.0), // Outside isolation window
            make_test_entry(3, 500.0, 20.0), // In window but wrong RT
        ];
        let mz_index = build_mz_index(&library);

        // With rt_tolerance=2.0, entries 0 and 1 should match (in isolation window and RT)
        let candidates = select_candidates(&spectrum, &library, &mz_index, 2.0, 100);

        assert!(
            candidates.contains(&0),
            "Entry 0 should match (in window, correct RT)"
        );
        assert!(
            candidates.contains(&1),
            "Entry 1 should match (in window, correct RT)"
        );
        assert!(
            !candidates.contains(&2),
            "Entry 2 should not match (outside window)"
        );
        assert!(
            !candidates.contains(&3),
            "Entry 3 should not match (wrong RT)"
        );
    }

    /// Verifies MzRTIndex bins each entry into 3 adjacent m/z bins (±1 Da).
    #[test]
    fn test_build_mz_index() {
        let library = vec![
            make_test_entry(0, 500.0, 10.0),
            make_test_entry(1, 500.3, 10.0), // same integer bin as entry 0
            make_test_entry(2, 600.0, 10.0), // different bin
        ];

        let index = build_mz_index(&library);

        // Helper: check if a bin contains a given library index
        let bin_contains = |bin: i32, idx: usize| -> bool {
            index
                .get(bin)
                .map_or(false, |entries| entries.iter().any(|&(_, i)| i == idx))
        };

        // Entry at 500.0 should appear in bins 499, 500, 501
        assert!(bin_contains(500, 0));
        assert!(bin_contains(499, 0));
        assert!(bin_contains(501, 0));

        // Entry at 500.3 also rounds to 500, so bins 499, 500, 501
        assert!(bin_contains(500, 1));

        // Entry at 600.0 should appear in bins 599, 600, 601
        assert!(bin_contains(600, 2));
        assert!(!bin_contains(600, 0));
    }

    /// Verifies MzRTIndex returns empty for an empty library.
    #[test]
    fn test_build_mz_index_empty() {
        let index = build_mz_index(&[]);
        assert!(index.is_empty());
    }

    /// Verifies that candidate selection respects the max_candidates limit.
    #[test]
    fn test_select_candidates_max_limit() {
        let spectrum = make_test_spectrum();
        // Create many entries that all match
        let library: Vec<LibraryEntry> = (0..20)
            .map(|i| make_test_entry(i, 500.0 + (i as f64 * 0.1), 10.0))
            .collect();
        let mz_index = build_mz_index(&library);

        let candidates = select_candidates(&spectrum, &library, &mz_index, 2.0, 5);
        assert!(candidates.len() <= 5, "Should respect max_candidates limit");
    }

    /// Verifies count_topn_fragment_overlap detects identical top-6 fragments.
    #[test]
    fn test_fragment_overlap_identical() {
        let frags: Vec<LibraryFragment> = [300.0, 400.0, 500.0, 600.0, 700.0, 800.0]
            .iter()
            .map(|&mz| LibraryFragment {
                mz,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            })
            .collect();

        // Same fragments should have 6/6 overlap
        let overlap = count_topn_fragment_overlap(&frags, &frags, 6, 0.5, ToleranceUnit::Mz);
        assert_eq!(overlap, 6, "Identical fragments should fully overlap");
    }

    /// Verifies count_topn_fragment_overlap returns 0 for completely disjoint fragments.
    #[test]
    fn test_fragment_overlap_disjoint() {
        let frags_a: Vec<LibraryFragment> = [300.0, 400.0, 500.0]
            .iter()
            .map(|&mz| LibraryFragment {
                mz,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            })
            .collect();
        let frags_b: Vec<LibraryFragment> = [350.0, 450.0, 550.0]
            .iter()
            .map(|&mz| LibraryFragment {
                mz,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            })
            .collect();

        let overlap = count_topn_fragment_overlap(&frags_a, &frags_b, 6, 0.5, ToleranceUnit::Mz);
        assert_eq!(overlap, 0, "Disjoint fragments should have zero overlap");
    }

    /// Verifies count_topn_fragment_overlap uses ppm tolerance correctly.
    #[test]
    fn test_fragment_overlap_ppm_tolerance() {
        let frags_a: Vec<LibraryFragment> = [500.0]
            .iter()
            .map(|&mz| LibraryFragment {
                mz,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            })
            .collect();
        // 500.005 is 10 ppm away from 500.0
        let frags_b: Vec<LibraryFragment> = [500.005]
            .iter()
            .map(|&mz| LibraryFragment {
                mz,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            })
            .collect();

        // Should match with 20 ppm tolerance
        let overlap_wide =
            count_topn_fragment_overlap(&frags_a, &frags_b, 6, 20.0, ToleranceUnit::Ppm);
        assert_eq!(overlap_wide, 1, "Should match within 20 ppm");

        // Should NOT match with 5 ppm tolerance
        let overlap_narrow =
            count_topn_fragment_overlap(&frags_a, &frags_b, 6, 5.0, ToleranceUnit::Ppm);
        assert_eq!(overlap_narrow, 0, "Should not match within 5 ppm");
    }

    /// Verifies top_n_fragment_mzs returns all fragments when n >= fragment count.
    #[test]
    fn test_top_n_fragments() {
        let frags: Vec<LibraryFragment> = [300.0, 400.0, 500.0]
            .iter()
            .enumerate()
            .map(|(i, &mz)| LibraryFragment {
                mz,
                relative_intensity: (100 - i as i32 * 30) as f32,
                annotation: FragmentAnnotation::default(),
            })
            .collect();

        let top3 = top_n_fragment_mzs(&frags, 3);
        assert_eq!(top3.len(), 3);

        // With n=2, should return the two most intense
        let top2 = top_n_fragment_mzs(&frags, 2);
        assert_eq!(top2.len(), 2);
        assert!(
            top2.contains(&300.0),
            "Most intense fragment should be included"
        );
        assert!(
            top2.contains(&400.0),
            "Second most intense should be included"
        );
    }

    /// Verifies extract_isolation_scheme detects DIA window cycle from spectra.
    #[test]
    fn test_extract_isolation_scheme() {
        // Create 2 cycles of 3 DIA windows
        let mut spectra = Vec::new();
        for cycle in 0..2 {
            for (i, center) in [500.0, 525.0, 550.0].iter().enumerate() {
                spectra.push(Spectrum {
                    scan_number: (cycle * 3 + i) as u32 + 1,
                    retention_time: (cycle * 3 + i) as f64 * 0.05,
                    precursor_mz: *center,
                    isolation_window: IsolationWindow::symmetric(*center, 12.5),
                    mzs: vec![300.0],
                    intensities: vec![100.0],
                });
            }
        }

        let scheme = extract_isolation_scheme(&spectra);
        assert!(scheme.is_some(), "Should detect isolation scheme");
        let scheme = scheme.unwrap();
        assert_eq!(scheme.windows.len(), 3, "Should find 3 unique windows");
    }

    /// Verifies extract_isolation_scheme returns None for empty spectra.
    #[test]
    fn test_extract_isolation_scheme_empty() {
        assert!(extract_isolation_scheme(&[]).is_none());
    }

    /// Verifies deduplication removes the lower-scoring entry when two entries
    /// share ≥50% of top-6 fragments, are in the same isolation window,
    /// and have similar apex RTs.
    #[test]
    fn test_deduplicate_removes_overlapping_entry() {
        // Two entries with identical fragments in the same isolation window
        let entry_a =
            make_entry_with_fragments(0, 500.0, 10.0, &[300.0, 400.0, 500.0, 600.0, 700.0, 800.0]);
        let entry_b =
            make_entry_with_fragments(1, 501.0, 10.0, &[300.0, 400.0, 500.0, 600.0, 700.0, 800.0]);
        let library = vec![entry_a, entry_b];

        // Create DIA spectra with a window covering both entries
        let spectra = make_dia_spectra(500.0, 25.0, 20);

        // Scored entries: A has higher peak_apex (survives), B has lower (removed)
        let scored = vec![
            ScoredEntry {
                lib_idx: 0,
                psm_id: "file_0".into(),
                file_name: "file".into(),
                apex_scan_number: 10,
                rt_coef_pairs: vec![(10.5, 100.0)],
                features: FeatureSet {
                    peak_apex: 100.0,
                    ..FeatureSet::default()
                },
                score: 0.9,
                run_qvalue: 0.01,
                experiment_qvalue: 0.01,
                pep: 0.01,
                is_decoy: false,
                fragment_peak_bounds: None,
            },
            ScoredEntry {
                lib_idx: 1,
                psm_id: "file_1".into(),
                file_name: "file".into(),
                apex_scan_number: 10,
                rt_coef_pairs: vec![(10.5, 50.0)],
                features: FeatureSet {
                    peak_apex: 50.0,
                    ..FeatureSet::default()
                },
                score: 0.7,
                run_qvalue: 0.01,
                experiment_qvalue: 0.01,
                pep: 0.01,
                is_decoy: false,
                fragment_peak_bounds: None,
            },
        ];

        let config = OspreyConfig::default();
        let result = deduplicate_double_counting(scored, &library, &spectra, None, &config);

        assert_eq!(
            result.len(),
            1,
            "One entry should be removed (shared fragments)"
        );
        assert_eq!(
            result[0].lib_idx, 0,
            "Higher-scoring entry A should survive"
        );
    }

    /// Verifies deduplication keeps both entries when fragments are disjoint.
    #[test]
    fn test_deduplicate_keeps_disjoint_entries() {
        // Two entries with completely different fragments
        let entry_a =
            make_entry_with_fragments(0, 500.0, 10.0, &[300.0, 400.0, 500.0, 600.0, 700.0, 800.0]);
        let entry_b =
            make_entry_with_fragments(1, 501.0, 10.0, &[310.0, 410.0, 510.0, 610.0, 710.0, 810.0]);
        let library = vec![entry_a, entry_b];

        let spectra = make_dia_spectra(500.0, 25.0, 20);

        let scored = vec![
            ScoredEntry {
                lib_idx: 0,
                psm_id: "file_0".into(),
                file_name: "file".into(),
                apex_scan_number: 10,
                rt_coef_pairs: vec![(10.5, 100.0)],
                features: FeatureSet {
                    peak_apex: 100.0,
                    ..FeatureSet::default()
                },
                score: 0.9,
                run_qvalue: 0.01,
                experiment_qvalue: 0.01,
                pep: 0.01,
                is_decoy: false,
                fragment_peak_bounds: None,
            },
            ScoredEntry {
                lib_idx: 1,
                psm_id: "file_1".into(),
                file_name: "file".into(),
                apex_scan_number: 10,
                rt_coef_pairs: vec![(10.5, 50.0)],
                features: FeatureSet {
                    peak_apex: 50.0,
                    ..FeatureSet::default()
                },
                score: 0.7,
                run_qvalue: 0.01,
                experiment_qvalue: 0.01,
                pep: 0.01,
                is_decoy: false,
                fragment_peak_bounds: None,
            },
        ];

        let config = OspreyConfig::default();
        let result = deduplicate_double_counting(scored, &library, &spectra, None, &config);

        assert_eq!(
            result.len(),
            2,
            "Both entries should survive (different fragments)"
        );
    }

    /// Verifies deduplication keeps entries that are distant in RT even with shared fragments.
    #[test]
    fn test_deduplicate_keeps_rt_distant_entries() {
        // Two entries with identical fragments but very different apex RTs
        let entry_a =
            make_entry_with_fragments(0, 500.0, 5.0, &[300.0, 400.0, 500.0, 600.0, 700.0, 800.0]);
        let entry_b =
            make_entry_with_fragments(1, 501.0, 25.0, &[300.0, 400.0, 500.0, 600.0, 700.0, 800.0]);
        let library = vec![entry_a, entry_b];

        let spectra = make_dia_spectra(500.0, 25.0, 20);

        // Entry A has apex near RT=10.5, entry B has apex near RT=11.5
        // With 20 spectra at 0.05 min intervals, 5*median_interval = 5*0.05 = 0.25 min
        // So entries >0.25 min apart should NOT be deduplicated
        let scored = vec![
            ScoredEntry {
                lib_idx: 0,
                psm_id: "file_0".into(),
                file_name: "file".into(),
                apex_scan_number: 1,
                rt_coef_pairs: vec![(10.0, 100.0)], // apex at RT=10.0
                features: FeatureSet {
                    peak_apex: 100.0,
                    ..FeatureSet::default()
                },
                score: 0.9,
                run_qvalue: 0.01,
                experiment_qvalue: 0.01,
                pep: 0.01,
                is_decoy: false,
                fragment_peak_bounds: None,
            },
            ScoredEntry {
                lib_idx: 1,
                psm_id: "file_1".into(),
                file_name: "file".into(),
                apex_scan_number: 20,
                rt_coef_pairs: vec![(11.0, 50.0)], // apex at RT=11.0 (1 min away)
                features: FeatureSet {
                    peak_apex: 50.0,
                    ..FeatureSet::default()
                },
                score: 0.7,
                run_qvalue: 0.01,
                experiment_qvalue: 0.01,
                pep: 0.01,
                is_decoy: false,
                fragment_peak_bounds: None,
            },
        ];

        let config = OspreyConfig::default();
        let result = deduplicate_double_counting(scored, &library, &spectra, None, &config);

        assert_eq!(
            result.len(),
            2,
            "Entries with distant apex RTs should both survive"
        );
    }

    /// Verifies deduplication with empty input returns empty output.
    #[test]
    fn test_deduplicate_empty() {
        let spectra = make_dia_spectra(500.0, 25.0, 10);
        let config = OspreyConfig::default();
        let result = deduplicate_double_counting(Vec::new(), &[], &spectra, None, &config);
        assert!(result.is_empty());
    }
}
