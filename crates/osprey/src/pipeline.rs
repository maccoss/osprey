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
    MzQCData, calculate_mz_calibration, save_calibration, calibration_filename,
};
use osprey_core::{
    BinConfig, DecoyMethod as CoreDecoyMethod, FeatureSet, LibraryEntry,
    MS1Spectrum, OspreyConfig, OspreyError, RegressionResult, ResolutionMode, Result, Spectrum,
};
use osprey_fdr::{FdrController, MokapotRunner, PsmFeatures};
use osprey_io::{load_library, load_ms1_spectra, BlibWriter, MS1Index, MzmlReader};
use osprey_regression::{Binner, DesignMatrixBuilder, RidgeSolver};
use osprey_scoring::{
    batch::{
        BatchScorer, MS1SpectrumLookup, PreprocessedLibrary, PreprocessedSpectra,
        run_libcosine_calibration_scoring_with_ms1, sample_library_for_calibration,
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
    let matrix_builder = DesignMatrixBuilder::new(binner.clone());

    // Get regularization lambda
    let lambda = match &config.regularization_lambda {
        osprey_core::RegularizationSetting::Fixed(l) => *l,
        _ => 1.0, // Default lambda for now
    };

    let solver = RidgeSolver::new(lambda);

    // Build lookup structure for candidates
    let library_by_mz = build_mz_index(&library);

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
    let mut shared_calibration: Option<(RTCalibration, f64)> = None; // (calibration, tolerance)

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
                &matrix_builder,
                &solver,
                &config,
                shared_calibration.as_ref(),
                preprocessed_library.as_ref(),
                &batch_scorer,
                use_windowed_scoring,
            )?;

        // Store calibration from first file for subsequent files
        if file_idx == 0 {
            if let Some((cal, tol)) = calibration_result {
                log::info!(
                    "Storing calibration from first file for subsequent files (tolerance: {:.2} min)",
                    tol
                );
                shared_calibration = Some((cal, tol));
            }
            // Store preprocessed spectra from first file for potential reuse
            _all_preprocessed_spectra = file_preprocessed_spectra;
        }

        all_results.extend(file_results);
        all_spectra.extend(file_spectra);
    }

    log::info!(
        "Analysis complete. {} total regression results",
        all_results.len()
    );

    // Compute FDR and filter results (using spectral similarity scores)
    log::info!("Computing FDR using target-decoy competition with spectral scoring");
    let scored_entries = compute_fdr_and_filter(&library, &all_results, &all_spectra, &config)?;
    log::info!(
        "FDR filtering: {} entries at {}% FDR",
        scored_entries.len(),
        config.run_fdr * 100.0
    );

    // Write blib output
    log::info!("Writing blib to {}", config.output_blib.display());
    write_blib_output_with_scores(&config, &library, &scored_entries, &config.input_files)?;

    // Write output report if specified
    if let Some(report_path) = &config.output_report {
        log::info!("Writing report to {}", report_path.display());
        write_scored_report(report_path, &library, &scored_entries)?;
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
    matrix_builder: &DesignMatrixBuilder,
    solver: &RidgeSolver,
    config: &OspreyConfig,
    existing_calibration: Option<&(RTCalibration, f64)>,
    preprocessed_library: Option<&PreprocessedLibrary>,
    batch_scorer: &BatchScorer,
    use_windowed_scoring: bool,
) -> Result<(
    Vec<RegressionResult>,
    Vec<Spectrum>,
    Option<PreprocessedSpectra>,
    Option<(RTCalibration, f64)>,
)> {
    // Open mzML file and collect spectra
    let reader = MzmlReader::open(path)?;
    let spectra: Vec<Spectrum> = reader.filter_map(|r| r.ok()).collect();
    log::info!("Read {} MS2 spectra", spectra.len());

    if spectra.is_empty() {
        return Ok((Vec::new(), Vec::new(), None, None));
    }

    // Check if RT calibration is enabled
    if !config.rt_calibration.enabled {
        log::info!(
            "RT calibration disabled, using fallback tolerance: {:.2} min",
            config.rt_calibration.fallback_rt_tolerance
        );
        let results = process_spectra(
            &spectra,
            library,
            library_by_mz,
            matrix_builder,
            solver,
            config.rt_calibration.fallback_rt_tolerance,
            config.max_candidates_per_spectrum,
            None,
        )?;
        return Ok((results, spectra, None, None));
    }

    // Check if we have an existing calibration from a previous file
    if let Some((cal, tol)) = existing_calibration {
        log::info!(
            "Using calibration from first file (tolerance: {:.2} min)",
            tol
        );

        let results = process_spectra(
            &spectra,
            library,
            library_by_mz,
            matrix_builder,
            solver,
            *tol,
            config.max_candidates_per_spectrum,
            Some(cal),
        )?;
        return Ok((results, spectra, None, None)); // Don't return calibration/preprocessing for subsequent files
    }

    // First file: Run calibration discovery with ALL peptides
    log::info!("First file: Calibration Discovery (using all peptides)");

    // Run calibration using appropriate method based on mode and library size
    #[cfg(feature = "streaming")]
    let calibration_result = if config.streaming {
        // Use streaming pipeline for memory-efficient calibration
        log::info!("Using streaming mode for calibration discovery");
        run_calibration_discovery_streaming(path, library, config)
    } else if use_windowed_scoring {
        // Use windowed scoring for large libraries (memory-efficient)
        run_calibration_discovery_windowed(path, library, &spectra, config)
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
        run_calibration_discovery_windowed(path, library, &spectra, config)
    };

    #[cfg(not(feature = "streaming"))]
    let calibration_result = if use_windowed_scoring {
        // Use windowed scoring for large libraries (memory-efficient)
        run_calibration_discovery_windowed(path, library, &spectra, config)
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
        run_calibration_discovery_windowed(path, library, &spectra, config)
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

            // Save calibration to JSON if output path available
            if let Some(output_dir) = config.output_blib.parent() {
                let cal_filename = calibration_filename(
                    config.output_blib.file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("osprey")
                );
                let cal_path = output_dir.join(&cal_filename);
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

    // Calibrated full search
    log::info!("Calibrated Full Search");

    let results = process_spectra(
        &spectra,
        library,
        library_by_mz,
        matrix_builder,
        solver,
        rt_tolerance,
        config.max_candidates_per_spectrum,
        calibration_opt.as_ref(),
    )?;

    // Return calibration for first file to be reused
    // Note: preprocessed_spectra is only available when not using windowed scoring
    let calibration_to_share = calibration_opt.map(|cal| (cal, rt_tolerance));
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

    // Calculate initial wide tolerance based on library RT range
    let library_rts: Vec<f64> = library
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.retention_time)
        .collect();
    let min_rt = library_rts.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_rt = library_rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let rt_range = max_rt - min_rt;

    // Wide tolerance = fraction of RT range (default 25%)
    let initial_tolerance = rt_range * rt_config.initial_tolerance_fraction;

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
    log::info!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min range)",
        initial_tolerance,
        rt_config.initial_tolerance_fraction * 100.0,
        rt_range
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

    // RT tolerance (scaled only if mapping was applied)
    let mapped_tolerance = initial_tolerance * rt_slope.abs();

    if !use_direct_rts {
        log::info!(
            "Mapped RT tolerance: {:.2} min (in measured RT scale)",
            mapped_tolerance
        );
    }

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
    input_path: &std::path::Path,
    library: &[LibraryEntry],
    spectra: &[Spectrum],
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

    // Wide tolerance = fraction of RT range (default 25%)
    // Use the mzML RT range for tolerance calculation since that's the measured scale
    let initial_tolerance = mzml_rt_range * rt_config.initial_tolerance_fraction;

    let n_targets = library.iter().filter(|e| !e.is_decoy && !e.fragments.is_empty()).count();
    let n_decoys = library.iter().filter(|e| e.is_decoy && !e.fragments.is_empty()).count();

    log::info!(
        "Calibration: library has {} targets + {} decoys (library RT: {:.1}-{:.1}, mzML RT: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        lib_min_rt,
        lib_max_rt,
        mzml_min_rt,
        mzml_max_rt
    );

    // Sample library for calibration (pyXcorrDIA strategy)
    // For large libraries, only score a subset of peptides for calibration
    let calibration_library = if rt_config.calibration_sample_size > 0 {
        let sampled = sample_library_for_calibration(
            library,
            rt_config.calibration_sample_size,
            42, // Fixed seed for reproducibility
        );
        log::info!(
            "Sampled {} entries for calibration (from {} total)",
            sampled.len(),
            library.len()
        );
        sampled
    } else {
        log::info!("Using full library for calibration (calibration_sample_size=0)");
        library.to_vec()
    };

    log::info!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min mzML range)",
        initial_tolerance,
        rt_config.initial_tolerance_fraction * 100.0,
        mzml_rt_range
    );
    log::info!(
        "Fragment tolerance: {} {}",
        config.fragment_tolerance.tolerance,
        match config.fragment_tolerance.unit {
            osprey_core::ToleranceUnit::Ppm => "ppm",
            osprey_core::ToleranceUnit::Da => "Da",
        }
    );

    // Load MS1 spectra for precursor mass calibration
    // Extracts the M+0 isotope peak from actual MS1 spectra instead of using isolation window center
    let ms1_index = match load_ms1_spectra(input_path) {
        Ok(index) => {
            log::info!(
                "Loaded {} MS1 spectra for precursor calibration",
                index.len()
            );
            Some(index)
        }
        Err(e) => {
            log::warn!(
                "Could not load MS1 spectra: {}. Falling back to isolation window center for precursor m/z.",
                e
            );
            None
        }
    };

    // Run LibCosine calibration scoring (ppm-based peak matching, NO binning)
    // Uses MS1 spectra to extract M+0 isotope peak for accurate precursor mass calibration
    let matches = if let Some(ref ms1) = ms1_index {
        run_libcosine_calibration_scoring_with_ms1(
            &calibration_library,
            spectra,
            Some(&MS1IndexWrapper(ms1)),
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
    let ms1_count_all = matches.iter().filter(|m| m.ms1_ppm_error.is_some()).count();
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

    // Run target-decoy competition (pyXcorrDIA style)
    // Each target competes with its paired decoy. Higher score wins. Ties go to decoy.
    // Then sort winners by score and compute FDR = cumulative_decoys / cumulative_targets
    let calibration_fdr = 0.01; // 1% FDR for calibration
    let fdr_controller = FdrController::new(calibration_fdr);

    // Convert matches to competition format: (item, score, is_decoy, entry_id)
    let competition_input = matches.iter().map(|m| {
        (m.clone(), m.score, m.is_decoy, m.entry_id)
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
    let mut mz_qc_data = MzQCData::new();

    for m in &competition_result.passing_targets {
        library_rts_detected.push(m.library_rt);
        measured_rts_detected.push(m.measured_rt);

        // Collect MS1 PPM error for mass calibration
        if let Some(ms1_error) = m.ms1_ppm_error {
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
        "MS1 calibration: mean={:.2} ppm, SD={:.2} ppm (from {} observations)",
        ms1_calibration.mean,
        ms1_calibration.sd,
        mz_qc_data.n_ms1()
    );

    // Build full CalibrationParams
    let calibration_params = CalibrationParams {
        metadata: CalibrationMetadata {
            num_confident_peptides,
            num_sampled_precursors: matches.len(),
            calibration_successful: true,
            timestamp: chrono::Utc::now().to_rfc3339(),
        },
        ms1_calibration,
        ms2_calibration,
        rt_calibration: RTCalibrationParams {
            method: RTCalibrationMethod::LOESS,
            residual_sd: rt_stats.residual_std,
            n_points: rt_stats.n_points,
            r_squared: rt_stats.r_squared,
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

    // Wide tolerance = fraction of RT range (default 25%)
    let initial_tolerance = rt_range * rt_config.initial_tolerance_fraction;

    let n_targets = library.iter().filter(|e| !e.is_decoy && !e.fragments.is_empty()).count();
    let n_decoys = library.iter().filter(|e| e.is_decoy && !e.fragments.is_empty()).count();

    log::info!(
        "RT calibration (cached): scoring {} targets + {} decoys (RT range: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        min_rt,
        max_rt
    );
    log::info!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min range)",
        initial_tolerance,
        rt_config.initial_tolerance_fraction * 100.0,
        rt_range
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

    // RT tolerance (scaled only if mapping was applied)
    let mapped_tolerance = initial_tolerance * rt_slope.abs();

    if !use_direct_rts {
        log::info!(
            "Mapped RT tolerance: {:.2} min (in measured RT scale)",
            mapped_tolerance
        );
    }

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
        },
        ms1_calibration: MzCalibration::uncalibrated(),
        ms2_calibration: MzCalibration::uncalibrated(),
        rt_calibration: RTCalibrationParams {
            method: RTCalibrationMethod::LOESS,
            residual_sd: rt_stats.residual_std,
            n_points: rt_stats.n_points,
            r_squared: rt_stats.r_squared,
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

/// Process spectra with optional RT calibration
fn process_spectra(
    spectra: &[Spectrum],
    library: &[LibraryEntry],
    library_by_mz: &HashMap<i32, Vec<usize>>,
    matrix_builder: &DesignMatrixBuilder,
    solver: &RidgeSolver,
    rt_tolerance: f64,
    max_candidates: usize,
    calibration: Option<&RTCalibration>,
) -> Result<Vec<RegressionResult>> {
    // Set up progress bar
    let pb = ProgressBar::new(spectra.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} spectra",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Process spectra in parallel
    let results: Vec<RegressionResult> = spectra
        .par_iter()
        .map(|spectrum| {
            let result = process_spectrum_with_filter(
                spectrum,
                library,
                library_by_mz,
                matrix_builder,
                solver,
                rt_tolerance,
                max_candidates,
                None, // No library filter
                calibration,
            );
            pb.inc(1);
            result
        })
        .filter_map(|r| r.ok())
        .filter(|r| !r.coefficients.is_empty())
        .collect();

    pb.finish_with_message("Done");

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
    /// Coefficient time series
    rt_coef_pairs: Vec<(f64, f64)>,
    /// Extracted features
    features: FeatureSet,
    /// Score (peak_apex for now)
    score: f64,
    /// Q-value from target-decoy competition
    q_value: f64,
    /// Is decoy
    is_decoy: bool,
}

/// Compute FDR using target-decoy competition and filter results
///
/// Uses spectral similarity scoring (LibCosine) as the primary score for FDR computation.
fn compute_fdr_and_filter(
    library: &[LibraryEntry],
    results: &[RegressionResult],
    spectra: &[Spectrum],
    config: &OspreyConfig,
) -> Result<Vec<ScoredEntry>> {
    use osprey_scoring::SpectralScorer;

    // Build ID-to-index map
    let id_to_index: HashMap<u32, usize> = library
        .iter()
        .enumerate()
        .map(|(idx, entry)| (entry.id, idx))
        .collect();

    // Aggregate results by library entry ID
    let mut entry_data: HashMap<u32, Vec<(f64, f64)>> = HashMap::new(); // (RT, coef)
    for result in results {
        for (lib_id, coef) in result.library_ids.iter().zip(result.coefficients.iter()) {
            entry_data
                .entry(*lib_id)
                .or_default()
                .push((result.retention_time, *coef));
        }
    }

    // Extract features and create scored entries (parallelized)
    let feature_extractor = FeatureExtractor::new();
    let spectral_scorer = SpectralScorer::new();
    let peak_detector = PeakDetector::new().with_min_height(0.05);

    let mut scored_entries: Vec<ScoredEntry> = entry_data
        .par_iter()
        .filter_map(|(lib_id, data_points)| {
            if data_points.len() < 3 {
                return None;
            }

            let lib_idx = id_to_index.get(lib_id).copied()?;
            let entry = &library[lib_idx];

            // Sort by RT
            let mut sorted_data = data_points.clone();
            sorted_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

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
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(rt, _)| *rt)
                .unwrap_or(entry.retention_time);

            // Get RT range from coefficient series for peak region filtering
            let rt_min = rt_coef_pairs.iter().map(|(rt, _)| *rt).fold(f64::INFINITY, f64::min);
            let rt_max = rt_coef_pairs.iter().map(|(rt, _)| *rt).fold(f64::NEG_INFINITY, f64::max);

            // Filter spectra to those in the peak region:
            // 1. Isolation window contains the peptide's precursor m/z
            // 2. RT is within the coefficient series range
            let peak_region_spectra: Vec<Spectrum> = spectra
                .iter()
                .filter(|s| {
                    s.contains_precursor(entry.precursor_mz)
                        && s.retention_time >= rt_min
                        && s.retention_time <= rt_max
                })
                .cloned()
                .collect();

            // Extract features using spectrum aggregation (FR-5.2.1)
            // This aggregates all spectra in the peak region weighted by coefficients
            let features = feature_extractor.extract_with_aggregation(
                entry,
                &rt_coef_pairs,
                &peak_region_spectra,
                Some(entry.retention_time),
            );

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

            Some(ScoredEntry {
                lib_idx,
                rt_coef_pairs,
                features,
                score,
                q_value: 1.0, // Will be computed below
                is_decoy,
            })
        })
        .collect();

    // Separate targets and decoys
    let target_scores: Vec<f64> = scored_entries
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.score)
        .collect();

    let decoy_scores: Vec<f64> = scored_entries
        .iter()
        .filter(|e| e.is_decoy)
        .map(|e| e.score)
        .collect();

    log::info!(
        "FDR computation: {} targets, {} decoys",
        target_scores.len(),
        decoy_scores.len()
    );

    // Compute q-values using target-decoy competition
    let fdr_controller = FdrController::new(config.run_fdr);
    let target_qvalues = fdr_controller.compute_qvalues(&target_scores, &decoy_scores)?;

    // Assign q-values back to entries
    let mut target_idx = 0;
    for entry in scored_entries.iter_mut() {
        if !entry.is_decoy {
            if target_idx < target_qvalues.len() {
                entry.q_value = target_qvalues[target_idx];
                target_idx += 1;
            }
        }
        // Decoys keep q_value = 1.0
    }

    // Report FDR statistics
    let passing_targets = scored_entries
        .iter()
        .filter(|e| !e.is_decoy && e.q_value <= config.run_fdr)
        .count();

    log::info!(
        "Targets passing {}% FDR: {}",
        config.run_fdr * 100.0,
        passing_targets
    );

    // Filter to passing targets only for output
    let filtered: Vec<ScoredEntry> = scored_entries
        .into_iter()
        .filter(|e| !e.is_decoy && e.q_value <= config.run_fdr)
        .collect();

    Ok(filtered)
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
            PsmFeatures {
                psm_id: format!("{}_{}", file_name, idx),
                peptide: lib_entry.modified_sequence.clone(),
                proteins: lib_entry.protein_ids.clone(),
                scan_number: entry.rt_coef_pairs
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .map(|(_, _)| 0u32)
                    .unwrap_or(0),
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

    log::info!("Wrote {} PSMs to PIN file for mokapot", psm_features.len());

    // Run mokapot
    let results = mokapot.run(&pin_file, &mokapot_dir)?;

    // Update q-values from mokapot results
    let result_map: std::collections::HashMap<String, f64> = results
        .into_iter()
        .map(|r| (r.psm_id, r.q_value))
        .collect();

    for (idx, entry) in scored_entries.iter_mut().enumerate() {
        let psm_id = format!("{}_{}", file_name, idx);
        if let Some(&q_value) = result_map.get(&psm_id) {
            entry.q_value = q_value;
        }
    }

    log::info!("Updated q-values from mokapot results");

    Ok(())
}

/// Run per-file and experiment-level FDR control
///
/// This implements the two-level FDR strategy:
/// 1. Per-file FDR: Each file is processed independently
/// 2. Experiment-level FDR: Combine results across all files
#[allow(dead_code)]
fn run_fdr_control(
    per_file_results: &mut [(String, Vec<ScoredEntry>)],
    library: &[LibraryEntry],
    run_fdr: f64,
    experiment_fdr: f64,
    use_mokapot: bool,
    output_dir: &std::path::Path,
) -> Result<Vec<ScoredEntry>> {
    // Step 1: Per-file FDR
    log::info!("Step 1: Per-file FDR control at {}%", run_fdr * 100.0);

    for (file_name, entries) in per_file_results.iter_mut() {
        if use_mokapot {
            compute_fdr_with_mokapot(
                library,
                entries,
                file_name,
                output_dir,
                run_fdr,
                run_fdr,
            )?;
        }

        // Filter to per-file FDR threshold
        let passing = entries.iter().filter(|e| !e.is_decoy && e.q_value <= run_fdr).count();
        log::info!("  {}: {} targets passing {}% FDR", file_name, passing, run_fdr * 100.0);
    }

    // Step 2: Experiment-level FDR
    log::info!("Step 2: Experiment-level FDR control at {}%", experiment_fdr * 100.0);

    // Combine all per-file results
    let mut all_entries: Vec<ScoredEntry> = per_file_results
        .iter()
        .flat_map(|(_, entries)| entries.iter().cloned())
        .collect();

    // Re-compute experiment-level q-values
    let fdr_controller = FdrController::new(experiment_fdr);

    let target_scores: Vec<f64> = all_entries
        .iter()
        .filter(|e| !e.is_decoy)
        .map(|e| e.score)
        .collect();

    let decoy_scores: Vec<f64> = all_entries
        .iter()
        .filter(|e| e.is_decoy)
        .map(|e| e.score)
        .collect();

    let experiment_qvalues = fdr_controller.compute_qvalues(&target_scores, &decoy_scores)?;

    // Update q-values with experiment-level values
    let mut target_idx = 0;
    for entry in all_entries.iter_mut() {
        if !entry.is_decoy {
            if target_idx < experiment_qvalues.len() {
                entry.q_value = experiment_qvalues[target_idx];
                target_idx += 1;
            }
        }
    }

    // Filter to experiment-level FDR
    let filtered: Vec<ScoredEntry> = all_entries
        .into_iter()
        .filter(|e| !e.is_decoy && e.q_value <= experiment_fdr)
        .collect();

    log::info!(
        "Experiment-level: {} targets passing {}% FDR",
        filtered.len(),
        experiment_fdr * 100.0
    );

    Ok(filtered)
}

/// Write scored results to blib format for Skyline
fn write_blib_output_with_scores(
    config: &OspreyConfig,
    library: &[LibraryEntry],
    scored_entries: &[ScoredEntry],
    input_files: &[std::path::PathBuf],
) -> Result<()> {
    let writer = BlibWriter::create(&config.output_blib)?;

    // Add metadata
    writer.add_metadata("osprey_version", env!("CARGO_PKG_VERSION"))?;
    writer.add_metadata("rt_calibration_enabled", &config.rt_calibration.enabled.to_string())?;
    writer.add_metadata("run_fdr", &config.run_fdr.to_string())?;
    writer.add_metadata("fdr_method", "target_decoy_competition")?;

    // Add source files
    let mut file_ids = std::collections::HashMap::new();
    for input_file in input_files {
        let file_name = input_file.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        let file_id = writer.add_source_file(file_name)?;
        file_ids.insert(input_file.clone(), file_id);
    }

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
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .map(|(rt, c)| (*rt, *c))
            .unwrap_or((0.0, 0.0));

        let file_id = *file_ids.values().next().unwrap_or(&1);

        // Get fragment peaks from library entry
        let mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
        let intensities: Vec<f32> = entry.fragments.iter().map(|f| f.relative_intensity).collect();

        // Use 1 - q_value as the score (higher is better, blib convention)
        let blib_score = 1.0 - scored.q_value;

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

        // Add peak boundaries
        let boundaries = osprey_core::PeakBoundaries {
            start_rt,
            end_rt,
            apex_rt,
            apex_coefficient: scored.score,
            integrated_area: scored.features.peak_area,
            peak_quality: osprey_core::PeakQuality::default(),
        };

        let file_name = input_files.first()
            .and_then(|p| p.file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        writer.add_peak_boundaries(ref_id, file_name, &boundaries)?;

        // Add run scores
        writer.add_run_scores(ref_id, file_name, blib_score, scored.score, scored.q_value)?;
    }

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
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
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
            scored.q_value
        )?;
    }

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
    let writer = BlibWriter::create(&config.output_blib)?;

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
        sorted_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

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
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
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

/// Process a single spectrum with optional library filtering and RT calibration
fn process_spectrum_with_filter(
    spectrum: &Spectrum,
    library: &[LibraryEntry],
    library_by_mz: &HashMap<i32, Vec<usize>>,
    matrix_builder: &DesignMatrixBuilder,
    solver: &RidgeSolver,
    rt_tolerance: f64,
    max_candidates: usize,
    library_filter: Option<&std::collections::HashSet<usize>>,
    calibration: Option<&RTCalibration>,
) -> Result<RegressionResult> {
    // Select candidates based on isolation window and RT
    let candidates = select_candidates_with_calibration(
        spectrum,
        library,
        library_by_mz,
        rt_tolerance,
        max_candidates,
        library_filter,
        calibration,
    );

    if candidates.is_empty() {
        return Ok(RegressionResult::new(
            spectrum.scan_number,
            spectrum.retention_time,
        ));
    }

    // Build design matrix
    let candidate_refs: Vec<&LibraryEntry> = candidates.iter().map(|&i| &library[i]).collect();
    let (design_matrix, library_ids) = matrix_builder.build_with_indices(&candidate_refs);

    // Bin observed spectrum
    let observed = matrix_builder.bin_observed(spectrum);

    // Solve regression
    let coefficients = solver.solve_nonnegative(&design_matrix, &observed, None)?;

    // Filter to non-zero coefficients
    let mut result = RegressionResult::new(spectrum.scan_number, spectrum.retention_time);
    for (id, coef) in library_ids.iter().zip(coefficients.iter()) {
        if *coef > 1e-6 {
            result.library_ids.push(*id);
            result.coefficients.push(*coef);
        }
    }

    // Compute residual
    let predicted = design_matrix.dot(&coefficients);
    let residual = (&observed - &predicted).mapv(|x| x * x).sum();
    result.residual = residual;

    Ok(result)
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

                if (expected_rt - spectrum.retention_time).abs() > rt_tolerance {
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
            rt_a.partial_cmp(&rt_b).unwrap()
        });
        candidates.truncate(max_candidates);
    }

    candidates
}

/// Filter library to only include peptides found in any file
///
/// Creates a filtered library containing only the union of peptides
/// that passed FDR in any of the input files.
#[allow(dead_code)]
fn filter_library_by_detections(
    library: &[LibraryEntry],
    detected_lib_indices: &[usize],
) -> Vec<LibraryEntry> {
    // Create unique set of detected indices
    let detected_set: std::collections::HashSet<usize> = detected_lib_indices.iter().copied().collect();

    // Filter library
    library
        .iter()
        .enumerate()
        .filter(|(idx, _)| detected_set.contains(idx))
        .map(|(_, entry)| entry.clone())
        .collect()
}

/// Run refinement search with filtered library
///
/// This is Step 6 of the workflow:
/// Repeat the search (without RT calibration) on the library filtered
/// by the union of peptides found from the per-file analysis.
///
/// Benefits:
/// - Faster because library is smaller
/// - May improve detection in files with poor initial results
/// - Can improve quantification by searching with focused candidates
#[allow(dead_code)]
fn run_refinement_search(
    input_files: &[std::path::PathBuf],
    filtered_library: &[LibraryEntry],
    config: &OspreyConfig,
) -> Result<Vec<RegressionResult>> {
    log::info!(
        "Starting refinement search with {} peptides (filtered library)",
        filtered_library.len()
    );

    if filtered_library.is_empty() {
        log::warn!("Filtered library is empty, skipping refinement search");
        return Ok(Vec::new());
    }

    // Set up binning based on resolution mode
    let bin_config = match config.resolution_mode {
        ResolutionMode::UnitResolution => BinConfig::unit_resolution(),
        ResolutionMode::HRAM => BinConfig::hram(),
        ResolutionMode::Auto => BinConfig::unit_resolution(),
    };
    let binner = Binner::new(bin_config);
    let matrix_builder = DesignMatrixBuilder::new(binner.clone());

    // Get regularization lambda
    let lambda = match &config.regularization_lambda {
        osprey_core::RegularizationSetting::Fixed(l) => *l,
        _ => 1.0,
    };
    let solver = RidgeSolver::new(lambda);

    // Build filtered library index
    let library_by_mz = build_mz_index(filtered_library);

    // Get RT tolerance from config (fixed, no calibration in refinement)
    let rt_tolerance = config.rt_calibration.rt_tolerance_factor * 2.0; // Use wider tolerance

    let mut all_results: Vec<RegressionResult> = Vec::new();

    // Process each file
    for path in input_files {
        let file_name = path.file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");

        log::info!("Refinement search: processing {}", file_name);

        // Open mzML file
        let reader = MzmlReader::open(path)?;
        let spectra: Vec<Spectrum> = reader.filter_map(|r| r.ok()).collect();

        if spectra.is_empty() {
            continue;
        }

        // Process spectra with fixed RT tolerance (no calibration)
        let results = process_spectra(
            &spectra,
            filtered_library,
            &library_by_mz,
            &matrix_builder,
            &solver,
            rt_tolerance,
            config.max_candidates_per_spectrum,
            None, // No calibration for refinement
        )?;

        log::info!(
            "Refinement search {}: {} results from {} spectra",
            file_name,
            results.len(),
            spectra.len()
        );

        all_results.extend(results);
    }

    log::info!(
        "Refinement search complete: {} total results",
        all_results.len()
    );

    Ok(all_results)
}

/// Get detected library indices from scored entries
#[allow(dead_code)]
fn get_detected_indices(scored_entries: &[ScoredEntry]) -> Vec<usize> {
    scored_entries.iter().map(|e| e.lib_idx).collect()
}

/// Run complete two-pass analysis with optional refinement
///
/// This implements the full workflow:
/// 1. Initial search with RT calibration
/// 2. FDR control (per-file and experiment-level)
/// 3. Optional refinement search with filtered library
#[allow(dead_code)]
fn run_analysis_with_refinement(
    config: &OspreyConfig,
    enable_refinement: bool,
) -> Result<Vec<RegressionResult>> {
    // Run initial analysis
    let initial_results = run_analysis(config.clone())?;

    if !enable_refinement {
        return Ok(initial_results);
    }

    // To do refinement, we need the scored entries
    // For now, log that refinement would happen
    log::info!("Refinement search is enabled but requires scored entries from initial pass");
    log::info!("Future: refinement would filter library and re-search");

    // The full refinement workflow would be:
    // 1. Get detected indices from scored_entries
    // 2. Filter library by detected peptides
    // 3. Run refinement search on filtered library
    // 4. Merge results and re-compute FDR

    Ok(initial_results)
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
    // Note: sorted by winning_libcosine descending to match FDR calculation order
    writeln!(
        writer,
        "target_entry_id,charge,target_sequence,decoy_sequence,\
         winning_libcosine,target_libcosine,decoy_libcosine,\
         winning_xcorr,target_xcorr,decoy_xcorr,\
         target_isotope_score,decoy_isotope_score,\
         target_precursor_error_ppm,decoy_precursor_error_ppm,\
         target_rt,decoy_rt,library_rt,expected_rt,target_delta_rt,decoy_delta_rt,\
         target_matched_frags,decoy_matched_frags,target_wins"
    ).map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

    // Write each paired result (already sorted by winning_libcosine descending)
    for p in &paired {
        writeln!(
            writer,
            "{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{},{},{},{},{:.2},{:.2},{:.2},{},{:.2},{:.2},{},{},{}",
            p.target_entry_id,
            p.charge,
            p.target_sequence,
            p.decoy_sequence,
            p.winning_libcosine,
            p.target_libcosine,
            p.decoy_libcosine,
            p.winning_xcorr,
            p.target_xcorr,
            p.decoy_xcorr,
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

    // Wide tolerance = fraction of RT range
    let initial_tolerance = rt_range * rt_config.initial_tolerance_fraction;

    let n_targets = library.iter().filter(|e| !e.is_decoy && !e.fragments.is_empty()).count();
    let n_decoys = library.iter().filter(|e| e.is_decoy && !e.fragments.is_empty()).count();

    log::info!(
        "Streaming calibration: scoring {} targets + {} decoys (RT range: {:.1}-{:.1} min)",
        n_targets,
        n_decoys,
        min_rt,
        max_rt
    );

    // Sample library for calibration if configured
    let calibration_library = if rt_config.calibration_sample_size > 0 {
        let sampled = sample_library_for_calibration(
            library,
            rt_config.calibration_sample_size,
            42,
        );
        log::info!(
            "Sampled {} entries for streaming calibration",
            sampled.len()
        );
        sampled
    } else {
        library.to_vec()
    };

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

    // Run target-decoy competition
    let calibration_fdr = 0.01;
    let fdr_controller = FdrController::new(calibration_fdr);

    let competition_input = matches.iter().map(|m| {
        (m.clone(), m.score, m.is_decoy, m.entry_id)
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
        },
        ms1_calibration: MzCalibration::uncalibrated(),
        ms2_calibration: MzCalibration::uncalibrated(),
        rt_calibration: RTCalibrationParams {
            method: RTCalibrationMethod::LOESS,
            residual_sd: rt_stats.residual_std,
            n_points: rt_stats.n_points,
            r_squared: rt_stats.r_squared,
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
