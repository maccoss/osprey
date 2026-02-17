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

use arrow::array::{
    ArrayRef, BinaryBuilder, BooleanBuilder, Float64Builder, StringBuilder, UInt32Builder,
    UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use indicatif::{ProgressBar, ProgressStyle};
use osprey_chromatography::{
    calculate_mz_calibration, calibration_path_for_input, load_calibration, save_calibration,
    trapezoidal_area, CalibrationMetadata, CalibrationParams, IsolationScheme, MzQCData,
    RTCalibration, RTCalibrationMethod, RTCalibrationParams, RTCalibrator, RTCalibratorConfig,
};
use osprey_core::{
    CoelutionFeatureSet, CoelutionScoredEntry, DecoyMethod as CoreDecoyMethod, FdrMethod,
    FragmentToleranceConfig, LibraryEntry, LibraryFragment, MS1Spectrum, OspreyConfig, OspreyError,
    Result, Spectrum, ToleranceUnit, XICPeakBounds,
};
use osprey_fdr::{
    get_pin_feature_names, percolator, pin_feature_value, MokapotRunner, NUM_PIN_FEATURES,
};
use osprey_io::{load_all_spectra, load_library, BlibWriter, MS1Index};
use osprey_scoring::{
    batch::{run_coelution_calibration_scoring, sample_library_for_calibration, MS1SpectrumLookup},
    DecoyGenerator, DecoyMethod, Enzyme, SpectralScorer,
};

/// Wrapper to implement MS1SpectrumLookup for MS1Index
struct MS1IndexWrapper<'a>(&'a MS1Index);

impl<'a> MS1SpectrumLookup for MS1IndexWrapper<'a> {
    fn find_nearest(&self, retention_time: f64) -> Option<&MS1Spectrum> {
        self.0.find_nearest(retention_time)
    }
}
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};

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

    // Set up XCorr scorer with appropriate binning for calibration LDA
    let is_hram = matches!(config.resolution_mode, osprey_core::ResolutionMode::HRAM);
    let xcorr_scorer = if is_hram {
        SpectralScorer::hram().with_tolerance_ppm(config.fragment_tolerance.tolerance)
    } else {
        SpectralScorer::new().with_tolerance_da(config.fragment_tolerance.tolerance)
    };

    // Calibration sampling with retry loop
    // Attempt 1: sample calibration_sample_size targets
    // Attempt 2: expand by retry_factor (default 2×)
    // Attempt 3: use ALL library entries (guaranteed fallback)
    // Matches accumulate across attempts — FDR runs on the combined set.
    // If the library has fewer targets than requested, we use all and skip retries.
    let sample_size = rt_config.calibration_sample_size;
    let retry_factor = rt_config.calibration_retry_factor;
    let n_total_targets = library.iter().filter(|e| !e.is_decoy).count();
    let max_attempts: usize =
        if sample_size > 0 && retry_factor > 1.0 && n_total_targets > sample_size {
            3
        } else {
            1
        };
    let mut current_sample_size = sample_size;

    log::info!(
        "Calibration: library has {} targets, requesting {} per attempt ({} attempt{} max)",
        n_total_targets,
        if current_sample_size == 0 || n_total_targets <= current_sample_size {
            "all".to_string()
        } else {
            format!("{}", current_sample_size)
        },
        max_attempts,
        if max_attempts == 1 { "" } else { "s" }
    );

    // Accumulate best match per entry across all attempts
    let mut accumulated_matches: HashMap<u32, CalibrationMatch> = HashMap::new();

    for attempt in 1..=max_attempts {
        let calibration_library =
            sample_library_for_calibration(library, current_sample_size, 42 + attempt as u64);

        let n_sampled_targets = calibration_library.iter().filter(|e| !e.is_decoy).count();
        let n_sampled_decoys = calibration_library.len() - n_sampled_targets;
        let used_all = calibration_library.len() == library.len();

        log::info!(
            "Calibration attempt {}/{}: {} targets + {} decoys{}",
            attempt,
            max_attempts,
            n_sampled_targets,
            n_sampled_decoys,
            if used_all {
                " (entire library)".to_string()
            } else {
                format!(
                    " ({:.0}% of library)",
                    n_sampled_targets as f64 / n_total_targets.max(1) as f64 * 100.0
                )
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
                None, // First pass: use library RT directly
                Some(&xcorr_scorer),
            )
        } else {
            run_coelution_calibration_scoring::<MS1IndexWrapper>(
                &calibration_library,
                spectra,
                None,
                config.fragment_tolerance,
                config.precursor_tolerance.tolerance,
                initial_tolerance,
                None, // First pass: use library RT directly
                Some(&xcorr_scorer),
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
        // Don't use isotope feature during calibration — MS1 isotope scoring
        // is only reliable after mass calibration, which hasn't happened yet.
        // Isotope features are used in the full scoring step after calibration.
        let calibration_fdr = 0.01;

        // Convert to Vec for LDA training
        let mut all_matches: Vec<CalibrationMatch> =
            accumulated_matches.values().cloned().collect();

        // Train LDA and score (4 features: correlation, libcosine, top6, snr)
        let _n_passing =
            osprey_scoring::calibration_ml::train_and_score_calibration(&mut all_matches, false)?;

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

        // Compute median peak width from confident matches for adaptive co-elution window
        // Log median peak width from calibration matches (diagnostic)
        {
            let mut widths: Vec<f64> = passing_targets
                .iter()
                .filter_map(|m| m.peak_width_minutes)
                .filter(|&w| w > 0.0 && w.is_finite())
                .collect();
            if widths.len() >= 10 {
                widths.sort_by(|a, b| a.total_cmp(b));
                let mid = widths.len() / 2;
                let median = if widths.len() % 2 == 0 {
                    (widths[mid - 1] + widths[mid]) / 2.0
                } else {
                    widths[mid]
                };
                log::info!(
                    "Measured peak width: median={:.2} min ({:.0} sec) from {} peaks",
                    median,
                    median * 60.0,
                    widths.len()
                );
            }
        }

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
            if attempt < max_attempts && !used_all {
                // Not the final attempt and haven't used all entries yet - retry with more
                // Determine next sample size
                if attempt + 1 == max_attempts {
                    // Final attempt always uses ALL library entries
                    current_sample_size = 0;
                } else {
                    let new_size = (current_sample_size as f64 * retry_factor) as usize;
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
            outlier_retention: 0.8, // Remove worst 20% of calibration points
        };
        let calibrator = RTCalibrator::with_config(calibrator_config);
        let mut rt_calibration = calibrator.fit(&library_rts_detected, &measured_rts_detected)?;
        let mut rt_stats = rt_calibration.stats();

        // === Iterative calibration refinement (2-pass, DIA-NN-style) ===
        // Compute first-pass tolerance and check if refinement is worthwhile
        let pass1_tolerance = if config.rt_calibration.use_percentile_tolerance {
            rt_calibration.percentile_tolerance(0.20)
        } else {
            rt_stats.residual_std * config.rt_calibration.rt_tolerance_factor
        };
        let pass1_tolerance = pass1_tolerance
            .max(config.rt_calibration.min_rt_tolerance)
            .min(config.rt_calibration.max_rt_tolerance);

        // Only refine if tolerance narrowed significantly (at least 2× tighter)
        if pass1_tolerance < initial_tolerance * 0.5 {
            log::info!(
                "Calibration refinement: re-scoring with {:.2} min tolerance (was {:.1} min)",
                pass1_tolerance,
                initial_tolerance
            );

            let predict_fn = |lib_rt: f64| -> f64 { rt_calibration.predict(lib_rt) };

            let refined_matches = if has_ms1 {
                run_coelution_calibration_scoring(
                    &calibration_library,
                    spectra,
                    Some(&MS1IndexWrapper(ms1_index)),
                    config.fragment_tolerance,
                    config.precursor_tolerance.tolerance,
                    pass1_tolerance,
                    Some(&predict_fn),
                    Some(&xcorr_scorer),
                )
            } else {
                run_coelution_calibration_scoring::<MS1IndexWrapper>(
                    &calibration_library,
                    spectra,
                    None,
                    config.fragment_tolerance,
                    config.precursor_tolerance.tolerance,
                    pass1_tolerance,
                    Some(&predict_fn),
                    Some(&xcorr_scorer),
                )
            };

            // Re-run LDA scoring on refined matches
            let mut refined_all: Vec<CalibrationMatch> = refined_matches;
            let _n_passing_refined = osprey_scoring::calibration_ml::train_and_score_calibration(
                &mut refined_all,
                false,
            )?;

            // Filter to passing targets
            let refined_passing: Vec<&CalibrationMatch> = refined_all
                .iter()
                .filter(|m| {
                    !m.is_decoy
                        && m.q_value <= calibration_fdr
                        && m.signal_to_noise >= MIN_SNR_FOR_RT_CAL
                })
                .collect();

            let n_refined = refined_passing.len();

            if n_refined >= ABSOLUTE_MIN_CALIBRATION_POINTS {
                // Refit LOESS on refined points
                let refined_lib_rts: Vec<f64> =
                    refined_passing.iter().map(|m| m.library_rt).collect();
                let refined_meas_rts: Vec<f64> =
                    refined_passing.iter().map(|m| m.measured_rt).collect();

                let rt_cal_refined = calibrator.fit(&refined_lib_rts, &refined_meas_rts)?;
                let rt_stats_refined = rt_cal_refined.stats();

                log::info!(
                    "Refined RT calibration: {} points, R²={:.4}, residual_SD={:.3} min (was {} points, R²={:.4}, {:.3} min)",
                    rt_stats_refined.n_points,
                    rt_stats_refined.r_squared,
                    rt_stats_refined.residual_std,
                    rt_stats.n_points,
                    rt_stats.r_squared,
                    rt_stats.residual_std
                );

                // Accept refined calibration if R² didn't degrade significantly
                if rt_stats_refined.r_squared >= rt_stats.r_squared * 0.99 {
                    rt_calibration = rt_cal_refined;
                    rt_stats = rt_stats_refined;

                    // Re-collect mass errors from refined matches
                    mz_qc_data = MzQCData::new(config.fragment_tolerance.unit);
                    for m in &refined_passing {
                        if let Some(ms1_error) = m.ms1_error {
                            mz_qc_data.add_ms1_error(ms1_error);
                        }
                        for &ms2_error in &m.ms2_mass_errors {
                            mz_qc_data.add_ms2_error(ms2_error);
                        }
                    }
                } else {
                    log::info!(
                        "Refined calibration not better (R² {:.4} vs {:.4}), keeping original",
                        rt_stats_refined.r_squared,
                        rt_stats.r_squared
                    );
                }
            } else {
                log::info!(
                    "Refinement pass: only {} points (need {}), keeping original calibration",
                    n_refined,
                    ABSOLUTE_MIN_CALIBRATION_POINTS
                );
            }
        }

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
                p20_abs_residual: Some(rt_stats.p20_abs_residual),
            },
        };

        return Ok((rt_calibration, calibration_params));
    }

    // Should not reach here (loop always returns or errors)
    unreachable!("Calibration retry loop exited without result")
}

/// Reverse m/z calibration: find where a theoretical m/z appears in uncalibrated data.
///
/// Since forward calibration corrects: corrected = observed - offset,
/// the reverse gives: observed ≈ theoretical + offset.
/// This avoids pre-calibrating all MS1 spectra (which would clone GBs of data).
fn reverse_calibrate_mz(
    theoretical_mz: f64,
    calibration: &osprey_chromatography::calibration::MzCalibration,
) -> f64 {
    if !calibration.calibrated {
        return theoretical_mz;
    }
    if calibration.unit == "Th" {
        theoretical_mz + calibration.mean
    } else {
        // PPM: offset_da = observed * mean / 1e6 ≈ theoretical * mean / 1e6
        theoretical_mz * (1.0 + calibration.mean / 1e6)
    }
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
/// Find the PIN feature that best separates targets from decoys using Cohen's d.
///
/// Returns `(feature_index, higher_is_better)` or `None` if separation cannot be computed.
/// Compute ROC AUC via the Mann-Whitney U statistic.
///
/// AUC > 0.5 means higher values are more target-like (targets rank higher).
/// AUC < 0.5 means lower values are more target-like (targets rank lower).
/// AUC = 0.5 means no discrimination.
fn compute_roc_auc(target_vals: &[f64], decoy_vals: &[f64]) -> f64 {
    let n_t = target_vals.len();
    let n_d = decoy_vals.len();
    if n_t == 0 || n_d == 0 {
        return 0.5;
    }

    // Combine all values with labels (true = target)
    let mut combined: Vec<(f64, bool)> = Vec::with_capacity(n_t + n_d);
    for &v in target_vals {
        combined.push((v, true));
    }
    for &v in decoy_vals {
        combined.push((v, false));
    }

    // Sort ascending by value; ties broken arbitrarily (doesn't affect average rank)
    combined.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Assign average ranks for tied groups
    let n = combined.len();
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j < n && combined[j].0 == combined[i].0 {
            j += 1;
        }
        // Average rank for positions i..j (1-indexed)
        let avg_rank = (i + 1 + j) as f64 / 2.0;
        for rank in ranks.iter_mut().take(j).skip(i) {
            *rank = avg_rank;
        }
        i = j;
    }

    // Sum of target ranks
    let target_rank_sum: f64 = ranks
        .iter()
        .zip(combined.iter())
        .filter(|(_, (_, is_target))| *is_target)
        .map(|(r, _)| *r)
        .sum();

    // U = target_rank_sum - n_t*(n_t+1)/2
    let u = target_rank_sum - (n_t * (n_t + 1)) as f64 / 2.0;
    u / (n_t as f64 * n_d as f64)
}
// =============================================================================
// Per-file score parquet caching
// =============================================================================

/// Derive the scores parquet path from an input mzML path.
/// e.g., `/data/sample1.mzML` → `/data/sample1.scores.parquet`
fn scores_path_for_input(input_path: &std::path::Path) -> std::path::PathBuf {
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    let parent = input_path.parent().unwrap_or(std::path::Path::new("."));
    parent.join(format!("{}.scores.parquet", stem))
}

/// Write per-file scored entries to parquet (pre-FDR, includes fragments and XIC for blib).
fn write_scores_parquet(path: &std::path::Path, entries: &[CoelutionScoredEntry]) -> Result<()> {
    use parquet::arrow::ArrowWriter;
    use parquet::basic::Compression;
    use parquet::file::properties::WriterProperties;

    if entries.is_empty() {
        return Ok(());
    }

    let n = entries.len();
    let feature_names = get_pin_feature_names();

    // Build schema
    let mut fields: Vec<Field> = vec![
        Field::new("entry_id", DataType::UInt32, false),
        Field::new("is_decoy", DataType::Boolean, false),
        Field::new("sequence", DataType::Utf8, false),
        Field::new("modified_sequence", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("precursor_mz", DataType::Float64, false),
        Field::new("protein_ids", DataType::Utf8, true),
        Field::new("scan_number", DataType::UInt32, false),
        Field::new("apex_rt", DataType::Float64, false),
        Field::new("start_rt", DataType::Float64, false),
        Field::new("end_rt", DataType::Float64, false),
        Field::new("bounds_area", DataType::Float64, false),
        Field::new("bounds_snr", DataType::Float64, false),
        Field::new("file_name", DataType::Utf8, false),
        // Variable-length arrays as binary (LE bytes)
        Field::new("fragment_mzs", DataType::Binary, true),
        Field::new("fragment_intensities", DataType::Binary, true),
        Field::new("reference_xic_rts", DataType::Binary, true),
        Field::new("reference_xic_intensities", DataType::Binary, true),
    ];
    for name in &feature_names {
        fields.push(Field::new(*name, DataType::Float64, false));
    }
    let schema = std::sync::Arc::new(Schema::new(fields));

    // Build column arrays
    let mut entry_id_b = UInt32Builder::with_capacity(n);
    let mut decoy_b = BooleanBuilder::with_capacity(n);
    let mut seq_b = StringBuilder::with_capacity(n, n * 20);
    let mut modseq_b = StringBuilder::with_capacity(n, n * 30);
    let mut charge_b = UInt8Builder::with_capacity(n);
    let mut mz_b = Float64Builder::with_capacity(n);
    let mut protein_b = StringBuilder::with_capacity(n, n * 15);
    let mut scan_b = UInt32Builder::with_capacity(n);
    let mut apex_rt_b = Float64Builder::with_capacity(n);
    let mut start_rt_b = Float64Builder::with_capacity(n);
    let mut end_rt_b = Float64Builder::with_capacity(n);
    let mut area_b = Float64Builder::with_capacity(n);
    let mut snr_b = Float64Builder::with_capacity(n);
    let mut fname_b = StringBuilder::with_capacity(n, n * 20);
    let mut frag_mz_b = BinaryBuilder::with_capacity(n, n * 48);
    let mut frag_int_b = BinaryBuilder::with_capacity(n, n * 24);
    let mut xic_rt_b = BinaryBuilder::with_capacity(n, n * 80);
    let mut xic_int_b = BinaryBuilder::with_capacity(n, n * 80);
    let mut feat_builders: Vec<Float64Builder> = (0..feature_names.len())
        .map(|_| Float64Builder::with_capacity(n))
        .collect();

    for entry in entries {
        entry_id_b.append_value(entry.entry_id);
        decoy_b.append_value(entry.is_decoy);
        seq_b.append_value(&entry.sequence);
        modseq_b.append_value(&entry.modified_sequence);
        charge_b.append_value(entry.charge);
        mz_b.append_value(entry.precursor_mz);
        protein_b.append_value(entry.protein_ids.join(";"));
        scan_b.append_value(entry.scan_number);
        apex_rt_b.append_value(entry.apex_rt);
        start_rt_b.append_value(entry.peak_bounds.start_rt);
        end_rt_b.append_value(entry.peak_bounds.end_rt);
        area_b.append_value(entry.peak_bounds.area);
        snr_b.append_value(entry.peak_bounds.signal_to_noise);
        fname_b.append_value(&entry.file_name);

        // Serialize variable-length arrays as contiguous LE bytes
        let mz_bytes: Vec<u8> = entry
            .fragment_mzs
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        frag_mz_b.append_value(&mz_bytes);

        let int_bytes: Vec<u8> = entry
            .fragment_intensities
            .iter()
            .flat_map(|v| v.to_le_bytes())
            .collect();
        frag_int_b.append_value(&int_bytes);

        let xic_rts: Vec<u8> = entry
            .reference_xic
            .iter()
            .flat_map(|(rt, _)| rt.to_le_bytes())
            .collect();
        xic_rt_b.append_value(&xic_rts);

        let xic_ints: Vec<u8> = entry
            .reference_xic
            .iter()
            .flat_map(|(_, intensity)| intensity.to_le_bytes())
            .collect();
        xic_int_b.append_value(&xic_ints);

        for (i, builder) in feat_builders.iter_mut().enumerate() {
            let v = pin_feature_value(&entry.features, i);
            builder.append_value(if v.is_finite() { v } else { 0.0 });
        }
    }

    // Assemble columns
    let mut columns: Vec<ArrayRef> = vec![
        std::sync::Arc::new(entry_id_b.finish()),
        std::sync::Arc::new(decoy_b.finish()),
        std::sync::Arc::new(seq_b.finish()),
        std::sync::Arc::new(modseq_b.finish()),
        std::sync::Arc::new(charge_b.finish()),
        std::sync::Arc::new(mz_b.finish()),
        std::sync::Arc::new(protein_b.finish()),
        std::sync::Arc::new(scan_b.finish()),
        std::sync::Arc::new(apex_rt_b.finish()),
        std::sync::Arc::new(start_rt_b.finish()),
        std::sync::Arc::new(end_rt_b.finish()),
        std::sync::Arc::new(area_b.finish()),
        std::sync::Arc::new(snr_b.finish()),
        std::sync::Arc::new(fname_b.finish()),
        std::sync::Arc::new(frag_mz_b.finish()),
        std::sync::Arc::new(frag_int_b.finish()),
        std::sync::Arc::new(xic_rt_b.finish()),
        std::sync::Arc::new(xic_int_b.finish()),
    ];
    for builder in &mut feat_builders {
        columns.push(std::sync::Arc::new(builder.finish()));
    }

    // Write with ZSTD compression
    let file = std::fs::File::create(path).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to create scores file {}: {}",
            path.display(),
            e
        ))
    })?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| OspreyError::OutputError(format!("Failed to create RecordBatch: {}", e)))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| OspreyError::OutputError(format!("Failed to create Parquet writer: {}", e)))?;
    writer
        .write(&batch)
        .map_err(|e| OspreyError::OutputError(format!("Failed to write Parquet batch: {}", e)))?;
    writer
        .close()
        .map_err(|e| OspreyError::OutputError(format!("Failed to close Parquet writer: {}", e)))?;

    log::info!("Wrote {} scores to {}", n, path.display());

    Ok(())
}

/// Load per-file scored entries from a cached parquet file.
fn load_scores_parquet(path: &std::path::Path) -> Result<Vec<CoelutionScoredEntry>> {
    use arrow::array::{
        Array, BinaryArray, BooleanArray, Float64Array, StringArray, UInt32Array, UInt8Array,
    };
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to open scores file {}: {}",
            path.display(),
            e
        ))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
        OspreyError::OutputError(format!("Failed to read parquet {}: {}", path.display(), e))
    })?;
    let reader = builder.build().map_err(|e| {
        OspreyError::OutputError(format!("Failed to build reader {}: {}", path.display(), e))
    })?;

    let feature_names = get_pin_feature_names();
    let mut entries = Vec::new();

    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| OspreyError::OutputError(format!("Failed to read batch: {}", e)))?;

        let entry_ids = batch
            .column_by_name("entry_id")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let is_decoys = batch
            .column_by_name("is_decoy")
            .unwrap()
            .as_any()
            .downcast_ref::<BooleanArray>()
            .unwrap();
        let sequences = batch
            .column_by_name("sequence")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let mod_sequences = batch
            .column_by_name("modified_sequence")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let charges = batch
            .column_by_name("charge")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt8Array>()
            .unwrap();
        let mzs = batch
            .column_by_name("precursor_mz")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let proteins = batch
            .column_by_name("protein_ids")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let scans = batch
            .column_by_name("scan_number")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let apex_rts = batch
            .column_by_name("apex_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let start_rts = batch
            .column_by_name("start_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let end_rts = batch
            .column_by_name("end_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let areas = batch
            .column_by_name("bounds_area")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let snrs = batch
            .column_by_name("bounds_snr")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let file_names = batch
            .column_by_name("file_name")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let frag_mz_col = batch
            .column_by_name("fragment_mzs")
            .unwrap()
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap();
        let frag_int_col = batch
            .column_by_name("fragment_intensities")
            .unwrap()
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap();
        let xic_rt_col = batch
            .column_by_name("reference_xic_rts")
            .unwrap()
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap();
        let xic_int_col = batch
            .column_by_name("reference_xic_intensities")
            .unwrap()
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap();

        // Read feature columns
        let feat_cols: Vec<&Float64Array> = feature_names
            .iter()
            .map(|name| {
                batch
                    .column_by_name(name)
                    .unwrap()
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .unwrap()
            })
            .collect();

        for row in 0..batch.num_rows() {
            // Deserialize binary arrays back to Vecs
            let frag_mzs: Vec<f64> = frag_mz_col
                .value(row)
                .chunks_exact(8)
                .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            let frag_ints: Vec<f32> = frag_int_col
                .value(row)
                .chunks_exact(4)
                .map(|chunk| f32::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            let xic_rts_raw: Vec<f64> = xic_rt_col
                .value(row)
                .chunks_exact(8)
                .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            let xic_ints_raw: Vec<f64> = xic_int_col
                .value(row)
                .chunks_exact(8)
                .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            let reference_xic: Vec<(f64, f64)> =
                xic_rts_raw.into_iter().zip(xic_ints_raw).collect();

            let protein_str = proteins.value(row);
            let protein_ids: Vec<String> = if protein_str.is_empty() {
                Vec::new()
            } else {
                protein_str.split(';').map(|s| s.to_string()).collect()
            };

            // Reconstruct features from columns
            let mut features = CoelutionFeatureSet::default();
            set_features_from_pin_values(&mut features, &feat_cols, row);

            let apex_rt_val = apex_rts.value(row);

            entries.push(CoelutionScoredEntry {
                entry_id: entry_ids.value(row),
                is_decoy: is_decoys.value(row),
                sequence: sequences.value(row).to_string(),
                modified_sequence: mod_sequences.value(row).to_string(),
                charge: charges.value(row),
                precursor_mz: mzs.value(row),
                protein_ids,
                scan_number: scans.value(row),
                apex_rt: apex_rt_val,
                peak_bounds: XICPeakBounds {
                    apex_rt: apex_rt_val,
                    apex_intensity: 0.0,
                    apex_index: 0,
                    start_rt: start_rts.value(row),
                    end_rt: end_rts.value(row),
                    start_index: 0,
                    end_index: 0,
                    area: areas.value(row),
                    signal_to_noise: snrs.value(row),
                },
                features,
                fragment_mzs: frag_mzs,
                fragment_intensities: frag_ints,
                reference_xic,
                file_name: file_names.value(row).to_string(),
                run_qvalue: 1.0,
                experiment_qvalue: 1.0,
                score: 0.0,
                pep: 1.0,
            });
        }
    }

    Ok(entries)
}

/// Set CoelutionFeatureSet fields from parquet feature columns.
///
/// This maps the 45 PIN feature columns back to the struct fields.
fn set_features_from_pin_values(
    features: &mut CoelutionFeatureSet,
    feat_cols: &[&arrow::array::Float64Array],
    row: usize,
) {
    // The PIN feature order matches the order in get_pin_feature_names()
    // and pin_feature_value(). We reconstruct by index.
    let v = |i: usize| -> f64 { feat_cols[i].value(row) };

    // Pairwise coelution (indices 0-10)
    features.coelution_sum = v(0);
    features.coelution_min = v(1);
    features.coelution_max = v(2);
    features.n_coeluting_fragments = v(3) as u8;
    features.n_fragment_pairs = v(4) as u8;
    features.fragment_corr = [v(5), v(6), v(7), v(8), v(9), v(10)];

    // Peak shape (indices 11-17)
    features.peak_apex = v(11);
    features.peak_area = v(12);
    features.peak_width = v(13);
    features.peak_symmetry = v(14);
    features.signal_to_noise = v(15);
    features.n_scans = v(16) as u16;
    features.peak_sharpness = v(17);

    // Spectral at apex (indices 18-33)
    features.hyperscore = v(18);
    features.xcorr = v(19);
    features.dot_product = v(20);
    features.dot_product_smz = v(21);
    features.dot_product_top6 = v(22);
    features.dot_product_top5 = v(23);
    features.dot_product_top4 = v(24);
    features.dot_product_smz_top6 = v(25);
    features.dot_product_smz_top5 = v(26);
    features.dot_product_smz_top4 = v(27);
    features.fragment_coverage = v(28);
    features.sequence_coverage = v(29);
    features.consecutive_ions = v(30) as u8;
    features.explained_intensity = v(31);
    features.elution_weighted_cosine = v(32);

    // Mass accuracy (indices 33-35)
    features.mass_accuracy_mean = v(33);
    features.abs_mass_accuracy_mean = v(34);
    features.mass_accuracy_std = v(35);

    // RT deviation (indices 36-37)
    features.rt_deviation = v(36);
    features.abs_rt_deviation = v(37);

    // MS1 features (indices 38-39)
    features.ms1_precursor_coelution = v(38);
    features.ms1_isotope_cosine = v(39);

    // Peptide properties (indices 40-42)
    features.peptide_length = v(40) as u8;
    features.missed_cleavages = v(41) as u8;

    // Median polish (indices 42-44)
    features.median_polish_cosine = v(42);
    features.median_polish_rsquared = v(43);
    features.median_polish_residual_ratio = v(44);
}

/// Write coelution scored entries to Parquet report (one row per precursor per run)
fn write_parquet_report(path: &std::path::Path, entries: &[CoelutionScoredEntry]) -> Result<()> {
    use parquet::arrow::ArrowWriter;
    use parquet::basic::Compression;
    use parquet::file::properties::WriterProperties;

    if entries.is_empty() {
        return Ok(());
    }

    let n = entries.len();
    let feature_names = get_pin_feature_names();

    // Build schema: identification + peak/RT + FDR + all features
    let mut fields: Vec<Field> = vec![
        Field::new("Run", DataType::Utf8, false),
        Field::new("Modified.Sequence", DataType::Utf8, false),
        Field::new("Stripped.Sequence", DataType::Utf8, false),
        Field::new("Precursor.Charge", DataType::UInt8, false),
        Field::new("Precursor.Mz", DataType::Float64, false),
        Field::new("Protein.Ids", DataType::Utf8, true),
        Field::new("Is.Decoy", DataType::Boolean, false),
        Field::new("RT", DataType::Float64, false),
        Field::new("RT.Start", DataType::Float64, false),
        Field::new("RT.Stop", DataType::Float64, false),
        Field::new("Peak.Width", DataType::Float64, false),
        Field::new("Library.RT", DataType::Float64, false),
        Field::new("Scan.Number", DataType::UInt32, false),
        Field::new("Score", DataType::Float64, false),
        Field::new("Q.Value", DataType::Float64, false),
        Field::new("Global.Q.Value", DataType::Float64, false),
        Field::new("PEP", DataType::Float64, false),
        Field::new("Search.Mode", DataType::Utf8, false),
    ];
    for name in &feature_names {
        fields.push(Field::new(*name, DataType::Float64, false));
    }
    let schema = std::sync::Arc::new(Schema::new(fields));

    // Build column arrays
    let mut run_b = StringBuilder::with_capacity(n, n * 20);
    let mut modseq_b = StringBuilder::with_capacity(n, n * 30);
    let mut seq_b = StringBuilder::with_capacity(n, n * 20);
    let mut charge_b = UInt8Builder::with_capacity(n);
    let mut mz_b = Float64Builder::with_capacity(n);
    let mut protein_b = StringBuilder::with_capacity(n, n * 15);
    let mut decoy_b = BooleanBuilder::with_capacity(n);
    let mut rt_b = Float64Builder::with_capacity(n);
    let mut rt_start_b = Float64Builder::with_capacity(n);
    let mut rt_stop_b = Float64Builder::with_capacity(n);
    let mut width_b = Float64Builder::with_capacity(n);
    let mut lib_rt_b = Float64Builder::with_capacity(n);
    let mut scan_b = UInt32Builder::with_capacity(n);
    let mut score_b = Float64Builder::with_capacity(n);
    let mut qval_b = Float64Builder::with_capacity(n);
    let mut gqval_b = Float64Builder::with_capacity(n);
    let mut pep_b = Float64Builder::with_capacity(n);
    let mut mode_b = StringBuilder::with_capacity(n, n * 10);
    let mut feat_builders: Vec<Float64Builder> = (0..feature_names.len())
        .map(|_| Float64Builder::with_capacity(n))
        .collect();

    for entry in entries {
        run_b.append_value(&entry.file_name);
        modseq_b.append_value(&entry.modified_sequence);
        seq_b.append_value(&entry.sequence);
        charge_b.append_value(entry.charge);
        mz_b.append_value(entry.precursor_mz);
        protein_b.append_value(entry.protein_ids.join(";"));
        decoy_b.append_value(entry.is_decoy);

        rt_b.append_value(entry.peak_bounds.apex_rt);
        rt_start_b.append_value(entry.peak_bounds.start_rt);
        rt_stop_b.append_value(entry.peak_bounds.end_rt);
        width_b.append_value(entry.peak_bounds.end_rt - entry.peak_bounds.start_rt);
        lib_rt_b.append_value(0.0); // CoelutionScoredEntry doesn't store library RT
        scan_b.append_value(entry.scan_number);
        score_b.append_value(entry.score);
        qval_b.append_value(entry.run_qvalue);
        gqval_b.append_value(entry.experiment_qvalue);
        pep_b.append_value(entry.pep);
        mode_b.append_value("coelution");

        // All 45 features
        for (i, builder) in feat_builders.iter_mut().enumerate() {
            let v = pin_feature_value(&entry.features, i);
            builder.append_value(if v.is_finite() { v } else { 0.0 });
        }
    }

    // Assemble columns
    let mut columns: Vec<ArrayRef> = vec![
        std::sync::Arc::new(run_b.finish()),
        std::sync::Arc::new(modseq_b.finish()),
        std::sync::Arc::new(seq_b.finish()),
        std::sync::Arc::new(charge_b.finish()),
        std::sync::Arc::new(mz_b.finish()),
        std::sync::Arc::new(protein_b.finish()),
        std::sync::Arc::new(decoy_b.finish()),
        std::sync::Arc::new(rt_b.finish()),
        std::sync::Arc::new(rt_start_b.finish()),
        std::sync::Arc::new(rt_stop_b.finish()),
        std::sync::Arc::new(width_b.finish()),
        std::sync::Arc::new(lib_rt_b.finish()),
        std::sync::Arc::new(scan_b.finish()),
        std::sync::Arc::new(score_b.finish()),
        std::sync::Arc::new(qval_b.finish()),
        std::sync::Arc::new(gqval_b.finish()),
        std::sync::Arc::new(pep_b.finish()),
        std::sync::Arc::new(mode_b.finish()),
    ];
    for builder in &mut feat_builders {
        columns.push(std::sync::Arc::new(builder.finish()));
    }

    // Write parquet with ZSTD compression
    let file = std::fs::File::create(path).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to create parquet file {}: {}",
            path.display(),
            e
        ))
    })?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| OspreyError::OutputError(format!("Failed to create RecordBatch: {}", e)))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| OspreyError::OutputError(format!("Failed to create Parquet writer: {}", e)))?;
    writer
        .write(&batch)
        .map_err(|e| OspreyError::OutputError(format!("Failed to write Parquet batch: {}", e)))?;
    writer
        .close()
        .map_err(|e| OspreyError::OutputError(format!("Failed to close Parquet writer: {}", e)))?;

    log::info!(
        "Wrote Parquet report ({} entries, {} columns) to {}",
        n,
        18 + feature_names.len(),
        path.display()
    );
    Ok(())
}

/// Write coelution scored entries to TSV report
fn write_scored_report(path: &std::path::Path, entries: &[CoelutionScoredEntry]) -> Result<()> {
    use std::io::Write;

    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "modified_sequence\tprecursor_mz\tcharge\tapex_rt\tpeak_apex\tpeak_area\tn_scans\tq_value"
    )?;

    for entry in entries {
        writeln!(
            file,
            "{}\t{:.4}\t{}\t{:.2}\t{:.4}\t{:.4}\t{}\t{:.6}",
            entry.modified_sequence,
            entry.precursor_mz,
            entry.charge,
            entry.peak_bounds.apex_rt,
            entry.features.peak_apex,
            entry.features.peak_area,
            entry.features.n_scans,
            entry.experiment_qvalue
        )?;
    }

    Ok(())
}

// =============================================================================
// Main analysis entry point
// =============================================================================

/// Run the complete Osprey analysis pipeline.
///
/// Uses DIA-NN-style fragment XIC extraction and pairwise correlation scoring.
/// Calibrates retention time per file, searches for peptide candidates, performs
/// FDR control, and writes results to blib format for Skyline integration.
pub fn run_analysis(config: OspreyConfig) -> Result<()> {
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

    // Generate decoys
    if !config.decoys_in_library {
        log::info!("Generating decoys using {:?} method", config.decoy_method);
        let decoy_method = match config.decoy_method {
            CoreDecoyMethod::Reverse => DecoyMethod::Reverse,
            CoreDecoyMethod::Shuffle => DecoyMethod::Shuffle,
            CoreDecoyMethod::FromLibrary => DecoyMethod::Reverse,
        };
        let generator = DecoyGenerator::with_enzyme(decoy_method, Enzyme::Trypsin);
        let targets: Vec<LibraryEntry> = library.iter().filter(|e| !e.is_decoy).cloned().collect();
        let n_original_targets = targets.len();
        let (valid_targets, decoys, _stats) =
            generator.generate_all_with_collision_detection(&targets);
        log::info!(
            "Generated {} decoys from {} targets ({} excluded due to collisions)",
            decoys.len(),
            n_original_targets,
            n_original_targets - valid_targets.len()
        );
        library = valid_targets;
        library.extend(decoys);
        log::info!("Total library size: {} (targets + decoys)", library.len());
    }

    // Set up output directories
    let output_dir = config
        .output_blib
        .parent()
        .unwrap_or(std::path::Path::new("."));
    let write_pin = config.write_pin || config.fdr_method == FdrMethod::Mokapot;
    let mokapot_dir = output_dir.join("mokapot");
    if write_pin {
        std::fs::create_dir_all(&mokapot_dir)?;
    }
    let mokapot = MokapotRunner::new().with_test_fdr(config.run_fdr);

    // Process each file
    let mut per_file_entries: Vec<(String, Vec<CoelutionScoredEntry>)> = Vec::new();
    let mut pin_files: HashMap<String, std::path::PathBuf> = HashMap::new();

    for (file_idx, input_file) in config.input_files.iter().enumerate() {
        log::info!("");
        log::info!(
            "===== Processing file {}/{}: {} =====",
            file_idx + 1,
            config.input_files.len(),
            input_file.display()
        );

        let file_name = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Check for cached scores parquet (skip search if already scored)
        let scores_path = scores_path_for_input(input_file);
        let entries = if scores_path.exists() {
            match load_scores_parquet(&scores_path) {
                Ok(entries) => {
                    log::info!(
                        "Loaded {} cached scores from {}",
                        entries.len(),
                        scores_path.display()
                    );
                    entries
                }
                Err(e) => {
                    log::warn!(
                        "Failed to load cached scores from {}: {}. Re-scoring.",
                        scores_path.display(),
                        e
                    );
                    // Fall through to scoring below
                    Vec::new()
                }
            }
        } else {
            Vec::new()
        };

        // If we got cached entries, skip loading spectra and scoring
        let entries = if !entries.is_empty() {
            entries
        } else {
            // Load spectra
            let (spectra, ms1_index) = load_all_spectra(input_file)?;
            if spectra.is_empty() {
                log::warn!("No spectra found in {}", input_file.display());
                continue;
            }

            // Run calibration (reuse the windowed calibration discovery)
            let (rt_cal, cal_params): (Option<RTCalibration>, Option<CalibrationParams>) = if config
                .rt_calibration
                .enabled
            {
                // Try to load cached calibration
                let cached = input_file.parent().and_then(|input_dir| {
                    let cal_path = calibration_path_for_input(input_file, input_dir);
                    if !cal_path.exists() {
                        return None;
                    }
                    let cal_params = load_calibration(&cal_path).ok()?;
                    if !cal_params.rt_calibration.has_model_data() || !cal_params.is_calibrated() {
                        return None;
                    }
                    let model_params = cal_params.rt_calibration.model_params.as_ref()?;
                    let rt_cal = RTCalibration::from_model_params(
                        model_params,
                        cal_params.rt_calibration.residual_sd,
                    )
                    .ok()?;
                    log::info!("Reusing cached calibration from {}", cal_path.display());
                    cal_params.log_summary();
                    Some((rt_cal, cal_params))
                });

                if let Some((rt_cal, cal_params)) = cached {
                    (Some(rt_cal), Some(cal_params))
                } else {
                    // Run calibration discovery
                    match run_calibration_discovery_windowed(
                        &library, &spectra, &ms1_index, &config,
                    ) {
                        Ok((rt_cal, cal_params)) => {
                            cal_params.log_summary();
                            if let Some(input_dir) = input_file.parent() {
                                let cal_path = calibration_path_for_input(input_file, input_dir);
                                if let Err(e) = save_calibration(&cal_params, &cal_path) {
                                    log::warn!("Failed to save calibration: {}", e);
                                }
                            }
                            (Some(rt_cal), Some(cal_params))
                        }
                        Err(e) => {
                            log::warn!("Calibration failed: {}. Using fallback tolerance.", e);
                            (None, None)
                        }
                    }
                }
            } else {
                (None, None)
            };

            // Run coelution search
            let entries = run_search(
                &library,
                &spectra,
                &ms1_index,
                cal_params.as_ref(),
                rt_cal.as_ref(),
                &config,
                &file_name,
            )?;

            log::info!(
                "Scored {} entries ({} targets, {} decoys) for {}",
                entries.len(),
                entries.iter().filter(|e| !e.is_decoy).count(),
                entries.iter().filter(|e| e.is_decoy).count(),
                file_name
            );

            // Remove double-counted entries (different precursors sharing fragment ions)
            let entries = deduplicate_double_counting(
                entries,
                &library,
                &spectra,
                cal_params.as_ref(),
                &config,
            );

            // Deduplicate: keep best target and best decoy per base_id.
            let entries = deduplicate_pairs(entries);

            // Save scores to parquet for future re-use
            if let Err(e) = write_scores_parquet(&scores_path, &entries) {
                log::warn!("Failed to save scores to {}: {}", scores_path.display(), e);
            }

            entries
        };

        // Write coelution PIN file (if mokapot or --write-pin)
        // Mokapot expects pre-competed data (winners only), so compete before writing PIN.
        // Native Percolator does internal paired competition, so keep both target and decoy.
        if write_pin {
            let pin_entries = if config.fdr_method == FdrMethod::Mokapot {
                compete_target_decoy_pairs(entries.clone())
            } else {
                entries.clone()
            };
            let pin_path = mokapot.write_pin_file(&file_name, &pin_entries, &mokapot_dir)?;
            pin_files.insert(file_name.clone(), pin_path);
        }

        per_file_entries.push((file_name, entries));
    }

    let total_scored: usize = per_file_entries.iter().map(|(_, s)| s.len()).sum();
    log::info!(
        "Coelution analysis complete. {} total scored entries across {} files",
        total_scored,
        config.input_files.len()
    );

    if per_file_entries.is_empty() || total_scored == 0 {
        log::warn!("No scored entries found. Cannot perform FDR control.");
        return Ok(());
    }

    // Dispatch FDR control based on method
    log::info!(
        "Running {} FDR control on coelution results...",
        config.fdr_method
    );

    match config.fdr_method {
        FdrMethod::Percolator => {
            run_percolator_fdr(&mut per_file_entries, &config)?;
        }
        FdrMethod::Mokapot => {
            run_mokapot_fdr(
                &mut per_file_entries,
                &mokapot,
                &pin_files,
                &mokapot_dir,
                &config,
            )?;
        }
        FdrMethod::Simple => {
            for (_, entries) in per_file_entries.iter_mut() {
                apply_simple_fdr(entries, config.run_fdr)?;
            }
        }
    }

    // Collect all entries (targets + decoys) with q-values assigned
    let all_entries: Vec<CoelutionScoredEntry> = per_file_entries
        .into_iter()
        .flat_map(|(_, entries)| entries)
        .collect();

    // Determine which precursors pass experiment-level FDR (using best observation)
    let passing_precursors: HashSet<(String, u8)> = all_entries
        .iter()
        .filter(|e| !e.is_decoy && e.experiment_qvalue <= config.experiment_fdr)
        .map(|e| (e.modified_sequence.clone(), e.charge))
        .collect();

    // Include ALL per-file target observations for passing precursors
    // (not just the winner — needed for per-file RT boundaries in blib)
    let mut passing_entries: Vec<CoelutionScoredEntry> = all_entries
        .iter()
        .filter(|e| {
            !e.is_decoy && passing_precursors.contains(&(e.modified_sequence.clone(), e.charge))
        })
        .cloned()
        .collect();

    // Propagate the best experiment_qvalue to all observations of each precursor
    let mut best_exp_q: HashMap<(String, u8), f64> = HashMap::new();
    for e in &passing_entries {
        best_exp_q
            .entry((e.modified_sequence.clone(), e.charge))
            .and_modify(|q| *q = q.min(e.experiment_qvalue))
            .or_insert(e.experiment_qvalue);
    }
    for e in passing_entries.iter_mut() {
        if let Some(&q) = best_exp_q.get(&(e.modified_sequence.clone(), e.charge)) {
            e.experiment_qvalue = q;
        }
    }

    // Write blib output
    if !passing_entries.is_empty() {
        log::info!("Writing blib to {}", config.output_blib.display());
        write_blib_output(&config, &library, &passing_entries)?;
    } else {
        log::warn!("No peptides passed FDR threshold, skipping blib output");
    }

    // Write output report if specified
    // Parquet includes all entries (targets + decoys) for feature analysis and SVM retraining;
    // TSV includes only passing entries.
    if let Some(report_path) = &config.output_report {
        let ext = report_path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("");
        if ext.eq_ignore_ascii_case("parquet") {
            log::info!("Writing Parquet report to {}", report_path.display());
            write_parquet_report(report_path, &all_entries)?;
        } else {
            log::info!("Writing TSV report to {}", report_path.display());
            write_scored_report(report_path, &passing_entries)?;
        }
    }

    Ok(())
}

/// Run native Percolator FDR on coelution entries
fn run_percolator_fdr(
    per_file_entries: &mut [(String, Vec<CoelutionScoredEntry>)],
    config: &OspreyConfig,
) -> Result<()> {
    log::info!("Running native Percolator FDR on coelution entries");

    // Convert CoelutionScoredEntry → PercolatorEntry
    let mut perc_entries = Vec::new();
    for (file_name, entries) in per_file_entries.iter() {
        for entry in entries.iter() {
            let features: Vec<f64> = (0..NUM_PIN_FEATURES)
                .map(|i| {
                    let v = pin_feature_value(&entry.features, i);
                    if v.is_finite() {
                        v
                    } else {
                        0.0
                    }
                })
                .collect();

            let psm_id = format!(
                "{}_{}_{}_{}",
                file_name, entry.modified_sequence, entry.charge, entry.scan_number
            );

            perc_entries.push(percolator::PercolatorEntry {
                id: psm_id,
                file_name: file_name.clone(),
                peptide: entry.modified_sequence.clone(),
                charge: entry.charge,
                is_decoy: entry.is_decoy,
                entry_id: entry.entry_id,
                features,
            });
        }
    }

    let coelution_feature_names = get_pin_feature_names();
    let perc_config = percolator::PercolatorConfig {
        train_fdr: config.run_fdr,
        test_fdr: config.run_fdr,
        feature_names: Some(
            coelution_feature_names
                .iter()
                .map(|s| s.to_string())
                .collect(),
        ),
        ..Default::default()
    };

    let results = percolator::run_percolator(&perc_entries, &perc_config)
        .map_err(|e| OspreyError::config(format!("Percolator failed: {}", e)))?;

    // Build result lookup
    let result_map: HashMap<&str, &percolator::PercolatorResult> =
        results.entries.iter().map(|r| (r.id.as_str(), r)).collect();

    // Map results back to CoelutionScoredEntries
    for (file_name, entries) in per_file_entries.iter_mut() {
        for entry in entries.iter_mut() {
            let psm_id = format!(
                "{}_{}_{}_{}",
                file_name, entry.modified_sequence, entry.charge, entry.scan_number
            );
            if let Some(result) = result_map.get(psm_id.as_str()) {
                // Enforce dual FDR: max of precursor and peptide q-values
                entry.run_qvalue = result.run_precursor_qvalue.max(result.run_peptide_qvalue);
                entry.experiment_qvalue = result
                    .experiment_precursor_qvalue
                    .max(result.experiment_peptide_qvalue);
                entry.score = result.score;
                entry.pep = result.pep;
            }
        }
    }

    // Per-file and experiment stats already logged by percolator.rs

    Ok(())
}

/// Run external mokapot FDR on coelution entries
fn run_mokapot_fdr(
    per_file_entries: &mut [(String, Vec<CoelutionScoredEntry>)],
    mokapot: &MokapotRunner,
    pin_files: &HashMap<String, std::path::PathBuf>,
    mokapot_dir: &std::path::Path,
    config: &OspreyConfig,
) -> Result<()> {
    if !mokapot.is_available() {
        log::warn!("Mokapot not available, falling back to native Percolator");
        return run_percolator_fdr(per_file_entries, config);
    }

    match mokapot.run_two_step_analysis(pin_files, mokapot_dir) {
        Ok((per_file_results, experiment_results)) => {
            let experiment_qmap: HashMap<String, f64> = experiment_results
                .iter()
                .map(|r| (r.psm_id.clone(), r.q_value))
                .collect();

            let experiment_score_map: HashMap<String, (f64, f64)> = experiment_results
                .iter()
                .map(|r| (r.psm_id.clone(), (r.score, r.pep)))
                .collect();

            // Parse run-level peptide q-values for dual FDR control
            let step1_dir = mokapot_dir.join("run_level");
            let mut run_peptide_qmaps: HashMap<String, HashMap<String, f64>> = HashMap::new();
            for (file_name, _) in per_file_entries.iter() {
                if let Some(pin_path) = pin_files.get(file_name) {
                    let pin_stem = pin_path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("unknown");
                    let peptides_file =
                        step1_dir.join(format!("{}.mokapot.peptides.txt", pin_stem));
                    if peptides_file.exists() {
                        if let Ok(pept_results) = mokapot.parse_results(&peptides_file) {
                            let pept_map: HashMap<String, f64> = pept_results
                                .into_iter()
                                .map(|r| {
                                    let stripped =
                                        osprey_fdr::mokapot::strip_flanking_chars(&r.peptide);
                                    (stripped, r.q_value)
                                })
                                .collect();
                            run_peptide_qmaps.insert(file_name.clone(), pept_map);
                        }
                    }
                }
            }

            // Parse experiment-level peptide q-values for dual FDR control
            let step2_dir = mokapot_dir.join("experiment_level");
            let exp_peptide_qmap: HashMap<String, f64> = {
                let peptides_file = step2_dir.join("mokapot.peptides.txt");
                if peptides_file.exists() {
                    mokapot
                        .parse_results(&peptides_file)
                        .unwrap_or_default()
                        .into_iter()
                        .map(|r| {
                            let stripped = osprey_fdr::mokapot::strip_flanking_chars(&r.peptide);
                            (stripped, r.q_value)
                        })
                        .collect()
                } else {
                    HashMap::new()
                }
            };

            for (file_name, entries) in per_file_entries.iter_mut() {
                if let Some(results) = per_file_results.get(file_name) {
                    let run_qmap: HashMap<String, (f64, f64, f64)> = results
                        .iter()
                        .map(|r| (r.psm_id.clone(), (r.q_value, r.score, r.pep)))
                        .collect();

                    let peptide_qmap = run_peptide_qmaps.get(file_name);

                    for entry in entries.iter_mut() {
                        let psm_id = format!(
                            "{}_{}_{}_{}",
                            file_name, entry.modified_sequence, entry.charge, entry.scan_number
                        );
                        if let Some(&(q, score, pep)) = run_qmap.get(&psm_id) {
                            // Enforce dual FDR: max of precursor and peptide q-values
                            let peptide_q = peptide_qmap
                                .and_then(|m| m.get(entry.modified_sequence.as_str()))
                                .copied()
                                .unwrap_or(1.0);
                            entry.run_qvalue = q.max(peptide_q);
                            entry.score = score;
                            entry.pep = pep;
                        }
                        if let Some(&q) = experiment_qmap.get(&psm_id) {
                            // Enforce dual FDR: max of precursor and peptide q-values
                            let peptide_q = exp_peptide_qmap
                                .get(entry.modified_sequence.as_str())
                                .copied()
                                .unwrap_or(1.0);
                            entry.experiment_qvalue = q.max(peptide_q);
                        }
                        if let Some(&(score, pep)) = experiment_score_map.get(&psm_id) {
                            entry.score = score;
                            entry.pep = pep;
                        }
                    }
                }
            }
        }
        Err(e) => {
            log::warn!(
                "Mokapot two-step FDR failed: {}. Falling back to Percolator.",
                e
            );
            run_percolator_fdr(per_file_entries, config)?;
        }
    }

    Ok(())
}

/// Target-decoy competition for coelution entries.
/// Each target competes with its paired decoy — higher coelution_sum wins.
/// Deduplicate coelution entries: keep the best target AND best decoy per base_id.
///
/// When a precursor appears in multiple isolation windows, we may have duplicate
/// entries. This keeps only the best-scoring target and the best-scoring decoy
/// for each base_id (paired by `entry_id & 0x7FFFFFFF`). Both are retained so
/// that Percolator/Mokapot can perform proper target-decoy competition using
/// ML-rescored values.
fn deduplicate_pairs(entries: Vec<CoelutionScoredEntry>) -> Vec<CoelutionScoredEntry> {
    let mut target_best: HashMap<u32, CoelutionScoredEntry> = HashMap::new();
    let mut decoy_best: HashMap<u32, CoelutionScoredEntry> = HashMap::new();

    for entry in entries {
        let base_id = entry.entry_id & 0x7FFFFFFF;
        let map = if entry.is_decoy {
            &mut decoy_best
        } else {
            &mut target_best
        };

        let is_better = map
            .get(&base_id)
            .map(|existing| entry.features.coelution_sum > existing.features.coelution_sum)
            .unwrap_or(true);
        if is_better {
            map.insert(base_id, entry);
        }
    }

    let n_targets = target_best.len();
    let n_decoys = decoy_best.len();
    let n_paired = target_best
        .keys()
        .filter(|k| decoy_best.contains_key(k))
        .count();

    // Collect all entries (both targets and decoys)
    let mut all_entries: Vec<CoelutionScoredEntry> =
        Vec::with_capacity(target_best.len() + decoy_best.len());
    all_entries.extend(target_best.into_values());
    all_entries.extend(decoy_best.into_values());

    log::info!(
        "Deduplicated to {} targets + {} decoys ({} paired, {} total)",
        n_targets,
        n_decoys,
        n_paired,
        all_entries.len()
    );

    all_entries
}

/// Select the best-separating coelution feature by ROC AUC.
fn select_best_separating_feature(entries: &[CoelutionScoredEntry]) -> Option<(usize, bool)> {
    let mut best_index = 0usize;
    let mut best_abs_deviation = 0.0f64;
    let mut best_auc = 0.5f64;

    let feature_names = get_pin_feature_names();

    for feat_idx in 0..NUM_PIN_FEATURES {
        let mut target_vals = Vec::new();
        let mut decoy_vals = Vec::new();

        for entry in entries {
            let val = pin_feature_value(&entry.features, feat_idx);
            if !val.is_finite() {
                continue;
            }
            if entry.is_decoy {
                decoy_vals.push(val);
            } else {
                target_vals.push(val);
            }
        }

        if target_vals.len() < 2 || decoy_vals.len() < 2 {
            continue;
        }

        let auc = compute_roc_auc(&target_vals, &decoy_vals);
        let abs_deviation = (auc - 0.5).abs();
        if abs_deviation > best_abs_deviation {
            best_abs_deviation = abs_deviation;
            best_auc = auc;
            best_index = feat_idx;
        }
    }

    if best_abs_deviation < 1e-10 {
        log::warn!("Coelution target-decoy competition: no feature separates targets from decoys");
        return None;
    }

    let higher_is_better = best_auc > 0.5;
    log::info!(
        "Coelution target-decoy competition: best feature = {} (ROC AUC = {:.3}, {})",
        feature_names[best_index],
        best_auc,
        if higher_is_better {
            "higher is better"
        } else {
            "lower is better"
        }
    );

    Some((best_index, higher_is_better))
}

/// Run target-decoy competition on coelution entries for mokapot PIN writing.
///
/// Mokapot expects pre-competed data (only winners). This function competes
/// paired target-decoy entries using the best-separating coelution feature
/// (selected by ROC AUC). Ties go to the decoy (conservative for FDR).
fn compete_target_decoy_pairs(entries: Vec<CoelutionScoredEntry>) -> Vec<CoelutionScoredEntry> {
    let (feat_idx, higher_is_better) = match select_best_separating_feature(&entries) {
        Some(result) => result,
        None => {
            // Fallback: coelution_sum (index 0, higher is better)
            log::warn!(
                "Coelution target-decoy competition: using coelution_sum as fallback feature"
            );
            (0, true)
        }
    };

    // Group entries by base_id
    let mut groups: HashMap<u32, (Option<CoelutionScoredEntry>, Option<CoelutionScoredEntry>)> =
        HashMap::new();

    for entry in entries {
        let base_id = entry.entry_id & 0x7FFFFFFF;
        let slot = groups.entry(base_id).or_insert((None, None));
        if entry.is_decoy {
            slot.1 = Some(entry);
        } else {
            slot.0 = Some(entry);
        }
    }

    let mut winners = Vec::with_capacity(groups.len());
    let mut n_target_wins = 0usize;
    let mut n_decoy_wins = 0usize;

    for (_base_id, (target_opt, decoy_opt)) in groups {
        match (target_opt, decoy_opt) {
            (Some(target), None) => {
                n_target_wins += 1;
                winners.push(target);
            }
            (None, Some(decoy)) => {
                n_decoy_wins += 1;
                winners.push(decoy);
            }
            (Some(target), Some(decoy)) => {
                let t_val = pin_feature_value(&target.features, feat_idx);
                let d_val = pin_feature_value(&decoy.features, feat_idx);

                let t_val = if t_val.is_finite() {
                    t_val
                } else {
                    f64::NEG_INFINITY
                };
                let d_val = if d_val.is_finite() {
                    d_val
                } else {
                    f64::NEG_INFINITY
                };

                let target_wins = if higher_is_better {
                    t_val > d_val // strict: ties go to decoy
                } else {
                    t_val < d_val
                };

                if target_wins {
                    n_target_wins += 1;
                    winners.push(target);
                } else {
                    n_decoy_wins += 1;
                    winners.push(decoy);
                }
            }
            (None, None) => {}
        }
    }

    log::info!(
        "Coelution target-decoy competition: {} target wins, {} decoy wins ({} total)",
        n_target_wins,
        n_decoy_wins,
        winners.len()
    );

    winners
}

/// Remove double-counted coelution entries where different precursors share fragment ions.
///
/// Two entries are considered double-counted if:
/// 1. Their precursor m/z falls within the same DIA isolation window
/// 2. Their apex RTs are within ±5 spectra of each other
/// 3. ≥50% of their top 6 fragment ions match within calibrated m/z tolerance
///
/// The entry with the lower coelution_sum is removed.
fn deduplicate_double_counting(
    entries: Vec<CoelutionScoredEntry>,
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    calibration_params: Option<&CalibrationParams>,
    config: &OspreyConfig,
) -> Vec<CoelutionScoredEntry> {
    let original_count = entries.len();

    // 1. Extract unique isolation windows
    let scheme = extract_isolation_scheme(spectra);
    let windows = match &scheme {
        Some(s) => &s.windows,
        None => {
            log::warn!("Could not extract isolation scheme, skipping double-counting check");
            return entries;
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
        let mut sorted_rts: Vec<f64> = spectra.iter().map(|s| s.retention_time).collect();
        sorted_rts.sort_by(|a, b| a.total_cmp(b));
        sorted_rts.dedup();
        let mut intervals: Vec<f64> = Vec::new();
        for w in sorted_rts.windows(2) {
            intervals.push(w[1] - w[0]);
        }
        let median_interval = if intervals.is_empty() {
            0.05
        } else {
            intervals.sort_by(|a, b| a.total_cmp(b));
            intervals[intervals.len() / 2]
        };
        5.0 * median_interval
    };

    // 4. Build ID-to-index map for library lookup (IDs may not match array indices)
    let lib_id_map: HashMap<u32, usize> =
        library.iter().enumerate().map(|(i, e)| (e.id, i)).collect();

    // 5. Pre-assign entries to isolation windows using a sorted m/z index.
    //    This replaces the O(n * w) per-window full scan with O(n log n + w log n).
    let mut mz_sorted_indices: Vec<usize> = (0..entries.len()).collect();
    mz_sorted_indices.sort_by(|&a, &b| entries[a].precursor_mz.total_cmp(&entries[b].precursor_mz));
    let mz_sorted_vals: Vec<f64> = mz_sorted_indices
        .iter()
        .map(|&i| entries[i].precursor_mz)
        .collect();

    // For each window, binary search for the range of entries that fall within it
    let window_entry_indices: Vec<Vec<usize>> = windows
        .iter()
        .map(|&(center, width)| {
            let win_lower = center - width / 2.0;
            let win_upper = center + width / 2.0;
            let lo = mz_sorted_vals.partition_point(|&mz| mz < win_lower);
            let hi = mz_sorted_vals.partition_point(|&mz| mz <= win_upper);
            mz_sorted_indices[lo..hi].to_vec()
        })
        .collect();

    // 6. Process each isolation window in parallel.
    //    Within each window, sort entries by apex_rt for early termination.
    //    Use AtomicBool for thread-safe removal marking.
    let removed: Vec<AtomicBool> = (0..entries.len()).map(|_| AtomicBool::new(false)).collect();

    window_entry_indices.par_iter().for_each(|indices| {
        if indices.is_empty() {
            return;
        }

        // Sort by apex_rt for locality-based early termination
        let mut rt_sorted: Vec<usize> = indices
            .iter()
            .copied()
            .filter(|&i| !removed[i].load(Ordering::Relaxed))
            .collect();
        rt_sorted.sort_by(|&a, &b| entries[a].apex_rt.total_cmp(&entries[b].apex_rt));

        for i_pos in 0..rt_sorted.len() {
            let idx_a = rt_sorted[i_pos];
            if removed[idx_a].load(Ordering::Relaxed) {
                continue;
            }
            let apex_a = entries[idx_a].apex_rt;

            // Scan forward only; break when RT distance exceeds neighborhood
            for &idx_b in &rt_sorted[(i_pos + 1)..] {
                let apex_b = entries[idx_b].apex_rt;

                // Early termination: entries are RT-sorted, so all subsequent are farther
                if apex_b - apex_a > rt_neighborhood {
                    break;
                }

                if removed[idx_b].load(Ordering::Relaxed) {
                    continue;
                }

                // Check top-6 fragment overlap (look up library entries by ID)
                let lib_idx_a = match lib_id_map.get(&entries[idx_a].entry_id) {
                    Some(&idx) => idx,
                    None => continue,
                };
                let lib_idx_b = match lib_id_map.get(&entries[idx_b].entry_id) {
                    Some(&idx) => idx,
                    None => continue,
                };
                let entry_a = &library[lib_idx_a];
                let entry_b = &library[lib_idx_b];
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
                    // Remove the entry with lower coelution_sum (keep the better one)
                    if entries[idx_a].features.coelution_sum
                        >= entries[idx_b].features.coelution_sum
                    {
                        removed[idx_b].store(true, Ordering::Relaxed);
                    } else {
                        removed[idx_a].store(true, Ordering::Relaxed);
                        break; // idx_a is removed, no need to check more pairs
                    }
                }
            }
        }
    });

    let removed_count = removed.iter().filter(|r| r.load(Ordering::Relaxed)).count();
    if removed_count > 0 {
        log::info!(
            "Double-counting deduplication: removed {} entries ({} remaining)",
            removed_count,
            original_count - removed_count
        );
    }

    entries
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !removed[*i].load(Ordering::Relaxed))
        .map(|(_, e)| e)
        .collect()
}

/// Simple target-decoy FDR for coelution entries (fallback when Mokapot unavailable).
fn apply_simple_fdr(entries: &mut [CoelutionScoredEntry], _fdr_threshold: f64) -> Result<()> {
    // Sort by coelution_sum descending
    entries.sort_by(|a, b| {
        b.features
            .coelution_sum
            .total_cmp(&a.features.coelution_sum)
    });

    let mut n_targets = 0usize;
    let mut n_decoys = 0usize;

    for entry in entries.iter_mut() {
        if entry.is_decoy {
            n_decoys += 1;
        } else {
            n_targets += 1;
        }

        let fdr = if n_targets > 0 {
            n_decoys as f64 / n_targets as f64
        } else {
            1.0
        };

        entry.run_qvalue = fdr;
    }

    // Convert to q-values (monotonic)
    let mut min_fdr = 1.0f64;
    for entry in entries.iter_mut().rev() {
        min_fdr = min_fdr.min(entry.run_qvalue);
        entry.run_qvalue = min_fdr;
        entry.experiment_qvalue = min_fdr; // No two-level distinction in fallback mode
    }

    Ok(())
}

/// Build shared peak boundaries per (modified_sequence, file_name) for coelution path.
///
/// When the same peptide is detected at multiple charge states in the same file,
/// they elute at the same time. This function selects boundaries from the charge
/// state with the lowest run q-value (best overall peak quality) for each peptide
/// per file. Returns the lookup and a count of peptide-file pairs where boundaries
/// were shared across charge states.
#[allow(clippy::type_complexity)]
fn build_shared_boundaries(
    entries: &[CoelutionScoredEntry],
) -> (HashMap<(String, String), (f64, f64, f64)>, usize) {
    // For each (modified_sequence, file_name), keep boundaries from the entry
    // with the lowest run q-value
    let mut best_by_key: HashMap<(String, String), (f64, f64, f64, f64)> = HashMap::new();
    // Track which keys have multiple charge states
    let mut charges_by_key: HashMap<(String, String), HashSet<u8>> = HashMap::new();

    for entry in entries {
        let key = (entry.modified_sequence.clone(), entry.file_name.clone());

        charges_by_key
            .entry(key.clone())
            .or_default()
            .insert(entry.charge);

        let should_update = best_by_key
            .get(&key)
            .map_or(true, |&(_, _, _, qval)| entry.run_qvalue < qval);

        if should_update {
            best_by_key.insert(
                key,
                (
                    entry.apex_rt,
                    entry.peak_bounds.start_rt,
                    entry.peak_bounds.end_rt,
                    entry.run_qvalue,
                ),
            );
        }
    }

    let n_shared = charges_by_key.values().filter(|c| c.len() > 1).count();

    let shared_bounds: HashMap<(String, String), (f64, f64, f64)> = best_by_key
        .into_iter()
        .map(|(k, (apex, start, end, _))| (k, (apex, start, end)))
        .collect();

    (shared_bounds, n_shared)
}

/// Write coelution-mode blib output for Skyline.
///
/// Implementation details:
/// - Groups entries by precursor (modified_sequence + charge)
/// - Uses the best run (lowest run q-value) for the RefSpectra entry
/// - Writes per-run RetentionTimes entries with run-level q-values
/// - Populates OspreyRunScores and OspreyExperimentScores tables
fn write_blib_output(
    config: &OspreyConfig,
    library: &[LibraryEntry],
    entries: &[CoelutionScoredEntry],
) -> Result<()> {
    let mut writer = BlibWriter::create(&config.output_blib)?;

    writer.add_metadata("osprey_version", env!("CARGO_PKG_VERSION"))?;
    writer.add_metadata("search_mode", "coelution")?;
    writer.add_metadata("run_fdr", &config.run_fdr.to_string())?;
    writer.add_metadata("experiment_fdr", &config.experiment_fdr.to_string())?;

    let library_name = config
        .library_source
        .path()
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("library")
        .to_string();

    // Add source files — build lookup by file stem for matching with entry.file_name
    let blib_dir = config.output_blib.parent();
    let mut file_stem_to_id: HashMap<String, i64> = HashMap::new();
    for input_file in &config.input_files {
        let file_name = input_file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
        let file_stem = input_file
            .file_stem()
            .and_then(|n| n.to_str())
            .unwrap_or("unknown");
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

    writer.begin_batch()?;

    // Build shared boundaries: for each peptide per file, use boundaries from
    // the charge state with the lowest run q-value (same peptide at different
    // charge states elutes at the same time)
    let (shared_bounds, n_shared) = build_shared_boundaries(entries);

    // Build library ID lookup for O(1) access instead of O(n) linear scan
    let lib_by_id: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    // Group entries by precursor (modified_sequence + charge) so we write one RefSpectra
    // per unique precursor, with per-run entries in RetentionTimes
    let mut precursor_groups: HashMap<(String, u8), Vec<&CoelutionScoredEntry>> = HashMap::new();
    for entry in entries {
        precursor_groups
            .entry((entry.modified_sequence.clone(), entry.charge))
            .or_default()
            .push(entry);
    }

    let mut n_written = 0usize;

    for ((_mod_seq, _charge), group) in &precursor_groups {
        // Find the best run (lowest run q-value) for this precursor
        let best = group
            .iter()
            .min_by(|a, b| a.run_qvalue.total_cmp(&b.run_qvalue))
            .unwrap();

        if best.fragment_mzs.is_empty() {
            continue;
        }

        let n_runs_detected = group.len() as i32;
        let file_id = *file_stem_to_id.get(&best.file_name).unwrap_or(&1);

        let tic: f64 = best.fragment_intensities.iter().map(|&x| x as f64).sum();

        // Use shared boundaries from the best charge state for this peptide in this file
        let shared_key = (best.modified_sequence.clone(), best.file_name.clone());
        let (shared_apex, shared_start, shared_end) =
            shared_bounds.get(&shared_key).copied().unwrap_or((
                best.apex_rt,
                best.peak_bounds.start_rt,
                best.peak_bounds.end_rt,
            ));

        // Score is the experiment q-value — Skyline GENERIC Q-VALUE convention
        let ref_id = writer.add_spectrum(
            &best.sequence,
            &best.modified_sequence,
            best.precursor_mz,
            best.charge as i32,
            shared_apex,
            shared_start,
            shared_end,
            &best.fragment_mzs,
            &best.fragment_intensities,
            best.experiment_qvalue,
            file_id,
            n_runs_detected,
            tic,
        )?;

        // Add modifications from library entry if present
        if let Some(lib_entry) = lib_by_id.get(&best.entry_id) {
            if !lib_entry.modifications.is_empty() {
                writer.add_modifications(ref_id, &lib_entry.modifications)?;
            }
            if !lib_entry.protein_ids.is_empty() {
                writer.add_protein_mapping(ref_id, &lib_entry.protein_ids)?;
            }
        }

        // Write RetentionTimes entries for every run where this precursor was detected
        for scored in group {
            let run_file_id = *file_stem_to_id.get(&scored.file_name).unwrap_or(&1);
            let is_best = std::ptr::eq(
                *scored as *const CoelutionScoredEntry,
                *best as *const CoelutionScoredEntry,
            );

            // Use shared boundaries for this peptide in this run's file
            let run_shared_key = (scored.modified_sequence.clone(), scored.file_name.clone());
            let (run_apex, run_start, run_end) =
                shared_bounds.get(&run_shared_key).copied().unwrap_or((
                    scored.apex_rt,
                    scored.peak_bounds.start_rt,
                    scored.peak_bounds.end_rt,
                ));

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

        // Add peak boundaries to Osprey extension table (shared across charge states)
        let boundaries = osprey_core::PeakBoundaries {
            start_rt: shared_start,
            end_rt: shared_end,
            apex_rt: shared_apex,
            apex_coefficient: best.features.dot_product,
            integrated_area: best.features.peak_area,
            peak_quality: osprey_core::PeakQuality::default(),
        };
        writer.add_peak_boundaries(ref_id, &best.file_name, &boundaries)?;

        // Add run-level scores for best run
        writer.add_run_scores(
            ref_id,
            &best.file_name,
            best.run_qvalue,
            best.features.dot_product,
            0.0, // No PEP from coelution scoring yet
        )?;

        // Add experiment-level scores
        writer.add_experiment_scores(
            ref_id,
            best.experiment_qvalue,
            n_runs_detected,
            config.input_files.len() as i32,
        )?;

        n_written += 1;
    }

    writer.commit()?;
    writer.finalize()?;

    log::info!(
        "Wrote {} total entries for {} precursors across {} files",
        entries.len(),
        n_written,
        config.input_files.len()
    );
    if n_shared > 0 {
        log::info!(
            "Shared peak boundaries across charge states for {} peptide-file pairs",
            n_shared
        );
    }

    Ok(())
}

// =============================================================================
// Coelution-based search (DIA-NN style)
// =============================================================================

/// Run coelution-based search for a single file.
///
/// For each precursor candidate in each isolation window:
/// 1. Extract top-6 fragment XICs
/// 2. Detect peak on reference XIC (DIA-NN-style valley detection)
/// 3. Compute pairwise fragment correlations within peak boundaries
/// 4. Score spectrum at apex (spectral similarity, mass accuracy)
/// 5. Compute peak shape, RT deviation, MS1, and median polish features
///
/// Returns `CoelutionScoredEntry` for each detected precursor (targets + decoys).
#[allow(clippy::too_many_arguments)]
fn run_search(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    ms1_index: &MS1Index,
    calibration: Option<&CalibrationParams>,
    rt_calibration: Option<&RTCalibration>,
    config: &OspreyConfig,
    file_name: &str,
) -> Result<Vec<CoelutionScoredEntry>> {
    use osprey_chromatography::calibration::calibrated_tolerance_ppm;
    use osprey_chromatography::{
        compute_snr, detect_all_xic_peaks, isolate_peak_region, EmgFitter,
    };
    use osprey_scoring::batch::{group_spectra_by_isolation_window, MIN_COELUTION_SPECTRA};
    use osprey_scoring::{
        compute_elution_weighted_cosine, compute_mass_accuracy, extract_fragment_xics,
        median_polish_libcosine, median_polish_residual_ratio, median_polish_rsquared,
        pearson_correlation_raw, tukey_median_polish, SpectralScorer,
    };

    let is_hram = matches!(config.resolution_mode, osprey_core::ResolutionMode::HRAM);

    // Set up fragment tolerance (calibrated if available)
    let fragment_tolerance = if let Some(cal) = calibration {
        let (tol_val, tol_unit) = osprey_chromatography::calibrated_tolerance(
            &cal.ms2_calibration,
            config.fragment_tolerance.tolerance,
            config.fragment_tolerance.unit,
        );
        log::info!(
            "Coelution search using calibrated fragment tolerance: {:.4} {}",
            tol_val,
            match tol_unit {
                ToleranceUnit::Ppm => "ppm",
                ToleranceUnit::Mz => "Th",
            }
        );
        FragmentToleranceConfig {
            tolerance: tol_val,
            unit: tol_unit,
        }
    } else {
        config.fragment_tolerance
    };

    let tol_da = match fragment_tolerance.unit {
        ToleranceUnit::Mz => fragment_tolerance.tolerance,
        ToleranceUnit::Ppm => 0.0,
    };
    let tol_ppm = match fragment_tolerance.unit {
        ToleranceUnit::Ppm => fragment_tolerance.tolerance,
        ToleranceUnit::Mz => 0.0,
    };

    // Apply MS2 calibration offset to spectra if available
    // This shifts each centroid by the mean offset so the error distribution is centered at 0
    let calibrated_spectra: Option<Vec<Spectrum>> = calibration.and_then(|cal| {
        if cal.ms2_calibration.calibrated {
            log::info!(
                "Applying MS2 calibration: mean error = {:.4} {} → correcting by {:+.4} {} (using 3×SD = {:.3} {} tolerance)",
                cal.ms2_calibration.mean,
                cal.ms2_calibration.unit,
                -cal.ms2_calibration.mean,
                cal.ms2_calibration.unit,
                3.0 * cal.ms2_calibration.sd,
                cal.ms2_calibration.unit
            );
            Some(
                spectra
                    .iter()
                    .map(|s| {
                        osprey_chromatography::apply_spectrum_calibration(s, &cal.ms2_calibration)
                    })
                    .collect(),
            )
        } else {
            None
        }
    });
    let spectra_ref = calibrated_spectra.as_deref().unwrap_or(spectra);

    // RT tolerance for candidate selection
    let rt_tolerance = if let Some(cal) = calibration {
        let (raw_tol, method_desc) = if config.rt_calibration.use_percentile_tolerance {
            // DIA-NN-style: max(2 × P20_abs_error, RT_range / 40)
            if let Some(p20) = cal.rt_calibration.p20_abs_residual {
                let rt_range = rt_calibration
                    .as_ref()
                    .map(|c| {
                        let (rmin, rmax) = c.measured_rt_range();
                        rmax - rmin
                    })
                    .unwrap_or_else(|| {
                        spectra
                            .iter()
                            .map(|s| s.retention_time)
                            .fold(f64::NEG_INFINITY, f64::max)
                            - spectra
                                .iter()
                                .map(|s| s.retention_time)
                                .fold(f64::INFINITY, f64::min)
                    });
                let from_percentile = 2.0 * p20;
                let from_range = rt_range / 40.0;
                let tol = from_percentile.max(from_range);
                (
                    tol,
                    format!(
                        "percentile: max(2×P20={:.2}, range/40={:.2})",
                        from_percentile, from_range
                    ),
                )
            } else {
                // Fallback for old calibration files without p20
                let tol =
                    cal.rt_calibration.residual_sd * config.rt_calibration.rt_tolerance_factor;
                (
                    tol,
                    format!(
                        "{}× residual SD {:.2} min (p20 unavailable)",
                        config.rt_calibration.rt_tolerance_factor, cal.rt_calibration.residual_sd
                    ),
                )
            }
        } else {
            // Classic SD-based tolerance
            let tol = cal.rt_calibration.residual_sd * config.rt_calibration.rt_tolerance_factor;
            (
                tol,
                format!(
                    "{}× residual SD {:.2} min",
                    config.rt_calibration.rt_tolerance_factor, cal.rt_calibration.residual_sd
                ),
            )
        };

        let tol = raw_tol.max(config.rt_calibration.min_rt_tolerance);
        let tol = tol.min(config.rt_calibration.max_rt_tolerance);

        if raw_tol > config.rt_calibration.max_rt_tolerance {
            log::warn!(
                "RT tolerance capped at {:.2} min (uncapped: {:.2} min, {}). \
                 Poor RT calibration may reduce search sensitivity.",
                tol,
                raw_tol,
                method_desc
            );
        } else {
            log::info!(
                "Using calibrated RT tolerance: {:.2} min ({})",
                tol,
                method_desc
            );
        }
        tol
    } else {
        log::info!(
            "Using fallback RT tolerance: {:.2} min",
            config.rt_calibration.fallback_rt_tolerance
        );
        config.rt_calibration.fallback_rt_tolerance
    };

    // Build MzRT index for candidate selection
    let mz_rt_index = MzRTIndex::build(library, rt_calibration);

    // Group spectra by isolation window
    let window_groups = group_spectra_by_isolation_window(spectra_ref);
    if window_groups.is_empty() {
        log::warn!("No spectra with isolation windows found");
        return Ok(Vec::new());
    }

    log::info!(
        "Coelution search: {} spectra in {} windows, {} library entries",
        spectra_ref.len(),
        window_groups.len(),
        library.len()
    );

    // Set up spectral scorer with appropriate XCorr binning
    let scorer = if is_hram {
        log::debug!(
            "HRAM XCorr bins (0.02 Th), fragment tolerance: {:.2} ppm",
            fragment_tolerance.tolerance
        );
        SpectralScorer::hram().with_tolerance_ppm(fragment_tolerance.tolerance)
    } else {
        log::debug!(
            "Unit resolution XCorr bins (1.0 Th), fragment tolerance: {:.4} Th",
            fragment_tolerance.tolerance
        );
        SpectralScorer::new().with_tolerance_da(fragment_tolerance.tolerance)
    };

    // Progress bar
    let total_candidates = library.len() as u64;
    let pb = ProgressBar::new(total_candidates);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} precursors ({per_sec})",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Process each isolation window — candidates within each window in parallel
    let all_entries: Vec<CoelutionScoredEntry> = window_groups
        .par_iter()
        .flat_map(|((lower, upper), spectrum_indices)| {
            // Gather spectra for this window, sorted by RT
            let mut window_pairs: Vec<(usize, &Spectrum)> = spectrum_indices
                .iter()
                .map(|&idx| (idx, &spectra_ref[idx]))
                .collect();
            window_pairs.sort_by(|a, b| {
                a.1.retention_time
                    .partial_cmp(&b.1.retention_time)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            if window_pairs.len() < MIN_COELUTION_SPECTRA {
                return Vec::new();
            }

            let window_spectra: Vec<&Spectrum> = window_pairs.iter().map(|(_, s)| *s).collect();
            let _global_indices: Vec<usize> = window_pairs.iter().map(|(idx, _)| *idx).collect();

            // Find candidate library entries for this window using MzRTIndex
            // Collect unique library indices whose precursor_mz falls in this window
            let _window_center = (*lower + *upper) / 2.0;
            let lower_bin = lower.floor() as i32;
            let upper_bin = upper.ceil() as i32;

            let mut candidate_indices: Vec<usize> = Vec::new();
            let mut seen = std::collections::HashSet::new();

            // Get median RT of spectra in this window for RT range filtering
            let window_rt_min = window_spectra
                .first()
                .map(|s| s.retention_time)
                .unwrap_or(0.0);
            let window_rt_max = window_spectra
                .last()
                .map(|s| s.retention_time)
                .unwrap_or(0.0);

            for bin in lower_bin..=upper_bin {
                if let Some(entries) = mz_rt_index.get(bin) {
                    for &(expected_rt, idx) in entries {
                        if seen.contains(&idx) {
                            continue;
                        }
                        // Check m/z falls within isolation window
                        let entry = &library[idx];
                        if entry.precursor_mz < *lower || entry.precursor_mz > *upper {
                            continue;
                        }
                        // Check RT is within a reasonable range of window spectra
                        if expected_rt < window_rt_min - rt_tolerance
                            || expected_rt > window_rt_max + rt_tolerance
                        {
                            continue;
                        }
                        seen.insert(idx);
                        candidate_indices.push(idx);
                    }
                }
            }

            // Score each candidate precursor
            let window_entries: Vec<CoelutionScoredEntry> = candidate_indices
                .iter()
                .filter_map(|&lib_idx| {
                    let entry = &library[lib_idx];

                    // Expected RT for this entry
                    let expected_rt = rt_calibration
                        .map(|cal| cal.predict(entry.retention_time))
                        .unwrap_or(entry.retention_time);

                    // Filter spectra within RT tolerance of expected RT
                    let rt_spec_pairs: Vec<(usize, &Spectrum)> = window_pairs
                        .iter()
                        .filter(|(_, spec)| {
                            (spec.retention_time - expected_rt).abs() <= rt_tolerance
                        })
                        .map(|(idx, spec)| (*idx, *spec))
                        .collect();

                    if rt_spec_pairs.len() < MIN_COELUTION_SPECTRA {
                        pb.inc(1);
                        return None;
                    }

                    let cand_spectra: Vec<&Spectrum> =
                        rt_spec_pairs.iter().map(|(_, s)| *s).collect();
                    let cand_global: Vec<usize> =
                        rt_spec_pairs.iter().map(|(idx, _)| *idx).collect();

                    // 1. Extract top-6 fragment XICs
                    let xics =
                        extract_fragment_xics(&entry.fragments, &cand_spectra, tol_da, tol_ppm, 6);

                    if xics.len() < 2 {
                        pb.inc(1);
                        return None;
                    }

                    // 2. Find reference XIC (highest total intensity) for scoring
                    let ref_idx = xics
                        .iter()
                        .enumerate()
                        .max_by(|a, b| {
                            let sum_a: f64 = a.1 .1.iter().map(|(_, v)| *v).sum();
                            let sum_b: f64 = b.1 .1.iter().map(|(_, v)| *v).sum();
                            sum_a.total_cmp(&sum_b)
                        })
                        .map(|(i, _)| i)?;

                    let ref_xic = &xics[ref_idx].1;

                    // 3. Detect candidate peaks and pick the best by coelution score.
                    // Find all candidate peaks from median polish elution profile (or
                    // reference XIC as fallback), then evaluate each candidate by computing
                    // pairwise fragment correlations within its boundaries. The candidate
                    // with the highest coelution score is selected — the real peak is
                    // where fragments co-elute.
                    let full_polish = tukey_median_polish(&xics, 10, 0.01);
                    let candidates = {
                        let mp_candidates = full_polish
                            .as_ref()
                            .map(|mp| detect_all_xic_peaks(&mp.elution_profile, 0.01, 5.0))
                            .unwrap_or_default();
                        if mp_candidates.is_empty() {
                            detect_all_xic_peaks(ref_xic, 0.01, 5.0)
                        } else {
                            mp_candidates
                        }
                    };

                    if candidates.is_empty() {
                        pb.inc(1);
                        return None;
                    }

                    // Score each candidate by mean pairwise fragment correlation
                    let peak = if candidates.len() == 1 {
                        let bp = &candidates[0];
                        let si = bp.start_index;
                        let ei = bp.end_index;
                        let (apex_idx, apex_val) = ref_xic[si..=ei]
                            .iter()
                            .enumerate()
                            .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                            .map(|(i, &(_, v))| (si + i, v))
                            .unwrap_or((bp.apex_index, 0.0));
                        let area = trapezoidal_area(&ref_xic[si..=ei]);
                        let raw_ints: Vec<f64> = ref_xic.iter().map(|(_, v)| *v).collect();
                        let snr = compute_snr(&raw_ints, apex_idx, si, ei);
                        XICPeakBounds {
                            apex_rt: ref_xic[apex_idx].0,
                            apex_intensity: apex_val,
                            apex_index: apex_idx,
                            start_rt: ref_xic[si].0,
                            end_rt: ref_xic[ei].0,
                            start_index: si,
                            end_index: ei,
                            area,
                            signal_to_noise: snr,
                        }
                    } else {
                        // Evaluate each candidate by coelution score
                        let best = candidates
                            .iter()
                            .map(|bp| {
                                let si = bp.start_index;
                                let ei = bp.end_index;
                                let peak_len = ei - si + 1;

                                let coelution_score = if peak_len >= 3 {
                                    let vals: Vec<Vec<f64>> = xics
                                        .iter()
                                        .map(|(_, xic_data)| {
                                            xic_data[si..=ei].iter().map(|(_, v)| *v).collect()
                                        })
                                        .collect();
                                    let mut sum = 0.0f64;
                                    let mut count = 0u32;
                                    for ii in 0..vals.len() {
                                        for jj in (ii + 1)..vals.len() {
                                            sum += pearson_correlation_raw(&vals[ii], &vals[jj]);
                                            count += 1;
                                        }
                                    }
                                    if count > 0 {
                                        sum / count as f64
                                    } else {
                                        0.0
                                    }
                                } else {
                                    0.0
                                };

                                (bp, coelution_score)
                            })
                            .max_by(|a, b| a.1.total_cmp(&b.1))
                            .unwrap();

                        let bp = best.0;
                        let si = bp.start_index;
                        let ei = bp.end_index;
                        let (apex_idx, apex_val) = ref_xic[si..=ei]
                            .iter()
                            .enumerate()
                            .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                            .map(|(i, &(_, v))| (si + i, v))
                            .unwrap_or((bp.apex_index, 0.0));
                        let area = trapezoidal_area(&ref_xic[si..=ei]);
                        let raw_ints: Vec<f64> = ref_xic.iter().map(|(_, v)| *v).collect();
                        let snr = compute_snr(&raw_ints, apex_idx, si, ei);
                        XICPeakBounds {
                            apex_rt: ref_xic[apex_idx].0,
                            apex_intensity: apex_val,
                            apex_index: apex_idx,
                            start_rt: ref_xic[si].0,
                            end_rt: ref_xic[ei].0,
                            start_index: si,
                            end_index: ei,
                            area,
                            signal_to_noise: snr,
                        }
                    };

                    // Refine boundaries via EMG fit to median polish consensus profile.
                    // The EMG (Exponentially Modified Gaussian) captures LC peak asymmetry
                    // and gives 95% CDF boundaries that are much tighter than intensity-based
                    // valley detection. Falls back to original boundaries if fit fails.
                    //
                    // CRITICAL: crop the profile to the single peak around the apex before
                    // fitting. The full profile may contain multiple peaks across the
                    // extraction window; the EMG is a single-peak model and would produce
                    // overly broad boundaries if fit to the entire multi-peak profile.
                    let peak = {
                        let full_profile = full_polish
                            .as_ref()
                            .map(|mp| &mp.elution_profile[..])
                            .unwrap_or(ref_xic);
                        let (crop_start, crop_end) =
                            isolate_peak_region(full_profile, peak.apex_index);
                        let emg_profile = &full_profile[crop_start..=crop_end];
                        let fitter = EmgFitter::new();
                        match fitter.fit(emg_profile) {
                            Ok(params) if params.r_squared > 0.5 => {
                                if let Some((left_rt, right_rt)) =
                                    EmgFitter::boundaries(&params, 0.99)
                                {
                                    let apex_rt = ref_xic[peak.apex_index].0;
                                    // Validate: EMG boundaries must contain the apex
                                    if left_rt < apex_rt && right_rt > apex_rt {
                                        // Map EMG CDF boundaries to XIC indices.
                                        // Don't constrain to original peak boundaries —
                                        // let the EMG determine its own boundaries.
                                        let n_xic = ref_xic.len();
                                        let si = ref_xic
                                            .partition_point(|(rt, _)| *rt < left_rt)
                                            .min(peak.apex_index); // safety: si <= apex
                                        let ei = ref_xic
                                            .partition_point(|(rt, _)| *rt <= right_rt)
                                            .saturating_sub(1)
                                            .min(n_xic - 1)
                                            .max(peak.apex_index); // safety: ei >= apex
                                        if ei - si >= 2 {
                                            let raw_ints: Vec<f64> =
                                                ref_xic.iter().map(|(_, v)| *v).collect();
                                            let snr =
                                                compute_snr(&raw_ints, peak.apex_index, si, ei);
                                            XICPeakBounds {
                                                start_index: si,
                                                end_index: ei,
                                                start_rt: ref_xic[si].0,
                                                end_rt: ref_xic[ei].0,
                                                area: trapezoidal_area(&ref_xic[si..=ei]),
                                                signal_to_noise: snr,
                                                ..peak
                                            }
                                        } else {
                                            peak // too narrow
                                        }
                                    } else {
                                        peak // EMG boundaries don't contain apex
                                    }
                                } else {
                                    peak
                                }
                            }
                            _ => peak,
                        }
                    };

                    if peak.end_index <= peak.start_index + 1 {
                        pb.inc(1);
                        return None;
                    }

                    // 4. Compute pairwise fragment correlations within peak boundaries
                    let n_xics = xics.len();
                    let peak_len = peak.end_index - peak.start_index + 1;

                    // Extract XIC values within peak boundaries
                    let xic_peak_values: Vec<Vec<f64>> = xics
                        .iter()
                        .map(|(_, xic)| {
                            xic[peak.start_index..=peak.end_index]
                                .iter()
                                .map(|(_, v)| *v)
                                .collect()
                        })
                        .collect();

                    // Pairwise correlations
                    let mut corr_sum = 0.0f64;
                    let mut corr_min = f64::INFINITY;
                    let mut corr_max = f64::NEG_INFINITY;
                    let mut n_pairs = 0u32;
                    let mut per_frag_corr_sum = vec![0.0f64; n_xics];
                    let mut per_frag_corr_count = vec![0u32; n_xics];

                    for i in 0..n_xics {
                        for j in (i + 1)..n_xics {
                            let r =
                                pearson_correlation_raw(&xic_peak_values[i], &xic_peak_values[j]);
                            corr_sum += r;
                            if r < corr_min {
                                corr_min = r;
                            }
                            if r > corr_max {
                                corr_max = r;
                            }
                            n_pairs += 1;
                            per_frag_corr_sum[i] += r;
                            per_frag_corr_count[i] += 1;
                            per_frag_corr_sum[j] += r;
                            per_frag_corr_count[j] += 1;
                        }
                    }

                    if n_pairs == 0 {
                        pb.inc(1);
                        return None;
                    }

                    let mut fragment_corr = [0.0f64; 6];
                    let mut n_coeluting = 0u8;
                    for i in 0..n_xics.min(6) {
                        if per_frag_corr_count[i] > 0 {
                            fragment_corr[i] = per_frag_corr_sum[i] / per_frag_corr_count[i] as f64;
                            if fragment_corr[i] > 0.0 {
                                n_coeluting += 1;
                            }
                        }
                    }

                    // 5. Find the spectrum closest to apex RT for spectral scoring
                    let apex_local_idx = cand_spectra
                        .iter()
                        .enumerate()
                        .min_by(|a, b| {
                            (a.1.retention_time - peak.apex_rt)
                                .abs()
                                .total_cmp(&(b.1.retention_time - peak.apex_rt).abs())
                        })
                        .map(|(i, _)| i)?;
                    let _apex_global_idx = cand_global[apex_local_idx];
                    let apex_spectrum = cand_spectra[apex_local_idx];

                    let spectral_score = scorer.xcorr(apex_spectrum, entry);

                    // 6. Elution-weighted cosine (LibCosine at each scan,
                    //    weighted by reference XIC intensity²)
                    let elution_weighted_cosine = compute_elution_weighted_cosine(
                        &entry.fragments,
                        ref_xic,
                        &cand_spectra,
                        tol_da,
                        tol_ppm,
                        peak.start_rt,
                        peak.end_rt,
                    );

                    // 7. Mass accuracy at apex
                    let frag_matches = scorer.match_fragments(apex_spectrum, entry);
                    let (mass_mean, mass_abs_mean, mass_std) =
                        compute_mass_accuracy(&frag_matches, fragment_tolerance.unit);

                    // 7. RT deviation
                    let rt_dev = peak.apex_rt - expected_rt;

                    // 8. MS1 features (HRAM only)
                    let (ms1_coelution, ms1_isotope_cos) = if is_hram && !ms1_index.is_empty() {
                        let ms1_tol_ppm = calibration
                            .map(|cal| calibrated_tolerance_ppm(&cal.ms1_calibration, 10.0))
                            .unwrap_or(10.0);
                        let search_mz = calibration
                            .map(|cal| {
                                reverse_calibrate_mz(entry.precursor_mz, &cal.ms1_calibration)
                            })
                            .unwrap_or(entry.precursor_mz);

                        // MS1 precursor coelution: correlate MS1 M+0 XIC with reference fragment XIC
                        let mut ms1_coelution_val = 0.0;
                        let peak_rts: Vec<(f64, f64)> =
                            ref_xic[peak.start_index..=peak.end_index].to_vec();
                        if peak_rts.len() >= 3 {
                            let mut ms1_intensities = Vec::with_capacity(peak_rts.len());
                            let mut ref_intensities = Vec::with_capacity(peak_rts.len());
                            for &(rt, ref_val) in &peak_rts {
                                if let Some(ms1) = ms1_index.find_nearest(rt) {
                                    let intensity = ms1
                                        .find_peak_ppm(search_mz, ms1_tol_ppm)
                                        .map(|(_, int)| int as f64)
                                        .unwrap_or(0.0);
                                    ms1_intensities.push(intensity);
                                    ref_intensities.push(ref_val);
                                }
                            }
                            if ms1_intensities.len() >= 3 {
                                ms1_coelution_val =
                                    pearson_correlation_raw(&ms1_intensities, &ref_intensities);
                            }
                        }

                        // MS1 isotope cosine at apex RT
                        let mut iso_cos = 0.0;
                        if let Some(apex_ms1) = ms1_index.find_nearest(peak.apex_rt) {
                            let envelope = osprey_core::IsotopeEnvelope::extract(
                                apex_ms1,
                                search_mz,
                                entry.charge,
                                ms1_tol_ppm,
                            );
                            if envelope.has_m0() {
                                iso_cos = osprey_core::peptide_isotope_cosine(
                                    &entry.sequence,
                                    &envelope.intensities,
                                )
                                .unwrap_or(0.0);
                            }
                        }

                        (ms1_coelution_val, iso_cos)
                    } else {
                        (0.0, 0.0)
                    };

                    // 9. Peak shape features
                    let peak_width = peak.end_rt - peak.start_rt;
                    let peak_symmetry = {
                        let left_slice = &ref_xic[peak.start_index..=peak.apex_index];
                        let right_slice = &ref_xic[peak.apex_index..=peak.end_index];
                        let left_area = trapezoidal_area(left_slice);
                        let right_area = trapezoidal_area(right_slice);
                        if right_area > 1e-10 {
                            (left_area / right_area).min(10.0)
                        } else {
                            1.0
                        }
                    };
                    let peak_sharpness = {
                        // Average edge steepness (intensity change per minute)
                        let left_edge = if peak.apex_index > peak.start_index {
                            let dt = peak.apex_rt - peak.start_rt;
                            if dt > 1e-10 {
                                (peak.apex_intensity - ref_xic[peak.start_index].1) / dt
                            } else {
                                0.0
                            }
                        } else {
                            0.0
                        };
                        let right_edge = if peak.end_index > peak.apex_index {
                            let dt = peak.end_rt - peak.apex_rt;
                            if dt > 1e-10 {
                                (peak.apex_intensity - ref_xic[peak.end_index].1) / dt
                            } else {
                                0.0
                            }
                        } else {
                            0.0
                        };
                        (left_edge + right_edge) / 2.0
                    };

                    // 10. Median polish features
                    let (mp_cosine, mp_rsquared, mp_residual) = {
                        // Build XIC matrix within peak boundaries
                        let peak_xics: Vec<(usize, Vec<(f64, f64)>)> = xics
                            .iter()
                            .map(|(frag_idx, xic)| {
                                let peak_slice: Vec<(f64, f64)> =
                                    xic[peak.start_index..=peak.end_index].to_vec();
                                (*frag_idx, peak_slice)
                            })
                            .collect();

                        if let Some(ref mp) = tukey_median_polish(&peak_xics, 10, 0.01) {
                            let cos = median_polish_libcosine(mp, &entry.fragments);
                            let rsq = median_polish_rsquared(mp);
                            let res = median_polish_residual_ratio(mp);
                            (cos, rsq, res)
                        } else {
                            (0.0, 0.0, 1.0)
                        }
                    };

                    // 11. Peptide properties
                    let mod_count = entry.modifications.len() as u8;
                    let pep_len = entry.sequence.len() as u8;
                    let missed_cleav = entry
                        .sequence
                        .chars()
                        .take(entry.sequence.len().saturating_sub(1))
                        .filter(|&c| c == 'K' || c == 'R')
                        .count() as u8;

                    // Fragment m/z and intensities from library for blib output
                    // Use library theoretical values (not observed peaks) so Skyline
                    // can build XICs for the correct b/y ions
                    let frag_mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
                    let frag_intensities: Vec<f32> = entry
                        .fragments
                        .iter()
                        .map(|f| f.relative_intensity)
                        .collect();

                    // Reference XIC within peak bounds for blib coefficient series
                    let ref_xic_peak: Vec<(f64, f64)> =
                        ref_xic[peak.start_index..=peak.end_index].to_vec();

                    let features = CoelutionFeatureSet {
                        // Pairwise coelution
                        coelution_sum: corr_sum,
                        coelution_min: if corr_min.is_finite() { corr_min } else { 0.0 },
                        coelution_max: if corr_max.is_finite() { corr_max } else { 0.0 },
                        n_coeluting_fragments: n_coeluting,
                        n_fragment_pairs: n_pairs.min(255) as u8,
                        fragment_corr,

                        // Peak shape
                        peak_apex: peak.apex_intensity,
                        peak_area: peak.area,
                        peak_width,
                        peak_symmetry,
                        signal_to_noise: peak.signal_to_noise,
                        n_scans: peak_len as u16,
                        peak_sharpness,

                        // Spectral at apex
                        hyperscore: spectral_score.hyperscore,
                        xcorr: spectral_score.xcorr,
                        dot_product: spectral_score.lib_cosine,
                        dot_product_smz: spectral_score.lib_cosine_smz,
                        dot_product_top6: spectral_score.dot_product_top6,
                        dot_product_top5: spectral_score.dot_product_top5,
                        dot_product_top4: spectral_score.dot_product_top4,
                        dot_product_smz_top6: spectral_score.dot_product_smz_top6,
                        dot_product_smz_top5: spectral_score.dot_product_smz_top5,
                        dot_product_smz_top4: spectral_score.dot_product_smz_top4,
                        fragment_coverage: spectral_score.fragment_coverage,
                        sequence_coverage: spectral_score.sequence_coverage,
                        consecutive_ions: spectral_score.consecutive_ions as u8,
                        base_peak_rank: spectral_score.base_peak_rank as u8,
                        top6_matches: spectral_score.top6_matches as u8,
                        explained_intensity: spectral_score.explained_intensity,
                        elution_weighted_cosine,

                        // Mass accuracy
                        mass_accuracy_mean: mass_mean,
                        abs_mass_accuracy_mean: mass_abs_mean,
                        mass_accuracy_std: mass_std,

                        // RT deviation
                        rt_deviation: rt_dev,
                        abs_rt_deviation: rt_dev.abs(),

                        // MS1
                        ms1_precursor_coelution: ms1_coelution,
                        ms1_isotope_cosine: ms1_isotope_cos,

                        // Peptide properties
                        modification_count: mod_count,
                        peptide_length: pep_len,
                        missed_cleavages: missed_cleav,

                        // Median polish
                        median_polish_cosine: mp_cosine,
                        median_polish_rsquared: mp_rsquared,
                        median_polish_residual_ratio: mp_residual,
                    };

                    pb.inc(1);

                    Some(CoelutionScoredEntry {
                        entry_id: entry.id,
                        is_decoy: entry.is_decoy,
                        sequence: entry.sequence.clone(),
                        modified_sequence: entry.modified_sequence.clone(),
                        charge: entry.charge,
                        precursor_mz: entry.precursor_mz,
                        protein_ids: entry.protein_ids.clone(),
                        scan_number: apex_spectrum.scan_number,
                        apex_rt: peak.apex_rt,
                        peak_bounds: peak,
                        features,
                        fragment_mzs: frag_mzs,
                        fragment_intensities: frag_intensities,
                        reference_xic: ref_xic_peak,
                        file_name: file_name.to_string(),
                        run_qvalue: 1.0,
                        experiment_qvalue: 1.0,
                        score: 0.0,
                        pep: 1.0,
                    })
                })
                .collect();

            pb.inc(candidate_indices.len().saturating_sub(window_entries.len()) as u64);
            window_entries
        })
        .collect();

    pb.finish_with_message("Done");

    log::info!(
        "Coelution search complete: {} scored entries ({} targets, {} decoys)",
        all_entries.len(),
        all_entries.iter().filter(|e| !e.is_decoy).count(),
        all_entries.iter().filter(|e| e.is_decoy).count(),
    );

    // Keep best score per precursor (a precursor can be scored in multiple
    // overlapping DIA windows; keep the one with highest coelution_sum)
    let pre_count = all_entries.len();
    let mut best_by_id: HashMap<u32, CoelutionScoredEntry> = HashMap::new();
    for entry in all_entries {
        let is_better = best_by_id
            .get(&entry.entry_id)
            .map(|existing| entry.features.coelution_sum > existing.features.coelution_sum)
            .unwrap_or(true);
        if is_better {
            best_by_id.insert(entry.entry_id, entry);
        }
    }

    let mut deduped: Vec<CoelutionScoredEntry> = best_by_id.into_values().collect();

    // Sort by (entry_id, scan_number) for deterministic output regardless of
    // Rayon thread scheduling or HashMap iteration order
    deduped.sort_by(|a, b| {
        a.entry_id
            .cmp(&b.entry_id)
            .then(a.scan_number.cmp(&b.scan_number))
    });

    if deduped.len() < pre_count {
        log::info!(
            "Best per precursor: {} entries (removed {} duplicate window hits)",
            deduped.len(),
            pre_count - deduped.len()
        );
    }

    Ok(deduped)
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

        // Safety check: NaN expected_rt values would bypass binary search
        // (fixed by handling duplicate library_rts in RTCalibration::predict)
        if cfg!(debug_assertions) {
            let nan_count: usize = bins
                .values()
                .flat_map(|v| v.iter())
                .filter(|(rt, _)| rt.is_nan())
                .count();
            if nan_count > 0 {
                log::error!(
                    "MzRTIndex has {} NaN expected_rt entries — calibration predict() bug!",
                    nan_count
                );
            }
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
/// 5. Optional top-N fragment pre-filter (requires at least 2 of top 6 peaks in spectrum)
///
/// The MzRTIndex stores entries sorted by expected_rt (calibrated if available),
/// enabling O(log n + k) candidate selection instead of O(n) per m/z bin.
#[cfg(test)]
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

                // Apply top-N fragment pre-filter if enabled
                if let Some(ref frag_tol) = fragment_tolerance {
                    if !osprey_scoring::has_topn_fragment_match(
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
                .is_some_and(|entries| entries.iter().any(|&(_, i)| i == idx))
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
}
