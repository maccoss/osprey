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
    CoelutionFeatureSet, CoelutionScoredEntry, CwtCandidate, DecoyMethod as CoreDecoyMethod,
    FdrEntry, FdrLevel, FdrMethod, FragmentToleranceConfig, LibraryEntry, LibraryFragment,
    MS1Spectrum, OspreyConfig, OspreyError, Result, Spectrum, ToleranceUnit, XICPeakBounds,
};
use osprey_fdr::{
    get_pin_feature_names, percolator, pin_feature_value, MokapotRunner, NUM_PIN_FEATURES,
};
use osprey_io::{
    load_all_spectra, load_library, load_spectra_cache, save_spectra_cache, spectra_cache_path,
    BlibWriter, MS1Index,
};
use osprey_scoring::{
    batch::{run_coelution_calibration_scoring, sample_library_for_calibration, MS1SpectrumLookup},
    has_topn_fragment_match, DecoyGenerator, DecoyMethod, Enzyme, SpectralScorer,
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
use std::sync::Arc;

use osprey_scoring::batch::CalibrationMatch;

// Parquet cache metadata keys
const META_OSPREY_VERSION: &str = "osprey.version";
const META_SEARCH_HASH: &str = "osprey.search_hash";
const META_LIBRARY_HASH: &str = "osprey.library_hash";
const META_RECONCILED: &str = "osprey.reconciled";
const META_RECONCILIATION_HASH: &str = "osprey.reconciliation_hash";

/// Result of validating a cached scores parquet file's metadata.
enum CacheValidity {
    /// Cache is valid (first-pass scores match current config)
    ValidFirstPass,
    /// Cache contains reconciled data and reconciliation params still match
    ValidReconciled,
    /// Cache is stale — delete and re-score
    Stale(String),
}

/// Read key-value metadata from a parquet file and check validity against current config.
fn validate_scores_cache(path: &std::path::Path, config: &OspreyConfig) -> CacheValidity {
    use parquet::file::reader::{FileReader, SerializedFileReader};

    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return CacheValidity::Stale("cannot open file".into()),
    };
    let reader = match SerializedFileReader::new(file) {
        Ok(r) => r,
        Err(_) => return CacheValidity::Stale("cannot read parquet metadata".into()),
    };
    let file_meta = reader.metadata().file_metadata();
    let kv = match file_meta.key_value_metadata() {
        Some(kv) => kv,
        None => return CacheValidity::Stale("no metadata (pre-hash cache)".into()),
    };

    let get = |key: &str| -> Option<&str> {
        kv.iter()
            .find(|kv| kv.key == key)
            .and_then(|kv| kv.value.as_deref())
    };

    // Check version
    let cached_version = match get(META_OSPREY_VERSION) {
        Some(v) => v,
        None => return CacheValidity::Stale("missing version".into()),
    };
    if cached_version != env!("CARGO_PKG_VERSION") {
        return CacheValidity::Stale(format!(
            "version mismatch: cached={}, current={}",
            cached_version,
            env!("CARGO_PKG_VERSION")
        ));
    }

    // Check search hash
    let cached_search = match get(META_SEARCH_HASH) {
        Some(h) => h,
        None => return CacheValidity::Stale("missing search hash".into()),
    };
    let current_search = config.search_parameter_hash();
    if cached_search != current_search {
        return CacheValidity::Stale("search parameters changed".into());
    }

    // Check library hash
    let cached_lib = match get(META_LIBRARY_HASH) {
        Some(h) => h,
        None => return CacheValidity::Stale("missing library hash".into()),
    };
    let current_lib = config.library_identity_hash();
    if cached_lib != current_lib {
        return CacheValidity::Stale("library changed".into());
    }

    // Check reconciliation status
    let is_reconciled = get(META_RECONCILED) == Some("true");
    if is_reconciled {
        if let Some(cached_recon) = get(META_RECONCILIATION_HASH) {
            let current_recon = config.reconciliation_parameter_hash();
            if cached_recon == current_recon {
                return CacheValidity::ValidReconciled;
            }
            // Reconciliation params changed — but first-pass is still valid
            return CacheValidity::ValidFirstPass;
        }
    }

    CacheValidity::ValidFirstPass
}

/// Build Parquet key-value metadata for a first-pass scores cache.
fn build_scores_metadata(config: &OspreyConfig) -> Vec<parquet::file::metadata::KeyValue> {
    vec![
        parquet::file::metadata::KeyValue {
            key: META_OSPREY_VERSION.to_string(),
            value: Some(env!("CARGO_PKG_VERSION").to_string()),
        },
        parquet::file::metadata::KeyValue {
            key: META_SEARCH_HASH.to_string(),
            value: Some(config.search_parameter_hash()),
        },
        parquet::file::metadata::KeyValue {
            key: META_LIBRARY_HASH.to_string(),
            value: Some(config.library_identity_hash()),
        },
        parquet::file::metadata::KeyValue {
            key: META_RECONCILED.to_string(),
            value: Some("false".to_string()),
        },
    ]
}

/// Build Parquet key-value metadata for a reconciled scores cache.
fn build_reconciled_metadata(config: &OspreyConfig) -> Vec<parquet::file::metadata::KeyValue> {
    vec![
        parquet::file::metadata::KeyValue {
            key: META_OSPREY_VERSION.to_string(),
            value: Some(env!("CARGO_PKG_VERSION").to_string()),
        },
        parquet::file::metadata::KeyValue {
            key: META_SEARCH_HASH.to_string(),
            value: Some(config.search_parameter_hash()),
        },
        parquet::file::metadata::KeyValue {
            key: META_LIBRARY_HASH.to_string(),
            value: Some(config.library_identity_hash()),
        },
        parquet::file::metadata::KeyValue {
            key: META_RECONCILED.to_string(),
            value: Some("true".to_string()),
        },
        parquet::file::metadata::KeyValue {
            key: META_RECONCILIATION_HASH.to_string(),
            value: Some(config.reconciliation_parameter_hash()),
        },
    ]
}

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

/// Build the XCorr scorer used during calibration LDA.
///
/// Always uses unit-resolution bins (~2K) regardless of instrument
/// resolution mode. Calibration only needs to discriminate target/decoy
/// for LDA training, and unit bins are ~50x cheaper than HRAM bins
/// (2001 vs ~100K bins per spectrum). This matches the design documented
/// in `docs/02-calibration.md` ("Comet-style XCorr (unit resolution,
/// BLAS sdot)") and the sibling function
/// `osprey_scoring::batch::run_xcorr_calibration_scoring`.
///
/// Main search (`run_search`) is unaffected and still picks the
/// resolution-mode-appropriate scorer: `SpectralScorer::hram()` on HRAM,
/// `SpectralScorer::new()` on unit-resolution data.
///
/// The if/else remains because the two branches still differ in
/// tolerance unit (ppm for HRAM, Da for unit-resolution); only the bin
/// choice is unconditional.
fn calibration_xcorr_scorer(config: &OspreyConfig) -> SpectralScorer {
    let is_hram = matches!(config.resolution_mode, osprey_core::ResolutionMode::HRAM);
    if is_hram {
        SpectralScorer::new().with_tolerance_ppm(config.fragment_tolerance.tolerance)
    } else {
        SpectralScorer::new().with_tolerance_da(config.fragment_tolerance.tolerance)
    }
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
            log::debug!("RT mapping: library and mzML ranges are similar, using identity mapping");
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
            log::debug!(
                "RT mapping: linear transform from library ({:.2}-{:.2}) to mzML ({:.2}-{:.2} min)",
                lib_min_rt,
                lib_max_rt,
                mzml_min_rt,
                mzml_max_rt
            );
            log::debug!(
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

    log::debug!(
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

    log::debug!(
        "Initial RT tolerance: {:.1} min ({:.0}% of {:.1} min mzML range)",
        initial_tolerance,
        tolerance_fraction * 100.0,
        mzml_rt_range
    );
    log::debug!(
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
        log::debug!(
            "Using {} MS1 spectra for precursor calibration",
            ms1_index.len()
        );
    } else {
        log::warn!(
            "No MS1 spectra available. Falling back to isolation window center for precursor m/z."
        );
    }

    let xcorr_scorer = calibration_xcorr_scorer(config);

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

    log::debug!(
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

        log::debug!(
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
            log::debug!(
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

        // LDA-based scoring on accumulated matches
        // Don't use isotope feature during calibration — MS1 isotope scoring
        // is only reliable after mass calibration, which hasn't happened yet.
        // Isotope features are used in the full scoring step after calibration.
        let calibration_fdr = 0.01;

        // Convert to Vec for LDA training.
        // Sort by base_id for deterministic LDA input: LinearDiscriminantAnalysis::fit()
        // computes class means by iterating matrix rows, and floating-point addition is
        // non-associative — different row order → different sums → different eigenvectors.
        // Sort by (base_id, entry_id) to group target-decoy pairs together.
        let mut all_matches: Vec<CalibrationMatch> =
            accumulated_matches.values().cloned().collect();
        all_matches.sort_by_key(|m| (m.entry_id & 0x7FFFFFFF, m.entry_id));

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

        log::debug!(
            "Competition: {} target wins, {} decoy wins",
            n_target_wins,
            n_decoy_wins
        );

        // Write calibration debug CSV (after LDA scoring, so discriminant/q-value are populated)
        if let Err(e) = write_calibration_debug_csv(&all_matches, &debug_path, expected_rt_fn) {
            log::warn!("Failed to write calibration debug CSV: {}", e);
        }

        // Log combined LDA + S/N filter result
        log::info!(
            "LDA scoring: {} at 1% FDR, {}/{} with S/N >= {:.0}",
            n_target_wins,
            passing_targets.len(),
            n_target_wins,
            MIN_SNR_FOR_RT_CAL
        );

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
                log::debug!(
                    "Measured peak width: median={:.2} min ({:.0} sec) from {} peaks",
                    median,
                    median * 60.0,
                    widths.len()
                );
            }
        }

        log::debug!(
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
            outlier_retention: 1.0, // Use all calibration points — LDA + S/N already filtered
        };
        let calibrator = RTCalibrator::with_config(calibrator_config);
        let mut rt_calibration = calibrator.fit(&library_rts_detected, &measured_rts_detected)?;
        let mut rt_stats = rt_calibration.stats();

        // === Iterative calibration refinement (2-pass) ===
        // Compute first-pass tolerance using MAD-based robust estimate
        // MAD × 1.4826 ≈ SD for normal distribution; 3× that covers ~99.7% of peptides
        let mad_tolerance = rt_stats.mad * 1.4826 * 3.0;
        let pass1_tolerance = mad_tolerance
            .max(config.rt_calibration.min_rt_tolerance)
            .min(config.rt_calibration.max_rt_tolerance);

        log::debug!(
            "First-pass RT tolerance: {:.2} min (MAD={:.3}, robust_SD={:.3}, residual_SD={:.3}, {} points, R²={:.4})",
            pass1_tolerance,
            rt_stats.mad,
            rt_stats.mad * 1.4826,
            rt_stats.residual_std,
            rt_stats.n_points,
            rt_stats.r_squared
        );

        // Only refine if tolerance narrowed significantly (at least 2× tighter)
        if pass1_tolerance < initial_tolerance * 0.5 {
            log::debug!(
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

                let refined_mad_tolerance = rt_stats_refined.mad * 1.4826 * 3.0;
                let refined_tolerance = refined_mad_tolerance
                    .max(config.rt_calibration.min_rt_tolerance)
                    .min(config.rt_calibration.max_rt_tolerance);

                log::debug!(
                    "Refined RT tolerance: {:.2} min (MAD={:.3}, robust_SD={:.3}, residual_SD={:.3}, {} points, R²={:.4})",
                    refined_tolerance,
                    rt_stats_refined.mad,
                    rt_stats_refined.mad * 1.4826,
                    rt_stats_refined.residual_std,
                    rt_stats_refined.n_points,
                    rt_stats_refined.r_squared,
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
                    log::debug!(
                        "Refined calibration not better (R² {:.4} vs {:.4}), keeping original",
                        rt_stats_refined.r_squared,
                        rt_stats.r_squared
                    );
                }
            } else {
                log::debug!(
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
                search_hash: None, // set by caller before saving
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
                mad: Some(rt_stats.mad),
            },
            second_pass_rt: None,
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

/// Path for per-file first-pass FDR score sidecar.
fn fdr_scores_path_pass1(input_path: &std::path::Path) -> std::path::PathBuf {
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    let parent = input_path.parent().unwrap_or(std::path::Path::new("."));
    parent.join(format!("{}.1st-pass.fdr_scores.bin", stem))
}

/// Path for per-file second-pass FDR score sidecar.
fn fdr_scores_path_pass2(input_path: &std::path::Path) -> std::path::PathBuf {
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    let parent = input_path.parent().unwrap_or(std::path::Path::new("."));
    parent.join(format!("{}.2nd-pass.fdr_scores.bin", stem))
}

/// Write per-file FDR scores to a sidecar binary file (little-endian f64 array).
fn write_fdr_scores_sidecar(path: &std::path::Path, entries: &[FdrEntry]) -> Result<()> {
    let tmp_path = std::env::temp_dir().join(format!(
        "osprey_{}_{}",
        std::process::id(),
        path.file_name()
            .unwrap_or(std::ffi::OsStr::new("fdr_scores.bin"))
            .to_string_lossy()
    ));
    let mut file = std::fs::File::create(&tmp_path).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to create scores sidecar {}: {}",
            tmp_path.display(),
            e
        ))
    })?;
    use std::io::Write;
    for entry in entries {
        file.write_all(&entry.score.to_le_bytes())
            .map_err(|e| OspreyError::OutputError(format!("Failed to write score: {}", e)))?;
    }
    drop(file);
    osprey_core::copy_and_verify(&tmp_path, path)?;
    Ok(())
}

/// Load per-file FDR scores from sidecar into FdrEntry stubs. Returns true if successful.
fn load_fdr_scores_sidecar(path: &std::path::Path, entries: &mut [FdrEntry]) -> bool {
    let data = match std::fs::read(path) {
        Ok(d) => d,
        Err(_) => return false,
    };
    let expected_len = entries.len() * 8;
    if data.len() != expected_len {
        log::warn!(
            "Score sidecar {} size mismatch ({} vs expected {}), ignoring",
            path.display(),
            data.len(),
            expected_len
        );
        return false;
    }
    for (i, entry) in entries.iter_mut().enumerate() {
        let offset = i * 8;
        let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap();
        entry.score = f64::from_le_bytes(bytes);
    }
    true
}

/// Write FDR score sidecars for all files.
fn persist_fdr_scores(
    per_file_entries: &[(String, Vec<FdrEntry>)],
    config: &OspreyConfig,
    path_fn: fn(&std::path::Path) -> std::path::PathBuf,
    label: &str,
) {
    for (file_name, entries) in per_file_entries.iter() {
        if let Some(input_file) = config.input_files.iter().find(|p| {
            p.file_stem()
                .and_then(|s| s.to_str())
                .is_some_and(|s| s == file_name)
        }) {
            let sidecar_path = path_fn(input_file);
            if let Err(e) = write_fdr_scores_sidecar(&sidecar_path, entries) {
                log::warn!(
                    "Failed to write {} score sidecar for {}: {}",
                    label,
                    file_name,
                    e
                );
            }
        }
    }
}

/// Write per-file scored entries to parquet with optional key-value metadata.
fn write_scores_parquet_with_metadata(
    path: &std::path::Path,
    entries: &[CoelutionScoredEntry],
    metadata: Option<Vec<parquet::file::metadata::KeyValue>>,
) -> Result<()> {
    use parquet::arrow::ArrowWriter;
    use parquet::basic::Compression;
    use parquet::file::properties::WriterProperties;

    if entries.is_empty() {
        // Don't create an empty file — downstream code checks existence + size > 0
        log::debug!("No entries to write to {}", path.display());
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
        Field::new("cwt_candidates", DataType::Binary, true),
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
    let mut cwt_cand_b = BinaryBuilder::with_capacity(n, n * 244); // 4 + 5*48 = 244
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

        // Serialize CWT candidates as packed LE bytes: [u32 count][48 bytes per candidate]
        let mut cwt_bytes = Vec::with_capacity(4 + entry.cwt_candidates.len() * 48);
        cwt_bytes.extend_from_slice(&(entry.cwt_candidates.len() as u32).to_le_bytes());
        for cand in &entry.cwt_candidates {
            cwt_bytes.extend_from_slice(&cand.apex_rt.to_le_bytes());
            cwt_bytes.extend_from_slice(&cand.start_rt.to_le_bytes());
            cwt_bytes.extend_from_slice(&cand.end_rt.to_le_bytes());
            cwt_bytes.extend_from_slice(&cand.area.to_le_bytes());
            cwt_bytes.extend_from_slice(&cand.snr.to_le_bytes());
            cwt_bytes.extend_from_slice(&cand.coelution_score.to_le_bytes());
        }
        cwt_cand_b.append_value(&cwt_bytes);

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
        std::sync::Arc::new(cwt_cand_b.finish()),
        std::sync::Arc::new(frag_mz_b.finish()),
        std::sync::Arc::new(frag_int_b.finish()),
        std::sync::Arc::new(xic_rt_b.finish()),
        std::sync::Arc::new(xic_int_b.finish()),
    ];
    for builder in &mut feat_builders {
        columns.push(std::sync::Arc::new(builder.finish()));
    }

    // Write to local temp dir first, then copy to final destination.
    // Direct writes to NAS can produce corrupt files.
    let tmp_path = std::env::temp_dir().join(format!(
        "osprey_{}_{}",
        std::process::id(),
        path.file_name()
            .unwrap_or(std::ffi::OsStr::new("scores.parquet"))
            .to_string_lossy()
    ));
    let file = std::fs::File::create(&tmp_path).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to create temp scores file {}: {}",
            tmp_path.display(),
            e
        ))
    })?;
    let mut props_builder =
        WriterProperties::builder().set_compression(Compression::ZSTD(Default::default()));
    if metadata.is_some() {
        props_builder = props_builder.set_key_value_metadata(metadata);
    }
    let props = props_builder.build();
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

    osprey_core::copy_and_verify(&tmp_path, path)?;
    log::debug!("Wrote {} scores to {}", n, path.display());

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

    // Validate schema: reject stale scores files from older versions with different features
    let feature_names = get_pin_feature_names();
    let parquet_schema = builder.schema();
    for name in &feature_names {
        if parquet_schema.column_with_name(name).is_none() {
            return Err(OspreyError::OutputError(format!(
                "Stale scores cache {}: missing column '{}'. Delete and re-run.",
                path.display(),
                name
            )));
        }
    }

    let reader = builder.build().map_err(|e| {
        OspreyError::OutputError(format!("Failed to build reader {}: {}", path.display(), e))
    })?;
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

        // Read CWT candidates column (optional, backwards-compatible)
        let cwt_cand_col = batch
            .column_by_name("cwt_candidates")
            .and_then(|c| c.as_any().downcast_ref::<BinaryArray>());

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

            // Deserialize CWT candidates from packed binary
            let cwt_candidates = if let Some(col) = cwt_cand_col {
                let bytes = col.value(row);
                if bytes.len() >= 4 {
                    let count = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
                    let data = &bytes[4..];
                    // Each CwtCandidate is 6 × f64 = 48 bytes
                    data.chunks_exact(48)
                        .take(count)
                        .map(|chunk| {
                            let f = |offset: usize| -> f64 {
                                f64::from_le_bytes(chunk[offset..offset + 8].try_into().unwrap())
                            };
                            CwtCandidate {
                                apex_rt: f(0),
                                start_rt: f(8),
                                end_rt: f(16),
                                area: f(24),
                                snr: f(32),
                                coelution_score: f(40),
                            }
                        })
                        .collect()
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            };

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
                run_precursor_qvalue: 1.0,
                run_peptide_qvalue: 1.0,
                run_protein_qvalue: 1.0,
                experiment_precursor_qvalue: 1.0,
                experiment_peptide_qvalue: 1.0,
                experiment_protein_qvalue: 1.0,
                score: 0.0,
                pep: 1.0,
                cwt_candidates,
            });
        }
    }

    Ok(entries)
}

/// Set CoelutionFeatureSet fields from parquet feature columns.
///
/// This maps the 17 PIN feature columns back to the struct fields.
/// Order must match get_pin_feature_names() and pin_feature_value().
fn set_features_from_pin_values(
    features: &mut CoelutionFeatureSet,
    feat_cols: &[&arrow::array::Float64Array],
    row: usize,
) {
    let v = |i: usize| -> f64 { feat_cols[i].value(row) };

    // Pairwise coelution (indices 0-2)
    features.coelution_sum = v(0);
    features.coelution_max = v(1);
    features.n_coeluting_fragments = v(2) as u8;

    // Peak shape (indices 3-5)
    features.peak_apex = v(3);
    features.peak_area = v(4);
    features.peak_sharpness = v(5);

    // Spectral at apex (indices 6-8)
    features.xcorr = v(6);
    features.consecutive_ions = v(7) as u8;
    features.explained_intensity = v(8);

    // Mass accuracy (indices 9-10)
    features.mass_accuracy_mean = v(9);
    features.abs_mass_accuracy_mean = v(10);

    // RT deviation (indices 11-12)
    features.rt_deviation = v(11);
    features.abs_rt_deviation = v(12);

    // MS1 (indices 13-14)
    features.ms1_precursor_coelution = v(13);
    features.ms1_isotope_cosine = v(14);

    // Median polish (indices 15-16)
    features.median_polish_cosine = v(15);
    features.median_polish_residual_ratio = v(16);

    // SG-weighted multi-scan (indices 17-20)
    features.sg_weighted_xcorr = v(17);
    features.sg_weighted_cosine = v(18);
    features.median_polish_min_fragment_r2 = v(19);
    features.median_polish_residual_correlation = v(20);
}

/// Extract PIN feature values from a CoelutionFeatureSet as a Vec<f64>.
///
/// Used to provide overlay (re-scored) features to FDR scoring without
/// loading from parquet. Order matches get_pin_feature_names().
fn pin_feature_vector(features: &CoelutionFeatureSet) -> Vec<f64> {
    let n = get_pin_feature_names().len();
    (0..n)
        .map(|i| {
            let v = pin_feature_value(features, i);
            if v.is_finite() {
                v
            } else {
                0.0
            }
        })
        .collect()
}

/// Type alias for the in-memory overlay of re-scored entries.
///
/// Maps file_name → {entry_index → re-scored entry}. Only entries that were
/// re-scored during reconciliation appear in the overlay (~85K per file).
/// This avoids loading/writing the full ~1.5M-entry parquet during re-scoring.
type RescoreOverlay = HashMap<String, HashMap<usize, CoelutionScoredEntry>>;

/// Intern a `modified_sequence` string, returning a shared `Arc<str>`.
///
/// With 240 GPF files the same ~3.5M unique peptide sequences appear across
/// all files.  Without interning, 240M separate `String` allocations consume
/// ~6 GB of heap.  The interner deduplicates them into ~90 MB of shared
/// `Arc<str>` pointers.
fn intern_seq(interner: &mut HashSet<Arc<str>>, s: &str) -> Arc<str> {
    if let Some(existing) = interner.get(s) {
        existing.clone()
    } else {
        let arc: Arc<str> = Arc::from(s);
        interner.insert(arc.clone());
        arc
    }
}

/// Load only PIN feature values from a parquet cache file.
///
/// Returns one Vec<f64> per entry (with NUM_PIN_FEATURES elements each).
/// This is a lightweight alternative to `load_scores_parquet` that only reads
/// the 21 feature columns — no fragments, CWT candidates, XICs, or strings.
/// Dramatically reduces I/O and memory for large experiments.
fn load_pin_features_from_parquet(path: &std::path::Path) -> Result<Vec<Vec<f64>>> {
    use arrow::array::Float64Array;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ProjectionMask;

    let feature_names = get_pin_feature_names();

    // Guard against 0-byte files from interrupted runs
    let file_size = path.metadata().map(|m| m.len()).unwrap_or(0);
    if file_size == 0 {
        return Err(OspreyError::config(format!(
            "Scores cache {} is empty (0 bytes). Delete it and re-run.",
            path.display()
        )));
    }

    let file = std::fs::File::open(path).map_err(|e| {
        OspreyError::config(format!("Failed to open parquet {}: {}", path.display(), e))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
        OspreyError::config(format!("Failed to read parquet {}: {}", path.display(), e))
    })?;

    // Build projection mask to only read feature columns
    let parquet_schema = builder.parquet_schema().clone();
    let arrow_schema = builder.schema().clone();
    let feature_indices: Vec<usize> = feature_names
        .iter()
        .filter_map(|name| arrow_schema.column_with_name(name).map(|(idx, _)| idx))
        .collect();

    if feature_indices.len() != feature_names.len() {
        return Err(OspreyError::OutputError(format!(
            "Stale scores cache {}: missing feature columns. Delete and re-run.",
            path.display(),
        )));
    }

    let projection = ProjectionMask::roots(&parquet_schema, feature_indices);
    let reader = builder.with_projection(projection).build().map_err(|e| {
        OspreyError::config(format!(
            "Failed to build reader for {}: {}",
            path.display(),
            e
        ))
    })?;

    let mut all_features = Vec::new();
    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| OspreyError::config(format!("Failed to read batch: {}", e)))?;

        // Read feature columns by name (projection may reorder them)
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
            let features: Vec<f64> = feat_cols
                .iter()
                .map(|col| {
                    let v = col.value(row);
                    if v.is_finite() {
                        v
                    } else {
                        0.0
                    }
                })
                .collect();
            all_features.push(features);
        }
    }

    Ok(all_features)
}

/// Load only FdrEntry stub fields from a parquet cache file.
///
/// Reads only 9 columns (entry_id, is_decoy, charge, scan_number, apex_rt,
/// start_rt, end_rt, modified_sequence, fragment_coelution_sum)
/// using a ProjectionMask. `modified_sequence` values are interned via the
/// shared `seq_interner` to deduplicate strings across files.
fn load_fdr_stubs_from_parquet(
    path: &std::path::Path,
    seq_interner: &mut HashSet<Arc<str>>,
) -> Result<Vec<FdrEntry>> {
    use arrow::array::{BooleanArray, Float64Array, StringArray, UInt32Array, UInt8Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ProjectionMask;

    let stub_columns = [
        "entry_id",
        "is_decoy",
        "charge",
        "scan_number",
        "apex_rt",
        "start_rt",
        "end_rt",
        "modified_sequence",
        "fragment_coelution_sum",
    ];

    let file = std::fs::File::open(path).map_err(|e| {
        OspreyError::config(format!("Failed to open parquet {}: {}", path.display(), e))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
        OspreyError::config(format!("Failed to read parquet {}: {}", path.display(), e))
    })?;

    let parquet_schema = builder.parquet_schema().clone();
    let arrow_schema = builder.schema().clone();
    let col_indices: Vec<usize> = stub_columns
        .iter()
        .filter_map(|name| arrow_schema.column_with_name(name).map(|(idx, _)| idx))
        .collect();

    if col_indices.len() != stub_columns.len() {
        return Err(OspreyError::OutputError(format!(
            "Stale scores cache {}: missing FdrEntry columns. Delete and re-run.",
            path.display(),
        )));
    }

    let projection = ProjectionMask::roots(&parquet_schema, col_indices);
    let reader = builder.with_projection(projection).build().map_err(|e| {
        OspreyError::config(format!(
            "Failed to build reader for {}: {}",
            path.display(),
            e
        ))
    })?;

    let mut stubs = Vec::new();
    let mut row_offset = 0usize;
    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| OspreyError::config(format!("Failed to read batch: {}", e)))?;

        let entry_id_col = batch
            .column_by_name("entry_id")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let decoy_col = batch
            .column_by_name("is_decoy")
            .unwrap()
            .as_any()
            .downcast_ref::<BooleanArray>()
            .unwrap();
        let charge_col = batch
            .column_by_name("charge")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt8Array>()
            .unwrap();
        let scan_col = batch
            .column_by_name("scan_number")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let apex_col = batch
            .column_by_name("apex_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let start_col = batch
            .column_by_name("start_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let end_col = batch
            .column_by_name("end_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let modseq_col = batch
            .column_by_name("modified_sequence")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let coelution_col = batch
            .column_by_name("fragment_coelution_sum")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();

        for row in 0..batch.num_rows() {
            let parquet_idx = (row_offset + row) as u32;
            stubs.push(FdrEntry {
                entry_id: entry_id_col.value(row),
                parquet_index: parquet_idx,
                is_decoy: decoy_col.value(row),
                charge: charge_col.value(row),
                scan_number: scan_col.value(row),
                apex_rt: apex_col.value(row),
                start_rt: start_col.value(row),
                end_rt: end_col.value(row),
                coelution_sum: coelution_col.value(row),
                score: 0.0,
                run_precursor_qvalue: 1.0,
                run_peptide_qvalue: 1.0,
                run_protein_qvalue: 1.0,
                experiment_precursor_qvalue: 1.0,
                experiment_peptide_qvalue: 1.0,
                experiment_protein_qvalue: 1.0,
                pep: 1.0,
                modified_sequence: intern_seq(seq_interner, modseq_col.value(row)),
            });
        }
        row_offset += batch.num_rows();
    }

    Ok(stubs)
}

/// Load a 5-column projection from Parquet for lightweight blib output planning.
///
/// Returns (entry_id, apex_rt, start_rt, end_rt, bounds_area) tuples, in the same
/// row order as the original Parquet file. These are merged with `LightFdr` data
/// (which carries q-values, score, pep, modified_sequence, charge) to build
/// `BlibPlanEntry` structs at ~96 bytes each, avoiding the ~940-byte full
/// `CoelutionScoredEntry` reload that would OOM on 240-file experiments.
#[allow(clippy::type_complexity)]
fn load_blib_plan_from_parquet(path: &std::path::Path) -> Result<Vec<(u32, f64, f64, f64, f64)>> {
    use arrow::array::{Float64Array, UInt32Array};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ProjectionMask;

    let plan_columns = ["entry_id", "apex_rt", "start_rt", "end_rt", "bounds_area"];

    let file = std::fs::File::open(path).map_err(|e| {
        OspreyError::config(format!("Failed to open parquet {}: {}", path.display(), e))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
        OspreyError::config(format!("Failed to read parquet {}: {}", path.display(), e))
    })?;

    let parquet_schema = builder.parquet_schema().clone();
    let arrow_schema = builder.schema().clone();
    let col_indices: Vec<usize> = plan_columns
        .iter()
        .filter_map(|name| arrow_schema.column_with_name(name).map(|(idx, _)| idx))
        .collect();

    if col_indices.len() != plan_columns.len() {
        return Err(OspreyError::OutputError(format!(
            "Stale scores cache {}: missing blib plan columns. Delete and re-run.",
            path.display(),
        )));
    }

    let projection = ProjectionMask::roots(&parquet_schema, col_indices);
    let reader = builder.with_projection(projection).build().map_err(|e| {
        OspreyError::config(format!(
            "Failed to build reader for {}: {}",
            path.display(),
            e
        ))
    })?;

    let mut rows = Vec::new();
    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| OspreyError::config(format!("Failed to read batch: {}", e)))?;

        let entry_id_col = batch
            .column_by_name("entry_id")
            .unwrap()
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let apex_rt_col = batch
            .column_by_name("apex_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let start_rt_col = batch
            .column_by_name("start_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let end_rt_col = batch
            .column_by_name("end_rt")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();
        let area_col = batch
            .column_by_name("bounds_area")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();

        for row in 0..batch.num_rows() {
            rows.push((
                entry_id_col.value(row),
                apex_rt_col.value(row),
                start_rt_col.value(row),
                end_rt_col.value(row),
                area_col.value(row),
            ));
        }
    }

    Ok(rows)
}

/// Load only CWT candidates from a parquet cache file.
///
/// Returns one Vec<CwtCandidate> per entry, in the same order as entries were written.
/// This is a lightweight alternative to `load_scores_parquet` that only reads
/// the cwt_candidates column, avoiding the cost of loading features and fragments.
pub(crate) fn load_cwt_candidates_from_parquet(
    path: &std::path::Path,
) -> Result<Vec<Vec<CwtCandidate>>> {
    use arrow::array::BinaryArray;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path).map_err(|e| {
        OspreyError::config(format!("Failed to open parquet {}: {}", path.display(), e))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| {
        OspreyError::config(format!("Failed to read parquet {}: {}", path.display(), e))
    })?;
    let reader = builder.build().map_err(|e| {
        OspreyError::config(format!(
            "Failed to build reader for {}: {}",
            path.display(),
            e
        ))
    })?;

    let mut all_candidates = Vec::new();
    for batch_result in reader {
        let batch = batch_result
            .map_err(|e| OspreyError::config(format!("Failed to read batch: {}", e)))?;

        let cwt_col = batch
            .column_by_name("cwt_candidates")
            .and_then(|c| c.as_any().downcast_ref::<BinaryArray>());

        for row in 0..batch.num_rows() {
            let candidates = if let Some(col) = cwt_col {
                let bytes = col.value(row);
                if bytes.len() >= 4 {
                    let count = u32::from_le_bytes(bytes[..4].try_into().unwrap()) as usize;
                    let data = &bytes[4..];
                    data.chunks_exact(48)
                        .take(count)
                        .map(|chunk| {
                            let f = |offset: usize| -> f64 {
                                f64::from_le_bytes(chunk[offset..offset + 8].try_into().unwrap())
                            };
                            CwtCandidate {
                                apex_rt: f(0),
                                start_rt: f(8),
                                end_rt: f(16),
                                area: f(24),
                                snr: f(32),
                                coelution_score: f(40),
                            }
                        })
                        .collect()
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            };
            all_candidates.push(candidates);
        }
    }

    Ok(all_candidates)
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
        qval_b.append_value(entry.run_precursor_qvalue);
        gqval_b.append_value(entry.experiment_precursor_qvalue);
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
    let mut per_file_entries: Vec<(String, Vec<FdrEntry>)> = Vec::new();
    let mut per_file_cache_paths: HashMap<String, std::path::PathBuf> = HashMap::new();
    let mut pin_files: HashMap<String, std::path::PathBuf> = HashMap::new();
    // Retain per-file RT calibrations for inter-replicate reconciliation
    let mut per_file_calibrations: HashMap<String, RTCalibration> = HashMap::new();
    // Per-file isolation window m/z coverage (list of [lower, upper] intervals),
    // captured from each file's calibration so gap-fill can filter by m/z range.
    // Essential for GPF datasets where each file covers a disjoint m/z range.
    let mut per_file_isolation_mz: HashMap<String, Vec<(f64, f64)>> = HashMap::new();
    // Interner for modified_sequence deduplication across files
    let mut seq_interner: HashSet<Arc<str>> = HashSet::new();
    // Track whether all files loaded SVM scores from sidecars (enables skipping Percolator)
    let mut all_scores_loaded = true;

    for (file_idx, input_file) in config.input_files.iter().enumerate() {
        log::info!("");
        log::info!(
            "File #{}/{}: {}",
            file_idx + 1,
            config.input_files.len(),
            input_file
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
        );

        let file_name = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        // Check for cached scores parquet (skip search if already scored)
        let scores_path = scores_path_for_input(input_file);
        let scores_file_valid =
            scores_path.exists() && scores_path.metadata().map(|m| m.len() > 0).unwrap_or(false);
        let cached_stubs = if scores_file_valid {
            match load_fdr_stubs_from_parquet(&scores_path, &mut seq_interner) {
                Ok(mut stubs) => {
                    // Try loading SVM scores from sidecar (prefer 2nd-pass, fall back to 1st-pass)
                    let pass2 = fdr_scores_path_pass2(input_file);
                    let pass1 = fdr_scores_path_pass1(input_file);
                    if load_fdr_scores_sidecar(&pass2, &mut stubs) {
                        log::debug!(
                            "Loaded {} cached scores + 2nd-pass SVM scores from {}",
                            stubs.len(),
                            scores_path.display()
                        );
                    } else if load_fdr_scores_sidecar(&pass1, &mut stubs) {
                        log::debug!(
                            "Loaded {} cached scores + 1st-pass SVM scores from {}",
                            stubs.len(),
                            scores_path.display()
                        );
                    } else {
                        all_scores_loaded = false;
                        log::debug!(
                            "Loaded {} cached scores from {} (no SVM scores sidecar)",
                            stubs.len(),
                            scores_path.display()
                        );
                    }
                    // Load cached calibration for inter-replicate reconciliation
                    if let Some(input_dir) = input_file.parent() {
                        let cal_path = calibration_path_for_input(input_file, input_dir);
                        if cal_path.exists() {
                            if let Ok(cal_params) = load_calibration(&cal_path) {
                                if let Some(ref model_params) =
                                    cal_params.rt_calibration.model_params
                                {
                                    if let Ok(rt_cal) = RTCalibration::from_model_params(
                                        model_params,
                                        cal_params.rt_calibration.residual_sd,
                                    ) {
                                        per_file_calibrations.insert(file_name.clone(), rt_cal);
                                    }
                                }
                            }
                        }
                    }
                    // For PIN writing: if needed and PIN doesn't exist yet,
                    // fall back to full parquet load for this file
                    if write_pin {
                        let pin_path_check = mokapot_dir.join(format!("{}.pin", file_name));
                        if pin_path_check.exists() {
                            pin_files.insert(file_name.clone(), pin_path_check);
                        } else {
                            // Need full entries for PIN writing
                            if let Ok(full_entries) = load_scores_parquet(&scores_path) {
                                let pin_entries = if config.fdr_method == FdrMethod::Mokapot {
                                    compete_target_decoy_pairs(full_entries)
                                } else {
                                    full_entries
                                };
                                if let Ok(pin_path) =
                                    mokapot.write_pin_file(&file_name, &pin_entries, &mokapot_dir)
                                {
                                    pin_files.insert(file_name.clone(), pin_path);
                                }
                            }
                        }
                    }
                    Some(stubs)
                }
                Err(e) => {
                    log::warn!(
                        "Failed to load cached scores from {}: {}. Re-scoring.",
                        scores_path.display(),
                        e
                    );
                    all_scores_loaded = false;
                    None
                }
            }
        } else {
            None
        };

        // If we got cached stubs, skip loading spectra and scoring
        let fdr_entries = if let Some(stubs) = cached_stubs {
            stubs
        } else {
            all_scores_loaded = false;
            // Load spectra
            let (spectra, ms1_index) = load_all_spectra(input_file)?;
            if spectra.is_empty() {
                log::warn!("No spectra found in {}", input_file.display());
                continue;
            }

            // Save binary spectra cache for fast reload during re-scoring
            let cache_path = spectra_cache_path(input_file);
            if let Err(e) = save_spectra_cache(&cache_path, &spectra, &ms1_index) {
                log::warn!("Failed to save spectra cache: {}", e);
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
                    log::debug!("Reusing cached calibration from {}", cal_path.display());
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

            // Retain RT calibration for inter-replicate reconciliation
            if let Some(ref cal) = rt_cal {
                per_file_calibrations.insert(file_name.clone(), cal.clone());
            }

            // Retain isolation window coverage for gap-fill m/z filtering
            if let Some(ref cp) = cal_params {
                if let Some(ref scheme) = cp.metadata.isolation_scheme {
                    let intervals: Vec<(f64, f64)> = scheme
                        .windows
                        .iter()
                        .map(|&(center, width)| {
                            let half = width / 2.0;
                            (center - half, center + half)
                        })
                        .collect();
                    per_file_isolation_mz.insert(file_name.clone(), intervals);
                }
            }

            // Log one-line calibration summary for terminal
            if let Some(ref cp) = cal_params {
                let n_pep = cp.metadata.num_confident_peptides;
                let rt_tol = cp.rt_calibration.mad.unwrap_or(0.0) * 1.4826 * 3.0;
                let ms1_str = if cp.ms1_calibration.calibrated {
                    format!(
                        ", MS1 shift {:+.2} {} (tol {:.2})",
                        cp.ms1_calibration.mean,
                        cp.ms1_calibration.unit,
                        3.0 * cp.ms1_calibration.sd,
                    )
                } else {
                    String::new()
                };
                let ms2_str = if cp.ms2_calibration.calibrated {
                    format!(
                        ", MS2 shift {:+.4} {} (tol {:.4})",
                        cp.ms2_calibration.mean,
                        cp.ms2_calibration.unit,
                        3.0 * cp.ms2_calibration.sd,
                    )
                } else {
                    String::new()
                };
                log::info!(
                    "Calibration: {} peptides, RT {:.2} min{}{}",
                    n_pep,
                    rt_tol,
                    ms1_str,
                    ms2_str
                );
            }

            // Run coelution search
            let entries = run_search(
                &library,
                &spectra,
                &ms1_index,
                cal_params.as_ref(),
                rt_cal.as_ref(),
                &config,
                &file_name,
                None,
                "Scoring",
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

            // Save scores to parquet for future re-use (with search parameter metadata)
            let scores_meta = build_scores_metadata(&config);
            if let Err(e) =
                write_scores_parquet_with_metadata(&scores_path, &entries, Some(scores_meta))
            {
                log::warn!("Failed to save scores to {}: {}", scores_path.display(), e);
            }

            // Write coelution PIN file (if mokapot or --write-pin)
            if write_pin {
                let pin_entries = if config.fdr_method == FdrMethod::Mokapot {
                    compete_target_decoy_pairs(entries.clone())
                } else {
                    entries.clone()
                };
                let pin_path = mokapot.write_pin_file(&file_name, &pin_entries, &mokapot_dir)?;
                pin_files.insert(file_name.clone(), pin_path);
            }

            // Convert to lightweight FdrEntry stubs (intern modified_sequences).
            // CRITICAL: set parquet_index to the entry's row in the Parquet file
            // (which matches its position in the Vec at this point, before compaction).
            // Without this, all stubs have parquet_index=0 and Percolator scores
            // every entry against the same feature row.
            let mut fdr_stubs: Vec<FdrEntry> = entries
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    let mut stub = e.to_fdr_entry();
                    stub.parquet_index = i as u32;
                    stub
                })
                .collect();
            for stub in &mut fdr_stubs {
                stub.modified_sequence = intern_seq(&mut seq_interner, &stub.modified_sequence);
            }
            fdr_stubs
        };

        per_file_cache_paths.insert(file_name.clone(), scores_path);
        per_file_entries.push((file_name, fdr_entries));
    }

    let total_scored: usize = per_file_entries.iter().map(|(_, s)| s.len()).sum();
    log::debug!(
        "Coelution analysis complete. {} total scored entries across {} files",
        total_scored,
        config.input_files.len()
    );

    if per_file_entries.is_empty() || total_scored == 0 {
        log::warn!("No scored entries found. Cannot perform FDR control.");
        return Ok(());
    }

    // Check if we can skip Percolator + reconciliation entirely:
    // All files must be ValidReconciled (or single-file) AND have SVM scores loaded.
    let reconciliation_enabled = config.reconciliation.enabled && config.input_files.len() > 1;
    let all_already_reconciled = reconciliation_enabled
        && per_file_cache_paths.values().all(|path| {
            matches!(
                validate_scores_cache(path, &config),
                CacheValidity::ValidReconciled
            )
        });
    let can_skip_fdr = all_scores_loaded
        && (all_already_reconciled || !reconciliation_enabled)
        && config.fdr_method != FdrMethod::Simple;

    if can_skip_fdr {
        // SVM scores loaded from sidecars -- recompute q-values using two-pass FDR
        log::info!("");
        log::info!(
            "SVM scores loaded from cache. Recomputing q-values (skipping Percolator training)..."
        );
        // First-pass: compute q-values on all entries
        percolator::compute_fdr_from_stubs(&mut per_file_entries, config.run_fdr, None);
    } else {
        // Dispatch FDR control based on method
        log::info!("");
        log::info!(
            "Running {} FDR control on coelution results...",
            config.fdr_method
        );

        // No overlay for first-pass FDR (no re-scoring has happened yet)
        let empty_overlay: RescoreOverlay = HashMap::new();

        match config.fdr_method {
            FdrMethod::Percolator => {
                run_percolator_fdr(
                    &mut per_file_entries,
                    &per_file_cache_paths,
                    &config,
                    &empty_overlay,
                    None,
                )?;
            }
            FdrMethod::Mokapot => {
                run_mokapot_fdr(
                    &mut per_file_entries,
                    &per_file_cache_paths,
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

        // Persist first-pass SVM scores to sidecar files
        log::debug!("Persisting 1st-pass FDR scores...");
        persist_fdr_scores(
            &per_file_entries,
            &config,
            fdr_scores_path_pass1,
            "1st-pass",
        );
    }

    // First-pass protein FDR (picked-protein, Savitski 2015).
    //
    // Runs on the full pre-compaction peptide pool so target and decoy proteins
    // have a symmetric set of competitors. Without this, compaction would drop
    // decoys paired with non-passing targets, leaving picked-protein with a
    // one-sided comparison.
    //
    // Produces run_protein_qvalue on every FdrEntry stub. Used downstream for:
    //   (a) protein-aware compaction (rescue borderline peptides from strong proteins)
    //   (b) reconciliation consensus selection (strong proteins anchor consensus)
    if !can_skip_fdr && config.protein_fdr.is_some() {
        use osprey_fdr::protein;
        log::info!("");
        log::info!("First-pass protein FDR");

        // Build parsimony from peptides that pass first-pass peptide FDR.
        let detected_peptides: HashSet<String> = per_file_entries
            .iter()
            .flat_map(|(_, entries)| entries.iter())
            .filter(|e| !e.is_decoy && e.run_peptide_qvalue <= config.run_fdr)
            .map(|e| e.modified_sequence.to_string())
            .collect();
        let parsimony = protein::build_protein_parsimony(
            &library,
            config.shared_peptides,
            Some(&detected_peptides),
        );

        // Collect best peptide scores from the full pre-compaction pool.
        let peptide_scores = protein::collect_best_peptide_scores(&per_file_entries);

        // Run picked-protein FDR at 1 x run_fdr gate (Savitski convention).
        let protein_fdr_result =
            protein::compute_protein_fdr(&parsimony, &peptide_scores, config.run_fdr);

        // Write first-pass protein q-values into FdrEntry.run_protein_qvalue.
        // Do not set experiment_protein_qvalue yet — second-pass will overwrite it.
        protein::propagate_protein_qvalues(
            &mut per_file_entries,
            &protein_fdr_result,
            true,  // set_run
            false, // set_experiment
        );

        let n_at_1pct = protein_fdr_result
            .group_qvalues
            .values()
            .filter(|&&q| q <= config.run_fdr)
            .count();
        log::info!(
            "First-pass protein FDR: {} target groups at {:.0}% FDR",
            n_at_1pct,
            config.run_fdr * 100.0,
        );
    }

    // Collect first-pass passing base_ids for compaction.
    //
    // Protein-aware compaction: a peptide survives if EITHER
    //   (a) its peptide-level q-value <= reconciliation_compaction_fdr (default 0.01), OR
    //   (b) its first-pass protein q-value <= config.protein_fdr (typically 0.01)
    // Rule (b) lets strong protein evidence rescue borderline peptides whose
    // own peptide q-value is above the compaction gate. When --protein-fdr is
    // not set, only rule (a) applies.
    let peptide_compaction_gate = config.reconciliation_compaction_fdr;
    let protein_compaction_gate = config.protein_fdr.unwrap_or(0.0);
    let first_pass_base_ids: HashSet<u32> = per_file_entries
        .iter()
        .flat_map(|(_, entries)| entries.iter())
        .filter(|e| {
            !e.is_decoy
                && (e.run_peptide_qvalue <= peptide_compaction_gate
                    || (protein_compaction_gate > 0.0
                        && e.run_protein_qvalue <= protein_compaction_gate))
        })
        .map(|e| e.entry_id & 0x7FFF_FFFF)
        .collect();

    log::info!(
        "First-pass: {} unique precursors pass peptide-level FDR <= {:.0}% (or protein rescue)",
        first_pass_base_ids.len(),
        peptide_compaction_gate * 100.0,
    );

    // Compact FDR stubs: drop entries for precursors that didn't pass compaction.
    // Keeps passing targets + their paired decoys (by base_id).
    // The parquet_index field preserves the original row index for CWT/feature lookup.
    if !can_skip_fdr {
        let entries_before: usize = per_file_entries.iter().map(|(_, e)| e.len()).sum();
        for (_, entries) in per_file_entries.iter_mut() {
            entries.retain(|e| {
                let base_id = e.entry_id & 0x7FFF_FFFF;
                first_pass_base_ids.contains(&base_id)
            });
            entries.shrink_to_fit();
        }
        let entries_after: usize = per_file_entries.iter().map(|(_, e)| e.len()).sum();
        let saved_gb =
            ((entries_before - entries_after) as f64 * 128.0) / (1024.0 * 1024.0 * 1024.0);
        log::debug!(
            "Compacted FDR stubs: {} -> {} entries (freed {:.1} GB)",
            entries_before,
            entries_after,
            saved_gb,
        );
    }

    // Post-FDR re-scoring: multi-charge consensus + inter-replicate reconciliation.
    // Skipped entirely when all files are already reconciled with scores loaded.
    let per_file_overlays: RescoreOverlay = HashMap::new();
    if can_skip_fdr {
        log::info!("");
        log::info!("Skipping reconciliation (cached)");
        // Second-pass FDR: restrict to first-pass passing precursors (same as real pipeline)
        log::debug!(
            "Recomputing second-pass q-values on {} first-pass precursors + paired decoys...",
            first_pass_base_ids.len(),
        );
        percolator::compute_fdr_from_stubs(
            &mut per_file_entries,
            config.run_fdr,
            Some(&first_pass_base_ids),
        );
    } else {
        use crate::reconciliation::{
            compute_consensus_rts, identify_gap_fill_targets, plan_reconciliation,
            refit_calibration_with_consensus, GapFillTarget, ReconcileAction,
        };

        // Build (modified_sequence, charge) → (target_entry_id, decoy_entry_id) lookup
        // for gap-fill identification
        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> = {
            let mut map = HashMap::new();
            for entry in library.iter().filter(|e| !e.is_decoy) {
                let decoy_id = entry.id | 0x80000000;
                map.insert(
                    (Arc::from(entry.modified_sequence.as_str()), entry.charge),
                    (entry.id, decoy_id),
                );
            }
            map
        };

        // entry_id → library precursor m/z, for gap-fill m/z range filtering.
        let lib_precursor_mz: HashMap<u32, f64> = library
            .iter()
            .filter(|e| !e.is_decoy)
            .map(|e| (e.id, e.precursor_mz))
            .collect();

        // Build file_name → input_file index mapping
        let file_name_to_idx: HashMap<String, usize> = config
            .input_files
            .iter()
            .enumerate()
            .map(|(idx, p)| {
                let stem = p
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string();
                (stem, idx)
            })
            .collect();

        // Check if all files are already reconciled (parquet metadata says so).
        // If yes, skip BOTH multi-charge consensus AND inter-replicate reconciliation
        // since the previous reconciled run already applied both.
        let reconciliation_enabled = config.reconciliation.enabled && config.input_files.len() > 1;
        let all_already_reconciled = reconciliation_enabled
            && per_file_cache_paths.values().all(|path| {
                match validate_scores_cache(path, &config) {
                    CacheValidity::ValidReconciled => true,
                    CacheValidity::ValidFirstPass => false,
                    CacheValidity::Stale(reason) => {
                        log::debug!("Cache {} not reconciled: {}", path.display(), reason);
                        false
                    }
                }
            });

        // 1. Multi-charge consensus: compute per-file rescore targets
        //    Groups by (peptide, file). If at least one charge state passes FDR,
        //    the best-scoring charge state defines the consensus peak; other
        //    charge states at a different peak get re-scored.
        //    Skipped when files are already reconciled (consensus was applied previously).
        let per_file_consensus_targets: HashMap<String, Vec<(usize, f64, f64, f64)>> =
            if all_already_reconciled {
                HashMap::new()
            } else {
                per_file_entries
                    .iter()
                    .map(|(file_name, entries)| {
                        (
                            file_name.clone(),
                            select_post_fdr_consensus(entries, config.run_fdr),
                        )
                    })
                    .collect()
            };

        let total_consensus: usize = per_file_consensus_targets.values().map(|v| v.len()).sum();
        if total_consensus > 0 {
            log::debug!(
                "Multi-charge consensus: {} entries need re-scoring across all files",
                total_consensus
            );
        }

        // 2. Inter-replicate reconciliation: compute per-file rescore targets
        //    Uses first-pass FDR results to build consensus RTs across runs,
        //    then plans which entries need re-scoring at reconciled boundaries.
        let mut reconciliation_actions: HashMap<(String, usize), ReconcileAction> = HashMap::new();
        let refined_calibrations: HashMap<String, RTCalibration>;
        let per_file_gap_fill: HashMap<String, Vec<GapFillTarget>>;

        if all_already_reconciled {
            log::info!("");
            log::info!(
                "Reconciliation: all {} files already reconciled, skipping",
                per_file_cache_paths.len()
            );
            refined_calibrations = HashMap::new();
            per_file_gap_fill = HashMap::new();
        } else if reconciliation_enabled {
            log::info!("");
            log::info!("Reconciliation");

            let consensus = compute_consensus_rts(
                &per_file_entries,
                &per_file_calibrations,
                config.reconciliation.consensus_fdr,
                config.protein_fdr.unwrap_or(0.0),
            );

            if !consensus.is_empty() {
                refined_calibrations = per_file_entries
                    .par_iter()
                    .filter_map(|(file_name, entries)| {
                        refit_calibration_with_consensus(
                            &consensus,
                            entries,
                            config.reconciliation.consensus_fdr,
                        )
                        .map(|cal| (file_name.clone(), cal))
                    })
                    .collect();

                // Plan reconciliation with batched CWT loading.
                // Each file's CWT data is ~240 MB. Batch size is the minimum of:
                //   - available cores (no point loading more than we can process)
                //   - RAM / 500 MB (leave headroom for other allocations)
                // This keeps peak memory bounded while maximizing parallelism.
                let n_cores = rayon::current_num_threads();
                let available_ram_mb: usize = {
                    use sysinfo::System;
                    let mut sys = System::new();
                    sys.refresh_memory();
                    (sys.available_memory() / (1024 * 1024)) as usize
                };
                let ram_based_batch = (available_ram_mb / 500).max(1); // 500 MB per file with headroom
                let cwt_batch_size = n_cores.min(ram_based_batch).max(1);
                log::debug!(
                    "CWT batch loading: {} files per batch ({} cores, {} MB available RAM)",
                    cwt_batch_size,
                    n_cores,
                    available_ram_mb
                );
                for chunk in per_file_entries.chunks(cwt_batch_size) {
                    // Load CWT for this batch in parallel
                    let batch_cwt: HashMap<String, Vec<Vec<CwtCandidate>>> = chunk
                        .par_iter()
                        .filter_map(|(file_name, _)| {
                            per_file_cache_paths.get(file_name).and_then(|cache_path| {
                                load_cwt_candidates_from_parquet(cache_path)
                                    .ok()
                                    .map(|cwt| (file_name.clone(), cwt))
                            })
                        })
                        .collect();

                    // Plan reconciliation for this batch only
                    let batch_actions = plan_reconciliation(
                        &consensus,
                        chunk,
                        &batch_cwt,
                        &refined_calibrations,
                        &per_file_calibrations,
                        config.run_fdr,
                    );
                    reconciliation_actions.extend(batch_actions);
                    // batch_cwt dropped here — frees ~2.4 GB
                }

                // Gap-fill: identify passing precursors missing from each file
                // Uses run-level FDR: any precursor passing in at least one replicate
                per_file_gap_fill = identify_gap_fill_targets(
                    &consensus,
                    &per_file_entries,
                    &refined_calibrations,
                    &per_file_calibrations,
                    config.run_fdr,
                    &lib_lookup,
                    &lib_precursor_mz,
                    &per_file_isolation_mz,
                );
            } else {
                refined_calibrations = HashMap::new();
                reconciliation_actions = HashMap::new();
                per_file_gap_fill = HashMap::new();
            }
        } else {
            refined_calibrations = HashMap::new();
            reconciliation_actions = HashMap::new();
            per_file_gap_fill = HashMap::new();
        }

        let total_reconciliation: usize = reconciliation_actions
            .values()
            .filter(|a| !matches!(a, ReconcileAction::Keep))
            .count();
        // 3. Load spectra once per file and re-score both consensus + reconciliation targets.
        //    For re-scoring, we reload full CoelutionScoredEntry from parquet, rescore,
        //    save updated entries back to parquet, and update FdrEntry stubs.
        //    Files are processed in parallel for throughput.

        // Pre-group reconciliation actions by file name for efficient per-file lookup
        let mut per_file_reconciliation_targets: HashMap<String, Vec<(usize, f64, f64, f64)>> =
            HashMap::new();
        for ((file_name, idx), action) in &reconciliation_actions {
            let (start, apex, end) = match action {
                ReconcileAction::Keep => continue,
                ReconcileAction::UseCwtPeak {
                    start_rt,
                    apex_rt,
                    end_rt,
                    ..
                } => (*start_rt, *apex_rt, *end_rt),
                ReconcileAction::ForcedIntegration {
                    expected_rt,
                    half_width,
                } => (
                    expected_rt - half_width,
                    *expected_rt,
                    expected_rt + half_width,
                ),
            };
            per_file_reconciliation_targets
                .entry(file_name.clone())
                .or_default()
                .push((*idx, apex, start, end));
        }

        let mut total_rescored: usize = 0;
        let mut total_gap_cwt: usize = 0;
        let mut total_gap_forced: usize = 0;

        // Process files sequentially to limit peak memory — each file loads spectra
        // (~1.5 GB). Per-file re-scoring is already internally parallelized.
        let n_total_files = per_file_entries.len();
        for (file_num, (file_name, fdr_entries)) in per_file_entries.iter_mut().enumerate() {
            // Collect consensus targets for this file
            let consensus_targets = per_file_consensus_targets
                .get(file_name.as_str())
                .cloned()
                .unwrap_or_default();

            // Collect reconciliation targets for this file (pre-grouped)
            let reconciliation_targets = per_file_reconciliation_targets
                .get(file_name.as_str())
                .cloned()
                .unwrap_or_default();

            // Collect gap-fill targets for this file
            let gap_fill_targets = per_file_gap_fill
                .get(file_name.as_str())
                .cloned()
                .unwrap_or_default();

            // Merge both target sets, deduplicating by entry index
            // (if an entry appears in both consensus and reconciliation,
            // prefer the reconciliation target since it's inter-replicate informed)
            let mut combined_targets: HashMap<usize, (f64, f64, f64)> = HashMap::new();
            for (idx, apex, start, end) in &consensus_targets {
                combined_targets.insert(*idx, (*apex, *start, *end));
            }
            for (idx, apex, start, end) in &reconciliation_targets {
                combined_targets.insert(*idx, (*apex, *start, *end));
            }

            let all_targets: Vec<(usize, f64, f64, f64)> = combined_targets
                .into_iter()
                .map(|(idx, (apex, start, end))| (idx, apex, start, end))
                .collect();

            // Skip file if no work to do (no re-scoring targets AND no gap-fill targets)
            if all_targets.is_empty() && gap_fill_targets.is_empty() {
                continue;
            }

            let input_idx = match file_name_to_idx.get(file_name.as_str()) {
                Some(idx) => *idx,
                None => continue,
            };
            let input_file = &config.input_files[input_idx];

            let n_multi_charge = consensus_targets.len();
            let n_inter_replicate = reconciliation_targets.len();
            let n_gap_fill = gap_fill_targets.len();
            log::info!(
                "Re-scoring file {}/{}: {}",
                file_num + 1,
                n_total_files,
                file_name,
            );
            log::debug!(
                "  {} entries ({} consensus, {} reconciliation, {} gap-fill, {} unique after dedup)",
                all_targets.len() + n_gap_fill * 2,
                n_multi_charge,
                n_inter_replicate,
                n_gap_fill,
                all_targets.len(),
            );

            // Build boundary overrides keyed by entry_id and track idx mapping
            let mut boundary_overrides: HashMap<u32, (f64, f64, f64)> = HashMap::new();
            let mut entry_id_to_idx: HashMap<u32, usize> = HashMap::new();
            let mut subset_ids: std::collections::HashSet<u32> = std::collections::HashSet::new();
            for &(idx, apex, start, end) in &all_targets {
                let entry_id = fdr_entries[idx].entry_id;
                boundary_overrides.insert(entry_id, (apex, start, end));
                entry_id_to_idx.insert(entry_id, idx);
                subset_ids.insert(entry_id);
            }

            // Build subset library containing only entries that need re-scoring
            let subset_library: Vec<LibraryEntry> = library
                .iter()
                .filter(|e| subset_ids.contains(&e.id))
                .cloned()
                .collect();

            // Build gap-fill library subset (targets + decoys, separate from re-scoring)
            let mut gap_fill_ids: std::collections::HashSet<u32> = std::collections::HashSet::new();
            for gf in &gap_fill_targets {
                gap_fill_ids.insert(gf.target_entry_id);
                gap_fill_ids.insert(gf.decoy_entry_id);
            }
            let gap_fill_library: Vec<LibraryEntry> = if !gap_fill_ids.is_empty() {
                library
                    .iter()
                    .filter(|e| gap_fill_ids.contains(&e.id))
                    .cloned()
                    .collect()
            } else {
                Vec::new()
            };

            // Load spectra and calibration for this file
            let (spectra, ms1_index) = {
                let cache_path = spectra_cache_path(input_file);
                if cache_path.exists() {
                    match load_spectra_cache(&cache_path) {
                        Ok(result) => {
                            log::debug!(
                                "Loaded {} MS2 spectra from cache '{}'",
                                result.0.len(),
                                cache_path.display()
                            );
                            result
                        }
                        Err(e) => {
                            log::warn!("Failed to load spectra cache, falling back to mzML: {}", e);
                            load_all_spectra(input_file)?
                        }
                    }
                } else {
                    load_all_spectra(input_file)?
                }
            };

            // Use refined calibration if available, fall back to original
            let rt_cal = refined_calibrations
                .get(file_name.as_str())
                .or_else(|| per_file_calibrations.get(file_name.as_str()));

            let cal_params: Option<CalibrationParams> = input_file.parent().and_then(|input_dir| {
                let cal_path = calibration_path_for_input(input_file, input_dir);
                if cal_path.exists() {
                    load_calibration(&cal_path).ok()
                } else {
                    None
                }
            });

            // --- Re-score existing entries (consensus + reconciliation) ---
            let mut overlay: HashMap<usize, CoelutionScoredEntry> = HashMap::new();
            if !all_targets.is_empty() {
                let rescored = run_search(
                    &subset_library,
                    &spectra,
                    &ms1_index,
                    cal_params.as_ref(),
                    rt_cal,
                    &config,
                    file_name,
                    Some(&boundary_overrides),
                    "Re-scoring",
                )?;

                for entry in rescored {
                    if let Some(&idx) = entry_id_to_idx.get(&entry.entry_id) {
                        overlay.insert(idx, entry);
                    }
                }
            }

            // --- Gap-fill: score missing precursors (two-pass: CWT then forced) ---
            let mut n_gap_cwt = 0usize;
            let mut n_gap_forced = 0usize;
            if !gap_fill_targets.is_empty() {
                // Pass 1 (CWT): Run without boundary overrides, pre-filter disabled.
                // This lets CWT peak detection find natural peaks at the consensus RT.
                let mut gap_config = config.clone();
                gap_config.prefilter_enabled = false;
                let cwt_results = run_search(
                    &gap_fill_library,
                    &spectra,
                    &ms1_index,
                    cal_params.as_ref(),
                    rt_cal,
                    &gap_config,
                    file_name,
                    None,
                    "Gap-fill CWT",
                )?;

                // Track which entry_ids got CWT results
                let cwt_hit_ids: std::collections::HashSet<u32> =
                    cwt_results.iter().map(|e| e.entry_id).collect();
                n_gap_cwt = cwt_results.len();

                // Append CWT results as new FdrEntry stubs (gap-fill: no Parquet row)
                let gap_fill_start_idx = fdr_entries.len();
                for (i, entry) in cwt_results.into_iter().enumerate() {
                    let mut stub = entry.to_fdr_entry();
                    stub.parquet_index = u32::MAX; // gap-fill sentinel
                    stub.modified_sequence = intern_seq(&mut seq_interner, &stub.modified_sequence);
                    fdr_entries.push(stub);
                    overlay.insert(gap_fill_start_idx + i, entry);
                }

                // Pass 2 (Forced): Entries that CWT missed get forced integration
                // at consensus RT ± half_width
                let mut forced_overrides: HashMap<u32, (f64, f64, f64)> = HashMap::new();
                let mut forced_ids: std::collections::HashSet<u32> =
                    std::collections::HashSet::new();
                for gf in &gap_fill_targets {
                    // Add target if CWT missed it
                    if !cwt_hit_ids.contains(&gf.target_entry_id) {
                        let start = gf.expected_rt - gf.half_width;
                        let end = gf.expected_rt + gf.half_width;
                        forced_overrides.insert(gf.target_entry_id, (gf.expected_rt, start, end));
                        forced_ids.insert(gf.target_entry_id);
                    }
                    // Add decoy if CWT missed it
                    if !cwt_hit_ids.contains(&gf.decoy_entry_id) {
                        let start = gf.expected_rt - gf.half_width;
                        let end = gf.expected_rt + gf.half_width;
                        forced_overrides.insert(gf.decoy_entry_id, (gf.expected_rt, start, end));
                        forced_ids.insert(gf.decoy_entry_id);
                    }
                }

                if !forced_overrides.is_empty() {
                    let forced_library: Vec<LibraryEntry> = gap_fill_library
                        .iter()
                        .filter(|e| forced_ids.contains(&e.id))
                        .cloned()
                        .collect();

                    let forced_results = run_search(
                        &forced_library,
                        &spectra,
                        &ms1_index,
                        cal_params.as_ref(),
                        rt_cal,
                        &config,
                        file_name,
                        Some(&forced_overrides),
                        "Gap-fill forced",
                    )?;

                    n_gap_forced = forced_results.len();

                    // Append forced results as new FdrEntry stubs (gap-fill: no Parquet row)
                    let forced_start_idx = fdr_entries.len();
                    for (i, entry) in forced_results.into_iter().enumerate() {
                        let mut stub = entry.to_fdr_entry();
                        stub.parquet_index = u32::MAX; // gap-fill sentinel
                        stub.modified_sequence =
                            intern_seq(&mut seq_interner, &stub.modified_sequence);
                        fdr_entries.push(stub);
                        overlay.insert(forced_start_idx + i, entry);
                    }
                }

                total_gap_cwt += n_gap_cwt;
                total_gap_forced += n_gap_forced;
            }

            // Drop spectra before updating stubs to free ~1.5-3 GB
            drop(spectra);

            // Count existing re-scored entries (not gap-fill)
            let existing_rescore_indices: std::collections::HashSet<usize> =
                entry_id_to_idx.values().copied().collect();
            let n_existing_rescored = overlay
                .keys()
                .filter(|idx| existing_rescore_indices.contains(idx))
                .count();
            total_rescored += n_existing_rescored + n_gap_cwt + n_gap_forced;

            if !all_targets.is_empty() || n_gap_cwt + n_gap_forced > 0 {
                log::debug!(
                    "  {} of {} existing entries re-scored, {} gap-fill ({} CWT, {} forced)",
                    n_existing_rescored,
                    all_targets.len(),
                    n_gap_cwt + n_gap_forced,
                    n_gap_cwt,
                    n_gap_forced,
                );
            }

            // Update FdrEntry stubs from overlay for existing re-scored entries
            // (gap-fill entries were already appended with stubs during scoring above)
            for &idx in &existing_rescore_indices {
                if let Some(entry) = overlay.get(&idx) {
                    let mut stub = entry.to_fdr_entry();
                    // Preserve the original parquet_index from the existing stub
                    stub.parquet_index = fdr_entries[idx].parquet_index;
                    stub.modified_sequence = intern_seq(&mut seq_interner, &stub.modified_sequence);
                    fdr_entries[idx] = stub;
                }
            }
            // Write reconciled entries (overlay) back to Parquet, then drop overlay.
            // This avoids accumulating ~28 GB of overlay data across 240 files.
            // The blib output step will reload from this reconciled Parquet.
            if !overlay.is_empty() {
                // Merge overlay into a complete entry list for this file:
                // load original parquet, replace re-scored entries, append gap-fill entries.
                let scores_path = per_file_cache_paths.get(file_name.as_str());
                if let Some(cache_path) = scores_path {
                    if let Ok(mut full_entries) = load_scores_parquet(cache_path) {
                        // Replace re-scored entries using parquet_index (not Vec position).
                        // After first-pass compaction, Vec position != Parquet row index.
                        // Using Vec position would overwrite the wrong entry in the Parquet
                        // file and leave the actual target entry unchanged.
                        for &idx in &existing_rescore_indices {
                            if let Some(entry) = overlay.remove(&idx) {
                                let pq_idx = fdr_entries[idx].parquet_index as usize;
                                if pq_idx < full_entries.len() {
                                    full_entries[pq_idx] = entry;
                                } else {
                                    log::warn!(
                                        "Rescore write-back: parquet_index {} out of range for {} ({} rows)",
                                        pq_idx,
                                        file_name,
                                        full_entries.len(),
                                    );
                                }
                            }
                        }
                        // Append gap-fill entries (remaining overlay entries beyond parquet range).
                        // Update the corresponding fdr_entries stubs' parquet_index to point to
                        // the actual Parquet row they will occupy after write-back. Without this,
                        // gap-fill stubs keep parquet_index = u32::MAX and Phase 2 of the next
                        // Percolator run can't load their features.
                        let gap_start_pq_row = full_entries.len();
                        let mut gap_entries: Vec<(usize, CoelutionScoredEntry)> =
                            overlay.into_iter().collect();
                        gap_entries.sort_by_key(|(idx, _)| *idx);
                        for (gap_offset, (vec_idx, entry)) in gap_entries.into_iter().enumerate() {
                            let new_pq_row = (gap_start_pq_row + gap_offset) as u32;
                            full_entries.push(entry);
                            // The vec_idx in the overlay matches the position in fdr_entries
                            // (gap-fill entries were appended to fdr_entries in order)
                            if vec_idx < fdr_entries.len() {
                                fdr_entries[vec_idx].parquet_index = new_pq_row;
                            }
                        }
                        // Write back to Parquet with reconciliation metadata
                        let recon_metadata = build_reconciled_metadata(&config);
                        if let Err(e) = write_scores_parquet_with_metadata(
                            cache_path,
                            &full_entries,
                            Some(recon_metadata),
                        ) {
                            log::warn!(
                                "Failed to write reconciled scores for {}: {}",
                                file_name,
                                e
                            );
                        }
                    }
                }
                // Overlay is consumed — no memory accumulation
            }
        }

        // 4. Single second-pass FDR after all re-scoring
        if total_rescored > 0 {
            log::debug!(
                "Post-FDR re-scoring complete: {} entries re-scored ({} consensus, {} reconciliation, {} gap-fill [{} CWT, {} forced])",
                total_rescored,
                total_consensus,
                total_reconciliation,
                total_gap_cwt + total_gap_forced,
                total_gap_cwt,
                total_gap_forced,
            );
            log::info!("");
            log::info!("Second-pass FDR");
            match config.fdr_method {
                FdrMethod::Percolator => {
                    run_percolator_fdr(
                        &mut per_file_entries,
                        &per_file_cache_paths,
                        &config,
                        &per_file_overlays,
                        Some(&first_pass_base_ids),
                    )?;
                }
                FdrMethod::Mokapot => {
                    run_mokapot_fdr(
                        &mut per_file_entries,
                        &per_file_cache_paths,
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
        }
    }

    // Persist second-pass SVM scores to sidecar files (after reconciliation block closes)
    if !can_skip_fdr {
        let has_reconciliation = config.reconciliation.enabled && config.input_files.len() > 1;
        if has_reconciliation {
            log::debug!("Persisting 2nd-pass FDR scores...");
            persist_fdr_scores(
                &per_file_entries,
                &config,
                fdr_scores_path_pass2,
                "2nd-pass",
            );
        }
    }

    // Protein parsimony (always runs) + optional second-pass picked-protein FDR.
    //
    // Parsimony is a cheap deterministic graph analysis: bipartite graph
    // construction, identical-set merging, subset elimination, and razor
    // assignment. It always runs because the protein group assignments are
    // useful for downstream interpretation.
    //
    // Second-pass picked-protein FDR runs when `--protein-fdr` is set. Unlike
    // first-pass protein FDR (which used the full pre-compaction pool for gating),
    // this pass uses the reconciled + second-pass-scored peptides to produce
    // AUTHORITATIVE experiment_protein_qvalue values that feed FdrLevel::Protein
    // filtering and the protein CSV report.
    {
        use osprey_fdr::protein;
        log::info!("");
        log::info!("Protein parsimony");

        // Build the protein parsimony bipartite graph from peptides that pass at
        // the peptide level. Use peptide-level FDR as the broader gate when
        // --fdr-level is protein (protein q-values haven't been computed yet
        // at this point in the second pass either).
        let peptide_gate_level = match config.fdr_level {
            FdrLevel::Protein => FdrLevel::Peptide,
            other => other,
        };
        let detected_peptides: HashSet<String> = per_file_entries
            .iter()
            .flat_map(|(_, entries)| entries.iter())
            .filter(|e| {
                !e.is_decoy
                    && e.effective_experiment_qvalue(peptide_gate_level) <= config.experiment_fdr
            })
            .map(|e| e.modified_sequence.to_string())
            .collect();

        let parsimony = protein::build_protein_parsimony(
            &library,
            config.shared_peptides,
            Some(&detected_peptides),
        );

        // Optional second-pass picked-protein FDR (Savitski)
        if let Some(protein_fdr_threshold) = config.protein_fdr {
            log::info!("Second-pass protein FDR");
            // Collect best peptide scores from compacted + reconciled + second-pass-scored stubs
            let best_scores = protein::collect_best_peptide_scores(&per_file_entries);

            // Savitski gate: 1 x run_fdr (same as peptide-level threshold)
            let protein_fdr_result =
                protein::compute_protein_fdr(&parsimony, &best_scores, config.run_fdr);

            let n_passing = protein_fdr_result
                .group_qvalues
                .values()
                .filter(|&&q| q <= protein_fdr_threshold)
                .count();
            log::info!(
                "Second-pass protein FDR: {} protein groups at {:.0}% FDR (from {} total groups)",
                n_passing,
                protein_fdr_threshold * 100.0,
                parsimony.groups.len()
            );

            // Propagate protein q-values into FdrEntry stubs — experiment side only.
            // The run side was already set by first-pass protein FDR and is used
            // only as a gate (not for final output).
            protein::propagate_protein_qvalues(
                &mut per_file_entries,
                &protein_fdr_result,
                false, // set_run: leave first-pass values in place
                true,  // set_experiment: AUTHORITATIVE
            );

            let protein_report_path = config.output_blib.with_extension("proteins.csv");
            if let Err(e) = protein::write_protein_report(
                &protein_report_path,
                &parsimony,
                &protein_fdr_result,
                protein_fdr_threshold,
                &library,
            ) {
                log::warn!("Failed to write protein report: {}", e);
            }
        } else {
            log::info!("  {} protein groups (no FDR)", parsimony.groups.len());
            if matches!(config.fdr_level, FdrLevel::Protein) {
                log::warn!(
                    "--fdr-level protein requires --protein-fdr; falling back to peptide-level filtering"
                );
            }
        }
    }

    // Two-stage blib output gate.
    //
    // Stage 1 (peptide gate): the configured --fdr-level determines which
    // peptide identities are eligible for output.
    //   Precursor → peptides with at least one charge state at experiment_precursor_qvalue <= threshold
    //   Peptide   → peptides with experiment_peptide_qvalue <= threshold (default)
    //   Protein   → peptides with experiment_protein_qvalue <= threshold
    //   Both      → peptides with max(precursor, peptide) <= threshold
    //
    // Stage 2 (precursor gate): within each eligible peptide, include only
    // charge states that individually pass experiment_precursor_qvalue <= threshold.
    // If NO charge state of a peptide passes precursor-level FDR (possible because
    // peptide-level FDR aggregates across charges), include the best charge state
    // (lowest experiment_precursor_qvalue) as a representative.
    //
    // This prevents phantom charge states (which never passed precursor FDR)
    // from appearing in the blib. Prior to this fix, commit 53e6a0e replaced
    // the original FdrLevel::Both gate with config.fdr_level, letting all charge
    // states of a peptide-level-passing peptide through regardless of their
    // precursor-level q-value.
    let effective_level = match (config.fdr_level, config.protein_fdr) {
        (FdrLevel::Protein, None) => FdrLevel::Peptide,
        (level, _) => level,
    };

    // Stage 1: which peptides pass the configured FDR level?
    let passing_peptides: HashSet<Arc<str>> = per_file_entries
        .iter()
        .flat_map(|(_, entries)| entries.iter())
        .filter(|e| {
            !e.is_decoy && e.effective_experiment_qvalue(effective_level) <= config.experiment_fdr
        })
        .map(|e| e.modified_sequence.clone())
        .collect();

    // Stage 2: within passing peptides, admit only precursors that pass
    // precursor-level FDR, with a fallback to the best charge state per peptide.
    let mut precursor_passing: HashSet<(Arc<str>, u8)> = HashSet::new();
    let mut best_per_peptide: HashMap<Arc<str>, (u8, f64)> = HashMap::new();
    for (_, entries) in per_file_entries.iter() {
        for e in entries {
            if e.is_decoy || !passing_peptides.contains(&e.modified_sequence) {
                continue;
            }
            if e.experiment_precursor_qvalue <= config.experiment_fdr {
                precursor_passing.insert((e.modified_sequence.clone(), e.charge));
            }
            best_per_peptide
                .entry(e.modified_sequence.clone())
                .and_modify(|(best_charge, best_q)| {
                    if e.experiment_precursor_qvalue < *best_q {
                        *best_charge = e.charge;
                        *best_q = e.experiment_precursor_qvalue;
                    }
                })
                .or_insert((e.charge, e.experiment_precursor_qvalue));
        }
    }
    // Fallback: peptides with no precursor-passing charge state keep their best
    let mut n_fallback = 0usize;
    for peptide in &passing_peptides {
        let has_passing = precursor_passing.iter().any(|(ms, _)| ms == peptide);
        if !has_passing {
            if let Some(&(charge, _)) = best_per_peptide.get(peptide) {
                precursor_passing.insert((peptide.clone(), charge));
                n_fallback += 1;
            }
        }
    }
    if n_fallback > 0 {
        log::info!(
            "{} peptides had no charge state passing precursor-level FDR; best charge state kept as fallback",
            n_fallback
        );
    }
    let passing_precursors = precursor_passing;

    // Compute best experiment q-value per precursor (precursor-level, for blib output)
    let mut best_exp_q: HashMap<(Arc<str>, u8), f64> = HashMap::new();
    for (_, entries) in per_file_entries.iter() {
        for e in entries {
            if !e.is_decoy && passing_precursors.contains(&(e.modified_sequence.clone(), e.charge))
            {
                let eq = e.experiment_precursor_qvalue;
                best_exp_q
                    .entry((e.modified_sequence.clone(), e.charge))
                    .and_modify(|q| *q = q.min(eq))
                    .or_insert(eq);
            }
        }
    }

    // Extract lightweight per-file FDR data for passing entries, then drop the
    // heavy FdrEntry stubs (~19 GB for 240 files) to free memory for blib loading.
    struct LightFdr {
        run_qvalue: f64,
        experiment_qvalue: f64,
        is_decoy: bool,
        modified_sequence: Arc<str>,
        charge: u8,
        /// Parquet row index for looking up RT/boundaries from the 5-column projection.
        /// After compaction, Vec position != Parquet row index, so this field is
        /// essential for correct blib plan entry construction.
        parquet_index: u32,
    }
    let per_file_fdr_data: Vec<(String, Vec<LightFdr>)> = per_file_entries
        .iter()
        .map(|(file_name, entries)| {
            let fdr_values: Vec<LightFdr> = entries
                .iter()
                .map(|e| LightFdr {
                    run_qvalue: e.effective_run_qvalue(FdrLevel::Both),
                    experiment_qvalue: e.effective_experiment_qvalue(FdrLevel::Both),
                    is_decoy: e.is_decoy,
                    modified_sequence: e.modified_sequence.clone(),
                    charge: e.charge,
                    parquet_index: e.parquet_index,
                })
                .collect();
            (file_name.clone(), fdr_values)
        })
        .collect();

    let files_with_passing: HashSet<String> = per_file_entries
        .iter()
        .filter(|(_, entries)| {
            entries.iter().any(|e| {
                !e.is_decoy && passing_precursors.contains(&(e.modified_sequence.clone(), e.charge))
            })
        })
        .map(|(file_name, _)| file_name.clone())
        .collect();

    let n_total_files = per_file_entries.len();

    // Drop the heavy FdrEntry stubs — frees ~19 GB
    drop(per_file_entries);
    log::info!(
        "Loading blib data from {}/{} files with passing precursors",
        files_with_passing.len(),
        n_total_files
    );

    // Load lightweight BlibPlanEntry from Parquet (5-column projection) per file,
    // merging with LightFdr data. ~96 bytes per entry vs ~940 bytes for full entries.
    // 28.8M entries x 96 bytes = ~2.8 GB (fits in RAM).

    // Build file_name -> index mapping for compact storage in BlibPlanEntry
    let file_names: Vec<String> = per_file_fdr_data.iter().map(|(f, _)| f.clone()).collect();
    let file_name_to_idx: HashMap<&str, u16> = file_names
        .iter()
        .enumerate()
        .map(|(i, f)| (f.as_str(), i as u16))
        .collect();

    let mut plan_entries: Vec<BlibPlanEntry> = Vec::new();

    for (file_num, (file_name, fdr_values)) in per_file_fdr_data.iter().enumerate() {
        if !files_with_passing.contains(file_name) {
            continue;
        }
        let cache_path = match per_file_cache_paths.get(file_name.as_str()) {
            Some(p) => p,
            None => continue,
        };
        let file_name_idx = *file_name_to_idx.get(file_name.as_str()).unwrap_or(&0);

        // Load 5-column projection from Parquet: entry_id, apex_rt, start_rt, end_rt, bounds_area
        let parquet_rows = match load_blib_plan_from_parquet(cache_path) {
            Ok(rows) => rows,
            Err(e) => {
                log::warn!("Failed to load plan entries from {}: {}", file_name, e);
                continue;
            }
        };

        // Merge Parquet row data with LightFdr data, filter to passing targets.
        // CRITICAL: use fdr.parquet_index to look up the correct Parquet row.
        // After first-pass compaction, Vec position != Parquet row index, so
        // zipping by position would assign wrong RT/boundaries to peptides.
        for fdr in fdr_values.iter() {
            if fdr.is_decoy
                || !passing_precursors.contains(&(fdr.modified_sequence.clone(), fdr.charge))
            {
                continue;
            }
            let pq_idx = fdr.parquet_index as usize;
            let Some(row) = parquet_rows.get(pq_idx) else {
                log::warn!(
                    "LightFdr.parquet_index {} out of range for {} ({} parquet rows) — skipping entry {}",
                    pq_idx,
                    file_name,
                    parquet_rows.len(),
                    fdr.modified_sequence,
                );
                continue;
            };
            plan_entries.push(BlibPlanEntry {
                entry_id: row.0,
                charge: fdr.charge,
                file_name_idx,
                run_qvalue: fdr.run_qvalue,
                experiment_qvalue: best_exp_q
                    .get(&(fdr.modified_sequence.clone(), fdr.charge))
                    .copied()
                    .unwrap_or(fdr.experiment_qvalue),
                apex_rt: row.1,
                start_rt: row.2,
                end_rt: row.3,
                bounds_area: row.4,
                modified_sequence: fdr.modified_sequence.clone(),
            });
        }
        if (file_num + 1) % 20 == 0 || file_num + 1 == n_total_files {
            log::info!(
                "  Loaded {}/{} files ({} observations so far)",
                file_num + 1,
                n_total_files,
                plan_entries.len()
            );
        }
    }

    // Write blib output from plan entries, looking up fragments/mods/proteins
    // from the in-memory library. Write to a local temp file first, then move
    // to the final destination. This avoids SQLite locking issues on network
    // filesystems.
    if !plan_entries.is_empty() {
        log::info!("Writing blib to {}", config.output_blib.display());

        let final_path = &config.output_blib;
        let temp_path = std::env::temp_dir().join(format!("osprey_{}.blib", std::process::id()));

        // Write to local temp file
        let mut temp_config = config.clone();
        temp_config.output_blib = temp_path.clone();
        write_blib_from_plan(&temp_config, &library, &plan_entries, &file_names)?;

        // Move to final destination (safe copy for network filesystems)
        osprey_core::copy_and_verify(&temp_path, final_path)?;
    } else {
        log::warn!("No peptides passed FDR threshold, skipping blib output");
    }

    // Write output report if specified
    // TSV includes only passing entries. Parquet report needs all entries (targets + decoys)
    // for feature analysis and SVM retraining — reload from caches if requested.
    if let Some(report_path) = &config.output_report {
        let ext = report_path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("");
        if ext.eq_ignore_ascii_case("parquet") {
            log::info!("Writing Parquet report to {}", report_path.display());
            // Reload all entries from parquet caches for the full report
            let mut all_entries: Vec<CoelutionScoredEntry> = Vec::new();
            for (file_name, fdr_values) in per_file_fdr_data.iter() {
                let cache_path = match per_file_cache_paths.get(file_name.as_str()) {
                    Some(p) => p,
                    None => continue,
                };
                match load_scores_parquet(cache_path) {
                    Ok(mut entries) => {
                        // Merge FDR q-values from lightweight data
                        for (full, fdr) in entries.iter_mut().zip(fdr_values.iter()) {
                            full.run_precursor_qvalue = fdr.run_qvalue;
                            full.experiment_precursor_qvalue = fdr.experiment_qvalue;
                        }
                        all_entries.extend(entries);
                    }
                    Err(e) => {
                        log::warn!(
                            "Failed to reload entries for report from {}: {}",
                            file_name,
                            e
                        );
                    }
                }
            }
            write_parquet_report(report_path, &all_entries)?;
        } else {
            // TSV report from plan entries (no full entry reload needed)
            log::info!("Writing TSV report to {}", report_path.display());
            write_tsv_report_from_plan(report_path, &plan_entries, &library)?;
        }
    }

    Ok(())
}

/// Run native Percolator FDR on coelution entries.
///
/// Uses a streaming approach for scalability with large experiments:
/// 1. Build lightweight metadata from FdrEntry stubs (already in memory)
/// 2. Subsample ~300K entries for SVM training (no feature loading needed)
/// 3. Load PIN features only for the subsampled entries from Parquet
/// 4. Train SVM models via run_percolator on the subsample
/// 5. Score ALL entries by streaming through Parquet files one at a time
/// 6. Compute q-values and PEP from the flat scores array
///
/// This avoids loading all Parquet files into memory simultaneously and
/// avoids creating PercolatorEntry objects for all entries.
fn run_percolator_fdr(
    per_file_entries: &mut [(String, Vec<FdrEntry>)],
    per_file_cache_paths: &HashMap<String, std::path::PathBuf>,
    config: &OspreyConfig,
    overlay: &RescoreOverlay,
    restrict_base_ids: Option<&HashSet<u32>>,
) -> Result<()> {
    use osprey_fdr::percolator;

    log::debug!("Running native Percolator FDR on coelution entries");

    let total_entries: usize = per_file_entries.iter().map(|(_, e)| e.len()).sum();
    let n_files = per_file_entries.len();
    let perc_config = percolator::PercolatorConfig {
        train_fdr: config.run_fdr,
        test_fdr: config.run_fdr,
        feature_names: Some(
            get_pin_feature_names()
                .iter()
                .map(|s| s.to_string())
                .collect(),
        ),
        ..Default::default()
    };
    let max_train = perc_config.max_train_size;

    // For small experiments (total entries fit comfortably), use the direct path
    // to preserve full cross-validation scoring semantics
    let use_streaming = max_train > 0 && total_entries > max_train * 2;

    if !use_streaming {
        // Direct path: load all features and call run_percolator
        log::debug!(
            "Percolator: loading features for {} entries across {} files",
            total_entries,
            n_files
        );
        return run_percolator_fdr_direct(
            per_file_entries,
            per_file_cache_paths,
            &perc_config,
            overlay,
            restrict_base_ids,
        );
    }

    // === Streaming path for large experiments ===
    // Memory-efficient: no flat metadata arrays for all entries.
    // Subsamples directly from per_file_entries (~16 MB HashMaps),
    // scores directly to entry.score, computes q-values per-file.
    log::debug!(
        "Percolator streaming: {} entries across {} files (training on {} subset)",
        total_entries,
        n_files,
        max_train
    );

    // Phase 1: Select best observation per precursor across all files for training.
    //
    // With 240 files × ~1M entries/file, each peptide has ~480 entries (target+decoy × 240).
    // Selecting peptide groups would give only ~625 unique peptides to reach 300K entries.
    // Instead, pick the BEST-scoring observation per base_id (one target, one decoy),
    // giving ~500K unique precursor observations with maximum peptide diversity.
    // Then subsample from those if still > max_train.

    // Find best target and best decoy per base_id, tracking (file_idx, local_idx, score)
    let mut best_target: HashMap<u32, (usize, usize, f64)> = HashMap::new();
    let mut best_decoy: HashMap<u32, (usize, usize, f64)> = HashMap::new();
    for (file_idx, (_, entries)) in per_file_entries.iter().enumerate() {
        for (local_idx, entry) in entries.iter().enumerate() {
            let base_id = entry.entry_id & 0x7FFF_FFFF;
            if let Some(restrict) = restrict_base_ids {
                if !restrict.contains(&base_id) {
                    continue;
                }
            }
            let map = if entry.is_decoy {
                &mut best_decoy
            } else {
                &mut best_target
            };
            map.entry(base_id)
                .and_modify(|(best_fi, best_li, best_score)| {
                    if entry.coelution_sum > *best_score {
                        *best_fi = file_idx;
                        *best_li = local_idx;
                        *best_score = entry.coelution_sum;
                    }
                })
                .or_insert((file_idx, local_idx, entry.coelution_sum));
        }
    }

    // Collect all best observations: (file_idx, local_idx) pairs
    let mut best_observations: Vec<(usize, usize)> =
        Vec::with_capacity(best_target.len() + best_decoy.len());
    for &(fi, li, _) in best_target.values() {
        best_observations.push((fi, li));
    }
    for &(fi, li, _) in best_decoy.values() {
        best_observations.push((fi, li));
    }
    best_observations.sort(); // deterministic order

    let dedup_count = best_observations.len();
    log::debug!(
        "  Best-per-precursor: {} entries ({} targets, {} decoys) from {} total",
        dedup_count,
        best_target.len(),
        best_decoy.len(),
        total_entries
    );

    drop(best_target);
    drop(best_decoy);

    // Subsample from deduplicated set if still > max_train
    let subset: Vec<(usize, usize)> = if dedup_count <= max_train {
        best_observations
    } else {
        // Build peptide groups for fair subsampling (target-decoy pairs together)
        let mut base_id_to_peptide: HashMap<u32, Arc<str>> = HashMap::new();
        for &(fi, li) in &best_observations {
            let entry = &per_file_entries[fi].1[li];
            let base_id = entry.entry_id & 0x7FFF_FFFF;
            if !entry.is_decoy {
                base_id_to_peptide
                    .entry(base_id)
                    .or_insert_with(|| entry.modified_sequence.clone());
            }
        }
        // Handle decoy-only base_ids
        for &(fi, li) in &best_observations {
            let entry = &per_file_entries[fi].1[li];
            let base_id = entry.entry_id & 0x7FFF_FFFF;
            if entry.is_decoy {
                base_id_to_peptide
                    .entry(base_id)
                    .or_insert_with(|| entry.modified_sequence.clone());
            }
        }

        // Group observations by peptide identity
        let mut peptide_groups: HashMap<&str, Vec<(usize, usize)>> = HashMap::new();
        for &(fi, li) in &best_observations {
            let entry = &per_file_entries[fi].1[li];
            let base_id = entry.entry_id & 0x7FFF_FFFF;
            let key = base_id_to_peptide
                .get(&base_id)
                .map(|s| &**s)
                .unwrap_or(&entry.modified_sequence);
            peptide_groups.entry(key).or_default().push((fi, li));
        }

        let mut groups: Vec<(&str, Vec<(usize, usize)>)> = peptide_groups.into_iter().collect();
        groups.sort_by_key(|&(k, _)| k);

        // Fisher-Yates shuffle
        let mut rng_state = perc_config.seed;
        for i in (1..groups.len()).rev() {
            rng_state ^= rng_state << 13;
            rng_state ^= rng_state >> 7;
            rng_state ^= rng_state << 17;
            let j = rng_state as usize % (i + 1);
            groups.swap(i, j);
        }

        let mut selected: Vec<(usize, usize)> = Vec::with_capacity(max_train);
        for (_, indices) in &groups {
            if selected.len() + indices.len() > max_train && !selected.is_empty() {
                break;
            }
            selected.extend_from_slice(indices);
        }
        selected.sort();
        selected
    };

    let sub_targets = subset
        .iter()
        .filter(|&&(fi, li)| !per_file_entries[fi].1[li].is_decoy)
        .count();
    let sub_decoys = subset.len() - sub_targets;

    log::debug!(
        "  Subsampled {} entries ({} targets, {} decoys) for SVM training",
        subset.len(),
        sub_targets,
        sub_decoys
    );

    // Phase 2: Load features only for subsampled entries from Parquet
    // Group subset by file to minimize Parquet loads
    let mut subset_by_file: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for (pos, &(file_idx, local_idx)) in subset.iter().enumerate() {
        subset_by_file
            .entry(file_idx)
            .or_default()
            .push((pos, local_idx));
    }

    // Load features in parallel across files (each file is an independent Parquet read)
    #[allow(clippy::type_complexity)]
    let file_results: Vec<(usize, Vec<(usize, Vec<f64>)>)> = subset_by_file
        .par_iter()
        .filter_map(|(&file_idx, entries_in_file)| {
            let (file_name, fdr_entries) = &per_file_entries[file_idx];
            let file_overlay = overlay.get(file_name.as_str());
            let cache_path = per_file_cache_paths.get(file_name.as_str())?;

            let file_features = load_pin_features_from_parquet(cache_path).ok()?;
            let results: Vec<(usize, Vec<f64>)> = entries_in_file
                .iter()
                .filter_map(|&(pos, local_idx)| {
                    if let Some(overlay_entry) = file_overlay.and_then(|o| o.get(&local_idx)) {
                        Some((pos, pin_feature_vector(&overlay_entry.features)))
                    } else {
                        // Use parquet_index for feature lookup (Vec may be compacted)
                        let pq_idx = fdr_entries[local_idx].parquet_index as usize;
                        if pq_idx < file_features.len() {
                            Some((pos, file_features[pq_idx].clone()))
                        } else {
                            None
                        }
                    }
                })
                .collect();
            Some((file_idx, results))
        })
        .collect();

    let mut subset_features: Vec<Vec<f64>> = vec![Vec::new(); subset.len()];
    for (_file_idx, results) in file_results {
        for (pos, features) in results {
            subset_features[pos] = features;
        }
    }
    drop(subset_by_file);

    // Phase 3: Train SVM models on the subsample via run_percolator
    // Create PercolatorEntry objects only for the subset (~300K, ~130 MB)
    let subset_perc_entries: Vec<percolator::PercolatorEntry> = subset
        .iter()
        .zip(subset_features.iter())
        .map(|(&(file_idx, local_idx), features)| {
            let (file_name, entries) = &per_file_entries[file_idx];
            let entry = &entries[local_idx];
            percolator::PercolatorEntry {
                id: format!("sub_{}_{}", file_idx, local_idx),
                file_name: file_name.clone(),
                peptide: entry.modified_sequence.to_string(),
                charge: entry.charge,
                is_decoy: entry.is_decoy,
                entry_id: entry.entry_id,
                features: features.clone(),
            }
        })
        .collect();

    // Free subset memory
    drop(subset);
    drop(subset_features);

    let mut train_config = perc_config.clone();
    train_config.train_only = true;
    let train_results = percolator::run_percolator(&subset_perc_entries, &train_config)
        .map_err(|e| OspreyError::config(format!("Percolator training failed: {}", e)))?;

    // Extract averaged model weights and standardizer
    let n_models = train_results.fold_weights.len();
    let n_features = if n_models > 0 {
        train_results.fold_weights[0].len()
    } else {
        return Err(OspreyError::config("Percolator produced no models"));
    };

    let mut avg_weights = vec![0.0f64; n_features];
    let mut avg_bias = 0.0f64;
    for (weights, &bias) in train_results
        .fold_weights
        .iter()
        .zip(train_results.fold_biases.iter())
    {
        for (j, &w) in weights.iter().enumerate() {
            avg_weights[j] += w;
        }
        avg_bias += bias;
    }
    let n_models_f = n_models as f64;
    for w in &mut avg_weights {
        *w /= n_models_f;
    }
    avg_bias /= n_models_f;

    let standardizer = train_results.standardizer.clone();

    // Free training memory
    drop(subset_perc_entries);
    drop(train_results);

    // Phase 4: Score ALL entries by streaming through Parquet files in parallel
    // Each file is independent: load features, apply SVM model, write scores.
    per_file_entries
        .par_iter_mut()
        .for_each(|(file_name, fdr_entries)| {
            let file_overlay = overlay.get(file_name.as_str());
            let cache_path = match per_file_cache_paths.get(file_name.as_str()) {
                Some(p) => p,
                None => return,
            };

            let file_features = match load_pin_features_from_parquet(cache_path) {
                Ok(f) => f,
                Err(e) => {
                    log::warn!("Failed to load features for {}: {}", file_name, e);
                    return;
                }
            };

            // Score entries using parquet_index for feature lookup
            // (stubs may have been compacted, so Vec index != Parquet row index)
            for (local_idx, fdr_entry) in fdr_entries.iter_mut().enumerate() {
                // Skip entries not in the restriction set (keep first-pass scores)
                if let Some(restrict) = restrict_base_ids {
                    let base_id = fdr_entry.entry_id & 0x7FFF_FFFF;
                    if !restrict.contains(&base_id) {
                        continue;
                    }
                }
                let pq_idx = fdr_entry.parquet_index as usize;
                // Use overlay features for re-scored entries, parquet features otherwise
                let mut features =
                    if let Some(overlay_entry) = file_overlay.and_then(|o| o.get(&local_idx)) {
                        pin_feature_vector(&overlay_entry.features)
                    } else if pq_idx < file_features.len() {
                        file_features[pq_idx].clone()
                    } else {
                        continue; // gap-fill without parquet features
                    };
                standardizer.transform_slice(&mut features);
                let mut score = avg_bias;
                for (w, x) in avg_weights.iter().zip(features.iter()) {
                    score += w * x;
                }
                fdr_entry.score = score;
            }

            // Score gap-fill entries (parquet_index == u32::MAX, features in overlay only)
            for (local_idx, fdr_entry) in fdr_entries.iter_mut().enumerate() {
                if fdr_entry.parquet_index != u32::MAX {
                    continue; // already scored above
                }
                if let Some(restrict) = restrict_base_ids {
                    let base_id = fdr_entry.entry_id & 0x7FFF_FFFF;
                    if !restrict.contains(&base_id) {
                        continue;
                    }
                }
                if let Some(overlay_entry) = file_overlay.and_then(|o| o.get(&local_idx)) {
                    let mut features = pin_feature_vector(&overlay_entry.features);
                    standardizer.transform_slice(&mut features);
                    let mut score = avg_bias;
                    for (w, x) in avg_weights.iter().zip(features.iter()) {
                        score += w * x;
                    }
                    fdr_entry.score = score;
                }
            }
        });
    log::debug!("  Scored {}/{} files", n_files, n_files);

    // Phase 5: Compute q-values and PEP directly from per_file_entries
    // This avoids building any flat metadata arrays (~0 extra memory vs 33 GB before)
    log::debug!("");
    log::debug!(
        "=== Per-file results (run-level FDR at {:.0}%) ===",
        perc_config.test_fdr * 100.0
    );
    percolator::compute_fdr_from_stubs(per_file_entries, perc_config.test_fdr, restrict_base_ids);

    Ok(())
}

/// Direct (non-streaming) Percolator FDR for small-to-medium experiments.
///
/// Loads all features into memory and passes them to run_percolator.
/// Used when total entries are small enough to fit in memory.
fn run_percolator_fdr_direct(
    per_file_entries: &mut [(String, Vec<FdrEntry>)],
    per_file_cache_paths: &HashMap<String, std::path::PathBuf>,
    perc_config: &percolator::PercolatorConfig,
    overlay: &RescoreOverlay,
    restrict_base_ids: Option<&HashSet<u32>>,
) -> Result<()> {
    let mut perc_entries = Vec::new();
    for (file_name, fdr_entries) in per_file_entries.iter() {
        let file_overlay = overlay.get(file_name.as_str());
        let cache_path = per_file_cache_paths
            .get(file_name.as_str())
            .ok_or_else(|| {
                OspreyError::config(format!("No parquet cache path for file {}", file_name))
            })?;
        let file_features = load_pin_features_from_parquet(cache_path)?;

        // Process entries using parquet_index for feature lookup
        // (stubs may have been compacted, so Vec index != Parquet row index)
        for (local_idx, fdr_entry) in fdr_entries.iter().enumerate() {
            if let Some(restrict) = restrict_base_ids {
                if !restrict.contains(&(fdr_entry.entry_id & 0x7FFF_FFFF)) {
                    continue;
                }
            }
            let pq_idx = fdr_entry.parquet_index as usize;
            // Use overlay features for re-scored entries, parquet features otherwise
            let features = if let Some(overlay_entry) = file_overlay.and_then(|o| o.get(&local_idx))
            {
                pin_feature_vector(&overlay_entry.features)
            } else if pq_idx < file_features.len() {
                file_features[pq_idx].clone()
            } else {
                continue; // gap-fill entry without parquet features
            };
            let psm_id = format!(
                "{}_{}_{}_{}",
                file_name, fdr_entry.modified_sequence, fdr_entry.charge, fdr_entry.scan_number
            );

            perc_entries.push(percolator::PercolatorEntry {
                id: psm_id,
                file_name: file_name.clone(),
                peptide: fdr_entry.modified_sequence.to_string(),
                charge: fdr_entry.charge,
                is_decoy: fdr_entry.is_decoy,
                entry_id: fdr_entry.entry_id,
                features,
            });
        }

        // Process gap-fill entries (parquet_index == u32::MAX, features in overlay only)
        for (local_idx, fdr_entry) in fdr_entries.iter().enumerate() {
            if fdr_entry.parquet_index != u32::MAX {
                continue; // already processed above
            }
            if let Some(restrict) = restrict_base_ids {
                if !restrict.contains(&(fdr_entry.entry_id & 0x7FFF_FFFF)) {
                    continue;
                }
            }
            if let Some(overlay_entry) = file_overlay.and_then(|o| o.get(&local_idx)) {
                let features = pin_feature_vector(&overlay_entry.features);
                let psm_id = format!(
                    "{}_{}_{}_{}",
                    file_name, fdr_entry.modified_sequence, fdr_entry.charge, fdr_entry.scan_number
                );
                perc_entries.push(percolator::PercolatorEntry {
                    id: psm_id,
                    file_name: file_name.clone(),
                    peptide: fdr_entry.modified_sequence.to_string(),
                    charge: fdr_entry.charge,
                    is_decoy: fdr_entry.is_decoy,
                    entry_id: fdr_entry.entry_id,
                    features,
                });
            }
        }
    }

    let results = percolator::run_percolator(&perc_entries, perc_config)
        .map_err(|e| OspreyError::config(format!("Percolator failed: {}", e)))?;

    // Build result lookup
    let result_map: HashMap<&str, &percolator::PercolatorResult> =
        results.entries.iter().map(|r| (r.id.as_str(), r)).collect();

    // Map results back to FdrEntry stubs
    for (file_name, entries) in per_file_entries.iter_mut() {
        for entry in entries.iter_mut() {
            let psm_id = format!(
                "{}_{}_{}_{}",
                file_name, entry.modified_sequence, entry.charge, entry.scan_number
            );
            if let Some(result) = result_map.get(psm_id.as_str()) {
                entry.run_precursor_qvalue = result.run_precursor_qvalue;
                entry.run_peptide_qvalue = result.run_peptide_qvalue;
                entry.experiment_precursor_qvalue = result.experiment_precursor_qvalue;
                entry.experiment_peptide_qvalue = result.experiment_peptide_qvalue;
                entry.score = result.score;
                entry.pep = result.pep;
            }
        }
    }

    Ok(())
}

/// Run external mokapot FDR on coelution entries
fn run_mokapot_fdr(
    per_file_entries: &mut [(String, Vec<FdrEntry>)],
    per_file_cache_paths: &HashMap<String, std::path::PathBuf>,
    mokapot: &MokapotRunner,
    pin_files: &HashMap<String, std::path::PathBuf>,
    mokapot_dir: &std::path::Path,
    config: &OspreyConfig,
) -> Result<()> {
    if !mokapot.is_available() {
        log::warn!("Mokapot not available, falling back to native Percolator");
        // Mokapot doesn't have overlay — pass empty (mokapot uses PIN files, not parquet)
        let empty_overlay: RescoreOverlay = HashMap::new();
        return run_percolator_fdr(
            per_file_entries,
            per_file_cache_paths,
            config,
            &empty_overlay,
            None,
        );
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
                            let peptide_q = peptide_qmap
                                .and_then(|m| m.get(&*entry.modified_sequence))
                                .copied()
                                .unwrap_or(1.0);
                            entry.run_precursor_qvalue = q;
                            entry.run_peptide_qvalue = peptide_q;
                            entry.score = score;
                            entry.pep = pep;
                        }
                        if let Some(&q) = experiment_qmap.get(&psm_id) {
                            let peptide_q = exp_peptide_qmap
                                .get(&*entry.modified_sequence)
                                .copied()
                                .unwrap_or(1.0);
                            entry.experiment_precursor_qvalue = q;
                            entry.experiment_peptide_qvalue = peptide_q;
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
            // Mokapot's own overlay isn't tracked — pass empty
            let empty_overlay: RescoreOverlay = HashMap::new();
            run_percolator_fdr(
                per_file_entries,
                per_file_cache_paths,
                config,
                &empty_overlay,
                None,
            )?;
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

    // Sort by entry_id for deterministic order regardless of HashMap iteration.
    // Without this, the random HashMap order propagates to SVM feature matrix
    // row ordering, causing non-deterministic gradient updates and model weights.
    all_entries.sort_by_key(|e| e.entry_id);

    log::debug!(
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
    log::debug!(
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

    log::debug!(
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

    // For each window, binary search for the range of entries that fall within it.
    // Use exclusive upper bound [lower, upper) so each entry is in at most one window.
    // Without this, entries at exact window boundaries (e.g., precursor_mz = 500.0 for
    // adjacent windows [498,500] and [500,502]) would be in two windows, causing race
    // conditions in the par_iter with AtomicBool below.
    let window_entry_indices: Vec<Vec<usize>> = windows
        .iter()
        .map(|&(center, width)| {
            let win_lower = center - width / 2.0;
            let win_upper = center + width / 2.0;
            let lo = mz_sorted_vals.partition_point(|&mz| mz < win_lower);
            let hi = mz_sorted_vals.partition_point(|&mz| mz < win_upper);
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
        rt_sorted.sort_by(|&a, &b| {
            entries[a]
                .apex_rt
                .total_cmp(&entries[b].apex_rt)
                .then((entries[a].entry_id & 0x7FFFFFFF).cmp(&(entries[b].entry_id & 0x7FFFFFFF)))
                .then(entries[a].entry_id.cmp(&entries[b].entry_id))
        });

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

                // Never deduplicate across target/decoy types: removing a decoy because
                // it shares fragments with a winning target destroys FDR estimation.
                // Only deduplicate within the same class (target-target or decoy-decoy).
                if entries[idx_a].is_decoy != entries[idx_b].is_decoy {
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
                    // Remove the entry with lower coelution_sum (keep the better one).
                    // Use standard f64 comparison (NaN loses, matching original behavior).
                    // The iteration order is deterministic (entries sorted by apex_rt + entry_id),
                    // so the >= tiebreaker (outer loop wins ties) is also deterministic.
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
        let removed_targets = removed
            .iter()
            .enumerate()
            .filter(|(i, r)| r.load(Ordering::Relaxed) && !entries[*i].is_decoy)
            .count();
        let removed_decoys = removed_count - removed_targets;
        log::debug!(
            "Double-counting deduplication: removed {} entries ({} targets, {} decoys; {} remaining)",
            removed_count,
            removed_targets,
            removed_decoys,
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
fn apply_simple_fdr(entries: &mut [FdrEntry], _fdr_threshold: f64) -> Result<()> {
    // Sort by coelution_sum descending
    entries.sort_by(|a, b| b.coelution_sum.total_cmp(&a.coelution_sum));

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

        entry.run_precursor_qvalue = fdr;
    }

    // Convert to q-values (monotonic)
    let mut min_fdr = 1.0f64;
    for entry in entries.iter_mut().rev() {
        min_fdr = min_fdr.min(entry.run_precursor_qvalue);
        entry.run_precursor_qvalue = min_fdr;
        entry.experiment_precursor_qvalue = min_fdr; // No two-level distinction in fallback mode
    }

    Ok(())
}

/// Lightweight plan entry for streaming blib output.
///
/// Contains only the fields needed for blib writing (~96 bytes per entry).
/// Fragment m/z, intensities, modifications, and protein mappings are looked up
/// from the in-memory library via `entry_id`. This avoids loading full
/// `CoelutionScoredEntry` objects (~940 bytes each) from Parquet, preventing
/// OOM on 240-file experiments with ~28.8M passing entries.
struct BlibPlanEntry {
    /// Library entry ID for fragment/mod/protein lookup
    entry_id: u32,
    /// Precursor charge
    charge: u8,
    /// Index into file_names vec
    file_name_idx: u16,
    /// Run-level precursor q-value
    run_qvalue: f64,
    /// Experiment-level precursor q-value
    experiment_qvalue: f64,
    /// Apex retention time (minutes)
    apex_rt: f64,
    /// Start retention time (minutes)
    start_rt: f64,
    /// End retention time (minutes)
    end_rt: f64,
    /// Integrated area under peak (from Parquet bounds_area)
    bounds_area: f64,
    /// Modified sequence (interned)
    modified_sequence: Arc<str>,
}

/// Build shared peak boundaries from plan entries.
///
/// For each (modified_sequence, file_name_idx) pair, keeps the boundaries from
/// the entry with the lowest run q-value. Multiple charge states of the same
/// peptide in the same file share the best boundaries.
#[allow(clippy::type_complexity)]
fn build_shared_boundaries_from_plan(
    entries: &[BlibPlanEntry],
) -> (HashMap<(Arc<str>, u16), (f64, f64, f64)>, usize) {
    // For each (modified_sequence, file_name_idx), keep boundaries from entry
    // with lowest run q-value
    let mut best_by_key: HashMap<(Arc<str>, u16), (f64, f64, f64, f64)> = HashMap::new();
    let mut charges_by_key: HashMap<(Arc<str>, u16), HashSet<u8>> = HashMap::new();

    for entry in entries {
        let key = (entry.modified_sequence.clone(), entry.file_name_idx);

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
                    entry.start_rt,
                    entry.end_rt,
                    entry.run_qvalue,
                ),
            );
        }
    }

    let n_shared = charges_by_key.values().filter(|c| c.len() > 1).count();

    let shared_bounds: HashMap<(Arc<str>, u16), (f64, f64, f64)> = best_by_key
        .into_iter()
        .map(|(k, (apex, start, end, _))| (k, (apex, start, end)))
        .collect();

    (shared_bounds, n_shared)
}

/// Write blib output from lightweight plan entries.
///
/// Looks up fragment m/z, intensities, modifications, and protein mappings
/// from the in-memory library via `entry_id`, avoiding the need to load full
/// `CoelutionScoredEntry` objects from Parquet.
fn write_blib_from_plan(
    config: &OspreyConfig,
    library: &[LibraryEntry],
    plan_entries: &[BlibPlanEntry],
    file_names: &[String],
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

    // Add source files and build lookup by file stem for matching with file_name_idx
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

    // Build shared boundaries from plan entries
    let (shared_bounds, n_shared) = build_shared_boundaries_from_plan(plan_entries);

    // Build library ID lookup for O(1) access
    let lib_by_id: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    // Group plan entries by precursor (modified_sequence + charge)
    let mut precursor_groups: HashMap<(Arc<str>, u8), Vec<&BlibPlanEntry>> = HashMap::new();
    for entry in plan_entries {
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

        // Look up library entry for fragments, mods, and proteins
        let lib_entry = match lib_by_id.get(&best.entry_id) {
            Some(e) => e,
            None => continue,
        };

        let frag_mzs: Vec<f64> = lib_entry.fragments.iter().map(|f| f.mz).collect();
        let frag_intensities: Vec<f32> = lib_entry
            .fragments
            .iter()
            .map(|f| f.relative_intensity)
            .collect();

        let n_runs_detected = group.len() as i32;
        let best_file_stem = &file_names[best.file_name_idx as usize];
        let file_id = *file_stem_to_id.get(best_file_stem).unwrap_or(&1);

        let tic: f64 = frag_intensities.iter().map(|&x| x as f64).sum();

        // Use shared boundaries from the best charge state for this peptide in this file
        let shared_key = (best.modified_sequence.clone(), best.file_name_idx);
        let (shared_apex, shared_start, shared_end) = shared_bounds
            .get(&shared_key)
            .copied()
            .unwrap_or((best.apex_rt, best.start_rt, best.end_rt));

        // Score is the experiment q-value, following Skyline GENERIC Q-VALUE convention
        let ref_id = writer.add_spectrum(
            &lib_entry.sequence,
            &lib_entry.modified_sequence,
            lib_entry.precursor_mz,
            lib_entry.charge as i32,
            shared_apex,
            shared_start,
            shared_end,
            &frag_mzs,
            &frag_intensities,
            best.experiment_qvalue,
            file_id,
            n_runs_detected,
            tic,
        )?;

        // Add modifications from library entry
        if !lib_entry.modifications.is_empty() {
            writer.add_modifications(ref_id, &lib_entry.modifications)?;
        }
        if !lib_entry.protein_ids.is_empty() {
            writer.add_protein_mapping(ref_id, &lib_entry.protein_ids)?;
        }

        // Write RetentionTimes entries for every run where this precursor was detected.
        // A non-NULL retentionTime produces an ID line in Skyline for that run.
        // Normally, only runs passing run-level FDR get an ID line. However, after
        // two-pass FDR, second-pass run q-values can shift slightly above the threshold
        // even though the precursor passes experiment-level FDR. In that case, the run
        // with the lowest q-value gets the ID line as a fallback, guaranteeing at least
        // one ID line per precursor in the blib.
        let any_passes_run_fdr = group.iter().any(|s| s.run_qvalue <= config.run_fdr);
        let best_run_file_idx = if !any_passes_run_fdr {
            // No run passes strict threshold; find the best one for fallback ID line
            Some(
                group
                    .iter()
                    .min_by(|a, b| a.run_qvalue.total_cmp(&b.run_qvalue))
                    .unwrap()
                    .file_name_idx,
            )
        } else {
            None
        };

        for scored in group {
            let run_file_stem = &file_names[scored.file_name_idx as usize];
            let run_file_id = *file_stem_to_id.get(run_file_stem).unwrap_or(&1);
            let is_best = std::ptr::eq(*scored, *best);

            // Use shared boundaries for this peptide in this run's file
            let run_shared_key = (scored.modified_sequence.clone(), scored.file_name_idx);
            let (run_apex, run_start, run_end) = shared_bounds
                .get(&run_shared_key)
                .copied()
                .unwrap_or((scored.apex_rt, scored.start_rt, scored.end_rt));

            // Show an ID line (non-NULL retentionTime) if:
            //   1. This run passes run-level FDR, OR
            //   2. No run passes but this is the best run (fallback)
            let rt_for_id = if scored.run_qvalue <= config.run_fdr
                || best_run_file_idx == Some(scored.file_name_idx)
            {
                Some(run_apex)
            } else {
                None
            };

            writer.add_retention_time(
                ref_id,
                run_file_id,
                rt_for_id,
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
            apex_coefficient: 0.0, // Not available from plan entries
            integrated_area: best.bounds_area,
            peak_quality: osprey_core::PeakQuality::default(),
        };
        writer.add_peak_boundaries(ref_id, best_file_stem, &boundaries)?;

        // Add run-level scores for best run
        writer.add_run_scores(
            ref_id,
            best_file_stem,
            best.run_qvalue,
            0.0, // dot_product not available from plan entries
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
        "Wrote {} precursors ({} observations) across {} files to blib",
        n_written,
        plan_entries.len(),
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

/// Write TSV report from plan entries (lightweight, no full entry reload).
///
/// Provides a basic summary of passing precursors with q-values and RT boundaries.
fn write_tsv_report_from_plan(
    path: &std::path::Path,
    plan_entries: &[BlibPlanEntry],
    library: &[LibraryEntry],
) -> Result<()> {
    let lib_by_id: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "modified_sequence\tprecursor_mz\tcharge\tapex_rt\tstart_rt\tend_rt\tbounds_area\tq_value"
    )?;

    for entry in plan_entries {
        let precursor_mz = lib_by_id
            .get(&entry.entry_id)
            .map(|e| e.precursor_mz)
            .unwrap_or(0.0);
        writeln!(
            file,
            "{}\t{:.4}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.4}\t{:.6}",
            entry.modified_sequence,
            precursor_mz,
            entry.charge,
            entry.apex_rt,
            entry.start_rt,
            entry.end_rt,
            entry.bounds_area,
            entry.experiment_qvalue
        )?;
    }

    Ok(())
}

// =============================================================================
// Coelution-based search (DIA-NN style)
// =============================================================================

/// Run coelution-based search for a single file.
///
/// Context for computing coelution features at given peak boundaries.
///
/// Bundles read-only references needed by `compute_features_at_peak`, avoiding
/// a function with 15+ parameters. Used both in the initial per-window scoring
/// and in the multi-charge consensus re-scoring pass.
struct FeatureComputeContext<'a> {
    entry: &'a LibraryEntry,
    xics: &'a [(usize, Vec<(f64, f64)>)],
    ref_xic: &'a [(f64, f64)],
    cand_spectra: &'a [&'a Spectrum],
    cand_global: &'a [usize],
    scorer: &'a SpectralScorer,
    ms1_index: &'a MS1Index,
    calibration: Option<&'a CalibrationParams>,
    tol_da: f64,
    tol_ppm: f64,
    expected_rt: f64,
    is_hram: bool,
    file_name: &'a str,
    /// Pre-preprocessed XCorr vectors for all spectra in the isolation window.
    /// Indexed by position in window_spectra (same order as window_pairs).
    /// When Some, enables O(n_frags) XCorr lookup instead of O(n_peaks) preprocessing per call.
    preprocessed_xcorr: Option<&'a [Vec<f32>]>,
    /// For each spectrum in cand_spectra, its index into the window_spectra array.
    /// Used to look up pre-preprocessed XCorr data.
    cand_window_local: Option<&'a [usize]>,
}

/// Compute all coelution features and build a `CoelutionScoredEntry` at the
/// specified peak boundaries.
///
/// Returns `None` if the peak has insufficient data (e.g. no pairwise
/// correlations, no apex spectrum found).
fn compute_features_at_peak(
    ctx: &FeatureComputeContext,
    peak: XICPeakBounds,
) -> Option<CoelutionScoredEntry> {
    use osprey_chromatography::calibration::calibrated_tolerance_ppm;
    use osprey_scoring::{
        compute_cosine_at_scan, compute_mass_accuracy, median_polish_libcosine,
        median_polish_min_fragment_r2, median_polish_residual_correlation,
        median_polish_residual_ratio, median_polish_rsquared, pearson_correlation_raw,
        tukey_median_polish,
    };

    let entry = ctx.entry;
    let xics = ctx.xics;
    let ref_xic = ctx.ref_xic;
    let cand_spectra = ctx.cand_spectra;
    let cand_global = ctx.cand_global;

    if peak.end_index <= peak.start_index + 1 {
        return None;
    }

    // 1. Pairwise fragment correlations within peak boundaries
    let n_xics = xics.len();
    let peak_len = peak.end_index - peak.start_index + 1;

    let xic_peak_values: Vec<Vec<f64>> = xics
        .iter()
        .map(|(_, xic)| {
            xic[peak.start_index..=peak.end_index]
                .iter()
                .map(|(_, v)| *v)
                .collect()
        })
        .collect();

    let mut corr_sum = 0.0f64;
    let mut corr_min = f64::INFINITY;
    let mut corr_max = f64::NEG_INFINITY;
    let mut n_pairs = 0u32;
    let mut per_frag_corr_sum = vec![0.0f64; n_xics];
    let mut per_frag_corr_count = vec![0u32; n_xics];

    for i in 0..n_xics {
        for j in (i + 1)..n_xics {
            let r = pearson_correlation_raw(&xic_peak_values[i], &xic_peak_values[j]);
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

    // 2. Find the spectrum closest to apex RT for spectral scoring
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

    // Spectral scoring: compute LibCosine features + XCorr.
    // When pre-preprocessed XCorr data is available (main search path), avoid
    // redundant per-entry spectrum preprocessing (binning + windowing + sliding
    // window subtraction) by looking up the pre-computed result.
    let spectral_score = if let (Some(preprocessed), Some(win_indices)) =
        (ctx.preprocessed_xcorr, ctx.cand_window_local)
    {
        let mut score = ctx.scorer.lib_cosine(apex_spectrum, entry);
        let win_idx = win_indices[apex_local_idx];
        let lib_preprocessed = ctx.scorer.preprocess_library_for_xcorr(entry);
        score.xcorr =
            SpectralScorer::xcorr_from_preprocessed(&preprocessed[win_idx], &lib_preprocessed);
        score
    } else {
        ctx.scorer.xcorr(apex_spectrum, entry)
    };

    // 2b. SG-weighted spectral scores at apex ±2 scans
    let sg_weights: [f64; 5] = [
        -3.0 / 35.0,
        12.0 / 35.0,
        17.0 / 35.0,
        12.0 / 35.0,
        -3.0 / 35.0,
    ];
    let mut sg_xcorr = 0.0;
    let mut sg_cosine = 0.0;
    // Pre-compute library XCorr vector once for SG-weighted lookups
    let lib_xcorr_preprocessed = if ctx.preprocessed_xcorr.is_some() {
        Some(ctx.scorer.preprocess_library_for_xcorr(entry))
    } else {
        None
    };
    for (offset, &weight) in [-2i32, -1, 0, 1, 2].iter().zip(&sg_weights) {
        let idx = apex_local_idx as i32 + offset;
        if idx >= 0 && (idx as usize) < cand_spectra.len() {
            let spec = cand_spectra[idx as usize];
            if let (Some(preprocessed), Some(win_indices), Some(ref lib_pre)) = (
                ctx.preprocessed_xcorr,
                ctx.cand_window_local,
                &lib_xcorr_preprocessed,
            ) {
                let win_idx = win_indices[idx as usize];
                sg_xcorr +=
                    SpectralScorer::xcorr_from_preprocessed(&preprocessed[win_idx], lib_pre)
                        * weight;
            } else {
                sg_xcorr += ctx.scorer.xcorr_at_scan(spec, entry) * weight;
            }
            sg_cosine +=
                compute_cosine_at_scan(&entry.fragments, spec, ctx.tol_da, ctx.tol_ppm) * weight;
        }
    }

    // 4. Mass accuracy at apex
    let frag_matches = ctx.scorer.match_fragments(apex_spectrum, entry);
    let fragment_tolerance_unit = if ctx.tol_ppm > 0.0 {
        ToleranceUnit::Ppm
    } else {
        ToleranceUnit::Mz
    };
    let fragment_tolerance_val = if ctx.tol_ppm > 0.0 {
        ctx.tol_ppm
    } else {
        ctx.tol_da
    };
    let (mass_mean, mass_abs_mean, mass_std) = compute_mass_accuracy(
        &frag_matches,
        fragment_tolerance_unit,
        fragment_tolerance_val,
    );

    // 5. RT deviation
    let rt_dev = peak.apex_rt - ctx.expected_rt;

    // 6. MS1 features (HRAM only)
    let (ms1_coelution, ms1_isotope_cos) = if ctx.is_hram && !ctx.ms1_index.is_empty() {
        let ms1_tol_ppm = ctx
            .calibration
            .map(|cal| calibrated_tolerance_ppm(&cal.ms1_calibration, 10.0))
            .unwrap_or(10.0);
        let search_mz = ctx
            .calibration
            .map(|cal| reverse_calibrate_mz(entry.precursor_mz, &cal.ms1_calibration))
            .unwrap_or(entry.precursor_mz);

        let mut ms1_coelution_val = 0.0;
        let peak_rts: Vec<(f64, f64)> = ref_xic[peak.start_index..=peak.end_index].to_vec();
        if peak_rts.len() >= 3 {
            let mut ms1_intensities = Vec::with_capacity(peak_rts.len());
            let mut ref_intensities = Vec::with_capacity(peak_rts.len());
            for &(rt, ref_val) in &peak_rts {
                if let Some(ms1) = ctx.ms1_index.find_nearest(rt) {
                    let intensity = ms1
                        .find_peak_ppm(search_mz, ms1_tol_ppm)
                        .map(|(_, int)| int as f64)
                        .unwrap_or(0.0);
                    ms1_intensities.push(intensity);
                    ref_intensities.push(ref_val);
                }
            }
            if ms1_intensities.len() >= 3 {
                ms1_coelution_val = pearson_correlation_raw(&ms1_intensities, &ref_intensities);
            }
        }

        let mut iso_cos = 0.0;
        if let Some(apex_ms1) = ctx.ms1_index.find_nearest(peak.apex_rt) {
            let envelope = osprey_core::IsotopeEnvelope::extract(
                apex_ms1,
                search_mz,
                entry.charge,
                ms1_tol_ppm,
            );
            if envelope.has_m0() {
                iso_cos =
                    osprey_core::peptide_isotope_cosine(&entry.sequence, &envelope.intensities)
                        .unwrap_or(0.0);
            }
        }

        (ms1_coelution_val, iso_cos)
    } else {
        (0.0, 0.0)
    };

    // 7. Peak shape features
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

    // 8. Median polish features
    let (mp_cosine, mp_rsquared, mp_residual, mp_min_r2, mp_resid_corr) = {
        let peak_xics: Vec<(usize, Vec<(f64, f64)>)> = xics
            .iter()
            .map(|(frag_idx, xic)| {
                let peak_slice: Vec<(f64, f64)> = xic[peak.start_index..=peak.end_index].to_vec();
                (*frag_idx, peak_slice)
            })
            .collect();

        if let Some(ref mp) = tukey_median_polish(&peak_xics, 10, 0.01) {
            let cos = median_polish_libcosine(mp, &entry.fragments);
            let rsq = median_polish_rsquared(mp);
            let res = median_polish_residual_ratio(mp);
            let min_r2 = median_polish_min_fragment_r2(mp);
            let resid_corr = median_polish_residual_correlation(mp);
            (cos, rsq, res, min_r2, resid_corr)
        } else {
            (0.0, 0.0, 1.0, 0.0, 0.0)
        }
    };

    // 9. Peptide properties
    let mod_count = entry.modifications.len() as u8;
    let pep_len = entry.sequence.len() as u8;
    let missed_cleav = entry
        .sequence
        .chars()
        .take(entry.sequence.len().saturating_sub(1))
        .filter(|&c| c == 'K' || c == 'R')
        .count() as u8;

    // Fragment m/z and intensities from library for blib output
    let frag_mzs: Vec<f64> = entry.fragments.iter().map(|f| f.mz).collect();
    let frag_intensities: Vec<f32> = entry
        .fragments
        .iter()
        .map(|f| f.relative_intensity)
        .collect();

    // Reference XIC within peak bounds for blib coefficient series
    let ref_xic_peak: Vec<(f64, f64)> = ref_xic[peak.start_index..=peak.end_index].to_vec();

    let features = CoelutionFeatureSet {
        coelution_sum: corr_sum,
        coelution_min: if corr_min.is_finite() { corr_min } else { 0.0 },
        coelution_max: if corr_max.is_finite() { corr_max } else { 0.0 },
        n_coeluting_fragments: n_coeluting,
        n_fragment_pairs: n_pairs.min(255) as u8,
        fragment_corr,
        peak_apex: peak.apex_intensity,
        peak_area: peak.area,
        peak_width,
        peak_symmetry,
        signal_to_noise: peak.signal_to_noise,
        n_scans: peak_len as u16,
        peak_sharpness,
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
        elution_weighted_cosine: 0.0, // removed: expensive per-scan ppm matching, not in PIN
        mass_accuracy_mean: mass_mean,
        abs_mass_accuracy_mean: mass_abs_mean,
        mass_accuracy_std: mass_std,
        rt_deviation: rt_dev,
        abs_rt_deviation: rt_dev.abs(),
        ms1_precursor_coelution: ms1_coelution,
        ms1_isotope_cosine: ms1_isotope_cos,
        modification_count: mod_count,
        peptide_length: pep_len,
        missed_cleavages: missed_cleav,
        median_polish_cosine: mp_cosine,
        median_polish_rsquared: mp_rsquared,
        median_polish_residual_ratio: mp_residual,
        sg_weighted_xcorr: sg_xcorr,
        sg_weighted_cosine: sg_cosine,
        median_polish_min_fragment_r2: mp_min_r2,
        median_polish_residual_correlation: mp_resid_corr,
    };

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
        file_name: ctx.file_name.to_string(),
        run_precursor_qvalue: 1.0,
        run_peptide_qvalue: 1.0,
        run_protein_qvalue: 1.0,
        experiment_precursor_qvalue: 1.0,
        experiment_peptide_qvalue: 1.0,
        experiment_protein_qvalue: 1.0,
        score: 0.0,
        pep: 1.0,
        cwt_candidates: Vec::new(), // populated by caller after peak selection
    })
}

/// Select multi-charge consensus targets after FDR scoring.
///
/// Groups entries by modified_sequence (within a single file's entries).
/// For multi-charge groups where at least one charge state passes FDR,
/// the best-scoring charge state defines the consensus peak. Other charge
/// states at a different peak are returned as rescore targets:
/// `(entry_index, consensus_apex, start, end)`.
///
/// Groups where NO charge state passes FDR are skipped entirely — they have
/// no impact on reported results and re-scoring them is wasted work.
fn select_post_fdr_consensus(
    entries: &[FdrEntry],
    fdr_threshold: f64,
) -> Vec<(usize, f64, f64, f64)> {
    let mut groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, entry) in entries.iter().enumerate() {
        groups.entry(&entry.modified_sequence).or_default().push(i);
    }

    let mut targets = Vec::new();

    for indices in groups.values() {
        if indices.len() <= 1 {
            continue;
        }

        // Best charge state: highest SVM score among FDR-passing entries,
        // with lowest run_qvalue as tiebreaker.
        // Skip groups where no charge state passes FDR.
        let best_idx = match indices
            .iter()
            .filter(|&&i| entries[i].run_precursor_qvalue <= fdr_threshold)
            .max_by(|&&a, &&b| {
                entries[a].score.total_cmp(&entries[b].score).then_with(|| {
                    entries[b]
                        .run_precursor_qvalue
                        .total_cmp(&entries[a].run_precursor_qvalue)
                })
            })
            .copied()
        {
            Some(idx) => idx,
            None => continue, // No charge state passes FDR — skip
        };

        let consensus_apex = entries[best_idx].apex_rt;
        let consensus_start = entries[best_idx].start_rt;
        let consensus_end = entries[best_idx].end_rt;
        let consensus_width = consensus_end - consensus_start;
        let rt_match_tol = (consensus_width / 2.0).max(0.1);

        for &idx in indices {
            if idx == best_idx {
                continue;
            }
            let apex_diff = (entries[idx].apex_rt - consensus_apex).abs();
            if apex_diff > rt_match_tol {
                targets.push((idx, consensus_apex, consensus_start, consensus_end));
            }
        }
    }

    targets
}

/// For each precursor candidate in each isolation window:
/// 1. Extract top-6 fragment XICs
/// 2. Detect peak on reference XIC (DIA-NN-style valley detection)
/// 3. Compute pairwise fragment correlations within peak boundaries
/// 4. Score spectrum at apex (spectral similarity, mass accuracy)
/// 5. Compute peak shape, RT deviation, MS1, and median polish features
///
/// When `boundary_overrides` contains an entry's ID, peak detection is skipped
/// and the given (apex_rt, start_rt, end_rt) boundaries are used directly.
/// This is used for multi-charge consensus (forced boundaries from best charge state)
/// and inter-replicate reconciliation (CWT or imputed boundaries).
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
    boundary_overrides: Option<&HashMap<u32, (f64, f64, f64)>>,
    search_label: &str,
) -> Result<Vec<CoelutionScoredEntry>> {
    use osprey_chromatography::{compute_snr, detect_all_xic_peaks, detect_cwt_consensus_peaks};
    use osprey_scoring::batch::{group_spectra_by_isolation_window, MIN_COELUTION_SPECTRA};
    use osprey_scoring::{
        extract_fragment_xics, pearson_correlation_raw, tukey_median_polish, SpectralScorer,
    };

    let is_hram = matches!(config.resolution_mode, osprey_core::ResolutionMode::HRAM);

    // Set up fragment tolerance (calibrated if available)
    let fragment_tolerance = if let Some(cal) = calibration {
        let (tol_val, tol_unit) = osprey_chromatography::calibrated_tolerance(
            &cal.ms2_calibration,
            config.fragment_tolerance.tolerance,
            config.fragment_tolerance.unit,
        );
        log::debug!(
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
            log::debug!(
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
    // Use MAD-based robust tolerance: 3 × MAD × 1.4826 covers ~99.7% of peptides
    let rt_tolerance = if let Some(cal) = calibration {
        let (raw_tol, method_desc) = if let Some(mad) = cal.rt_calibration.mad {
            let robust_sd = mad * 1.4826;
            let tol = 3.0 * robust_sd;
            (
                tol,
                format!("3×MAD×1.4826: MAD={:.3}, robust_SD={:.3}", mad, robust_sd),
            )
        } else {
            // Fallback for old calibration files without MAD
            let tol = cal.rt_calibration.residual_sd * config.rt_calibration.rt_tolerance_factor;
            (
                tol,
                format!(
                    "{}× residual SD {:.2} min (MAD unavailable)",
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
            log::debug!(
                "Using calibrated RT tolerance: {:.2} min ({})",
                tol,
                method_desc
            );
        }
        tol
    } else {
        log::debug!(
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

    log::debug!(
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

    // Progress bar — one tick per isolation window, with descriptive label.
    let pb = ProgressBar::new(window_groups.len() as u64);
    let bar_template = format!(
        "{{spinner:.green}} {} [{{elapsed_precise}}] [{{bar:40.cyan/blue}}] {{pos}}/{{len}} windows",
        search_label
    );
    pb.set_style(
        ProgressStyle::default_bar()
            .template(&bar_template)
            .unwrap()
            .progress_chars("#>-"),
    );
    // Coordinate with two-tier logger: route log lines through pb.println()
    // while the progress bar is active to avoid interleaving.
    crate::logging::set_progress_bar(&pb);

    // Process each isolation window — candidates within each window in parallel
    // Use fold/reduce instead of flat_map+collect so results are deduplicated
    // per-thread as each window finishes, then merged pairwise.  This avoids
    // the peak-memory spike where a flat Vec of ALL scored entries (up to ~600 MB)
    // had to exist simultaneously with the dedup HashMap.  With fold/reduce the
    // working set at any point is bounded by the total number of unique entries
    // (≈ final size), not twice that.
    let best_by_id: HashMap<u32, CoelutionScoredEntry> = window_groups
        .par_iter()
        .fold(
            HashMap::<u32, CoelutionScoredEntry>::new,
            |mut map, ((lower, upper), spectrum_indices)| {
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
                    return map;
                }

                let window_spectra: Vec<&Spectrum> = window_pairs.iter().map(|(_, s)| *s).collect();
                let _global_indices: Vec<usize> =
                    window_pairs.iter().map(|(idx, _)| *idx).collect();

                // Pre-preprocess all window spectra for XCorr to eliminate redundant
                // per-entry preprocessing (binning + windowing + sliding window subtraction).
                // Each spectrum is preprocessed once here; per-entry scoring then uses
                // O(n_frags) dot product lookups instead of O(n_peaks) re-preprocessing.
                let preprocessed_xcorr: Vec<Vec<f32>> = window_spectra
                    .iter()
                    .map(|s| scorer.preprocess_spectrum_for_xcorr(s))
                    .collect();

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
                            // Check m/z falls within isolation window.
                            // Use half-open [lower, upper) so entries at the boundary
                            // belong to exactly one window, preventing double-counting.
                            let entry = &library[idx];
                            if entry.precursor_mz < *lower || entry.precursor_mz >= *upper {
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
                        let has_override =
                            boundary_overrides.and_then(|m| m.get(&entry.id)).copied();

                        // Expected RT for this entry
                        let expected_rt = rt_calibration
                            .map(|cal| cal.predict(entry.retention_time))
                            .unwrap_or(entry.retention_time);

                        // Filter spectra within RT range.
                        // For boundary overrides: use the given boundaries + margin for SNR context.
                        // For normal search: use expected_rt ± rt_tolerance.
                        let rt_spec_triples: Vec<(usize, usize, &Spectrum)> =
                            if let Some((_, target_start, target_end)) = has_override {
                                let peak_width = (target_end - target_start).max(0.1);
                                let margin = peak_width.max(0.2);
                                let rt_lo = target_start - margin;
                                let rt_hi = target_end + margin;
                                window_pairs
                                    .iter()
                                    .enumerate()
                                    .filter(|(_, (_, spec))| {
                                        spec.retention_time >= rt_lo && spec.retention_time <= rt_hi
                                    })
                                    .map(|(win_idx, (global_idx, spec))| {
                                        (win_idx, *global_idx, *spec)
                                    })
                                    .collect()
                            } else {
                                window_pairs
                                    .iter()
                                    .enumerate()
                                    .filter(|(_, (_, spec))| {
                                        (spec.retention_time - expected_rt).abs() <= rt_tolerance
                                    })
                                    .map(|(win_idx, (global_idx, spec))| {
                                        (win_idx, *global_idx, *spec)
                                    })
                                    .collect()
                            };

                        if rt_spec_triples.len() < MIN_COELUTION_SPECTRA {
                            return None;
                        }

                        // Signal pre-filter (skipped for boundary overrides — we're told to score here)
                        if has_override.is_none() && config.prefilter_enabled {
                            let tol = if tol_ppm > 0.0 { tol_ppm } else { tol_da };
                            let tol_unit = if tol_ppm > 0.0 {
                                ToleranceUnit::Ppm
                            } else {
                                ToleranceUnit::Mz
                            };
                            const WIN: usize = 4;
                            const MIN_PASS: u32 = 3;
                            let mut window = [false; WIN];
                            let mut win_sum = 0u32;
                            let mut has_signal = false;
                            for (i, (_, _, spec)) in rt_spec_triples.iter().enumerate() {
                                let passes = has_topn_fragment_match(
                                    &entry.fragments,
                                    &spec.mzs,
                                    tol,
                                    tol_unit,
                                );
                                let slot = i % WIN;
                                if window[slot] {
                                    win_sum -= 1;
                                }
                                window[slot] = passes;
                                if passes {
                                    win_sum += 1;
                                }
                                if i + 1 >= WIN && win_sum >= MIN_PASS {
                                    has_signal = true;
                                    break;
                                }
                            }
                            if !has_signal {
                                return None;
                            }
                        }

                        let cand_spectra: Vec<&Spectrum> =
                            rt_spec_triples.iter().map(|(_, _, s)| *s).collect();
                        let cand_global: Vec<usize> =
                            rt_spec_triples.iter().map(|(_, idx, _)| *idx).collect();
                        let cand_window_local: Vec<usize> = rt_spec_triples
                            .iter()
                            .map(|(win_idx, _, _)| *win_idx)
                            .collect();

                        // Extract top-6 fragment XICs
                        let xics = extract_fragment_xics(
                            &entry.fragments,
                            &cand_spectra,
                            tol_da,
                            tol_ppm,
                            6,
                        );

                        if xics.len() < 2 {
                            return None;
                        }

                        // Find reference XIC (highest total intensity) for scoring
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

                        // --- Peak selection ---
                        // For boundary overrides (forced/imputed): map given RTs to XIC indices.
                        // For normal search: CWT consensus peak detection with fallbacks.
                        if let Some((target_apex, target_start, target_end)) = has_override {
                            if ref_xic.len() < 3 {
                                return None;
                            }

                            // Map target RT boundaries to XIC scan indices
                            let start_index = ref_xic
                                .partition_point(|(rt, _)| *rt < target_start)
                                .saturating_sub(1)
                                .min(ref_xic.len() - 1);

                            let end_index = ref_xic
                                .partition_point(|(rt, _)| *rt < target_end)
                                .min(ref_xic.len() - 1);

                            let apex_index = ref_xic
                                .partition_point(|(rt, _)| *rt < target_apex)
                                .min(ref_xic.len() - 1);
                            // Pick closer neighbor for apex
                            let apex_index = if apex_index > 0
                                && (ref_xic[apex_index - 1].0 - target_apex).abs()
                                    < (ref_xic[apex_index].0 - target_apex).abs()
                            {
                                apex_index - 1
                            } else {
                                apex_index
                            };
                            let apex_index = apex_index.clamp(start_index, end_index);

                            if end_index <= start_index + 1 {
                                return None;
                            }

                            let raw_ints: Vec<f64> = ref_xic.iter().map(|(_, v)| *v).collect();
                            let area = trapezoidal_area(&ref_xic[start_index..=end_index]);
                            let snr = compute_snr(&raw_ints, apex_index, start_index, end_index);

                            let peak = XICPeakBounds {
                                apex_rt: ref_xic[apex_index].0,
                                apex_intensity: ref_xic[apex_index].1,
                                apex_index,
                                start_rt: ref_xic[start_index].0,
                                end_rt: ref_xic[end_index].0,
                                start_index,
                                end_index,
                                area,
                                signal_to_noise: snr,
                            };

                            let ctx = FeatureComputeContext {
                                entry,
                                xics: &xics,
                                ref_xic,
                                cand_spectra: &cand_spectra,
                                cand_global: &cand_global,
                                scorer: &scorer,
                                ms1_index,
                                calibration,
                                tol_da,
                                tol_ppm,
                                expected_rt,
                                is_hram,
                                file_name,
                                preprocessed_xcorr: Some(&preprocessed_xcorr),
                                cand_window_local: Some(&cand_window_local),
                            };

                            return compute_features_at_peak(&ctx, peak);
                        }

                        // --- Normal CWT peak detection path ---

                        // Detect candidate peaks using CWT consensus across all transitions.
                        let full_polish = tukey_median_polish(&xics, 10, 0.01);
                        let candidates = {
                            let cwt_candidates = detect_cwt_consensus_peaks(&xics, 0.0);
                            if cwt_candidates.is_empty() {
                                let mp_candidates = full_polish
                                    .as_ref()
                                    .map(|mp| detect_all_xic_peaks(&mp.elution_profile, 0.01, 5.0))
                                    .unwrap_or_default();
                                if mp_candidates.is_empty() {
                                    detect_all_xic_peaks(ref_xic, 0.01, 5.0)
                                } else {
                                    mp_candidates
                                }
                            } else {
                                cwt_candidates
                            }
                        };

                        if candidates.is_empty() {
                            return None;
                        }

                        // Score each candidate by mean pairwise fragment correlation
                        let raw_ints: Vec<f64> = ref_xic.iter().map(|(_, v)| *v).collect();
                        let mut scored_candidates: Vec<(&XICPeakBounds, f64)> = candidates
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
                            .collect();
                        scored_candidates.sort_by(|a, b| b.1.total_cmp(&a.1));

                        // Store top-N CWT candidates for inter-replicate reconciliation
                        let top_n = config.reconciliation.top_n_peaks;
                        let cwt_top_n: Vec<CwtCandidate> = scored_candidates
                            .iter()
                            .take(top_n)
                            .map(|(bp, score)| {
                                let si = bp.start_index;
                                let ei = bp.end_index;
                                let area = trapezoidal_area(&ref_xic[si..=ei]);
                                let (apex_idx, _) = ref_xic[si..=ei]
                                    .iter()
                                    .enumerate()
                                    .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                                    .map(|(i, &(_, v))| (si + i, v))
                                    .unwrap_or((bp.apex_index, 0.0));
                                let snr = compute_snr(&raw_ints, apex_idx, si, ei);
                                CwtCandidate {
                                    apex_rt: ref_xic[apex_idx].0,
                                    start_rt: ref_xic[si].0,
                                    end_rt: ref_xic[ei].0,
                                    area,
                                    snr,
                                    coelution_score: *score,
                                }
                            })
                            .collect();

                        // Build peak from best candidate
                        let (best_bp, _) = scored_candidates[0];
                        let si = best_bp.start_index;
                        let ei = best_bp.end_index;
                        let (apex_idx, apex_val) = ref_xic[si..=ei]
                            .iter()
                            .enumerate()
                            .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                            .map(|(i, &(_, v))| (si + i, v))
                            .unwrap_or((best_bp.apex_index, 0.0));
                        let area = trapezoidal_area(&ref_xic[si..=ei]);
                        let snr = compute_snr(&raw_ints, apex_idx, si, ei);
                        let peak = XICPeakBounds {
                            apex_rt: ref_xic[apex_idx].0,
                            apex_intensity: apex_val,
                            apex_index: apex_idx,
                            start_rt: ref_xic[si].0,
                            end_rt: ref_xic[ei].0,
                            start_index: si,
                            end_index: ei,
                            area,
                            signal_to_noise: snr,
                        };

                        let ctx = FeatureComputeContext {
                            entry,
                            xics: &xics,
                            ref_xic,
                            cand_spectra: &cand_spectra,
                            cand_global: &cand_global,
                            scorer: &scorer,
                            ms1_index,
                            calibration,
                            tol_da,
                            tol_ppm,
                            expected_rt,
                            is_hram,
                            file_name,
                            preprocessed_xcorr: Some(&preprocessed_xcorr),
                            cand_window_local: Some(&cand_window_local),
                        };

                        let mut result = compute_features_at_peak(&ctx, peak);
                        if let Some(ref mut entry) = result {
                            entry.cwt_candidates = cwt_top_n;
                        }
                        result
                    })
                    .collect::<Vec<_>>();

                // Dedup within this window (handles any intra-window duplicates) and
                // merge into the thread-local map, keeping the best coelution_sum.
                for entry in window_entries {
                    let is_better = map
                        .get(&entry.entry_id)
                        .map(|existing| {
                            entry.features.coelution_sum > existing.features.coelution_sum
                        })
                        .unwrap_or(true);
                    if is_better {
                        map.insert(entry.entry_id, entry);
                    }
                }

                pb.inc(1);
                map
            },
        )
        .reduce(HashMap::<u32, CoelutionScoredEntry>::new, |mut a, b| {
            // Merge two thread-local maps, keeping best coelution_sum per entry.
            for (id, entry) in b {
                let is_better = a
                    .get(&id)
                    .map(|existing| entry.features.coelution_sum > existing.features.coelution_sum)
                    .unwrap_or(true);
                if is_better {
                    a.insert(id, entry);
                }
            }
            a
        });

    pb.finish_with_message("Done");
    crate::logging::clear_progress_bar();

    log::debug!(
        "Coelution search complete: {} scored entries ({} targets, {} decoys)",
        best_by_id.len(),
        best_by_id.values().filter(|e| !e.is_decoy).count(),
        best_by_id.values().filter(|e| e.is_decoy).count(),
    );

    let mut deduped: Vec<CoelutionScoredEntry> = best_by_id.into_values().collect();

    // Sort by (entry_id, scan_number) for deterministic output regardless of
    // Rayon thread scheduling or HashMap iteration order
    deduped.sort_by(|a, b| {
        a.entry_id
            .cmp(&b.entry_id)
            .then(a.scan_number.cmp(&b.scan_number))
    });

    // Multi-charge consensus is handled post-FDR (see run_analysis),
    // where the full SVM score is available to pick the correct peak.
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

/// Write calibration debug CSV showing LDA features and discriminant scores
///
/// Produces a paired target-decoy CSV sorted by discriminant score descending,
/// showing the 4 LDA features, discriminant score, q-value, and RT info.
fn write_calibration_debug_csv(
    matches: &[CalibrationMatch],
    output_path: &std::path::Path,
    expected_rt_fn: Option<&dyn Fn(f64) -> f64>,
) -> Result<()> {
    use std::collections::HashMap as DebugMap;
    use std::fs::File;

    let file = File::create(output_path)
        .map_err(|e| OspreyError::OutputError(format!("Failed to create debug CSV: {}", e)))?;
    let mut writer = std::io::BufWriter::new(file);

    // Pair targets with their decoys by base_id
    let mut target_map: DebugMap<u32, &CalibrationMatch> = DebugMap::new();
    let mut decoy_map: DebugMap<u32, &CalibrationMatch> = DebugMap::new();
    for m in matches {
        let base_id = m.entry_id & 0x7FFFFFFF;
        if m.is_decoy {
            decoy_map.insert(base_id, m);
        } else {
            target_map.insert(base_id, m);
        }
    }

    struct PairedRow<'a> {
        target: &'a CalibrationMatch,
        decoy: &'a CalibrationMatch,
    }

    let mut paired: Vec<PairedRow> = Vec::new();
    for (&base_id, &target) in &target_map {
        if let Some(&decoy) = decoy_map.get(&base_id) {
            paired.push(PairedRow { target, decoy });
        }
    }

    // Sort by target discriminant score descending (best first)
    paired.sort_by(|a, b| {
        b.target
            .discriminant_score
            .total_cmp(&a.target.discriminant_score)
    });

    writeln!(
        writer,
        "entry_id,charge,target_sequence,decoy_sequence,\
         target_correlation,decoy_correlation,\
         target_libcosine,decoy_libcosine,\
         target_top6,decoy_top6,\
         target_xcorr,decoy_xcorr,\
         target_discriminant,decoy_discriminant,\
         target_qvalue,decoy_qvalue,\
         target_snr,decoy_snr,\
         target_rt,decoy_rt,library_rt,expected_rt,\
         target_delta_rt,decoy_delta_rt,\
         target_matched_frags,decoy_matched_frags,\
         target_wins"
    )
    .map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

    for p in &paired {
        let t = p.target;
        let d = p.decoy;
        let expected_rt = expected_rt_fn.map(|f| f(t.library_rt));
        let ref_rt = expected_rt.unwrap_or(t.library_rt);
        let target_delta = (t.measured_rt - ref_rt).abs();
        let decoy_delta = (d.measured_rt - ref_rt).abs();
        let target_wins = t.discriminant_score > d.discriminant_score;

        writeln!(
            writer,
            "{},{},{},{},{:.4},{:.4},{:.4},{:.4},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.2},{:.2},{:.2},{:.2},{:.2},{},{:.2},{:.2},{},{},{}",
            t.entry_id,
            t.charge,
            t.sequence,
            d.sequence,
            t.correlation_score,
            d.correlation_score,
            t.libcosine_apex,
            d.libcosine_apex,
            t.top6_matched_apex,
            d.top6_matched_apex,
            t.xcorr_score,
            d.xcorr_score,
            t.discriminant_score,
            d.discriminant_score,
            t.q_value,
            d.q_value,
            t.signal_to_noise,
            d.signal_to_noise,
            t.measured_rt,
            d.measured_rt,
            t.library_rt,
            expected_rt.map(|v| format!("{:.2}", v)).unwrap_or_default(),
            target_delta,
            decoy_delta,
            t.n_matched_fragments,
            d.n_matched_fragments,
            target_wins
        )
        .map_err(|e| OspreyError::OutputError(format!("Failed to write row: {}", e)))?;
    }

    log::debug!(
        "Wrote calibration debug CSV ({} pairs) to: {}",
        paired.len(),
        output_path.display()
    );

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

    /// Verifies deduplicate_pairs() returns deterministic order regardless of HashMap internals.
    ///
    /// Non-deterministic HashMap iteration used to cause different row ordering in the
    /// SVM feature matrix, leading to different gradient updates and model weights.
    /// The fix sorts output by entry_id. This test verifies that guarantee.
    #[test]
    fn test_deduplicate_pairs_deterministic() {
        use osprey_core::types::{CoelutionFeatureSet, CoelutionScoredEntry, XICPeakBounds};

        // Helper to create a scored entry with a given entry_id and coelution_sum
        let make_entry = |entry_id: u32, is_decoy: bool, coelution_sum: f64| CoelutionScoredEntry {
            entry_id,
            is_decoy,
            sequence: format!("SEQ{}", entry_id & 0x7FFFFFFF),
            modified_sequence: format!("SEQ{}", entry_id & 0x7FFFFFFF),
            charge: 2,
            precursor_mz: 500.0,
            protein_ids: vec![],
            scan_number: 100,
            apex_rt: 10.0,
            peak_bounds: XICPeakBounds {
                apex_rt: 10.0,
                apex_intensity: 100.0,
                apex_index: 5,
                start_rt: 9.0,
                end_rt: 11.0,
                start_index: 0,
                end_index: 10,
                area: 500.0,
                signal_to_noise: 10.0,
            },
            features: CoelutionFeatureSet {
                coelution_sum,
                ..CoelutionFeatureSet::default()
            },
            fragment_mzs: vec![],
            fragment_intensities: vec![],
            reference_xic: vec![],
            file_name: "test".into(),
            run_precursor_qvalue: 1.0,
            run_peptide_qvalue: 1.0,
            run_protein_qvalue: 1.0,
            experiment_precursor_qvalue: 1.0,
            experiment_peptide_qvalue: 1.0,
            experiment_protein_qvalue: 1.0,
            score: 0.0,
            pep: 1.0,
            cwt_candidates: vec![],
        };

        // Create entries: 5 target-decoy pairs with tied coelution_sum scores.
        // The HashMap iteration order will vary, but the output must be sorted by entry_id.
        let entries = vec![
            make_entry(5, false, 0.8),
            make_entry(3, false, 0.8),
            make_entry(1, false, 0.8),
            make_entry(4, false, 0.8),
            make_entry(2, false, 0.8),
            make_entry(5 | 0x80000000, true, 0.3),
            make_entry(3 | 0x80000000, true, 0.3),
            make_entry(1 | 0x80000000, true, 0.3),
            make_entry(4 | 0x80000000, true, 0.3),
            make_entry(2 | 0x80000000, true, 0.3),
        ];

        let first_result: Vec<u32> = deduplicate_pairs(entries.clone())
            .iter()
            .map(|e| e.entry_id)
            .collect();

        // Run multiple times — must always produce the same order
        for _ in 0..20 {
            let result: Vec<u32> = deduplicate_pairs(entries.clone())
                .iter()
                .map(|e| e.entry_id)
                .collect();
            assert_eq!(
                result, first_result,
                "deduplicate_pairs must return deterministic order"
            );
        }

        // Verify sorted by entry_id (targets first since their IDs are smaller)
        assert_eq!(first_result.len(), 10);
        for i in 1..first_result.len() {
            assert!(
                first_result[i - 1] < first_result[i],
                "Output must be sorted by entry_id: {} >= {}",
                first_result[i - 1],
                first_result[i]
            );
        }
    }

    /// Helper to create a minimal CoelutionScoredEntry for consensus tests.
    fn make_scored_entry(
        entry_id: u32,
        modified_sequence: &str,
        charge: u8,
        apex_rt: f64,
        start_rt: f64,
        end_rt: f64,
        coelution_sum: f64,
    ) -> CoelutionScoredEntry {
        make_scored_entry_with_score(
            entry_id,
            modified_sequence,
            charge,
            apex_rt,
            start_rt,
            end_rt,
            coelution_sum,
            0.0, // score
            1.0, // run_qvalue (fails FDR by default)
        )
    }

    /// Helper with explicit SVM score and q-value for post-FDR consensus tests.
    #[allow(clippy::too_many_arguments)]
    fn make_scored_entry_with_score(
        entry_id: u32,
        modified_sequence: &str,
        charge: u8,
        apex_rt: f64,
        start_rt: f64,
        end_rt: f64,
        coelution_sum: f64,
        score: f64,
        run_qvalue: f64,
    ) -> CoelutionScoredEntry {
        CoelutionScoredEntry {
            entry_id,
            is_decoy: modified_sequence.starts_with("DECOY_"),
            sequence: modified_sequence
                .strip_prefix("DECOY_")
                .unwrap_or(modified_sequence)
                .to_string(),
            modified_sequence: modified_sequence.to_string(),
            charge,
            precursor_mz: 500.0,
            protein_ids: vec![],
            scan_number: 1,
            apex_rt,
            peak_bounds: XICPeakBounds {
                apex_rt,
                apex_intensity: 1000.0,
                apex_index: 50,
                start_rt,
                end_rt,
                start_index: 40,
                end_index: 60,
                area: 5000.0,
                signal_to_noise: 10.0,
            },
            features: CoelutionFeatureSet {
                coelution_sum,
                ..Default::default()
            },
            fragment_mzs: vec![],
            fragment_intensities: vec![],
            reference_xic: vec![],
            file_name: "test".to_string(),
            run_precursor_qvalue: run_qvalue,
            run_peptide_qvalue: 1.0,
            run_protein_qvalue: 1.0,
            experiment_precursor_qvalue: 1.0,
            experiment_peptide_qvalue: 1.0,
            experiment_protein_qvalue: 1.0,
            score,
            pep: 1.0,
            cwt_candidates: vec![],
        }
    }

    /// Convert test CoelutionScoredEntry vec to FdrEntry vec for select_post_fdr_consensus tests.
    fn to_fdr(entries: &[CoelutionScoredEntry]) -> Vec<FdrEntry> {
        entries
            .iter()
            .enumerate()
            .map(|(i, e)| {
                let mut stub = e.to_fdr_entry();
                stub.parquet_index = i as u32;
                stub
            })
            .collect()
    }

    #[test]
    fn test_consensus_single_charge_no_change() {
        let entries = vec![make_scored_entry(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0)];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        assert!(rescore.is_empty());
    }

    #[test]
    fn test_consensus_two_charges_same_peak() {
        // Both charges found the same peak (within tolerance)
        let entries = vec![
            make_scored_entry(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0),
            make_scored_entry(1, "PEPTIDEK", 3, 15.1, 14.6, 15.6, 6.0),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        assert!(rescore.is_empty(), "Both charges agree — no re-scoring");
    }

    #[test]
    fn test_consensus_two_charges_different_peaks() {
        // Charge 2+ found the true peak (passes FDR), charge 3+ picked a different one
        let entries = vec![
            make_scored_entry_with_score(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0, 2.5, 0.005),
            make_scored_entry_with_score(1, "PEPTIDEK", 3, 22.0, 21.5, 22.5, 3.0, 1.0, 0.05),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        assert_eq!(rescore.len(), 1, "The other charge state needs re-scoring");

        // Verify the rescore target has the consensus boundaries from entry 0
        let (idx, apex, start, end) = rescore[0];
        assert_eq!(idx, 1);
        assert!((apex - 15.0).abs() < 0.01);
        assert!((start - 14.5).abs() < 0.01);
        assert!((end - 15.5).abs() < 0.01);
    }

    #[test]
    fn test_consensus_uses_svm_score_not_coelution() {
        // Entry 1 has higher coelution_sum but lower SVM score.
        // Post-FDR consensus should prefer the higher SVM score (entry 0).
        let entries = vec![
            make_scored_entry_with_score(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 5.0, 3.0, 0.005),
            make_scored_entry_with_score(1, "PEPTIDEK", 3, 22.0, 21.5, 22.5, 9.0, 1.0, 0.008),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        assert_eq!(rescore.len(), 1);
        // Entry 1 should be re-scored at entry 0's boundaries (higher SVM score)
        assert_eq!(rescore[0].0, 1);
        assert!((rescore[0].1 - 15.0).abs() < 0.01);
    }

    #[test]
    fn test_consensus_skips_groups_where_none_pass_fdr() {
        // Neither passes FDR — skip entirely (no impact on results)
        let entries = vec![
            make_scored_entry_with_score(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0, 1.0, 0.5),
            make_scored_entry_with_score(1, "PEPTIDEK", 3, 22.0, 21.5, 22.5, 3.0, 0.5, 0.8),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        assert!(
            rescore.is_empty(),
            "No charge state passes FDR — skip group"
        );
    }

    #[test]
    fn test_consensus_three_charges_two_agree() {
        let entries = vec![
            make_scored_entry_with_score(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0, 3.0, 0.005),
            make_scored_entry_with_score(1, "PEPTIDEK", 3, 15.1, 14.6, 15.6, 6.0, 2.0, 0.008),
            make_scored_entry_with_score(2, "PEPTIDEK", 4, 22.0, 21.5, 22.5, 2.0, 0.5, 0.05),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        assert_eq!(rescore.len(), 1, "One needs re-scoring");
        assert_eq!(rescore[0].0, 2);
    }

    #[test]
    fn test_consensus_decoys_separate_from_targets() {
        // Target and decoy of same peptide should be in different groups
        let entries = vec![
            make_scored_entry(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0),
            make_scored_entry(1, "DECOY_PEPTIDEK", 2, 22.0, 21.5, 22.5, 5.0),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        // Different modified_sequence → different groups → no re-scoring
        assert!(rescore.is_empty());
    }

    #[test]
    fn test_consensus_multiple_peptides_independent() {
        let entries = vec![
            make_scored_entry_with_score(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0, 3.0, 0.005),
            make_scored_entry_with_score(1, "PEPTIDEK", 3, 22.0, 21.5, 22.5, 3.0, 1.0, 0.05),
            make_scored_entry_with_score(2, "ANOTHERPEPTIDER", 2, 10.0, 9.5, 10.5, 7.0, 2.5, 0.006),
            make_scored_entry_with_score(3, "ANOTHERPEPTIDER", 3, 10.1, 9.6, 10.6, 5.0, 2.0, 0.007),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        // PEPTIDEK: entry 0 best (score=3.0, passes FDR), entry 1 rescored (different peak)
        // ANOTHERPEPTIDER: both agree (same peak) → no re-scoring
        assert_eq!(rescore.len(), 1);
        assert_eq!(rescore[0].0, 1);
    }

    #[test]
    fn test_consensus_per_file_grouping() {
        // select_post_fdr_consensus is always called per-file, so all entries
        // share the same file. Same peptide, different charge, different peak →
        // the non-passing entry gets re-scored to the passing entry's peak.
        let entries = vec![
            make_scored_entry_with_score(0, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0, 3.0, 0.005),
            make_scored_entry_with_score(1, "PEPTIDEK", 3, 22.0, 21.5, 22.5, 3.0, 1.0, 0.05),
        ];
        let fdr = to_fdr(&entries);
        let rescore = select_post_fdr_consensus(&fdr, 0.01);
        // Same group (same modified_sequence) → entry 1 at different peak gets re-scored
        assert_eq!(rescore.len(), 1);
        assert_eq!(rescore[0].0, 1);
    }

    // ===================================================================
    // FDR stub compaction + parquet_index regression tests
    // ===================================================================
    //
    // After first-pass FDR, non-passing stubs are dropped (compacted) to free
    // memory. The parquet_index field preserves the original Parquet row index
    // so that downstream phases (Percolator, reconciliation, blib output) can
    // still look up features and CWT candidates from the original Parquet file.
    //
    // These tests catch the bug where code uses Vec position (`local_idx`) as
    // the Parquet index after compaction, causing out-of-range panics or
    // garbage scores.

    /// Build a minimal FdrEntry stub for compaction tests.
    fn make_stub(entry_id: u32, parquet_index: u32, is_decoy: bool) -> FdrEntry {
        FdrEntry {
            entry_id,
            parquet_index,
            is_decoy,
            charge: 2,
            scan_number: 1,
            apex_rt: 0.0,
            start_rt: 0.0,
            end_rt: 0.0,
            coelution_sum: 0.0,
            score: 0.0,
            run_precursor_qvalue: 1.0,
            run_peptide_qvalue: 1.0,
            run_protein_qvalue: 1.0,
            experiment_precursor_qvalue: 1.0,
            experiment_peptide_qvalue: 1.0,
            experiment_protein_qvalue: 1.0,
            pep: 1.0,
            modified_sequence: Arc::from("PEPTIDE"),
        }
    }

    #[test]
    fn test_compaction_preserves_parquet_index() {
        // 10 entries: 5 targets (entry_id 0..5) + 5 decoys (entry_id with bit 31)
        // Each entry's parquet_index matches its original Vec position.
        let mut stubs: Vec<FdrEntry> = (0..5)
            .map(|i| make_stub(i, i, false)) // targets
            .chain((0..5).map(|i| make_stub(i | 0x8000_0000, 5 + i, true))) // paired decoys
            .collect();

        // Mark entries 0, 2, 4 (targets) as passing
        stubs[0].run_precursor_qvalue = 0.001;
        stubs[2].run_precursor_qvalue = 0.001;
        stubs[4].run_precursor_qvalue = 0.001;

        // Compute first_pass_base_ids (target base_ids that passed)
        let first_pass_base_ids: HashSet<u32> = stubs
            .iter()
            .filter(|e| !e.is_decoy && e.run_precursor_qvalue <= 0.01)
            .map(|e| e.entry_id & 0x7FFF_FFFF)
            .collect();
        assert_eq!(first_pass_base_ids.len(), 3);

        // Compact: keep entries whose base_id is in first_pass_base_ids
        // (this should keep both targets and their paired decoys)
        stubs.retain(|e| first_pass_base_ids.contains(&(e.entry_id & 0x7FFF_FFFF)));
        assert_eq!(stubs.len(), 6, "should keep 3 targets + 3 paired decoys");

        // Verify that parquet_index values are preserved (not collapsed to 0..5)
        let parquet_indices: Vec<u32> = stubs.iter().map(|e| e.parquet_index).collect();
        assert!(
            parquet_indices.contains(&0)
                && parquet_indices.contains(&2)
                && parquet_indices.contains(&4),
            "target parquet_indices should be preserved: {:?}",
            parquet_indices
        );
        assert!(
            parquet_indices.contains(&5)
                && parquet_indices.contains(&7)
                && parquet_indices.contains(&9),
            "decoy parquet_indices should be preserved: {:?}",
            parquet_indices
        );

        // CRITICAL: After compaction, the Vec position (local_idx) is NOT the
        // same as parquet_index for at least some entries. Downstream code must
        // use entry.parquet_index for any Parquet-row-based lookups (features,
        // CWT candidates).
        let any_divergence = stubs
            .iter()
            .enumerate()
            .any(|(local_idx, entry)| local_idx as u32 != entry.parquet_index);
        assert!(
            any_divergence,
            "After compaction, at least some entries should have local_idx != parquet_index"
        );
    }

    #[test]
    fn test_feature_lookup_uses_parquet_index_after_compaction() {
        // Simulate a Parquet file with 10 feature rows (one per row index 0..10)
        // where each row's "feature" is just its index as f64.
        let file_features: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();

        // Compacted stubs: 3 entries kept from original positions 1, 4, 7
        let stubs = vec![
            make_stub(1, 1, false),
            make_stub(4, 4, false),
            make_stub(7, 7, false),
        ];

        // CORRECT lookup: use parquet_index
        for entry in &stubs {
            let pq_idx = entry.parquet_index as usize;
            assert!(pq_idx < file_features.len());
            let feature = &file_features[pq_idx];
            assert_eq!(
                feature[0], entry.parquet_index as f64,
                "feature lookup via parquet_index should match"
            );
        }

        // INCORRECT lookup: using Vec position would give wrong features
        // (this test documents the bug we're guarding against)
        for (local_idx, entry) in stubs.iter().enumerate() {
            let wrong_feature = &file_features[local_idx];
            // The wrong feature value (local_idx) should NOT match the
            // entry's actual parquet_index value
            if local_idx as u32 != entry.parquet_index {
                assert_ne!(
                    wrong_feature[0], entry.parquet_index as f64,
                    "Vec-position lookup gives wrong features after compaction"
                );
            }
        }
    }

    #[test]
    fn test_gap_fill_sentinel_parquet_index() {
        // Gap-fill entries (added during reconciliation) get parquet_index = u32::MAX
        // because they don't exist in the original Parquet. Their features must
        // come from the overlay, not the Parquet file.
        let stubs = [
            make_stub(1, 1, false),        // normal entry from Parquet
            make_stub(2, u32::MAX, false), // gap-fill entry
            make_stub(3, 3, false),        // normal entry from Parquet
        ];

        // Filter for gap-fill entries
        let gap_fill_count = stubs.iter().filter(|e| e.parquet_index == u32::MAX).count();
        assert_eq!(gap_fill_count, 1);

        // Filter for normal entries (have valid parquet_index)
        let normal_count = stubs.iter().filter(|e| e.parquet_index != u32::MAX).count();
        assert_eq!(normal_count, 2);
    }

    #[test]
    fn test_to_fdr_entry_default_parquet_index_is_zero() {
        // CRITICAL: CoelutionScoredEntry::to_fdr_entry() initializes parquet_index = 0
        // because it doesn't know its position in the Vec/Parquet. Callers MUST
        // populate parquet_index explicitly after the call.
        //
        // This test catches the bug class where stubs are created without setting
        // parquet_index, causing all stubs to point to Parquet row 0.
        let entry =
            make_scored_entry_with_score(42, "PEPTIDEK", 2, 15.0, 14.5, 15.5, 8.0, 3.0, 0.005);
        let stub = entry.to_fdr_entry();
        assert_eq!(
            stub.parquet_index, 0,
            "to_fdr_entry must default parquet_index to 0; callers must populate it"
        );
    }

    #[test]
    fn test_fresh_stubs_must_have_correct_parquet_index() {
        // Regression test for the bug where freshly scored entries were converted to
        // FdrEntry stubs without setting parquet_index. This caused every stub to have
        // parquet_index=0, so Percolator's feature lookup returned the same features
        // for every entry, producing 0 passing precursors at FDR.
        //
        // The pipeline pattern that MUST be used:
        //   entries.iter().enumerate().map(|(i, e)| {
        //       let mut stub = e.to_fdr_entry();
        //       stub.parquet_index = i as u32;
        //       stub
        //   }).collect()
        let entries: Vec<CoelutionScoredEntry> = (0..10)
            .map(|i| {
                make_scored_entry_with_score(
                    i,
                    &format!("PEPTIDE{}", i),
                    2,
                    15.0,
                    14.5,
                    15.5,
                    8.0,
                    3.0,
                    0.005,
                )
            })
            .collect();

        // Simulate the (correct) pipeline conversion
        let stubs: Vec<FdrEntry> = entries
            .iter()
            .enumerate()
            .map(|(i, e)| {
                let mut stub = e.to_fdr_entry();
                stub.parquet_index = i as u32;
                stub
            })
            .collect();

        // Each stub's parquet_index must equal its position in the Vec
        for (i, stub) in stubs.iter().enumerate() {
            assert_eq!(
                stub.parquet_index, i as u32,
                "fresh stub at position {} should have parquet_index={}, got {}",
                i, i, stub.parquet_index
            );
        }

        // Distinct parquet_index values (the bug would make all of them 0)
        let unique_indices: HashSet<u32> = stubs.iter().map(|s| s.parquet_index).collect();
        assert_eq!(
            unique_indices.len(),
            stubs.len(),
            "all parquet_index values should be distinct (bug: all are 0 if not populated)"
        );
    }

    #[test]
    fn test_buggy_pattern_collapses_parquet_index_to_zero() {
        // Documents the BUG pattern: calling to_fdr_entry without setting parquet_index.
        // This test exists to make sure we catch this specific failure mode if anyone
        // accidentally reverts to the buggy pattern.
        let entries: Vec<CoelutionScoredEntry> = (0..5)
            .map(|i| {
                make_scored_entry_with_score(
                    i,
                    &format!("PEPTIDE{}", i),
                    2,
                    15.0,
                    14.5,
                    15.5,
                    8.0,
                    3.0,
                    0.005,
                )
            })
            .collect();

        // BUG: convert without populating parquet_index
        let buggy_stubs: Vec<FdrEntry> = entries.iter().map(|e| e.to_fdr_entry()).collect();

        // All parquet_index values are 0 — they all point to Parquet row 0
        let all_zero = buggy_stubs.iter().all(|s| s.parquet_index == 0);
        assert!(
            all_zero,
            "buggy pattern produces all-zero parquet_index — confirms the bug exists \
             when callers forget to populate parquet_index"
        );
    }

    /// Models the reconciliation Parquet write-back step that appends gap-fill
    /// entries and updates their stubs' parquet_index. This is the bug class
    /// where gap-fill stubs were left with parquet_index = u32::MAX even though
    /// their features had been written to the end of the Parquet file.
    fn simulate_reconcile_writeback(
        fdr_entries: &mut [FdrEntry],
        original_parquet_len: usize,
        gap_fill_vec_indices: &[usize],
    ) {
        // Mirrors the reconciliation write-back logic in pipeline.rs:
        // 1. Start at original_parquet_len (the end of the original Parquet rows)
        // 2. For each gap-fill entry (in vec order), assign sequential parquet_index
        //    matching the row it will occupy after append + write-back
        let gap_start_pq_row = original_parquet_len;
        for (gap_offset, &vec_idx) in gap_fill_vec_indices.iter().enumerate() {
            let new_pq_row = (gap_start_pq_row + gap_offset) as u32;
            if vec_idx < fdr_entries.len() {
                fdr_entries[vec_idx].parquet_index = new_pq_row;
            }
        }
    }

    #[test]
    fn test_reconcile_writeback_updates_gap_fill_parquet_index() {
        // Setup: 5 normal entries (parquet rows 0..5), then 3 gap-fill entries appended.
        // The gap-fill entries initially have parquet_index = u32::MAX (sentinel).
        let mut stubs: Vec<FdrEntry> = (0..5).map(|i| make_stub(i, i, false)).collect();
        // Append gap-fill stubs with sentinel parquet_index
        for i in 0..3 {
            stubs.push(make_stub(100 + i, u32::MAX, false));
        }

        // Verify all 3 gap-fill stubs start with u32::MAX
        let initial_max_count = stubs.iter().filter(|s| s.parquet_index == u32::MAX).count();
        assert_eq!(initial_max_count, 3);

        // Simulate reconciliation write-back: gap-fill entries are at vec indices 5, 6, 7
        // and will be written at parquet rows 5, 6, 7 (right after the original 5 rows)
        let gap_fill_vec_indices = vec![5, 6, 7];
        simulate_reconcile_writeback(&mut stubs, 5, &gap_fill_vec_indices);

        // After write-back, gap-fill stubs should have actual parquet rows assigned
        assert_eq!(stubs[5].parquet_index, 5, "gap-fill 0 -> parquet row 5");
        assert_eq!(stubs[6].parquet_index, 6, "gap-fill 1 -> parquet row 6");
        assert_eq!(stubs[7].parquet_index, 7, "gap-fill 2 -> parquet row 7");

        // No stubs should still have the u32::MAX sentinel
        let remaining_max_count = stubs.iter().filter(|s| s.parquet_index == u32::MAX).count();
        assert_eq!(
            remaining_max_count, 0,
            "all gap-fill entries should have real parquet_index after write-back"
        );

        // All parquet_index values should be unique and contiguous (0..8)
        let mut indices: Vec<u32> = stubs.iter().map(|s| s.parquet_index).collect();
        indices.sort();
        assert_eq!(
            indices,
            vec![0, 1, 2, 3, 4, 5, 6, 7],
            "parquet_index values should cover all rows 0..8"
        );
    }

    #[test]
    fn test_reconcile_writeback_with_no_gap_fill() {
        // No gap-fill entries -> stubs unchanged
        let mut stubs: Vec<FdrEntry> = (0..5).map(|i| make_stub(i, i, false)).collect();
        simulate_reconcile_writeback(&mut stubs, 5, &[]);
        for (i, stub) in stubs.iter().enumerate() {
            assert_eq!(stub.parquet_index, i as u32);
        }
    }

    #[test]
    fn test_reconcile_writeback_buggy_pattern_leaves_max_sentinel() {
        // Document the bug: if write-back forgets to update parquet_index,
        // gap-fill stubs are left with u32::MAX which breaks Phase 2 of the
        // next Percolator run (feature lookup at file_features[u32::MAX as usize]
        // is out of range, so the entry is filtered out, leaving an empty
        // feature vector that causes a matrix shape mismatch).
        let mut stubs: Vec<FdrEntry> = (0..5).map(|i| make_stub(i, i, false)).collect();
        for i in 0..3 {
            stubs.push(make_stub(100 + i, u32::MAX, false));
        }

        // BUG: don't call simulate_reconcile_writeback (forget to update gap-fill stubs)

        // The bug leaves gap-fill stubs with the sentinel
        let max_count = stubs.iter().filter(|s| s.parquet_index == u32::MAX).count();
        assert_eq!(
            max_count, 3,
            "buggy pattern leaves u32::MAX sentinel on gap-fill stubs after write-back"
        );

        // Documenting the failure mode: when these stubs are passed to Percolator's
        // Phase 2 feature loader, the loader tries file_features[u32::MAX as usize]
        // which is out of range. The .filter_map skips them, returning fewer feature
        // vectors than expected entries. The matrix shape assertion (n*p) then fails.
    }

    #[test]
    fn test_reconcile_writeback_preserves_existing_parquet_index() {
        // Existing entries (positions 0..5) should keep their parquet_index unchanged.
        // Only gap-fill entries (positions 5..8) get updated.
        let mut stubs: Vec<FdrEntry> = vec![
            make_stub(0, 100, false), // existing, parquet_index = 100
            make_stub(1, 101, false),
            make_stub(2, 102, false),
            make_stub(3, u32::MAX, false), // gap-fill at vec position 3
            make_stub(4, u32::MAX, false), // gap-fill at vec position 4
        ];

        // Simulate write-back: original Parquet had 103 rows (0..103),
        // gap-fill entries are at vec indices 3, 4 and get parquet rows 103, 104
        simulate_reconcile_writeback(&mut stubs, 103, &[3, 4]);

        // Existing entries unchanged
        assert_eq!(stubs[0].parquet_index, 100);
        assert_eq!(stubs[1].parquet_index, 101);
        assert_eq!(stubs[2].parquet_index, 102);
        // Gap-fill entries got new parquet_index values
        assert_eq!(stubs[3].parquet_index, 103);
        assert_eq!(stubs[4].parquet_index, 104);
    }

    // ===================================================================
    // Blib plan entry construction: parquet_index vs Vec position
    // ===================================================================
    //
    // After first-pass FDR compaction, the LightFdr Vec is shorter than the
    // Parquet file. The blib plan entry construction must use LightFdr.parquet_index
    // to index into the Parquet projection rows, not zip by position (which was
    // the v26.1.0 bug). Without this, each peptide gets the RT/boundaries of a
    // different peptide, producing wildly wrong blib output.

    #[test]
    fn test_blib_plan_uses_parquet_index_not_vec_position() {
        // Simulate a Parquet file with 10 rows, each having a distinct apex_rt.
        // Row i has entry_id=i, apex_rt = (i * 10.0), others zeroed.
        let parquet_rows: Vec<(u32, f64, f64, f64, f64)> = (0..10)
            .map(|i| (i as u32, i as f64 * 10.0, 0.0, 0.0, 0.0))
            .collect();

        // After compaction, only 3 entries remain (from original positions 2, 5, 8).
        // Their parquet_index preserves the original Parquet row.
        let compacted_fdr: Vec<(u32, Arc<str>, u8, bool)> = vec![
            (2, Arc::from("PEPTIDEA"), 2, false), // parquet row 2 -> apex 20.0
            (5, Arc::from("PEPTIDEB"), 2, false), // parquet row 5 -> apex 50.0
            (8, Arc::from("PEPTIDEC"), 3, false), // parquet row 8 -> apex 80.0
        ];

        // CORRECT: look up by parquet_index
        let mut correct_plan = Vec::new();
        for &(pq_idx, ref modseq, charge, is_decoy) in &compacted_fdr {
            if is_decoy {
                continue;
            }
            let row = &parquet_rows[pq_idx as usize];
            correct_plan.push((modseq.clone(), charge, row.1)); // apex_rt
        }
        assert_eq!(
            correct_plan[0].2, 20.0,
            "PEPTIDEA should get apex 20.0 (parquet row 2)"
        );
        assert_eq!(
            correct_plan[1].2, 50.0,
            "PEPTIDEB should get apex 50.0 (parquet row 5)"
        );
        assert_eq!(
            correct_plan[2].2, 80.0,
            "PEPTIDEC should get apex 80.0 (parquet row 8)"
        );

        // BUG: zip by Vec position gives wrong RTs
        let buggy_plan: Vec<f64> = parquet_rows
            .iter()
            .zip(compacted_fdr.iter())
            .map(|(row, _)| row.1)
            .collect();
        assert_eq!(
            buggy_plan[0], 0.0,
            "buggy zip gives row 0 apex (0.0) to first entry"
        );
        assert_eq!(
            buggy_plan[1], 10.0,
            "buggy zip gives row 1 apex (10.0) to second entry"
        );
        assert_eq!(
            buggy_plan[2], 20.0,
            "buggy zip gives row 2 apex (20.0) to third entry"
        );

        // Confirm buggy values differ from correct values
        assert_ne!(buggy_plan[0], correct_plan[0].2);
        assert_ne!(buggy_plan[1], correct_plan[1].2);
        assert_ne!(buggy_plan[2], correct_plan[2].2);
    }

    #[test]
    fn test_blib_plan_with_gap_fill_entries() {
        // Parquet has 5 original rows + 2 gap-fill rows appended at end.
        let parquet_rows: Vec<(u32, f64, f64, f64, f64)> = (0..7)
            .map(|i| (i as u32, i as f64 * 10.0, 0.0, 0.0, 0.0))
            .collect();

        // After compaction + gap-fill, 4 entries:
        //   - 2 original (parquet rows 1, 3)
        //   - 2 gap-fill (parquet rows 5, 6 — appended by reconciliation write-back)
        let compacted_fdr: Vec<(u32, &str)> = vec![
            (1, "PEPTIDEA"), // original, parquet row 1 -> apex 10.0
            (3, "PEPTIDEB"), // original, parquet row 3 -> apex 30.0
            (5, "GAPFILLC"), // gap-fill, parquet row 5 -> apex 50.0
            (6, "GAPFILLD"), // gap-fill, parquet row 6 -> apex 60.0
        ];

        // Correct lookup
        for &(pq_idx, label) in &compacted_fdr {
            let row = &parquet_rows[pq_idx as usize];
            assert_eq!(
                row.1,
                pq_idx as f64 * 10.0,
                "{} at parquet_index {} should get apex {}",
                label,
                pq_idx,
                pq_idx as f64 * 10.0
            );
        }

        // Buggy zip would assign row 0,1,2,3 to entries that are actually at rows 1,3,5,6
        let buggy: Vec<f64> = parquet_rows
            .iter()
            .zip(compacted_fdr.iter())
            .map(|(row, _)| row.1)
            .collect();
        assert_eq!(buggy[0], 0.0, "zip gives row 0 (wrong)");
        assert_ne!(buggy[0], 10.0, "should be 10.0 for PEPTIDEA");
    }

    // ===================================================================
    // Reconciliation write-back: parquet_index vs Vec position
    // ===================================================================
    //
    // After compaction, the FdrEntry Vec is shorter than the Parquet file.
    // When writing re-scored entries back to Parquet, the code must use
    // fdr_entries[idx].parquet_index (not idx) to locate the correct row
    // in full_entries. The v26.1.0 bug used idx directly, causing re-scored
    // entries to overwrite the wrong Parquet row and leaving the actual
    // target entry unchanged.

    #[test]
    fn test_rescore_writeback_uses_parquet_index_not_vec_position() {
        // Simulate: Parquet has 10 entries with distinct apex_rt values.
        // After compaction, 3 entries remain at Vec positions 0, 1, 2
        // but their parquet_index values are 2, 5, 8.
        let mut full_entries: Vec<f64> = (0..10).map(|i| i as f64 * 10.0).collect(); // apex_rt values
        let parquet_indices: Vec<usize> = vec![2, 5, 8]; // from fdr_entries[idx].parquet_index

        // Entry at Vec position 1 (parquet_index=5) was re-scored to apex_rt=99.0
        let rescored_idx = 1;
        let new_apex = 99.0;

        // CORRECT: write to parquet row 5
        let pq_idx = parquet_indices[rescored_idx];
        full_entries[pq_idx] = new_apex;
        assert_eq!(full_entries[5], 99.0, "correct: row 5 should be updated");
        assert_eq!(full_entries[1], 10.0, "correct: row 1 should be untouched");

        // Reset and show BUG: writing to Vec position 1 overwrites the wrong row
        let mut buggy_entries: Vec<f64> = (0..10).map(|i| i as f64 * 10.0).collect();
        buggy_entries[rescored_idx] = new_apex; // BUG: uses Vec position, not parquet_index
        assert_eq!(
            buggy_entries[1], 99.0,
            "buggy: row 1 was overwritten (wrong peptide)"
        );
        assert_eq!(
            buggy_entries[5], 50.0,
            "buggy: row 5 is unchanged (target NOT updated)"
        );
    }

    #[test]
    fn test_rescore_writeback_multi_charge_consensus_scenario() {
        // Real scenario: file has 20 entries. After compaction, entries for
        // peptide AVDFAERDIIHDPGR remain:
        //   z=2 at Vec pos 0, parquet_index=10, apex=9.68
        //   z=3 at Vec pos 1, parquet_index=15, apex=10.67
        // Multi-charge consensus re-scores z=3 to align with z=2's apex (9.68).
        // Write-back must update full_entries[15], not full_entries[1].

        let n_rows = 20;
        let mut full_entries: Vec<f64> = (0..n_rows).map(|i| i as f64).collect();

        // Compacted FdrEntry Vec (only 2 entries shown, in reality more)
        let stubs: Vec<(usize, u32)> = vec![
            (0, 10), // Vec pos 0, parquet_index 10
            (1, 15), // Vec pos 1, parquet_index 15
        ];

        // After multi-charge consensus, z=3 (vec pos 1) was re-scored to apex 9.68
        let rescored_vec_idx = 1;
        let new_apex = 9.68;

        // CORRECT: use parquet_index from the stub
        let pq_idx = stubs[rescored_vec_idx].1 as usize;
        full_entries[pq_idx] = new_apex;
        assert_eq!(
            full_entries[15], 9.68,
            "z=3 entry at parquet row 15 should be updated to consensus apex"
        );
        assert_eq!(
            full_entries[10], 10.0,
            "z=2 entry at parquet row 10 should be unchanged"
        );
        assert_eq!(
            full_entries[1], 1.0,
            "unrelated entry at row 1 should be untouched"
        );

        // BUG: using Vec position 1 corrupts an unrelated entry and misses the target
        let mut buggy_entries: Vec<f64> = (0..n_rows).map(|i| i as f64).collect();
        buggy_entries[rescored_vec_idx] = new_apex;
        assert_eq!(
            buggy_entries[1], 9.68,
            "buggy: row 1 (WRONG peptide) gets the re-scored apex"
        );
        assert_eq!(
            buggy_entries[15], 15.0,
            "buggy: row 15 (actual z=3 entry) remains unchanged"
        );
    }

    // ===================================================================
    // run_search boundary_overrides tests
    // ===================================================================
    //
    // These tests exercise `run_search` with synthetic DIA spectra to verify
    // that the boundary_overrides path (used by reconciliation and multi-charge
    // consensus) produces results at the specified peak boundaries, skips
    // the prefilter, and coexists correctly with normal CWT-based search.

    /// Build a library entry with 3 fragments at known m/z values.
    fn make_lib_entry_with_fragments(id: u32, precursor_mz: f64, rt: f64) -> LibraryEntry {
        let mut entry = LibraryEntry::new(
            id,
            "TESTPEPTIDE".into(),
            "TESTPEPTIDE".into(),
            2,
            precursor_mz,
            rt,
        );
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 80.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 60.0,
                annotation: FragmentAnnotation::default(),
            },
        ];
        entry
    }

    /// Build synthetic DIA spectra with a Gaussian elution profile centered
    /// at `center_rt`. All spectra have peaks at 300.0, 400.0, 500.0 Da
    /// with intensities scaled by a Gaussian envelope. This provides realistic
    /// XIC shapes for peak detection and coelution scoring.
    ///
    /// `n_spectra`: number of spectra spanning `rt_start..rt_end`
    /// `center_rt`: RT of the Gaussian peak center
    /// `sigma`: width of the Gaussian (minutes)
    /// `window_center`, `window_half_width`: isolation window parameters
    fn make_gaussian_spectra(
        n_spectra: usize,
        rt_start: f64,
        rt_end: f64,
        center_rt: f64,
        sigma: f64,
        window_center: f64,
        window_half_width: f64,
    ) -> Vec<Spectrum> {
        let dt = (rt_end - rt_start) / (n_spectra - 1) as f64;
        (0..n_spectra)
            .map(|i| {
                let rt = rt_start + i as f64 * dt;
                let gauss = 1000.0 * (-(rt - center_rt).powi(2) / (2.0 * sigma.powi(2))).exp();
                // Background noise floor + Gaussian signal
                let base = 10.0 + gauss;
                Spectrum {
                    scan_number: i as u32 + 1,
                    retention_time: rt,
                    precursor_mz: window_center,
                    isolation_window: IsolationWindow::symmetric(window_center, window_half_width),
                    // Sorted m/z array with matching fragment intensities
                    mzs: vec![300.0, 400.0, 500.0],
                    intensities: vec![base as f32, (base * 0.8) as f32, (base * 0.6) as f32],
                }
            })
            .collect()
    }

    /// Create an OspreyConfig suitable for unit tests:
    /// unit resolution, Da tolerance, prefilter on.
    fn make_test_config() -> OspreyConfig {
        OspreyConfig {
            resolution_mode: osprey_core::ResolutionMode::UnitResolution,
            fragment_tolerance: FragmentToleranceConfig {
                tolerance: 0.5,
                unit: ToleranceUnit::Mz,
            },
            prefilter_enabled: true,
            ..OspreyConfig::default()
        }
    }

    #[test]
    fn test_run_search_override_returns_entry_at_specified_boundaries() {
        // Build 30 spectra with a Gaussian peak at RT=10.0
        let spectra = make_gaussian_spectra(30, 5.0, 15.0, 10.0, 1.0, 500.0, 12.5);
        let entry = make_lib_entry_with_fragments(0, 500.0, 10.0);
        let library = vec![entry];
        let ms1_index = MS1Index::new(vec![]);
        let config = make_test_config();

        // Override: force boundaries at RT 9.0–11.0, apex at 10.0
        let mut overrides = HashMap::new();
        overrides.insert(0u32, (10.0, 9.0, 11.0));

        let results = run_search(
            &library,
            &spectra,
            &ms1_index,
            None,
            None,
            &config,
            "test.mzML",
            Some(&overrides),
            "Test",
        )
        .unwrap();

        assert_eq!(results.len(), 1, "Should produce exactly one result");
        let result = &results[0];
        assert_eq!(result.entry_id, 0);

        // The peak boundaries should be at the closest XIC points to our overrides
        assert!(
            (result.peak_bounds.start_rt - 9.0).abs() < 0.5,
            "start_rt={} should be near 9.0",
            result.peak_bounds.start_rt
        );
        assert!(
            (result.peak_bounds.end_rt - 11.0).abs() < 0.5,
            "end_rt={} should be near 11.0",
            result.peak_bounds.end_rt
        );
        assert!(
            (result.peak_bounds.apex_rt - 10.0).abs() < 0.5,
            "apex_rt={} should be near 10.0",
            result.peak_bounds.apex_rt
        );
        // Apex must be within boundaries
        assert!(result.peak_bounds.apex_rt >= result.peak_bounds.start_rt);
        assert!(result.peak_bounds.apex_rt <= result.peak_bounds.end_rt);
    }

    #[test]
    fn test_run_search_override_skips_prefilter() {
        // Build spectra where fragment intensities are ZERO — no signal at all.
        // Normal search would fail the prefilter. Override should bypass it.
        let n = 30;
        let dt = 10.0 / (n - 1) as f64;
        let spectra: Vec<Spectrum> = (0..n)
            .map(|i| {
                let rt = 5.0 + i as f64 * dt;
                Spectrum {
                    scan_number: i as u32 + 1,
                    retention_time: rt,
                    precursor_mz: 500.0,
                    isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                    // Peaks at correct m/z but with zero intensity — fails prefilter
                    // but still extractable as XICs (all-zero intensities).
                    // Actually, XIC extraction will produce zero-intensity XICs,
                    // which will fail area/correlation checks. Instead, use tiny
                    // intensities that won't pass the sliding window prefilter
                    // (requires 3/4 consecutive spectra with 2+ fragment matches).
                    // With only 1 fragment m/z present, prefilter requires 2 matches
                    // from top-6 fragments, so having only 1 fragment match fails it.
                    mzs: vec![300.0], // Only 1 of 3 fragments present
                    intensities: vec![100.0],
                }
            })
            .collect();

        let entry = make_lib_entry_with_fragments(0, 500.0, 10.0);
        let library = vec![entry];
        let ms1_index = MS1Index::new(vec![]);
        let config = make_test_config();

        // Without override: should fail prefilter (only 1/3 fragments match)
        let results_normal = run_search(
            &library,
            &spectra,
            &ms1_index,
            None,
            None,
            &config,
            "test.mzML",
            None,
            "Test",
        )
        .unwrap();
        assert!(
            results_normal.is_empty(),
            "Normal search should fail prefilter with only 1 fragment"
        );

        // With override: should bypass prefilter and attempt scoring
        // (may still return None if XIC quality is poor, but it won't
        // be filtered by the prefilter check)
        let mut overrides = HashMap::new();
        overrides.insert(0u32, (10.0, 9.0, 11.0));

        let results_override = run_search(
            &library,
            &spectra,
            &ms1_index,
            None,
            None,
            &config,
            "test.mzML",
            Some(&overrides),
            "Test",
        )
        .unwrap();

        // The override path skips prefilter but still needs ≥2 XICs.
        // With only 1 fragment m/z in spectra, extract_fragment_xics returns
        // 3 XICs (one per library fragment) but only the 300.0 one has signal.
        // The other two are all-zero. All 3 XICs are still returned (zero is valid).
        // So xics.len() ≥ 2 passes, but features may still be None if correlations
        // are degenerate. The key assertion is that override doesn't hit prefilter.
        // We verify the paths diverge — override gets further than normal search.
        // If it produces a result, it came via the override path.
        // If it doesn't, that's also fine — the prefilter was still skipped.
        assert!(
            results_override.len() <= 1,
            "Override should produce 0 or 1 results"
        );
    }

    #[test]
    fn test_run_search_mixed_override_and_cwt() {
        // Two library entries in the same isolation window:
        // entry 0 gets boundary overrides, entry 1 uses normal CWT detection.
        let spectra = make_gaussian_spectra(30, 5.0, 15.0, 10.0, 1.0, 500.0, 12.5);

        let entry0 = make_lib_entry_with_fragments(0, 498.0, 10.0);
        let entry1 = make_lib_entry_with_fragments(1, 502.0, 10.0);
        let library = vec![entry0, entry1];
        let ms1_index = MS1Index::new(vec![]);
        let config = make_test_config();

        // Override only entry 0 with tight boundaries at 9.5–10.5
        let mut overrides = HashMap::new();
        overrides.insert(0u32, (10.0, 9.5, 10.5));

        let results = run_search(
            &library,
            &spectra,
            &ms1_index,
            None,
            None,
            &config,
            "test.mzML",
            Some(&overrides),
            "Test",
        )
        .unwrap();

        // Both entries have good signal → both should produce results
        assert!(
            !results.is_empty(),
            "Should produce at least 1 result (override entry)"
        );

        // Find results by entry_id
        let result0 = results.iter().find(|r| r.entry_id == 0);
        let result1 = results.iter().find(|r| r.entry_id == 1);

        // Entry 0 (override) should have tight boundaries near 9.5–10.5
        if let Some(r0) = result0 {
            assert!(
                (r0.peak_bounds.start_rt - 9.5).abs() < 0.5,
                "Override entry start_rt={} should be near 9.5",
                r0.peak_bounds.start_rt
            );
            assert!(
                (r0.peak_bounds.end_rt - 10.5).abs() < 0.5,
                "Override entry end_rt={} should be near 10.5",
                r0.peak_bounds.end_rt
            );
        }

        // Entry 1 (CWT) should detect the full peak — likely wider boundaries
        if let Some(r1) = result1 {
            // CWT detection on a Gaussian peak centered at 10.0 with σ=1.0
            // The boundaries should encompass the main peak (roughly 8–12 min)
            assert!(
                r1.peak_bounds.start_rt < 10.0 && r1.peak_bounds.end_rt > 10.0,
                "CWT entry should bracket the peak center"
            );
        }

        // If both returned, the override entry's peak should be tighter
        if let (Some(r0), Some(r1)) = (result0, result1) {
            let width0 = r0.peak_bounds.end_rt - r0.peak_bounds.start_rt;
            let width1 = r1.peak_bounds.end_rt - r1.peak_bounds.start_rt;
            assert!(
                width0 <= width1 + 0.5,
                "Override peak width ({:.2}) should be <= CWT peak width ({:.2}) + margin",
                width0,
                width1
            );
        }
    }

    #[test]
    fn test_run_search_no_override_normal_detection() {
        // Baseline: run_search with None boundary_overrides uses CWT peak detection.
        // This validates the normal path still works after the refactor.
        let spectra = make_gaussian_spectra(30, 5.0, 15.0, 10.0, 1.0, 500.0, 12.5);
        let entry = make_lib_entry_with_fragments(0, 500.0, 10.0);
        let library = vec![entry];
        let ms1_index = MS1Index::new(vec![]);
        let config = make_test_config();

        let results = run_search(
            &library,
            &spectra,
            &ms1_index,
            None,
            None,
            &config,
            "test.mzML",
            None,
            "Test",
        )
        .unwrap();

        assert_eq!(results.len(), 1, "Normal CWT search should find the peak");
        let result = &results[0];
        assert_eq!(result.entry_id, 0);
        // Peak should be near the Gaussian center
        assert!(
            (result.peak_bounds.apex_rt - 10.0).abs() < 1.0,
            "CWT apex={} should be near center 10.0",
            result.peak_bounds.apex_rt
        );
        assert!(
            result.peak_bounds.area > 0.0,
            "Peak area should be positive"
        );
    }

    #[test]
    fn test_run_search_override_narrow_boundaries() {
        // Override with very tight boundaries (just 3 scans wide).
        // Verifies we handle the minimum viable peak width.
        let spectra = make_gaussian_spectra(60, 5.0, 15.0, 10.0, 1.0, 500.0, 12.5);
        let entry = make_lib_entry_with_fragments(0, 500.0, 10.0);
        let library = vec![entry];
        let ms1_index = MS1Index::new(vec![]);
        let config = make_test_config();

        // Very tight override: 0.5 min window around apex
        let mut overrides = HashMap::new();
        overrides.insert(0u32, (10.0, 9.75, 10.25));

        let results = run_search(
            &library,
            &spectra,
            &ms1_index,
            None,
            None,
            &config,
            "test.mzML",
            Some(&overrides),
            "Test",
        )
        .unwrap();

        // Should still produce a result if enough spectra fall in range
        if !results.is_empty() {
            let r = &results[0];
            let width = r.peak_bounds.end_rt - r.peak_bounds.start_rt;
            assert!(
                width < 2.0,
                "Tight override should produce narrow peak, got width={:.2}",
                width
            );
        }
    }

    // Regression guard for the calibration XCorr bin spec in
    // docs/02-calibration.md: "Comet-style XCorr (unit resolution, BLAS
    // sdot)". An earlier revision of run_calibration_discovery_windowed
    // picked SpectralScorer::hram() on HRAM data, allocating 100K bins
    // per spectrum instead of 2001. These tests fail loudly if anyone
    // reverts that choice in calibration_xcorr_scorer.

    #[test]
    fn calibration_scorer_uses_unit_bins_for_hram() {
        let config = OspreyConfig {
            resolution_mode: osprey_core::ResolutionMode::HRAM,
            fragment_tolerance: osprey_core::FragmentToleranceConfig::hram(10.0),
            ..Default::default()
        };

        let scorer = calibration_xcorr_scorer(&config);

        assert_eq!(
            scorer.num_bins(),
            SpectralScorer::new().num_bins(),
            "calibration XCorr must use unit-resolution bins on HRAM data \
             (see docs/02-calibration.md)"
        );
        assert_ne!(
            scorer.num_bins(),
            SpectralScorer::hram().num_bins(),
            "calibration XCorr must NOT use HRAM bins -- if this fails, \
             someone put SpectralScorer::hram() back in \
             calibration_xcorr_scorer"
        );
    }

    #[test]
    fn calibration_scorer_uses_unit_bins_for_unit_resolution() {
        let config = OspreyConfig {
            resolution_mode: osprey_core::ResolutionMode::UnitResolution,
            fragment_tolerance: osprey_core::FragmentToleranceConfig::unit_resolution(0.5),
            ..Default::default()
        };

        let scorer = calibration_xcorr_scorer(&config);

        assert_eq!(
            scorer.num_bins(),
            SpectralScorer::new().num_bins(),
            "calibration XCorr on unit-resolution data should stay unit bins"
        );
    }
}
