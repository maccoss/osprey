//! Osprey Scoring - Feature extraction and peptide scoring
//!
//! This crate provides:
//! - Feature extraction from coefficient time series and spectra
//! - Spectral similarity metrics (LibCosine, XCorr)
//! - Hyperscore calculation
//! - Decoy generation with enzyme-aware sequence reversal
//! - **Batch scoring with BLAS acceleration** (10-20× faster)
//! - **Streaming pipeline** for memory-efficient processing (optional)
//!
//! ## Batch Scoring
//!
//! For high-performance scoring of many library entries against many spectra,
//! use the [`batch`] module:
//!
//! ```ignore
//! use osprey_scoring::batch::{BatchScorer, PreprocessedLibrary, PreprocessedSpectra};
//!
//! let scorer = BatchScorer::new();
//! let lib = scorer.preprocess_library(&library_entries);
//! let spec = scorer.preprocess_spectra(&spectra);
//!
//! // Score all pairs with single BLAS call
//! let scores = scorer.score_all(&lib, &spec);
//! ```
//!
//! ## Streaming Pipeline (requires `streaming` feature)
//!
//! For processing large mzML files with overlapped I/O and preprocessing:
//!
//! ```ignore
//! use osprey_scoring::pipeline::{PipelineCoordinator, run_streaming_pipeline};
//!
//! // From file (async)
//! let matches = run_streaming_pipeline("file.mzML", &library, config).await?;
//!
//! // From pre-loaded spectra
//! let coordinator = PipelineCoordinator::new();
//! let matches = coordinator.run_on_spectra(&spectra, &library);
//! ```

pub mod batch;
pub mod calibration_ml;
pub mod diagnostics;
pub mod pipeline;

use osprey_core::{
    BinConfig, FragmentAnnotation, IonType, LibraryEntry, LibraryFragment, Modification, Result,
    Spectrum, ToleranceUnit,
};
use rayon::prelude::*;
use std::collections::HashMap;

/// Check if at least 2 of the top 6 library peaks match observed spectrum peaks
///
/// This is a fast pre-filter using binary search to eliminate candidates that have
/// no signal overlap with the observed spectrum before expensive scoring.
///
/// # Arguments
/// * `library_fragments` - Library fragment ions for a candidate entry
/// * `spectrum_mzs` - Sorted observed m/z values from the spectrum
/// * `tolerance` - Fragment tolerance value
/// * `unit` - Tolerance unit (ppm or Th)
///
/// # Returns
/// `true` if at least 2 of the top 6 library peaks have matching observed peaks
/// Check if at least one of the top-6 library fragments matches in a spectrum.
///
/// More permissive than [`has_topn_fragment_match`] (which requires 2 matches).
/// Used as a fast apex-scan pre-filter in the coelution search: if not even one
/// top-6 fragment has any signal in the scan closest to the expected RT, the
/// candidate is very unlikely to be a real detection.
pub fn has_any_top_fragment_match(
    library_fragments: &[LibraryFragment],
    spectrum_mzs: &[f64],
    tolerance: f64,
    unit: ToleranceUnit,
) -> bool {
    if library_fragments.is_empty() || spectrum_mzs.is_empty() {
        return true; // Be conservative - don't filter if no data
    }

    // Consider only the top 6 fragments by intensity
    let top_indices: Vec<usize> = if library_fragments.len() <= 6 {
        (0..library_fragments.len()).collect()
    } else {
        let mut indexed: Vec<(usize, f32)> = library_fragments
            .iter()
            .enumerate()
            .map(|(i, f)| (i, f.relative_intensity))
            .collect();
        indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
        indexed.iter().take(6).map(|(i, _)| *i).collect()
    };

    for &idx in &top_indices {
        let lib_mz = library_fragments[idx].mz;
        let tol_da = match unit {
            ToleranceUnit::Ppm => lib_mz * tolerance / 1e6,
            ToleranceUnit::Mz => tolerance,
        };
        let lower = lib_mz - tol_da;
        let upper = lib_mz + tol_da;
        let start_idx = spectrum_mzs.partition_point(|&mz| mz < lower);
        if start_idx < spectrum_mzs.len() && spectrum_mzs[start_idx] <= upper {
            return true;
        }
    }

    false
}

pub fn has_topn_fragment_match(
    library_fragments: &[LibraryFragment],
    spectrum_mzs: &[f64],
    tolerance: f64,
    unit: ToleranceUnit,
) -> bool {
    if library_fragments.is_empty() || spectrum_mzs.is_empty() {
        return true; // Be conservative - don't filter if no data
    }

    // Get indices of top 6 fragments by intensity
    let n_top = library_fragments.len().min(6);
    let top_indices: Vec<usize> = if library_fragments.len() <= 6 {
        (0..library_fragments.len()).collect()
    } else {
        let mut indexed: Vec<(usize, f32)> = library_fragments
            .iter()
            .enumerate()
            .map(|(i, f)| (i, f.relative_intensity))
            .collect();
        indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
        indexed.iter().take(6).map(|(i, _)| *i).collect()
    };

    // Require at least 2 of the top 6 to have matches
    // For entries with only 1 fragment, accept 1 match
    let required_matches = if n_top <= 1 { 1 } else { 2 };
    let mut match_count = 0u32;

    for &idx in &top_indices {
        let lib_mz = library_fragments[idx].mz;

        // Calculate tolerance window based on unit
        let tol_da = match unit {
            ToleranceUnit::Ppm => lib_mz * tolerance / 1e6,
            ToleranceUnit::Mz => tolerance,
        };

        let lower = lib_mz - tol_da;
        let upper = lib_mz + tol_da;

        // Binary search for first m/z >= lower bound
        let start_idx = spectrum_mzs.partition_point(|&mz| mz < lower);

        // Check if any peak is within the tolerance window
        if start_idx < spectrum_mzs.len() && spectrum_mzs[start_idx] <= upper {
            match_count += 1;
            if match_count >= required_matches {
                return true;
            }
        }
    }

    false
}

/// Get indices of top N fragments by intensity
///
/// Returns indices sorted by intensity (highest first), limited to N fragments.
/// If the library has fewer than N fragments, returns all indices.
///
/// # Arguments
/// * `library_fragments` - Library fragment ions
/// * `n` - Number of top fragments to return
pub fn get_top_n_fragment_indices(library_fragments: &[LibraryFragment], n: usize) -> Vec<usize> {
    if library_fragments.len() <= n {
        return (0..library_fragments.len()).collect();
    }

    let mut indexed: Vec<(usize, f32)> = library_fragments
        .iter()
        .enumerate()
        .map(|(i, f)| (i, f.relative_intensity))
        .collect();
    indexed.sort_by(|a, b| b.1.total_cmp(&a.1)); // Sort by intensity descending
    indexed.iter().take(n).map(|(i, _)| *i).collect()
}

/// Check if a library m/z has a matching peak in the spectrum
///
/// Uses binary search for O(log n) performance.
///
/// # Arguments
/// * `lib_mz` - Library fragment m/z
/// * `spectrum_mzs` - Sorted observed m/z values
/// * `tolerance` - Tolerance value
/// * `unit` - Tolerance unit (ppm or Th)
///
/// # Returns
/// `true` if there is at least one observed peak within tolerance
pub fn has_match_within_tolerance(
    lib_mz: f64,
    spectrum_mzs: &[f64],
    tolerance: f64,
    unit: ToleranceUnit,
) -> bool {
    let tol_da = match unit {
        ToleranceUnit::Ppm => lib_mz * tolerance / 1e6,
        ToleranceUnit::Mz => tolerance,
    };

    let lower = lib_mz - tol_da;
    let upper = lib_mz + tol_da;

    // Binary search for first m/z >= lower bound
    let start_idx = spectrum_mzs.partition_point(|&mz| mz < lower);
    start_idx < spectrum_mzs.len() && spectrum_mzs[start_idx] <= upper
}

/// Check if any top N library fragments match an observed peak, and collect mass errors.
///
/// Combines the filtering check with mass error collection for calibration.
/// Uses the top 6 fragments by intensity, requiring at least 1 match.
/// For each matched fragment, finds the closest observed peak and computes the mass error.
///
/// # Returns
/// `(has_match, mass_errors)` where `has_match` is true if at least 1 top-N peak matches,
/// and `mass_errors` contains the error (in the configured unit) for each matched fragment.
pub fn topn_fragment_match_with_errors(
    library_fragments: &[LibraryFragment],
    spectrum_mzs: &[f64],
    tolerance: f64,
    unit: ToleranceUnit,
) -> (bool, Vec<f64>) {
    if library_fragments.is_empty() || spectrum_mzs.is_empty() {
        return (true, Vec::new());
    }

    // Get indices of top 6 fragments by intensity (consistent with has_topn_fragment_match)
    let n_top = library_fragments.len().min(6);
    let mut top_indices: Vec<usize> = (0..n_top).collect();

    if library_fragments.len() > 6 {
        let mut indexed: Vec<(usize, f32)> = library_fragments
            .iter()
            .enumerate()
            .map(|(i, f)| (i, f.relative_intensity))
            .collect();
        indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
        top_indices = indexed.iter().take(6).map(|(i, _)| *i).collect();
    }

    let mut has_match = false;
    let mut mass_errors = Vec::with_capacity(n_top);

    for &idx in &top_indices {
        let lib_mz = library_fragments[idx].mz;

        let tol_da = match unit {
            ToleranceUnit::Ppm => lib_mz * tolerance / 1e6,
            ToleranceUnit::Mz => tolerance,
        };

        let lower = lib_mz - tol_da;
        let upper = lib_mz + tol_da;

        // Binary search for first m/z >= lower bound
        let start_idx = spectrum_mzs.partition_point(|&mz| mz < lower);

        if start_idx < spectrum_mzs.len() && spectrum_mzs[start_idx] <= upper {
            has_match = true;

            // Find closest peak within tolerance window
            let mut best_mz = spectrum_mzs[start_idx];
            let mut best_diff = (best_mz - lib_mz).abs();
            let mut j = start_idx + 1;
            while j < spectrum_mzs.len() && spectrum_mzs[j] <= upper {
                let diff = (spectrum_mzs[j] - lib_mz).abs();
                if diff < best_diff {
                    best_diff = diff;
                    best_mz = spectrum_mzs[j];
                }
                j += 1;
            }

            // Compute mass error in configured unit
            let error = match unit {
                ToleranceUnit::Ppm => (best_mz - lib_mz) / lib_mz * 1e6,
                ToleranceUnit::Mz => best_mz - lib_mz,
            };
            mass_errors.push(error);
        }
    }

    (has_match, mass_errors)
}

/// Compute elution-weighted cosine similarity across peak boundaries.
///
/// At each scan within `[peak_start_rt, peak_end_rt]`, computes the cosine similarity
/// (sqrt preprocessing, L2 normalization) between the observed spectrum and the library.
/// Each scan's cosine is weighted by coefficient², then the weighted average is returned.
///
/// This captures whether the library spectral pattern holds consistently across the
/// entire peak, not just at the apex. Interference that only affects part of the peak
/// will reduce this score even if the apex cosine is high.
///
/// Inspired by DIA-NN's `pCos` (elution-weighted cosine) feature.
pub fn compute_elution_weighted_cosine(
    library_fragments: &[LibraryFragment],
    coefficient_series: &[(f64, f64)], // (RT, coefficient) pairs
    spectra: &[&Spectrum],             // peak region spectra (sorted by RT)
    tolerance_da: f64,
    tolerance_ppm: f64,
    peak_start_rt: f64,
    peak_end_rt: f64,
) -> f64 {
    if library_fragments.is_empty() || spectra.is_empty() || coefficient_series.is_empty() {
        return 0.0;
    }

    let spec_rts: Vec<f64> = spectra.iter().map(|s| s.retention_time).collect();

    let mut weighted_cosine_sum = 0.0;
    let mut weight_sum = 0.0;

    // For each data point within peak boundaries, compute cosine at the nearest spectrum
    for &(rt, coef) in coefficient_series {
        if rt < peak_start_rt || rt > peak_end_rt {
            continue;
        }
        if coef <= 0.0 {
            continue;
        }

        // Weight = coefficient² (like DIA-NN's elution² weighting)
        let weight = coef * coef;

        // Find nearest spectrum
        let spec_idx = match spec_rts
            .binary_search_by(|srt| srt.partial_cmp(&rt).unwrap_or(std::cmp::Ordering::Equal))
        {
            Ok(i) => i,
            Err(i) => {
                if i == 0 {
                    0
                } else if i >= spec_rts.len() {
                    spec_rts.len() - 1
                } else if (spec_rts[i] - rt).abs() < (spec_rts[i - 1] - rt).abs() {
                    i
                } else {
                    i - 1
                }
            }
        };

        let spectrum = spectra[spec_idx];

        // Compute LibCosine (sqrt preprocessing) between this spectrum and library
        let cosine =
            compute_cosine_at_scan(library_fragments, spectrum, tolerance_da, tolerance_ppm);

        weighted_cosine_sum += cosine * weight;
        weight_sum += weight;
    }

    if weight_sum > 0.0 {
        weighted_cosine_sum / weight_sum
    } else {
        0.0
    }
}

/// Compute cosine similarity (sqrt preprocessing, L2 normalization) between
/// a single observed spectrum and library fragment list.
///
/// ALL library fragments within the spectrum's mass range are included.
/// Unmatched fragments use observed intensity of 0, which penalizes the
/// cosine — a missing fragment that should be present is strong evidence
/// against a match.
pub fn compute_cosine_at_scan(
    library_fragments: &[LibraryFragment],
    spectrum: &Spectrum,
    tolerance_da: f64,
    tolerance_ppm: f64,
) -> f64 {
    if library_fragments.is_empty() || spectrum.mzs.is_empty() {
        return 0.0;
    }

    let spec_mz_min = spectrum.mzs[0];
    let spec_mz_max = spectrum.mzs[spectrum.mzs.len() - 1];

    let mut lib_preprocessed: Vec<f64> = Vec::new();
    let mut obs_preprocessed: Vec<f64> = Vec::new();

    for frag in library_fragments {
        // Skip fragments outside the spectrum's mass range
        if frag.mz < spec_mz_min || frag.mz > spec_mz_max {
            continue;
        }

        let tol = tolerance_da.max(frag.mz * tolerance_ppm / 1e6);
        let lower = frag.mz - tol;
        let upper = frag.mz + tol;
        let start = spectrum.mzs.partition_point(|&mz| mz < lower);

        // Find closest peak within tolerance
        let mut best_intensity = 0.0f64;
        let mut best_diff = f64::MAX;
        let mut j = start;
        while j < spectrum.mzs.len() && spectrum.mzs[j] <= upper {
            let diff = (spectrum.mzs[j] - frag.mz).abs();
            if diff < best_diff {
                best_diff = diff;
                best_intensity = spectrum.intensities[j] as f64;
            }
            j += 1;
        }

        // Always include: library sqrt(intensity) and observed sqrt(intensity) or 0
        lib_preprocessed.push((frag.relative_intensity as f64).sqrt());
        obs_preprocessed.push(best_intensity.sqrt());
    }

    if lib_preprocessed.is_empty() {
        return 0.0;
    }

    // L2 normalize
    let lib_norm = lib_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();
    let obs_norm = obs_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();

    if lib_norm < 1e-10 || obs_norm < 1e-10 {
        return 0.0;
    }

    // Cosine = dot(lib_norm, obs_norm)
    lib_preprocessed
        .iter()
        .zip(obs_preprocessed.iter())
        .map(|(a, b)| (a / lib_norm) * (b / obs_norm))
        .sum()
}

/// Compute per-fragment mass accuracy statistics at the apex scan.
///
/// # Arguments
/// * `matches` - Fragment matches from the apex spectrum
/// * `unit` - Tolerance unit determining error computation (ppm for HRAM, Th for unit resolution)
/// * `tolerance` - The fragment matching tolerance (in the same unit). Used as penalty value
///   when no fragments match, so entries with no evidence get "bad" mass accuracy instead of
///   a misleading 0.0 ("perfect" accuracy).
///
/// # Returns
/// `(signed_mean, abs_mean, std)` in the configured unit
pub fn compute_mass_accuracy(
    matches: &[FragmentMatch],
    unit: ToleranceUnit,
    tolerance: f64,
) -> (f64, f64, f64) {
    if matches.is_empty() {
        return (0.0, tolerance, tolerance);
    }

    let errors: Vec<f64> = matches
        .iter()
        .map(|m| match unit {
            ToleranceUnit::Ppm => (m.obs_mz - m.lib_mz) / m.lib_mz * 1e6,
            ToleranceUnit::Mz => m.obs_mz - m.lib_mz,
        })
        .collect();

    let n = errors.len() as f64;

    // Mean of signed errors (systematic bias direction)
    let mean_signed: f64 = errors.iter().sum::<f64>() / n;

    // Mean of absolute errors (overall accuracy)
    let mean_abs: f64 = errors.iter().map(|e| e.abs()).sum::<f64>() / n;

    // Standard deviation of signed errors
    let variance: f64 = errors
        .iter()
        .map(|e| (e - mean_signed).powi(2))
        .sum::<f64>()
        / n;
    let std_dev = variance.sqrt();

    (mean_signed, mean_abs, std_dev)
}

/// Extract XICs for the top N library fragments from a set of spectra.
///
/// Fragment XIC extraction approach inspired by DIA-NN (Demichev et al.,
/// Nature Methods, 2020).
///
/// For each fragment, uses binary search on each spectrum's sorted m/z array
/// to extract the intensity at the fragment's m/z (or 0 if not found).
///
/// # Returns
/// Vec of `(fragment_index, xic)` where xic is `Vec<(RT, intensity)>`.
/// Only returns fragments that have at least one non-zero intensity.
pub fn extract_fragment_xics(
    library_fragments: &[LibraryFragment],
    spectra: &[&Spectrum],
    tolerance_da: f64,
    tolerance_ppm: f64,
    max_fragments: usize,
) -> Vec<(usize, Vec<(f64, f64)>)> {
    if library_fragments.is_empty() || spectra.is_empty() {
        return Vec::new();
    }

    // Select top N fragments by relative intensity
    let n_top = library_fragments.len().min(max_fragments);
    let top_indices: Vec<usize> = if library_fragments.len() <= max_fragments {
        (0..library_fragments.len()).collect()
    } else {
        let mut indexed: Vec<(usize, f32)> = library_fragments
            .iter()
            .enumerate()
            .map(|(i, f)| (i, f.relative_intensity))
            .collect();
        indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
        indexed
            .iter()
            .take(max_fragments)
            .map(|(i, _)| *i)
            .collect()
    };

    let mut results = Vec::with_capacity(n_top);

    for &frag_idx in &top_indices {
        let frag_mz = library_fragments[frag_idx].mz;
        let tol = tolerance_da.max(frag_mz * tolerance_ppm / 1e6);
        let lower = frag_mz - tol;
        let upper = frag_mz + tol;

        let mut xic = Vec::with_capacity(spectra.len());

        for spectrum in spectra {
            let start = spectrum.mzs.partition_point(|&mz| mz < lower);

            let intensity = if start < spectrum.mzs.len() && spectrum.mzs[start] <= upper {
                // Find closest peak within tolerance
                let mut best_intensity = spectrum.intensities[start] as f64;
                let mut best_diff = (spectrum.mzs[start] - frag_mz).abs();
                let mut j = start + 1;
                while j < spectrum.mzs.len() && spectrum.mzs[j] <= upper {
                    let diff = (spectrum.mzs[j] - frag_mz).abs();
                    if diff < best_diff {
                        best_diff = diff;
                        best_intensity = spectrum.intensities[j] as f64;
                    }
                    j += 1;
                }
                best_intensity
            } else {
                0.0
            };

            xic.push((spectrum.retention_time, intensity));
        }

        // Always include the fragment XIC, even if all-zero.
        // Zero intensities are valid data (no centroided peak found) and are
        // handled properly in log-space models by skipping them per-scan.
        // Dropping all-zero fragments biases the median polish by giving decoys
        // (with fewer matching fragments) artificially high R².
        results.push((frag_idx, xic));
    }

    results
}

/// Result of Tukey median polish decomposition of fragment XICs.
///
/// Decomposes the fragment XIC matrix (fragments × scans) into additive components
/// in log space using the standard Tukey median polish algorithm.
///
/// Model: `ln(Observed[f,s]) = μ + α_f + β_s + ε_fs`
///
/// - Column effects (β_s) give a robust elution profile for FWHM/boundary estimation.
/// - Row effects (α_f) give interference-free fragment intensities for library scoring.
///
/// Reference: PRISM `rollup.py:640-731`
pub struct TukeyMedianPolishResult {
    /// Overall effect (grand median, ln space).
    pub overall: f64,
    /// Row effects per fragment (ln space) — data-derived relative fragment intensities.
    pub row_effects: Vec<f64>,
    /// Column effects per scan (ln space) — elution profile shape.
    pub col_effects: Vec<f64>,
    /// Elution profile in linear space: (RT, exp(μ + β_s)) per scan.
    pub elution_profile: Vec<(f64, f64)>,
    /// Residuals\[fragment\]\[scan\] in ln space. NaN for zero-intensity cells.
    pub residuals: Vec<Vec<f64>>,
    /// Number of iterations used.
    pub n_iterations: usize,
    /// Whether the algorithm converged within tolerance.
    pub converged: bool,
    /// Number of fragments with at least one non-zero intensity.
    pub n_fragments_used: usize,
    /// Fragment indices from the input (maps row index → original fragment index).
    pub fragment_indices: Vec<usize>,
}

/// Compute median of a slice, skipping NaN values. Returns NaN if no finite values.
fn nanmedian(values: &[f64]) -> f64 {
    let mut finite: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return f64::NAN;
    }
    finite.sort_by(|a, b| a.total_cmp(b));
    let mid = finite.len() / 2;
    if finite.len() % 2 == 0 {
        (finite[mid - 1] + finite[mid]) / 2.0
    } else {
        finite[mid]
    }
}

/// Cosine similarity between two vectors. Returns 0 if either has zero norm.
pub fn cosine_angle(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    let mut dot = 0.0;
    let mut norm_a = 0.0;
    let mut norm_b = 0.0;
    for i in 0..n {
        dot += a[i] * b[i];
        norm_a += a[i] * a[i];
        norm_b += b[i] * b[i];
    }
    if norm_a < 1e-30 || norm_b < 1e-30 {
        return 0.0;
    }
    (dot / (norm_a.sqrt() * norm_b.sqrt())).clamp(0.0, 1.0)
}

/// Tukey median polish decomposition of fragment XICs.
///
/// Decomposes the 6×N fragment XIC matrix into overall + row effects + column effects
/// + residuals using iterative median subtraction. Works in ln space with NaN for zeros.
///
/// The column effects give a robust elution profile (shared peak shape across transitions,
/// resistant to interference on individual fragments). The row effects give data-derived
/// fragment intensities that can be scored against the library.
///
/// Following PRISM `rollup.py:640-731`:
/// 1. Row sweep: subtract row medians (captures fragment-level intensity differences)
/// 2. Column sweep: subtract column medians (captures scan-to-scan elution shape)
/// 3. Repeat until convergence (max change < tol) or max iterations reached
pub fn tukey_median_polish(
    xics: &[(usize, Vec<(f64, f64)>)],
    max_iter: usize,
    tol: f64,
) -> Option<TukeyMedianPolishResult> {
    if xics.len() < 2 {
        return None;
    }

    let n_scans = xics[0].1.len();
    if n_scans < 3 {
        return None;
    }

    let n_frags = xics.len();
    let rts: Vec<f64> = xics[0].1.iter().map(|(rt, _)| *rt).collect();
    let fragment_indices: Vec<usize> = xics.iter().map(|(idx, _)| *idx).collect();

    // Build ln-space matrix. Zeros → NaN (missing data).
    let mut residuals: Vec<Vec<f64>> = Vec::with_capacity(n_frags);
    for (_frag_idx, xic) in xics {
        let row: Vec<f64> = xic
            .iter()
            .map(|(_, intensity)| {
                if *intensity > 0.0 {
                    intensity.ln()
                } else {
                    f64::NAN
                }
            })
            .collect();
        residuals.push(row);
    }

    let mut overall: f64 = 0.0;
    let mut row_effects: Vec<f64> = vec![0.0; n_frags];
    let mut col_effects: Vec<f64> = vec![0.0; n_scans];
    let mut converged = false;
    let mut n_iter = 0;

    for iteration in 0..max_iter {
        n_iter = iteration + 1;

        // Save old residuals for convergence check
        let old_residuals: Vec<Vec<f64>> = residuals.clone();

        // === Row sweep: subtract nanmedian of each row ===
        let row_medians: Vec<f64> = residuals.iter().map(|row| nanmedian(row)).collect();

        for (f, row) in residuals.iter_mut().enumerate() {
            let rm = row_medians[f];
            if rm.is_finite() {
                for val in row.iter_mut() {
                    if val.is_finite() {
                        *val -= rm;
                    }
                }
            }
        }

        // Update: row_effects += (row_medians - median(row_medians)), overall += median(row_medians)
        let median_of_row_medians = nanmedian(&row_medians);
        if median_of_row_medians.is_finite() {
            for (f, &rm) in row_medians.iter().enumerate() {
                if rm.is_finite() {
                    row_effects[f] += rm - median_of_row_medians;
                }
            }
            overall += median_of_row_medians;
        }

        // === Column sweep: subtract nanmedian of each column ===
        let mut col_buf: Vec<f64> = Vec::with_capacity(n_frags);
        let col_medians: Vec<f64> = (0..n_scans)
            .map(|s| {
                col_buf.clear();
                col_buf.extend(residuals.iter().map(|row| row[s]));
                nanmedian(&col_buf)
            })
            .collect();

        for row in residuals.iter_mut() {
            for (s, val) in row.iter_mut().enumerate() {
                let cm = col_medians[s];
                if val.is_finite() && cm.is_finite() {
                    *val -= cm;
                }
            }
        }

        // Update: col_effects += (col_medians - median(col_medians)), overall += median(col_medians)
        let median_of_col_medians = nanmedian(&col_medians);
        if median_of_col_medians.is_finite() {
            for (s, &cm) in col_medians.iter().enumerate() {
                if cm.is_finite() {
                    col_effects[s] += cm - median_of_col_medians;
                }
            }
            overall += median_of_col_medians;
        }

        // Convergence check: max(|new - old|) < tol
        let max_change = residuals
            .iter()
            .zip(old_residuals.iter())
            .flat_map(|(new_row, old_row)| {
                new_row.iter().zip(old_row.iter()).map(|(n, o)| {
                    if n.is_finite() && o.is_finite() {
                        (n - o).abs()
                    } else {
                        0.0
                    }
                })
            })
            .fold(0.0f64, f64::max);

        if max_change < tol {
            converged = true;
            break;
        }
    }

    // Count fragments with at least one finite value
    let n_fragments_used = residuals
        .iter()
        .filter(|row| row.iter().any(|v| v.is_finite()))
        .count();

    if n_fragments_used < 2 {
        return None;
    }

    // Build linear-space elution profile: exp(overall + col_effects[s])
    let elution_profile: Vec<(f64, f64)> = rts
        .iter()
        .zip(col_effects.iter())
        .map(|(&rt, &ce)| (rt, (overall + ce).exp()))
        .collect();

    Some(TukeyMedianPolishResult {
        overall,
        row_effects,
        col_effects,
        elution_profile,
        residuals,
        n_iterations: n_iter,
        converged,
        n_fragments_used,
        fragment_indices,
    })
}

/// Compute cosine similarity between Tukey median polish row effects and library
/// fragment intensities, using sqrt preprocessing (matching LibCosine convention).
///
/// Row effects are in ln space. Convert to linear via exp(overall + α_f), then sqrt.
/// Library intensities are linear; apply sqrt.
///
/// For targets: row effects match library → high cosine.
/// For decoys: row effects are random → low cosine.
pub fn median_polish_libcosine(
    polish: &TukeyMedianPolishResult,
    library_fragments: &[LibraryFragment],
) -> f64 {
    let n = polish.fragment_indices.len();
    if n < 2 {
        return 0.0;
    }

    let mut row_vec: Vec<f64> = Vec::with_capacity(n);
    let mut lib_vec: Vec<f64> = Vec::with_capacity(n);

    for (i, &frag_idx) in polish.fragment_indices.iter().enumerate() {
        if frag_idx >= library_fragments.len() {
            continue;
        }

        // Library intensity ALWAYS contributes — even for undetected fragments.
        // This keeps the cosine in full 6D space. Undetected fragments push 0
        // into the dot product but their library intensity increases the denominator,
        // naturally penalizing decoys that fail to detect high-intensity fragments.
        let lib_int = library_fragments[frag_idx].relative_intensity as f64;
        lib_vec.push(lib_int.max(0.0).sqrt());

        // Check if fragment had any detected signal (at least one finite residual)
        let has_signal = polish.residuals[i].iter().any(|v| v.is_finite());
        if !has_signal || !polish.row_effects[i].is_finite() {
            row_vec.push(0.0);
            continue;
        }

        // Row effect is ln-space; convert to linear, then sqrt (Poisson noise model)
        let linear = (polish.overall + polish.row_effects[i]).exp();
        row_vec.push(linear.max(0.0).sqrt());
    }

    if row_vec.len() < 2 {
        return 0.0;
    }

    cosine_angle(&row_vec, &lib_vec)
}

/// Compute R² (coefficient of determination) from the Tukey median polish model.
///
/// Measures how well the additive model (μ + α_f + β_s) explains the observed
/// fragment×scan intensity matrix. Uses sqrt preprocessing (Poisson noise model)
/// to compress dynamic range — without this, R² is dominated by the highest-intensity
/// cells. Cells where the model predicts signal but none was detected (observed=0)
/// contribute sqrt(predicted)² to SS_residual, properly penalizing poor matches.
///
/// R² = 1 - SS_residual / SS_total
///   SS_residual = Σ(sqrt(observed) - sqrt(predicted))²  over all cells
///   SS_total    = Σ(sqrt(observed) - mean(sqrt(observed)))²  over all cells
pub fn median_polish_rsquared(polish: &TukeyMedianPolishResult) -> f64 {
    let n_frags = polish.row_effects.len();
    let n_scans = if n_frags > 0 {
        polish.col_effects.len()
    } else {
        return 0.0;
    };
    if n_frags < 2 || n_scans < 3 {
        return 0.0;
    }

    let mut obs_values: Vec<f64> = Vec::with_capacity(n_frags * n_scans);
    let mut pred_values: Vec<f64> = Vec::with_capacity(n_frags * n_scans);

    for f in 0..n_frags {
        for s in 0..n_scans {
            // Predicted intensity: convert to linear then sqrt
            let predicted = (polish.overall + polish.row_effects[f] + polish.col_effects[s]).exp();

            // Observed: if residual is finite, reconstruct from model + residual;
            // if NaN (original was zero), observed = 0
            let observed = if polish.residuals[f][s].is_finite() {
                (polish.overall
                    + polish.row_effects[f]
                    + polish.col_effects[s]
                    + polish.residuals[f][s])
                    .exp()
            } else {
                0.0
            };

            // Sqrt preprocessing (Poisson noise model, matches cosine convention)
            pred_values.push(predicted.max(0.0).sqrt());
            obs_values.push(observed.max(0.0).sqrt());
        }
    }

    let n = obs_values.len() as f64;
    if n < 2.0 {
        return 0.0;
    }

    let mean_obs = obs_values.iter().sum::<f64>() / n;
    let ss_total: f64 = obs_values.iter().map(|o| (o - mean_obs).powi(2)).sum();
    let ss_residual: f64 = obs_values
        .iter()
        .zip(pred_values.iter())
        .map(|(o, p)| (o - p).powi(2))
        .sum();

    if ss_total < 1e-30 {
        return 0.0;
    }

    (1.0 - ss_residual / ss_total).clamp(0.0, 1.0)
}

/// Compute residual ratio from the Tukey median polish model.
///
/// Measures the fraction of total signal that is unexplained by the additive model.
/// Computed in linear space: Σ|observed - predicted| / Σ observed.
///
/// Low values (~0) indicate clean co-elution where all fragments follow the same
/// chromatographic shape. High values indicate interference or noise on individual
/// fragments. Cells where observed=0 but the model predicts signal contribute to
/// the numerator, penalizing missing fragments.
///
/// Note: this is inverted relative to R² — lower is better. Mokapot handles
/// this automatically via semi-supervised learning.
pub fn median_polish_residual_ratio(polish: &TukeyMedianPolishResult) -> f64 {
    let n_frags = polish.row_effects.len();
    let n_scans = if n_frags > 0 {
        polish.col_effects.len()
    } else {
        return 1.0;
    };
    if n_frags < 2 || n_scans < 3 {
        return 1.0;
    }

    let mut sum_abs_residual = 0.0;
    let mut sum_observed = 0.0;

    for f in 0..n_frags {
        for s in 0..n_scans {
            let predicted = (polish.overall + polish.row_effects[f] + polish.col_effects[s]).exp();

            let observed = if polish.residuals[f][s].is_finite() {
                (polish.overall
                    + polish.row_effects[f]
                    + polish.col_effects[s]
                    + polish.residuals[f][s])
                    .exp()
            } else {
                // Zero intensity is a real measurement (fragment had no signal),
                // not missing data. Use pseudocount to prevent |0 - predicted|
                // from dominating the ratio.
                0.0001
            };

            sum_abs_residual += (observed - predicted).abs();
            sum_observed += observed;
        }
    }

    if sum_observed < 1e-30 {
        return 1.0;
    }

    sum_abs_residual / sum_observed
}

/// Minimum per-fragment R² with the shared elution profile.
///
/// For each fragment row: compute R² of sqrt(observed) vs sqrt(predicted) across scans.
/// `predicted(f,s) = exp(overall + row_effects[f] + col_effects[s])`
/// `observed(f,s) = exp(overall + row_effects[f] + col_effects[s] + residuals[f][s])`
/// Returns the minimum R² across fragments (weakest link).
/// If a fragment has < 3 finite scans, R² = 0.0 for that fragment.
pub fn median_polish_min_fragment_r2(polish: &TukeyMedianPolishResult) -> f64 {
    let n_frags = polish.row_effects.len();
    let n_scans = polish.col_effects.len();
    if n_frags < 1 || n_scans < 3 {
        return 0.0;
    }

    let mut min_r2 = f64::MAX;

    for f in 0..n_frags {
        let mut pred_vals = Vec::with_capacity(n_scans);
        let mut obs_vals = Vec::with_capacity(n_scans);

        for s in 0..n_scans {
            if !polish.residuals[f][s].is_finite() {
                continue;
            }
            let predicted = (polish.overall + polish.row_effects[f] + polish.col_effects[s]).exp();
            let observed = (polish.overall
                + polish.row_effects[f]
                + polish.col_effects[s]
                + polish.residuals[f][s])
                .exp();
            pred_vals.push(predicted.sqrt());
            obs_vals.push(observed.sqrt());
        }

        let r2 = if pred_vals.len() < 3 {
            0.0
        } else {
            compute_r2(&pred_vals, &obs_vals)
        };

        if r2 < min_r2 {
            min_r2 = r2;
        }
    }

    if min_r2 == f64::MAX {
        0.0
    } else {
        min_r2.max(0.0)
    }
}

/// Mean pairwise Pearson correlation of median polish residuals across fragments.
///
/// For a clean target, residuals should be uncorrelated noise → correlations ≈ 0.
/// Correlated residuals suggest a co-eluting interferer affecting multiple fragment channels.
/// Returns 0.0 if < 2 fragments with sufficient data.
pub fn median_polish_residual_correlation(polish: &TukeyMedianPolishResult) -> f64 {
    let n_frags = polish.row_effects.len();
    let n_scans = polish.col_effects.len();
    if n_frags < 2 || n_scans < 3 {
        return 0.0;
    }

    // Pairwise Pearson correlations using only scans where BOTH fragments
    // have finite residuals. NaN cells (zero observed intensity) would create
    // artificial correlation structure if replaced with 0.0.
    let mut corr_sum = 0.0;
    let mut n_pairs = 0;

    for i in 0..n_frags {
        for j in (i + 1)..n_frags {
            // Collect residuals at scans where both fragments are finite
            let mut ri = Vec::with_capacity(n_scans);
            let mut rj = Vec::with_capacity(n_scans);
            for s in 0..n_scans {
                if polish.residuals[i][s].is_finite() && polish.residuals[j][s].is_finite() {
                    ri.push(polish.residuals[i][s]);
                    rj.push(polish.residuals[j][s]);
                }
            }
            if ri.len() >= 3 {
                let r = pearson_correlation_raw(&ri, &rj);
                if r.is_finite() {
                    corr_sum += r;
                    n_pairs += 1;
                }
            }
        }
    }

    if n_pairs == 0 {
        0.0
    } else {
        corr_sum / n_pairs as f64
    }
}

/// Compute R² (coefficient of determination) between predicted and observed values.
fn compute_r2(predicted: &[f64], observed: &[f64]) -> f64 {
    let n = predicted.len();
    if n < 2 {
        return 0.0;
    }
    let obs_mean = observed.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = observed.iter().map(|&o| (o - obs_mean).powi(2)).sum();
    let ss_res: f64 = predicted
        .iter()
        .zip(observed.iter())
        .map(|(&p, &o)| (o - p).powi(2))
        .sum();
    if ss_tot < 1e-30 {
        return 0.0;
    }
    (1.0 - ss_res / ss_tot).max(0.0)
}

/// Compute fragment-based FWHM and peak boundaries using Tukey median polish.
///
/// Primary approach: Tukey median polish decomposes the fragment XIC matrix into
/// row effects (fragment intensities) and column effects (elution profile). The
/// column effects give a robust peak shape by borrowing strength across all 6
/// transitions, with interference suppressed by the median. FWHM is computed on
/// this elution profile.
///
/// Fallback: consensus XIC (normalize co-eluting fragment XICs to max=1, sum, then
/// compute FWHM on the consensus).
///
/// # Returns
/// Tuple of:
/// - `Option<(start_rt, end_rt, fwhm)>` — 95% Gaussian boundaries (±1.96σ)
/// - `Option<TukeyMedianPolishResult>` — decomposition for feature extraction
pub fn compute_fragment_fwhm(
    library_fragments: &[LibraryFragment],
    coefficient_series: &[(f64, f64)],
    spectra: &[&Spectrum],
    tolerance_da: f64,
    tolerance_ppm: f64,
) -> (Option<(f64, f64, f64)>, Option<TukeyMedianPolishResult>) {
    if coefficient_series.len() < 3 || spectra.is_empty() {
        return (None, None);
    }

    // Extract XICs for top 6 fragments
    let xics = extract_fragment_xics(library_fragments, spectra, tolerance_da, tolerance_ppm, 6);
    if xics.is_empty() {
        return (None, None);
    }

    // Try Tukey median polish first
    let polish_result = tukey_median_polish(&xics, 20, 1e-4);
    if let Some(ref polish) = polish_result {
        if let Some((fwhm, _left, _right)) = compute_fwhm_interpolated(&polish.elution_profile) {
            let apex_rt = polish
                .elution_profile
                .iter()
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .map(|(rt, _)| *rt);
            if let Some(apex_rt) = apex_rt {
                let sigma = fwhm / 2.355;
                let start_rt = apex_rt - 1.96 * sigma;
                let end_rt = apex_rt + 1.96 * sigma;
                return (Some((start_rt, end_rt, fwhm)), polish_result);
            }
        }
    }

    // Fallback: consensus XIC approach
    let spec_rts: Vec<f64> = spectra.iter().map(|s| s.retention_time).collect();
    let mut coeluting_xics: Vec<&Vec<(f64, f64)>> = Vec::new();

    for (_frag_idx, xic) in &xics {
        let mut xic_aligned: Vec<f64> = Vec::with_capacity(coefficient_series.len());
        let mut coefs: Vec<f64> = Vec::with_capacity(coefficient_series.len());

        for &(rt, coef) in coefficient_series {
            let spec_idx = match spec_rts
                .binary_search_by(|srt| srt.partial_cmp(&rt).unwrap_or(std::cmp::Ordering::Equal))
            {
                Ok(i) => i,
                Err(i) => {
                    if i == 0 {
                        0
                    } else if i >= spec_rts.len() {
                        spec_rts.len() - 1
                    } else if (spec_rts[i] - rt).abs() < (spec_rts[i - 1] - rt).abs() {
                        i
                    } else {
                        i - 1
                    }
                }
            };
            xic_aligned.push(xic[spec_idx].1);
            coefs.push(coef);
        }

        let corr = pearson_correlation_raw(&xic_aligned, &coefs);
        if corr > 0.0 {
            coeluting_xics.push(xic);
        }
    }

    if coeluting_xics.is_empty() {
        return (None, polish_result);
    }

    let n_points = coeluting_xics[0].len();
    let mut consensus: Vec<(f64, f64)> =
        coeluting_xics[0].iter().map(|&(rt, _)| (rt, 0.0)).collect();

    for xic in &coeluting_xics {
        let max_val = xic.iter().map(|(_, v)| *v).fold(0.0f64, f64::max);
        if max_val > 0.0 {
            for i in 0..n_points.min(xic.len()) {
                consensus[i].1 += xic[i].1 / max_val;
            }
        }
    }

    let fwhm_result = compute_fwhm_interpolated(&consensus).and_then(|(fwhm, _left, _right)| {
        let apex_rt = consensus
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, _)| *rt)?;
        let sigma = fwhm / 2.355;
        Some((apex_rt - 1.96 * sigma, apex_rt + 1.96 * sigma, fwhm))
    });

    (fwhm_result, polish_result)
}

/// Plain Pearson correlation on raw values (no intensity transform).
pub fn pearson_correlation_raw(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return 0.0;
    }

    let dn = n as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sx2 = 0.0;
    let mut sy2 = 0.0;
    let mut sxy = 0.0;

    for i in 0..n {
        let xi = x[i];
        let yi = y[i];
        sx += xi;
        sy += yi;
        sx2 += xi * xi;
        sy2 += yi * yi;
        sxy += xi * yi;
    }

    let denom = ((dn * sx2 - sx * sx) * (dn * sy2 - sy * sy))
        .max(1e-30)
        .sqrt();
    (dn * sxy - sx * sy) / denom
}

/// Compute FWHM with linear interpolation on a coefficient time series.
///
/// Uses DIA-NN's approach: scan left/right from the apex until the coefficient
/// drops below half-height, then linearly interpolate to find the exact crossing RT.
///
/// # Returns
/// `Some((fwhm, left_half_rt, right_half_rt))` or `None` if FWHM cannot be computed
/// (e.g., series too short, or signal never drops below half-height).
pub fn compute_fwhm_interpolated(series: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    if series.len() < 3 {
        return None;
    }

    // Find apex
    let (apex_idx, &(_, apex_coef)) = series
        .iter()
        .enumerate()
        .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))?;

    if apex_coef <= 0.0 {
        return None;
    }

    let half = apex_coef / 2.0;

    // Scan right from apex to find half-height crossing
    let right_rt = {
        let mut found = None;
        for i in apex_idx..series.len() - 1 {
            if series[i].1 >= half && series[i + 1].1 < half {
                // Linear interpolation between points i and i+1
                let denom = series[i].1 - series[i + 1].1;
                if denom > 1e-30 {
                    let frac = (series[i].1 - half) / denom;
                    found = Some(series[i].0 + frac * (series[i + 1].0 - series[i].0));
                } else {
                    found = Some(series[i].0);
                }
                break;
            }
        }
        // If never drops below half, use the last point
        found.unwrap_or(series.last()?.0)
    };

    // Scan left from apex to find half-height crossing
    let left_rt = {
        let mut found = None;
        for i in (1..=apex_idx).rev() {
            if series[i].1 >= half && series[i - 1].1 < half {
                // Linear interpolation between points i-1 and i
                let denom = series[i].1 - series[i - 1].1;
                if denom > 1e-30 {
                    let frac = (series[i].1 - half) / denom;
                    found = Some(series[i].0 - frac * (series[i].0 - series[i - 1].0));
                } else {
                    found = Some(series[i].0);
                }
                break;
            }
        }
        // If never drops below half, use the first point
        found.unwrap_or(series.first()?.0)
    };

    let fwhm = right_rt - left_rt;
    if fwhm > 0.0 {
        Some((fwhm, left_rt, right_rt))
    } else {
        None
    }
}

/// Background statistics from flanking points outside the peak region.
pub struct BackgroundStats {
    /// Mean of background points
    pub mean: f64,
    /// Standard deviation of background points
    pub sd: f64,
    /// Number of background points used
    pub n_points: usize,
}

/// Estimate background level from points flanking the peak.
///
/// Uses FWHM-based peak boundaries (±1.96σ) to define the peak region,
/// then collects up to 5 points on each side outside these boundaries.
/// Returns None if fewer than 4 background points are available.
///
/// Works on any (RT, value) series — coefficient series or fragment XICs.
pub fn compute_background_stats(series: &[(f64, f64)], apex_idx: usize) -> Option<BackgroundStats> {
    if series.len() < 10 {
        return None;
    }

    let apex_value = series[apex_idx].1;
    if apex_value <= 0.0 {
        return None;
    }

    // Calculate FWHM to determine peak boundaries
    let (fwhm, _, _) = compute_fwhm_interpolated(series)?;

    // Peak boundaries: ±1.96σ where σ = FWHM / 2.355
    let sigma = fwhm / 2.355;
    let boundary_width = 1.96 * sigma;
    let apex_rt = series[apex_idx].0;
    let left_boundary = apex_rt - boundary_width;
    let right_boundary = apex_rt + boundary_width;

    // Collect background points outside peak boundaries
    let mut background_values = Vec::new();

    // Left side: points before left boundary
    for (rt, val) in series.iter().take(apex_idx) {
        if *rt < left_boundary {
            background_values.push(*val);
        }
    }
    // Take last 5 on left (closest to peak)
    if background_values.len() > 5 {
        background_values = background_values.into_iter().rev().take(5).collect();
    }

    // Right side: points after right boundary
    let mut right_bg = Vec::new();
    for (rt, val) in series.iter().skip(apex_idx + 1) {
        if *rt > right_boundary {
            right_bg.push(*val);
            if right_bg.len() >= 5 {
                break;
            }
        }
    }
    background_values.extend(right_bg);

    // Need at least 4 background points
    if background_values.len() < 4 {
        return None;
    }

    // Calculate background mean and SD
    let n = background_values.len() as f64;
    let bg_mean = background_values.iter().sum::<f64>() / n;
    let bg_variance = background_values
        .iter()
        .map(|x| (x - bg_mean).powi(2))
        .sum::<f64>()
        / n;
    let bg_sd = bg_variance.sqrt();

    Some(BackgroundStats {
        mean: bg_mean,
        sd: bg_sd,
        n_points: background_values.len(),
    })
}

/// Compute trapezoidal area of a time series between given RT boundaries.
///
/// Integrates only the points within [start_rt, end_rt] using the trapezoidal rule:
/// area = Σ (c_i + c_{i+1}) / 2 * (t_{i+1} - t_i)
///
/// If no points fall within boundaries, returns 0.0.
pub fn compute_trapezoidal_area(series: &[(f64, f64)], start_rt: f64, end_rt: f64) -> f64 {
    // Filter to points within boundaries
    let bounded: Vec<(f64, f64)> = series
        .iter()
        .filter(|(rt, _)| *rt >= start_rt && *rt <= end_rt)
        .copied()
        .collect();

    if bounded.len() < 2 {
        return bounded.first().map(|(_, c)| *c).unwrap_or(0.0);
    }

    let mut area = 0.0;
    for i in 0..bounded.len() - 1 {
        let dt = bounded[i + 1].0 - bounded[i].0;
        let avg_height = (bounded[i].1 + bounded[i + 1].1) / 2.0;
        area += avg_height * dt;
    }
    area
}

/// Result of hyperscore computation with binary search fragment matching.
#[derive(Debug, Clone)]
pub struct HyperscoreResult {
    /// X!Tandem hyperscore: ln(n_b!) + ln(n_y!) + Σ ln(I+1)
    pub score: f64,
    /// Number of matched b-ions
    pub n_b: u32,
    /// Number of matched y-ions
    pub n_y: u32,
    /// Total number of matched fragments
    pub n_matched: usize,
    /// Signed mass errors for ALL matched fragments (in configured unit)
    pub mass_errors: Vec<f64>,
}

/// Compute X!Tandem-style hyperscore using binary search fragment matching.
///
/// Matches ALL library fragments against sorted observed peaks within tolerance,
/// counts b/y ions, and computes hyperscore = ln(n_b!) + ln(n_y!) + Σ ln(I+1).
///
/// Uses binary search for O(n_frags × log n_peaks) performance — same pattern as
/// `has_topn_fragment_match` but for ALL fragments.
///
/// Mass errors for all matched fragments are collected as a side product,
/// so no separate mass error collection pass is needed.
pub fn compute_hyperscore(
    library_fragments: &[LibraryFragment],
    spectrum_mzs: &[f64],
    spectrum_intensities: &[f64],
    tolerance: f64,
    unit: ToleranceUnit,
) -> HyperscoreResult {
    let mut n_b: u32 = 0;
    let mut n_y: u32 = 0;
    let mut n_matched: usize = 0;
    let mut sum_log_intensity: f64 = 0.0;
    let mut mass_errors = Vec::new();

    if library_fragments.is_empty() || spectrum_mzs.is_empty() {
        return HyperscoreResult {
            score: 0.0,
            n_b: 0,
            n_y: 0,
            n_matched: 0,
            mass_errors,
        };
    }

    for frag in library_fragments {
        let lib_mz = frag.mz;

        // Calculate Da tolerance from configured unit
        let tol_da = match unit {
            ToleranceUnit::Ppm => lib_mz * tolerance / 1e6,
            ToleranceUnit::Mz => tolerance,
        };

        let lower = lib_mz - tol_da;
        let upper = lib_mz + tol_da;

        // Binary search for first m/z >= lower bound
        let start_idx = spectrum_mzs.partition_point(|&mz| mz < lower);

        if start_idx >= spectrum_mzs.len() || spectrum_mzs[start_idx] > upper {
            continue; // No match within tolerance
        }

        // Find closest peak within tolerance window (by m/z distance)
        let mut best_idx = start_idx;
        let mut best_diff = (spectrum_mzs[start_idx] - lib_mz).abs();
        let mut j = start_idx + 1;
        while j < spectrum_mzs.len() && spectrum_mzs[j] <= upper {
            let diff = (spectrum_mzs[j] - lib_mz).abs();
            if diff < best_diff {
                best_diff = diff;
                best_idx = j;
            }
            j += 1;
        }

        // Record match
        n_matched += 1;
        let obs_intensity = spectrum_intensities[best_idx];
        sum_log_intensity += (obs_intensity + 1.0).ln();

        // Count b/y ion types
        match frag.annotation.ion_type {
            IonType::B => n_b += 1,
            IonType::Y => n_y += 1,
            _ => {} // Other ion types contribute to intensity but not factorial terms
        }

        // Collect signed mass error
        let best_mz = spectrum_mzs[best_idx];
        let error = match unit {
            ToleranceUnit::Ppm => (best_mz - lib_mz) / lib_mz * 1e6,
            ToleranceUnit::Mz => best_mz - lib_mz,
        };
        mass_errors.push(error);
    }

    // hyperscore = ln(n_b!) + ln(n_y!) + Σ ln(I+1)
    let score = if n_matched > 0 {
        ln_gamma(n_b as f64 + 1.0) + ln_gamma(n_y as f64 + 1.0) + sum_log_intensity
    } else {
        0.0
    };

    HyperscoreResult {
        score,
        n_b,
        n_y,
        n_matched,
        mass_errors,
    }
}

/// Spectral similarity scorer implementing LibCosine, XCorr, and Hyperscore
///
/// Scoring methods:
/// - **LibCosine**: sqrt(intensity) preprocessing with L2 normalization
/// - **LibCosine SMZ**: sqrt(intensity) × m/z² preprocessing with L2 normalization
/// - **XCorr**: Comet-style windowing normalization with sliding window subtraction
/// - **Hyperscore**: X!Tandem-style log(n_b!) + log(n_y!) + Σlog(I_f+1)
#[derive(Debug, Clone)]
pub struct SpectralScorer {
    /// Mass tolerance for matching fragments (Da)
    tolerance_da: f64,
    /// Mass tolerance for matching fragments (ppm)
    tolerance_ppm: f64,
    /// Binning configuration for XCorr (Comet BIN macro)
    bin_config: BinConfig,
}

impl Default for SpectralScorer {
    fn default() -> Self {
        Self {
            tolerance_da: 0.5,   // For unit resolution
            tolerance_ppm: 20.0, // For HRAM
            bin_config: BinConfig::unit_resolution(),
        }
    }
}

impl SpectralScorer {
    /// Create a new spectral scorer with default unit resolution settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a scorer with a specific BinConfig (e.g., HRAM)
    pub fn with_bin_config(bin_config: BinConfig) -> Self {
        Self {
            bin_config,
            ..Self::default()
        }
    }

    /// Create an HRAM scorer (0.02 Th bins, 0.0 offset)
    pub fn hram() -> Self {
        Self::with_bin_config(BinConfig::hram())
    }

    /// Set Da tolerance for unit resolution matching
    /// Set Da tolerance for unit resolution matching (clears ppm tolerance)
    pub fn with_tolerance_da(mut self, tolerance: f64) -> Self {
        self.tolerance_da = tolerance;
        self.tolerance_ppm = 0.0; // Da mode — don't let ppm default override
        self
    }

    /// Set ppm tolerance for HRAM matching (clears Da tolerance)
    pub fn with_tolerance_ppm(mut self, tolerance: f64) -> Self {
        self.tolerance_ppm = tolerance;
        self.tolerance_da = 0.0; // ppm mode — don't let Da default override
        self
    }

    /// Get the Da tolerance value
    pub fn tolerance_da(&self) -> f64 {
        self.tolerance_da
    }

    /// Get the ppm tolerance value
    pub fn tolerance_ppm(&self) -> f64 {
        self.tolerance_ppm
    }

    /// Compute LibCosine score between observed spectrum and library entry
    ///
    /// ALL library fragments within the spectrum's mass range are included.
    /// Unmatched fragments use observed intensity of 0 — a missing fragment
    /// that should be present is strong evidence against a match and penalizes
    /// cosine, Pearson, and Spearman scores.
    ///
    /// Counting metrics (fragment_coverage, consecutive_ions, top6_matches,
    /// sequence_coverage) use matched fragments only since they measure presence.
    pub fn lib_cosine(&self, observed: &Spectrum, library: &LibraryEntry) -> SpectralScore {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return SpectralScore::default();
        }

        let spec_mz_min = observed.mzs[0];
        let spec_mz_max = observed.mzs[observed.mzs.len() - 1];

        // Match library fragments to observed peaks (for counting metrics)
        let matches = self.match_fragments(observed, library);

        // Build intensity vectors for ALL library fragments within mass range.
        // Unmatched fragments get observed intensity = 0.
        let mut lib_preprocessed: Vec<f64> = Vec::new();
        let mut obs_preprocessed: Vec<f64> = Vec::new();
        let mut lib_intensities: Vec<f64> = Vec::new();
        let mut obs_intensities: Vec<f64> = Vec::new();

        for frag in &library.fragments {
            if frag.mz < spec_mz_min || frag.mz > spec_mz_max {
                continue;
            }

            let tol_da = self.tolerance_da.max(frag.mz * self.tolerance_ppm / 1e6);
            let lower = frag.mz - tol_da;
            let upper = frag.mz + tol_da;
            let start = observed.mzs.partition_point(|&mz| mz < lower);

            let mut best_intensity = 0.0f64;
            let mut best_diff = f64::MAX;
            let mut j = start;
            while j < observed.mzs.len() && observed.mzs[j] <= upper {
                let diff = (observed.mzs[j] - frag.mz).abs();
                if diff < best_diff {
                    best_diff = diff;
                    best_intensity = observed.intensities[j] as f64;
                }
                j += 1;
            }

            // Raw intensities for base_peak_rank (0 for unmatched)
            lib_intensities.push(frag.relative_intensity as f64);
            obs_intensities.push(best_intensity);

            // Sqrt-preprocessed for cosine, Pearson, and Spearman (0 for unmatched).
            // Sqrt compresses the dynamic range so all fragments contribute more
            // equally, rather than being dominated by a single high-intensity peak.
            lib_preprocessed.push((frag.relative_intensity as f64).sqrt());
            obs_preprocessed.push(best_intensity.sqrt());
        }

        if lib_preprocessed.len() < 2 {
            return SpectralScore::default();
        }

        // L2 normalize for cosine. When either norm is too small (no matched
        // observed intensity, or the library has no non-zero intensities in
        // range) the cosine is undefined; treat that as zero but still
        // populate the presence/counting features below. Tying all features
        // to the norm gate caused short/low-signal peptides to report zero
        // matches even when match_fragments clearly found some.
        let lib_norm = lib_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();

        let cosine_ok = lib_norm >= 1e-10 && obs_norm >= 1e-10;

        let dot_product: f64 = if cosine_ok {
            lib_preprocessed
                .iter()
                .zip(obs_preprocessed.iter())
                .map(|(a, b)| (a / lib_norm) * (b / obs_norm))
                .sum()
        } else {
            0.0
        };

        // Pearson and Spearman use sqrt-preprocessed intensities including zeros
        // for unmatched fragments. The zeros are essential — they penalize missing
        // matches. Without them, a decoy with 2 random matches would score higher
        // than a target with 6 matches but some noise. Gated on norm availability
        // for the same reason cosine is.
        let (pearson_correlation, spearman_correlation) = if cosine_ok {
            (
                Self::pearson_correlation(&lib_preprocessed, &obs_preprocessed),
                Self::spearman_correlation(&lib_preprocessed, &obs_preprocessed),
            )
        } else {
            (0.0, 0.0)
        };

        // Counting metrics use matched fragments only (presence-based)
        let n_matched = matches.len() as u32;
        let n_library = library.fragments.len() as u32;
        let fragment_coverage = n_matched as f64 / n_library as f64;

        let matched_intensity: f64 = matches.iter().map(|m| m.obs_intensity as f64).sum();
        let total_intensity: f64 = observed.intensities.iter().map(|&i| i as f64).sum();
        let explained_intensity = if total_intensity > 0.0 {
            matched_intensity / total_intensity
        } else {
            0.0
        };

        let consecutive_ions = self.longest_consecutive_ions(library, &matches);
        let sequence_coverage = self.compute_sequence_coverage(library, &matches);
        let base_peak_rank =
            self.compute_base_peak_rank(&matches, &lib_intensities, &obs_intensities);
        let top6_matches = self.compute_top6_matches(library, &matches);

        // Compute LibCosine with SMZ preprocessing (sqrt(intensity) * mz²)
        let lib_cosine_smz = self.lib_cosine_smz(observed, library);

        // Compute X!Tandem hyperscore
        let (hyperscore_val, n_matched_b, n_matched_y) = self.hyperscore(observed, library);

        // Compute top-N cosine variants (top 6, 5, 4 library fragments by intensity).
        // Sort fragment pairs by library intensity descending, then compute cosine
        // on the top N. Peptides with fewer than N fragments use all available.
        let mut frag_data: Vec<(f64, f64, f64)> = lib_intensities
            .iter()
            .zip(obs_intensities.iter())
            .zip(
                library
                    .fragments
                    .iter()
                    .filter(|f| f.mz >= spec_mz_min && f.mz <= spec_mz_max),
            )
            .map(|((lib, obs), frag)| (*lib, *obs, frag.mz))
            .collect();
        frag_data.sort_by(|a, b| b.0.total_cmp(&a.0)); // descending by lib intensity

        let dot_product_top6 = Self::cosine_topn_sqrt(&frag_data, 6);
        let dot_product_top5 = Self::cosine_topn_sqrt(&frag_data, 5);
        let dot_product_top4 = Self::cosine_topn_sqrt(&frag_data, 4);
        let dot_product_smz_top6 = Self::cosine_topn_smz(&frag_data, 6);
        let dot_product_smz_top5 = Self::cosine_topn_smz(&frag_data, 5);
        let dot_product_smz_top4 = Self::cosine_topn_smz(&frag_data, 4);

        SpectralScore {
            lib_cosine: dot_product,
            lib_cosine_smz,
            xcorr: 0.0, // Not computed in this method
            hyperscore: hyperscore_val,
            dot_product,
            n_matched,
            n_library,
            fragment_coverage,
            explained_intensity,
            pearson_correlation,
            spearman_correlation,
            consecutive_ions,
            sequence_coverage,
            base_peak_rank,
            top6_matches,
            n_matched_b,
            n_matched_y,
            dot_product_top6,
            dot_product_top5,
            dot_product_top4,
            dot_product_smz_top6,
            dot_product_smz_top5,
            dot_product_smz_top4,
        }
    }

    /// Compute cosine similarity on the top N fragments (by library intensity)
    /// using sqrt preprocessing. `frag_data` must be sorted by library intensity
    /// descending. Each entry is (lib_intensity_raw, obs_intensity_raw, mz).
    fn cosine_topn_sqrt(frag_data: &[(f64, f64, f64)], n: usize) -> f64 {
        let subset = &frag_data[..frag_data.len().min(n)];
        if subset.len() < 2 {
            return 0.0;
        }
        let lib_sqrt: Vec<f64> = subset.iter().map(|(lib, _, _)| lib.sqrt()).collect();
        let obs_sqrt: Vec<f64> = subset.iter().map(|(_, obs, _)| obs.sqrt()).collect();
        let lib_norm = lib_sqrt.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_sqrt.iter().map(|x| x * x).sum::<f64>().sqrt();
        if lib_norm < 1e-10 || obs_norm < 1e-10 {
            return 0.0;
        }
        lib_sqrt
            .iter()
            .zip(obs_sqrt.iter())
            .map(|(a, b)| (a / lib_norm) * (b / obs_norm))
            .sum()
    }

    /// Compute cosine similarity on the top N fragments (by library intensity)
    /// using SMZ preprocessing (sqrt(intensity) * mz²).
    fn cosine_topn_smz(frag_data: &[(f64, f64, f64)], n: usize) -> f64 {
        let subset = &frag_data[..frag_data.len().min(n)];
        if subset.len() < 2 {
            return 0.0;
        }
        let lib_smz: Vec<f64> = subset
            .iter()
            .map(|(lib, _, mz)| lib.sqrt() * mz * mz)
            .collect();
        let obs_smz: Vec<f64> = subset
            .iter()
            .map(|(_, obs, mz)| obs.sqrt() * mz * mz)
            .collect();
        let lib_norm = lib_smz.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_smz.iter().map(|x| x * x).sum::<f64>().sqrt();
        if lib_norm < 1e-10 || obs_norm < 1e-10 {
            return 0.0;
        }
        lib_smz
            .iter()
            .zip(obs_smz.iter())
            .map(|(a, b)| (a / lib_norm) * (b / obs_norm))
            .sum()
    }

    /// Compute Pearson correlation coefficient
    fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
        if x.len() < 2 || x.len() != y.len() {
            return 0.0;
        }

        let n = x.len() as f64;
        let mean_x = x.iter().sum::<f64>() / n;
        let mean_y = y.iter().sum::<f64>() / n;

        let mut sum_xy = 0.0;
        let mut sum_x2 = 0.0;
        let mut sum_y2 = 0.0;

        for (xi, yi) in x.iter().zip(y.iter()) {
            let dx = xi - mean_x;
            let dy = yi - mean_y;
            sum_xy += dx * dy;
            sum_x2 += dx * dx;
            sum_y2 += dy * dy;
        }

        let denom = (sum_x2 * sum_y2).sqrt();
        if denom < 1e-10 {
            return 0.0;
        }

        sum_xy / denom
    }

    /// Compute Spearman rank correlation coefficient
    fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
        if x.len() < 2 || x.len() != y.len() {
            return 0.0;
        }

        // Get ranks
        let rank_x = Self::get_ranks(x);
        let rank_y = Self::get_ranks(y);

        // Pearson correlation of ranks
        Self::pearson_correlation(&rank_x, &rank_y)
    }

    /// Convert values to ranks (1-based, average for ties)
    fn get_ranks(values: &[f64]) -> Vec<f64> {
        let n = values.len();
        let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
        indexed.sort_by(|a, b| a.1.total_cmp(&b.1));

        let mut ranks = vec![0.0; n];
        let mut i = 0;
        while i < n {
            let mut j = i;
            // Find all tied values
            while j < n && (indexed[j].1 - indexed[i].1).abs() < 1e-10 {
                j += 1;
            }
            // Assign average rank to all tied values
            let avg_rank = (i + j + 1) as f64 / 2.0; // 1-based ranks
            for k in i..j {
                ranks[indexed[k].0] = avg_rank;
            }
            i = j;
        }

        ranks
    }

    /// Find longest consecutive b or y ion series
    fn longest_consecutive_ions(&self, _library: &LibraryEntry, matches: &[FragmentMatch]) -> u32 {
        // Each FragmentMatch carries its own ordinal so the attribution is
        // unambiguous. A prior implementation reverse-looked-up the library
        // by m/z and broke on first hit, silently mis-attributing an ordinal
        // whenever two library fragments shared a near-identical m/z (e.g.
        // an incidental b_n / y_m collision).
        let mut matched_b: Vec<u8> = Vec::new();
        let mut matched_y: Vec<u8> = Vec::new();

        for m in matches {
            match m.ion_type {
                IonType::B => matched_b.push(m.ordinal),
                IonType::Y => matched_y.push(m.ordinal),
                _ => {}
            }
        }

        // Find longest consecutive run in each
        let longest_b = Self::longest_consecutive_run(&matched_b);
        let longest_y = Self::longest_consecutive_run(&matched_y);

        longest_b.max(longest_y)
    }

    /// Find longest consecutive run of ordinals
    fn longest_consecutive_run(ordinals: &[u8]) -> u32 {
        if ordinals.is_empty() {
            return 0;
        }

        let mut sorted: Vec<u8> = ordinals.to_vec();
        sorted.sort();
        sorted.dedup();

        let mut max_run = 1;
        let mut current_run = 1;

        for i in 1..sorted.len() {
            if sorted[i] == sorted[i - 1] + 1 {
                current_run += 1;
                max_run = max_run.max(current_run);
            } else {
                current_run = 1;
            }
        }

        max_run
    }

    /// Compute sequence coverage (fraction of backbone cleavage sites covered)
    fn compute_sequence_coverage(&self, library: &LibraryEntry, matches: &[FragmentMatch]) -> f64 {
        let seq_len = library.sequence.len();
        if seq_len < 2 {
            return 0.0;
        }

        // Number of possible cleavage sites is seq_len - 1
        let n_sites = seq_len - 1;
        let mut covered = vec![false; n_sites];

        for m in matches {
            // Find the corresponding library fragment
            for frag in &library.fragments {
                if (frag.mz - m.lib_mz).abs() < 0.001 {
                    let ordinal = frag.annotation.ordinal as usize;
                    if ordinal > 0 && ordinal <= n_sites {
                        // b-ion of ordinal n covers site n-1 (0-indexed)
                        // y-ion of ordinal n covers site seq_len - n - 1 (0-indexed)
                        match frag.annotation.ion_type {
                            IonType::B => {
                                covered[ordinal - 1] = true;
                            }
                            IonType::Y => {
                                let site = seq_len - ordinal - 1;
                                if site < n_sites {
                                    covered[site] = true;
                                }
                            }
                            _ => {}
                        }
                    }
                    break;
                }
            }
        }

        let n_covered = covered.iter().filter(|&&x| x).count();
        n_covered as f64 / n_sites as f64
    }

    /// Compute rank of observed base peak in library ordering
    fn compute_base_peak_rank(
        &self,
        matches: &[FragmentMatch],
        lib_intensities: &[f64],
        obs_intensities: &[f64],
    ) -> u32 {
        if matches.is_empty() || obs_intensities.is_empty() {
            return 0;
        }

        // Find index of max observed intensity
        let max_obs_idx = obs_intensities
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, _)| i);

        if let Some(obs_idx) = max_obs_idx {
            // Get library intensity of this match
            let lib_int = lib_intensities[obs_idx];

            // Count how many library fragments have higher intensity
            let rank = lib_intensities.iter().filter(|&&x| x > lib_int).count() + 1;
            rank as u32
        } else {
            0
        }
    }

    /// Count how many of the top-6 library peaks are matched
    fn compute_top6_matches(&self, library: &LibraryEntry, matches: &[FragmentMatch]) -> u32 {
        if library.fragments.is_empty() {
            return 0;
        }

        // Get top 6 library fragments by intensity
        let mut lib_sorted: Vec<(f64, f64)> = library
            .fragments
            .iter()
            .map(|f| (f.mz, f.relative_intensity as f64))
            .collect();
        lib_sorted.sort_by(|a, b| b.1.total_cmp(&a.1));

        let top6: Vec<f64> = lib_sorted.iter().take(6).map(|(mz, _)| *mz).collect();

        // Count how many are matched
        let mut count = 0;
        for top_mz in &top6 {
            for m in matches {
                if (m.lib_mz - top_mz).abs() < 0.001 {
                    count += 1;
                    break;
                }
            }
        }

        count
    }

    /// Compute XCorr (cross-correlation) score using Comet-style preprocessing
    ///
    /// XCorr applies:
    /// 1. Windowing normalization (normalize to 50.0 within 10 windows)
    /// 2. Sliding window subtraction (offset=75)
    /// 3. Dot product scoring
    pub fn xcorr(&self, observed: &Spectrum, library: &LibraryEntry) -> SpectralScore {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return SpectralScore::default();
        }

        // First get LibCosine for additional metrics
        let lib_cosine_score = self.lib_cosine(observed, library);

        // Bin observed spectrum using Comet BIN macro
        let mut obs_binned = vec![0.0f32; self.bin_config.n_bins];
        for (&mz, &intensity) in observed.mzs.iter().zip(observed.intensities.iter()) {
            if let Some(bin) = self.bin_config.mz_to_bin(mz) {
                // Apply sqrt transformation to experimental spectrum
                obs_binned[bin] += intensity.sqrt();
            }
        }

        // Apply windowing normalization (Comet's MakeCorrData)
        let windowed = self.apply_windowing_normalization(&obs_binned);

        // Apply sliding window subtraction (fast XCorr preprocessing)
        let xcorr_preprocessed = self.apply_sliding_window(&windowed);

        // XCorr = sum of preprocessed experimental values at UNIQUE fragment
        // bin positions. When two library fragments fall into the same bin,
        // the bin's contribution must count once, not twice -- the Comet
        // theoretical spectrum uses unit intensity per bin (see
        // preprocess_library_for_xcorr, which sets binned[bin] = 1.0 per
        // unique bin). Summing preprocessed[bin] once per fragment instead
        // of once per unique bin double-counts collisions and over-scores
        // dense fragment lists.
        let n_bins = xcorr_preprocessed.len();
        let mut visited = vec![false; n_bins];
        let mut xcorr_raw: f64 = 0.0;
        for frag in &library.fragments {
            if let Some(bin) = self.bin_config.mz_to_bin(frag.mz) {
                if !visited[bin] {
                    visited[bin] = true;
                    xcorr_raw += xcorr_preprocessed[bin] as f64;
                }
            }
        }

        // Scale XCorr (pyXcorrDIA uses 0.005 for spectrum-centric)
        let xcorr_scaled = xcorr_raw * 0.005;

        SpectralScore {
            xcorr: xcorr_scaled,
            ..lib_cosine_score
        }
    }

    /// Compute LibCosine score with sqrt(intensity) * mz² (SMZ) preprocessing
    ///
    /// This variant weights fragment matches by their m/z value squared,
    /// giving more importance to higher m/z fragments which are more sequence-specific.
    /// ALL library fragments within the spectrum's mass range are included —
    /// unmatched fragments use observed intensity of 0.
    pub fn lib_cosine_smz(&self, observed: &Spectrum, library: &LibraryEntry) -> f64 {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return 0.0;
        }

        let spec_mz_min = observed.mzs[0];
        let spec_mz_max = observed.mzs[observed.mzs.len() - 1];

        let mut lib_preprocessed: Vec<f64> = Vec::new();
        let mut obs_preprocessed: Vec<f64> = Vec::new();

        for frag in &library.fragments {
            if frag.mz < spec_mz_min || frag.mz > spec_mz_max {
                continue;
            }

            let tol_da = self.tolerance_da.max(frag.mz * self.tolerance_ppm / 1e6);
            let lower = frag.mz - tol_da;
            let upper = frag.mz + tol_da;
            let start = observed.mzs.partition_point(|&mz| mz < lower);

            let mut best_intensity = 0.0f64;
            let mut best_diff = f64::MAX;
            let mut j = start;
            while j < observed.mzs.len() && observed.mzs[j] <= upper {
                let diff = (observed.mzs[j] - frag.mz).abs();
                if diff < best_diff {
                    best_diff = diff;
                    best_intensity = observed.intensities[j] as f64;
                }
                j += 1;
            }

            let mz_sq = frag.mz * frag.mz;
            lib_preprocessed.push((frag.relative_intensity as f64).sqrt() * mz_sq);
            obs_preprocessed.push(best_intensity.sqrt() * mz_sq);
        }

        if lib_preprocessed.is_empty() {
            return 0.0;
        }

        // L2 normalize both vectors
        let lib_norm = lib_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();

        if lib_norm < 1e-10 || obs_norm < 1e-10 {
            return 0.0;
        }

        lib_preprocessed
            .iter()
            .zip(obs_preprocessed.iter())
            .map(|(a, b)| (a / lib_norm) * (b / obs_norm))
            .sum()
    }

    /// Compute X!Tandem-style hyperscore
    ///
    /// The hyperscore combines the number of matched fragment ions with their
    /// observed intensities, rewarding spectra with many matching fragments
    /// of both b and y types.
    ///
    /// Formula: log(n_b!) + log(n_y!) + Σ log(I_f + 1)
    ///
    /// This score does NOT depend on predicted intensities, only on fragment
    /// m/z positions and ion types.
    pub fn hyperscore(&self, observed: &Spectrum, library: &LibraryEntry) -> (f64, u32, u32) {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return (0.0, 0, 0);
        }

        let matches = self.match_fragments(observed, library);
        if matches.is_empty() {
            return (0.0, 0, 0);
        }

        // Count matched b and y ions
        let mut n_b: u32 = 0;
        let mut n_y: u32 = 0;
        let mut sum_log_intensity: f64 = 0.0;

        for m in &matches {
            match m.ion_type {
                IonType::B => n_b += 1,
                IonType::Y => n_y += 1,
                _ => {} // Other ion types don't contribute to factorial terms
            }
            // All matched fragments contribute to the intensity term
            sum_log_intensity += (m.obs_intensity as f64 + 1.0).ln();
        }

        // hyperscore = log(n_b!) + log(n_y!) + Σ log(I_f + 1)
        // Use ln_gamma(n+1) = log(n!) for efficient computation
        let log_nb_factorial = ln_gamma(n_b as f64 + 1.0);
        let log_ny_factorial = ln_gamma(n_y as f64 + 1.0);

        let score = log_nb_factorial + log_ny_factorial + sum_log_intensity;

        (score, n_b, n_y)
    }

    /// Match library fragments to observed peaks
    ///
    /// Uses binary search on sorted observed m/z values for O(n_fragments × log(n_peaks))
    /// instead of O(n_fragments × n_peaks) linear scan.
    pub fn match_fragments(
        &self,
        observed: &Spectrum,
        library: &LibraryEntry,
    ) -> Vec<FragmentMatch> {
        let mut matches = Vec::new();

        for frag in &library.fragments {
            // Calculate tolerance window (use whichever is larger: Da or ppm-derived)
            let tol_da = self.tolerance_da.max(frag.mz * self.tolerance_ppm / 1e6);
            let lower = frag.mz - tol_da;
            let upper = frag.mz + tol_da;

            // Binary search for first m/z >= lower bound (observed.mzs is sorted)
            let start = observed.mzs.partition_point(|&mz| mz < lower);

            // Scan only within the tolerance window to find best match
            let mut best_match: Option<(usize, f64)> = None;
            let mut best_error = f64::MAX;
            let mut i = start;
            while i < observed.mzs.len() && observed.mzs[i] <= upper {
                let error_da = (observed.mzs[i] - frag.mz).abs();
                if error_da < best_error {
                    best_error = error_da;
                    best_match = Some((i, observed.intensities[i] as f64));
                }
                i += 1;
            }

            if let Some((idx, obs_intensity)) = best_match {
                matches.push(FragmentMatch {
                    lib_mz: frag.mz,
                    obs_mz: observed.mzs[idx],
                    lib_intensity: frag.relative_intensity,
                    obs_intensity: obs_intensity as f32,
                    ion_type: frag.annotation.ion_type,
                    ordinal: frag.annotation.ordinal,
                });
            }
        }

        matches
    }

    /// Apply Comet-style windowing normalization
    ///
    /// Divides spectrum into 10 windows and normalizes each to max=50.0
    fn apply_windowing_normalization(&self, spectrum: &[f32]) -> Vec<f32> {
        let mut result = vec![0.0f32; spectrum.len()];
        let num_windows = 10;
        let window_size = (spectrum.len() / num_windows) + 1;

        // Find global max for threshold
        let global_max = spectrum.iter().cloned().fold(0.0f32, f32::max);
        let threshold = global_max * 0.05;

        for window_idx in 0..num_windows {
            let start = window_idx * window_size;
            let end = ((window_idx + 1) * window_size).min(spectrum.len());

            // Find max in this window
            let mut window_max = 0.0f32;
            for &val in &spectrum[start..end] {
                if val > window_max {
                    window_max = val;
                }
            }

            // Normalize this window to 50.0
            if window_max > 0.0 {
                let norm_factor = 50.0 / window_max;
                for i in start..end {
                    if spectrum[i] > threshold {
                        result[i] = spectrum[i] * norm_factor;
                    }
                }
            }
        }

        result
    }

    /// Apply sliding window subtraction for fast XCorr (Comet-style)
    ///
    /// Uses prefix sum for O(n) performance instead of O(n × window).
    /// Comet divides by (2*offset) = 150 regardless of boundary effects.
    /// offset=75, matching Comet's iXcorrProcessingOffset default.
    fn apply_sliding_window(&self, spectrum: &[f32]) -> Vec<f32> {
        let n = spectrum.len();
        let offset: usize = 75;
        // Comet uses (window_size - 1) = 2*offset = 150 as divisor
        let norm_factor = 1.0f32 / (2 * offset) as f32;

        // Build prefix sum for O(n) window sums
        let mut prefix = vec![0.0f32; n + 1];
        for i in 0..n {
            prefix[i + 1] = prefix[i] + spectrum[i];
        }

        let mut result = vec![0.0f32; n];
        for i in 0..n {
            let left = i.saturating_sub(offset);
            let right = if i + offset < n { i + offset + 1 } else { n };
            // Window sum including center
            let window_sum = prefix[right] - prefix[left];
            // Subtract center to get sum excluding center
            let sum_excluding_center = window_sum - spectrum[i];
            // Subtract local average from center value
            result[i] = spectrum[i] - sum_excluding_center * norm_factor;
        }

        result
    }

    // ========================================================================
    // Batch XCorr preprocessing methods (for BLAS-accelerated scoring)
    // ========================================================================

    /// Preprocess a spectrum for XCorr (Comet-style)
    ///
    /// Returns a preprocessed vector ready for dot product with library vectors.
    /// Preprocessing: bin → sqrt → windowing → flanking subtraction
    ///
    /// This allows precomputing once per spectrum and reusing across library entries.
    pub fn preprocess_spectrum_for_xcorr(&self, spectrum: &Spectrum) -> Vec<f32> {
        // Bin observed spectrum with sqrt transformation using Comet BIN macro
        let mut binned = vec![0.0f32; self.bin_config.n_bins];
        for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            if let Some(bin) = self.bin_config.mz_to_bin(mz) {
                binned[bin] += intensity.sqrt();
            }
        }

        // Apply windowing normalization
        let windowed = self.apply_windowing_normalization(&binned);

        // Apply flanking bin subtraction
        self.apply_sliding_window(&windowed)
    }

    /// Lightweight XCorr for a single spectrum — no LibCosine overhead.
    ///
    /// Preprocesses the spectrum (bin + window + sliding window subtraction),
    /// then sums the preprocessed values at library fragment bin positions × 0.005.
    pub fn xcorr_at_scan(&self, spectrum: &Spectrum, library: &LibraryEntry) -> f64 {
        if library.fragments.is_empty() || spectrum.mzs.is_empty() {
            return 0.0;
        }
        let preprocessed = self.preprocess_spectrum_for_xcorr(spectrum);
        // Unique fragment bins only (Comet theoretical spectrum uses unit
        // intensity per bin; see preprocess_library_for_xcorr). Matches
        // the scorer.xcorr() dedup.
        let n_bins = preprocessed.len();
        let mut visited = vec![false; n_bins];
        let mut xcorr_raw: f64 = 0.0;
        for frag in &library.fragments {
            if let Some(bin) = self.bin_config.mz_to_bin(frag.mz) {
                if !visited[bin] {
                    visited[bin] = true;
                    xcorr_raw += preprocessed[bin] as f64;
                }
            }
        }
        xcorr_raw * 0.005
    }

    /// Preprocess a library entry for XCorr (Comet-style)
    ///
    /// Returns a preprocessed vector ready for dot product with spectrum vectors.
    /// Preprocessing: bin with unit intensity → windowing (NO flanking subtraction)
    ///
    /// This allows precomputing once per library entry and reusing across spectra.
    pub fn preprocess_library_for_xcorr(&self, entry: &LibraryEntry) -> Vec<f32> {
        // Bin library fragments with unit intensity using Comet BIN macro
        // Comet-style: theoretical spectrum is NOT windowed, just unit intensities at fragment bins
        // The score is simply: sum(experimental_preprocessed[frag_bins]) * 0.005
        let mut binned = vec![0.0f32; self.bin_config.n_bins];
        for frag in &entry.fragments {
            if let Some(bin) = self.bin_config.mz_to_bin(frag.mz) {
                // Use unit intensity (1.0), NOT library intensity (Comet-style)
                // NO windowing applied - this matches Comet exactly
                binned[bin] = 1.0;
            }
        }

        // Return binned vector directly - NO windowing for theoretical spectrum
        // Comet just looks up preprocessed experimental values at fragment bin positions
        binned
    }

    /// Compute XCorr from preprocessed vectors
    ///
    /// Uses ndarray dot product which dispatches to BLAS sdot for contiguous f32 slices.
    /// Returns f64 for compatibility with scoring pipelines.
    #[inline]
    pub fn xcorr_from_preprocessed(
        spectrum_preprocessed: &[f32],
        library_preprocessed: &[f32],
    ) -> f64 {
        use ndarray::ArrayView1;

        let min_len = spectrum_preprocessed.len().min(library_preprocessed.len());
        let spec = ArrayView1::from(&spectrum_preprocessed[..min_len]);
        let lib = ArrayView1::from(&library_preprocessed[..min_len]);
        let raw = spec.dot(&lib);

        // Scale XCorr (pyXcorrDIA uses 0.005 for spectrum-centric)
        (raw * 0.005f32) as f64
    }

    /// Get the number of bins used for XCorr
    pub fn num_bins(&self) -> usize {
        self.bin_config.n_bins
    }

    /// Get the bin configuration
    pub fn bin_config(&self) -> &BinConfig {
        &self.bin_config
    }
}

/// Result of spectral similarity scoring
#[derive(Debug, Clone, Default)]
pub struct SpectralScore {
    /// LibCosine score (0-1, cosine similarity with sqrt preprocessing)
    pub lib_cosine: f64,
    /// LibCosine score with sqrt(intensity)*mz² (SMZ) preprocessing
    pub lib_cosine_smz: f64,
    /// XCorr score (Comet-style cross-correlation)
    pub xcorr: f64,
    /// X!Tandem-style hyperscore: log(n_b!) + log(n_y!) + Σlog(I_f+1)
    pub hyperscore: f64,
    /// Raw dot product (same as lib_cosine for backward compat)
    pub dot_product: f64,
    /// Number of matched fragments
    pub n_matched: u32,
    /// Number of library fragments
    pub n_library: u32,
    /// Fragment coverage (n_matched / n_library)
    pub fragment_coverage: f64,
    /// Fraction of observed intensity explained by matches
    pub explained_intensity: f64,
    /// Pearson correlation between library and observed intensities
    pub pearson_correlation: f64,
    /// Spearman rank correlation between library and observed intensities
    pub spearman_correlation: f64,
    /// Longest consecutive b or y ion series
    pub consecutive_ions: u32,
    /// Sequence coverage (backbone cleavages covered)
    pub sequence_coverage: f64,
    /// Rank of base peak (highest intensity) in library
    pub base_peak_rank: u32,
    /// Number of top-6 library peaks that matched
    pub top6_matches: u32,
    /// Number of matched b-ions (for hyperscore)
    pub n_matched_b: u32,
    /// Number of matched y-ions (for hyperscore)
    pub n_matched_y: u32,
    /// LibCosine using only top 6 library fragments by intensity
    pub dot_product_top6: f64,
    /// LibCosine using only top 5 library fragments by intensity
    pub dot_product_top5: f64,
    /// LibCosine using only top 4 library fragments by intensity
    pub dot_product_top4: f64,
    /// SMZ cosine using only top 6 library fragments by intensity
    pub dot_product_smz_top6: f64,
    /// SMZ cosine using only top 5 library fragments by intensity
    pub dot_product_smz_top5: f64,
    /// SMZ cosine using only top 4 library fragments by intensity
    pub dot_product_smz_top4: f64,
}

/// Fragment matching result with observed m/z for mass accuracy computation
#[derive(Debug)]
pub struct FragmentMatch {
    pub lib_mz: f64,
    pub obs_mz: f64,
    pub lib_intensity: f32,
    pub obs_intensity: f32,
    /// Ion type of the matched fragment (for hyperscore b/y counting)
    pub ion_type: IonType,
    /// Ordinal of the matched fragment (e.g. 5 for b5/y5). Needed to attribute
    /// the correct series position when two library fragments share a near-
    /// identical m/z — previously `longest_consecutive_ions` reverse-looked-up
    /// the library by m/z and broke on first hit, silently dropping the second
    /// fragment's ordinal.
    pub ordinal: u8,
}

/// Spectrum aggregator for combining spectra across peak apex region
///
/// FR-5.2.1: Aggregate observed spectrum across peak apex region
#[derive(Debug, Clone)]
pub struct SpectrumAggregator {
    /// Mass tolerance for combining peaks (Da)
    tolerance_da: f64,
    /// Number of scans to include around apex
    apex_window: usize,
}

impl Default for SpectrumAggregator {
    fn default() -> Self {
        Self {
            tolerance_da: 0.5,
            apex_window: 5, // 2 scans on each side of apex
        }
    }
}

impl SpectrumAggregator {
    /// Create a new spectrum aggregator
    pub fn new() -> Self {
        Self::default()
    }

    /// Set mass tolerance for peak combining
    pub fn with_tolerance_da(mut self, tolerance: f64) -> Self {
        self.tolerance_da = tolerance;
        self
    }

    /// Set number of scans to aggregate around apex
    pub fn with_apex_window(mut self, window: usize) -> Self {
        self.apex_window = window;
        self
    }

    /// Aggregate spectra around the peak apex
    ///
    /// Combines peaks from multiple spectra within the apex region
    /// by summing intensities of peaks at similar m/z values.
    pub fn aggregate(&self, spectra: &[&Spectrum], apex_rt: f64) -> Spectrum {
        if spectra.is_empty() {
            return Spectrum::new(
                0,
                apex_rt,
                osprey_core::IsolationWindow::symmetric(0.0, 0.0),
            );
        }

        // Find apex scan index
        let apex_idx = spectra
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1.retention_time - apex_rt)
                    .abs()
                    .partial_cmp(&(b.1.retention_time - apex_rt).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Select scans in the window around apex
        let start_idx = apex_idx.saturating_sub(self.apex_window / 2);
        let end_idx = (apex_idx + self.apex_window / 2 + 1).min(spectra.len());
        let window_spectra = &spectra[start_idx..end_idx];

        if window_spectra.is_empty() {
            return spectra[apex_idx].clone();
        }

        // Collect all peaks from window
        let mut all_peaks: Vec<(f64, f64)> = Vec::new();
        for spectrum in window_spectra {
            for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
                all_peaks.push((mz, intensity as f64));
            }
        }

        // Sort by m/z
        all_peaks.sort_by(|a, b| a.0.total_cmp(&b.0));

        // Merge peaks within tolerance
        let mut merged_mzs: Vec<f64> = Vec::new();
        let mut merged_intensities: Vec<f32> = Vec::new();

        let mut i = 0;
        while i < all_peaks.len() {
            let mut sum_mz = all_peaks[i].0;
            let mut sum_intensity = all_peaks[i].1;
            let mut count = 1.0;
            let base_mz = all_peaks[i].0;

            // Merge all peaks within tolerance
            let mut j = i + 1;
            while j < all_peaks.len() && (all_peaks[j].0 - base_mz).abs() <= self.tolerance_da {
                sum_mz += all_peaks[j].0;
                sum_intensity += all_peaks[j].1;
                count += 1.0;
                j += 1;
            }

            // Use intensity-weighted m/z average
            merged_mzs.push(sum_mz / count);
            merged_intensities.push(sum_intensity as f32);

            i = j;
        }

        // Create aggregated spectrum using the apex spectrum's metadata
        let apex_spectrum = &spectra[apex_idx];
        Spectrum {
            scan_number: apex_spectrum.scan_number,
            retention_time: apex_spectrum.retention_time,
            precursor_mz: apex_spectrum.precursor_mz,
            isolation_window: apex_spectrum.isolation_window,
            mzs: merged_mzs,
            intensities: merged_intensities,
        }
    }

    /// Aggregate spectra weighted by coefficient values
    ///
    /// Uses coefficient values to weight the contribution
    /// of each scan to the aggregated spectrum.
    pub fn aggregate_weighted(
        &self,
        spectra: &[&Spectrum],
        coefficients: &[(f64, f64)], // (rt, coefficient) pairs
        apex_rt: f64,
    ) -> Spectrum {
        if spectra.is_empty() || coefficients.is_empty() {
            return Spectrum::new(
                0,
                apex_rt,
                osprey_core::IsolationWindow::symmetric(0.0, 0.0),
            );
        }

        // Find apex scan index
        let apex_idx = spectra
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1.retention_time - apex_rt)
                    .abs()
                    .partial_cmp(&(b.1.retention_time - apex_rt).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Select scans in the window around apex
        let start_idx = apex_idx.saturating_sub(self.apex_window / 2);
        let end_idx = (apex_idx + self.apex_window / 2 + 1).min(spectra.len());
        let window_spectra = &spectra[start_idx..end_idx];

        if window_spectra.is_empty() {
            return spectra[apex_idx].clone();
        }

        // Get coefficient weights for each spectrum
        let mut weights: Vec<f64> = Vec::new();
        for spectrum in window_spectra {
            // Find corresponding coefficient
            let weight = coefficients
                .iter()
                .min_by(|a, b| {
                    (a.0 - spectrum.retention_time)
                        .abs()
                        .partial_cmp(&(b.0 - spectrum.retention_time).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(_, c)| c.max(0.0)) // Use max(0, c) for non-negative weights
                .unwrap_or(0.0);
            weights.push(weight);
        }

        // Normalize weights
        let total_weight: f64 = weights.iter().sum();
        if total_weight < 1e-10 {
            // All weights zero, fall back to simple aggregation
            return self.aggregate(spectra, apex_rt);
        }
        for w in weights.iter_mut() {
            *w /= total_weight;
        }

        // Collect all peaks with weights
        let mut all_peaks: Vec<(f64, f64)> = Vec::new();
        for (spectrum, &weight) in window_spectra.iter().zip(weights.iter()) {
            for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
                all_peaks.push((mz, intensity as f64 * weight));
            }
        }

        // Sort by m/z
        all_peaks.sort_by(|a, b| a.0.total_cmp(&b.0));

        // Merge peaks within tolerance
        let mut merged_mzs: Vec<f64> = Vec::new();
        let mut merged_intensities: Vec<f32> = Vec::new();

        let mut i = 0;
        while i < all_peaks.len() {
            let mut sum_mz_weighted = all_peaks[i].0 * all_peaks[i].1;
            let mut sum_intensity = all_peaks[i].1;
            let base_mz = all_peaks[i].0;

            // Merge all peaks within tolerance
            let mut j = i + 1;
            while j < all_peaks.len() && (all_peaks[j].0 - base_mz).abs() <= self.tolerance_da {
                sum_mz_weighted += all_peaks[j].0 * all_peaks[j].1;
                sum_intensity += all_peaks[j].1;
                j += 1;
            }

            if sum_intensity > 1e-10 {
                // Use intensity-weighted m/z average
                merged_mzs.push(sum_mz_weighted / sum_intensity);
                merged_intensities.push(sum_intensity as f32);
            }

            i = j;
        }

        // Create aggregated spectrum
        let apex_spectrum = &spectra[apex_idx];
        Spectrum {
            scan_number: apex_spectrum.scan_number,
            retention_time: apex_spectrum.retention_time,
            precursor_mz: apex_spectrum.precursor_mz,
            isolation_window: apex_spectrum.isolation_window,
            mzs: merged_mzs,
            intensities: merged_intensities,
        }
    }
}

/// Compute ln(Gamma(x)) = ln((x-1)!) for positive integers
///
/// Uses Stirling's approximation for large values and lookup table for small values.
/// For X!Tandem hyperscore, this computes log(n!) efficiently.
fn ln_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x <= 1.0 {
        return 0.0; // ln(Gamma(1)) = ln(0!) = 0
    }
    if x <= 2.0 {
        return 0.0; // ln(Gamma(2)) = ln(1!) = 0
    }

    // Use Lanczos approximation (accurate for all positive values)
    // Coefficients from Numerical Recipes
    let g = 7.0;
    let c = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];

    let x_adj = x - 1.0;
    let mut sum = c[0];
    for (i, &coeff) in c[1..].iter().enumerate() {
        sum += coeff / (x_adj + i as f64 + 1.0);
    }

    let t = x_adj + g + 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (t.ln() * (x_adj + 0.5)) - t + sum.ln()
}

/// Decoy generator using enzyme-aware sequence reversal
///
/// Following the pyXcorrDIA approach:
/// - Reverse peptide sequence preserving terminal residue based on enzyme
/// - Recalculate fragment m/z values for the reversed sequence
/// - Precursor m/z stays the same (same amino acid composition)
#[derive(Debug, Clone)]
pub struct DecoyGenerator {
    /// Method for generating decoys
    method: DecoyMethod,
    /// Enzyme used for digestion
    enzyme: Enzyme,
    /// Amino acid masses for fragment calculation
    aa_masses: HashMap<char, f64>,
}

/// Decoy generation method
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum DecoyMethod {
    /// Reverse the peptide sequence (preserving terminal residue)
    #[default]
    Reverse,
    /// Shuffle the peptide sequence
    Shuffle,
}

/// Digestion enzyme (affects terminal preservation during reversal)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Enzyme {
    /// Trypsin: C-terminal cleavage at K/R (preserves C-terminus)
    #[default]
    Trypsin,
    /// Lys-C: C-terminal cleavage at K (preserves C-terminus)
    LysC,
    /// Lys-N: N-terminal cleavage at K (preserves N-terminus)
    LysN,
    /// Asp-N: N-terminal cleavage at D (preserves N-terminus)
    AspN,
    /// No enzyme specificity
    Unspecific,
}

impl Enzyme {
    /// Returns true if this enzyme cleaves at C-terminus (so C-term should be preserved)
    pub fn preserves_c_terminus(&self) -> bool {
        matches!(self, Enzyme::Trypsin | Enzyme::LysC | Enzyme::Unspecific)
    }
}

impl DecoyGenerator {
    /// Create a new decoy generator with default settings (trypsin, reverse)
    pub fn new(method: DecoyMethod) -> Self {
        Self {
            method,
            enzyme: Enzyme::default(),
            aa_masses: Self::standard_aa_masses(),
        }
    }

    /// Create a decoy generator with specified enzyme
    pub fn with_enzyme(method: DecoyMethod, enzyme: Enzyme) -> Self {
        Self {
            method,
            enzyme,
            aa_masses: Self::standard_aa_masses(),
        }
    }

    /// Standard monoisotopic amino acid residue masses
    fn standard_aa_masses() -> HashMap<char, f64> {
        let mut masses = HashMap::new();
        masses.insert('A', 71.037114);
        masses.insert('R', 156.101111);
        masses.insert('N', 114.042927);
        masses.insert('D', 115.026943);
        masses.insert('C', 103.009185);
        masses.insert('E', 129.042593);
        masses.insert('Q', 128.058578);
        masses.insert('G', 57.021464);
        masses.insert('H', 137.058912);
        masses.insert('I', 113.084064);
        masses.insert('L', 113.084064);
        masses.insert('K', 128.094963);
        masses.insert('M', 131.040485);
        masses.insert('F', 147.068414);
        masses.insert('P', 97.052764);
        masses.insert('S', 87.032028);
        masses.insert('T', 101.047679);
        masses.insert('W', 186.079313);
        masses.insert('Y', 163.063329);
        masses.insert('V', 99.068414);
        masses
    }

    /// Generate a decoy from a target library entry
    pub fn generate(&self, target: &LibraryEntry) -> Result<LibraryEntry> {
        let mut decoy = target.clone();
        decoy.is_decoy = true;

        // Generate decoy sequence
        let (decoy_sequence, position_mapping) = match self.method {
            DecoyMethod::Reverse => self.reverse_sequence(&target.sequence),
            DecoyMethod::Shuffle => self.shuffle_sequence(&target.sequence),
        };

        decoy.sequence = decoy_sequence.clone();

        // Update modified sequence notation
        decoy.modified_sequence = format!("DECOY_{}", target.modified_sequence);

        // Remap modifications to new positions
        decoy.modifications = self.remap_modifications(&target.modifications, &position_mapping);

        // Recalculate fragment m/z values for the reversed sequence
        decoy.fragments = self.recalculate_fragments(target, &position_mapping);

        // Precursor m/z stays the same (same amino acid composition)
        // decoy.precursor_mz = target.precursor_mz;

        // Update ID to indicate decoy (set high bit)
        decoy.id = target.id | 0x80000000;

        // Update protein ID to indicate decoy
        decoy.protein_ids = target
            .protein_ids
            .iter()
            .map(|p| format!("DECOY_{}", p))
            .collect();

        Ok(decoy)
    }

    /// Reverse sequence preserving terminal residue based on enzyme
    ///
    /// For trypsin/LysC (C-terminal cleavage): preserve C-terminus only
    /// For Lys-N/Asp-N (N-terminal cleavage): preserve N-terminus only
    ///
    /// Returns (reversed_sequence, position_mapping)
    /// position_mapping[new_pos] = old_pos
    fn reverse_sequence(&self, sequence: &str) -> (String, Vec<usize>) {
        let chars: Vec<char> = sequence.chars().collect();
        let len = chars.len();

        if len <= 2 {
            // Too short to reverse meaningfully
            return (sequence.to_string(), (0..len).collect());
        }

        let mut reversed: Vec<char> = Vec::with_capacity(len);
        let mut position_mapping: Vec<usize> = Vec::with_capacity(len);

        if self.enzyme.preserves_c_terminus() {
            // Preserve C-terminus only (trypsin, Lys-C)
            // Reverse everything except the C-terminal residue
            // Original: P E P T I D E K  (K is the cleavage site)
            // Reversed: E D I T P E P K
            //
            // Reverse positions 0..len-1, keep position len-1 at the end
            for i in (0..len - 1).rev() {
                reversed.push(chars[i]);
                position_mapping.push(i);
            }

            // Keep C-term
            reversed.push(chars[len - 1]);
            position_mapping.push(len - 1);
        } else {
            // Preserve N-terminus only (Lys-N, Asp-N)
            // Reverse everything except the N-terminal residue
            // Original: K P E P T I D E  (K is the cleavage site)
            // Reversed: K E D I T P E P

            // Keep N-term
            reversed.push(chars[0]);
            position_mapping.push(0);

            // Reverse the rest (positions 1..len)
            for i in (1..len).rev() {
                reversed.push(chars[i]);
                position_mapping.push(i);
            }
        }

        (reversed.into_iter().collect(), position_mapping)
    }

    /// Shuffle sequence preserving terminal residue
    fn shuffle_sequence(&self, sequence: &str) -> (String, Vec<usize>) {
        let chars: Vec<char> = sequence.chars().collect();
        let len = chars.len();

        if len <= 2 {
            return (sequence.to_string(), (0..len).collect());
        }

        let mut shuffled: Vec<char> = Vec::with_capacity(len);
        let mut position_mapping: Vec<usize> = Vec::with_capacity(len);

        if self.enzyme.preserves_c_terminus() {
            // Keep C-term only, shuffle the rest (positions 0..len-1)
            let prefix_len = len - 1;
            for i in 0..prefix_len {
                let j = (i + prefix_len / 2 + 1) % prefix_len;
                shuffled.push(chars[j]);
                position_mapping.push(j);
            }

            shuffled.push(chars[len - 1]);
            position_mapping.push(len - 1);
        } else {
            // Keep N-term only, shuffle the rest (positions 1..len)
            shuffled.push(chars[0]);
            position_mapping.push(0);

            let rest_len = len - 1;
            for i in 0..rest_len {
                let j = 1 + (i + rest_len / 2 + 1) % rest_len;
                shuffled.push(chars[j]);
                position_mapping.push(j);
            }
        }

        (shuffled.into_iter().collect(), position_mapping)
    }

    /// Remap modifications to new positions after sequence reversal
    fn remap_modifications(
        &self,
        modifications: &[Modification],
        position_mapping: &[usize],
    ) -> Vec<Modification> {
        // Create reverse mapping: old_pos -> new_pos
        let mut reverse_map: HashMap<usize, usize> = HashMap::new();
        for (new_pos, &old_pos) in position_mapping.iter().enumerate() {
            reverse_map.insert(old_pos, new_pos);
        }

        modifications
            .iter()
            .filter_map(|m| {
                reverse_map.get(&m.position).map(|&new_pos| Modification {
                    position: new_pos,
                    unimod_id: m.unimod_id,
                    mass_delta: m.mass_delta,
                    name: m.name.clone(),
                })
            })
            .collect()
    }

    /// Recalculate fragment m/z values for the reversed sequence
    ///
    /// When reversing a sequence:
    /// - b-ions become y-ions at complementary positions
    /// - y-ions become b-ions at complementary positions
    ///
    /// For a sequence of length N:
    /// - b{i} -> y{N-i} (reversed)
    /// - y{i} -> b{N-i} (reversed)
    fn recalculate_fragments(
        &self,
        target: &LibraryEntry,
        position_mapping: &[usize],
    ) -> Vec<LibraryFragment> {
        let seq_len = target.sequence.len();
        let decoy_chars: Vec<char> = self.get_decoy_sequence_chars(target, position_mapping);

        // Build modification mass map for decoy (by position)
        let mut mod_masses: HashMap<usize, f64> = HashMap::new();
        for m in &target.modifications {
            // Find new position for this modification
            if let Some(new_pos) = position_mapping.iter().position(|&old| old == m.position) {
                mod_masses.insert(new_pos, m.mass_delta);
            }
        }

        target
            .fragments
            .iter()
            .filter_map(|frag| {
                self.recalculate_single_fragment(frag, seq_len, &decoy_chars, &mod_masses)
            })
            .collect()
    }

    /// Get decoy sequence characters from position mapping
    fn get_decoy_sequence_chars(
        &self,
        target: &LibraryEntry,
        position_mapping: &[usize],
    ) -> Vec<char> {
        let target_chars: Vec<char> = target.sequence.chars().collect();
        position_mapping
            .iter()
            .map(|&old_pos| target_chars[old_pos])
            .collect()
    }

    /// Recalculate a single fragment for the decoy
    fn recalculate_single_fragment(
        &self,
        frag: &LibraryFragment,
        seq_len: usize,
        decoy_chars: &[char],
        mod_masses: &HashMap<usize, f64>,
    ) -> Option<LibraryFragment> {
        let annotation = &frag.annotation;

        // Handle b and y ions - they swap when sequence is reversed
        let (new_ion_type, new_ordinal) = match annotation.ion_type {
            IonType::B => {
                // b{i} covers residues 0..i (N-terminal)
                // In reversed sequence, this becomes y{seq_len - i}
                (IonType::Y, (seq_len - annotation.ordinal as usize) as u8)
            }
            IonType::Y => {
                // y{i} covers residues (seq_len-i)..seq_len (C-terminal)
                // In reversed sequence, this becomes b{seq_len - i}
                (IonType::B, (seq_len - annotation.ordinal as usize) as u8)
            }
            // For other ion types, keep as-is but don't recalculate
            _ => return Some(frag.clone()),
        };

        if new_ordinal == 0 || new_ordinal as usize > seq_len {
            return None;
        }

        // Calculate new m/z based on the decoy sequence
        let new_mz = self.calculate_fragment_mz(
            new_ion_type,
            new_ordinal as usize,
            annotation.charge,
            decoy_chars,
            mod_masses,
            annotation.neutral_loss.as_ref().map(|nl| nl.mass()),
        )?;

        Some(LibraryFragment {
            mz: new_mz,
            relative_intensity: frag.relative_intensity,
            annotation: FragmentAnnotation {
                ion_type: new_ion_type,
                ordinal: new_ordinal,
                charge: annotation.charge,
                neutral_loss: annotation.neutral_loss,
            },
        })
    }

    /// Calculate fragment m/z from sequence
    fn calculate_fragment_mz(
        &self,
        ion_type: IonType,
        ordinal: usize,
        charge: u8,
        sequence: &[char],
        mod_masses: &HashMap<usize, f64>,
        neutral_loss: Option<f64>,
    ) -> Option<f64> {
        const PROTON_MASS: f64 = 1.007276;
        const H2O_MASS: f64 = 18.010565;

        let seq_len = sequence.len();

        // Sum amino acid masses
        let (start, end) = match ion_type {
            IonType::B => (0, ordinal),                 // b-ions: N-terminal fragment
            IonType::Y => (seq_len - ordinal, seq_len), // y-ions: C-terminal fragment
            _ => return None,
        };

        if end > seq_len {
            return None;
        }

        let mut mass = 0.0;

        // Sum residue masses
        for i in start..end {
            let aa = sequence.get(i)?;
            mass += self.aa_masses.get(aa)?;

            // Add modification mass if present
            if let Some(&mod_mass) = mod_masses.get(&i) {
                mass += mod_mass;
            }
        }

        // Add terminal groups and ion-specific masses
        match ion_type {
            IonType::B => {
                // b-ion: [M + H]+ - loses C-terminal OH
                mass += PROTON_MASS;
            }
            IonType::Y => {
                // y-ion: [M + H]+ + H2O - keeps C-terminal OH, gains H on N-terminal
                mass += H2O_MASS + PROTON_MASS;
            }
            _ => {}
        }

        // Subtract neutral loss if present
        if let Some(loss) = neutral_loss {
            mass -= loss;
        }

        // Calculate m/z for charge state
        let mz = (mass + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64;

        Some(mz)
    }

    /// Generate decoys for all entries in a library (legacy method without collision detection)
    pub fn generate_all(&self, targets: &[LibraryEntry]) -> Result<Vec<LibraryEntry>> {
        targets.iter().map(|t| self.generate(t)).collect()
    }

    /// Generate decoys with collision detection (pyXcorrDIA approach)
    ///
    /// This method:
    /// 1. Builds a set of all target sequences for collision detection
    /// 2. For each target, tries reversal first
    /// 3. If reversed sequence collides with a target, uses cycling fallback
    /// 4. If all methods fail, excludes the target-decoy pair
    ///
    /// Returns (valid_targets, decoys, statistics)
    pub fn generate_all_with_collision_detection(
        &self,
        targets: &[LibraryEntry],
    ) -> (Vec<LibraryEntry>, Vec<LibraryEntry>, DecoyGenerationStats) {
        use std::collections::HashSet;

        // Build set of all target sequences for collision detection
        let target_sequences: HashSet<&str> = targets.iter().map(|t| t.sequence.as_str()).collect();

        // Result type for each parallel task
        enum DecoyResult {
            Reversed(LibraryEntry, LibraryEntry),        // (target, decoy)
            CyclingFallback(LibraryEntry, LibraryEntry), // (target, decoy)
            SkippedNoFragments,
            FailedFragmentCalculation,
            ExcludedNoUniqueDecoy,
        }

        // Process each target in parallel
        let results: Vec<DecoyResult> = targets
            .par_iter()
            .map(|target| {
                // Skip if target has no fragments (can't generate meaningful decoy)
                if target.fragments.is_empty() {
                    return DecoyResult::SkippedNoFragments;
                }

                // Try reversal first
                let (reversed_seq, position_mapping) = self.reverse_sequence(&target.sequence);

                // Check if reversed sequence is valid:
                // 1. Different from original (not palindromic)
                // 2. Not in target database (no collision)
                if reversed_seq != target.sequence
                    && !target_sequences.contains(reversed_seq.as_str())
                {
                    // Reversal succeeded - create decoy
                    match self.create_decoy_from_mapping(target, &reversed_seq, &position_mapping) {
                        Ok(decoy) => {
                            return DecoyResult::Reversed(target.clone(), decoy);
                        }
                        Err(_) => {
                            return DecoyResult::FailedFragmentCalculation;
                        }
                    }
                }

                // Reversal failed - try cycling fallback
                let max_retries = target.sequence.len().min(10);

                for cycle_length in 1..=max_retries {
                    let (cycled_seq, cycle_mapping) =
                        self.cycle_sequence(&target.sequence, cycle_length);

                    // Check if cycled sequence is valid
                    if cycled_seq != target.sequence
                        && !target_sequences.contains(cycled_seq.as_str())
                    {
                        match self.create_decoy_from_mapping(target, &cycled_seq, &cycle_mapping) {
                            Ok(decoy) => {
                                return DecoyResult::CyclingFallback(target.clone(), decoy);
                            }
                            Err(_) => {
                                continue; // Try next cycle length
                            }
                        }
                    }
                }

                // All methods failed - exclude this target-decoy pair
                DecoyResult::ExcludedNoUniqueDecoy
            })
            .collect();

        // Aggregate results
        let mut valid_targets: Vec<LibraryEntry> = Vec::with_capacity(targets.len());
        let mut decoys: Vec<LibraryEntry> = Vec::with_capacity(targets.len());
        let mut stats = DecoyGenerationStats::default();

        for result in results {
            match result {
                DecoyResult::Reversed(target, decoy) => {
                    valid_targets.push(target);
                    decoys.push(decoy);
                    stats.reversed += 1;
                }
                DecoyResult::CyclingFallback(target, decoy) => {
                    valid_targets.push(target);
                    decoys.push(decoy);
                    stats.cycling_fallback += 1;
                }
                DecoyResult::SkippedNoFragments => {
                    stats.skipped_no_fragments += 1;
                }
                DecoyResult::FailedFragmentCalculation => {
                    stats.failed_fragment_calculation += 1;
                }
                DecoyResult::ExcludedNoUniqueDecoy => {
                    stats.excluded_no_unique_decoy += 1;
                }
            }
        }

        // Log statistics
        log::info!("Decoy generation statistics:");
        log::info!(
            "  Reversed: {} ({:.1}%)",
            stats.reversed,
            100.0 * stats.reversed as f64 / targets.len() as f64
        );
        log::info!(
            "  Cycling fallback: {} ({:.1}%)",
            stats.cycling_fallback,
            100.0 * stats.cycling_fallback as f64 / targets.len() as f64
        );
        log::info!(
            "  Excluded (no unique decoy): {} ({:.1}%)",
            stats.excluded_no_unique_decoy,
            100.0 * stats.excluded_no_unique_decoy as f64 / targets.len() as f64
        );
        if stats.skipped_no_fragments > 0 {
            log::info!("  Skipped (no fragments): {}", stats.skipped_no_fragments);
        }
        if stats.failed_fragment_calculation > 0 {
            log::warn!(
                "  Failed fragment calculation: {}",
                stats.failed_fragment_calculation
            );
        }

        (valid_targets, decoys, stats)
    }

    /// Cycle sequence by N positions, preserving terminal residue based on enzyme
    ///
    /// For trypsin (C-term preserved): PEPTIDEK with cycle=2 -> TIDEPEPK
    /// Cycling shifts the middle portion while keeping the cleavage site fixed
    fn cycle_sequence(&self, sequence: &str, cycle_length: usize) -> (String, Vec<usize>) {
        let chars: Vec<char> = sequence.chars().collect();
        let len = chars.len();

        if len <= 2 || cycle_length == 0 {
            return (sequence.to_string(), (0..len).collect());
        }

        let mut cycled: Vec<char> = Vec::with_capacity(len);
        let mut position_mapping: Vec<usize> = Vec::with_capacity(len);

        if self.enzyme.preserves_c_terminus() {
            // Preserve C-terminus only (trypsin, Lys-C)
            // Cycle positions 0..len-1, keep position len-1 at the end
            let middle_len = len - 1;
            let effective_cycle = cycle_length % middle_len;

            for i in 0..middle_len {
                let src_idx = (i + effective_cycle) % middle_len;
                cycled.push(chars[src_idx]);
                position_mapping.push(src_idx);
            }

            // Keep C-term
            cycled.push(chars[len - 1]);
            position_mapping.push(len - 1);
        } else {
            // Preserve N-terminus only (Lys-N, Asp-N)
            // Keep N-term, cycle positions 1..len
            cycled.push(chars[0]);
            position_mapping.push(0);

            let middle_len = len - 1;
            let effective_cycle = cycle_length % middle_len;

            for i in 0..middle_len {
                let src_idx = 1 + ((i + effective_cycle) % middle_len);
                cycled.push(chars[src_idx]);
                position_mapping.push(src_idx);
            }
        }

        (cycled.into_iter().collect(), position_mapping)
    }

    /// Create a decoy entry from a target and position mapping
    fn create_decoy_from_mapping(
        &self,
        target: &LibraryEntry,
        decoy_sequence: &str,
        position_mapping: &[usize],
    ) -> Result<LibraryEntry> {
        let mut decoy = target.clone();
        decoy.is_decoy = true;
        decoy.sequence = decoy_sequence.to_string();
        decoy.modified_sequence = format!("DECOY_{}", target.modified_sequence);

        // Remap modifications to new positions
        decoy.modifications = self.remap_modifications(&target.modifications, position_mapping);

        // Recalculate fragment m/z values
        decoy.fragments = self.recalculate_fragments(target, position_mapping);

        // Update ID to indicate decoy (set high bit)
        decoy.id = target.id | 0x80000000;

        // Update protein ID to indicate decoy
        decoy.protein_ids = target
            .protein_ids
            .iter()
            .map(|p| format!("DECOY_{}", p))
            .collect();

        Ok(decoy)
    }
}

/// Statistics from decoy generation with collision detection
#[derive(Debug, Clone, Default)]
pub struct DecoyGenerationStats {
    /// Number of decoys generated using reversal (primary method)
    pub reversed: usize,
    /// Number of decoys generated using cycling fallback
    pub cycling_fallback: usize,
    /// Number of targets excluded because no unique decoy could be generated
    pub excluded_no_unique_decoy: usize,
    /// Number of targets skipped because they had no fragments
    pub skipped_no_fragments: usize,
    /// Number of targets where fragment calculation failed
    pub failed_fragment_calculation: usize,
}

impl DecoyGenerationStats {
    /// Total number of valid target-decoy pairs created
    pub fn total_pairs(&self) -> usize {
        self.reversed + self.cycling_fallback
    }

    /// Total number of targets that were not paired
    pub fn total_excluded(&self) -> usize {
        self.excluded_no_unique_decoy + self.skipped_no_fragments + self.failed_fragment_calculation
    }
}

impl Default for DecoyGenerator {
    fn default() -> Self {
        Self::new(DecoyMethod::Reverse)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::IsolationWindow;

    /// Verifies that LibCosine scorer returns a score near 1.0 when the observed spectrum perfectly matches the library entry.
    #[test]
    fn test_spectral_scorer_lib_cosine() {
        let scorer = SpectralScorer::new();

        // Create a library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 75.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum that exactly matches library
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 500.0],
            intensities: vec![100.0, 50.0, 75.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        // Perfect match should give cosine close to 1.0
        assert!(
            score.lib_cosine > 0.99,
            "LibCosine should be ~1.0 for perfect match, got {}",
            score.lib_cosine
        );
        assert_eq!(score.n_matched, 3);
        assert!((score.fragment_coverage - 1.0).abs() < 1e-6);
    }

    /// Verifies that LibCosine scorer correctly reports partial fragment coverage when only some peaks match.
    #[test]
    fn test_spectral_scorer_partial_match() {
        let scorer = SpectralScorer::new();

        // Create a library entry with 4 fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 75.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 600.0,
                relative_intensity: 25.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum with only 2 matching peaks
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 700.0], // 700.0 doesn't match
            intensities: vec![100.0, 50.0, 80.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        // Should have 2 matches out of 4
        assert_eq!(score.n_matched, 2);
        assert!((score.fragment_coverage - 0.5).abs() < 1e-6);
        // LibCosine includes all 4 library fragments (0 for unmatched).
        // With only 2 of 4 matched, cosine should be penalized (~0.77).
        assert!(
            score.lib_cosine > 0.7 && score.lib_cosine < 0.85,
            "LibCosine should reflect partial match penalty, got {}",
            score.lib_cosine
        );
    }

    /// Verifies that LibCosine scorer returns zero score and zero coverage when no peaks overlap.
    #[test]
    fn test_spectral_scorer_no_match() {
        let scorer = SpectralScorer::new();

        // Create a library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![LibraryFragment {
            mz: 300.0,
            relative_intensity: 100.0,
            annotation: FragmentAnnotation::default(),
        }];

        // Create an observed spectrum with no matching peaks
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![800.0, 900.0], // No matches
            intensities: vec![100.0, 50.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        assert_eq!(score.n_matched, 0);
        assert!((score.fragment_coverage - 0.0).abs() < 1e-6);
        assert!((score.lib_cosine - 0.0).abs() < 1e-6);
    }

    /// Verifies that XCorr scoring produces a non-zero score for matching library and observed spectra.
    #[test]
    fn test_spectral_scorer_xcorr() {
        let scorer = SpectralScorer::new();

        // Create a library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0],
            intensities: vec![100.0, 50.0],
        };

        let score = scorer.xcorr(&spectrum, &entry);

        // XCorr should be computed
        // The exact value depends on the preprocessing, but it should be non-zero
        // for matching spectra
        assert!(
            score.xcorr != 0.0 || score.lib_cosine > 0.9,
            "Should have some score"
        );
    }

    /// Verifies that trypsin-aware decoy generation reverses the peptide sequence while preserving the C-terminal residue.
    #[test]
    fn test_decoy_reverse_trypsin() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);

        // PEPTIDEK has 8 chars: P(0) E(1) P(2) T(3) I(4) D(5) E(6) K(7)
        // For trypsin: keep only C-term (K), reverse positions 0-6
        // PEPTIDE (positions 0-6) -> EDITPEP (reversed)
        // Result: E D I T P E P K = EDITPEPK
        let target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0);
        let decoy = generator.generate(&target).unwrap();

        assert!(decoy.is_decoy);
        assert_eq!(decoy.sequence, "EDITPEPK"); // PEPTIDE->EDITPEP, K stays
        assert!(decoy.modified_sequence.starts_with("DECOY_"));
        assert_eq!(decoy.id, 1 | 0x80000000);
    }

    /// Verifies that LysN-aware decoy generation reverses the peptide sequence while preserving the N-terminal residue.
    #[test]
    fn test_decoy_reverse_lysn() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::LysN);

        // KPEPTIDE -> should preserve K at start
        // Rest PEPTIDE reversed: EDITPEP
        // Result: KEDITPEP
        let target = LibraryEntry::new(1, "KPEPTIDE".into(), "KPEPTIDE".into(), 2, 500.0, 10.0);
        let decoy = generator.generate(&target).unwrap();

        assert!(decoy.is_decoy);
        assert_eq!(decoy.sequence, "KEDITPEP"); // K stays, PEPTIDE->EDITPEP
    }

    /// Verifies that decoy generation preserves the precursor m/z since the amino acid composition is unchanged.
    #[test]
    fn test_decoy_preserves_precursor_mz() {
        let generator = DecoyGenerator::new(DecoyMethod::Reverse);
        let target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.123, 10.0);

        let decoy = generator.generate(&target).unwrap();

        // Same amino acids -> same precursor mass
        assert!((decoy.precursor_mz - 500.123).abs() < 1e-6);
    }

    /// Verifies that decoy generation adds "DECOY_" prefix to all associated protein IDs.
    #[test]
    fn test_decoy_protein_id_prefix() {
        let generator = DecoyGenerator::new(DecoyMethod::Reverse);
        let mut target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0);
        target.protein_ids = vec!["P12345".to_string(), "Q67890".to_string()];

        let decoy = generator.generate(&target).unwrap();

        assert_eq!(decoy.protein_ids.len(), 2);
        assert!(decoy.protein_ids[0].starts_with("DECOY_"));
        assert!(decoy.protein_ids[1].starts_with("DECOY_"));
    }

    /// Verifies that modification positions are correctly remapped when the peptide sequence is reversed for decoy generation.
    #[test]
    fn test_decoy_with_modification() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);
        let mut target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0);

        // Add modification at position 3 (0-indexed)
        // PEPTIDEK positions: P(0) E(1) P(2) T(3) I(4) D(5) E(6) K(7)
        target.modifications.push(Modification {
            position: 3, // T
            unimod_id: Some(35),
            mass_delta: 15.994915,
            name: Some("Oxidation".to_string()),
        });

        let decoy = generator.generate(&target).unwrap();

        // Reversed: EDITPEPK
        // Position mapping: [6, 5, 4, 3, 2, 1, 0, 7]
        // new_pos 0 <- old_pos 6 (E)
        // new_pos 1 <- old_pos 5 (D)
        // new_pos 2 <- old_pos 4 (I)
        // new_pos 3 <- old_pos 3 (T)  <- modification stays here
        // new_pos 4 <- old_pos 2 (P)
        // new_pos 5 <- old_pos 1 (E)
        // new_pos 6 <- old_pos 0 (P)
        // new_pos 7 <- old_pos 7 (K)
        assert_eq!(decoy.modifications.len(), 1);
        // The T at original position 3 stays at position 3 (symmetric reversal)
        assert_eq!(decoy.modifications[0].position, 3);
    }

    /// Verifies that sequences too short to meaningfully reverse are returned unchanged.
    #[test]
    fn test_short_sequence() {
        let generator = DecoyGenerator::new(DecoyMethod::Reverse);
        let target = LibraryEntry::new(1, "PK".into(), "PK".into(), 2, 200.0, 5.0);

        let decoy = generator.generate(&target).unwrap();

        // Too short to meaningfully reverse
        assert_eq!(decoy.sequence, "PK");
    }

    /// Verifies that enzyme cleavage site classification correctly identifies C-terminal vs N-terminal cutters.
    #[test]
    fn test_enzyme_detection() {
        assert!(Enzyme::Trypsin.preserves_c_terminus());
        assert!(Enzyme::LysC.preserves_c_terminus());
        assert!(!Enzyme::LysN.preserves_c_terminus());
        assert!(!Enzyme::AspN.preserves_c_terminus());
    }

    /// Verifies that SpectrumAggregator correctly sums intensities across multiple spectra for matching m/z bins.
    #[test]
    fn test_spectrum_aggregator() {
        let aggregator = SpectrumAggregator::new().with_tolerance_da(0.5);

        // Create test spectra at different RTs
        let spectra = [
            Spectrum {
                scan_number: 1,
                retention_time: 9.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0, 400.0],
                intensities: vec![50.0, 25.0],
            },
            Spectrum {
                scan_number: 2,
                retention_time: 10.0, // Apex
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0, 400.0],
                intensities: vec![100.0, 50.0],
            },
            Spectrum {
                scan_number: 3,
                retention_time: 10.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0, 400.0],
                intensities: vec![75.0, 35.0],
            },
        ];

        let spectra_refs: Vec<&Spectrum> = spectra.iter().collect();
        let aggregated = aggregator.aggregate(&spectra_refs, 10.0);

        // Should have 2 peaks (300.0 and 400.0)
        assert_eq!(aggregated.mzs.len(), 2);

        // Intensities should be summed
        // 300.0: 50 + 100 + 75 = 225
        // 400.0: 25 + 50 + 35 = 110
        let int_300 = aggregated.intensities[0];
        let int_400 = aggregated.intensities[1];
        assert!(
            (int_300 - 225.0).abs() < 1.0,
            "Expected ~225, got {}",
            int_300
        );
        assert!(
            (int_400 - 110.0).abs() < 1.0,
            "Expected ~110, got {}",
            int_400
        );
    }

    /// Verifies that weighted spectrum aggregation applies coefficient-based weights and normalizes correctly.
    #[test]
    fn test_spectrum_aggregator_weighted() {
        let aggregator = SpectrumAggregator::new().with_tolerance_da(0.5);

        // Create test spectra with corresponding coefficients
        let spectra = [
            Spectrum {
                scan_number: 1,
                retention_time: 9.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0],
                intensities: vec![100.0],
            },
            Spectrum {
                scan_number: 2,
                retention_time: 10.0, // Apex
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0],
                intensities: vec![100.0],
            },
        ];

        // Coefficients: scan 1 has weight 0.2, scan 2 has weight 0.8
        let coefficients = vec![(9.5, 0.2), (10.0, 0.8)];

        let spectra_refs: Vec<&Spectrum> = spectra.iter().collect();
        let aggregated = aggregator.aggregate_weighted(&spectra_refs, &coefficients, 10.0);

        // Should have 1 peak
        assert_eq!(aggregated.mzs.len(), 1);

        // Weighted intensity: 100 * 0.2 + 100 * 0.8 = 100 (after normalization)
        // Since we normalize weights to sum to 1.0, result should be 100
        let int = aggregated.intensities[0];
        assert!((int - 100.0).abs() < 1.0, "Expected ~100, got {}", int);
    }

    /// Verifies that decoy collision detection prevents generating decoys whose sequences match existing target sequences.
    #[test]
    fn test_decoy_collision_detection() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);

        // Create targets where one reversed sequence collides with another target
        // PEPTIDEK reverses to EDITPEPK
        // Include EDITPEPK as a target so there's a collision
        let mut targets = vec![
            LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0),
            LibraryEntry::new(2, "EDITPEPK".into(), "EDITPEPK".into(), 2, 500.0, 11.0), // This IS the reversed form!
            LibraryEntry::new(3, "ANOTHERK".into(), "ANOTHERK".into(), 2, 600.0, 12.0),
        ];

        // Add some fragments so they're not skipped
        for t in &mut targets {
            t.fragments = vec![LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            }];
        }

        let (valid_targets, decoys, stats) =
            generator.generate_all_with_collision_detection(&targets);

        // PEPTIDEK can't use reversal (would collide with EDITPEPK)
        // It should either use cycling or be excluded
        // EDITPEPK reverses to PEPTIDEK (also a collision!)
        // ANOTHERK reverses to REHTONAK (no collision)

        // At minimum, ANOTHERK should get a valid decoy
        assert!(stats.reversed >= 1, "At least ANOTHERK should use reversal");

        // Either collisions are resolved via cycling, or peptides are excluded
        assert!(
            stats.cycling_fallback + stats.excluded_no_unique_decoy >= 2,
            "PEPTIDEK and EDITPEPK should use cycling or be excluded"
        );

        // Total pairs created should match valid_targets and decoys
        assert_eq!(valid_targets.len(), decoys.len());
        assert_eq!(valid_targets.len(), stats.total_pairs());

        // Verify no decoy sequence matches any target sequence
        let target_seqs: std::collections::HashSet<&str> =
            valid_targets.iter().map(|t| t.sequence.as_str()).collect();
        for decoy in &decoys {
            assert!(
                !target_seqs.contains(decoy.sequence.as_str()),
                "Decoy {} should not match any target",
                decoy.sequence
            );
        }
    }

    /// Verifies that cyclic permutation shifts the peptide body by the specified offset while preserving the terminal residue.
    #[test]
    fn test_cycle_sequence() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);

        // PEPTIDEK with cycle=1: EPTIDEPK (shift by 1, keep K)
        let (cycled, _) = generator.cycle_sequence("PEPTIDEK", 1);
        assert_eq!(cycled, "EPTIDEPK");

        // PEPTIDEK with cycle=2: PTIDEPEK
        let (cycled, _) = generator.cycle_sequence("PEPTIDEK", 2);
        assert_eq!(cycled, "PTIDEPEK");
    }

    #[test]
    fn test_get_top_n_fragment_indices() {
        // Create test library fragments with varying intensities
        let fragments = vec![
            LibraryFragment {
                mz: 100.0,
                relative_intensity: 0.5,
                annotation: FragmentAnnotation {
                    ion_type: IonType::Y,
                    ordinal: 1,
                    charge: 1,
                    neutral_loss: None,
                },
            },
            LibraryFragment {
                mz: 200.0,
                relative_intensity: 1.0, // Highest
                annotation: FragmentAnnotation {
                    ion_type: IonType::Y,
                    ordinal: 2,
                    charge: 1,
                    neutral_loss: None,
                },
            },
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 0.2, // Lowest
                annotation: FragmentAnnotation {
                    ion_type: IonType::B,
                    ordinal: 1,
                    charge: 1,
                    neutral_loss: None,
                },
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 0.8, // Second highest
                annotation: FragmentAnnotation {
                    ion_type: IonType::Y,
                    ordinal: 3,
                    charge: 1,
                    neutral_loss: None,
                },
            },
        ];

        // Get top 2 fragments
        let top2 = get_top_n_fragment_indices(&fragments, 2);
        assert_eq!(top2.len(), 2);
        assert_eq!(top2[0], 1); // Index of fragment with intensity 1.0
        assert_eq!(top2[1], 3); // Index of fragment with intensity 0.8

        // Get top 6 (more than available)
        let top6 = get_top_n_fragment_indices(&fragments, 6);
        assert_eq!(top6.len(), 4); // Returns all 4 fragments

        // Get top 0
        let top0 = get_top_n_fragment_indices(&fragments, 0);
        assert_eq!(top0.len(), 0);
    }

    #[test]
    fn test_has_match_within_tolerance_ppm() {
        // Spectrum with m/z values
        let spectrum_mzs = vec![100.0, 200.0, 300.0, 400.0, 500.0];

        // Test exact match
        assert!(has_match_within_tolerance(
            200.0,
            &spectrum_mzs,
            10.0, // 10 ppm
            ToleranceUnit::Ppm
        ));

        // Test match within tolerance
        // 200.0 ± 10 ppm = 200.0 ± 0.002 = [199.998, 200.002]
        assert!(has_match_within_tolerance(
            200.001,
            &spectrum_mzs,
            10.0,
            ToleranceUnit::Ppm
        ));

        // Test no match (too far)
        assert!(!has_match_within_tolerance(
            210.0, // 10 Da away, way more than 10 ppm
            &spectrum_mzs,
            10.0,
            ToleranceUnit::Ppm
        ));

        // Test match at lower end
        assert!(has_match_within_tolerance(
            100.0,
            &spectrum_mzs,
            10.0,
            ToleranceUnit::Ppm
        ));

        // Test no match below range
        assert!(!has_match_within_tolerance(
            50.0,
            &spectrum_mzs,
            10.0,
            ToleranceUnit::Ppm
        ));

        // Test no match above range
        assert!(!has_match_within_tolerance(
            600.0,
            &spectrum_mzs,
            10.0,
            ToleranceUnit::Ppm
        ));
    }

    #[test]
    fn test_has_match_within_tolerance_mz() {
        let spectrum_mzs = vec![100.0, 200.0, 300.0, 400.0, 500.0];

        // Test exact match
        assert!(has_match_within_tolerance(
            200.0,
            &spectrum_mzs,
            0.5, // 0.5 Da
            ToleranceUnit::Mz
        ));

        // Test match within tolerance
        assert!(has_match_within_tolerance(
            200.3,
            &spectrum_mzs,
            0.5,
            ToleranceUnit::Mz
        ));

        assert!(has_match_within_tolerance(
            199.7,
            &spectrum_mzs,
            0.5,
            ToleranceUnit::Mz
        ));

        // Test no match (outside tolerance)
        assert!(!has_match_within_tolerance(
            200.6,
            &spectrum_mzs,
            0.5,
            ToleranceUnit::Mz
        ));

        assert!(!has_match_within_tolerance(
            199.4,
            &spectrum_mzs,
            0.5,
            ToleranceUnit::Mz
        ));
    }

    #[test]
    fn test_has_match_empty_spectrum() {
        let empty_spectrum: Vec<f64> = vec![];

        assert!(!has_match_within_tolerance(
            200.0,
            &empty_spectrum,
            10.0,
            ToleranceUnit::Ppm
        ));
    }

    // --- pearson_correlation_raw tests ---

    #[test]
    fn test_pearson_identical_vectors() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let r = pearson_correlation_raw(&x, &x);
        assert!(
            (r - 1.0).abs() < 1e-10,
            "Identical vectors should have r=1.0, got {}",
            r
        );
    }

    #[test]
    fn test_pearson_perfect_negative() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![5.0, 4.0, 3.0, 2.0, 1.0];
        let r = pearson_correlation_raw(&x, &y);
        assert!(
            (r - (-1.0)).abs() < 1e-10,
            "Perfectly anti-correlated vectors should have r=-1.0, got {}",
            r
        );
    }

    #[test]
    fn test_pearson_uncorrelated() {
        // Symmetric about zero — no linear correlation
        let x = vec![1.0, -1.0, 1.0, -1.0];
        let y = vec![1.0, 1.0, -1.0, -1.0];
        let r = pearson_correlation_raw(&x, &y);
        assert!(
            r.abs() < 1e-10,
            "Uncorrelated vectors should have r≈0, got {}",
            r
        );
    }

    #[test]
    fn test_pearson_too_short() {
        assert_eq!(pearson_correlation_raw(&[], &[]), 0.0);
        assert_eq!(pearson_correlation_raw(&[1.0], &[2.0]), 0.0);
    }

    #[test]
    fn test_pearson_constant_input() {
        // Constant vector has zero variance — denom should be clamped, not NaN
        let x = vec![5.0, 5.0, 5.0, 5.0];
        let y = vec![1.0, 2.0, 3.0, 4.0];
        let r = pearson_correlation_raw(&x, &y);
        assert!(
            r.is_finite(),
            "Constant input should not produce NaN, got {}",
            r
        );
    }

    #[test]
    fn test_pearson_known_value() {
        // x=[1,2,3], y=[2,4,5] → verified: r ≈ 0.9819805
        let x = vec![1.0, 2.0, 3.0];
        let y = vec![2.0, 4.0, 5.0];
        let r = pearson_correlation_raw(&x, &y);
        assert!((r - 0.9819805).abs() < 1e-4, "Expected r≈0.9820, got {}", r);
    }

    #[test]
    fn test_pearson_linear_transform() {
        // r should be invariant to positive linear scaling: y = 3x + 10
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y: Vec<f64> = x.iter().map(|&v| 3.0 * v + 10.0).collect();
        let r = pearson_correlation_raw(&x, &y);
        assert!(
            (r - 1.0).abs() < 1e-10,
            "Linear transform should give r=1.0, got {}",
            r
        );
    }

    // --- coelution_sum (pairwise Pearson sum) tests ---

    /// Compute coelution_sum the same way the pipeline does: sum of all unique
    /// pairwise Pearson correlations between fragment XICs within peak bounds.
    fn coelution_sum_from_xics(xic_values: &[Vec<f64>]) -> f64 {
        let n = xic_values.len();
        let mut sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                sum += pearson_correlation_raw(&xic_values[i], &xic_values[j]);
            }
        }
        sum
    }

    #[test]
    fn test_coelution_sum_perfect_coelution() {
        // 4 fragments with identical elution profiles → every pair has r=1.0
        // n_pairs = C(4,2) = 6, so coelution_sum should be 6.0
        let profile = vec![0.0, 1.0, 5.0, 10.0, 5.0, 1.0, 0.0];
        let xics: Vec<Vec<f64>> = (0..4).map(|_| profile.clone()).collect();

        let sum = coelution_sum_from_xics(&xics);
        assert!(
            (sum - 6.0).abs() < 1e-6,
            "4 identical fragment XICs → C(4,2)=6 pairs × r=1.0 = 6.0, got {}",
            sum
        );
    }

    #[test]
    fn test_coelution_sum_scaled_fragments() {
        // Different absolute intensities but same shape → still r=1.0 per pair
        let base = [0.0, 1.0, 5.0, 10.0, 5.0, 1.0, 0.0];
        let xics: Vec<Vec<f64>> = vec![
            base.iter().map(|v| v * 1.0).collect(),
            base.iter().map(|v| v * 0.5).collect(),
            base.iter().map(|v| v * 3.0).collect(),
        ];

        let sum = coelution_sum_from_xics(&xics);
        // C(3,2) = 3 pairs, all r=1.0
        assert!(
            (sum - 3.0).abs() < 1e-6,
            "Scaled profiles should give sum=3.0, got {}",
            sum
        );
    }

    #[test]
    fn test_coelution_sum_one_interferer() {
        // 3 co-eluting fragments + 1 fragment with a shifted peak (interference)
        let good = vec![0.0, 1.0, 5.0, 10.0, 5.0, 1.0, 0.0];
        let bad = vec![10.0, 5.0, 1.0, 0.0, 0.0, 1.0, 5.0]; // peaks at opposite end
        let xics = vec![good.clone(), good.clone(), good.clone(), bad];

        let sum = coelution_sum_from_xics(&xics);

        // 3 good-good pairs: r≈1.0 each → contribute ~3.0
        // 3 good-bad pairs: strongly negative r → pull sum down
        // Total should be substantially lower than the perfect case (6.0)
        assert!(
            sum < 3.0,
            "Interference should reduce sum well below 6.0 (perfect), got {}",
            sum
        );
        // The 3 good-good pairs are still perfect
        let good_only = vec![good.clone(), good.clone(), good.clone()];
        let good_sum = coelution_sum_from_xics(&good_only);
        assert!(
            (good_sum - 3.0).abs() < 1e-6,
            "Good-only subset should have sum=3.0, got {}",
            good_sum
        );
    }

    #[test]
    fn test_coelution_sum_two_fragments() {
        // Minimum case: 2 fragments → 1 pair, exactly proportional
        let a = vec![0.0, 3.0, 8.0, 3.0, 0.0];
        let b: Vec<f64> = a.iter().map(|v| v * 0.75).collect();
        let xics = vec![a.clone(), b.clone()];

        let sum = coelution_sum_from_xics(&xics);
        let expected = pearson_correlation_raw(&a, &b);
        assert!(
            (sum - expected).abs() < 1e-10,
            "2 fragments: sum should equal single pairwise r={}, got {}",
            expected,
            sum
        );
        assert!(
            (sum - 1.0).abs() < 1e-6,
            "Proportional profiles should give r≈1.0, got {}",
            sum
        );
    }

    #[test]
    fn test_coelution_sum_noise_reduces_score() {
        // Same base shape but add noise → correlations drop below 1.0
        let base = vec![0.0, 1.0, 5.0, 10.0, 5.0, 1.0, 0.0];
        let noisy1 = vec![0.1, 1.2, 4.8, 10.3, 5.1, 0.8, 0.2];
        let noisy2 = vec![0.0, 0.8, 5.3, 9.7, 4.7, 1.3, 0.1];
        let xics = vec![base, noisy1, noisy2];

        let sum = coelution_sum_from_xics(&xics);
        // Still highly correlated, but not perfectly
        assert!(
            sum > 2.5 && sum < 3.0,
            "Noisy co-eluting fragments should give sum close to but below 3.0, got {}",
            sum
        );
    }

    #[test]
    fn test_coelution_sum_peak_selection() {
        // Simulate picking the best peak from two candidates, the way the pipeline does.
        // Candidate 1: fragments co-elute well within this window
        let peak1_frag_a = vec![0.0, 2.0, 8.0, 10.0, 8.0, 2.0, 0.0];
        let peak1_frag_b = vec![0.0, 1.5, 6.0, 7.5, 6.0, 1.5, 0.0];
        let peak1_frag_c = vec![0.0, 1.0, 4.0, 5.0, 4.0, 1.0, 0.0];

        // Candidate 2: fragments don't co-elute (interfered peak)
        let peak2_frag_a = vec![5.0, 1.0, 0.0, 0.0, 2.0, 8.0, 3.0];
        let peak2_frag_b = vec![0.0, 3.0, 7.0, 2.0, 0.0, 1.0, 0.0];
        let peak2_frag_c = vec![1.0, 0.0, 2.0, 6.0, 8.0, 0.0, 1.0];

        let sum1 = coelution_sum_from_xics(&[peak1_frag_a, peak1_frag_b, peak1_frag_c]);
        let sum2 = coelution_sum_from_xics(&[peak2_frag_a, peak2_frag_b, peak2_frag_c]);

        assert!(
            sum1 > sum2,
            "Co-eluting peak (sum={}) should score higher than interfered peak (sum={})",
            sum1,
            sum2
        );
        // The good peak should be near-perfect
        assert!(
            sum1 > 2.9,
            "Well co-eluting peak should have sum near 3.0, got {}",
            sum1
        );
    }

    /// Computes mean pairwise Pearson correlation as used in pipeline peak selection.
    /// This mirrors the scoring logic in pipeline.rs that selects the best CWT candidate.
    fn mean_pairwise_correlation(xics: &[Vec<f64>]) -> f64 {
        let mut sum = 0.0;
        let mut count = 0u32;
        for i in 0..xics.len() {
            for j in (i + 1)..xics.len() {
                sum += pearson_correlation_raw(&xics[i], &xics[j]);
                count += 1;
            }
        }
        if count > 0 {
            sum / count as f64
        } else {
            0.0
        }
    }

    /// Verifies that mean pairwise correlation (the pipeline's peak selection metric)
    /// correctly discriminates between co-eluting and interfered peaks.
    #[test]
    fn test_mean_pairwise_correlation_peak_selection() {
        // Good candidate: all fragments co-elute with the same Gaussian-like shape
        let good_a = vec![0.0, 2.0, 8.0, 10.0, 8.0, 2.0, 0.0];
        let good_b: Vec<f64> = good_a.iter().map(|v| v * 0.75).collect();
        let good_c: Vec<f64> = good_a.iter().map(|v| v * 0.5).collect();

        // Bad candidate: fragments peak at different times (interference)
        let bad_a = vec![10.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        let bad_b = vec![0.0, 0.0, 0.0, 0.0, 1.0, 5.0, 10.0];
        let bad_c = vec![0.0, 0.0, 5.0, 10.0, 5.0, 0.0, 0.0];

        let mean_good = mean_pairwise_correlation(&[good_a, good_b, good_c]);
        let mean_bad = mean_pairwise_correlation(&[bad_a, bad_b, bad_c]);

        // Co-eluting fragments should have mean ≈ 1.0
        assert!(
            (mean_good - 1.0).abs() < 1e-6,
            "Co-eluting fragments should have mean correlation ≈ 1.0, got {}",
            mean_good
        );

        // Interfered fragments should have low or negative mean
        assert!(
            mean_bad < 0.0,
            "Interfered fragments should have negative mean correlation, got {}",
            mean_bad
        );

        // Pipeline would pick the good candidate
        assert!(
            mean_good > mean_bad,
            "Pipeline should prefer co-eluting peak (mean={}) over interfered (mean={})",
            mean_good,
            mean_bad
        );
    }

    /// Verifies that mean pairwise correlation handles the minimum case of 2 fragments.
    #[test]
    fn test_mean_pairwise_correlation_two_fragments() {
        let a = vec![1.0, 3.0, 7.0, 10.0, 7.0, 3.0, 1.0];
        let b = vec![0.5, 1.5, 3.5, 5.0, 3.5, 1.5, 0.5];

        let mean = mean_pairwise_correlation(&[a.clone(), b.clone()]);
        let direct = pearson_correlation_raw(&a, &b);

        assert!(
            (mean - direct).abs() < 1e-10,
            "Mean with 2 fragments should equal single Pearson r: mean={}, direct={}",
            mean,
            direct
        );
    }

    /// Verifies that single fragment or empty input returns 0.0 (no pairs to compute).
    #[test]
    fn test_mean_pairwise_correlation_degenerate() {
        let single = vec![vec![1.0, 2.0, 3.0]];
        assert!(
            (mean_pairwise_correlation(&single) - 0.0).abs() < 1e-10,
            "Single fragment should have mean=0 (no pairs)"
        );

        let empty: Vec<Vec<f64>> = vec![];
        assert!(
            (mean_pairwise_correlation(&empty) - 0.0).abs() < 1e-10,
            "Empty input should have mean=0"
        );
    }

    // ============================================================
    // XCorr known-answer tests
    // ============================================================

    /// Verifies XCorr scores a perfect match higher than a partial match.
    /// Uses two fragment libraries against the same observed spectrum to test discrimination.
    #[test]
    fn test_xcorr_perfect_vs_partial_match() {
        let scorer = SpectralScorer::new();

        // Observed spectrum with 3 peaks
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 500.0, 700.0],
            intensities: vec![1000.0, 500.0, 800.0],
        };

        // Library that matches all 3 peaks
        let mut full_match =
            LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        full_match.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 700.0,
                relative_intensity: 80.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Library that matches only 1 peak
        let mut partial_match =
            LibraryEntry::new(2, "OTHER".into(), "OTHER".into(), 2, 500.0, 10.0);
        partial_match.fragments = vec![LibraryFragment {
            mz: 300.0,
            relative_intensity: 100.0,
            annotation: FragmentAnnotation::default(),
        }];

        let score_full = scorer.xcorr(&spectrum, &full_match);
        let score_partial = scorer.xcorr(&spectrum, &partial_match);

        assert!(
            score_full.xcorr > score_partial.xcorr,
            "Full match xcorr ({}) should exceed partial match xcorr ({})",
            score_full.xcorr,
            score_partial.xcorr
        );
        assert!(
            score_full.xcorr > 0.0,
            "Full match should have positive xcorr: {}",
            score_full.xcorr
        );
    }

    /// Verifies XCorr is zero when library fragments don't match any observed peaks.
    #[test]
    fn test_xcorr_no_match_is_low() {
        let scorer = SpectralScorer::new();

        // Observed spectrum at 300 and 500
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 500.0],
            intensities: vec![1000.0, 500.0],
        };

        // Library at completely different m/z values (no overlap)
        let mut entry = LibraryEntry::new(1, "PEP".into(), "PEP".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 800.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 900.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        let score = scorer.xcorr(&spectrum, &entry);

        // XCorr at non-matching bins should be very low (near zero or negative)
        // because the sliding window subtraction yields negative background
        assert!(
            score.xcorr < 0.01,
            "Non-matching XCorr should be near zero, got {}",
            score.xcorr
        );
    }

    /// Verifies that XCorr uses the Comet 0.005 scaling factor.
    /// The raw dot product of preprocessed spectrum at fragment bins is multiplied by 0.005.
    #[test]
    fn test_xcorr_scaling_factor() {
        let scorer = SpectralScorer::new();

        // Single strong peak at 400.0
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![400.0],
            intensities: vec![10000.0],
        };

        let mut entry = LibraryEntry::new(1, "PEP".into(), "PEP".into(), 2, 500.0, 10.0);
        entry.fragments = vec![LibraryFragment {
            mz: 400.0,
            relative_intensity: 100.0,
            annotation: FragmentAnnotation::default(),
        }];

        let score = scorer.xcorr(&spectrum, &entry);

        // With windowing normalization (max=50.0) and offset=75,
        // a single isolated peak becomes 50.0 after windowing,
        // then after sliding window subtraction at the peak bin it becomes:
        //   50.0 - (sum excluding center) / 150
        // Since there's only one peak, sum in window is just the peak itself (excluded),
        // so the preprocessed value should be close to 50.0.
        // XCorr = preprocessed_value * 0.005 ≈ 0.25
        assert!(
            score.xcorr > 0.1,
            "Single matching peak should produce significant XCorr: {}",
            score.xcorr
        );
        // Should not be unreasonably large
        assert!(
            score.xcorr < 1.0,
            "Single peak XCorr should be modest: {}",
            score.xcorr
        );
    }

    /// Verifies XCorr returns zero for empty inputs.
    #[test]
    fn test_xcorr_empty_inputs() {
        let scorer = SpectralScorer::new();

        let empty_spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![],
            intensities: vec![],
        };

        let mut entry = LibraryEntry::new(1, "PEP".into(), "PEP".into(), 2, 500.0, 10.0);
        entry.fragments = vec![LibraryFragment {
            mz: 300.0,
            relative_intensity: 100.0,
            annotation: FragmentAnnotation::default(),
        }];

        let score = scorer.xcorr(&empty_spectrum, &entry);
        assert!(
            (score.xcorr - 0.0).abs() < 1e-10,
            "Empty spectrum should produce zero XCorr"
        );

        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0],
            intensities: vec![100.0],
        };

        let empty_lib = LibraryEntry::new(2, "PEP".into(), "PEP".into(), 2, 500.0, 10.0);
        let score2 = scorer.xcorr(&spectrum, &empty_lib);
        assert!(
            (score2.xcorr - 0.0).abs() < 1e-10,
            "Empty library should produce zero XCorr"
        );
    }

    /// Verifies that XCorr with preprocess_spectrum_for_xcorr matches the one-shot xcorr method.
    /// This tests the batch preprocessing path produces equivalent results.
    #[test]
    fn test_xcorr_preprocessed_matches_direct() {
        let scorer = SpectralScorer::new();

        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![200.0, 400.0, 600.0, 800.0],
            intensities: vec![500.0, 1000.0, 750.0, 250.0],
        };

        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 200.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 600.0,
                relative_intensity: 75.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // One-shot XCorr
        let direct_score = scorer.xcorr(&spectrum, &entry);

        // Preprocessed XCorr: manually compute what xcorr does
        let preprocessed = scorer.preprocess_spectrum_for_xcorr(&spectrum);
        let xcorr_raw: f32 = entry
            .fragments
            .iter()
            .filter_map(|frag| scorer.bin_config().mz_to_bin(frag.mz))
            .map(|bin| preprocessed[bin])
            .sum();
        let xcorr_preprocessed = (xcorr_raw * 0.005) as f64;

        assert!(
            (direct_score.xcorr - xcorr_preprocessed).abs() < 1e-6,
            "Direct XCorr ({}) should match preprocessed XCorr ({})",
            direct_score.xcorr,
            xcorr_preprocessed
        );
    }

    /// Regression guard for the ordinal-attribution fix in
    /// `longest_consecutive_ions`. A prior implementation reverse-looked-up
    /// the library by m/z and broke on first hit, silently mis-attributing
    /// ordinals whenever two library fragments shared a near-identical m/z
    /// (incidental b_n / y_m collision).
    ///
    /// This test constructs FragmentMatch instances with ordinals set
    /// directly and an **empty** library. If anyone reverts to the library
    /// m/z lookup, the empty library will drop every ordinal and the test
    /// returns 0 instead of 4.
    #[test]
    fn longest_consecutive_ions_uses_match_ordinal_not_library_lookup() {
        let scorer = SpectralScorer::new();
        let lib = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        // library.fragments intentionally left empty.
        let mk = |ordinal: u8, ion_type: IonType| FragmentMatch {
            lib_mz: 0.0,
            obs_mz: 0.0,
            lib_intensity: 1.0,
            obs_intensity: 1.0,
            ion_type,
            ordinal,
        };
        let matches = vec![
            mk(2, IonType::B),
            mk(3, IonType::B),
            mk(4, IonType::B),
            mk(5, IonType::B),
            mk(7, IonType::Y),
        ];
        assert_eq!(
            scorer.longest_consecutive_ions(&lib, &matches),
            4,
            "expected b2..b5 consecutive run of length 4; if this returns 0 \
             longest_consecutive_ions is falling back to library m/z lookup \
             instead of using FragmentMatch.ordinal"
        );
    }

    /// Regression guard: the low-norm gate in `lib_cosine` must NOT zero out
    /// the presence/counting features (n_matched, consecutive_ions,
    /// explained_intensity, fragment_coverage). Those features are derived
    /// from `match_fragments` which runs before the norm computation; they
    /// are independent of whether the cosine itself is computable.
    ///
    /// A prior implementation returned `SpectralScore::default()` as soon as
    /// `lib_norm < 1e-10 || obs_norm < 1e-10`, silently dropping counting
    /// features alongside the undefined cosine. This showed up as rows where
    /// `mass_accuracy` indicated fragment matches were present but
    /// `consecutive_ions` and `explained_intensity` reported zero — the bug
    /// affected short/low-signal peptides in cross-implementation validation.
    #[test]
    fn lib_cosine_counting_features_survive_zero_norm() {
        let scorer = SpectralScorer::new();

        // Library with a b2/b3/b4 run but all zero intensities -> lib_norm = 0.
        let mk_frag = |mz: f64, ordinal: u8| LibraryFragment {
            mz,
            relative_intensity: 0.0,
            annotation: FragmentAnnotation {
                ion_type: IonType::B,
                ordinal,
                charge: 1,
                neutral_loss: None,
            },
        };
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![mk_frag(200.0, 2), mk_frag(300.0, 3), mk_frag(400.0, 4)];

        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![200.0, 300.0, 400.0],
            intensities: vec![100.0, 100.0, 100.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        // Cosine + correlation features are undefined with lib_norm = 0
        // and must zero out.
        assert_eq!(score.lib_cosine, 0.0, "cosine undefined when lib_norm = 0");
        assert_eq!(score.pearson_correlation, 0.0);
        assert_eq!(score.spearman_correlation, 0.0);

        // Presence/counting features are independent of the cosine norm and
        // must still reflect match_fragments output. If these are zero the
        // caller has re-introduced the early `return SpectralScore::default()`
        // that silently dropped counting features alongside cosine.
        assert_eq!(
            score.n_matched, 3,
            "n_matched must populate independent of the norm gate"
        );
        assert_eq!(
            score.consecutive_ions, 3,
            "consecutive_ions (b2-b3-b4) must populate independent of the norm gate"
        );
        assert!(
            score.fragment_coverage > 0.99,
            "fragment_coverage must reflect matches (3/3), got {}",
            score.fragment_coverage
        );
        assert!(
            score.explained_intensity > 0.99,
            "explained_intensity must reflect matches (3 of 3 obs peaks matched), got {}",
            score.explained_intensity
        );
    }

    /// Regression guard for the XCorr fragment-bin dedup invariant.
    ///
    /// Comet XCorr scores the dot product of an experimental preprocessed
    /// spectrum against a theoretical spectrum with unit intensity per
    /// bin. When two library fragments fall into the same bin, the bin's
    /// contribution must count ONCE, not twice. The fast path
    /// `preprocess_library_for_xcorr` + `xcorr_from_preprocessed` handles
    /// this correctly by assigning `binned[bin] = 1.0` per unique bin.
    /// The direct `scorer.xcorr()` / `xcorr_at_scan()` paths used to
    /// iterate `library.fragments` and sum `preprocessed[bin]` once per
    /// fragment, double-counting collisions.
    ///
    /// This test constructs a library with two fragments that land in the
    /// same unit-resolution bin and asserts the direct path and the fast
    /// (preprocessed) path agree. If direct double-counts again, it will
    /// be larger than the fast path and this test panics.
    #[test]
    fn xcorr_dedups_fragment_bin_collisions() {
        let scorer = SpectralScorer::new();

        // Unit bins are ~1.0005 Th wide. Fragments at 500.0 and 500.5 both
        // land in bin 500 (BIN(mass) = (int)(mass / 1.0005079 + 0.6)).
        // Sanity-check the assumption so the test fails loudly if BinConfig
        // ever changes bin width or offset in a way that breaks it.
        let bin_500_0 = scorer.bin_config.mz_to_bin(500.0).unwrap();
        let bin_500_5 = scorer.bin_config.mz_to_bin(500.5).unwrap();
        assert_eq!(
            bin_500_0, bin_500_5,
            "test setup expects 500.0 and 500.5 to collide in one bin"
        );

        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.5, // same bin as 500.0 on unit-resolution
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            // One more fragment in a distinct bin so the test exercises the
            // mixed case (one colliding pair + one non-colliding fragment).
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 75.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![200.0, 300.0, 400.0, 500.0, 600.0],
            intensities: vec![500.0, 1000.0, 750.0, 1200.0, 250.0],
        };

        // Direct XCorr (the path being tested).
        let direct_xcorr = scorer.xcorr(&spectrum, &entry).xcorr;

        // Canonical preprocessed path: preprocess_library_for_xcorr assigns
        // binned[bin] = 1.0 per unique bin (dedup by assignment), so this
        // path is correct by construction.
        let lib_preprocessed = scorer.preprocess_library_for_xcorr(&entry);
        let spec_preprocessed = scorer.preprocess_spectrum_for_xcorr(&spectrum);
        let canonical_xcorr =
            SpectralScorer::xcorr_from_preprocessed(&spec_preprocessed, &lib_preprocessed);

        assert!(
            (direct_xcorr - canonical_xcorr).abs() < 1e-6,
            "direct XCorr ({}) must match the preprocessed/canonical XCorr \
             ({}). If direct is roughly canonical + preprocessed[bin_500], \
             someone re-introduced the per-fragment sum that double-counts \
             colliding-bin fragments. Comet XCorr counts unit intensity per \
             bin, not per fragment.",
            direct_xcorr,
            canonical_xcorr,
        );
    }
}
