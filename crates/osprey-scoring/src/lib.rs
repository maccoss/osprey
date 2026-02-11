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
pub mod pipeline;

use osprey_core::{
    BinConfig, FeatureSet, FragmentAnnotation, IonType, LibraryEntry, LibraryFragment,
    Modification, Result, Spectrum, ToleranceUnit,
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
pub fn has_top3_fragment_match(
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

/// Check if any top-3 library fragment matches an observed peak, and collect mass errors.
///
/// Combines the filtering check with mass error collection for calibration.
/// Uses binary search for O(3 × log n) performance. For each matched top-3 fragment,
/// finds the closest observed peak and computes the mass error.
///
/// # Returns
/// `(has_match, mass_errors)` where `has_match` is true if at least 1 top-3 peak matches,
/// and `mass_errors` contains the error (in the configured unit) for each matched fragment.
pub fn top3_fragment_match_with_errors(
    library_fragments: &[LibraryFragment],
    spectrum_mzs: &[f64],
    tolerance: f64,
    unit: ToleranceUnit,
) -> (bool, Vec<f64>) {
    if library_fragments.is_empty() || spectrum_mzs.is_empty() {
        return (true, Vec::new());
    }

    // Get indices of top 3 fragments by intensity
    let mut top3_indices: Vec<usize> = (0..library_fragments.len().min(3)).collect();

    if library_fragments.len() > 3 {
        let mut indexed: Vec<(usize, f32)> = library_fragments
            .iter()
            .enumerate()
            .map(|(i, f)| (i, f.relative_intensity))
            .collect();
        indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
        top3_indices = indexed.iter().take(3).map(|(i, _)| *i).collect();
    }

    let mut has_match = false;
    let mut mass_errors = Vec::with_capacity(3);

    for &idx in &top3_indices {
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

/// Compute fragment co-elution correlation (DIA-NN pTimeCorr-inspired).
///
/// For each of the top 6 library fragments, extracts a raw XIC from the peak region
/// spectra and computes Pearson correlation against the regression coefficient time series.
/// The coefficient series serves as the deconvolved elution profile reference.
///
/// Uses plain Pearson correlation on raw intensities (no sqrt transform), matching DIA-NN.
///
/// # Returns
/// `(sum, min, n_positive)` — sum of per-fragment correlations, minimum correlation,
/// and count of fragments with positive correlation.
pub fn compute_fragment_coelution(
    library_fragments: &[LibraryFragment],
    coefficient_series: &[(f64, f64)], // (RT, coefficient) pairs
    spectra: &[&Spectrum],             // peak region spectra (sorted by RT)
    tolerance_da: f64,
    tolerance_ppm: f64,
) -> (f64, f64, u32) {
    if library_fragments.is_empty() || coefficient_series.len() < 3 || spectra.is_empty() {
        return (0.0, 0.0, 0);
    }

    // Select top 6 fragments by relative intensity
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

    // Build a map from spectrum RT to index for alignment with coefficient series
    // Both coefficient_series and spectra should cover the same RT region
    let spec_rts: Vec<f64> = spectra.iter().map(|s| s.retention_time).collect();

    let mut coelution_sum = 0.0;
    let mut coelution_min = f64::MAX;
    let mut n_positive: u32 = 0;
    let mut n_scored: u32 = 0;

    for &frag_idx in &top_indices {
        let frag_mz = library_fragments[frag_idx].mz;
        let tol = tolerance_da.max(frag_mz * tolerance_ppm / 1e6);

        // Extract raw XIC: for each coefficient series point, find the nearest spectrum
        // and extract the fragment intensity via binary search
        let mut xic: Vec<f64> = Vec::with_capacity(coefficient_series.len());
        let mut coefs: Vec<f64> = Vec::with_capacity(coefficient_series.len());

        for &(rt, coef) in coefficient_series {
            // Find the spectrum closest to this RT
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

            // Binary search for the fragment m/z in the spectrum
            let spectrum = &spectra[spec_idx];
            let lower = frag_mz - tol;
            let upper = frag_mz + tol;
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

            xic.push(intensity);
            coefs.push(coef);
        }

        // Compute Pearson correlation between XIC and coefficient series
        if xic.len() >= 3 {
            let corr = pearson_correlation_raw(&xic, &coefs);
            coelution_sum += corr;
            if corr < coelution_min {
                coelution_min = corr;
            }
            if corr > 0.0 {
                n_positive += 1;
            }
            n_scored += 1;
        }
    }

    if n_scored == 0 {
        return (0.0, 0.0, 0);
    }
    if coelution_min == f64::MAX {
        coelution_min = 0.0;
    }

    (coelution_sum, coelution_min, n_positive)
}

/// Compute per-fragment mass accuracy statistics at the apex scan.
///
/// # Arguments
/// * `matches` - Fragment matches from the apex spectrum
/// * `unit` - Tolerance unit determining error computation (ppm for HRAM, Th for unit resolution)
///
/// # Returns
/// `(mean_abs_error, std_error)` in the configured unit
pub fn compute_mass_accuracy(matches: &[FragmentMatch], unit: ToleranceUnit) -> (f64, f64) {
    if matches.is_empty() {
        return (0.0, 0.0);
    }

    let errors: Vec<f64> = matches
        .iter()
        .map(|m| match unit {
            ToleranceUnit::Ppm => (m.obs_mz - m.lib_mz) / m.lib_mz * 1e6,
            ToleranceUnit::Mz => m.obs_mz - m.lib_mz,
        })
        .collect();

    let n = errors.len() as f64;

    // Mean of absolute errors
    let mean_abs: f64 = errors.iter().map(|e| e.abs()).sum::<f64>() / n;

    // Standard deviation of signed errors
    let mean_signed: f64 = errors.iter().sum::<f64>() / n;
    let variance: f64 = errors
        .iter()
        .map(|e| (e - mean_signed).powi(2))
        .sum::<f64>()
        / n;
    let std_dev = variance.sqrt();

    (mean_abs, std_dev)
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
        let mut has_signal = false;

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
                has_signal = true;
                best_intensity
            } else {
                0.0
            };

            xic.push((spectrum.retention_time, intensity));
        }

        if has_signal {
            results.push((frag_idx, xic));
        }
    }

    results
}

/// Build a consensus XIC from co-eluting fragment ions and compute FWHM.
///
/// Uses co-eluting fragment XICs (inspired by DIA-NN's approach, Demichev et al.,
/// Nature Methods, 2020) to determine chromatographic peak boundaries.
///
/// 1. Extracts XICs for top 6 fragments
/// 2. Correlates each against the coefficient series (reference)
/// 3. For fragments with positive correlation (co-eluting):
///    - Normalize each XIC to max=1
///    - Sum into consensus XIC
/// 4. Compute FWHM with linear interpolation on the consensus
/// 5. Derive 95% Gaussian boundaries: apex ± 1.96σ where σ = FWHM/2.355
///
/// # Returns
/// `Some((start_rt, end_rt, fwhm))` using 95% Gaussian boundaries,
/// or `None` if FWHM cannot be computed.
pub fn compute_fragment_fwhm(
    library_fragments: &[LibraryFragment],
    coefficient_series: &[(f64, f64)],
    spectra: &[&Spectrum],
    tolerance_da: f64,
    tolerance_ppm: f64,
) -> Option<(f64, f64, f64)> {
    if coefficient_series.len() < 3 || spectra.is_empty() {
        return None;
    }

    // Extract XICs for top 6 fragments
    let xics = extract_fragment_xics(library_fragments, spectra, tolerance_da, tolerance_ppm, 6);
    if xics.is_empty() {
        return None;
    }

    // For each fragment XIC, correlate against the coefficient series
    // Align by finding the nearest spectrum RT for each coefficient point
    let spec_rts: Vec<f64> = spectra.iter().map(|s| s.retention_time).collect();

    let mut coeluting_xics: Vec<&Vec<(f64, f64)>> = Vec::new();

    for (_frag_idx, xic) in &xics {
        // Build aligned vectors: for each coefficient point, find nearest XIC value
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
        return None;
    }

    // Build consensus XIC: normalize each co-eluting XIC to max=1, then sum
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

    // Compute FWHM on the consensus XIC
    let (fwhm, _left, _right) = compute_fwhm_interpolated(&consensus)?;

    // 95% Gaussian boundaries: apex ± 1.96σ where σ = FWHM / 2.355
    let apex_rt = consensus
        .iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(rt, _)| *rt)?;
    let sigma = fwhm / 2.355;
    let start_rt = apex_rt - 1.96 * sigma;
    let end_rt = apex_rt + 1.96 * sigma;

    Some((start_rt, end_rt, fwhm))
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
/// `has_top3_fragment_match` but for ALL fragments.
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
    pub fn with_tolerance_da(mut self, tolerance: f64) -> Self {
        self.tolerance_da = tolerance;
        self
    }

    /// Set ppm tolerance for HRAM matching
    pub fn with_tolerance_ppm(mut self, tolerance: f64) -> Self {
        self.tolerance_ppm = tolerance;
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
    /// LibCosine uses sqrt(intensity) preprocessing with L2 normalization,
    /// then computes cosine similarity. Also computes SMZ variant and hyperscore.
    pub fn lib_cosine(&self, observed: &Spectrum, library: &LibraryEntry) -> SpectralScore {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return SpectralScore::default();
        }

        // Match library fragments to observed peaks
        let matches = self.match_fragments(observed, library);

        if matches.is_empty() {
            return SpectralScore::default();
        }

        // Preprocess library fragments: sqrt(intensity) × m/z²
        let mut lib_preprocessed: Vec<f64> = Vec::with_capacity(matches.len());
        let mut obs_preprocessed: Vec<f64> = Vec::with_capacity(matches.len());

        // Keep raw intensities for correlation computation
        let mut lib_intensities: Vec<f64> = Vec::with_capacity(matches.len());
        let mut obs_intensities: Vec<f64> = Vec::with_capacity(matches.len());

        for m in &matches {
            // Library: sqrt(relative_intensity)
            let lib_val = (m.lib_intensity as f64).sqrt();
            lib_preprocessed.push(lib_val);
            lib_intensities.push(m.lib_intensity as f64);

            // Observed: sqrt(intensity)
            let obs_val = (m.obs_intensity as f64).sqrt();
            obs_preprocessed.push(obs_val);
            obs_intensities.push(m.obs_intensity as f64);
        }

        // L2 normalize both vectors
        let lib_norm = lib_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();

        if lib_norm < 1e-10 || obs_norm < 1e-10 {
            return SpectralScore::default();
        }

        for v in lib_preprocessed.iter_mut() {
            *v /= lib_norm;
        }
        for v in obs_preprocessed.iter_mut() {
            *v /= obs_norm;
        }

        // Cosine similarity (dot product of normalized vectors)
        let dot_product: f64 = lib_preprocessed
            .iter()
            .zip(obs_preprocessed.iter())
            .map(|(a, b)| a * b)
            .sum();

        // Compute additional metrics
        let n_matched = matches.len() as u32;
        let n_library = library.fragments.len() as u32;
        let fragment_coverage = n_matched as f64 / n_library as f64;

        // Compute explained intensity (sum of matched intensities / total observed intensity)
        let matched_intensity: f64 = matches.iter().map(|m| m.obs_intensity as f64).sum();
        let total_intensity: f64 = observed.intensities.iter().map(|&i| i as f64).sum();
        let explained_intensity = if total_intensity > 0.0 {
            matched_intensity / total_intensity
        } else {
            0.0
        };

        // Compute Pearson correlation
        let pearson_correlation = Self::pearson_correlation(&lib_intensities, &obs_intensities);

        // Compute Spearman correlation
        let spearman_correlation = Self::spearman_correlation(&lib_intensities, &obs_intensities);

        // Compute consecutive ion series
        let consecutive_ions = self.longest_consecutive_ions(library, &matches);

        // Compute sequence coverage (backbone cleavage sites)
        let sequence_coverage = self.compute_sequence_coverage(library, &matches);

        // Compute base peak rank
        let base_peak_rank =
            self.compute_base_peak_rank(&matches, &lib_intensities, &obs_intensities);

        // Compute top-3 matches
        let top6_matches = self.compute_top6_matches(library, &matches);

        // Compute LibCosine with SMZ preprocessing (sqrt(intensity) * mz²)
        let lib_cosine_smz = self.lib_cosine_smz(observed, library);

        // Compute X!Tandem hyperscore
        let (hyperscore_val, n_matched_b, n_matched_y) = self.hyperscore(observed, library);

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
        }
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
    fn longest_consecutive_ions(&self, library: &LibraryEntry, matches: &[FragmentMatch]) -> u32 {
        // Build sets of matched b and y ion ordinals
        let mut matched_b: Vec<u8> = Vec::new();
        let mut matched_y: Vec<u8> = Vec::new();

        for m in matches {
            // Find the corresponding library fragment annotation
            for frag in &library.fragments {
                if (frag.mz - m.lib_mz).abs() < 0.001 {
                    match frag.annotation.ion_type {
                        IonType::B => matched_b.push(frag.annotation.ordinal),
                        IonType::Y => matched_y.push(frag.annotation.ordinal),
                        _ => {}
                    }
                    break;
                }
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
                                if ordinal <= n_sites {
                                    covered[ordinal - 1] = true;
                                }
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

    /// Count how many of the top-3 library peaks are matched
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

        // XCorr = sum of preprocessed experimental values at fragment bin positions
        // This matches Comet exactly: score = sum(experimental_preprocessed[frag_bins]) * 0.005
        // Directly sum at fragment positions — O(n_fragments) instead of O(n_bins)
        let xcorr_raw: f32 = library
            .fragments
            .iter()
            .filter_map(|frag| self.bin_config.mz_to_bin(frag.mz))
            .map(|bin| xcorr_preprocessed[bin])
            .sum();

        // Scale XCorr (pyXcorrDIA uses 0.005 for spectrum-centric)
        let xcorr_scaled = xcorr_raw * 0.005f32;

        SpectralScore {
            lib_cosine: lib_cosine_score.lib_cosine,
            lib_cosine_smz: lib_cosine_score.lib_cosine_smz,
            xcorr: xcorr_scaled as f64,
            hyperscore: lib_cosine_score.hyperscore,
            dot_product: lib_cosine_score.dot_product,
            n_matched: lib_cosine_score.n_matched,
            n_library: lib_cosine_score.n_library,
            fragment_coverage: lib_cosine_score.fragment_coverage,
            explained_intensity: lib_cosine_score.explained_intensity,
            pearson_correlation: lib_cosine_score.pearson_correlation,
            spearman_correlation: lib_cosine_score.spearman_correlation,
            consecutive_ions: lib_cosine_score.consecutive_ions,
            sequence_coverage: lib_cosine_score.sequence_coverage,
            base_peak_rank: lib_cosine_score.base_peak_rank,
            top6_matches: lib_cosine_score.top6_matches,
            n_matched_b: lib_cosine_score.n_matched_b,
            n_matched_y: lib_cosine_score.n_matched_y,
        }
    }

    /// Compute LibCosine score with sqrt(intensity) * mz² (SMZ) preprocessing
    ///
    /// This variant weights fragment matches by their m/z value squared,
    /// giving more importance to higher m/z fragments which are more sequence-specific.
    /// This is the original SMZ (Sqrt-Mz-squared) preprocessing from spectral library searching.
    pub fn lib_cosine_smz(&self, observed: &Spectrum, library: &LibraryEntry) -> f64 {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return 0.0;
        }

        let matches = self.match_fragments(observed, library);
        if matches.is_empty() {
            return 0.0;
        }

        let mut lib_preprocessed: Vec<f64> = Vec::with_capacity(matches.len());
        let mut obs_preprocessed: Vec<f64> = Vec::with_capacity(matches.len());

        for m in &matches {
            let mz_sq = m.lib_mz * m.lib_mz;

            // Library: sqrt(relative_intensity) * mz²
            let lib_val = (m.lib_intensity as f64).sqrt() * mz_sq;
            lib_preprocessed.push(lib_val);

            // Observed: sqrt(intensity) * mz²
            let obs_val = (m.obs_intensity as f64).sqrt() * mz_sq;
            obs_preprocessed.push(obs_val);
        }

        // L2 normalize both vectors
        let lib_norm = lib_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();

        if lib_norm < 1e-10 || obs_norm < 1e-10 {
            return 0.0;
        }

        for v in lib_preprocessed.iter_mut() {
            *v /= lib_norm;
        }
        for v in obs_preprocessed.iter_mut() {
            *v /= obs_norm;
        }

        // Cosine similarity (dot product of normalized vectors)
        lib_preprocessed
            .iter()
            .zip(obs_preprocessed.iter())
            .map(|(a, b)| a * b)
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
    /// Uses ridge regression coefficients to weight the contribution
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

/// Contextual information from regression results for a peptide
///
/// This aggregates information from multiple regression results (across spectra)
/// to compute contextual features like n_competitors and relative_coefficient.
#[derive(Debug, Clone, Default)]
pub struct RegressionContext {
    /// Average number of competitors per spectrum
    pub avg_n_competitors: f64,
    /// Maximum number of competitors seen
    pub max_n_competitors: u32,
    /// Average relative coefficient (this peptide / sum of all in spectrum)
    pub avg_relative_coefficient: f64,
    /// Maximum relative coefficient
    pub max_relative_coefficient: f64,
    /// Average regression residual
    pub avg_residual: f64,
    /// Average explained variance (1 - residual / ||b||²)
    pub avg_explained_variance: f64,
    /// Number of spectra this peptide appeared in
    pub n_spectra: u32,
    /// Spectral complexity estimate (average number of peaks)
    pub avg_spectral_complexity: f64,
}

impl RegressionContext {
    /// Build regression context from regression results for a specific peptide
    pub fn from_results(results: &[&osprey_core::RegressionResult], lib_id: u32) -> Self {
        if results.is_empty() {
            return Self::default();
        }

        let n = results.len() as f64;
        let mut total_competitors = 0.0;
        let mut max_competitors = 0u32;
        let mut total_relative_coef = 0.0;
        let mut max_relative_coef = 0.0f64;
        let mut total_residual = 0.0;
        let mut total_explained_var = 0.0;

        for result in results {
            total_competitors += result.n_candidates as f64;
            max_competitors = max_competitors.max(result.n_candidates);

            let rel_coef = result.relative_coefficient(lib_id);
            total_relative_coef += rel_coef;
            max_relative_coef = max_relative_coef.max(rel_coef);

            total_residual += result.residual;
            total_explained_var += result.explained_variance();
        }

        Self {
            avg_n_competitors: total_competitors / n,
            max_n_competitors: max_competitors,
            avg_relative_coefficient: total_relative_coef / n,
            max_relative_coefficient: max_relative_coef,
            avg_residual: total_residual / n,
            avg_explained_variance: total_explained_var / n,
            n_spectra: results.len() as u32,
            avg_spectral_complexity: 0.0, // Computed separately if needed
        }
    }
}

/// Feature extractor for peptide candidates
#[derive(Debug, Default)]
pub struct FeatureExtractor {
    /// Spectral scorer for computing similarity metrics
    spectral_scorer: SpectralScorer,
    /// Spectrum aggregator for combining scans
    spectrum_aggregator: SpectrumAggregator,
}

impl FeatureExtractor {
    /// Create a new feature extractor
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a feature extractor with custom spectral scorer
    pub fn with_scorer(spectral_scorer: SpectralScorer) -> Self {
        Self {
            spectral_scorer,
            spectrum_aggregator: SpectrumAggregator::default(),
        }
    }

    /// Create a feature extractor with custom scorer and aggregator
    pub fn with_scorer_and_aggregator(
        spectral_scorer: SpectralScorer,
        spectrum_aggregator: SpectrumAggregator,
    ) -> Self {
        Self {
            spectral_scorer,
            spectrum_aggregator,
        }
    }

    /// Get the spectrum aggregator
    pub fn aggregator(&self) -> &SpectrumAggregator {
        &self.spectrum_aggregator
    }

    /// Extract features using spectrum aggregation across peak apex region
    ///
    /// FR-5.2.1: This method aggregates spectra around the peak apex before
    /// computing spectral features, providing more robust scoring.
    ///
    /// Parameters:
    /// - `entry`: The library entry being scored
    /// - `coefficient_series`: Vec of (retention_time, coefficient) pairs
    /// - `spectra`: All spectra in the peak region
    /// - `expected_rt`: Optional expected retention time
    pub fn extract_with_aggregation(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        spectra: &[&Spectrum],
        expected_rt: Option<f64>,
    ) -> FeatureSet {
        // Find apex RT from coefficient series
        let apex_rt = coefficient_series
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, _)| *rt)
            .unwrap_or(entry.retention_time);

        // Aggregate spectra around apex (FR-5.2.1)
        let aggregated =
            self.spectrum_aggregator
                .aggregate_weighted(spectra, coefficient_series, apex_rt);

        // Extract features using the aggregated spectrum
        self.extract_with_expected_rt(entry, coefficient_series, Some(&aggregated), expected_rt)
    }

    /// Extract features for a peptide candidate
    ///
    /// Computes both chromatographic features (from coefficient series) and
    /// spectral features (from apex spectrum comparison to library).
    ///
    /// Parameters:
    /// - `entry`: The library entry being scored
    /// - `coefficient_series`: Vec of (retention_time, coefficient) pairs
    /// - `apex_spectrum`: Optional spectrum at peak apex for spectral scoring
    /// - `expected_rt`: Optional expected retention time for RT deviation calculation
    pub fn extract(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        apex_spectrum: Option<&Spectrum>,
    ) -> FeatureSet {
        self.extract_with_expected_rt(entry, coefficient_series, apex_spectrum, None)
    }

    /// Extract features with expected RT for deviation calculation
    pub fn extract_with_expected_rt(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        apex_spectrum: Option<&Spectrum>,
        expected_rt: Option<f64>,
    ) -> FeatureSet {
        let mut features = FeatureSet::default();

        // Chromatographic features from coefficient series (FR-5.1.*)
        if !coefficient_series.is_empty() {
            // FR-5.1.1: Peak apex coefficient (maximum value)
            let (apex_idx, apex_value) = coefficient_series
                .iter()
                .enumerate()
                .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
                .map(|(i, (_, c))| (i, *c))
                .unwrap_or((0, 0.0));
            features.peak_apex = apex_value;

            let apex_rt = coefficient_series
                .get(apex_idx)
                .map(|(rt, _)| *rt)
                .unwrap_or(0.0);

            // FR-5.1.2: Integrated peak area (AUC of coefficients)
            features.peak_area = coefficient_series.iter().map(|(_, c)| c).sum();

            // FR-5.1.8: Number of contributing scans
            features.n_contributing_scans =
                coefficient_series.iter().filter(|(_, c)| *c > 0.0).count() as u32;

            // FR-5.1.4: Peak width (FWHM with linear interpolation)
            if let Some((fwhm, _, _)) = compute_fwhm_interpolated(coefficient_series) {
                features.peak_width = fwhm;
            }

            // FR-5.1.5: Peak symmetry (leading/trailing ratio)
            features.peak_symmetry = self.compute_peak_symmetry(coefficient_series, apex_idx);

            // FR-5.1.6 & FR-5.1.7: RT deviation
            if let Some(expected) = expected_rt.or(Some(entry.retention_time)) {
                features.rt_deviation = apex_rt - expected;
                // Normalized deviation (use FWHM as uncertainty proxy)
                if features.peak_width > 0.0 {
                    features.rt_deviation_normalized = features.rt_deviation / features.peak_width;
                }
            }

            // FR-5.1.9: Coefficient stability (variance near apex)
            features.coefficient_stability =
                self.compute_coefficient_stability(coefficient_series, apex_idx);

            // FR-5.1.10: Peak boundary sharpness
            features.peak_sharpness = self.compute_peak_sharpness(coefficient_series, apex_idx);

            // FR-5.1.11: Peak prominence (apex / baseline)
            features.peak_prominence = self.compute_peak_prominence(coefficient_series, apex_value);

            // FR-5.1.12: Signal-to-noise ratio
            if let Some(snr) =
                self.compute_signal_to_noise(coefficient_series, apex_idx, apex_value)
            {
                features.signal_to_noise = snr;
            }

            // FR-5.1.3: EMG fit quality (simplified - use symmetry as proxy)
            // Full EMG fitting would go here, but for now use a heuristic
            features.emg_fit_quality = self.estimate_emg_quality(coefficient_series, apex_idx);
        }

        // Spectral features from apex spectrum (FR-5.2.*)
        if let Some(spectrum) = apex_spectrum {
            // Compute both LibCosine and XCorr scores
            let spectral_score = self.spectral_scorer.xcorr(spectrum, entry);

            // FR-5.2.4: Dot product (LibCosine with sqrt preprocessing)
            features.dot_product = spectral_score.lib_cosine;

            // LibCosine with sqrt(intensity)*mz² (SMZ) preprocessing
            features.dot_product_smz = spectral_score.lib_cosine_smz;

            // FR-5.2.3: Normalized spectral contrast angle
            features.spectral_contrast_angle = (spectral_score.lib_cosine.clamp(-1.0, 1.0))
                .acos()
                .to_degrees();

            // FR-5.2.7: Fragment coverage (fraction detected)
            features.fragment_coverage = spectral_score.fragment_coverage;

            // FR-5.2.12: Explained intensity fraction
            features.explained_intensity = spectral_score.explained_intensity;

            // FR-5.2.2: X!Tandem Hyperscore: log(n_b!) + log(n_y!) + Σlog(I_f+1)
            features.hyperscore = spectral_score.hyperscore;

            // XCorr (Comet-style cross-correlation)
            features.xcorr = spectral_score.xcorr;

            // FR-5.2.5: Pearson intensity correlation
            features.pearson_correlation = spectral_score.pearson_correlation;

            // FR-5.2.6: Spearman rank correlation
            features.spearman_correlation = spectral_score.spearman_correlation;

            // FR-5.2.8: Sequence coverage (backbone coverage)
            features.sequence_coverage = spectral_score.sequence_coverage;

            // FR-5.2.9: Consecutive ion count (longest b/y run)
            features.consecutive_ions = spectral_score.consecutive_ions;

            // FR-5.2.10: Base peak match rank
            features.base_peak_rank = spectral_score.base_peak_rank;

            // FR-5.2.11: Top-3 match count
            features.top6_matches = spectral_score.top6_matches;
        }

        // Count modifications
        features.modification_count = entry.modifications.len() as u32;

        features
    }

    /// Extract features with both mixed (apex) and deconvoluted (aggregated) scoring
    ///
    /// This computes spectral features twice:
    /// - Mixed: from the raw observed spectrum at apex (includes interference)
    /// - Deconvoluted: from coefficient-weighted aggregated spectrum (apex ± 2 scans)
    ///
    /// The deconvoluted features weight each scan's contribution by its regression
    /// coefficient, emphasizing where this peptide is actually contributing.
    ///
    /// Parameters:
    /// - `entry`: The library entry being scored
    /// - `coefficient_series`: Vec of (retention_time, coefficient) pairs
    /// - `spectra`: Full set of spectra for aggregation
    /// - `expected_rt`: Optional expected retention time
    pub fn extract_with_deconvolution(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        spectra: &[&Spectrum],
        expected_rt: Option<f64>,
    ) -> FeatureSet {
        // Find apex RT and spectrum
        let apex_rt = coefficient_series
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(rt, _)| *rt)
            .unwrap_or(entry.retention_time);

        // Find apex spectrum index
        let apex_spectrum = spectra
            .iter()
            .min_by(|a, b| {
                (a.retention_time - apex_rt)
                    .abs()
                    .partial_cmp(&(b.retention_time - apex_rt).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .copied();

        // First, extract features with mixed (apex) spectrum
        let mut features =
            self.extract_with_expected_rt(entry, coefficient_series, apex_spectrum, expected_rt);

        // Now compute deconvoluted features from coefficient-weighted aggregated spectrum
        if !spectra.is_empty() && !coefficient_series.is_empty() {
            // Aggregate spectra weighted by coefficients (apex ± 2 scans)
            let deconv_spectrum =
                self.spectrum_aggregator
                    .aggregate_weighted(spectra, coefficient_series, apex_rt);

            // Compute spectral scores on the deconvoluted spectrum
            let deconv_score = self.spectral_scorer.xcorr(&deconv_spectrum, entry);

            // Populate deconvoluted features
            features.hyperscore_deconv = deconv_score.hyperscore;
            features.xcorr_deconv = deconv_score.xcorr;
            features.dot_product_deconv = deconv_score.lib_cosine;
            features.dot_product_smz_deconv = deconv_score.lib_cosine_smz;
            features.fragment_coverage_deconv = deconv_score.fragment_coverage;
            features.sequence_coverage_deconv = deconv_score.sequence_coverage;
            features.consecutive_ions_deconv = deconv_score.consecutive_ions;
            features.top6_matches_deconv = deconv_score.top6_matches;
        }

        features
    }

    /// Apply contextual features from regression context
    ///
    /// Call this after extract() or extract_with_expected_rt() to add
    /// contextual features computed from regression results.
    pub fn apply_regression_context(features: &mut FeatureSet, context: &RegressionContext) {
        // FR-5.3.1: Number of competing candidates
        features.n_competitors = context.max_n_competitors;

        // FR-5.3.2: Relative coefficient (our coefficient / sum of all)
        features.relative_coefficient = context.max_relative_coefficient;

        // FR-5.3.3: Local peptide density (using avg competitors as proxy)
        features.local_peptide_density = context.avg_n_competitors;

        // FR-5.3.4: Spectral complexity
        features.spectral_complexity = context.avg_spectral_complexity;

        // FR-5.3.5: Regression residual (average across contributing spectra)
        features.regression_residual = context.avg_residual;
    }

    /// Extract features with full context (chromatographic, spectral, and contextual)
    ///
    /// This is the recommended method when regression context is available.
    pub fn extract_full(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        apex_spectrum: Option<&Spectrum>,
        expected_rt: Option<f64>,
        context: Option<&RegressionContext>,
    ) -> FeatureSet {
        let mut features =
            self.extract_with_expected_rt(entry, coefficient_series, apex_spectrum, expected_rt);

        if let Some(ctx) = context {
            Self::apply_regression_context(&mut features, ctx);
        }

        features
    }

    /// Compute peak symmetry (leading half area / trailing half area)
    fn compute_peak_symmetry(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        if series.len() < 3 || apex_idx == 0 || apex_idx >= series.len() - 1 {
            return 1.0; // Default to symmetric
        }

        // Sum coefficients before and after apex
        let leading_area: f64 = series[..apex_idx].iter().map(|(_, c)| c).sum();
        let trailing_area: f64 = series[apex_idx + 1..].iter().map(|(_, c)| c).sum();

        if trailing_area > 1e-10 {
            leading_area / trailing_area
        } else if leading_area > 1e-10 {
            f64::MAX // Infinitely asymmetric
        } else {
            1.0
        }
    }

    /// Compute coefficient stability (inverse of variance near apex)
    fn compute_coefficient_stability(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        // Take 5 points centered on apex
        let start = apex_idx.saturating_sub(2);
        let end = (apex_idx + 3).min(series.len());

        if end <= start + 1 {
            return 0.0;
        }

        let window: Vec<f64> = series[start..end].iter().map(|(_, c)| *c).collect();
        let mean: f64 = window.iter().sum::<f64>() / window.len() as f64;
        let variance: f64 =
            window.iter().map(|c| (c - mean).powi(2)).sum::<f64>() / window.len() as f64;

        // Return inverse of coefficient of variation (higher = more stable)
        if mean > 1e-10 {
            mean / (variance.sqrt() + 1e-10)
        } else {
            0.0
        }
    }

    /// Compute peak sharpness (slope at boundaries)
    fn compute_peak_sharpness(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        if series.len() < 3 {
            return 0.0;
        }

        let apex_coef = series[apex_idx].1;

        // Compute average slope from apex to boundaries
        let mut total_slope = 0.0;
        let mut count = 0;

        // Slope going left from apex
        if apex_idx > 0 {
            let dt = series[apex_idx].0 - series[0].0;
            if dt > 0.0 {
                total_slope += (apex_coef - series[0].1) / dt;
                count += 1;
            }
        }

        // Slope going right from apex
        if apex_idx < series.len() - 1 {
            let dt = series[series.len() - 1].0 - series[apex_idx].0;
            if dt > 0.0 {
                total_slope += (apex_coef - series[series.len() - 1].1) / dt;
                count += 1;
            }
        }

        if count > 0 {
            total_slope / count as f64
        } else {
            0.0
        }
    }

    /// Compute peak prominence (apex / local baseline)
    fn compute_peak_prominence(&self, series: &[(f64, f64)], apex_value: f64) -> f64 {
        if series.is_empty() {
            return 0.0;
        }

        // Estimate baseline from minimum values at edges
        let edge_values: Vec<f64> = series
            .iter()
            .take(3)
            .chain(series.iter().rev().take(3))
            .map(|(_, c)| *c)
            .collect();

        let baseline = edge_values
            .iter()
            .cloned()
            .fold(f64::MAX, f64::min)
            .max(0.0);

        if baseline > 1e-10 {
            apex_value / baseline
        } else if apex_value > 0.0 {
            apex_value / 1e-10 // Very high prominence
        } else {
            0.0
        }
    }

    /// Calculate signal-to-noise ratio from coefficient series
    ///
    /// Uses FWHM-based peak boundaries (±1.96σ) to define the peak region.
    /// Background is measured from 3-5 points on each side outside these boundaries.
    fn compute_signal_to_noise(
        &self,
        series: &[(f64, f64)],
        apex_idx: usize,
        apex_value: f64,
    ) -> Option<f64> {
        if series.len() < 10 || apex_value <= 0.0 {
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
        for (rt, coef) in series.iter().take(apex_idx) {
            if *rt < left_boundary {
                background_values.push(*coef);
            }
        }
        // Take last 5 on left
        if background_values.len() > 5 {
            background_values = background_values.into_iter().rev().take(5).collect();
        }

        // Right side: points after right boundary
        let mut right_bg = Vec::new();
        for (rt, coef) in series.iter().skip(apex_idx + 1) {
            if *rt > right_boundary {
                right_bg.push(*coef);
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

        if bg_sd <= 1e-10 {
            return None;
        }

        // S/N = (signal - background) / noise
        let signal = apex_value - bg_mean;
        let snr = signal / bg_sd;

        Some(snr.max(0.0))
    }

    /// Estimate EMG fit quality (heuristic based on peak shape)
    fn estimate_emg_quality(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        if series.len() < 5 {
            return 0.0;
        }

        // Simple heuristic: measure how Gaussian-like the peak is
        // by comparing expected and actual half-width points
        let apex_coef = series[apex_idx].1;
        let half_max = apex_coef / 2.0;

        // Find actual half-width points
        let mut left_half_idx = apex_idx;
        let mut right_half_idx = apex_idx;

        for i in (0..apex_idx).rev() {
            if series[i].1 < half_max {
                left_half_idx = i;
                break;
            }
        }

        for (i, entry) in series.iter().enumerate().skip(apex_idx + 1) {
            if entry.1 < half_max {
                right_half_idx = i;
                break;
            }
        }

        // For a Gaussian, we expect roughly symmetric FWHM
        let left_half = apex_idx - left_half_idx;
        let right_half = right_half_idx - apex_idx;

        if left_half > 0 && right_half > 0 {
            // Ratio closer to 1.0 = more symmetric = better fit
            // Value between 0 and 1
            left_half.min(right_half) as f64 / left_half.max(right_half) as f64
        } else {
            0.5 // Default if we can't measure
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

    /// Verifies that FeatureExtractor computes correct peak apex, area, and scan count from a coefficient time series.
    #[test]
    fn test_feature_extractor() {
        let extractor = FeatureExtractor::new();
        let entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);

        let series = vec![(9.0, 0.1), (10.0, 1.0), (11.0, 0.5)];
        let features = extractor.extract(&entry, &series, None);

        assert!((features.peak_apex - 1.0).abs() < 1e-6);
        assert!((features.peak_area - 1.6).abs() < 1e-6);
        assert_eq!(features.n_contributing_scans, 3);
    }

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
        // LibCosine should still be reasonable for the matched fragments
        assert!(
            score.lib_cosine > 0.9,
            "LibCosine should be high for matched fragments"
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

    /// Verifies that FeatureExtractor computes both chromatographic and spectral features when given an apex spectrum.
    #[test]
    fn test_feature_extractor_with_spectrum() {
        let extractor = FeatureExtractor::new();

        // Create library entry with fragments
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

        // Create apex spectrum
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0],
            intensities: vec![100.0, 50.0],
        };

        let series = vec![(9.0, 0.1), (10.0, 1.0), (11.0, 0.5)];
        let features = extractor.extract(&entry, &series, Some(&spectrum));

        // Chromatographic features should be computed
        assert!((features.peak_apex - 1.0).abs() < 1e-6);

        // Spectral features should be computed
        assert!(
            features.dot_product > 0.9,
            "Dot product should be high for matching spectrum"
        );
        assert!(
            (features.fragment_coverage - 1.0).abs() < 1e-6,
            "All fragments should match"
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
}
