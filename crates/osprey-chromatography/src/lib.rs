//! Osprey Chromatography - Peak detection and chromatographic analysis
//!
//! This crate provides:
//! - Peak detection in coefficient time series
//! - EMG (exponentially modified Gaussian) fitting
//! - Peak boundary determination
//! - Chromatographic feature extraction
//! - RT calibration using LOESS regression
//! - Mass calibration for MS1/MS2 measurements

pub mod calibration;

use osprey_core::{PeakBoundaries, PeakQuality, Result, XICPeakBounds};

// Re-export calibration types
pub use calibration::{
    apply_mz_calibration,
    apply_spectrum_calibration,
    calculate_mz_calibration,
    calculate_ppm_error,
    calibrated_tolerance,
    calibration_filename,
    calibration_filename_for_input,
    calibration_path_for_input,
    load_calibration,
    // I/O
    save_calibration,
    CalibrationMetadata,
    // Core types
    CalibrationParams,
    IsolationScheme,
    MzCalibration,
    // Mass calibration
    MzQCData,
    // RT calibration
    RTCalibration,
    RTCalibrationMethod,
    RTCalibrationParams,
    RTCalibrationStats,
    RTCalibrator,
    RTCalibratorConfig,
    RTModelParams,
    RTStratifiedSampler,
};

/// Compute trapezoidal area of a (RT, value) time series slice.
///
/// area = Σ (v_i + v_{i+1}) / 2 × (t_{i+1} - t_i)
pub fn trapezoidal_area(series: &[(f64, f64)]) -> f64 {
    if series.len() < 2 {
        return series.first().map(|(_, c)| *c).unwrap_or(0.0);
    }
    let mut area = 0.0;
    for i in 0..series.len() - 1 {
        let dt = series[i + 1].0 - series[i].0;
        let avg_height = (series[i].1 + series[i + 1].1) / 2.0;
        area += avg_height * dt;
    }
    area
}

/// Simple peak detector for coefficient time series
#[derive(Debug, Default)]
pub struct PeakDetector {
    /// Minimum peak height threshold
    min_height: f64,
    /// Minimum peak width in scans
    min_width: usize,
}

impl PeakDetector {
    /// Create a new peak detector
    pub fn new() -> Self {
        Self {
            min_height: 0.01,
            min_width: 3,
        }
    }

    /// Set minimum peak height
    pub fn with_min_height(mut self, height: f64) -> Self {
        self.min_height = height;
        self
    }

    /// Detect peaks in a coefficient time series
    ///
    /// Input: vector of (retention_time, coefficient) pairs
    /// Output: vector of peak boundaries
    pub fn detect(&self, series: &[(f64, f64)]) -> Vec<PeakBoundaries> {
        if series.is_empty() {
            return Vec::new();
        }

        let mut peaks = Vec::new();

        // Find local maxima above threshold
        let mut in_peak = false;
        let mut peak_start_idx = 0;
        let mut apex_idx = 0;
        let mut apex_value = 0.0;

        for (i, (_, coef)) in series.iter().enumerate() {
            if *coef >= self.min_height {
                if !in_peak {
                    // Start of new peak
                    in_peak = true;
                    peak_start_idx = i;
                    apex_idx = i;
                    apex_value = *coef;
                } else if *coef > apex_value {
                    // New apex
                    apex_idx = i;
                    apex_value = *coef;
                }
            } else if in_peak {
                // End of peak
                let peak_end_idx = i.saturating_sub(1);
                if peak_end_idx - peak_start_idx + 1 >= self.min_width {
                    // Calculate peak area using trapezoidal integration
                    let peak_slice = &series[peak_start_idx..=peak_end_idx];
                    let area = trapezoidal_area(peak_slice);

                    peaks.push(PeakBoundaries {
                        start_rt: series[peak_start_idx].0,
                        end_rt: series[peak_end_idx].0,
                        apex_rt: series[apex_idx].0,
                        apex_coefficient: apex_value,
                        integrated_area: area,
                        peak_quality: PeakQuality::default(),
                    });
                }
                in_peak = false;
            }
        }

        // Handle peak at end of series
        if in_peak {
            let peak_end_idx = series.len() - 1;
            if peak_end_idx - peak_start_idx + 1 >= self.min_width {
                let peak_slice = &series[peak_start_idx..=peak_end_idx];
                let area = trapezoidal_area(peak_slice);

                peaks.push(PeakBoundaries {
                    start_rt: series[peak_start_idx].0,
                    end_rt: series[peak_end_idx].0,
                    apex_rt: series[apex_idx].0,
                    apex_coefficient: apex_value,
                    integrated_area: area,
                    peak_quality: PeakQuality::default(),
                });
            }
        }

        peaks
    }

    /// Find the best peak near the expected retention time
    pub fn find_best_peak<'a>(
        &self,
        peaks: &'a [PeakBoundaries],
        expected_rt: f64,
        rt_tolerance: f64,
    ) -> Option<&'a PeakBoundaries> {
        peaks
            .iter()
            .filter(|p| (p.apex_rt - expected_rt).abs() <= rt_tolerance)
            .max_by(|a, b| {
                a.apex_coefficient
                    .partial_cmp(&b.apex_coefficient)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    }
}

// =============================================================================
// DIA-NN-style peak detection with valley-based boundary determination
// =============================================================================

/// Smooth a time series with 5-point Savitzky-Golay quadratic filter.
///
/// Coefficients: [-3, 12, 17, 12, -3] / 35
/// Preserves peak position and shape better than triangular kernels by fitting
/// a local quadratic polynomial at each point. Endpoints (first 2 and last 2
/// points) are left unsmoothed. Negative values from the filter are clamped
/// to zero. Series shorter than 5 points are returned unsmoothed.
pub fn smooth_savitzky_golay(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    if n < 5 {
        return values.to_vec();
    }
    let mut out = vec![0.0; n];
    out[0] = values[0];
    out[1] = values[1];
    out[n - 2] = values[n - 2];
    out[n - 1] = values[n - 1];
    for i in 2..n - 2 {
        out[i] = ((-3.0 * values[i - 2]
            + 12.0 * values[i - 1]
            + 17.0 * values[i]
            + 12.0 * values[i + 1]
            + -3.0 * values[i + 2])
            / 35.0)
            .max(0.0);
    }
    out
}

/// Detect a peak in a fragment XIC using DIA-NN-style adaptive boundary detection.
///
/// Algorithm:
/// 1. Smooth the XIC with weighted moving average
/// 2. Find the apex (local maximum nearest to expected_rt, or global max)
/// 3. Walk left/right from apex to find boundaries using:
///    - Intensity falloff: stop when intensity < apex / peak_boundary
///    - Valley detection: stop at local minimum that is <50% of apex AND <50% of neighbor
///
/// # Arguments
/// * `xic` - (RT, intensity) pairs, must be sorted by RT
/// * `min_height` - Minimum apex intensity after smoothing (default: 0.01)
/// * `peak_boundary` - Intensity divisor for boundary threshold (default: 5.0 → 20% of apex)
/// * `expected_rt` - If provided, prefer peaks nearest to this RT
pub fn detect_xic_peak(
    xic: &[(f64, f64)],
    min_height: f64,
    peak_boundary: f64,
    expected_rt: Option<f64>,
) -> Option<XICPeakBounds> {
    if xic.len() < 3 {
        return None;
    }

    // Extract intensity values and smooth
    let raw_intensities: Vec<f64> = xic.iter().map(|(_, v)| *v).collect();
    let smoothed = smooth_savitzky_golay(&raw_intensities);

    // Find apex: local maximum above min_height
    let apex_idx = find_apex(&smoothed, min_height, expected_rt, xic)?;
    let apex_intensity = smoothed[apex_idx];

    // Estimate local background from minimum of smoothed intensity.
    // For low S/N peaks, background can be 20-30% of apex intensity, so computing
    // the boundary threshold relative to the signal above background prevents
    // boundaries from extending into baseline noise.
    let background = smoothed
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(0.0);
    let signal_above_bg = (apex_intensity - background).max(0.0);
    let boundary_threshold = background + signal_above_bg / peak_boundary;

    // Walk left to find start boundary
    let mut start_idx = walk_boundary_left(&smoothed, apex_idx, apex_intensity, boundary_threshold);

    // Walk right to find end boundary
    let mut end_idx = walk_boundary_right(&smoothed, apex_idx, apex_intensity, boundary_threshold);

    // Asymmetric FWHM-based boundary capping: prevents boundaries from extending
    // too far on slowly decaying tails. Cap at apex +/- cap_factor * half_width.
    // Factor 2.0 with HWHM ≈ 1.177σ gives cap at ~2.35σ (~98% Gaussian area).
    if let Some((left_hw, right_hw)) = compute_asymmetric_half_widths(&smoothed, xic, apex_idx) {
        let cap_factor = 2.0;
        let apex_rt = xic[apex_idx].0;

        // Cap left boundary
        let min_start_rt = apex_rt - cap_factor * left_hw;
        if xic[start_idx].0 < min_start_rt {
            if let Some((i, _)) = xic[start_idx..apex_idx]
                .iter()
                .enumerate()
                .find(|(_, (rt, _))| *rt >= min_start_rt)
            {
                start_idx += i;
            }
        }

        // Cap right boundary
        let max_end_rt = apex_rt + cap_factor * right_hw;
        if xic[end_idx].0 > max_end_rt {
            if let Some((i, _)) = xic[apex_idx..=end_idx]
                .iter()
                .enumerate()
                .rev()
                .find(|(_, (rt, _))| *rt <= max_end_rt)
            {
                end_idx = apex_idx + i;
            }
        }
    }

    // Need at least 3 points in the peak
    if end_idx - start_idx + 1 < 3 {
        return None;
    }

    // Compute area via trapezoidal integration
    let area = trapezoidal_area(&xic[start_idx..=end_idx]);

    // Compute SNR using raw intensities: apex vs background outside boundaries
    let snr = compute_snr(&raw_intensities, apex_idx, start_idx, end_idx);

    Some(XICPeakBounds {
        apex_rt: xic[apex_idx].0,
        apex_intensity,
        apex_index: apex_idx,
        start_rt: xic[start_idx].0,
        end_rt: xic[end_idx].0,
        start_index: start_idx,
        end_index: end_idx,
        area,
        signal_to_noise: snr,
    })
}

/// Compute asymmetric half-widths at half-maximum on a smoothed intensity series.
///
/// Unlike symmetric FWHM, this returns separate left and right half-widths,
/// which naturally captures chromatographic peak tailing. Uses linear
/// interpolation to find exact half-height crossing points.
///
/// # Returns
/// `Some((left_half_width, right_half_width))` in RT units (minutes), measured
/// from the apex. Returns `None` if either crossing cannot be found (e.g.,
/// signal never drops below 50% of apex on one side).
fn compute_asymmetric_half_widths(
    smoothed: &[f64],
    xic: &[(f64, f64)],
    apex_idx: usize,
) -> Option<(f64, f64)> {
    let apex_val = smoothed[apex_idx];
    if apex_val <= 0.0 {
        return None;
    }
    let half = apex_val / 2.0;
    let apex_rt = xic[apex_idx].0;

    // Scan left from apex to find half-height crossing
    let left_hw = {
        let mut found = None;
        for i in (1..=apex_idx).rev() {
            if smoothed[i] >= half && smoothed[i - 1] < half {
                let denom = smoothed[i] - smoothed[i - 1];
                let crossing_rt = if denom > 1e-30 {
                    let frac = (smoothed[i] - half) / denom;
                    xic[i].0 - frac * (xic[i].0 - xic[i - 1].0)
                } else {
                    xic[i].0
                };
                found = Some(apex_rt - crossing_rt);
                break;
            }
        }
        found
    };

    // Scan right from apex to find half-height crossing
    let right_hw = {
        let mut found = None;
        for i in apex_idx..smoothed.len() - 1 {
            if smoothed[i] >= half && smoothed[i + 1] < half {
                let denom = smoothed[i] - smoothed[i + 1];
                let crossing_rt = if denom > 1e-30 {
                    let frac = (smoothed[i] - half) / denom;
                    xic[i].0 + frac * (xic[i + 1].0 - xic[i].0)
                } else {
                    xic[i].0
                };
                found = Some(crossing_rt - apex_rt);
                break;
            }
        }
        found
    };

    match (left_hw, right_hw) {
        (Some(l), Some(r)) if l > 0.0 && r > 0.0 => Some((l, r)),
        // One side found: use it for both (assume roughly symmetric)
        (Some(l), None) if l > 0.0 => Some((l, l)),
        (None, Some(r)) if r > 0.0 => Some((r, r)),
        _ => None,
    }
}

/// Find the best apex in a smoothed intensity series.
///
/// Returns the global maximum among local maxima above `min_height`.
/// The `expected_rt` parameter is accepted for API compatibility but not used
/// for apex selection — the caller (pipeline) should use `detect_all_xic_peaks`
/// and pick the best candidate by coelution score instead.
fn find_apex(
    smoothed: &[f64],
    min_height: f64,
    _expected_rt: Option<f64>,
    _xic: &[(f64, f64)],
) -> Option<usize> {
    let n = smoothed.len();

    // Collect all local maxima above threshold
    let mut candidates: Vec<(usize, f64)> = Vec::new();
    for i in 1..n - 1 {
        if smoothed[i] >= min_height
            && smoothed[i] >= smoothed[i - 1]
            && smoothed[i] >= smoothed[i + 1]
        {
            candidates.push((i, smoothed[i]));
        }
    }
    // Also check endpoints
    if n >= 2 && smoothed[0] >= min_height && smoothed[0] >= smoothed[1] {
        candidates.push((0, smoothed[0]));
    }
    if n >= 2 && smoothed[n - 1] >= min_height && smoothed[n - 1] >= smoothed[n - 2] {
        candidates.push((n - 1, smoothed[n - 1]));
    }

    if candidates.is_empty() {
        return None;
    }

    // Pick global maximum
    candidates
        .iter()
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|&(idx, _)| idx)
}

/// Detect ALL candidate peaks in a XIC, returning them sorted by apex intensity (descending).
///
/// Each local maximum above `min_height` produces a candidate peak with its own
/// boundaries (valley detection + FWHM capping). The caller should evaluate each
/// candidate using coelution scoring and pick the best one.
pub fn detect_all_xic_peaks(
    xic: &[(f64, f64)],
    min_height: f64,
    peak_boundary: f64,
) -> Vec<XICPeakBounds> {
    if xic.len() < 3 {
        return Vec::new();
    }

    let raw_intensities: Vec<f64> = xic.iter().map(|(_, v)| *v).collect();
    let smoothed = smooth_savitzky_golay(&raw_intensities);
    let n = smoothed.len();

    // Find all local maxima above threshold
    let mut apex_indices: Vec<(usize, f64)> = Vec::new();
    for i in 1..n - 1 {
        if smoothed[i] >= min_height
            && smoothed[i] >= smoothed[i - 1]
            && smoothed[i] >= smoothed[i + 1]
        {
            apex_indices.push((i, smoothed[i]));
        }
    }
    if n >= 2 && smoothed[0] >= min_height && smoothed[0] >= smoothed[1] {
        apex_indices.push((0, smoothed[0]));
    }
    if n >= 2 && smoothed[n - 1] >= min_height && smoothed[n - 1] >= smoothed[n - 2] {
        apex_indices.push((n - 1, smoothed[n - 1]));
    }

    // Sort by intensity descending so best candidates come first
    apex_indices.sort_by(|a, b| b.1.total_cmp(&a.1));

    let background = smoothed
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(0.0);

    let mut peaks = Vec::new();
    for &(apex_idx, apex_intensity) in &apex_indices {
        let signal_above_bg = (apex_intensity - background).max(0.0);
        let boundary_threshold = background + signal_above_bg / peak_boundary;

        let mut start_idx =
            walk_boundary_left(&smoothed, apex_idx, apex_intensity, boundary_threshold);
        let mut end_idx =
            walk_boundary_right(&smoothed, apex_idx, apex_intensity, boundary_threshold);

        // FWHM capping
        if let Some((left_hw, right_hw)) = compute_asymmetric_half_widths(&smoothed, xic, apex_idx)
        {
            let cap_factor = 2.0;
            let apex_rt = xic[apex_idx].0;

            let min_start_rt = apex_rt - cap_factor * left_hw;
            if xic[start_idx].0 < min_start_rt {
                if let Some((i, _)) = xic[start_idx..apex_idx]
                    .iter()
                    .enumerate()
                    .find(|(_, (rt, _))| *rt >= min_start_rt)
                {
                    start_idx += i;
                }
            }

            let max_end_rt = apex_rt + cap_factor * right_hw;
            if xic[end_idx].0 > max_end_rt {
                if let Some((i, _)) = xic[apex_idx..=end_idx]
                    .iter()
                    .enumerate()
                    .rev()
                    .find(|(_, (rt, _))| *rt <= max_end_rt)
                {
                    end_idx = apex_idx + i;
                }
            }
        }

        if end_idx - start_idx + 1 < 3 {
            continue;
        }

        let area = trapezoidal_area(&xic[start_idx..=end_idx]);
        let snr = compute_snr(&raw_intensities, apex_idx, start_idx, end_idx);

        peaks.push(XICPeakBounds {
            apex_rt: xic[apex_idx].0,
            apex_intensity,
            apex_index: apex_idx,
            start_rt: xic[start_idx].0,
            end_rt: xic[end_idx].0,
            start_index: start_idx,
            end_index: end_idx,
            area,
            signal_to_noise: snr,
        });
    }

    peaks
}

/// Walk left from apex to find start boundary using DIA-NN-style valley detection.
///
/// Stops at whichever comes first:
/// - Intensity drops below `boundary_threshold` (apex / peak_boundary)
/// - Valley detected: local minimum is <50% of apex AND <50% of the next rising neighbor
fn walk_boundary_left(
    smoothed: &[f64],
    apex_idx: usize,
    apex_intensity: f64,
    boundary_threshold: f64,
) -> usize {
    let mut valley = smoothed[apex_idx];
    let mut valley_pos = apex_idx;

    for i in (0..apex_idx).rev() {
        if smoothed[i] < valley {
            valley = smoothed[i];
            valley_pos = i;
        } else if valley < apex_intensity / 2.0 && valley < smoothed[i] / 2.0 {
            // Valley detected: intensity dropped to <50% of apex and is <50% of rising neighbor
            return valley_pos;
        }

        if smoothed[i] < boundary_threshold {
            return i;
        }
    }

    0 // Reached beginning of series
}

/// Walk right from apex to find end boundary using DIA-NN-style valley detection.
fn walk_boundary_right(
    smoothed: &[f64],
    apex_idx: usize,
    apex_intensity: f64,
    boundary_threshold: f64,
) -> usize {
    let n = smoothed.len();
    let mut valley = smoothed[apex_idx];
    let mut valley_pos = apex_idx;

    for (i, &val) in smoothed.iter().enumerate().skip(apex_idx + 1) {
        if val < valley {
            valley = val;
            valley_pos = i;
        } else if valley < apex_intensity / 2.0 && valley < val / 2.0 {
            // Valley detected
            return valley_pos;
        }

        if val < boundary_threshold {
            return i;
        }
    }

    n - 1 // Reached end of series
}

/// Compute signal-to-noise ratio for a peak using raw intensities.
///
/// Uses 5 points immediately before the peak start and 5 points immediately
/// after the peak end as background. This captures the local baseline the peak
/// sits on without reaching into distant interferences.
///
/// Signal = apex_intensity - mean(background)
/// Noise  = SD(background)
/// SNR    = Signal / Noise
///
/// This is the canonical SNR function used by coelution and calibration scoring.
pub fn compute_snr(intensities: &[f64], apex_idx: usize, start_idx: usize, end_idx: usize) -> f64 {
    // Collect 5 background points immediately before peak start
    let mut bg_points: Vec<f64> = Vec::new();
    let left_start = start_idx.saturating_sub(5);
    for &val in &intensities[left_start..start_idx] {
        bg_points.push(val);
    }

    // Collect 5 background points immediately after peak end
    let right_end = (end_idx + 6).min(intensities.len());
    for &val in &intensities[(end_idx + 1)..right_end] {
        bg_points.push(val);
    }

    if bg_points.is_empty() {
        // No background points at all — can't estimate noise
        return 0.0;
    }

    let n = bg_points.len() as f64;
    let bg_mean = bg_points.iter().sum::<f64>() / n;
    let bg_sd = (bg_points.iter().map(|x| (x - bg_mean).powi(2)).sum::<f64>() / n).sqrt();

    if bg_sd > 1e-10 {
        ((intensities[apex_idx] - bg_mean) / bg_sd).max(0.0)
    } else {
        // Background is flat (no noise variance) — report signal above baseline
        // relative to the apex intensity as a rough SNR
        let signal = intensities[apex_idx] - bg_mean;
        if signal > 0.0 && bg_mean > 1e-10 {
            signal / bg_mean
        } else if signal > 0.0 {
            signal * 100.0
        } else {
            0.0
        }
    }
}

/// Isolate a single peak from a profile by walking outward from the apex index
/// and stopping when a significant rise is detected (indicating another peak).
///
/// This is critical for EMG fitting: the fitter assumes a single-peak model,
/// so passing a multi-peak profile (common in wide extraction windows) causes
/// either a broad fit or convergence failure. This function crops to the
/// valley between peaks while preserving the full tails of the target peak.
///
/// Algorithm: track the running minimum as we walk away from the apex. If the
/// profile rises more than 5% of the apex intensity above that minimum, we've
/// entered another peak — crop at the minimum. Small noise bumps on the tail
/// won't trigger this threshold, so single-peak tails are fully preserved.
///
/// Returns `(start_index, end_index)` of the isolated peak region.
pub fn isolate_peak_region(profile: &[(f64, f64)], apex_index: usize) -> (usize, usize) {
    let n = profile.len();
    if n < 3 || apex_index >= n {
        return (0, n.saturating_sub(1));
    }

    let apex_val = profile[apex_index].1;
    let rise_threshold = apex_val * 0.05; // 5% of apex = significant rise

    // Walk left from apex: detect another peak rising
    let mut left = 0;
    let mut running_min = apex_val;
    let mut min_idx = apex_index;
    for i in (0..apex_index).rev() {
        if profile[i].1 < running_min {
            running_min = profile[i].1;
            min_idx = i;
        }
        // Intensity has risen significantly above the running minimum →
        // we've entered another peak. Crop at the minimum (valley).
        if profile[i].1 - running_min > rise_threshold {
            left = min_idx;
            break;
        }
    }

    // Walk right from apex: detect another peak rising
    let mut right = n - 1;
    let mut running_min = apex_val;
    let mut min_idx = apex_index;
    for (i, &(_, intensity)) in profile.iter().enumerate().skip(apex_index + 1) {
        if intensity < running_min {
            running_min = intensity;
            min_idx = i;
        }
        if intensity - running_min > rise_threshold {
            right = min_idx;
            break;
        }
    }

    (left, right)
}

// =============================================================================
// EMG (Exponentially Modified Gaussian) Peak Fitting
// =============================================================================

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
/// Maximum error < 1.5e-7.
fn erfc_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let result = poly * (-x * x).exp();
    if x >= 0.0 {
        result
    } else {
        2.0 - result
    }
}

/// Standard normal CDF: Φ(x) = P(Z ≤ x).
fn normal_cdf(x: f64) -> f64 {
    0.5 * erfc_approx(-x / std::f64::consts::SQRT_2)
}

/// EMG probability density function (unnormalized, with amplitude A).
///
/// f(t) = (A/(2τ)) * exp(σ²/(2τ²) - (t-μ)/τ) * erfc((1/√2)(σ/τ - (t-μ)/σ))
///
/// Models an LC chromatographic peak as the convolution of a Gaussian (μ, σ)
/// with an exponential decay (τ). τ > 0 gives right-tailing peaks.
pub fn emg_pdf(t: f64, mu: f64, sigma: f64, tau: f64, amplitude: f64) -> f64 {
    if sigma <= 0.0 || tau <= 0.0 {
        return 0.0;
    }
    let s_over_tau = sigma / tau;
    let z = (t - mu) / sigma;
    let exponent = 0.5 * s_over_tau * s_over_tau - z * s_over_tau;
    // Prevent overflow: exp(709) ≈ f64::MAX
    if exponent > 500.0 {
        return 0.0;
    }
    let erfc_arg = (s_over_tau - z) / std::f64::consts::SQRT_2;
    amplitude / (2.0 * tau) * exponent.exp() * erfc_approx(erfc_arg)
}

/// EMG cumulative distribution function (normalized, amplitude-independent).
///
/// F(t) = Φ((t-μ)/σ) - exp(σ²/(2τ²) - (t-μ)/τ) * Φ((t-μ)/σ - σ/τ)
pub fn emg_cdf(t: f64, mu: f64, sigma: f64, tau: f64) -> f64 {
    if sigma <= 0.0 || tau <= 0.0 {
        return 0.0;
    }
    let z = (t - mu) / sigma;
    let s_over_tau = sigma / tau;
    let exponent = 0.5 * s_over_tau * s_over_tau - z * s_over_tau;
    if exponent > 500.0 {
        return normal_cdf(z);
    }
    normal_cdf(z) - exponent.exp() * normal_cdf(z - s_over_tau)
}

/// EMG (Exponentially Modified Gaussian) peak fitter.
///
/// Fits an EMG model to a chromatographic profile using Levenberg-Marquardt
/// optimization. The EMG is the gold standard for LC peaks since they're
/// typically asymmetric (Gaussian convolved with exponential decay).
///
/// Workflow:
/// 1. Initialize from intensity-weighted moments (M₁, M₂, M₃ → μ, σ, τ)
/// 2. Refine via Levenberg-Marquardt (4 params: μ, σ, τ, A)
/// 3. Derive peak boundaries from the fitted CDF quantiles
pub struct EmgFitter {
    max_iterations: usize,
}

impl EmgFitter {
    pub fn new() -> Self {
        Self { max_iterations: 50 }
    }

    /// Fit EMG to a chromatographic profile (RT, intensity pairs).
    ///
    /// The profile should be the consensus elution shape (e.g., from Tukey
    /// median polish). Points with zero intensity are included in the fit.
    pub fn fit(&self, profile: &[(f64, f64)]) -> Result<EmgParameters> {
        if profile.len() < 5 {
            return Err(osprey_core::OspreyError::PeakDetectionError(
                "EMG fit requires at least 5 data points".to_string(),
            ));
        }

        let rts: Vec<f64> = profile.iter().map(|(rt, _)| *rt).collect();
        let intensities: Vec<f64> = profile.iter().map(|(_, v)| *v).collect();

        // Find apex and its value for amplitude initialization
        let (apex_idx, apex_val) = intensities
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, &v)| (i, v))
            .unwrap_or((0, 0.0));

        if apex_val <= 0.0 {
            return Err(osprey_core::OspreyError::PeakDetectionError(
                "EMG fit: no positive signal in profile".to_string(),
            ));
        }

        // Initialize from intensity-weighted moments
        let (mu_init, sigma_init, tau_init) =
            Self::initialize_from_moments(&rts, &intensities, rts[apex_idx]);

        let mut params = [mu_init, sigma_init, tau_init, apex_val];

        // Levenberg-Marquardt optimization
        let mut lambda = 0.001_f64;
        let mut best_ssr = Self::sum_squared_residuals(&rts, &intensities, &params);

        for _iter in 0..self.max_iterations {
            // Compute Jacobian via forward finite differences
            let eps = [1e-6, 1e-6, 1e-6, apex_val * 1e-6];
            let n = rts.len();
            let mut jac = vec![0.0; n * 4]; // column-major: jac[row + col*n]
            let residuals = Self::residuals(&rts, &intensities, &params);

            for p in 0..4 {
                let mut params_plus = params;
                params_plus[p] += eps[p];
                Self::clamp_params(&mut params_plus);
                let res_plus = Self::residuals(&rts, &intensities, &params_plus);
                for i in 0..n {
                    jac[i + p * n] = (res_plus[i] - residuals[i]) / eps[p];
                }
            }

            // Solve (JᵀJ + λI)δ = Jᵀr using 4x4 system
            let mut jtj = [0.0; 16]; // 4x4 row-major
            let mut jtr = [0.0; 4];
            for p in 0..4 {
                for q in 0..4 {
                    let mut sum = 0.0;
                    for i in 0..n {
                        sum += jac[i + p * n] * jac[i + q * n];
                    }
                    jtj[p * 4 + q] = sum;
                }
                let mut sum = 0.0;
                for i in 0..n {
                    sum += jac[i + p * n] * residuals[i];
                }
                jtr[p] = sum;
            }

            // Add damping: JᵀJ + λ * diag(JᵀJ)
            for p in 0..4 {
                jtj[p * 4 + p] *= 1.0 + lambda;
            }

            // Solve 4x4 system via Gaussian elimination
            let delta = match Self::solve_4x4(&mut jtj, &mut jtr) {
                Some(d) => d,
                None => break, // Singular matrix, stop
            };

            // Trial step
            let mut new_params = [
                params[0] - delta[0],
                params[1] - delta[1],
                params[2] - delta[2],
                params[3] - delta[3],
            ];
            Self::clamp_params(&mut new_params);

            let new_ssr = Self::sum_squared_residuals(&rts, &intensities, &new_params);
            if new_ssr < best_ssr {
                params = new_params;
                best_ssr = new_ssr;
                lambda *= 0.1;
                lambda = lambda.max(1e-10);

                // Check convergence
                let max_delta = delta.iter().map(|d| d.abs()).fold(0.0_f64, f64::max);
                if max_delta < 1e-8 {
                    break;
                }
            } else {
                lambda *= 10.0;
                if lambda > 1e10 {
                    break;
                }
            }
        }

        // Compute R²
        let mean_int = intensities.iter().sum::<f64>() / intensities.len() as f64;
        let ss_tot: f64 = intensities.iter().map(|v| (v - mean_int).powi(2)).sum();
        let r_squared = if ss_tot > 1e-30 {
            1.0 - best_ssr / ss_tot
        } else {
            0.0
        };

        Ok(EmgParameters {
            mu: params[0],
            sigma: params[1],
            tau: params[2],
            amplitude: params[3],
            r_squared,
        })
    }

    /// Compute peak boundaries from EMG CDF quantiles.
    ///
    /// Returns (left_rt, right_rt) such that the CDF covers `coverage` of the
    /// peak area (e.g., 0.95 → 2.5th to 97.5th percentile).
    pub fn boundaries(params: &EmgParameters, coverage: f64) -> Option<(f64, f64)> {
        if params.sigma <= 0.0 || params.tau <= 0.0 {
            return None;
        }
        let lower_q = (1.0 - coverage) / 2.0;
        let upper_q = 1.0 - lower_q;

        // EMG mean and variance for search range
        let emg_mean = params.mu + params.tau;
        let emg_var = params.sigma * params.sigma + params.tau * params.tau;
        let emg_sd = emg_var.sqrt();
        let search_range = 10.0 * emg_sd; // generous search range

        let left = Self::find_cdf_quantile(
            params.mu,
            params.sigma,
            params.tau,
            lower_q,
            emg_mean - search_range,
            emg_mean,
        )?;
        let right = Self::find_cdf_quantile(
            params.mu,
            params.sigma,
            params.tau,
            upper_q,
            emg_mean,
            emg_mean + search_range,
        )?;

        Some((left, right))
    }

    /// Initialize EMG parameters from intensity-weighted moments.
    ///
    /// M₁ = μ + τ, M₂ = σ² + τ², M₃ = 2τ³
    fn initialize_from_moments(rts: &[f64], intensities: &[f64], apex_rt: f64) -> (f64, f64, f64) {
        let total: f64 = intensities.iter().sum();
        if total <= 0.0 {
            let dt = if rts.len() >= 2 {
                (rts.last().unwrap() - rts.first().unwrap()) / rts.len() as f64
            } else {
                0.1
            };
            return (apex_rt, dt, dt * 0.1);
        }

        // Intensity-weighted moments
        let m1: f64 = rts.iter().zip(intensities).map(|(r, v)| r * v).sum::<f64>() / total;
        let m2: f64 = rts
            .iter()
            .zip(intensities)
            .map(|(r, v)| (r - m1).powi(2) * v)
            .sum::<f64>()
            / total;
        let m3: f64 = rts
            .iter()
            .zip(intensities)
            .map(|(r, v)| (r - m1).powi(3) * v)
            .sum::<f64>()
            / total;

        // From moments: τ = (M₃/2)^(1/3), σ² = M₂ - τ², μ = M₁ - τ
        let tau = if m3 > 0.0 {
            (m3 / 2.0).cbrt()
        } else {
            // Symmetric or left-tailing: use small τ
            m2.sqrt() * 0.1
        };
        let sigma_sq = (m2 - tau * tau).max(1e-6);
        let sigma = sigma_sq.sqrt();
        let mu = m1 - tau;

        // Sanity check: μ should be near apex
        let dt = if rts.len() >= 2 {
            (rts.last().unwrap() - rts.first().unwrap()) / rts.len() as f64
        } else {
            0.1
        };
        let mu = if (mu - apex_rt).abs() > 5.0 * m2.sqrt() {
            apex_rt // Fall back to apex if moment estimate is wild
        } else {
            mu
        };
        let sigma = sigma.max(dt * 0.5);
        let tau = tau.max(dt * 0.01);

        (mu, sigma, tau)
    }

    /// Compute residuals: observed - predicted
    fn residuals(rts: &[f64], observed: &[f64], params: &[f64; 4]) -> Vec<f64> {
        rts.iter()
            .zip(observed)
            .map(|(&t, &obs)| obs - emg_pdf(t, params[0], params[1], params[2], params[3]))
            .collect()
    }

    /// Sum of squared residuals
    fn sum_squared_residuals(rts: &[f64], observed: &[f64], params: &[f64; 4]) -> f64 {
        rts.iter()
            .zip(observed)
            .map(|(&t, &obs)| {
                let pred = emg_pdf(t, params[0], params[1], params[2], params[3]);
                (obs - pred).powi(2)
            })
            .sum()
    }

    /// Clamp parameters to valid ranges
    fn clamp_params(params: &mut [f64; 4]) {
        params[1] = params[1].max(1e-6); // sigma > 0
        params[2] = params[2].max(1e-6); // tau > 0
        params[3] = params[3].max(0.0); // amplitude >= 0
    }

    /// Solve a 4x4 linear system Ax=b via Gaussian elimination with partial pivoting.
    /// Modifies a and b in place. Returns None if singular.
    fn solve_4x4(a: &mut [f64; 16], b: &mut [f64; 4]) -> Option<[f64; 4]> {
        for col in 0..4 {
            // Partial pivoting
            let mut max_val = a[col * 4 + col].abs();
            let mut max_row = col;
            for row in (col + 1)..4 {
                let val = a[row * 4 + col].abs();
                if val > max_val {
                    max_val = val;
                    max_row = row;
                }
            }
            if max_val < 1e-30 {
                return None;
            }
            if max_row != col {
                for k in 0..4 {
                    a.swap(col * 4 + k, max_row * 4 + k);
                }
                b.swap(col, max_row);
            }
            // Eliminate below
            for row in (col + 1)..4 {
                let factor = a[row * 4 + col] / a[col * 4 + col];
                for k in col..4 {
                    a[row * 4 + k] -= factor * a[col * 4 + k];
                }
                b[row] -= factor * b[col];
            }
        }
        // Back-substitute
        let mut x = [0.0; 4];
        for col in (0..4).rev() {
            let mut sum = b[col];
            for k in (col + 1)..4 {
                sum -= a[col * 4 + k] * x[k];
            }
            if a[col * 4 + col].abs() < 1e-30 {
                return None;
            }
            x[col] = sum / a[col * 4 + col];
        }
        Some(x)
    }

    /// Find the RT where emg_cdf = target_quantile using bisection.
    fn find_cdf_quantile(
        mu: f64,
        sigma: f64,
        tau: f64,
        target: f64,
        mut lo: f64,
        mut hi: f64,
    ) -> Option<f64> {
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            let cdf_val = emg_cdf(mid, mu, sigma, tau);
            if (cdf_val - target).abs() < 1e-8 {
                return Some(mid);
            }
            if cdf_val < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some((lo + hi) / 2.0)
    }
}

impl Default for EmgFitter {
    fn default() -> Self {
        Self::new()
    }
}

/// EMG (Exponentially Modified Gaussian) parameters.
///
/// Models an LC chromatographic peak as the convolution of a Gaussian (μ, σ)
/// with an exponential decay (τ). For right-tailing peaks (typical in LC-MS),
/// τ > 0 shifts the peak rightward and adds an exponential tail.
#[derive(Debug, Clone, Copy, Default)]
pub struct EmgParameters {
    /// Gaussian center (mu) — RT of the underlying Gaussian component
    pub mu: f64,
    /// Gaussian width (sigma) — standard deviation of the Gaussian component
    pub sigma: f64,
    /// Exponential decay constant (tau) — tailing factor
    pub tau: f64,
    /// Peak amplitude
    pub amplitude: f64,
    /// Fit quality (R²) — fraction of variance explained
    pub r_squared: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies peak detection finds a single Gaussian-like peak with the correct apex near RT 10.
    #[test]
    fn test_peak_detector() {
        let detector = PeakDetector::new().with_min_height(0.1);

        // Create a simple Gaussian-like peak
        let series: Vec<(f64, f64)> = (0..20)
            .map(|i| {
                let rt = i as f64;
                let coef = (-((rt - 10.0).powi(2) / 8.0)).exp();
                (rt, coef)
            })
            .collect();

        let peaks = detector.detect(&series);
        assert_eq!(peaks.len(), 1);
        assert!((peaks[0].apex_rt - 10.0).abs() < 1.0);
    }

    /// Verifies find_best_peak selects the highest-coefficient peak within the given RT tolerance window.
    #[test]
    fn test_find_best_peak() {
        let detector = PeakDetector::new();
        let peaks = vec![
            PeakBoundaries {
                start_rt: 5.0,
                end_rt: 8.0,
                apex_rt: 6.5,
                apex_coefficient: 0.5,
                integrated_area: 1.0,
                peak_quality: PeakQuality::default(),
            },
            PeakBoundaries {
                start_rt: 10.0,
                end_rt: 15.0,
                apex_rt: 12.0,
                apex_coefficient: 1.0,
                integrated_area: 3.0,
                peak_quality: PeakQuality::default(),
            },
        ];

        let best = detector.find_best_peak(&peaks, 11.0, 2.0);
        assert!(best.is_some());
        assert!((best.unwrap().apex_rt - 12.0).abs() < 0.1);
    }

    /// Helper: create a Gaussian XIC with given parameters
    fn make_gaussian_xic(center: f64, sigma: f64, amplitude: f64, n: usize) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let rt = i as f64 * 0.1; // 0.1 min spacing (6 sec cycle time)
                let intensity =
                    amplitude * (-((rt - center).powi(2)) / (2.0 * sigma * sigma)).exp();
                (rt, intensity)
            })
            .collect()
    }

    /// Helper: create a tailing peak (Gaussian rise, exponential decay tail).
    /// This models real chromatographic tailing where the trailing edge decays
    /// much slower than Gaussian, causing valley detection to extend too far.
    fn make_tailing_xic(
        center: f64,
        sigma_left: f64,
        tau: f64,
        amplitude: f64,
        n: usize,
    ) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let rt = i as f64 * 0.1;
                let intensity = if rt <= center {
                    // Gaussian leading edge
                    amplitude * (-((rt - center).powi(2)) / (2.0 * sigma_left * sigma_left)).exp()
                } else {
                    // Exponential decay trailing edge
                    amplitude * (-(rt - center) / tau).exp()
                };
                (rt, intensity)
            })
            .collect()
    }

    /// Symmetric Gaussian peak: FWHM capping should not significantly change boundaries
    #[test]
    fn test_fwhm_cap_symmetric_peak() {
        let xic = make_gaussian_xic(5.0, 0.3, 1000.0, 100);
        let peak = detect_xic_peak(&xic, 0.01, 5.0, Some(5.0)).unwrap();

        // For a Gaussian with σ=0.3, FWHM = 2.355 * 0.3 = 0.707 min
        // Valley detection at 20% stops around ±1.2σ = ±0.36 min
        // FWHM cap at 2.0 * HWHM = 2.0 * 0.353 = 0.706 min from apex
        // The cap is wider than the valley boundary, so no tightening expected
        assert!((peak.apex_rt - 5.0).abs() < 0.15);
        assert!(peak.end_rt - peak.start_rt > 0.4); // Peak is at least 0.4 min wide
        assert!(peak.end_rt - peak.start_rt < 2.0); // But not excessively wide
    }

    /// Tailing peak: FWHM capping should tighten the right boundary.
    /// Uses Gaussian rise + exponential decay to model chromatographic tailing.
    /// Exponential tails decay slower than Gaussian, so valley detection (20% threshold)
    /// extends too far, but FWHM capping (2× half-width) is tighter.
    #[test]
    fn test_fwhm_cap_tailing_peak() {
        // Gaussian leading edge (σ=0.3) + exponential trailing edge (τ=0.8)
        // The exponential tail is much slower than Gaussian decay
        let xic = make_tailing_xic(5.0, 0.3, 0.8, 1000.0, 150);

        // Get boundaries without FWHM capping (using raw valley detection)
        let smoothed: Vec<f64> =
            smooth_savitzky_golay(&xic.iter().map(|(_, v)| *v).collect::<Vec<_>>());
        let apex_idx = find_apex(&smoothed, 0.01, Some(5.0), &xic).unwrap();
        let apex_intensity = smoothed[apex_idx];
        let threshold = apex_intensity / 5.0;
        let valley_end = walk_boundary_right(&smoothed, apex_idx, apex_intensity, threshold);
        let valley_end_rt = xic[valley_end].0;

        // Get boundaries with FWHM capping
        let peak = detect_xic_peak(&xic, 0.01, 5.0, Some(5.0)).unwrap();

        // The capped right boundary should be tighter than valley detection
        // For exponential tail: valley at ln(5)*τ = 1.29 min, cap at 2*ln(2)*τ = 1.11 min
        assert!(
            peak.end_rt < valley_end_rt,
            "FWHM cap should tighten trailing edge: capped={:.3} vs valley={:.3}",
            peak.end_rt,
            valley_end_rt
        );

        // Left boundary should still be close to the peak (σ=0.3 is narrow)
        assert!(peak.apex_rt - peak.start_rt < 1.5);
    }

    /// Adjacent peaks with valley: valley boundary should be preserved (not widened)
    #[test]
    fn test_fwhm_cap_preserves_valley() {
        // Two peaks: one at RT 3.0, one at RT 5.0, with a valley between them
        let xic: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                let rt = i as f64 * 0.1;
                let peak1 = 1000.0 * (-((rt - 3.0).powi(2)) / (2.0 * 0.3 * 0.3)).exp();
                let peak2 = 800.0 * (-((rt - 5.0).powi(2)) / (2.0 * 0.3 * 0.3)).exp();
                (rt, peak1 + peak2)
            })
            .collect();

        // Detect the first peak
        let peak = detect_xic_peak(&xic, 0.01, 5.0, Some(3.0)).unwrap();

        // The right boundary should stop before the second peak's apex
        assert!(
            peak.end_rt < 5.0,
            "Boundary should not extend into second peak: end_rt={:.3}",
            peak.end_rt
        );
    }

    /// FWHM computation fails gracefully when signal never drops below 50%
    #[test]
    fn test_fwhm_cap_fallback() {
        // Very short XIC where half-height crossing may not be found
        let xic = vec![
            (0.0, 100.0),
            (0.1, 500.0),
            (0.2, 1000.0),
            (0.3, 600.0),
            (0.4, 100.0),
        ];
        let result = detect_xic_peak(&xic, 0.01, 5.0, None);
        // Should still detect a peak (falls back to valley detection)
        assert!(result.is_some());
    }

    /// Verify asymmetric half-widths are computed correctly for a tailing peak.
    /// Uses Gaussian rise + exponential decay to model chromatographic tailing.
    #[test]
    fn test_asymmetric_half_widths() {
        // Gaussian leading edge (σ=0.3) + exponential trailing edge (τ=0.8)
        let xic = make_tailing_xic(5.0, 0.3, 0.8, 1000.0, 150);
        let smoothed: Vec<f64> =
            smooth_savitzky_golay(&xic.iter().map(|(_, v)| *v).collect::<Vec<_>>());
        let apex_idx = find_apex(&smoothed, 0.01, Some(5.0), &xic).unwrap();

        let (left_hw, right_hw) =
            compute_asymmetric_half_widths(&smoothed, &xic, apex_idx).unwrap();

        // Right half-width should be larger than left (tailing peak)
        assert!(
            right_hw > left_hw,
            "Right half-width ({:.4}) should exceed left ({:.4}) for tailing peak",
            right_hw,
            left_hw
        );

        // Left HWHM ≈ σ * √(2ln2) = 0.3 * 1.177 = 0.353 (broadened by smoothing)
        // Right HWHM ≈ τ * ln(2) = 0.8 * 0.693 = 0.554
        assert!(left_hw > 0.2, "Left half-width too small: {:.4}", left_hw);
        assert!(left_hw < 0.6, "Left half-width too large: {:.4}", left_hw);
        assert!(
            right_hw > 0.4,
            "Right half-width too small: {:.4}",
            right_hw
        );
        assert!(
            right_hw < 1.0,
            "Right half-width too large: {:.4}",
            right_hw
        );
    }

    /// detect_all_xic_peaks returns multiple candidate peaks sorted by intensity.
    #[test]
    fn test_detect_all_xic_peaks_finds_multiple() {
        // Two well-separated peaks: large at RT 5.8, small at RT 3.0
        let xic: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                let rt = i as f64 * 0.1;
                let peak1 = 8000.0 * (-((rt - 5.8).powi(2)) / (2.0 * 0.2 * 0.2)).exp();
                let peak2 = 1000.0 * (-((rt - 3.0).powi(2)) / (2.0 * 0.15 * 0.15)).exp();
                (rt, peak1 + peak2)
            })
            .collect();

        let peaks = detect_all_xic_peaks(&xic, 0.01, 5.0);
        assert!(
            peaks.len() >= 2,
            "Should find at least 2 peaks, got {}",
            peaks.len()
        );
        // First peak should be the tallest (sorted by intensity descending)
        assert!(
            (peaks[0].apex_rt - 5.8).abs() < 0.2,
            "First peak should be at 5.8, got {:.2}",
            peaks[0].apex_rt
        );
        assert!(
            peaks[0].apex_intensity > peaks[1].apex_intensity,
            "Peaks should be sorted by intensity descending"
        );
    }

    /// detect_all_xic_peaks returns a single peak for a clean Gaussian.
    #[test]
    fn test_detect_all_xic_peaks_single_peak() {
        let xic = make_gaussian_xic(5.0, 0.3, 1000.0, 100);
        let peaks = detect_all_xic_peaks(&xic, 0.01, 5.0);
        assert_eq!(peaks.len(), 1);
        assert!((peaks[0].apex_rt - 5.0).abs() < 0.2);
    }

    /// detect_all_xic_peaks with noise bump + real peak — caller can pick best by score.
    /// This models the real-world File 49 case where the old code picked noise near
    /// expected_rt over the real 8x-taller peak 0.3 min further away.
    #[test]
    fn test_detect_all_xic_peaks_noise_plus_real() {
        let xic: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                let rt = i as f64 * 0.1;
                let noise = 1000.0 * (-((rt - 5.4).powi(2)) / (2.0 * 0.15 * 0.15)).exp();
                let real_peak = 8000.0 * (-((rt - 5.8).powi(2)) / (2.0 * 0.2 * 0.2)).exp();
                (rt, noise + real_peak)
            })
            .collect();

        let peaks = detect_all_xic_peaks(&xic, 0.01, 5.0);
        assert!(!peaks.is_empty());
        // The tallest candidate (first) should be the real peak
        assert!(
            (peaks[0].apex_rt - 5.8).abs() < 0.2,
            "Tallest candidate should be real peak at 5.8, got {:.2}",
            peaks[0].apex_rt
        );
    }

    // =========================================================================
    // EMG fitting tests
    // =========================================================================

    /// erfc_approx matches known values
    #[test]
    fn test_erfc_approx() {
        assert!((erfc_approx(0.0) - 1.0).abs() < 1e-6);
        assert!((erfc_approx(1.0) - 0.157299).abs() < 1e-4);
        assert!((erfc_approx(-1.0) - 1.842701).abs() < 1e-4);
        assert!(erfc_approx(5.0) < 1e-10);
        assert!((erfc_approx(-5.0) - 2.0).abs() < 1e-10);
    }

    /// normal_cdf matches known values
    #[test]
    fn test_normal_cdf() {
        assert!((normal_cdf(0.0) - 0.5).abs() < 1e-6);
        assert!((normal_cdf(1.96) - 0.975).abs() < 1e-3);
        assert!((normal_cdf(-1.96) - 0.025).abs() < 1e-3);
    }

    /// EMG PDF produces a right-tailing peak and integrates to ~A
    #[test]
    fn test_emg_pdf_basic() {
        let mu = 10.0;
        let sigma = 0.1;
        let tau = 0.05;
        let amplitude = 1000.0;

        // Peak should be near mu + tau
        let apex_approx = mu + tau;
        let val_at_apex = emg_pdf(apex_approx, mu, sigma, tau, amplitude);
        assert!(val_at_apex > 0.0, "EMG should be positive at apex");

        // Should be near-zero far from peak
        let far_left = emg_pdf(mu - 5.0 * sigma, mu, sigma, tau, amplitude);
        let far_right = emg_pdf(mu + 10.0 * tau, mu, sigma, tau, amplitude);
        assert!(far_left < val_at_apex * 0.01);
        assert!(far_right < val_at_apex * 0.01);
    }

    /// EMG CDF is monotonically increasing from 0 to 1
    #[test]
    fn test_emg_cdf_monotonic() {
        let mu = 10.0;
        let sigma = 0.1;
        let tau = 0.05;

        let mut prev = 0.0;
        for i in 0..200 {
            let t = 8.0 + i as f64 * 0.02;
            let cdf = emg_cdf(t, mu, sigma, tau);
            assert!(
                cdf >= prev - 1e-10,
                "CDF should be non-decreasing: at t={:.2}, cdf={:.6} < prev={:.6}",
                t,
                cdf,
                prev
            );
            prev = cdf;
        }
        // Should approach 1.0 far right
        assert!(emg_cdf(mu + 20.0 * sigma, mu, sigma, tau) > 0.99);
        // Should approach 0.0 far left
        assert!(emg_cdf(mu - 10.0 * sigma, mu, sigma, tau) < 0.01);
    }

    /// Fit EMG to a synthetic Gaussian peak (τ→0 limit)
    #[test]
    fn test_emg_fit_gaussian() {
        let true_mu = 10.0;
        let true_sigma = 0.15;
        let true_amp = 5000.0;

        let profile: Vec<(f64, f64)> = (0..200)
            .map(|i| {
                let rt = 9.0 + i as f64 * 0.01;
                let v =
                    true_amp * (-((rt - true_mu).powi(2)) / (2.0 * true_sigma * true_sigma)).exp();
                (rt, v)
            })
            .collect();

        let fitter = EmgFitter::new();
        let params = fitter.fit(&profile).unwrap();

        assert!(
            params.r_squared > 0.95,
            "R² should be high for clean Gaussian: {:.4}",
            params.r_squared
        );
        assert!(
            (params.mu - true_mu).abs() < 0.1,
            "μ should be near {}: got {:.4}",
            true_mu,
            params.mu
        );
        assert!(
            (params.sigma - true_sigma).abs() < 0.1,
            "σ should be near {}: got {:.4}",
            true_sigma,
            params.sigma
        );
    }

    /// Fit EMG to a synthetic EMG peak and recover parameters
    #[test]
    fn test_emg_fit_tailing_peak() {
        let true_mu = 10.0;
        let true_sigma = 0.08;
        let true_tau = 0.04;
        let true_amp = 5000.0;

        let profile: Vec<(f64, f64)> = (0..300)
            .map(|i| {
                let rt = 9.5 + i as f64 * 0.005;
                (rt, emg_pdf(rt, true_mu, true_sigma, true_tau, true_amp))
            })
            .collect();

        let fitter = EmgFitter::new();
        let params = fitter.fit(&profile).unwrap();

        assert!(
            params.r_squared > 0.99,
            "R² should be very high for exact EMG data: {:.4}",
            params.r_squared
        );
        assert!(
            (params.mu - true_mu).abs() < 0.02,
            "μ: expected {}, got {:.4}",
            true_mu,
            params.mu
        );
        assert!(
            (params.sigma - true_sigma).abs() < 0.02,
            "σ: expected {}, got {:.4}",
            true_sigma,
            params.sigma
        );
        assert!(
            (params.tau - true_tau).abs() < 0.02,
            "τ: expected {}, got {:.4}",
            true_tau,
            params.tau
        );
    }

    /// EMG boundaries capture 95% of the peak area
    #[test]
    fn test_emg_boundaries() {
        let params = EmgParameters {
            mu: 10.0,
            sigma: 0.1,
            tau: 0.05,
            amplitude: 1000.0,
            r_squared: 0.99,
        };

        let (left, right) = EmgFitter::boundaries(&params, 0.95).unwrap();

        // Boundaries should bracket the peak
        assert!(left < params.mu, "Left boundary should be before μ");
        assert!(
            right > params.mu + params.tau,
            "Right boundary should be after μ+τ"
        );

        // Width should be reasonable (not too narrow, not too wide)
        let width = right - left;
        assert!(width > 0.2, "95% boundaries too narrow: {:.4} min", width);
        assert!(width < 2.0, "95% boundaries too wide: {:.4} min", width);

        // Verify CDF at boundaries
        let cdf_left = emg_cdf(left, params.mu, params.sigma, params.tau);
        let cdf_right = emg_cdf(right, params.mu, params.sigma, params.tau);
        assert!(
            (cdf_left - 0.025).abs() < 0.01,
            "CDF at left boundary: {:.4}, expected 0.025",
            cdf_left
        );
        assert!(
            (cdf_right - 0.975).abs() < 0.01,
            "CDF at right boundary: {:.4}, expected 0.975",
            cdf_right
        );
    }

    /// Full pipeline: fit EMG to synthetic peak and verify boundaries are tight
    #[test]
    fn test_emg_fit_and_boundaries() {
        let true_mu = 13.8;
        let true_sigma = 0.06;
        let true_tau = 0.03;
        let true_amp = 100_000.0;

        // Simulate a typical Astral peak: ~0.15 min wide, 0.01 min spacing
        let profile: Vec<(f64, f64)> = (0..400)
            .map(|i| {
                let rt = 12.0 + i as f64 * 0.01;
                (rt, emg_pdf(rt, true_mu, true_sigma, true_tau, true_amp))
            })
            .collect();

        let fitter = EmgFitter::new();
        let params = fitter.fit(&profile).unwrap();
        let (left, right) = EmgFitter::boundaries(&params, 0.95).unwrap();
        let width = right - left;

        // For σ=0.06, τ=0.03, expect ~0.3-0.5 min boundaries (not 3+ min!)
        assert!(
            width < 1.0,
            "EMG boundaries should be <1 min for a narrow peak: got {:.3} min",
            width
        );
        assert!(
            width > 0.1,
            "EMG boundaries should be >0.1 min: got {:.3} min",
            width
        );

        // Center should be near the peak
        let center = (left + right) / 2.0;
        assert!(
            (center - true_mu).abs() < 0.2,
            "Boundary center should be near μ={}: got {:.3}",
            true_mu,
            center
        );
    }
}
