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
/// Scans for local maxima above `min_height`. If `expected_rt` is provided,
/// returns the qualifying apex nearest to it. Otherwise returns the global maximum.
fn find_apex(
    smoothed: &[f64],
    min_height: f64,
    expected_rt: Option<f64>,
    xic: &[(f64, f64)],
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

    if let Some(expected) = expected_rt {
        // Pick the local maximum nearest to expected RT
        candidates
            .iter()
            .min_by(|a, b| {
                let dist_a = (xic[a.0].0 - expected).abs();
                let dist_b = (xic[b.0].0 - expected).abs();
                dist_a.total_cmp(&dist_b)
            })
            .map(|&(idx, _)| idx)
    } else {
        // Pick global maximum
        candidates
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .map(|&(idx, _)| idx)
    }
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
/// This is the canonical SNR function used by coelution, regression, and
/// calibration scoring.
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

/// Placeholder for EMG fitting
/// TODO: Implement in Phase 2
pub struct EmgFitter;

impl EmgFitter {
    pub fn new() -> Self {
        Self
    }

    /// Fit EMG to coefficient series (placeholder)
    pub fn fit(&self, _series: &[(f64, f64)]) -> Result<EmgParameters> {
        // TODO: Implement Levenberg-Marquardt fitting
        Ok(EmgParameters::default())
    }
}

impl Default for EmgFitter {
    fn default() -> Self {
        Self::new()
    }
}

/// EMG (Exponentially Modified Gaussian) parameters
#[derive(Debug, Clone, Copy, Default)]
pub struct EmgParameters {
    /// Center (mu)
    pub mu: f64,
    /// Width (sigma)
    pub sigma: f64,
    /// Tailing factor (tau)
    pub tau: f64,
    /// Amplitude
    pub amplitude: f64,
    /// Fit quality (R²)
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
}
