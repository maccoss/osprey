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

/// Smooth a time series with weighted moving average [0.25, 0.5, 0.25].
///
/// Endpoints use asymmetric weights: [2/3, 1/3] for first and [1/3, 2/3] for last.
/// This matches DIA-NN's smoothing kernel.
pub fn smooth_weighted_avg(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![values[0]];
    }
    let mut out = vec![0.0; n];
    out[0] = (2.0 / 3.0) * values[0] + (1.0 / 3.0) * values[1];
    out[n - 1] = (2.0 / 3.0) * values[n - 1] + (1.0 / 3.0) * values[n - 2];
    for i in 1..n - 1 {
        out[i] = 0.25 * values[i - 1] + 0.5 * values[i] + 0.25 * values[i + 1];
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
///    - Valley detection: stop at local minimum that is <33% of apex AND <50% of neighbor
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
    let smoothed = smooth_weighted_avg(&raw_intensities);

    // Find apex: local maximum above min_height
    let apex_idx = find_apex(&smoothed, min_height, expected_rt, xic)?;
    let apex_intensity = smoothed[apex_idx];

    // Walk left to find start boundary
    let boundary_threshold = apex_intensity / peak_boundary;
    let start_idx = walk_boundary_left(&smoothed, apex_idx, apex_intensity, boundary_threshold);

    // Walk right to find end boundary
    let end_idx = walk_boundary_right(&smoothed, apex_idx, apex_intensity, boundary_threshold);

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
/// - Valley detected: local minimum is <33% of apex AND <50% of the next rising neighbor
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
        } else if valley < apex_intensity / 3.0 && valley < smoothed[i] / 2.0 {
            // Valley detected: intensity dropped to <33% of apex and is <50% of rising neighbor
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
        } else if valley < apex_intensity / 3.0 && valley < val / 2.0 {
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
}
