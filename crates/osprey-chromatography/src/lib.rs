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

use osprey_core::{PeakBoundaries, PeakQuality, Result};

// Re-export calibration types
pub use calibration::{
    apply_mz_calibration,
    calculate_mz_calibration,
    calculate_ppm_error,
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
fn trapezoidal_area(series: &[(f64, f64)]) -> f64 {
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
