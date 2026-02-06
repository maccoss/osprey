//! Mass and retention time calibration for DIA searches
//!
//! This module implements auto-calibration to correct systematic errors in:
//! - m/z measurements (MS1 and MS2)
//! - Retention time predictions
//!
//! The calibration workflow:
//! 1. Sample high-quality library peptides
//! 2. Run calibration search with target-decoy competition
//! 3. Filter to 1% FDR
//! 4. Calculate calibration parameters from confident matches
//! 5. Apply corrections during main search

pub mod mass;
pub mod rt;
pub mod io;

use serde::{Deserialize, Serialize};

/// Complete calibration parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationParams {
    /// Metadata about the calibration process
    pub metadata: CalibrationMetadata,
    /// MS1 (precursor) m/z calibration
    pub ms1_calibration: MzCalibration,
    /// MS2 (fragment) m/z calibration
    pub ms2_calibration: MzCalibration,
    /// Retention time calibration
    pub rt_calibration: RTCalibrationParams,
}

impl CalibrationParams {
    /// Create default uncalibrated parameters
    pub fn uncalibrated() -> Self {
        Self {
            metadata: CalibrationMetadata {
                num_confident_peptides: 0,
                num_sampled_precursors: 0,
                calibration_successful: false,
                timestamp: chrono::Utc::now().to_rfc3339(),
                isolation_scheme: None,
            },
            ms1_calibration: MzCalibration::uncalibrated(),
            ms2_calibration: MzCalibration::uncalibrated(),
            rt_calibration: RTCalibrationParams::uncalibrated(),
        }
    }

    /// Check if calibration was successful
    pub fn is_calibrated(&self) -> bool {
        self.metadata.calibration_successful
    }

    /// Log calibration summary
    pub fn log_summary(&self) {
        log::info!("=== Calibration Summary ===");
        log::info!(
            "Confident peptides: {} (from {} sampled)",
            self.metadata.num_confident_peptides,
            self.metadata.num_sampled_precursors
        );
        log::info!(
            "MS1 calibration: mean={:.2} {}, median={:.2} {}, SD={:.2} {} (n={})",
            self.ms1_calibration.mean,
            self.ms1_calibration.unit,
            self.ms1_calibration.median,
            self.ms1_calibration.unit,
            self.ms1_calibration.sd,
            self.ms1_calibration.unit,
            self.ms1_calibration.count
        );
        log::info!(
            "MS2 calibration: mean={:.4} {}, median={:.4} {}, SD={:.4} {} (n={})",
            self.ms2_calibration.mean,
            self.ms2_calibration.unit,
            self.ms2_calibration.median,
            self.ms2_calibration.unit,
            self.ms2_calibration.sd,
            self.ms2_calibration.unit,
            self.ms2_calibration.count
        );
        log::info!(
            "RT calibration: method={:?}, residual_SD={:.2} min, R²={:.4}",
            self.rt_calibration.method,
            self.rt_calibration.residual_sd,
            self.rt_calibration.r_squared
        );

        // Display ASCII histograms if available
        if self.ms1_calibration.calibrated {
            self.ms1_calibration.log_histogram("MS1 (precursor)");
        }
        if self.ms2_calibration.calibrated {
            self.ms2_calibration.log_histogram("MS2 (fragment)");
        }
    }
}

/// Calibration metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibrationMetadata {
    /// Number of confident peptides used for calibration (after FDR filtering)
    pub num_confident_peptides: usize,
    /// Number of precursors sampled for calibration discovery
    pub num_sampled_precursors: usize,
    /// Whether calibration was successful
    pub calibration_successful: bool,
    /// Timestamp when calibration was performed
    pub timestamp: String,
    /// DIA isolation window scheme (from first cycle of MS2 spectra)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub isolation_scheme: Option<IsolationScheme>,
}

/// DIA isolation window scheme extracted from mzML
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsolationScheme {
    /// Number of isolation windows per cycle
    pub num_windows: usize,
    /// Minimum isolation window center m/z
    pub mz_min: f64,
    /// Maximum isolation window center m/z
    pub mz_max: f64,
    /// Typical isolation window width (Da) - may vary across range
    pub typical_width: f64,
    /// Whether all windows have the same width
    pub uniform_width: bool,
    /// Individual windows: (center, width) for each window in cycle
    pub windows: Vec<(f64, f64)>,
}

/// m/z calibration parameters
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MzCalibration {
    /// Mean PPM error (systematic error)
    pub mean: f64,

    /// Median PPM error
    pub median: f64,

    /// Standard deviation of PPM errors
    pub sd: f64,

    /// Number of observations
    pub count: usize,

    /// Unit (always "ppm")
    pub unit: String,

    /// Adjusted tolerance: |mean| + 3*sd
    #[serde(skip_serializing_if = "Option::is_none")]
    pub adjusted_tolerance: Option<f64>,

    /// Window halfwidth multiplier (default: 3.0)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub window_halfwidth_multiplier: Option<f64>,

    /// Histogram of PPM errors (bin edges and counts)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub histogram: Option<MzHistogram>,

    /// Whether calibration was successfully performed
    pub calibrated: bool,
}

/// Histogram data for m/z error distribution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MzHistogram {
    /// Bin edges (N+1 edges for N bins)
    pub bin_edges: Vec<f64>,
    /// Counts in each bin
    pub counts: Vec<usize>,
    /// Bin width in ppm
    pub bin_width: f64,
}

impl MzCalibration {
    /// Create uncalibrated parameters
    pub fn uncalibrated() -> Self {
        Self {
            mean: 0.0,
            median: 0.0,
            sd: 0.0,
            count: 0,
            unit: "ppm".to_string(),
            adjusted_tolerance: None,
            window_halfwidth_multiplier: None,
            histogram: None,
            calibrated: false,
        }
    }

    /// Get the effective tolerance to use for matching
    ///
    /// Returns the calibrated tolerance if available, otherwise the base tolerance
    pub fn effective_tolerance(&self, base_tolerance_ppm: f64) -> f64 {
        if self.calibrated {
            self.adjusted_tolerance.unwrap_or(base_tolerance_ppm)
        } else {
            base_tolerance_ppm
        }
    }

    /// Log ASCII histogram of mass errors
    pub fn log_histogram(&self, label: &str) {
        if let Some(ref hist) = self.histogram {
            log::info!("--- {} Mass Error Histogram ({}) ---", label, self.unit);

            // Find max count for scaling
            let max_count = hist.counts.iter().max().copied().unwrap_or(1);
            let max_bar_width = 40;

            // Use appropriate format precision based on unit
            let is_th = self.unit == "Th";

            for (i, &count) in hist.counts.iter().enumerate() {
                let bin_start = hist.bin_edges[i];
                let bin_end = hist.bin_edges[i + 1];
                let bin_center = (bin_start + bin_end) / 2.0;

                // Scale bar width
                let bar_width = if max_count > 0 {
                    (count as f64 / max_count as f64 * max_bar_width as f64).round() as usize
                } else {
                    0
                };

                let bar: String = "█".repeat(bar_width);
                if is_th {
                    log::info!(
                        "{:>7.3} | {:6} | {}",
                        bin_center,
                        count,
                        bar
                    );
                } else {
                    log::info!(
                        "{:>7.1} | {:6} | {}",
                        bin_center,
                        count,
                        bar
                    );
                }
            }
            if is_th {
                log::info!(
                    "        | mean={:.4}, median={:.4}, SD={:.4} {} (n={})",
                    self.mean,
                    self.median,
                    self.sd,
                    self.unit,
                    self.count
                );
            } else {
                log::info!(
                    "        | mean={:.2}, median={:.2}, SD={:.2} {} (n={})",
                    self.mean,
                    self.median,
                    self.sd,
                    self.unit,
                    self.count
                );
            }
        }
    }
}

/// RT calibration parameters (serializable summary)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RTCalibrationParams {
    /// Calibration method used
    pub method: RTCalibrationMethod,
    /// Residual standard deviation (for RT window sizing)
    pub residual_sd: f64,
    /// Number of calibration points
    pub n_points: usize,
    /// R² (coefficient of determination)
    pub r_squared: f64,
    /// LOESS model parameters for reconstruction (optional, for reuse)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub model_params: Option<RTModelParams>,
}

/// LOESS model parameters for serialization
///
/// Stores the library RTs and fitted measured RTs which allows
/// the calibration curve to be reconstructed via interpolation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RTModelParams {
    /// Library retention times (sorted, used as x-values for interpolation)
    pub library_rts: Vec<f64>,
    /// Fitted measured retention times (corresponding to library_rts)
    pub fitted_rts: Vec<f64>,
    /// Absolute residuals at each calibration point (for local RT tolerance)
    /// This field is optional for backwards compatibility with older calibration files
    #[serde(default)]
    pub abs_residuals: Vec<f64>,
}

impl RTCalibrationParams {
    /// Create uncalibrated parameters
    pub fn uncalibrated() -> Self {
        Self {
            method: RTCalibrationMethod::None,
            residual_sd: 0.0,
            n_points: 0,
            r_squared: 0.0,
            model_params: None,
        }
    }

    /// Check if RT was calibrated
    pub fn is_calibrated(&self) -> bool {
        self.method != RTCalibrationMethod::None && self.n_points > 0
    }

    /// Check if this calibration has model data for reconstruction
    pub fn has_model_data(&self) -> bool {
        self.model_params.as_ref()
            .map(|m| !m.library_rts.is_empty() && m.library_rts.len() == m.fitted_rts.len())
            .unwrap_or(false)
    }
}

/// RT calibration method
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RTCalibrationMethod {
    /// No calibration performed
    None,
    /// LOESS (locally weighted regression)
    LOESS,
    /// Linear regression (fallback)
    Linear,
}

// Re-export key types and functions from submodules
pub use mass::{calculate_mz_calibration, apply_mz_calibration, apply_spectrum_calibration, calibrated_tolerance_ppm, calibrated_tolerance, calculate_ppm_error, MzQCData};
pub use rt::{RTCalibration, RTCalibrationStats, RTCalibrator, RTCalibratorConfig, RTStratifiedSampler};
pub use io::{save_calibration, load_calibration, calibration_filename, calibration_filename_for_input, calibration_path_for_input};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uncalibrated_params() {
        let params = CalibrationParams::uncalibrated();
        assert!(!params.is_calibrated());
        assert!(!params.ms1_calibration.calibrated);
        assert!(!params.ms2_calibration.calibrated);
        assert!(!params.rt_calibration.is_calibrated());
    }

    #[test]
    fn test_effective_tolerance() {
        let uncalibrated = MzCalibration::uncalibrated();
        assert_eq!(uncalibrated.effective_tolerance(10.0), 10.0);

        let calibrated = MzCalibration {
            mean: -2.5,
            median: -2.4,
            sd: 1.0,
            count: 100,
            unit: "ppm".to_string(),
            adjusted_tolerance: Some(5.5),
            window_halfwidth_multiplier: Some(3.0),
            histogram: None,
            calibrated: true,
        };
        assert_eq!(calibrated.effective_tolerance(10.0), 5.5);
    }

    #[test]
    fn test_serialization() {
        let params = CalibrationParams {
            metadata: CalibrationMetadata {
                num_confident_peptides: 100,
                num_sampled_precursors: 1000,
                calibration_successful: true,
                timestamp: "2024-01-15T10:30:00Z".to_string(),
                isolation_scheme: None,
            },
            ms1_calibration: MzCalibration {
                mean: -2.5,
                median: -2.4,
                sd: 0.8,
                count: 100,
                unit: "ppm".to_string(),
                adjusted_tolerance: Some(4.9),
                window_halfwidth_multiplier: Some(3.0),
                histogram: Some(MzHistogram {
                    bin_edges: vec![-5.0, -4.0, -3.0, -2.0, -1.0, 0.0],
                    counts: vec![5, 15, 50, 25, 5],
                    bin_width: 1.0,
                }),
                calibrated: true,
            },
            ms2_calibration: MzCalibration {
                mean: 1.2,
                median: 1.1,
                sd: 1.0,
                count: 500,
                unit: "ppm".to_string(),
                adjusted_tolerance: Some(4.2),
                window_halfwidth_multiplier: Some(3.0),
                histogram: None,
                calibrated: true,
            },
            rt_calibration: RTCalibrationParams {
                method: RTCalibrationMethod::LOESS,
                residual_sd: 0.8,
                n_points: 100,
                r_squared: 0.98,
                model_params: Some(RTModelParams {
                    library_rts: vec![0.0, 10.0, 20.0, 30.0],
                    fitted_rts: vec![1.0, 11.0, 21.0, 31.0],
                    abs_residuals: vec![0.1, 0.2, 0.3, 0.4],
                }),
            },
        };

        // Serialize to JSON
        let json = serde_json::to_string_pretty(&params).unwrap();
        assert!(json.contains("\"ms1_calibration\""));
        assert!(json.contains("\"mean\": -2.5"));
        assert!(json.contains("\"median\": -2.4"));
        assert!(json.contains("\"histogram\""));

        // Deserialize back
        let loaded: CalibrationParams = serde_json::from_str(&json).unwrap();
        assert!(loaded.is_calibrated());
        assert!((loaded.ms1_calibration.mean + 2.5).abs() < 0.001);
        assert!(loaded.ms1_calibration.histogram.is_some());
    }
}
