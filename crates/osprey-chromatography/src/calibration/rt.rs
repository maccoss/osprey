//! RT Calibration using LOESS (Local Polynomial Regression)
//!
//! This module provides retention time calibration to convert library RTs
//! (which may be iRT, predicted, or uncalibrated) to measured RTs.
//!
//! The calibration uses LOESS (Locally Estimated Scatterplot Smoothing)
//! which fits local polynomial regressions weighted by distance.

use super::RTModelParams;
use osprey_core::{OspreyError, Result};

/// RT calibration configuration
#[derive(Debug, Clone)]
pub struct RTCalibratorConfig {
    /// Bandwidth parameter for LOESS (fraction of data to use for each local fit)
    /// Typical values: 0.2-0.5
    pub bandwidth: f64,
    /// Polynomial degree for local fits (1 = linear, 2 = quadratic)
    pub degree: usize,
    /// Minimum number of points required for calibration
    pub min_points: usize,
    /// Number of robustness iterations (0 = no robustness weighting)
    pub robustness_iter: usize,
    /// Fraction of points to retain after outlier removal (0.8 = keep best 80%).
    /// After LOESS fitting, points with the largest absolute residuals are removed
    /// and the fit is repeated on the clean data (mirrors DIA-NN's map_RT approach).
    /// Set to 1.0 to disable outlier removal.
    pub outlier_retention: f64,
    /// When false (default), `residuals` is captured from the initial fit BEFORE
    /// the robustness loop and reused across all iterations. This produces a
    /// single refinement regardless of `robustness_iter >= 1` — the historical
    /// behavior of this crate. When true, each iteration recomputes residuals
    /// from the current fit, which is the classical Cleveland (1979) robust
    /// LOESS algorithm. Default false preserves backward compatibility; the
    /// `OSPREY_LOESS_CLASSICAL_ROBUST=1` env var flips to classical mode for
    /// cross-implementation validation (OspreySharp honors the same env var).
    pub classical_robust_iterations: bool,
}

impl Default for RTCalibratorConfig {
    fn default() -> Self {
        Self {
            bandwidth: 0.3,
            degree: 1, // Linear local fits
            min_points: 20,
            robustness_iter: 2,
            outlier_retention: 0.8, // Keep best 80%, remove worst 20% (DIA-NN default)
            classical_robust_iterations: false,
        }
    }
}

/// RT Calibrator using LOESS regression
#[derive(Debug, Clone)]
pub struct RTCalibrator {
    config: RTCalibratorConfig,
}

impl RTCalibrator {
    /// Create a new RT calibrator with default settings
    pub fn new() -> Self {
        Self {
            config: RTCalibratorConfig::default(),
        }
    }

    /// Create a new RT calibrator with custom settings
    pub fn with_config(config: RTCalibratorConfig) -> Self {
        Self { config }
    }

    /// Set the bandwidth parameter
    pub fn with_bandwidth(mut self, bandwidth: f64) -> Self {
        self.config.bandwidth = bandwidth;
        self
    }

    /// Set the polynomial degree
    pub fn with_degree(mut self, degree: usize) -> Self {
        self.config.degree = degree;
        self
    }

    /// Set the outlier retention fraction (0.0-1.0). Set to 1.0 to disable outlier removal.
    pub fn with_outlier_retention(mut self, outlier_retention: f64) -> Self {
        self.config.outlier_retention = outlier_retention;
        self
    }

    /// Fit a calibration curve from (library_rt, measured_rt) pairs
    ///
    /// # Arguments
    /// * `library_rts` - Library retention times (iRT or predicted)
    /// * `measured_rts` - Corresponding measured retention times
    ///
    /// # Returns
    /// A fitted `RTCalibration` that can predict measured RT from library RT
    pub fn fit(&self, library_rts: &[f64], measured_rts: &[f64]) -> Result<RTCalibration> {
        if library_rts.len() != measured_rts.len() {
            return Err(OspreyError::ConfigError(format!(
                "RT arrays must have same length: {} vs {}",
                library_rts.len(),
                measured_rts.len()
            )));
        }

        if library_rts.len() < self.config.min_points {
            return Err(OspreyError::ConfigError(format!(
                "Need at least {} calibration points, got {}",
                self.config.min_points,
                library_rts.len()
            )));
        }

        // Sort by library RT
        let mut pairs: Vec<(f64, f64)> = library_rts
            .iter()
            .zip(measured_rts.iter())
            .map(|(&x, &y)| (x, y))
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

        let x: Vec<f64> = pairs.iter().map(|(x, _)| *x).collect();
        let y: Vec<f64> = pairs.iter().map(|(_, y)| *y).collect();

        // Compute initial LOESS fit
        let fitted = self.loess_fit(&x, &y, None)?;

        // Compute residuals from the INITIAL fit. In the default
        // (non-classical) mode these residuals are reused unchanged across
        // every iteration of the robustness loop — the historical behavior
        // of this crate, which produces a single refinement regardless of
        // `robustness_iter >= 1`. In classical mode (Cleveland 1979),
        // residuals are refreshed from the current fit at the top of each
        // iteration, producing the expected multi-pass robust LOESS.
        // OspreySharp mirrors both modes via the same
        // OSPREY_LOESS_CLASSICAL_ROBUST env var so the two implementations
        // stay bit-identical in either configuration.
        let mut abs_residuals: Vec<f64> = y
            .iter()
            .zip(fitted.iter())
            .map(|(obs, pred)| (obs - pred).abs())
            .collect();

        // Robustness iterations
        let mut weights: Vec<f64> = vec![1.0; x.len()];
        let mut final_fitted = fitted;

        for iter in 0..self.config.robustness_iter {
            if self.config.classical_robust_iterations && iter > 0 {
                // Refresh residuals from the current fit (classical Cleveland 1979)
                for (i, (yv, fv)) in y.iter().zip(final_fitted.iter()).enumerate() {
                    abs_residuals[i] = (yv - fv).abs();
                }
            }

            // Compute robust weights using bisquare function
            let median_abs_residual = median(&abs_residuals);
            let s = 6.0 * median_abs_residual; // MAD scale estimate

            if s > 1e-10 {
                weights = abs_residuals
                    .iter()
                    .map(|r| {
                        let u = r / s;
                        if u.abs() < 1.0 {
                            (1.0 - u * u).powi(2)
                        } else {
                            0.0
                        }
                    })
                    .collect();
            }

            // Refit with weights
            final_fitted = self.loess_fit(&x, &y, Some(&weights))?;
        }

        // Outlier removal: remove points with largest absolute residuals, then refit
        // This mirrors DIA-NN's map_RT() approach — fit, remove worst 20%, refit
        let (x, y, final_fitted) = if self.config.outlier_retention < 1.0 {
            let abs_resid: Vec<f64> = y
                .iter()
                .zip(final_fitted.iter())
                .map(|(obs, pred)| (obs - pred).abs())
                .collect();
            let threshold = percentile_value(&abs_resid, self.config.outlier_retention);

            let keep: Vec<bool> = abs_resid.iter().map(|&r| r <= threshold + 1e-12).collect();
            let n_keep = keep.iter().filter(|&&k| k).count();

            if n_keep >= self.config.min_points {
                let x_filtered: Vec<f64> = x
                    .iter()
                    .zip(&keep)
                    .filter(|(_, &k)| k)
                    .map(|(&v, _)| v)
                    .collect();
                let y_filtered: Vec<f64> = y
                    .iter()
                    .zip(&keep)
                    .filter(|(_, &k)| k)
                    .map(|(&v, _)| v)
                    .collect();

                log::info!(
                    "RT outlier removal: kept {}/{} points (threshold={:.3} min)",
                    n_keep,
                    x.len(),
                    threshold
                );

                // Refit LOESS on clean data (no robustness weights needed)
                let refitted = self.loess_fit(&x_filtered, &y_filtered, None)?;
                (x_filtered, y_filtered, refitted)
            } else {
                log::warn!(
                    "Outlier removal would leave too few points ({} < {}), skipping",
                    n_keep,
                    self.config.min_points
                );
                (x, y, final_fitted)
            }
        } else {
            (x, y, final_fitted)
        };

        // Compute residuals and absolute residuals on (possibly filtered) data
        let residuals: Vec<f64> = y
            .iter()
            .zip(final_fitted.iter())
            .map(|(obs, pred)| obs - pred)
            .collect();
        let abs_residuals: Vec<f64> = residuals.iter().map(|r| r.abs()).collect();
        let residual_std = std_dev(&residuals);

        Ok(RTCalibration {
            library_rts: x,
            measured_rts: y,
            fitted_values: final_fitted,
            abs_residuals,
            bandwidth: self.config.bandwidth,
            degree: self.config.degree,
            residual_std,
        })
    }

    /// Perform LOESS fit
    fn loess_fit(&self, x: &[f64], y: &[f64], weights: Option<&[f64]>) -> Result<Vec<f64>> {
        let n = x.len();
        let k = (self.config.bandwidth * n as f64).ceil() as usize;
        let k = k.max(self.config.degree + 2).min(n); // Ensure enough points for polynomial

        let mut fitted = Vec::with_capacity(n);

        for i in 0..n {
            let xi = x[i];

            // Find k nearest neighbors as contiguous range [lo, hi).
            // x is sorted, so the k nearest neighbors are always a contiguous subrange.
            // Two-pointer expansion from index i: O(k) instead of O(n log n) sort.
            let (lo, hi) = find_k_nearest_sorted(x, i, k);
            let max_dist = (x[lo] - xi).abs().max((x[hi - 1] - xi).abs());

            // Accumulate weighted sums in a single pass (zero allocation)
            let mut sum_w = 0.0f64;
            let mut sum_wx = 0.0f64;
            let mut sum_wy = 0.0f64;
            let mut sum_wxx = 0.0f64;
            let mut sum_wxy = 0.0f64;

            for j in lo..hi {
                let dist = (x[j] - xi).abs();
                let u = if max_dist > 1e-10 {
                    dist / max_dist
                } else {
                    0.0
                };
                let tw = if u < 1.0 {
                    (1.0 - u.powi(3)).powi(3)
                } else {
                    0.0
                };
                let w = if let Some(wts) = weights {
                    tw * wts[j]
                } else {
                    tw
                };

                sum_w += w;
                sum_wx += w * x[j];
                sum_wy += w * y[j];
                sum_wxx += w * x[j] * x[j];
                sum_wxy += w * x[j] * y[j];
            }

            // Fit local polynomial (degree 0 = weighted mean, else weighted linear)
            let yi = if self.config.degree == 0 {
                if sum_w > 1e-10 {
                    sum_wy / sum_w
                } else {
                    y[i]
                }
            } else {
                // Weighted linear regression (degree 1 or 2, quadratic falls back to linear)
                let det = sum_w * sum_wxx - sum_wx * sum_wx;
                if det.abs() < 1e-10 {
                    if sum_w > 1e-10 {
                        sum_wy / sum_w
                    } else {
                        y[i]
                    }
                } else {
                    let b0 = (sum_wxx * sum_wy - sum_wx * sum_wxy) / det;
                    let b1 = (sum_w * sum_wxy - sum_wx * sum_wy) / det;
                    b0 + b1 * xi
                }
            };

            fitted.push(yi);
        }

        Ok(fitted)
    }
}

impl Default for RTCalibrator {
    fn default() -> Self {
        Self::new()
    }
}

/// Find the contiguous range [lo, hi) of k nearest neighbors of x[center]
/// in a sorted array. Returns (lo, hi) where hi - lo == k.
///
/// Since x is sorted, the k nearest neighbors are always a contiguous subrange.
/// Uses two-pointer expansion from center: O(k) instead of O(n log n) sort.
fn find_k_nearest_sorted(x: &[f64], center: usize, k: usize) -> (usize, usize) {
    let n = x.len();
    debug_assert!(k <= n);
    debug_assert!(center < n);

    let mut lo = center;
    let mut hi = center + 1;

    while hi - lo < k {
        let can_go_left = lo > 0;
        let can_go_right = hi < n;

        if can_go_left && can_go_right {
            // Expand toward the closer neighbor (ties go left for stability)
            if (x[center] - x[lo - 1]) <= (x[hi] - x[center]) {
                lo -= 1;
            } else {
                hi += 1;
            }
        } else if can_go_left {
            lo -= 1;
        } else {
            hi += 1;
        }
    }

    (lo, hi)
}

/// Fitted RT calibration curve
#[derive(Debug, Clone)]
pub struct RTCalibration {
    /// Library RTs used for fitting (sorted)
    library_rts: Vec<f64>,
    /// Corresponding measured RTs
    measured_rts: Vec<f64>,
    /// Fitted values at training points
    fitted_values: Vec<f64>,
    /// Absolute residuals at each calibration point (for local tolerance)
    abs_residuals: Vec<f64>,
    /// Bandwidth used for fitting (stored for stats/debugging)
    #[allow(dead_code)]
    bandwidth: f64,
    /// Polynomial degree used (stored for stats/debugging)
    #[allow(dead_code)]
    degree: usize,
    /// Residual standard deviation
    residual_std: f64,
}

impl RTCalibration {
    /// Predict measured RT from library RT
    ///
    /// Uses local interpolation for points within the training range
    /// and linear extrapolation for points outside.
    pub fn predict(&self, library_rt: f64) -> f64 {
        let n = self.library_rts.len();
        if n == 0 {
            return library_rt;
        }

        // Find position in sorted array
        let idx = self
            .library_rts
            .binary_search_by(|&x| x.total_cmp(&library_rt))
            .unwrap_or_else(|i| i);

        if idx == 0 {
            // Extrapolate below range — find first pair with non-zero denominator
            self.extrapolate_below(library_rt)
        } else if idx >= n {
            // Extrapolate above range — find last pair with non-zero denominator
            self.extrapolate_above(library_rt)
        } else {
            // Interpolate between idx-1 and idx
            let x0 = self.library_rts[idx - 1];
            let x1 = self.library_rts[idx];

            if (x1 - x0).abs() < 1e-12 {
                // Duplicate library_rts — average the fitted values to avoid NaN
                let y0 = self.fitted_values[idx - 1];
                let y1 = self.fitted_values[idx];
                (y0 + y1) / 2.0
            } else {
                let y0 = self.fitted_values[idx - 1];
                let y1 = self.fitted_values[idx];
                let t = (library_rt - x0) / (x1 - x0);
                y0 + t * (y1 - y0)
            }
        }
    }

    /// Inverse predict: convert measured RT back to library RT space.
    ///
    /// Searches the `fitted_values` array and interpolates back to `library_rts`.
    /// This is the exact inverse of `predict()` when the LOESS fit is monotonic
    /// (which it virtually always is for LC-MS RT calibration).
    ///
    /// For non-monotonic fits (rare), finds the nearest fitted value and returns
    /// the corresponding library RT.
    pub fn inverse_predict(&self, measured_rt: f64) -> f64 {
        let n = self.fitted_values.len();
        if n == 0 {
            return measured_rt;
        }

        // Check monotonicity of fitted_values (almost always true for RT calibration)
        let is_monotonic = self.fitted_values.windows(2).all(|w| w[1] >= w[0] - 1e-12);

        if is_monotonic {
            // Binary search on fitted_values (they're sorted if monotonic)
            let idx = self
                .fitted_values
                .binary_search_by(|&y| y.total_cmp(&measured_rt))
                .unwrap_or_else(|i| i);

            if idx == 0 {
                // Extrapolate below range
                self.inverse_extrapolate_below(measured_rt)
            } else if idx >= n {
                // Extrapolate above range
                self.inverse_extrapolate_above(measured_rt)
            } else {
                // Interpolate between idx-1 and idx
                let y0 = self.fitted_values[idx - 1];
                let y1 = self.fitted_values[idx];

                if (y1 - y0).abs() < 1e-12 {
                    // Duplicate fitted_values — average the library_rts
                    (self.library_rts[idx - 1] + self.library_rts[idx]) / 2.0
                } else {
                    let x0 = self.library_rts[idx - 1];
                    let x1 = self.library_rts[idx];
                    let t = (measured_rt - y0) / (y1 - y0);
                    x0 + t * (x1 - x0)
                }
            }
        } else {
            // Non-monotonic fallback: find nearest fitted value
            let nearest_idx = self
                .fitted_values
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    (measured_rt - *a)
                        .abs()
                        .total_cmp(&(measured_rt - *b).abs())
                })
                .map(|(i, _)| i)
                .unwrap_or(0);
            self.library_rts[nearest_idx]
        }
    }

    /// Extrapolate below fitted range for inverse prediction
    fn inverse_extrapolate_below(&self, measured_rt: f64) -> f64 {
        let n = self.fitted_values.len();
        if n < 2 {
            return self.library_rts[0];
        }
        for i in 0..n - 1 {
            let dy = self.fitted_values[i + 1] - self.fitted_values[i];
            if dy.abs() > 1e-12 {
                let slope = (self.library_rts[i + 1] - self.library_rts[i]) / dy;
                return self.library_rts[i] + slope * (measured_rt - self.fitted_values[i]);
            }
        }
        self.library_rts[0]
    }

    /// Extrapolate above fitted range for inverse prediction
    fn inverse_extrapolate_above(&self, measured_rt: f64) -> f64 {
        let n = self.fitted_values.len();
        if n < 2 {
            return self.library_rts[n - 1];
        }
        for i in (0..n - 1).rev() {
            let dy = self.fitted_values[i + 1] - self.fitted_values[i];
            if dy.abs() > 1e-12 {
                let slope = (self.library_rts[i + 1] - self.library_rts[i]) / dy;
                return self.library_rts[i + 1] + slope * (measured_rt - self.fitted_values[i + 1]);
            }
        }
        self.library_rts[n - 1]
    }

    /// Extrapolate below the calibration range, skipping duplicate points
    fn extrapolate_below(&self, library_rt: f64) -> f64 {
        let n = self.library_rts.len();
        if n < 2 {
            return self.fitted_values[0];
        }
        // Find first pair with distinct library_rts for a valid slope
        for i in 0..n - 1 {
            let dx = self.library_rts[i + 1] - self.library_rts[i];
            if dx.abs() > 1e-12 {
                let slope = (self.fitted_values[i + 1] - self.fitted_values[i]) / dx;
                return self.fitted_values[i] + slope * (library_rt - self.library_rts[i]);
            }
        }
        // All points have the same library_rt — return average fitted value
        self.fitted_values[0]
    }

    /// Extrapolate above the calibration range, skipping duplicate points
    fn extrapolate_above(&self, library_rt: f64) -> f64 {
        let n = self.library_rts.len();
        if n < 2 {
            return self.fitted_values[n - 1];
        }
        // Find last pair with distinct library_rts for a valid slope
        for i in (0..n - 1).rev() {
            let dx = self.library_rts[i + 1] - self.library_rts[i];
            if dx.abs() > 1e-12 {
                let slope = (self.fitted_values[i + 1] - self.fitted_values[i]) / dx;
                return self.fitted_values[i + 1] + slope * (library_rt - self.library_rts[i + 1]);
            }
        }
        // All points have the same library_rt — return average fitted value
        self.fitted_values[n - 1]
    }

    /// Get the residual standard deviation
    ///
    /// This can be used to set adaptive RT tolerance (e.g., 3× residual_std)
    pub fn residual_std(&self) -> f64 {
        self.residual_std
    }

    /// Get local RT tolerance at a given library RT
    ///
    /// Interpolates absolute residuals from calibration points and applies
    /// smoothing to reduce noise from individual outliers.
    ///
    /// # Arguments
    /// * `library_rt` - The library retention time to query
    /// * `factor` - Multiplier for the local residual (typically 3.0)
    /// * `min_tolerance` - Minimum tolerance floor in minutes
    ///
    /// # Returns
    /// The local RT tolerance in minutes
    pub fn local_tolerance(&self, library_rt: f64, factor: f64, min_tolerance: f64) -> f64 {
        let local_residual = self.interpolate_abs_residual(library_rt);
        (local_residual * factor).max(min_tolerance)
    }

    /// Interpolate absolute residual at a query library RT
    ///
    /// Uses the same interpolation logic as predict(), but operates on
    /// smoothed absolute residuals instead of fitted values.
    fn interpolate_abs_residual(&self, library_rt: f64) -> f64 {
        let n = self.library_rts.len();
        if n == 0 {
            return self.residual_std;
        }

        // Find position in sorted array
        let idx = self
            .library_rts
            .binary_search_by(|&x| x.total_cmp(&library_rt))
            .unwrap_or_else(|i| i);

        if idx == 0 {
            // Extrapolate below range - use smoothed residual at first point
            self.smoothed_abs_residual(0)
        } else if idx >= n {
            // Extrapolate above range - use smoothed residual at last point
            self.smoothed_abs_residual(n - 1)
        } else {
            // Interpolate between idx-1 and idx using smoothed residuals
            let x0 = self.library_rts[idx - 1];
            let x1 = self.library_rts[idx];
            let r0 = self.smoothed_abs_residual(idx - 1);
            let r1 = self.smoothed_abs_residual(idx);

            if (x1 - x0).abs() < 1e-12 {
                // Duplicate library_rts — average the residuals
                (r0 + r1) / 2.0
            } else {
                let t = (library_rt - x0) / (x1 - x0);
                r0 + t * (r1 - r0)
            }
        }
    }

    /// Get smoothed absolute residual at a given index
    ///
    /// Averages residuals over ±2 neighbors to reduce noise from outliers.
    fn smoothed_abs_residual(&self, idx: usize) -> f64 {
        let start = idx.saturating_sub(2);
        let end = (idx + 3).min(self.abs_residuals.len());
        let sum: f64 = self.abs_residuals[start..end].iter().sum();
        sum / (end - start) as f64
    }

    /// Get the range of library RTs used for calibration
    pub fn library_rt_range(&self) -> (f64, f64) {
        if self.library_rts.is_empty() {
            return (0.0, 0.0);
        }
        (
            self.library_rts[0],
            self.library_rts[self.library_rts.len() - 1],
        )
    }

    /// Get the range of measured RTs
    pub fn measured_rt_range(&self) -> (f64, f64) {
        if self.measured_rts.is_empty() {
            return (0.0, 0.0);
        }
        let min = self
            .measured_rts
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max = self
            .measured_rts
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        (min, max)
    }

    /// Get the absolute residuals from the calibration fit.
    pub fn abs_residuals(&self) -> &[f64] {
        &self.abs_residuals
    }

    /// Get calibration statistics
    pub fn stats(&self) -> RTCalibrationStats {
        let n = self.library_rts.len();
        let residuals: Vec<f64> = self
            .measured_rts
            .iter()
            .zip(self.fitted_values.iter())
            .map(|(obs, pred)| obs - pred)
            .collect();

        let mean_residual = if n > 0 {
            residuals.iter().sum::<f64>() / n as f64
        } else {
            0.0
        };

        let max_residual = residuals.iter().fold(0.0f64, |m, &r| m.max(r.abs()));

        // Compute R²
        let y_mean: f64 = self.measured_rts.iter().sum::<f64>() / n as f64;
        let ss_tot: f64 = self
            .measured_rts
            .iter()
            .map(|&y| (y - y_mean).powi(2))
            .sum();
        let ss_res: f64 = residuals.iter().map(|r| r.powi(2)).sum();
        let r_squared = if ss_tot > 1e-10 {
            1.0 - ss_res / ss_tot
        } else {
            0.0
        };

        // Compute percentile stats from absolute residuals
        let p20_abs_residual = percentile_value(&self.abs_residuals, 0.20);
        let p80_abs_residual = percentile_value(&self.abs_residuals, 0.80);

        // MAD = median(|residual_i|) — robust measure of spread
        // Since median of residuals ≈ 0 for LOESS, MAD ≈ P50 of absolute residuals
        let mad = percentile_value(&self.abs_residuals, 0.50);

        RTCalibrationStats {
            n_points: n,
            residual_std: self.residual_std,
            mean_residual,
            max_residual,
            r_squared,
            p20_abs_residual,
            p80_abs_residual,
            mad,
        }
    }

    /// Compute percentile-based RT tolerance (DIA-NN approach).
    ///
    /// Returns `max(2.0 * abs_error_at_percentile, measured_rt_range / 40.0)`.
    /// This is much more robust to outliers than SD-based tolerance.
    pub fn percentile_tolerance(&self, percentile: f64) -> f64 {
        if self.abs_residuals.is_empty() {
            return self.residual_std * 3.0;
        }
        let p_val = percentile_value(&self.abs_residuals, percentile);
        let (rt_min, rt_max) = self.measured_rt_range();
        let rt_range = rt_max - rt_min;
        (2.0 * p_val).max(rt_range / 40.0)
    }

    /// Check if a library RT is within the calibration range
    pub fn is_within_range(&self, library_rt: f64) -> bool {
        if self.library_rts.is_empty() {
            return false;
        }
        library_rt >= self.library_rts[0] && library_rt <= *self.library_rts.last().unwrap()
    }

    /// Export model parameters for serialization
    ///
    /// Returns the library RTs, fitted values, and absolute residuals which
    /// can be used to reconstruct the calibration curve via interpolation.
    pub fn export_model_params(&self) -> RTModelParams {
        RTModelParams {
            library_rts: self.library_rts.clone(),
            fitted_rts: self.fitted_values.clone(),
            abs_residuals: self.abs_residuals.clone(),
        }
    }

    /// Reconstruct an RTCalibration from saved model parameters
    ///
    /// This creates a calibration that can predict using interpolation
    /// between the saved calibration points. The residual_std must be
    /// provided separately (from RTCalibrationParams).
    ///
    /// # Arguments
    /// * `params` - Model parameters containing library_rts and fitted_rts
    /// * `residual_std` - The residual standard deviation from the original fit
    ///
    /// # Returns
    /// A reconstructed RTCalibration that can be used for prediction
    pub fn from_model_params(params: &RTModelParams, residual_std: f64) -> Result<Self> {
        if params.library_rts.len() != params.fitted_rts.len() {
            return Err(OspreyError::ConfigError(format!(
                "Model params have mismatched lengths: {} library_rts vs {} fitted_rts",
                params.library_rts.len(),
                params.fitted_rts.len()
            )));
        }

        if params.library_rts.is_empty() {
            return Err(OspreyError::ConfigError(
                "Model params have no calibration points".to_string(),
            ));
        }

        // Verify library_rts are sorted
        let mut sorted_library = params.library_rts.clone();
        sorted_library.sort_by(|a, b| a.total_cmp(b));

        if sorted_library != params.library_rts {
            return Err(OspreyError::ConfigError(
                "Model params library_rts are not sorted".to_string(),
            ));
        }

        // Handle backwards compatibility: if abs_residuals not present or wrong length,
        // use uniform residual_std as fallback
        let abs_residuals = if params.abs_residuals.len() == params.library_rts.len() {
            params.abs_residuals.clone()
        } else {
            vec![residual_std; params.library_rts.len()]
        };

        Ok(Self {
            library_rts: params.library_rts.clone(),
            // measured_rts not available from params, use fitted values as placeholder
            measured_rts: params.fitted_rts.clone(),
            fitted_values: params.fitted_rts.clone(),
            abs_residuals,
            bandwidth: 0.3, // Default value, not used for prediction
            degree: 1,      // Default value, not used for prediction
            residual_std,
        })
    }
}

/// Statistics for RT calibration quality
#[derive(Debug, Clone, Copy)]
pub struct RTCalibrationStats {
    /// Number of calibration points
    pub n_points: usize,
    /// Residual standard deviation
    pub residual_std: f64,
    /// Mean residual (should be near 0)
    pub mean_residual: f64,
    /// Maximum absolute residual
    pub max_residual: f64,
    /// R² (coefficient of determination)
    pub r_squared: f64,
    /// 20th percentile of absolute residuals (for DIA-NN-style tolerance)
    pub p20_abs_residual: f64,
    /// 80th percentile of absolute residuals
    pub p80_abs_residual: f64,
    /// Median absolute deviation (MAD) of residuals — robust measure of spread
    /// For normal distribution: SD ≈ MAD × 1.4826
    pub mad: f64,
}

/// Stratified sampler for RT calibration discovery
#[derive(Debug, Clone)]
pub struct RTStratifiedSampler {
    /// Number of RT bins for stratification
    n_bins: usize,
    /// Target peptides per bin
    peptides_per_bin: usize,
    /// Minimum total calibration points
    min_calibration_points: usize,
}

impl RTStratifiedSampler {
    /// Create a new stratified sampler with default settings
    pub fn new() -> Self {
        Self {
            n_bins: 10,
            peptides_per_bin: 100,
            min_calibration_points: 50,
        }
    }

    /// Set the number of RT bins
    pub fn with_n_bins(mut self, n_bins: usize) -> Self {
        self.n_bins = n_bins;
        self
    }

    /// Set the target peptides per bin
    pub fn with_peptides_per_bin(mut self, peptides_per_bin: usize) -> Self {
        self.peptides_per_bin = peptides_per_bin;
        self
    }

    /// Set the minimum calibration points
    pub fn with_min_points(mut self, min_points: usize) -> Self {
        self.min_calibration_points = min_points;
        self
    }

    /// Sample library entries stratified by RT
    ///
    /// Returns indices of library entries to use for calibration discovery
    pub fn sample(&self, retention_times: &[f64]) -> Vec<usize> {
        if retention_times.is_empty() {
            return Vec::new();
        }

        // Find RT range
        let (min_rt, max_rt) = retention_times
            .iter()
            .fold((f64::MAX, f64::MIN), |(min, max), &rt| {
                (min.min(rt), max.max(rt))
            });

        if (max_rt - min_rt) < 1e-10 {
            // All same RT, just return first peptides_per_bin entries
            return (0..retention_times.len().min(self.peptides_per_bin)).collect();
        }

        let bin_width = (max_rt - min_rt) / self.n_bins as f64;

        // Assign entries to bins
        let mut bins: Vec<Vec<usize>> = vec![Vec::new(); self.n_bins];
        for (i, &rt) in retention_times.iter().enumerate() {
            let bin = ((rt - min_rt) / bin_width) as usize;
            let bin = bin.min(self.n_bins - 1); // Handle edge case at max_rt
            bins[bin].push(i);
        }

        // Sample from each bin using deterministic selection
        // (first peptides_per_bin entries per bin for reproducibility)
        let mut sampled = Vec::new();
        for bin in bins.iter() {
            let n_take = bin.len().min(self.peptides_per_bin);
            sampled.extend(bin.iter().take(n_take).cloned());
        }

        sampled
    }

    /// Get the minimum calibration points setting
    pub fn min_calibration_points(&self) -> usize {
        self.min_calibration_points
    }
}

impl Default for RTStratifiedSampler {
    fn default() -> Self {
        Self::new()
    }
}

// Helper functions

/// Compute median of a slice
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

/// Compute the value at a given percentile (0.0 to 1.0) of a slice
pub fn percentile_value(values: &[f64], p: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let idx = (p * (sorted.len() - 1) as f64).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

/// Compute standard deviation of a slice
pub fn std_dev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies LOESS calibration fits a linear RT relationship with high R-squared and accurate prediction/extrapolation.
    #[test]
    fn test_rt_calibration_linear() {
        // Linear relationship: measured_rt = 2 * library_rt + 5
        let library_rts: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let measured_rts: Vec<f64> = library_rts.iter().map(|&x| 2.0 * x + 5.0).collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Test prediction
        let pred = calibration.predict(25.0);
        assert!(
            (pred - 55.0).abs() < 0.5,
            "Prediction at 25 should be ~55, got {}",
            pred
        );

        // Test extrapolation
        let pred_extrap = calibration.predict(60.0);
        assert!(
            (pred_extrap - 125.0).abs() < 5.0,
            "Extrapolation at 60 should be ~125, got {}",
            pred_extrap
        );

        // Test stats
        let stats = calibration.stats();
        assert!(
            stats.r_squared > 0.99,
            "R² should be > 0.99 for linear data"
        );
    }

    /// Verifies LOESS calibration captures a sinusoidal nonlinear RT relationship with R-squared above 0.9.
    #[test]
    fn test_rt_calibration_nonlinear() {
        // Nonlinear relationship with some noise
        let library_rts: Vec<f64> = (0..100).map(|i| i as f64 * 0.5).collect();
        let measured_rts: Vec<f64> = library_rts
            .iter()
            .map(|&x| x + 10.0 * (x / 50.0).sin())
            .collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.2);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // LOESS should capture nonlinearity
        let stats = calibration.stats();
        assert!(
            stats.r_squared > 0.9,
            "R² should be > 0.9 for smooth nonlinear data, got {}",
            stats.r_squared
        );
    }

    /// Verifies the RT stratified sampler distributes samples evenly across all bins.
    #[test]
    fn test_stratified_sampler() {
        // Create retention times spanning 0-100
        let rts: Vec<f64> = (0..1000).map(|i| i as f64 / 10.0).collect();

        let sampler = RTStratifiedSampler::new()
            .with_n_bins(10)
            .with_peptides_per_bin(50);

        let sampled = sampler.sample(&rts);

        // Should sample from all bins
        assert!(
            sampled.len() <= 500,
            "Should sample at most 500 (10 bins × 50 each)"
        );
        assert!(sampled.len() >= 100, "Should sample at least 100 entries");

        // Check distribution across bins
        let mut bin_counts = vec![0; 10];
        for &idx in &sampled {
            let rt = rts[idx];
            let bin = (rt / 10.0) as usize;
            let bin = bin.min(9);
            bin_counts[bin] += 1;
        }

        // Each bin should have roughly equal samples
        for count in bin_counts {
            assert!(count > 0, "Each bin should have samples");
        }
    }

    /// Verifies calibration succeeds with exactly the minimum number of points and fails with fewer.
    #[test]
    fn test_calibration_edge_cases() {
        // Minimum points
        let library_rts: Vec<f64> = (0..20).map(|i| i as f64).collect();
        let measured_rts: Vec<f64> = library_rts.iter().map(|&x| x + 1.0).collect();

        let calibrator = RTCalibrator::new();
        let result = calibrator.fit(&library_rts, &measured_rts);
        assert!(result.is_ok(), "Should succeed with exactly min_points");

        // Too few points
        let library_rts_short: Vec<f64> = vec![1.0, 2.0, 3.0];
        let measured_rts_short: Vec<f64> = vec![2.0, 3.0, 4.0];
        let result = calibrator.fit(&library_rts_short, &measured_rts_short);
        assert!(result.is_err(), "Should fail with too few points");
    }

    /// Verifies median computation for odd-length, even-length, and empty slices.
    #[test]
    fn test_median() {
        assert!((median(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-10);
        assert!((median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
        assert!((median(&[]) - 0.0).abs() < 1e-10);
    }

    /// Verifies sample standard deviation calculation against a known reference value.
    #[test]
    fn test_std_dev() {
        let vals = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let sd = std_dev(&vals);
        assert!(
            (sd - 2.138).abs() < 0.01,
            "std dev should be ~2.138, got {}",
            sd
        );
    }

    /// Verifies that local RT tolerance values are computed and respect the minimum floor across the RT range.
    #[test]
    fn test_local_tolerance_varies_with_rt() {
        // Create calibration with varying residuals: small at early RT, large at late RT
        // Simulate a case where early RTs have better prediction accuracy
        let library_rts: Vec<f64> = (0..50).map(|i| i as f64).collect();
        // Add noise that increases with RT
        let measured_rts: Vec<f64> = library_rts
            .iter()
            .map(|&x| x + 0.1 * (x / 10.0).sin() * (x / 50.0 + 0.1))
            .collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Check that local tolerances are computed
        let tol_early = calibration.local_tolerance(5.0, 3.0, 0.25);
        let tol_late = calibration.local_tolerance(45.0, 3.0, 0.25);

        // Both should be at least the minimum floor
        assert!(tol_early >= 0.25, "Early tolerance should be at least 0.25");
        assert!(tol_late >= 0.25, "Late tolerance should be at least 0.25");
    }

    /// Verifies that local tolerance returns the minimum floor when residuals are near zero.
    #[test]
    fn test_local_tolerance_minimum_floor() {
        // Create calibration with very small residuals
        let library_rts: Vec<f64> = (0..50).map(|i| i as f64).collect();
        // Nearly perfect fit - almost no residuals
        let measured_rts: Vec<f64> = library_rts.iter().map(|&x| x + 0.001).collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Local tolerance should respect minimum floor
        let tol = calibration.local_tolerance(25.0, 3.0, 0.25);
        assert!(
            (tol - 0.25).abs() < 0.01,
            "Tolerance should be at minimum floor 0.25, got {}",
            tol
        );
    }

    /// Verifies that local tolerance returns valid values when queried outside the calibration RT range.
    #[test]
    fn test_local_tolerance_extrapolation() {
        // Test extrapolation outside calibration range
        let library_rts: Vec<f64> = (10..40).map(|i| i as f64).collect();
        let measured_rts: Vec<f64> = library_rts.iter().map(|&x| x + 0.5).collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Query outside range - should use endpoint residual
        let tol_below = calibration.local_tolerance(5.0, 3.0, 0.25);
        let tol_above = calibration.local_tolerance(50.0, 3.0, 0.25);

        // Both should be valid (at least minimum floor)
        assert!(
            tol_below >= 0.25,
            "Below-range tolerance should be at least 0.25"
        );
        assert!(
            tol_above >= 0.25,
            "Above-range tolerance should be at least 0.25"
        );
    }

    /// Verifies that absolute residuals survive an export/import roundtrip and produce matching local tolerances.
    #[test]
    fn test_model_params_roundtrip_with_abs_residuals() {
        // Test that abs_residuals are preserved through export/import
        let library_rts: Vec<f64> = (0..30).map(|i| i as f64).collect();
        let measured_rts: Vec<f64> = library_rts
            .iter()
            .map(|&x| x + 0.3 * (x / 15.0).sin())
            .collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Export model params
        let params = calibration.export_model_params();
        assert_eq!(
            params.abs_residuals.len(),
            params.library_rts.len(),
            "abs_residuals should have same length as library_rts"
        );

        // Reconstruct calibration
        let reconstructed =
            RTCalibration::from_model_params(&params, calibration.residual_std()).unwrap();

        // Check that local_tolerance works on reconstructed calibration
        let orig_tol = calibration.local_tolerance(15.0, 3.0, 0.25);
        let recon_tol = reconstructed.local_tolerance(15.0, 3.0, 0.25);
        assert!(
            (orig_tol - recon_tol).abs() < 0.01,
            "Reconstructed calibration should give same local tolerance"
        );
    }

    /// Verifies that loading old model params without abs_residuals falls back to uniform residual_std for tolerance.
    #[test]
    fn test_backwards_compatibility_no_abs_residuals() {
        // Test loading old model params without abs_residuals
        let params = RTModelParams {
            library_rts: vec![0.0, 10.0, 20.0, 30.0],
            fitted_rts: vec![1.0, 11.0, 21.0, 31.0],
            abs_residuals: vec![], // Empty - simulating old format
        };

        let residual_std = 0.5;
        let calibration = RTCalibration::from_model_params(&params, residual_std).unwrap();

        // Should fall back to uniform residual_std for local tolerance
        let tol = calibration.local_tolerance(15.0, 3.0, 0.25);
        // Expected: residual_std * factor = 0.5 * 3.0 = 1.5
        assert!(
            (tol - 1.5).abs() < 0.01,
            "Should use uniform residual_std when abs_residuals empty, got {}",
            tol
        );
    }

    #[test]
    fn test_predict_with_duplicate_library_rts() {
        // Regression test: duplicate library_rts caused division by zero (NaN)
        // in predict(), which silently bypassed binary search in MzRTIndex
        let library_rts = vec![5.0, 10.0, 10.0, 15.0, 20.0];
        let fitted_values = vec![5.1, 10.2, 10.3, 15.1, 20.0];
        let abs_residuals = vec![0.1, 0.2, 0.3, 0.1, 0.0];

        let cal = RTCalibration {
            library_rts: library_rts.clone(),
            measured_rts: fitted_values.clone(), // For test, measured ≈ fitted
            fitted_values: fitted_values.clone(),
            abs_residuals,
            residual_std: 0.5,
            bandwidth: 0.3,
            degree: 1,
        };

        // Test prediction at duplicate point — should NOT return NaN
        let pred = cal.predict(10.0);
        assert!(
            !pred.is_nan(),
            "predict() must not return NaN for duplicate library_rts"
        );
        // Should average the two fitted values at x=10.0: (10.2 + 10.3) / 2
        assert!(
            (pred - 10.25).abs() < 0.01,
            "Expected ~10.25 for duplicate point, got {}",
            pred
        );

        // Test prediction between duplicate and next point — normal interpolation
        let pred_12 = cal.predict(12.0);
        assert!(
            !pred_12.is_nan(),
            "predict() must not return NaN near duplicate points"
        );

        // Test prediction near duplicate point
        let pred_10_1 = cal.predict(10.1);
        assert!(
            !pred_10_1.is_nan(),
            "predict() must not return NaN near duplicate points"
        );

        // Test local_tolerance at duplicate point — should not return NaN
        let tol = cal.local_tolerance(10.0, 3.0, 0.1);
        assert!(
            !tol.is_nan(),
            "local_tolerance() must not return NaN for duplicate library_rts"
        );
        assert!(
            tol >= 0.1,
            "local_tolerance should be at least min_tolerance"
        );

        // Test extrapolation with duplicates at boundaries
        let boundary_rts = vec![5.0, 5.0, 10.0, 15.0];
        let boundary_fitted = vec![5.1, 5.2, 10.1, 15.0];
        let boundary_residuals = vec![0.1, 0.2, 0.1, 0.0];
        let cal2 = RTCalibration {
            library_rts: boundary_rts,
            measured_rts: boundary_fitted.clone(),
            fitted_values: boundary_fitted,
            abs_residuals: boundary_residuals,
            residual_std: 0.5,
            bandwidth: 0.3,
            degree: 1,
        };

        let pred_below = cal2.predict(3.0);
        assert!(
            !pred_below.is_nan(),
            "Extrapolation below range with duplicate first points must not return NaN"
        );

        // Test extrapolation above with duplicates at end
        let end_rts = vec![5.0, 10.0, 15.0, 15.0];
        let end_fitted = vec![5.1, 10.1, 15.0, 15.1];
        let end_residuals = vec![0.1, 0.1, 0.0, 0.1];
        let cal3 = RTCalibration {
            library_rts: end_rts,
            measured_rts: end_fitted.clone(),
            fitted_values: end_fitted,
            abs_residuals: end_residuals,
            residual_std: 0.5,
            bandwidth: 0.3,
            degree: 1,
        };

        let pred_above = cal3.predict(20.0);
        assert!(
            !pred_above.is_nan(),
            "Extrapolation above range with duplicate last points must not return NaN"
        );
    }

    /// Verifies inverse_predict roundtrips with predict for a linear calibration.
    #[test]
    fn test_inverse_predict_roundtrip_linear() {
        let library_rts: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let measured_rts: Vec<f64> = library_rts.iter().map(|&x| 2.0 * x + 5.0).collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Test roundtrip at several points within range
        for &lib_rt in &[5.0, 15.0, 25.0, 35.0, 45.0] {
            let measured = calibration.predict(lib_rt);
            let recovered = calibration.inverse_predict(measured);
            assert!(
                (recovered - lib_rt).abs() < 0.5,
                "inverse_predict(predict({})) = {}, expected ~{}",
                lib_rt,
                recovered,
                lib_rt
            );
        }

        // Test extrapolation roundtrip
        let measured_extrap = calibration.predict(55.0);
        let recovered_extrap = calibration.inverse_predict(measured_extrap);
        assert!(
            (recovered_extrap - 55.0).abs() < 2.0,
            "Extrapolation roundtrip: inverse_predict(predict(55)) = {}, expected ~55",
            recovered_extrap
        );
    }

    /// Verifies inverse_predict roundtrips for a nonlinear (sinusoidal) calibration.
    #[test]
    fn test_inverse_predict_roundtrip_nonlinear() {
        let library_rts: Vec<f64> = (0..100).map(|i| i as f64 * 0.5).collect();
        let measured_rts: Vec<f64> = library_rts
            .iter()
            .map(|&x| x + 2.0 * (x / 20.0).sin())
            .collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.2);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        for &lib_rt in &[5.0, 15.0, 25.0, 35.0, 45.0] {
            let measured = calibration.predict(lib_rt);
            let recovered = calibration.inverse_predict(measured);
            assert!(
                (recovered - lib_rt).abs() < 1.0,
                "Nonlinear roundtrip: inverse_predict(predict({})) = {}, expected ~{}",
                lib_rt,
                recovered,
                lib_rt
            );
        }
    }

    /// Verifies inverse_predict handles duplicate fitted values without NaN.
    #[test]
    fn test_inverse_predict_duplicate_fitted() {
        let cal = RTCalibration {
            library_rts: vec![5.0, 10.0, 15.0, 20.0],
            measured_rts: vec![5.5, 10.5, 10.5, 20.5],
            fitted_values: vec![5.5, 10.5, 10.5, 20.5],
            abs_residuals: vec![0.1, 0.1, 0.1, 0.1],
            residual_std: 0.5,
            bandwidth: 0.3,
            degree: 1,
        };

        let result = cal.inverse_predict(10.5);
        assert!(
            !result.is_nan(),
            "inverse_predict must not return NaN for duplicate fitted_values"
        );
        // Should average the library_rts: (10 + 15) / 2 = 12.5
        assert!(
            (result - 12.5).abs() < 0.01,
            "Expected ~12.5 for duplicate fitted_values, got {}",
            result
        );
    }
}
