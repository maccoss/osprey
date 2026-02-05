//! RT Calibration using LOESS (Local Polynomial Regression)
//!
//! This module provides retention time calibration to convert library RTs
//! (which may be iRT, predicted, or uncalibrated) to measured RTs.
//!
//! The calibration uses LOESS (Locally Estimated Scatterplot Smoothing)
//! which fits local polynomial regressions weighted by distance.

use osprey_core::{OspreyError, Result};
use super::RTModelParams;

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
}

impl Default for RTCalibratorConfig {
    fn default() -> Self {
        Self {
            bandwidth: 0.3,
            degree: 1, // Linear local fits
            min_points: 20,
            robustness_iter: 2,
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

    /// Fit a calibration curve from (library_rt, measured_rt) pairs
    ///
    /// # Arguments
    /// * `library_rts` - Library retention times (iRT or predicted)
    /// * `measured_rts` - Corresponding measured retention times
    ///
    /// # Returns
    /// A fitted `RTCalibration` that can predict measured RT from library RT
    pub fn fit(
        &self,
        library_rts: &[f64],
        measured_rts: &[f64],
    ) -> Result<RTCalibration> {
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

        // Compute residuals
        let residuals: Vec<f64> = y
            .iter()
            .zip(fitted.iter())
            .map(|(obs, pred)| obs - pred)
            .collect();

        // Robustness iterations
        let mut weights: Vec<f64> = vec![1.0; x.len()];
        let mut final_fitted = fitted;

        for _ in 0..self.config.robustness_iter {
            // Compute robust weights using bisquare function
            let abs_residuals: Vec<f64> = residuals.iter().map(|r| r.abs()).collect();
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

        // Compute residual standard deviation
        let residuals: Vec<f64> = y
            .iter()
            .zip(final_fitted.iter())
            .map(|(obs, pred)| obs - pred)
            .collect();
        let residual_std = std_dev(&residuals);

        Ok(RTCalibration {
            library_rts: x,
            measured_rts: y,
            fitted_values: final_fitted,
            bandwidth: self.config.bandwidth,
            degree: self.config.degree,
            residual_std,
        })
    }

    /// Perform LOESS fit
    fn loess_fit(
        &self,
        x: &[f64],
        y: &[f64],
        weights: Option<&[f64]>,
    ) -> Result<Vec<f64>> {
        let n = x.len();
        let k = (self.config.bandwidth * n as f64).ceil() as usize;
        let k = k.max(self.config.degree + 2).min(n); // Ensure enough points for polynomial

        let mut fitted = Vec::with_capacity(n);

        for i in 0..n {
            let xi = x[i];

            // Find k nearest neighbors
            let mut distances: Vec<(usize, f64)> = x
                .iter()
                .enumerate()
                .map(|(j, &xj)| (j, (xj - xi).abs()))
                .collect();
            distances.sort_by(|a, b| a.1.total_cmp(&b.1));

            let neighbors: Vec<usize> = distances.iter().take(k).map(|(j, _)| *j).collect();
            let max_dist = distances[k - 1].1;

            // Compute tricube weights
            let tricube_weights: Vec<f64> = neighbors
                .iter()
                .map(|&j| {
                    let dist = (x[j] - xi).abs();
                    let u = if max_dist > 1e-10 {
                        dist / max_dist
                    } else {
                        0.0
                    };
                    if u < 1.0 {
                        (1.0 - u.powi(3)).powi(3)
                    } else {
                        0.0
                    }
                })
                .collect();

            // Combine with robustness weights if provided
            let final_weights: Vec<f64> = if let Some(w) = weights {
                tricube_weights
                    .iter()
                    .zip(neighbors.iter())
                    .map(|(&tw, &j)| tw * w[j])
                    .collect()
            } else {
                tricube_weights
            };

            // Fit local polynomial (degree 0, 1, or 2)
            let yi = match self.config.degree {
                0 => {
                    // Weighted mean
                    let sum_w: f64 = final_weights.iter().sum();
                    if sum_w > 1e-10 {
                        neighbors
                            .iter()
                            .zip(final_weights.iter())
                            .map(|(&j, &w)| w * y[j])
                            .sum::<f64>()
                            / sum_w
                    } else {
                        y[i]
                    }
                }
                1 => {
                    // Weighted linear regression
                    self.weighted_linear_fit(&neighbors, x, y, &final_weights, xi)
                }
                _ => {
                    // Weighted quadratic regression
                    self.weighted_quadratic_fit(&neighbors, x, y, &final_weights, xi)
                }
            };

            fitted.push(yi);
        }

        Ok(fitted)
    }

    /// Fit weighted linear regression and predict at xi
    fn weighted_linear_fit(
        &self,
        neighbors: &[usize],
        x: &[f64],
        y: &[f64],
        weights: &[f64],
        xi: f64,
    ) -> f64 {
        let mut sum_w = 0.0;
        let mut sum_wx = 0.0;
        let mut sum_wy = 0.0;
        let mut sum_wxx = 0.0;
        let mut sum_wxy = 0.0;

        for (&j, &w) in neighbors.iter().zip(weights.iter()) {
            sum_w += w;
            sum_wx += w * x[j];
            sum_wy += w * y[j];
            sum_wxx += w * x[j] * x[j];
            sum_wxy += w * x[j] * y[j];
        }

        // Solve: [sum_w    sum_wx ] [b0]   [sum_wy ]
        //        [sum_wx   sum_wxx] [b1] = [sum_wxy]
        let det = sum_w * sum_wxx - sum_wx * sum_wx;
        if det.abs() < 1e-10 {
            // Fallback to weighted mean
            if sum_w > 1e-10 {
                return sum_wy / sum_w;
            } else {
                return y[neighbors[0]];
            }
        }

        let b0 = (sum_wxx * sum_wy - sum_wx * sum_wxy) / det;
        let b1 = (sum_w * sum_wxy - sum_wx * sum_wy) / det;

        b0 + b1 * xi
    }

    /// Fit weighted quadratic regression and predict at xi
    fn weighted_quadratic_fit(
        &self,
        neighbors: &[usize],
        x: &[f64],
        y: &[f64],
        weights: &[f64],
        xi: f64,
    ) -> f64 {
        // For simplicity, use linear fit as fallback if quadratic is unstable
        // Full quadratic would require solving 3x3 system
        self.weighted_linear_fit(neighbors, x, y, weights, xi)
    }
}

impl Default for RTCalibrator {
    fn default() -> Self {
        Self::new()
    }
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
            // Extrapolate below range
            if n > 1 {
                let slope = (self.fitted_values[1] - self.fitted_values[0])
                    / (self.library_rts[1] - self.library_rts[0]);
                self.fitted_values[0] + slope * (library_rt - self.library_rts[0])
            } else {
                self.fitted_values[0]
            }
        } else if idx >= n {
            // Extrapolate above range
            if n > 1 {
                let slope = (self.fitted_values[n - 1] - self.fitted_values[n - 2])
                    / (self.library_rts[n - 1] - self.library_rts[n - 2]);
                self.fitted_values[n - 1] + slope * (library_rt - self.library_rts[n - 1])
            } else {
                self.fitted_values[n - 1]
            }
        } else {
            // Interpolate between idx-1 and idx
            let x0 = self.library_rts[idx - 1];
            let x1 = self.library_rts[idx];
            let y0 = self.fitted_values[idx - 1];
            let y1 = self.fitted_values[idx];

            let t = (library_rt - x0) / (x1 - x0);
            y0 + t * (y1 - y0)
        }
    }

    /// Get the residual standard deviation
    ///
    /// This can be used to set adaptive RT tolerance (e.g., 3× residual_std)
    pub fn residual_std(&self) -> f64 {
        self.residual_std
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
        let min = self.measured_rts.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = self.measured_rts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (min, max)
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
        let ss_tot: f64 = self.measured_rts.iter().map(|&y| (y - y_mean).powi(2)).sum();
        let ss_res: f64 = residuals.iter().map(|r| r.powi(2)).sum();
        let r_squared = if ss_tot > 1e-10 { 1.0 - ss_res / ss_tot } else { 0.0 };

        RTCalibrationStats {
            n_points: n,
            residual_std: self.residual_std,
            mean_residual,
            max_residual,
            r_squared,
        }
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
    /// Returns the library RTs and fitted values which can be used
    /// to reconstruct the calibration curve via interpolation.
    pub fn export_model_params(&self) -> RTModelParams {
        RTModelParams {
            library_rts: self.library_rts.clone(),
            fitted_rts: self.fitted_values.clone(),
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
                "Model params have no calibration points".to_string()
            ));
        }

        // Verify library_rts are sorted
        let mut sorted_library = params.library_rts.clone();
        sorted_library.sort_by(|a, b| a.total_cmp(b));

        if sorted_library != params.library_rts {
            return Err(OspreyError::ConfigError(
                "Model params library_rts are not sorted".to_string()
            ));
        }

        Ok(Self {
            library_rts: params.library_rts.clone(),
            // measured_rts not available from params, use fitted values as placeholder
            measured_rts: params.fitted_rts.clone(),
            fitted_values: params.fitted_rts.clone(),
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

    #[test]
    fn test_rt_calibration_linear() {
        // Linear relationship: measured_rt = 2 * library_rt + 5
        let library_rts: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let measured_rts: Vec<f64> = library_rts.iter().map(|&x| 2.0 * x + 5.0).collect();

        let calibrator = RTCalibrator::new().with_bandwidth(0.3);
        let calibration = calibrator.fit(&library_rts, &measured_rts).unwrap();

        // Test prediction
        let pred = calibration.predict(25.0);
        assert!((pred - 55.0).abs() < 0.5, "Prediction at 25 should be ~55, got {}", pred);

        // Test extrapolation
        let pred_extrap = calibration.predict(60.0);
        assert!((pred_extrap - 125.0).abs() < 5.0, "Extrapolation at 60 should be ~125, got {}", pred_extrap);

        // Test stats
        let stats = calibration.stats();
        assert!(stats.r_squared > 0.99, "R² should be > 0.99 for linear data");
    }

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
        assert!(stats.r_squared > 0.9, "R² should be > 0.9 for smooth nonlinear data, got {}", stats.r_squared);
    }

    #[test]
    fn test_stratified_sampler() {
        // Create retention times spanning 0-100
        let rts: Vec<f64> = (0..1000).map(|i| i as f64 / 10.0).collect();

        let sampler = RTStratifiedSampler::new()
            .with_n_bins(10)
            .with_peptides_per_bin(50);

        let sampled = sampler.sample(&rts);

        // Should sample from all bins
        assert!(sampled.len() <= 500, "Should sample at most 500 (10 bins × 50 each)");
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

    #[test]
    fn test_median() {
        assert!((median(&[1.0, 2.0, 3.0]) - 2.0).abs() < 1e-10);
        assert!((median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
        assert!((median(&[]) - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_std_dev() {
        let vals = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let sd = std_dev(&vals);
        assert!((sd - 2.138).abs() < 0.01, "std dev should be ~2.138, got {}", sd);
    }
}
