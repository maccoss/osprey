//! m/z calibration for MS1 and MS2 measurements
//!
//! Corrects systematic mass measurement errors using mean PPM offset
//! and provides adaptive tolerance based on measurement variability.

use super::{MzCalibration, MzHistogram};

/// QC data for m/z calibration
#[derive(Debug, Clone, Default)]
pub struct MzQCData {
    /// MS1 m/z errors in PPM
    pub ms1_errors_ppm: Vec<f64>,

    /// MS2 m/z errors in PPM
    pub ms2_errors_ppm: Vec<f64>,
}

impl MzQCData {
    /// Create new empty QC data
    pub fn new() -> Self {
        Self::default()
    }

    /// Add MS1 error measurement
    pub fn add_ms1_error(&mut self, error_ppm: f64) {
        self.ms1_errors_ppm.push(error_ppm);
    }

    /// Add MS2 error measurement
    pub fn add_ms2_error(&mut self, error_ppm: f64) {
        self.ms2_errors_ppm.push(error_ppm);
    }

    /// Number of MS1 observations
    pub fn n_ms1(&self) -> usize {
        self.ms1_errors_ppm.len()
    }

    /// Number of MS2 observations
    pub fn n_ms2(&self) -> usize {
        self.ms2_errors_ppm.len()
    }
}

/// Calculate m/z calibration parameters from QC data
///
/// Uses mean PPM offset to correct systematic mass measurement errors,
/// and calculates adaptive tolerance as mean + 3*SD.
///
/// # Arguments
/// * `qc_data` - QC data from confident calibration matches
///
/// # Returns
/// Calibration parameters for MS1 and MS2
///
/// # Example
/// ```
/// use osprey_chromatography::calibration::mass::{calculate_mz_calibration, MzQCData};
///
/// let qc_data = MzQCData {
///     ms1_errors_ppm: vec![-2.5, -2.3, -2.7, -2.4],
///     ms2_errors_ppm: vec![1.2, 1.1, 1.3, 1.0],
/// };
///
/// let (ms1_cal, ms2_cal) = calculate_mz_calibration(&qc_data);
/// println!("MS1 offset: {:.2} ppm (SD: {:.2})", ms1_cal.mean, ms1_cal.sd);
/// println!("MS2 offset: {:.2} ppm (SD: {:.2})", ms2_cal.mean, ms2_cal.sd);
/// ```
pub fn calculate_mz_calibration(qc_data: &MzQCData) -> (MzCalibration, MzCalibration) {
    let ms1_calibration = calculate_single_calibration(&qc_data.ms1_errors_ppm);
    let ms2_calibration = calculate_single_calibration(&qc_data.ms2_errors_ppm);
    (ms1_calibration, ms2_calibration)
}

/// Calculate calibration parameters for a single error vector
fn calculate_single_calibration(errors: &[f64]) -> MzCalibration {
    if errors.is_empty() {
        return MzCalibration {
            mean: 0.0,
            median: 0.0,
            sd: 0.0,
            count: 0,
            unit: "ppm".to_string(),
            adjusted_tolerance: None,
            window_halfwidth_multiplier: None,
            histogram: None,
            calibrated: false,
        };
    }

    // Calculate mean
    let mean = errors.iter().sum::<f64>() / errors.len() as f64;

    // Calculate median
    let mut sorted = errors.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    // Calculate standard deviation
    let variance = if errors.len() > 1 {
        errors.iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>() / (errors.len() - 1) as f64
    } else {
        0.0
    };
    let sd = variance.sqrt();

    // Generate histogram with 1 ppm bins
    let histogram = generate_histogram(&sorted, 1.0);

    MzCalibration {
        mean,
        median,
        sd,
        count: errors.len(),
        unit: "ppm".to_string(),
        adjusted_tolerance: Some(mean.abs() + 3.0 * sd),
        window_halfwidth_multiplier: Some(3.0),
        histogram: Some(histogram),
        calibrated: true,
    }
}

/// Generate histogram with specified bin width
fn generate_histogram(sorted_values: &[f64], bin_width: f64) -> MzHistogram {
    if sorted_values.is_empty() {
        return MzHistogram {
            bin_edges: vec![],
            counts: vec![],
            bin_width,
        };
    }

    let min_val = sorted_values[0];
    let max_val = sorted_values[sorted_values.len() - 1];

    // Round to nearest bin boundary
    let bin_start = (min_val / bin_width).floor() * bin_width;
    let bin_end = (max_val / bin_width).ceil() * bin_width;

    let num_bins = ((bin_end - bin_start) / bin_width).ceil() as usize;
    let num_bins = num_bins.max(1);

    let mut bin_edges = Vec::with_capacity(num_bins + 1);
    let mut counts = vec![0usize; num_bins];

    for i in 0..=num_bins {
        bin_edges.push(bin_start + i as f64 * bin_width);
    }

    // Count values in each bin
    for &val in sorted_values {
        let bin_idx = ((val - bin_start) / bin_width).floor() as usize;
        let bin_idx = bin_idx.min(num_bins - 1); // Clamp to last bin
        counts[bin_idx] += 1;
    }

    MzHistogram {
        bin_edges,
        counts,
        bin_width,
    }
}

/// Apply m/z calibration to correct an observed m/z value
///
/// # Arguments
/// * `observed_mz` - Measured m/z value
/// * `calibration` - Calibration parameters
///
/// # Returns
/// Corrected m/z value
///
/// # Example
/// ```
/// use osprey_chromatography::calibration::{MzCalibration, mass::apply_mz_calibration};
///
/// let calibration = MzCalibration {
///     mean: -2.5,  // Systematic error of -2.5 ppm
///     median: -2.5,
///     sd: 0.5,
///     count: 100,
///     unit: "ppm".to_string(),
///     adjusted_tolerance: Some(3.0),
///     window_halfwidth_multiplier: Some(3.0),
///     histogram: None,
///     calibrated: true,
/// };
///
/// let observed = 500.0;
/// let corrected = apply_mz_calibration(observed, &calibration);
///
/// // Correction: 500.0 - (500.0 * -2.5 / 1e6) = 500.00125
/// assert!((corrected - 500.00125).abs() < 0.00001);
/// ```
pub fn apply_mz_calibration(observed_mz: f64, calibration: &MzCalibration) -> f64 {
    if !calibration.calibrated {
        return observed_mz;
    }

    // Correct for systematic offset using mean
    let correction_da = observed_mz * calibration.mean / 1e6;
    observed_mz - correction_da
}

/// Calculate m/z error in PPM
///
/// # Arguments
/// * `observed` - Observed m/z
/// * `theoretical` - Theoretical m/z
///
/// # Returns
/// Error in PPM (positive = observed > theoretical)
pub fn calculate_ppm_error(observed: f64, theoretical: f64) -> f64 {
    ((observed - theoretical) / theoretical) * 1e6
}

/// Check if an m/z match is within calibrated tolerance
///
/// # Arguments
/// * `observed_mz` - Observed m/z value
/// * `theoretical_mz` - Theoretical m/z value
/// * `calibration` - Calibration parameters
/// * `base_tolerance_ppm` - Base tolerance to use if not calibrated
///
/// # Returns
/// True if the observed m/z is within tolerance of theoretical
pub fn is_within_calibrated_tolerance(
    observed_mz: f64,
    theoretical_mz: f64,
    calibration: &MzCalibration,
    base_tolerance_ppm: f64,
) -> bool {
    let error_ppm = calculate_ppm_error(observed_mz, theoretical_mz);

    if calibration.calibrated {
        // Use calibrated tolerance centered on mean
        let tolerance = calibration.adjusted_tolerance.unwrap_or(base_tolerance_ppm);
        (error_ppm - calibration.mean).abs() <= tolerance
    } else {
        // Use base tolerance centered on 0
        error_ppm.abs() <= base_tolerance_ppm
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ppm_error_calculation() {
        // Observed = 500.001, Theoretical = 500.0
        // Error = (500.001 - 500.0) / 500.0 * 1e6 = 2.0 ppm
        let error = calculate_ppm_error(500.001, 500.0);
        assert!((error - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_apply_calibration_negative_offset() {
        let calibration = MzCalibration {
            mean: -2.5,
            median: -2.5,
            sd: 0.1,
            count: 100,
            unit: "ppm".to_string(),
            adjusted_tolerance: Some(2.5 + 3.0 * 0.1),
            window_halfwidth_multiplier: Some(3.0),
            histogram: None,
            calibrated: true,
        };

        let observed = 500.0;
        let corrected = apply_mz_calibration(observed, &calibration);

        // Correction: 500.0 - (500.0 * -2.5 / 1e6) = 500.0 + 0.00125 = 500.00125
        assert!((corrected - 500.00125).abs() < 0.00001);
    }

    #[test]
    fn test_apply_calibration_positive_offset() {
        let calibration = MzCalibration {
            mean: 2.5,
            median: 2.5,
            sd: 0.1,
            count: 100,
            unit: "ppm".to_string(),
            adjusted_tolerance: Some(2.5 + 3.0 * 0.1),
            window_halfwidth_multiplier: Some(3.0),
            histogram: None,
            calibrated: true,
        };

        let observed = 500.0;
        let corrected = apply_mz_calibration(observed, &calibration);

        // Correction: 500.0 - (500.0 * 2.5 / 1e6) = 500.0 - 0.00125 = 499.99875
        assert!((corrected - 499.99875).abs() < 0.00001);
    }

    #[test]
    fn test_calculate_mz_calibration() {
        let qc_data = MzQCData {
            ms1_errors_ppm: vec![-2.5, -2.3, -2.7, -2.4, -2.6],
            ms2_errors_ppm: vec![1.2, 1.1, 1.3, 1.0, 1.4],
        };

        let (ms1_cal, ms2_cal) = calculate_mz_calibration(&qc_data);

        // MS1 mean should be -2.5, median should be -2.5
        assert!((ms1_cal.mean + 2.5).abs() < 0.01);
        assert!((ms1_cal.median + 2.5).abs() < 0.01);
        assert!(ms1_cal.sd > 0.0);
        assert_eq!(ms1_cal.count, 5);
        assert!(ms1_cal.calibrated);
        assert!(ms1_cal.histogram.is_some());

        // MS2 mean should be 1.2, median should be 1.2
        assert!((ms2_cal.mean - 1.2).abs() < 0.01);
        assert!((ms2_cal.median - 1.2).abs() < 0.01);
        assert!(ms2_cal.sd > 0.0);
        assert_eq!(ms2_cal.count, 5);
        assert!(ms2_cal.calibrated);
        assert!(ms2_cal.histogram.is_some());
    }

    #[test]
    fn test_empty_qc_data() {
        let qc_data = MzQCData {
            ms1_errors_ppm: vec![],
            ms2_errors_ppm: vec![],
        };

        let (ms1_cal, ms2_cal) = calculate_mz_calibration(&qc_data);

        assert_eq!(ms1_cal.mean, 0.0);
        assert_eq!(ms1_cal.sd, 0.0);
        assert!(!ms1_cal.calibrated);

        assert_eq!(ms2_cal.mean, 0.0);
        assert_eq!(ms2_cal.sd, 0.0);
        assert!(!ms2_cal.calibrated);
    }

    #[test]
    fn test_uncalibrated_passthrough() {
        let calibration = MzCalibration {
            mean: 0.0,
            median: 0.0,
            sd: 0.0,
            count: 0,
            unit: "ppm".to_string(),
            adjusted_tolerance: None,
            window_halfwidth_multiplier: None,
            histogram: None,
            calibrated: false,
        };

        let observed = 500.12345;
        let corrected = apply_mz_calibration(observed, &calibration);
        assert_eq!(corrected, observed);
    }

    #[test]
    fn test_within_calibrated_tolerance() {
        let calibration = MzCalibration {
            mean: -2.5,
            median: -2.5,
            sd: 1.0,
            count: 100,
            unit: "ppm".to_string(),
            adjusted_tolerance: Some(5.5), // -2.5 + 3*1.0
            window_halfwidth_multiplier: Some(3.0),
            histogram: None,
            calibrated: true,
        };

        // Observed with -2.5 ppm error (exactly at mean) should be within tolerance
        let theoretical = 500.0;
        let observed_at_mean = theoretical * (1.0 - 2.5 / 1e6);
        assert!(is_within_calibrated_tolerance(observed_at_mean, theoretical, &calibration, 10.0));

        // Observed with -6.0 ppm error (3.5 ppm from mean) should be within tolerance
        let observed_within = theoretical * (1.0 - 6.0 / 1e6);
        assert!(is_within_calibrated_tolerance(observed_within, theoretical, &calibration, 10.0));

        // Observed with +5.0 ppm error (7.5 ppm from mean) should be outside tolerance
        let observed_outside = theoretical * (1.0 + 5.0 / 1e6);
        assert!(!is_within_calibrated_tolerance(observed_outside, theoretical, &calibration, 10.0));
    }

    #[test]
    fn test_histogram_generation() {
        let qc_data = MzQCData {
            ms1_errors_ppm: vec![-2.5, -2.3, -2.7, -1.8, -3.2, 0.5, -2.1],
            ms2_errors_ppm: vec![],
        };

        let (ms1_cal, _) = calculate_mz_calibration(&qc_data);
        let hist = ms1_cal.histogram.unwrap();

        // Check histogram properties
        assert!(!hist.bin_edges.is_empty());
        assert!(!hist.counts.is_empty());
        assert_eq!(hist.bin_edges.len(), hist.counts.len() + 1);

        // Total counts should equal number of observations
        let total_count: usize = hist.counts.iter().sum();
        assert_eq!(total_count, 7);
    }

    #[test]
    fn test_median_calculation() {
        // Odd number of elements
        let qc_data = MzQCData {
            ms1_errors_ppm: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            ms2_errors_ppm: vec![],
        };
        let (ms1_cal, _) = calculate_mz_calibration(&qc_data);
        assert!((ms1_cal.median - 3.0).abs() < 0.001);

        // Even number of elements
        let qc_data = MzQCData {
            ms1_errors_ppm: vec![1.0, 2.0, 3.0, 4.0],
            ms2_errors_ppm: vec![],
        };
        let (ms1_cal, _) = calculate_mz_calibration(&qc_data);
        assert!((ms1_cal.median - 2.5).abs() < 0.001);
    }
}
