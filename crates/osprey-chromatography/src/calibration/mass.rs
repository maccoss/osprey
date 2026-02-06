//! m/z calibration for MS1 and MS2 measurements
//!
//! Corrects systematic mass measurement errors using mean offset
//! and provides adaptive tolerance based on measurement variability.
//!
//! Supports both ppm-based calibration (HRAM) and absolute m/z (Th) calibration
//! (unit resolution).

use super::{MzCalibration, MzHistogram};
use osprey_core::{Spectrum, ToleranceUnit};

/// QC data for m/z calibration
///
/// For unit resolution instruments (e.g., Stellar), both MS1 and MS2
/// errors are measured in Th (Thomson). For HRAM instruments, both
/// are measured in ppm.
#[derive(Debug, Clone)]
pub struct MzQCData {
    /// MS1 m/z errors (in configured unit - ppm or Th)
    pub ms1_errors: Vec<f64>,

    /// MS2 m/z errors (in configured unit - ppm or Th)
    pub ms2_errors: Vec<f64>,

    /// Unit for mass errors (same for both MS1 and MS2)
    pub unit: ToleranceUnit,
}

impl Default for MzQCData {
    fn default() -> Self {
        Self {
            ms1_errors: Vec::new(),
            ms2_errors: Vec::new(),
            unit: ToleranceUnit::Ppm,
        }
    }
}

impl MzQCData {
    /// Create new empty QC data with specified unit
    ///
    /// The unit applies to both MS1 and MS2 errors:
    /// - `ToleranceUnit::Ppm` for HRAM instruments
    /// - `ToleranceUnit::Mz` (Th) for unit resolution instruments
    pub fn new(unit: ToleranceUnit) -> Self {
        Self {
            ms1_errors: Vec::new(),
            ms2_errors: Vec::new(),
            unit,
        }
    }

    /// Add MS1 error measurement (in configured unit)
    pub fn add_ms1_error(&mut self, error: f64) {
        self.ms1_errors.push(error);
    }

    /// Add MS2 error measurement (in configured unit)
    pub fn add_ms2_error(&mut self, error: f64) {
        self.ms2_errors.push(error);
    }

    /// Number of MS1 observations
    pub fn n_ms1(&self) -> usize {
        self.ms1_errors.len()
    }

    /// Number of MS2 observations
    pub fn n_ms2(&self) -> usize {
        self.ms2_errors.len()
    }
}

/// Calculate m/z calibration parameters from QC data
///
/// Uses mean offset to correct systematic mass measurement errors,
/// and calculates adaptive tolerance as |mean| + 3*SD.
///
/// # Arguments
/// * `qc_data` - QC data from confident calibration matches
///
/// # Returns
/// Calibration parameters for MS1 and MS2 (both use same unit from qc_data)
pub fn calculate_mz_calibration(qc_data: &MzQCData) -> (MzCalibration, MzCalibration) {
    // Both MS1 and MS2 use the same unit (ppm for HRAM, Th for unit resolution)
    let ms1_calibration = calculate_single_calibration(&qc_data.ms1_errors, qc_data.unit);
    let ms2_calibration = calculate_single_calibration(&qc_data.ms2_errors, qc_data.unit);
    (ms1_calibration, ms2_calibration)
}

/// Calculate calibration parameters for a single error vector
fn calculate_single_calibration(errors: &[f64], unit: ToleranceUnit) -> MzCalibration {
    let unit_str = match unit {
        ToleranceUnit::Ppm => "ppm",
        ToleranceUnit::Mz => "Th",
    };

    if errors.is_empty() {
        return MzCalibration {
            mean: 0.0,
            median: 0.0,
            sd: 0.0,
            count: 0,
            unit: unit_str.to_string(),
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

    // Generate histogram with appropriate bin width
    // ppm: 1.0 bins, Th: 0.01 bins
    let bin_width = match unit {
        ToleranceUnit::Ppm => 1.0,
        ToleranceUnit::Mz => 0.01,
    };
    let histogram = generate_histogram(&sorted, bin_width);

    MzCalibration {
        mean,
        median,
        sd,
        count: errors.len(),
        unit: unit_str.to_string(),
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
/// # Unit-aware correction
/// - For ppm: `corrected = observed - observed * mean / 1e6`
/// - For Th: `corrected = observed - mean`
pub fn apply_mz_calibration(observed_mz: f64, calibration: &MzCalibration) -> f64 {
    if !calibration.calibrated {
        return observed_mz;
    }

    // Correct for systematic offset using mean
    // Unit is determined by calibration.unit field
    if calibration.unit == "Th" {
        // Absolute m/z correction (unit resolution)
        observed_mz - calibration.mean
    } else {
        // PPM correction (HRAM)
        let correction_da = observed_mz * calibration.mean / 1e6;
        observed_mz - correction_da
    }
}

/// Apply m/z calibration to all peaks in a spectrum
///
/// Creates a new spectrum with corrected m/z values. This shifts all
/// measured m/z values by the systematic offset determined during calibration.
///
/// # Arguments
/// * `spectrum` - Observed spectrum with uncorrected m/z values
/// * `calibration` - MS2 calibration parameters
///
/// # Returns
/// New spectrum with calibration-corrected m/z values
pub fn apply_spectrum_calibration(spectrum: &Spectrum, calibration: &MzCalibration) -> Spectrum {
    if !calibration.calibrated {
        return spectrum.clone();
    }

    // Correct each m/z value
    let corrected_mzs: Vec<f64> = spectrum.mzs
        .iter()
        .map(|&mz| apply_mz_calibration(mz, calibration))
        .collect();

    Spectrum {
        scan_number: spectrum.scan_number,
        retention_time: spectrum.retention_time,
        precursor_mz: spectrum.precursor_mz,
        isolation_window: spectrum.isolation_window.clone(),
        mzs: corrected_mzs,
        intensities: spectrum.intensities.clone(),
    }
}

/// Get calibrated tolerance for fragment matching
///
/// Returns 3×SD from calibration if calibrated, otherwise returns base tolerance.
/// This provides the appropriate tolerance for ppm-based fragment matching
/// after calibration has been applied.
///
/// # Arguments
/// * `calibration` - MS2 calibration parameters
/// * `base_tolerance_ppm` - Default tolerance to use if not calibrated
///
/// # Returns
/// Tolerance in ppm to use for fragment matching
pub fn calibrated_tolerance_ppm(calibration: &MzCalibration, base_tolerance_ppm: f64) -> f64 {
    if calibration.calibrated {
        // Use 3×SD as the tolerance (mean offset is already corrected)
        let tolerance_3sd = 3.0 * calibration.sd;
        // Ensure we have at least a minimum tolerance
        tolerance_3sd.max(1.0)
    } else {
        base_tolerance_ppm
    }
}

/// Get calibrated tolerance for fragment matching (unit-aware)
///
/// Returns 3×SD from calibration in the appropriate unit (ppm or Th).
/// This is the preferred function for unit resolution data where calibration
/// is done in Th instead of ppm.
///
/// # Arguments
/// * `calibration` - MS2 calibration parameters (with unit in `calibration.unit`)
/// * `base_tolerance` - Default tolerance to use if not calibrated
/// * `base_unit` - Unit for the base tolerance
///
/// # Returns
/// Tuple of (tolerance value, unit) for fragment matching
pub fn calibrated_tolerance(
    calibration: &MzCalibration,
    base_tolerance: f64,
    base_unit: ToleranceUnit,
) -> (f64, ToleranceUnit) {
    if calibration.calibrated {
        // Use 3×SD as the tolerance (mean offset is already corrected)
        let tolerance_3sd = 3.0 * calibration.sd;

        // Determine unit from calibration
        let unit = if calibration.unit == "Th" {
            ToleranceUnit::Mz
        } else {
            ToleranceUnit::Ppm
        };

        // Apply appropriate minimum based on unit
        let min_tolerance = if unit == ToleranceUnit::Mz {
            0.05 // Minimum 0.05 Th for unit resolution
        } else {
            1.0 // Minimum 1.0 ppm for HRAM
        };

        (tolerance_3sd.max(min_tolerance), unit)
    } else {
        (base_tolerance, base_unit)
    }
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
        let mut qc_data = MzQCData::new(ToleranceUnit::Ppm);
        for &e in &[-2.5, -2.3, -2.7, -2.4, -2.6] {
            qc_data.add_ms1_error(e);
        }
        for &e in &[1.2, 1.1, 1.3, 1.0, 1.4] {
            qc_data.add_ms2_error(e);
        }

        let (ms1_cal, ms2_cal) = calculate_mz_calibration(&qc_data);

        // MS1 mean should be -2.5, median should be -2.5
        assert!((ms1_cal.mean + 2.5).abs() < 0.01);
        assert!((ms1_cal.median + 2.5).abs() < 0.01);
        assert!(ms1_cal.sd > 0.0);
        assert_eq!(ms1_cal.count, 5);
        assert!(ms1_cal.calibrated);
        assert!(ms1_cal.histogram.is_some());
        assert_eq!(ms1_cal.unit, "ppm");

        // MS2 mean should be 1.2, median should be 1.2
        assert!((ms2_cal.mean - 1.2).abs() < 0.01);
        assert!((ms2_cal.median - 1.2).abs() < 0.01);
        assert!(ms2_cal.sd > 0.0);
        assert_eq!(ms2_cal.count, 5);
        assert!(ms2_cal.calibrated);
        assert!(ms2_cal.histogram.is_some());
        assert_eq!(ms2_cal.unit, "ppm");
    }

    #[test]
    fn test_empty_qc_data() {
        let qc_data = MzQCData::new(ToleranceUnit::Ppm);

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
        let mut qc_data = MzQCData::new(ToleranceUnit::Ppm);
        for &e in &[-2.5, -2.3, -2.7, -1.8, -3.2, 0.5, -2.1] {
            qc_data.add_ms1_error(e);
        }

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
        let mut qc_data = MzQCData::new(ToleranceUnit::Ppm);
        for &e in &[1.0, 2.0, 3.0, 4.0, 5.0] {
            qc_data.add_ms1_error(e);
        }
        let (ms1_cal, _) = calculate_mz_calibration(&qc_data);
        assert!((ms1_cal.median - 3.0).abs() < 0.001);

        // Even number of elements
        let mut qc_data = MzQCData::new(ToleranceUnit::Ppm);
        for &e in &[1.0, 2.0, 3.0, 4.0] {
            qc_data.add_ms1_error(e);
        }
        let (ms1_cal, _) = calculate_mz_calibration(&qc_data);
        assert!((ms1_cal.median - 2.5).abs() < 0.001);
    }

    #[test]
    fn test_unit_resolution_calibration() {
        // Test with Th units for unit resolution instruments
        let mut qc_data = MzQCData::new(ToleranceUnit::Mz);
        for &e in &[-0.01, 0.02, -0.005, 0.015, 0.0] {
            qc_data.add_ms1_error(e);
            qc_data.add_ms2_error(e * 1.1); // Slightly different for MS2
        }

        let (ms1_cal, ms2_cal) = calculate_mz_calibration(&qc_data);

        // Both should be in Th units
        assert_eq!(ms1_cal.unit, "Th");
        assert_eq!(ms2_cal.unit, "Th");

        // Both should be calibrated
        assert!(ms1_cal.calibrated);
        assert!(ms2_cal.calibrated);
    }
}
