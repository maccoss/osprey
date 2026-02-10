//! Calibration parameter serialization and deserialization
//!
//! Saves/loads calibration parameters to/from JSON files for reuse across searches.

use super::CalibrationParams;
use osprey_core::Result;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

/// Save calibration parameters to JSON file
///
/// # Arguments
/// * `calibration` - Calibration parameters to save
/// * `output_path` - Path to output JSON file
///
/// # Returns
/// Ok(()) on success
///
/// # Example
/// ```no_run
/// use osprey_chromatography::calibration::{CalibrationParams, save_calibration};
/// use std::path::Path;
///
/// # let calibration: CalibrationParams = unimplemented!();
/// save_calibration(&calibration, Path::new("output.calibration.json"))
///     .expect("Failed to save calibration");
/// ```
pub fn save_calibration(
    calibration: &CalibrationParams,
    output_path: &Path,
) -> Result<()> {
    let file = File::create(output_path)
        .map_err(|e| osprey_core::OspreyError::OutputError(format!(
            "Failed to create calibration file {}: {}",
            output_path.display(),
            e
        )))?;

    let writer = BufWriter::new(file);

    serde_json::to_writer_pretty(writer, calibration)
        .map_err(|e| osprey_core::OspreyError::OutputError(format!(
            "Failed to serialize calibration parameters: {}",
            e
        )))?;

    log::info!("Saved calibration to: {}", output_path.display());

    Ok(())
}

/// Load calibration parameters from JSON file
///
/// # Arguments
/// * `input_path` - Path to calibration JSON file
///
/// # Returns
/// Loaded calibration parameters
///
/// # Example
/// ```no_run
/// use osprey_chromatography::calibration::load_calibration;
/// use std::path::Path;
///
/// let calibration = load_calibration(Path::new("output.calibration.json"))
///     .expect("Failed to load calibration");
///
/// println!("MS1 offset: {:.2} ppm", calibration.ms1_calibration.mean);
/// ```
pub fn load_calibration(input_path: &Path) -> Result<CalibrationParams> {
    let file = File::open(input_path)
        .map_err(|e| osprey_core::OspreyError::FileNotFound(format!(
            "Failed to open calibration file {}: {}",
            input_path.display(),
            e
        )))?;

    let reader = BufReader::new(file);

    let calibration: CalibrationParams = serde_json::from_reader(reader)
        .map_err(|e| osprey_core::OspreyError::ConfigError(format!(
            "Failed to deserialize calibration parameters: {}",
            e
        )))?;

    log::info!("Loaded calibration from: {}", input_path.display());

    Ok(calibration)
}

/// Generate calibration filename from base output name
///
/// # Arguments
/// * `output_base` - Base name for output files (e.g., "results")
///
/// # Returns
/// Calibration filename (e.g., "results.calibration.json")
///
/// # Example
/// ```
/// use osprey_chromatography::calibration::calibration_filename;
///
/// let filename = calibration_filename("my_search");
/// assert_eq!(filename, "my_search.calibration.json");
/// ```
pub fn calibration_filename(output_base: &str) -> String {
    format!("{}.calibration.json", output_base)
}

/// Generate calibration filename from input mzML file path
///
/// The calibration file is named after the input file, not the output,
/// so that calibration can be reused when reprocessing the same data file.
///
/// # Arguments
/// * `input_path` - Path to the input mzML file
///
/// # Returns
/// Calibration filename based on input file (e.g., "sample.mzML" -> "sample.calibration.json")
///
/// # Example
/// ```
/// use osprey_chromatography::calibration::calibration_filename_for_input;
/// use std::path::Path;
///
/// let filename = calibration_filename_for_input(Path::new("/data/sample.mzML"));
/// assert_eq!(filename, "sample.calibration.json");
///
/// let filename = calibration_filename_for_input(Path::new("test.dia.mzML"));
/// assert_eq!(filename, "test.dia.calibration.json");
/// ```
pub fn calibration_filename_for_input(input_path: &Path) -> String {
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    format!("{}.calibration.json", stem)
}

/// Get the full path for a calibration file based on input mzML and output directory
///
/// # Arguments
/// * `input_path` - Path to the input mzML file
/// * `output_dir` - Directory where calibration file should be stored
///
/// # Returns
/// Full path to calibration file
pub fn calibration_path_for_input(input_path: &Path, output_dir: &Path) -> std::path::PathBuf {
    let filename = calibration_filename_for_input(input_path);
    output_dir.join(filename)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::calibration::{CalibrationMetadata, MzCalibration, RTCalibrationParams, RTCalibrationMethod};

    #[allow(dead_code)]
    fn create_test_calibration() -> CalibrationParams {
        CalibrationParams {
            metadata: CalibrationMetadata {
                num_confident_peptides: 150,
                num_sampled_precursors: 2000,
                calibration_successful: true,
                timestamp: "2024-01-15T10:30:00Z".to_string(),
                isolation_scheme: None,
            },
            ms1_calibration: MzCalibration {
                mean: -2.5,
                median: -2.4,
                sd: 0.8,
                count: 150,
                unit: "ppm".to_string(),
                adjusted_tolerance: Some(4.9),
                window_halfwidth_multiplier: Some(3.0),
                histogram: None,
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
                n_points: 150,
                r_squared: 0.98,
                model_params: None,
            },
        }
    }

    /// Verifies that calibration_filename appends ".calibration.json" to the base name.
    #[test]
    fn test_calibration_filename() {
        assert_eq!(calibration_filename("results"), "results.calibration.json");
        assert_eq!(calibration_filename("my_search"), "my_search.calibration.json");
        assert_eq!(calibration_filename("test"), "test.calibration.json");
    }

    /// Verifies calibration filename derivation from input mzML paths, including multi-extension filenames.
    #[test]
    fn test_calibration_filename_for_input() {
        use std::path::Path;

        // Basic mzML file
        assert_eq!(
            calibration_filename_for_input(Path::new("/data/sample.mzML")),
            "sample.calibration.json"
        );

        // File with multiple extensions (e.g., sample.dia.mzML)
        assert_eq!(
            calibration_filename_for_input(Path::new("test.dia.mzML")),
            "test.dia.calibration.json"
        );

        // Simple filename only (no path)
        assert_eq!(
            calibration_filename_for_input(Path::new("experiment.mzML")),
            "experiment.calibration.json"
        );
    }

    /// Verifies that calibration_path_for_input joins the output directory with the input-derived calibration filename.
    #[test]
    fn test_calibration_path_for_input() {
        use std::path::Path;

        let input = Path::new("/data/raw/sample.mzML");
        let output_dir = Path::new("/output/results");

        let cal_path = calibration_path_for_input(input, output_dir);
        assert_eq!(
            cal_path,
            Path::new("/output/results/sample.calibration.json")
        );
    }

    // Note: File I/O tests require tempfile crate which isn't available here
    // These would be integration tests in a real setup
}
