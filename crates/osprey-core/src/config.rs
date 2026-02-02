//! Configuration structures for Osprey analysis
//!
//! This module provides all configuration options for the analysis pipeline,
//! including input/output settings, resolution modes, and algorithm parameters.
//!
//! Configuration can be loaded from YAML files for reproducibility.

use crate::{OspreyError, Result};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};

/// Main configuration structure for Osprey analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OspreyConfig {
    // Input/Output
    /// Input mzML file paths
    pub input_files: Vec<PathBuf>,
    /// Spectral library source
    pub library_source: LibrarySource,
    /// Primary output: blib for Skyline
    pub output_blib: PathBuf,
    /// Optional: TSV report
    pub output_report: Option<PathBuf>,

    // Resolution settings
    /// Resolution mode for binning
    pub resolution_mode: ResolutionMode,
    /// Custom bin width override
    pub custom_bin_width: Option<f64>,

    // Candidate selection
    /// Maximum candidates per spectrum
    pub max_candidates_per_spectrum: usize,

    // RT Calibration
    /// RT calibration configuration
    pub rt_calibration: RTCalibrationConfig,

    // Regression
    /// Regularization parameter selection
    pub regularization_lambda: RegularizationSetting,
    /// Maximum iterations for solver
    pub max_iterations: usize,
    /// Convergence threshold
    pub convergence_threshold: f64,

    // Background correction
    /// RT offsets for background estimation
    pub background_rt_offsets: Vec<f64>,

    // Two-step search
    /// Two-step search configuration
    pub two_step_search: TwoStepConfig,

    // FDR control
    /// Run-level FDR threshold
    pub run_fdr: f64,
    /// Experiment-level FDR threshold
    pub experiment_fdr: f64,
    /// Decoy generation method
    pub decoy_method: DecoyMethod,
    /// Whether library already contains decoys
    pub decoys_in_library: bool,

    // Performance
    /// Number of threads to use
    pub n_threads: usize,
    /// Memory limit in GB
    pub memory_limit_gb: Option<f64>,

    // Output options
    /// Include coefficient time series in blib
    pub export_coefficients: bool,
    /// Write features to separate TSV
    pub export_features: bool,
}

impl Default for OspreyConfig {
    fn default() -> Self {
        Self {
            input_files: Vec::new(),
            library_source: LibrarySource::DiannTsv(PathBuf::new()),
            output_blib: PathBuf::from("results.blib"),
            output_report: None,
            resolution_mode: ResolutionMode::Auto,
            custom_bin_width: None,
            max_candidates_per_spectrum: 200,
            rt_calibration: RTCalibrationConfig::default(),
            regularization_lambda: RegularizationSetting::CrossValidated,
            max_iterations: 1000,
            convergence_threshold: 1e-6,
            background_rt_offsets: vec![5.0, 8.0],
            two_step_search: TwoStepConfig::default(),
            run_fdr: 0.01,
            experiment_fdr: 0.01,
            decoy_method: DecoyMethod::Reverse,
            decoys_in_library: false,
            n_threads: num_cpus(),
            memory_limit_gb: None,
            export_coefficients: false,
            export_features: false,
        }
    }
}

/// Get the number of CPUs available
fn num_cpus() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

impl OspreyConfig {
    /// Load configuration from a YAML file
    ///
    /// # Example
    ///
    /// ```yaml
    /// # osprey_config.yaml
    /// input_files:
    ///   - sample1.mzML
    ///   - sample2.mzML
    /// library_source:
    ///   DiannTsv: library.tsv
    /// output_blib: results.blib
    /// rt_tolerance: 2.0
    /// run_fdr: 0.01
    /// ```
    pub fn from_yaml<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let content = fs::read_to_string(path).map_err(|e| {
            OspreyError::ConfigError(format!(
                "Failed to read config file '{}': {}",
                path.display(),
                e
            ))
        })?;

        serde_yaml::from_str(&content).map_err(|e| {
            OspreyError::ConfigError(format!(
                "Failed to parse config file '{}': {}",
                path.display(),
                e
            ))
        })
    }

    /// Save configuration to a YAML file
    pub fn to_yaml<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let content = serde_yaml::to_string(self).map_err(|e| {
            OspreyError::ConfigError(format!("Failed to serialize config: {}", e))
        })?;

        fs::write(path.as_ref(), content).map_err(|e| {
            OspreyError::ConfigError(format!(
                "Failed to write config file '{}': {}",
                path.as_ref().display(),
                e
            ))
        })
    }

    /// Convert configuration to YAML string
    pub fn to_yaml_string(&self) -> Result<String> {
        serde_yaml::to_string(self).map_err(|e| {
            OspreyError::ConfigError(format!("Failed to serialize config: {}", e))
        })
    }

    /// Create a template configuration file with comments
    pub fn create_template<P: AsRef<Path>>(path: P) -> Result<()> {
        let template = r#"# Osprey Configuration File
# See documentation for full options

# Input/Output
input_files:
  - sample1.mzML
  - sample2.mzML
  # Add more files as needed

library_source:
  DiannTsv: library.tsv
  # Or: Blib: library.blib
  # Or: Elib: library.elib

output_blib: results.blib
output_report: results.tsv  # Optional TSV report

# Resolution mode
resolution_mode: Auto
# Options: Auto, UnitResolution, or HRAM with tolerance_ppm
# resolution_mode:
#   HRAM:
#     tolerance_ppm: 20.0

# Candidate selection (precursor filtering uses isolation window from mzML)
max_candidates_per_spectrum: 200

# RT Calibration (multi-file strategy)
# First file: calibrates with ALL peptides using wide tolerance (25% of RT range)
# Subsequent files: reuse calibration from first file
rt_calibration:
  enabled: true
  loess_bandwidth: 0.3            # Fraction of data for local fits (0.2-0.5)
  min_calibration_points: 50      # Minimum detections required
  rt_tolerance_factor: 3.0        # Multiplier for residual SD (for calibrated search)
  initial_tolerance_fraction: 0.25  # Initial tolerance as fraction of RT range (first file)

# Regression
regularization_lambda: CrossValidated
# Options: CrossValidated, Adaptive, or Fixed with value
# regularization_lambda:
#   Fixed: 0.5

# Two-step search (recommended)
two_step_search:
  enabled: true
  step1_fdr: 0.01
  min_runs_for_step2: 1

# FDR control
run_fdr: 0.01
experiment_fdr: 0.01
decoy_method: Reverse  # Options: Reverse, Shuffle, FromLibrary
decoys_in_library: false

# Performance
n_threads: 0  # 0 = auto-detect
# memory_limit_gb: 16.0  # Optional memory limit

# Output options
export_coefficients: false  # Include coefficient time series
export_features: false      # Write features to TSV
"#;

        fs::write(path.as_ref(), template).map_err(|e| {
            OspreyError::ConfigError(format!(
                "Failed to write template config: {}",
                e
            ))
        })
    }

    /// Merge command-line arguments over file configuration
    ///
    /// This allows CLI args to override YAML config values
    pub fn merge_with_args(&mut self, args: &ConfigOverrides) {
        if let Some(ref files) = args.input_files {
            self.input_files = files.clone();
        }
        if let Some(ref library) = args.library {
            self.library_source = LibrarySource::from_path(library.clone());
        }
        if let Some(ref output) = args.output {
            self.output_blib = output.clone();
        }
        if let Some(ref report) = args.report {
            self.output_report = Some(report.clone());
        }
        if let Some(rt_tol) = args.rt_tolerance {
            self.rt_calibration.fallback_rt_tolerance = rt_tol;
        }
        if args.disable_rt_calibration {
            self.rt_calibration.enabled = false;
        }
        if let Some(fdr) = args.run_fdr {
            self.run_fdr = fdr;
        }
        if let Some(threads) = args.n_threads {
            self.n_threads = threads;
        }
        if let Some(lambda) = args.lambda {
            self.regularization_lambda = RegularizationSetting::Fixed(lambda);
        }
        if args.verbose {
            // Logging level would be handled separately
        }
    }
}

/// Command-line argument overrides for config
#[derive(Debug, Default, Clone)]
pub struct ConfigOverrides {
    pub input_files: Option<Vec<PathBuf>>,
    pub library: Option<PathBuf>,
    pub output: Option<PathBuf>,
    pub report: Option<PathBuf>,
    pub rt_tolerance: Option<f64>,
    pub run_fdr: Option<f64>,
    pub n_threads: Option<usize>,
    pub lambda: Option<f64>,
    pub verbose: bool,
    pub disable_rt_calibration: bool,
}

/// Two-step search configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TwoStepConfig {
    /// Enable two-step search
    pub enabled: bool,
    /// Step 1 FDR threshold
    pub step1_fdr: f64,
    /// Minimum runs for peptide inclusion in Step 2
    pub min_runs_for_step2: usize,
    /// Include library peptides even if not detected
    pub include_high_confidence: bool,
}

impl Default for TwoStepConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            step1_fdr: 0.01,
            min_runs_for_step2: 1,
            include_high_confidence: false,
        }
    }
}

impl TwoStepConfig {
    /// Create enabled two-step configuration
    pub fn enabled() -> Self {
        Self::default()
    }

    /// Create disabled two-step configuration
    pub fn disabled() -> Self {
        Self {
            enabled: false,
            ..Self::default()
        }
    }
}

/// RT Calibration configuration
///
/// Controls how library retention times are calibrated against measured RTs.
/// Uses LOESS (Local Polynomial Regression) for calibration.
///
/// ## Multi-File Strategy
///
/// - **First file**: Uses ALL library peptides with wide tolerance (initial_tolerance_fraction × RT range)
/// - **Subsequent files**: Reuse calibration from first file with tight tolerance (residual SD × factor)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RTCalibrationConfig {
    /// Enable RT calibration
    pub enabled: bool,
    /// LOESS bandwidth (0.0-1.0, fraction of data to use for local fits)
    pub loess_bandwidth: f64,
    /// Minimum detections required for calibration
    pub min_calibration_points: usize,
    /// RT tolerance multiplier (× residual SD) for calibrated search
    pub rt_tolerance_factor: f64,
    /// Fallback RT tolerance (minutes) if calibration fails
    pub fallback_rt_tolerance: f64,
    /// Initial tolerance as fraction of library RT range for first file calibration (default: 0.25 = 25%)
    pub initial_tolerance_fraction: f64,
}

impl Default for RTCalibrationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            loess_bandwidth: 0.3,
            min_calibration_points: 50,
            rt_tolerance_factor: 3.0,
            fallback_rt_tolerance: 2.0,
            initial_tolerance_fraction: 0.25, // 25% of RT range
        }
    }
}

impl RTCalibrationConfig {
    /// Create enabled calibration with default settings
    pub fn enabled() -> Self {
        Self::default()
    }

    /// Create disabled calibration (uses fallback_rt_tolerance)
    pub fn disabled() -> Self {
        Self {
            enabled: false,
            ..Self::default()
        }
    }

    /// Set the fallback RT tolerance
    pub fn with_fallback_tolerance(mut self, tolerance: f64) -> Self {
        self.fallback_rt_tolerance = tolerance;
        self
    }
}

/// Resolution mode for spectrum binning
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ResolutionMode {
    /// Unit resolution (~1 Th bins)
    UnitResolution,
    /// High-resolution accurate mass
    HRAM {
        /// Tolerance in ppm
        tolerance_ppm: f64,
    },
    /// Auto-detect from data
    Auto,
}

impl Default for ResolutionMode {
    fn default() -> Self {
        ResolutionMode::Auto
    }
}

/// Regularization parameter selection strategy
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RegularizationSetting {
    /// Fixed lambda value
    Fixed(f64),
    /// Select via cross-validation
    CrossValidated,
    /// Adaptive per spectrum
    Adaptive,
}

impl Default for RegularizationSetting {
    fn default() -> Self {
        RegularizationSetting::CrossValidated
    }
}

/// Decoy generation method
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DecoyMethod {
    /// Reverse peptide sequence
    Reverse,
    /// Shuffle peptide sequence
    Shuffle,
    /// Use decoys already in library
    FromLibrary,
}

impl Default for DecoyMethod {
    fn default() -> Self {
        DecoyMethod::Reverse
    }
}

/// Spectral library source
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LibrarySource {
    /// DIA-NN TSV format library
    DiannTsv(PathBuf),
    /// Skyline BiblioSpec library (.blib)
    Blib(PathBuf),
    /// EncyclopeDIA chromatogram library (.elib)
    Elib(PathBuf),
    /// Skyline document (extracts library from document)
    SkylineDocument(PathBuf),
}

impl LibrarySource {
    /// Get the file path for the library source
    pub fn path(&self) -> &PathBuf {
        match self {
            LibrarySource::DiannTsv(p) => p,
            LibrarySource::Blib(p) => p,
            LibrarySource::Elib(p) => p,
            LibrarySource::SkylineDocument(p) => p,
        }
    }

    /// Detect library format from file extension
    pub fn from_path(path: PathBuf) -> Self {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .map(|e| e.to_lowercase());

        match ext.as_deref() {
            Some("blib") => LibrarySource::Blib(path),
            Some("elib") => LibrarySource::Elib(path),
            Some("sky") => LibrarySource::SkylineDocument(path),
            _ => LibrarySource::DiannTsv(path), // Default to TSV
        }
    }
}

/// Library metadata
#[derive(Debug, Clone)]
pub struct LibraryInfo {
    /// Library source
    pub source: LibrarySource,
    /// Number of target entries
    pub n_targets: usize,
    /// Number of decoy entries
    pub n_decoys: usize,
    /// Retention time type
    pub rt_type: RtType,
    /// Unique modifications in library
    pub modifications: Vec<String>,
    /// Proteome coverage if known
    pub proteome_coverage: Option<f64>,
}

/// Retention time type in library
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum RtType {
    /// From predictor, not calibrated
    Predicted,
    /// iRT or similar normalized scale
    Normalized,
    /// Calibrated to specific LC conditions
    CalibratedMinutes,
}

/// Binning configuration
#[derive(Debug, Clone, Copy)]
pub struct BinConfig {
    /// Bin width in Th
    pub bin_width: f64,
    /// Bin offset in Th
    pub bin_offset: f64,
    /// Minimum m/z
    pub min_mz: f64,
    /// Maximum m/z
    pub max_mz: f64,
    /// Total number of bins
    pub n_bins: usize,
}

impl BinConfig {
    /// Create unit resolution binning configuration (Comet-style)
    pub fn unit_resolution() -> Self {
        Self {
            bin_width: 1.0005079,
            bin_offset: 0.4,
            min_mz: 100.0,
            max_mz: 2000.0,
            n_bins: 1899,
        }
    }

    /// Create HRAM binning configuration
    pub fn hram(tolerance_th: f64) -> Self {
        let min_mz = 100.0;
        let max_mz = 2000.0;
        Self {
            bin_width: tolerance_th,
            bin_offset: 0.0,
            min_mz,
            max_mz,
            n_bins: ((max_mz - min_mz) / tolerance_th).ceil() as usize,
        }
    }

    /// Create HRAM binning configuration from ppm tolerance
    pub fn hram_ppm(tolerance_ppm: f64) -> Self {
        // Use tolerance at a reference m/z (e.g., 500 Th)
        let ref_mz = 500.0;
        let tolerance_th = ref_mz * tolerance_ppm / 1_000_000.0;
        Self::hram(tolerance_th)
    }

    /// Convert m/z to bin index
    pub fn mz_to_bin(&self, mz: f64) -> Option<usize> {
        if mz < self.min_mz || mz > self.max_mz {
            return None;
        }
        let bin = ((mz - self.min_mz + self.bin_offset) / self.bin_width).floor() as usize;
        if bin < self.n_bins {
            Some(bin)
        } else {
            None
        }
    }

    /// Convert bin index to m/z (bin center)
    pub fn bin_to_mz(&self, bin: usize) -> f64 {
        self.min_mz - self.bin_offset + (bin as f64 + 0.5) * self.bin_width
    }
}

impl Default for BinConfig {
    fn default() -> Self {
        Self::unit_resolution()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bin_config_unit_resolution() {
        let config = BinConfig::unit_resolution();
        assert!((config.bin_width - 1.0005079).abs() < 1e-6);

        // Test bin assignment
        let bin = config.mz_to_bin(500.0).unwrap();
        assert!(bin > 0 && bin < config.n_bins);

        // Test round-trip approximately
        let mz_back = config.bin_to_mz(bin);
        assert!((mz_back - 500.0).abs() < config.bin_width);
    }

    #[test]
    fn test_library_source_from_path() {
        let tsv = LibrarySource::from_path(PathBuf::from("library.tsv"));
        assert!(matches!(tsv, LibrarySource::DiannTsv(_)));

        let blib = LibrarySource::from_path(PathBuf::from("library.blib"));
        assert!(matches!(blib, LibrarySource::Blib(_)));

        let elib = LibrarySource::from_path(PathBuf::from("library.elib"));
        assert!(matches!(elib, LibrarySource::Elib(_)));
    }

    #[test]
    fn test_default_config() {
        let config = OspreyConfig::default();
        assert_eq!(config.rt_calibration.fallback_rt_tolerance, 2.0);
        assert_eq!(config.run_fdr, 0.01);
        assert!(config.two_step_search.enabled);
        assert!(config.rt_calibration.enabled);
    }

    #[test]
    fn test_yaml_roundtrip() {
        let config = OspreyConfig::default();
        let yaml = config.to_yaml_string().unwrap();

        // Verify it's valid YAML that can be parsed back
        let parsed: OspreyConfig = serde_yaml::from_str(&yaml).unwrap();
        assert_eq!(parsed.rt_calibration.fallback_rt_tolerance, config.rt_calibration.fallback_rt_tolerance);
        assert_eq!(parsed.run_fdr, config.run_fdr);
    }

    #[test]
    fn test_config_merge() {
        let mut config = OspreyConfig::default();
        let overrides = ConfigOverrides {
            rt_tolerance: Some(3.0),
            run_fdr: Some(0.05),
            ..Default::default()
        };

        config.merge_with_args(&overrides);
        assert_eq!(config.rt_calibration.fallback_rt_tolerance, 3.0);
        assert_eq!(config.run_fdr, 0.05);
    }

    #[test]
    fn test_rt_calibration_config() {
        let config = RTCalibrationConfig::default();
        assert!(config.enabled);
        assert_eq!(config.loess_bandwidth, 0.3);
        assert_eq!(config.initial_tolerance_fraction, 0.25);

        let disabled = RTCalibrationConfig::disabled();
        assert!(!disabled.enabled);
    }
}
