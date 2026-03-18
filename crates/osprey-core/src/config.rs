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
    /// Fragment tolerance for LibCosine scoring (ppm-based matching, no binning)
    pub fragment_tolerance: FragmentToleranceConfig,
    /// Precursor tolerance for MS1 matching
    pub precursor_tolerance: FragmentToleranceConfig,

    // Candidate selection
    /// Maximum candidates per spectrum
    pub max_candidates_per_spectrum: usize,

    // RT Calibration
    /// RT calibration configuration
    pub rt_calibration: RTCalibrationConfig,

    // FDR control
    /// Run-level FDR threshold
    pub run_fdr: f64,
    /// Experiment-level FDR threshold
    pub experiment_fdr: f64,
    /// Decoy generation method
    pub decoy_method: DecoyMethod,
    /// Whether library already contains decoys
    pub decoys_in_library: bool,

    // FDR method
    /// FDR method: native Percolator (default), external mokapot, or simple target-decoy
    #[serde(default)]
    pub fdr_method: FdrMethod,

    /// Write PIN files for external tools (default: false)
    #[serde(default)]
    pub write_pin: bool,

    // Inter-replicate reconciliation
    /// Inter-replicate peak reconciliation settings
    #[serde(default)]
    pub reconciliation: ReconciliationConfig,

    // Search pre-filter
    /// Enable the coelution signal pre-filter (3-of-4 consecutive scans with ≥2 top-6 fragments).
    /// Speeds up searches ~30% with minimal sensitivity loss. Disable with --no-prefilter.
    #[serde(default = "default_true")]
    pub prefilter_enabled: bool,

    // Performance
    /// Number of threads to use
    pub n_threads: usize,
    /// Memory limit in GB
    pub memory_limit_gb: Option<f64>,
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
            fragment_tolerance: FragmentToleranceConfig::default(), // 10 ppm for HRAM
            precursor_tolerance: FragmentToleranceConfig::hram(10.0), // 10 ppm for precursor
            max_candidates_per_spectrum: 5250,
            rt_calibration: RTCalibrationConfig::default(),
            fdr_method: FdrMethod::default(),
            write_pin: false,
            run_fdr: 0.01,
            experiment_fdr: 0.01,
            decoy_method: DecoyMethod::Reverse,
            decoys_in_library: false,
            reconciliation: ReconciliationConfig::default(),
            prefilter_enabled: true,
            n_threads: num_cpus(),
            memory_limit_gb: None,
        }
    }
}

fn default_true() -> bool {
    true
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
        let content = serde_yaml::to_string(self)
            .map_err(|e| OspreyError::ConfigError(format!("Failed to serialize config: {}", e)))?;

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
        serde_yaml::to_string(self)
            .map_err(|e| OspreyError::ConfigError(format!("Failed to serialize config: {}", e)))
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
# Options: Auto, UnitResolution, or HRAM
# Unit resolution: 1.0005079 Th bins, 0.4 offset (for Stellar/low-res instruments)
# HRAM: 0.02 Th bins, 0 offset (for Astral/high-res instruments)

# Fragment tolerance for LibCosine scoring (ppm-based peak matching)
# This is used for calibration scoring and RT/m/z calibration
fragment_tolerance:
  tolerance: 10.0   # Tolerance value
  unit: Ppm         # Unit: Ppm or Da (use Ppm for HRAM, Da for unit resolution)

# Precursor tolerance for MS1 matching
precursor_tolerance:
  tolerance: 10.0
  unit: Ppm

# Candidate selection (precursor filtering uses isolation window from mzML)
max_candidates_per_spectrum: 500

# RT Calibration
# Samples a subset of peptides for fast calibration
rt_calibration:
  enabled: true
  loess_bandwidth: 0.3            # Fraction of data for local fits (0.2-0.5)
  min_calibration_points: 200     # Minimum detections required for robust LOESS fit
  rt_tolerance_factor: 3.0        # Multiplier for residual SD (for calibrated search)
  initial_tolerance_fraction: 1.0   # Initial tolerance as fraction of RT range (1.0 = no RT filtering during calibration)
  calibration_sample_size: 100000 # Target peptides to sample for calibration (0 = all)
  calibration_retry_factor: 2.0   # Multiply sample size on retry if too few calibration points

# FDR control
run_fdr: 0.01
experiment_fdr: 0.01
decoy_method: Reverse  # Options: Reverse, Shuffle, FromLibrary
decoys_in_library: false
fdr_method: Percolator  # Options: Percolator (native SVM), Mokapot (external Python), Simple (no ML)
write_pin: false  # Write PIN files for external tools

# Performance
n_threads: 0  # 0 = auto-detect
# memory_limit_gb: 16.0  # Optional memory limit

"#;

        fs::write(path.as_ref(), template).map_err(|e| {
            OspreyError::ConfigError(format!("Failed to write template config: {}", e))
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
        if args.verbose {
            // Logging level would be handled separately
        }
        if let Some(tol) = args.fragment_tolerance {
            self.fragment_tolerance.tolerance = tol;
        }
        if let Some(unit) = args.fragment_unit {
            self.fragment_tolerance.unit = unit;
        }
        if let Some(tol) = args.precursor_tolerance {
            self.precursor_tolerance.tolerance = tol;
        }
        if let Some(unit) = args.precursor_unit {
            self.precursor_tolerance.unit = unit;
        }
        if let Some(method) = args.fdr_method {
            self.fdr_method = method;
        }
        if args.write_pin {
            self.write_pin = true;
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
    pub verbose: bool,
    pub disable_rt_calibration: bool,
    pub fragment_tolerance: Option<f64>,
    pub fragment_unit: Option<ToleranceUnit>,
    pub precursor_tolerance: Option<f64>,
    pub precursor_unit: Option<ToleranceUnit>,
    pub fdr_method: Option<FdrMethod>,
    pub write_pin: bool,
}

/// RT Calibration configuration
///
/// Controls how library retention times are calibrated against measured RTs.
/// Uses LOESS (Local Polynomial Regression) for calibration.
///
/// ## Calibration Sampling Strategy (pyXcorrDIA-style)
///
/// For large libraries (millions of entries), calibration samples a subset:
/// - **Attempt 1: 100,000 peptides** sampled for calibration scoring
/// - **Attempt 2: 200,000 peptides** (2× retry) if first attempt doesn't yield enough calibration points
/// - **Attempt 3: ALL peptides** as final fallback
/// - FDR is cumulative across attempts (matches accumulate)
///
/// ## Multi-File Strategy
///
/// - **First file**: Sample peptides with wide RT tolerance (initial_tolerance_fraction × RT range)
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
    /// Minimum RT tolerance (minutes) floor for local tolerance calculation
    /// This prevents over-filtering in regions with very tight calibration
    #[serde(default = "default_min_rt_tolerance")]
    pub min_rt_tolerance: f64,
    /// Initial tolerance as fraction of library RT range for first file calibration (default: 1.0 = 100%)
    /// Set to 1.0 to match pyXcorrDIA's approach of no RT filtering during calibration
    pub initial_tolerance_fraction: f64,
    /// Number of target peptides to sample for calibration (0 = use all). Default: 100,000.
    /// Sampling uses 2D stratification across RT and m/z for uniform coverage.
    #[serde(default = "default_calibration_sample_size")]
    pub calibration_sample_size: usize,
    /// Multiplier for expanding sample on retry if too few calibration points. Default: 2.0.
    /// On each retry, sample_size is multiplied by this factor. Final attempt uses all entries.
    #[serde(default = "default_calibration_retry_factor")]
    pub calibration_retry_factor: f64,
    /// Maximum RT tolerance (minutes). Hard cap to prevent catastrophic widening
    /// when calibration is poor. Default: 3.0 min.
    #[serde(default = "default_max_rt_tolerance")]
    pub max_rt_tolerance: f64,
    /// Use percentile-based RT tolerance instead of SD-based (DIA-NN approach).
    /// When true: tolerance = max(2 × P20_abs_error, RT_range / 40)
    /// When false: tolerance = rt_tolerance_factor × residual_SD
    #[serde(default = "default_use_percentile_tolerance")]
    pub use_percentile_tolerance: bool,
}

fn default_min_rt_tolerance() -> f64 {
    0.5
}

fn default_calibration_sample_size() -> usize {
    100000
}

fn default_calibration_retry_factor() -> f64 {
    2.0
}

fn default_max_rt_tolerance() -> f64 {
    3.0
}

fn default_use_percentile_tolerance() -> bool {
    false
}

impl Default for RTCalibrationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            loess_bandwidth: 0.3,
            min_calibration_points: 200,
            rt_tolerance_factor: 3.0,
            fallback_rt_tolerance: 2.0,
            min_rt_tolerance: 0.5,           // 0.5 minute minimum
            initial_tolerance_fraction: 1.0, // 100% of RT range (no RT filtering, like pyXcorrDIA)
            calibration_sample_size: 100000,
            calibration_retry_factor: 2.0,
            max_rt_tolerance: 3.0,           // Hard cap at 3 minutes
            use_percentile_tolerance: false, // Use 3×SD tolerance (covers >99% of peptides)
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
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize, Default)]
pub enum ResolutionMode {
    /// Unit resolution (~1.0005079 Th bins, 0.4 offset)
    UnitResolution,
    /// High-resolution accurate mass (0.02 Th bins, 0 offset)
    HRAM,
    /// Auto-detect from data
    #[default]
    Auto,
}

/// Decoy generation method
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum DecoyMethod {
    /// Reverse peptide sequence
    #[default]
    Reverse,
    /// Shuffle peptide sequence
    Shuffle,
    /// Use decoys already in library
    FromLibrary,
}

/// FDR control method
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum FdrMethod {
    /// Native Rust Percolator implementation (linear SVM, default)
    #[default]
    Percolator,
    /// External Python mokapot
    Mokapot,
    /// Simple target-decoy competition (no ML rescoring)
    Simple,
}

impl std::fmt::Display for FdrMethod {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FdrMethod::Percolator => write!(f, "percolator"),
            FdrMethod::Mokapot => write!(f, "mokapot"),
            FdrMethod::Simple => write!(f, "simple"),
        }
    }
}

/// Inter-replicate peak reconciliation configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReconciliationConfig {
    /// Enable inter-replicate peak reconciliation (default: true for multi-file)
    pub enabled: bool,
    /// Number of CWT candidate peaks to store per precursor (default: 5)
    pub top_n_peaks: usize,
    /// FDR threshold for selecting consensus peptides (default: 0.05)
    /// More liberal than final FDR to cast a wider net and rescue peptides
    /// with supporting evidence across replicates.
    pub consensus_fdr: f64,
}

impl Default for ReconciliationConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            top_n_peaks: 5,
            consensus_fdr: 0.05,
        }
    }
}

/// Fragment tolerance configuration for LibCosine scoring
///
/// LibCosine uses direct peak matching with a tolerance window, NOT binning.
/// This matches pyXcorrDIA's implementation exactly.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct FragmentToleranceConfig {
    /// Fragment m/z tolerance value
    pub tolerance: f64,
    /// Tolerance unit (ppm or Da)
    pub unit: ToleranceUnit,
}

impl Default for FragmentToleranceConfig {
    fn default() -> Self {
        // Default: 10 ppm for HRAM data (matches pyXcorrDIA default)
        Self {
            tolerance: 10.0,
            unit: ToleranceUnit::Ppm,
        }
    }
}

impl FragmentToleranceConfig {
    /// Create configuration for HRAM data (ppm-based)
    pub fn hram(ppm: f64) -> Self {
        Self {
            tolerance: ppm,
            unit: ToleranceUnit::Ppm,
        }
    }

    /// Create configuration for unit resolution data (Da-based)
    pub fn unit_resolution(da: f64) -> Self {
        Self {
            tolerance: da,
            unit: ToleranceUnit::Mz,
        }
    }

    /// Calculate tolerance in Da for a given m/z
    pub fn tolerance_da(&self, mz: f64) -> f64 {
        match self.unit {
            ToleranceUnit::Ppm => mz * self.tolerance / 1e6,
            ToleranceUnit::Mz => self.tolerance,
        }
    }

    /// Check if a given m/z delta is within tolerance
    pub fn within_tolerance(&self, lib_mz: f64, obs_mz: f64) -> bool {
        let delta = (obs_mz - lib_mz).abs();
        match self.unit {
            ToleranceUnit::Ppm => {
                let ppm_error = delta / lib_mz * 1e6;
                ppm_error <= self.tolerance
            }
            ToleranceUnit::Mz => delta <= self.tolerance,
        }
    }

    /// Calculate mass error in the configured unit
    pub fn mass_error(&self, lib_mz: f64, obs_mz: f64) -> f64 {
        let delta = obs_mz - lib_mz;
        match self.unit {
            ToleranceUnit::Ppm => delta / lib_mz * 1e6,
            ToleranceUnit::Mz => delta,
        }
    }
}

/// Tolerance unit for m/z matching
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ToleranceUnit {
    /// Parts per million (relative to m/z)
    Ppm,
    /// m/z units (Thomson, Th) - absolute tolerance in mass-to-charge
    #[serde(alias = "Da")] // Accept "Da" for backwards compatibility
    Mz,
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

/// Binning configuration using Comet's BIN macro:
///   BIN(mass) = (int)(mass / bin_width + (1 - offset))
///
/// Unit resolution: bin_width=1.0005079, offset=0.4 → BIN(m) = (int)(m/1.0005 + 0.6)
/// HRAM:            bin_width=0.02,      offset=0.0 → BIN(m) = (int)(m/0.02 + 1.0)
#[derive(Debug, Clone, Copy)]
pub struct BinConfig {
    /// Bin width in Th
    pub bin_width: f64,
    /// Bin offset (Comet fragment_bin_offset): 0.4 for unit resolution, 0.0 for HRAM
    pub bin_offset: f64,
    /// Precomputed: 1.0 / bin_width (for fast BIN macro)
    pub inverse_bin_width: f64,
    /// Precomputed: 1.0 - bin_offset (for fast BIN macro)
    pub one_minus_offset: f64,
    /// Maximum m/z
    pub max_mz: f64,
    /// Total number of bins
    pub n_bins: usize,
}

impl BinConfig {
    /// Create unit resolution binning configuration (Comet-style)
    ///
    /// BIN(mass) = (int)(mass / 1.0005079 + 0.6)
    pub fn unit_resolution() -> Self {
        let bin_width = 1.0005079;
        let bin_offset = 0.4;
        let max_mz = 2000.0;
        let inverse_bin_width = 1.0 / bin_width;
        let one_minus_offset = 1.0 - bin_offset;
        let n_bins = (max_mz * inverse_bin_width + one_minus_offset) as usize + 1;
        Self {
            bin_width,
            bin_offset,
            inverse_bin_width,
            one_minus_offset,
            max_mz,
            n_bins,
        }
    }

    /// Create HRAM binning configuration (Comet-style, 0.02 Th bins, 0 offset)
    ///
    /// BIN(mass) = (int)(mass / 0.02 + 1.0)
    pub fn hram() -> Self {
        let bin_width = 0.02;
        let bin_offset = 0.0;
        let max_mz = 2000.0;
        let inverse_bin_width = 1.0 / bin_width;
        let one_minus_offset = 1.0 - bin_offset;
        let n_bins = (max_mz * inverse_bin_width + one_minus_offset) as usize + 1;
        Self {
            bin_width,
            bin_offset,
            inverse_bin_width,
            one_minus_offset,
            max_mz,
            n_bins,
        }
    }

    /// Convert m/z to bin index using Comet BIN macro:
    ///   BIN(mass) = (int)(mass * inverse_bin_width + one_minus_offset)
    #[inline]
    pub fn mz_to_bin(&self, mz: f64) -> Option<usize> {
        if mz < 0.0 || mz > self.max_mz {
            return None;
        }
        let bin = (mz * self.inverse_bin_width + self.one_minus_offset) as usize;
        if bin < self.n_bins {
            Some(bin)
        } else {
            None
        }
    }

    /// Convert bin index to approximate m/z (bin center)
    pub fn bin_to_mz(&self, bin: usize) -> f64 {
        (bin as f64 + 0.5 - self.one_minus_offset) * self.bin_width
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

    /// Verifies that unit resolution BinConfig creates bins with correct width and m/z round-trip.
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

    /// Verifies that LibrarySource auto-detects format from file extension.
    #[test]
    fn test_library_source_from_path() {
        let tsv = LibrarySource::from_path(PathBuf::from("library.tsv"));
        assert!(matches!(tsv, LibrarySource::DiannTsv(_)));

        let blib = LibrarySource::from_path(PathBuf::from("library.blib"));
        assert!(matches!(blib, LibrarySource::Blib(_)));

        let elib = LibrarySource::from_path(PathBuf::from("library.elib"));
        assert!(matches!(elib, LibrarySource::Elib(_)));
    }

    /// Verifies that default OspreyConfig has expected RT tolerance, FDR, and calibration values.
    #[test]
    fn test_default_config() {
        let config = OspreyConfig::default();
        assert_eq!(config.rt_calibration.fallback_rt_tolerance, 2.0);
        assert_eq!(config.run_fdr, 0.01);
        assert!(config.rt_calibration.enabled);
    }

    /// Verifies that YAML serialization and deserialization preserves config field values.
    #[test]
    fn test_yaml_roundtrip() {
        let config = OspreyConfig::default();
        let yaml = config.to_yaml_string().unwrap();

        // Verify it's valid YAML that can be parsed back
        let parsed: OspreyConfig = serde_yaml::from_str(&yaml).unwrap();
        assert_eq!(
            parsed.rt_calibration.fallback_rt_tolerance,
            config.rt_calibration.fallback_rt_tolerance
        );
        assert_eq!(parsed.run_fdr, config.run_fdr);
    }

    /// Verifies that CLI ConfigOverrides correctly override config file values.
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

    /// Verifies RTCalibrationConfig defaults and that disabled() turns off calibration.
    #[test]
    fn test_rt_calibration_config() {
        let config = RTCalibrationConfig::default();
        assert!(config.enabled);
        assert_eq!(config.loess_bandwidth, 0.3);
        assert_eq!(config.initial_tolerance_fraction, 1.0);

        let disabled = RTCalibrationConfig::disabled();
        assert!(!disabled.enabled);
    }

    /// Verifies ppm-based fragment tolerance Da conversion, within_tolerance, and mass_error.
    #[test]
    fn test_fragment_tolerance_ppm() {
        let config = FragmentToleranceConfig::hram(10.0);
        assert_eq!(config.tolerance, 10.0);
        assert_eq!(config.unit, ToleranceUnit::Ppm);

        // Test tolerance_da conversion at 500 m/z
        // 10 ppm of 500 = 0.005 Da
        let tol_da = config.tolerance_da(500.0);
        assert!((tol_da - 0.005).abs() < 1e-10);

        // Test within_tolerance at exactly 10 ppm
        assert!(config.within_tolerance(500.0, 500.005));
        assert!(!config.within_tolerance(500.0, 500.006)); // >10 ppm

        // Test mass_error calculation
        // (500.005 - 500.0) / 500.0 * 1e6 = 10 ppm
        let error = config.mass_error(500.0, 500.005);
        assert!((error - 10.0).abs() < 0.1);
    }

    /// Verifies Da-based fragment tolerance is constant across m/z values.
    #[test]
    fn test_fragment_tolerance_da() {
        let config = FragmentToleranceConfig::unit_resolution(0.3);
        assert_eq!(config.tolerance, 0.3);
        assert_eq!(config.unit, ToleranceUnit::Mz);

        // Test tolerance_da (should be constant)
        assert!((config.tolerance_da(500.0) - 0.3).abs() < 1e-10);
        assert!((config.tolerance_da(1000.0) - 0.3).abs() < 1e-10);

        // Test within_tolerance
        assert!(config.within_tolerance(500.0, 500.29)); // 0.29 Da
        assert!(!config.within_tolerance(500.0, 500.31)); // 0.31 Da

        // Test mass_error calculation
        let error = config.mass_error(500.0, 500.25);
        assert!((error - 0.25).abs() < 1e-10);
    }

    /// Verifies that default FragmentToleranceConfig is 10 ppm for HRAM data.
    #[test]
    fn test_fragment_tolerance_default() {
        let config = FragmentToleranceConfig::default();
        // Default should be 10 ppm for HRAM
        assert_eq!(config.tolerance, 10.0);
        assert_eq!(config.unit, ToleranceUnit::Ppm);
    }
}
