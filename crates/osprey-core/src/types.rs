//! Core data types for Osprey DIA analysis
//!
//! These types represent the fundamental data structures used throughout
//! the analysis pipeline, from spectral library entries to detection results.

use serde::{Deserialize, Serialize};

/// Library entry representing a peptide precursor with spectral information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibraryEntry {
    /// Unique identifier
    pub id: u32,
    /// Unmodified amino acid sequence
    pub sequence: String,
    /// Sequence with modification notation
    pub modified_sequence: String,
    /// List of modifications
    pub modifications: Vec<Modification>,
    /// Precursor charge state
    pub charge: u8,
    /// Precursor m/z
    pub precursor_mz: f64,
    /// Normalized or calibrated retention time
    pub retention_time: f64,
    /// Whether RT is calibrated to this run
    pub rt_calibrated: bool,
    /// Fragment ions with intensities
    pub fragments: Vec<LibraryFragment>,
    /// Protein accession(s)
    pub protein_ids: Vec<String>,
    /// Gene name(s)
    pub gene_names: Vec<String>,
    /// Whether this is a decoy entry
    pub is_decoy: bool,
}

impl LibraryEntry {
    /// Create a new library entry with minimal required fields
    pub fn new(
        id: u32,
        sequence: String,
        modified_sequence: String,
        charge: u8,
        precursor_mz: f64,
        retention_time: f64,
    ) -> Self {
        Self {
            id,
            sequence,
            modified_sequence,
            modifications: Vec::new(),
            charge,
            precursor_mz,
            retention_time,
            rt_calibrated: false,
            fragments: Vec::new(),
            protein_ids: Vec::new(),
            gene_names: Vec::new(),
            is_decoy: false,
        }
    }
}

/// Fragment ion from library
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibraryFragment {
    /// Fragment m/z
    pub mz: f64,
    /// Normalized relative intensity (0-1 or 0-100)
    pub relative_intensity: f32,
    /// Fragment annotation
    pub annotation: FragmentAnnotation,
}

/// Modification on a peptide
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Modification {
    /// 0-indexed position in sequence
    pub position: usize,
    /// UniMod accession if known
    pub unimod_id: Option<u32>,
    /// Monoisotopic mass change
    pub mass_delta: f64,
    /// Human-readable name (e.g., "Oxidation")
    pub name: Option<String>,
}

/// Fragment ion annotation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FragmentAnnotation {
    /// Ion type (b, y, etc.)
    pub ion_type: IonType,
    /// Ion number (ordinal)
    pub ordinal: u8,
    /// Fragment charge
    pub charge: u8,
    /// Neutral loss if any
    pub neutral_loss: Option<NeutralLoss>,
}

impl Default for FragmentAnnotation {
    fn default() -> Self {
        Self {
            ion_type: IonType::Unknown,
            ordinal: 0,
            charge: 1,
            neutral_loss: None,
        }
    }
}

/// Ion types for fragment annotation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IonType {
    B,
    Y,
    A,
    C,
    X,
    Z,
    Precursor,
    Internal,
    Immonium,
    Unknown,
}

impl IonType {
    /// Parse ion type from character
    pub fn from_char(c: char) -> Self {
        match c.to_ascii_lowercase() {
            'b' => IonType::B,
            'y' => IonType::Y,
            'a' => IonType::A,
            'c' => IonType::C,
            'x' => IonType::X,
            'z' => IonType::Z,
            'p' => IonType::Precursor,
            'm' => IonType::Precursor, // M for precursor in some formats
            _ => IonType::Unknown,
        }
    }
}

/// Neutral losses
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum NeutralLoss {
    /// Water loss (-18.0106)
    H2O,
    /// Ammonia loss (-17.0265)
    NH3,
    /// Phosphoric acid loss (-97.9769)
    H3PO4,
    /// Custom neutral loss with mass
    Custom(f64),
}

impl NeutralLoss {
    /// Get the mass of the neutral loss
    pub fn mass(&self) -> f64 {
        match self {
            NeutralLoss::H2O => 18.010565,
            NeutralLoss::NH3 => 17.026549,
            NeutralLoss::H3PO4 => 97.976896,
            NeutralLoss::Custom(m) => *m,
        }
    }

    /// Parse neutral loss from string
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_uppercase().as_str() {
            "" | "NOLOSS" => None,
            "H2O" | "WATER" => Some(NeutralLoss::H2O),
            "NH3" | "AMMONIA" => Some(NeutralLoss::NH3),
            "H3PO4" | "PHOSPHO" => Some(NeutralLoss::H3PO4),
            _ => {
                // Try to parse as a number
                s.parse::<f64>().ok().map(NeutralLoss::Custom)
            }
        }
    }
}

/// MS/MS spectrum from data file
#[derive(Debug, Clone)]
pub struct Spectrum {
    /// Scan number in the file
    pub scan_number: u32,
    /// Retention time in minutes
    pub retention_time: f64,
    /// Precursor m/z (center of isolation window)
    pub precursor_mz: f64,
    /// DIA isolation window
    pub isolation_window: IsolationWindow,
    /// Fragment m/z values
    pub mzs: Vec<f64>,
    /// Fragment intensities
    pub intensities: Vec<f32>,
}

impl Spectrum {
    /// Create a new empty spectrum
    pub fn new(scan_number: u32, retention_time: f64, isolation_window: IsolationWindow) -> Self {
        Self {
            scan_number,
            retention_time,
            precursor_mz: isolation_window.center,
            isolation_window,
            mzs: Vec::new(),
            intensities: Vec::new(),
        }
    }

    /// Check if a precursor m/z falls within this spectrum's isolation window
    pub fn contains_precursor(&self, mz: f64) -> bool {
        self.isolation_window.contains(mz)
    }

    /// Get the number of peaks in the spectrum
    pub fn len(&self) -> usize {
        self.mzs.len()
    }

    /// Check if spectrum is empty
    pub fn is_empty(&self) -> bool {
        self.mzs.is_empty()
    }
}

/// DIA isolation window
#[derive(Debug, Clone, Copy)]
pub struct IsolationWindow {
    /// Center m/z
    pub center: f64,
    /// Lower offset from center (positive value)
    pub lower_offset: f64,
    /// Upper offset from center (positive value)
    pub upper_offset: f64,
}

impl IsolationWindow {
    /// Create a symmetric isolation window
    pub fn symmetric(center: f64, half_width: f64) -> Self {
        Self {
            center,
            lower_offset: half_width,
            upper_offset: half_width,
        }
    }

    /// Create an asymmetric isolation window
    pub fn new(center: f64, lower_offset: f64, upper_offset: f64) -> Self {
        Self {
            center,
            lower_offset,
            upper_offset,
        }
    }

    /// Get the lower bound of the window
    pub fn lower_bound(&self) -> f64 {
        self.center - self.lower_offset
    }

    /// Get the upper bound of the window
    pub fn upper_bound(&self) -> f64 {
        self.center + self.upper_offset
    }

    /// Get the total width of the window
    pub fn width(&self) -> f64 {
        self.lower_offset + self.upper_offset
    }

    /// Check if an m/z value falls within this window
    pub fn contains(&self, mz: f64) -> bool {
        mz >= self.lower_bound() && mz <= self.upper_bound()
    }
}

/// Binned spectrum for regression
#[derive(Debug, Clone)]
pub struct BinnedSpectrum {
    /// Which bins have signal
    pub bin_indices: Vec<u32>,
    /// Intensity in each bin
    pub intensities: Vec<f32>,
}

impl BinnedSpectrum {
    /// Create a new binned spectrum
    pub fn new(bin_indices: Vec<u32>, intensities: Vec<f32>) -> Self {
        debug_assert_eq!(bin_indices.len(), intensities.len());
        Self {
            bin_indices,
            intensities,
        }
    }

    /// Create an empty binned spectrum
    pub fn empty() -> Self {
        Self {
            bin_indices: Vec::new(),
            intensities: Vec::new(),
        }
    }

    /// Get the number of non-zero bins
    pub fn len(&self) -> usize {
        self.bin_indices.len()
    }

    /// Check if the binned spectrum is empty
    pub fn is_empty(&self) -> bool {
        self.bin_indices.is_empty()
    }
}

/// Regression result for one spectrum
#[derive(Debug, Clone)]
pub struct RegressionResult {
    /// Scan number
    pub scan_number: u32,
    /// Retention time
    pub retention_time: f64,
    /// Indices into library for candidates with non-zero coefficients
    pub library_ids: Vec<u32>,
    /// Corresponding coefficients
    pub coefficients: Vec<f64>,
    /// Unexplained intensity (residual)
    pub residual: f64,
}

impl RegressionResult {
    /// Create a new regression result
    pub fn new(scan_number: u32, retention_time: f64) -> Self {
        Self {
            scan_number,
            retention_time,
            library_ids: Vec::new(),
            coefficients: Vec::new(),
            residual: 0.0,
        }
    }
}

/// Peptide detection result across the experiment
#[derive(Debug, Clone)]
pub struct PeptideDetection {
    /// Reference to library entry ID
    pub library_entry_id: u32,
    /// Results per run
    pub run_results: Vec<RunDetection>,
    /// Experiment-level q-value
    pub experiment_q_value: f64,
    /// Number of runs with detection
    pub n_runs_detected: u32,
    /// Total runs searched
    pub n_runs_searched: u32,
    /// Was this detected in Step 1 of two-step search?
    pub detected_step1: bool,
}

/// Per-run detection result
#[derive(Debug, Clone)]
pub struct RunDetection {
    /// Run file name
    pub file_name: String,
    /// Whether peptide was detected in this run
    pub detected: bool,
    /// Peak boundaries if detected
    pub peak_boundaries: Option<PeakBoundaries>,
    /// Run-level q-value
    pub run_q_value: f64,
    /// Discriminant score from ML model
    pub discriminant_score: f64,
    /// Full feature set
    pub features: FeatureSet,
}

/// Peak boundary information for Skyline
#[derive(Debug, Clone, Copy)]
pub struct PeakBoundaries {
    /// Peak start time (minutes)
    pub start_rt: f64,
    /// Peak end time (minutes)
    pub end_rt: f64,
    /// Peak apex time (minutes)
    pub apex_rt: f64,
    /// Coefficient at apex
    pub apex_coefficient: f64,
    /// Sum of coefficients (for QC, not quant)
    pub integrated_area: f64,
    /// Peak quality metrics
    pub peak_quality: PeakQuality,
}

/// Peak quality flags
#[derive(Debug, Clone, Copy, Default)]
pub struct PeakQuality {
    /// R² of EMG fit
    pub emg_fit_r2: f64,
    /// Peak appears split
    pub is_split: bool,
    /// Peak truncated at gradient edge
    pub is_truncated: bool,
    /// Shoulder detected
    pub has_shoulder: bool,
    /// Width relative to other peaks (0-100 percentile)
    pub width_percentile: f64,
}

/// Complete feature set for scoring
#[derive(Debug, Clone, Default)]
pub struct FeatureSet {
    // Chromatographic features
    /// Peak apex coefficient (maximum value)
    pub peak_apex: f64,
    /// Integrated peak area (AUC of coefficients)
    pub peak_area: f64,
    /// EMG fit quality (R²)
    pub emg_fit_quality: f64,
    /// Peak width (FWHM in minutes)
    pub peak_width: f64,
    /// Peak symmetry (leading/trailing ratio)
    pub peak_symmetry: f64,
    /// RT deviation from prediction (minutes)
    pub rt_deviation: f64,
    /// Normalized RT deviation
    pub rt_deviation_normalized: f64,
    /// Number of scans contributing to peak
    pub n_contributing_scans: u32,
    /// Coefficient variance near apex
    pub coefficient_stability: f64,
    /// Peak boundary sharpness
    pub peak_sharpness: f64,
    /// Peak prominence (apex / baseline)
    pub peak_prominence: f64,

    // Spectral features
    /// X!Tandem-style hyperscore
    pub hyperscore: f64,
    /// Normalized spectral contrast angle
    pub spectral_contrast_angle: f64,
    /// Dot product
    pub dot_product: f64,
    /// Pearson intensity correlation
    pub pearson_correlation: f64,
    /// Spearman rank correlation
    pub spearman_correlation: f64,
    /// Fraction of predicted fragments detected
    pub fragment_coverage: f64,
    /// Backbone coverage
    pub sequence_coverage: f64,
    /// Longest consecutive b/y ion run
    pub consecutive_ions: u32,
    /// Rank of base peak in predicted
    pub base_peak_rank: u32,
    /// Number of top-3 predicted fragments matched
    pub top3_matches: u32,
    /// Fraction of observed intensity explained
    pub explained_intensity: f64,

    // Contextual features
    /// Number of competing candidates
    pub n_competitors: u32,
    /// Coefficient relative to sum of all
    pub relative_coefficient: f64,
    /// Local peptide density
    pub local_peptide_density: f64,
    /// Spectral complexity estimate
    pub spectral_complexity: f64,
    /// Regression residual
    pub regression_residual: f64,
    /// Precursor intensity if MS1 available
    pub precursor_intensity: Option<f64>,
    /// Number of modifications
    pub modification_count: u32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isolation_window_contains() {
        let window = IsolationWindow::symmetric(500.0, 12.5);
        assert!(window.contains(500.0));
        assert!(window.contains(487.5));
        assert!(window.contains(512.5));
        assert!(!window.contains(487.4));
        assert!(!window.contains(512.6));
    }

    #[test]
    fn test_neutral_loss_mass() {
        assert!((NeutralLoss::H2O.mass() - 18.010565).abs() < 1e-6);
        assert!((NeutralLoss::NH3.mass() - 17.026549).abs() < 1e-6);
    }

    #[test]
    fn test_ion_type_from_char() {
        assert_eq!(IonType::from_char('b'), IonType::B);
        assert_eq!(IonType::from_char('Y'), IonType::Y);
        assert_eq!(IonType::from_char('z'), IonType::Z);
        assert_eq!(IonType::from_char('?'), IonType::Unknown);
    }
}
