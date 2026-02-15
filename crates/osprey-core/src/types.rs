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
    ///
    /// Returns `None` for empty string or "NOLOSS" (indicating no neutral loss)
    pub fn parse(s: &str) -> Option<Self> {
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

/// MS1 (survey) spectrum from data file
///
/// Used for precursor isotope envelope extraction and MS1 mass calibration.
/// pyXcorrDIA extracts M+0 peak from MS1 spectra for accurate mass calibration.
#[derive(Debug, Clone)]
pub struct MS1Spectrum {
    /// Scan number in the file
    pub scan_number: u32,
    /// Retention time in minutes
    pub retention_time: f64,
    /// m/z values
    pub mzs: Vec<f64>,
    /// Intensities
    pub intensities: Vec<f32>,
}

impl MS1Spectrum {
    /// Create a new MS1 spectrum
    pub fn new(scan_number: u32, retention_time: f64) -> Self {
        Self {
            scan_number,
            retention_time,
            mzs: Vec::new(),
            intensities: Vec::new(),
        }
    }

    /// Get the number of peaks
    pub fn len(&self) -> usize {
        self.mzs.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.mzs.is_empty()
    }

    /// Find peak within tolerance and return (observed_mz, intensity)
    ///
    /// Returns the most intense peak within the tolerance window.
    pub fn find_peak_ppm(&self, target_mz: f64, tolerance_ppm: f64) -> Option<(f64, f32)> {
        let mut best_intensity = 0.0f32;
        let mut best_mz = None;

        for (&mz, &intensity) in self.mzs.iter().zip(self.intensities.iter()) {
            let ppm_error = ((mz - target_mz) / target_mz).abs() * 1e6;
            if ppm_error <= tolerance_ppm && intensity > best_intensity {
                best_intensity = intensity;
                best_mz = Some(mz);
            }
        }

        best_mz.map(|mz| (mz, best_intensity))
    }
}

/// Isotope envelope extracted from MS1 spectrum (pyXcorrDIA-compatible)
///
/// Contains intensities for M-1, M+0, M+1, M+2, M+3 isotope peaks.
/// The M+0 (monoisotopic) peak is used for MS1 mass calibration.
#[derive(Debug, Clone)]
pub struct IsotopeEnvelope {
    /// Extracted intensities [M-1, M+0, M+1, M+2, M+3]
    pub intensities: [f64; 5],
    /// Observed M+0 m/z (for mass calibration)
    pub m0_observed_mz: Option<f64>,
    /// Mass error of M+0 peak (observed - theoretical) in Da
    pub m0_mass_error: Option<f64>,
}

impl IsotopeEnvelope {
    /// Neutron mass for isotope spacing calculation
    ///
    /// This value (1.002868 Da) comes from Monocle:
    /// Rad R, et al. "Improved Monoisotopic Mass Estimation for Deeper Proteome Coverage"
    /// https://pmc.ncbi.nlm.nih.gov/articles/PMC12204116/
    ///
    /// Also used by pyXcorrDIA for isotope envelope extraction.
    pub const NEUTRON_MASS: f64 = 1.002868;

    /// Calculate expected isotope m/z values for a precursor
    ///
    /// Returns [M-1, M+0, M+1, M+2, M+3] m/z values
    pub fn calculate_isotope_mzs(precursor_mz: f64, charge: u8) -> [f64; 5] {
        let isotope_gap = Self::NEUTRON_MASS / charge as f64;
        [
            precursor_mz - isotope_gap,       // M-1
            precursor_mz,                     // M+0 (monoisotopic)
            precursor_mz + isotope_gap,       // M+1
            precursor_mz + 2.0 * isotope_gap, // M+2
            precursor_mz + 3.0 * isotope_gap, // M+3
        ]
    }

    /// Extract isotope envelope from MS1 spectrum (pyXcorrDIA-compatible)
    ///
    /// For each expected isotope m/z, finds the best matching peak within tolerance.
    /// Returns extracted intensities and M+0 mass error for calibration.
    ///
    /// # Arguments
    /// * `ms1` - MS1 spectrum to extract from
    /// * `precursor_mz` - Theoretical monoisotopic precursor m/z
    /// * `charge` - Precursor charge state
    /// * `tolerance_ppm` - m/z matching tolerance in ppm (default: 10)
    pub fn extract(ms1: &MS1Spectrum, precursor_mz: f64, charge: u8, tolerance_ppm: f64) -> Self {
        let expected_mzs = Self::calculate_isotope_mzs(precursor_mz, charge);
        let mut intensities = [0.0f64; 5];
        let mut m0_observed_mz = None;
        let mut m0_mass_error = None;

        for (i, &target_mz) in expected_mzs.iter().enumerate() {
            if let Some((obs_mz, intensity)) = ms1.find_peak_ppm(target_mz, tolerance_ppm) {
                intensities[i] = intensity as f64;

                // Track M+0 peak (index 1) for mass calibration
                if i == 1 {
                    m0_observed_mz = Some(obs_mz);
                    m0_mass_error = Some(obs_mz - target_mz);
                }
            }
        }

        Self {
            intensities,
            m0_observed_mz,
            m0_mass_error,
        }
    }

    /// Calculate M+0 mass error in ppm
    pub fn m0_ppm_error(&self, precursor_mz: f64) -> Option<f64> {
        self.m0_mass_error.map(|delta| delta / precursor_mz * 1e6)
    }

    /// Get M+0 intensity
    pub fn m0_intensity(&self) -> f64 {
        self.intensities[1]
    }

    /// Check if M+0 was detected
    pub fn has_m0(&self) -> bool {
        self.m0_observed_mz.is_some()
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
    /// Unexplained intensity (residual = ||Ax - b||²)
    pub residual: f64,
    /// Total number of candidates considered in this spectrum
    pub n_candidates: u32,
    /// Sum of all coefficients (for relative_coefficient calculation)
    pub coefficient_sum: f64,
    /// Observed spectrum norm (||b||² for explained variance calculation)
    pub observed_norm: f64,
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
            n_candidates: 0,
            coefficient_sum: 0.0,
            observed_norm: 0.0,
        }
    }

    /// Compute explained variance: 1 - residual / ||b||²
    pub fn explained_variance(&self) -> f64 {
        if self.observed_norm > 1e-10 {
            1.0 - self.residual / self.observed_norm
        } else {
            0.0
        }
    }

    /// Get relative coefficient for a given library ID
    pub fn relative_coefficient(&self, lib_id: u32) -> f64 {
        if self.coefficient_sum > 1e-10 {
            self.library_ids
                .iter()
                .zip(self.coefficients.iter())
                .find(|(id, _)| **id == lib_id)
                .map(|(_, coef)| *coef / self.coefficient_sum)
                .unwrap_or(0.0)
        } else {
            0.0
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
    /// Signal-to-noise ratio on coefficient series (peak apex vs background)
    pub signal_to_noise: f64,
    /// Signal-to-noise ratio on best fragment XIC (highest correlation to coefficient series)
    pub xic_signal_to_noise: f64,

    // Spectral features (from mixed/observed spectrum at apex)
    /// X!Tandem-style hyperscore: log(n_b!) + log(n_y!) + Σlog(I_f+1)
    pub hyperscore: f64,
    /// XCorr score (Comet-style cross-correlation)
    pub xcorr: f64,
    /// Dot product (LibCosine with sqrt preprocessing)
    pub dot_product: f64,
    /// LibCosine with sqrt(intensity)*mz² (SMZ) preprocessing
    pub dot_product_smz: f64,
    /// Fraction of predicted fragments detected
    pub fragment_coverage: f64,
    /// Backbone coverage
    pub sequence_coverage: f64,
    /// Longest consecutive b/y ion run
    pub consecutive_ions: u32,
    /// Rank of base peak in predicted
    pub base_peak_rank: u32,
    /// Number of top-6 predicted fragments matched
    pub top6_matches: u32,
    /// Fraction of observed intensity explained (at apex spectrum)
    pub explained_intensity: f64,

    // Deconvoluted spectral features (from summed deconvoluted signal, apex ± 2 scans)
    // These compare the deconvoluted (coefficient-scaled) contribution to the library,
    // removing interference from co-eluting peptides.
    /// Hyperscore from deconvoluted spectrum (apex ± 2 scans)
    pub hyperscore_deconv: f64,
    /// XCorr from deconvoluted spectrum
    pub xcorr_deconv: f64,
    /// Dot product from deconvoluted spectrum
    pub dot_product_deconv: f64,
    /// Dot product SMZ from deconvoluted spectrum
    pub dot_product_smz_deconv: f64,
    /// Fragment coverage from deconvoluted spectrum
    pub fragment_coverage_deconv: f64,
    /// Sequence coverage from deconvoluted spectrum
    pub sequence_coverage_deconv: f64,
    /// Consecutive ions from deconvoluted spectrum
    pub consecutive_ions_deconv: u32,
    /// Top-6 matches from deconvoluted spectrum
    pub top6_matches_deconv: u32,

    // Contextual features (computed at apex regression result)
    /// Number of competing candidates at apex spectrum
    pub n_competitors: u32,
    /// Coefficient relative to sum of all at apex spectrum
    pub relative_coefficient: f64,
    /// Mean relative coefficient over apex ± 2 spectra
    pub relative_coefficient_mean: f64,
    /// Local peptide density
    pub local_peptide_density: f64,
    /// Spectral complexity estimate
    pub spectral_complexity: f64,
    /// Regression residual
    pub regression_residual: f64,
    /// Number of modifications
    pub modification_count: u32,

    // Fragment co-elution features (DIA-NN pTimeCorr-inspired)
    // Computed within peak integration boundaries only (start_rt..end_rt)
    /// Sum of Pearson correlations between each fragment XIC and the coefficient time series
    pub fragment_coelution_sum: f64,
    /// Minimum per-fragment correlation with coefficient time series
    pub fragment_coelution_min: f64,
    /// Number of fragments with positive correlation to coefficient series
    pub n_coeluting_fragments: u32,

    // Per-fragment correlation array (top 6 fragments by library intensity)
    // Correlations computed within peak integration boundaries only.
    // Fragments not present in library are 0-padded.
    /// Correlation of fragment 0 (highest library intensity) XIC vs coefficient series
    pub fragment_corr_0: f64,
    /// Correlation of fragment 1 XIC vs coefficient series
    pub fragment_corr_1: f64,
    /// Correlation of fragment 2 XIC vs coefficient series
    pub fragment_corr_2: f64,
    /// Correlation of fragment 3 XIC vs coefficient series
    pub fragment_corr_3: f64,
    /// Correlation of fragment 4 XIC vs coefficient series
    pub fragment_corr_4: f64,
    /// Correlation of fragment 5 XIC vs coefficient series
    pub fragment_corr_5: f64,

    // Elution-weighted spectral similarity
    /// Cosine similarity computed at each scan within peak boundaries,
    /// weighted by coefficient², then averaged. Captures whether the library
    /// pattern holds across the entire peak, not just at the apex.
    pub elution_weighted_cosine: f64,

    // Per-fragment mass accuracy at apex (ppm for HRAM, Th for unit resolution)
    /// Signed mean mass error across matched fragments at apex (systematic bias)
    pub mass_accuracy_deviation_mean: f64,
    /// Absolute mean mass error across matched fragments at apex
    pub abs_mass_accuracy_deviation_mean: f64,
    /// Standard deviation of mass errors across matched fragments at apex
    pub mass_accuracy_std: f64,

    // Percolator-style features
    /// Absolute RT deviation (|apex_rt - expected_rt|, minutes)
    pub abs_rt_deviation: f64,
    /// Peptide sequence length (unmodified)
    pub peptide_length: u32,
    /// Number of missed tryptic cleavages (internal K/R)
    pub missed_cleavages: u32,
    /// ln(number of candidates in regression)
    pub ln_num_candidates: f64,
    /// Z-score of this peptide's coefficient at apex spectrum
    /// (how many SDs above the mean of all coefficients in the spectrum)
    pub coef_zscore: f64,
    /// Mean z-score across apex ± 2 spectra (smoothed version)
    pub coef_zscore_mean: f64,

    // MS1-based features (HRAM only, 0.0 for unit resolution or missing MS1)
    /// Pearson correlation between coefficient time series and MS1 monoisotopic
    /// precursor XIC within peak integration boundaries.
    pub ms1_precursor_coelution: f64,
    /// Cosine similarity between observed MS1 isotope envelope at apex RT and
    /// theoretical isotope distribution from peptide sequence.
    pub ms1_isotope_cosine: f64,

    // Tukey median polish features (fragment XIC decomposition)
    /// Cosine similarity between median-polish row effects (data-derived fragment intensities)
    /// and library fragment intensities in linear space with sqrt preprocessing.
    /// Higher = better agreement between observed and expected fragment pattern.
    pub median_polish_cosine: f64,
    /// R² of the additive model (μ + α_f + β_s) in linear space.
    /// Measures how well the shared elution profile + fragment intensities explain
    /// the observed XIC matrix. Zero-intensity cells included (penalizes missing signal).
    pub median_polish_rsquared: f64,
    /// Fraction of total signal unexplained: Σ|obs-pred| / Σobs in linear space.
    /// Lower = cleaner co-elution. Zero-intensity cells included.
    pub median_polish_residual_ratio: f64,
}

// =============================================================================
// Coelution search types (DIA-NN-style fragment coelution without regression)
// =============================================================================

/// Peak boundaries detected from fragment XICs using DIA-NN-style valley detection.
///
/// Unlike coefficient-based `PeakBoundaries`, these are derived from fragment
/// chromatograms with adaptive boundary detection that handles overlapping peaks.
#[derive(Debug, Clone)]
pub struct XICPeakBounds {
    /// Retention time at peak apex (minutes)
    pub apex_rt: f64,
    /// Intensity at apex (from smoothed reference XIC)
    pub apex_intensity: f64,
    /// Index into the spectrum/XIC array at peak apex
    pub apex_index: usize,
    /// Start retention time of peak (minutes)
    pub start_rt: f64,
    /// End retention time of peak (minutes)
    pub end_rt: f64,
    /// Index into the spectrum/XIC array at peak start
    pub start_index: usize,
    /// Index into the spectrum/XIC array at peak end
    pub end_index: usize,
    /// Integrated area (trapezoidal) within boundaries
    pub area: f64,
    /// Signal-to-noise: (apex - bg_mean) / bg_sd
    pub signal_to_noise: f64,
}

/// Feature set for coelution-based scoring (no ridge regression).
///
/// ~45 features computed from fragment XICs, pairwise correlation,
/// spectral matching at apex, and peptide properties.
#[derive(Debug, Clone, Default)]
pub struct CoelutionFeatureSet {
    // --- Pairwise coelution features (from fragment XIC correlation matrix) ---
    /// Sum of all pairwise Pearson correlations between fragment XICs
    pub coelution_sum: f64,
    /// Minimum pairwise correlation (worst pair)
    pub coelution_min: f64,
    /// Maximum pairwise correlation (best pair)
    pub coelution_max: f64,
    /// Number of fragments with positive mean correlation
    pub n_coeluting_fragments: u8,
    /// Number of fragment pairs in correlation matrix
    pub n_fragment_pairs: u8,
    /// Per-fragment average correlation with all other fragments (top 6, 0-padded)
    pub fragment_corr: [f64; 6],

    // --- Peak shape features (from reference XIC) ---
    /// Peak apex intensity (background-subtracted)
    pub peak_apex: f64,
    /// Integrated peak area within boundaries
    pub peak_area: f64,
    /// Peak width (FWHM in minutes)
    pub peak_width: f64,
    /// Peak symmetry (leading/trailing area ratio around apex)
    pub peak_symmetry: f64,
    /// Signal-to-noise ratio
    pub signal_to_noise: f64,
    /// Number of scans within peak boundaries
    pub n_scans: u16,
    /// Peak sharpness (steepness of edges)
    pub peak_sharpness: f64,

    // --- Spectral features at apex (from SpectralScorer) ---
    /// X!Tandem-style hyperscore
    pub hyperscore: f64,
    /// Comet-style cross-correlation score
    pub xcorr: f64,
    /// Library cosine (dot product with sqrt intensity)
    pub dot_product: f64,
    /// Library cosine with sqrt(intensity) * mz² preprocessing
    pub dot_product_smz: f64,
    /// Library cosine using only top 6 library fragments by intensity
    pub dot_product_top6: f64,
    /// Library cosine using only top 5 library fragments by intensity
    pub dot_product_top5: f64,
    /// Library cosine using only top 4 library fragments by intensity
    pub dot_product_top4: f64,
    /// SMZ cosine using only top 6 library fragments by intensity
    pub dot_product_smz_top6: f64,
    /// SMZ cosine using only top 5 library fragments by intensity
    pub dot_product_smz_top5: f64,
    /// SMZ cosine using only top 4 library fragments by intensity
    pub dot_product_smz_top4: f64,
    /// Fraction of library fragments matched in observed spectrum
    pub fragment_coverage: f64,
    /// Fraction of peptide backbone covered by b/y ions
    pub sequence_coverage: f64,
    /// Longest consecutive b or y ion series
    pub consecutive_ions: u8,
    /// Rank of observed base peak among library fragments
    pub base_peak_rank: u8,
    /// Count of top 6 library fragments matched
    pub top6_matches: u8,
    /// Fraction of total observed intensity explained by matches
    pub explained_intensity: f64,
    /// LibCosine at each scan within peak, weighted by reference XIC intensity²
    pub elution_weighted_cosine: f64,

    // --- Mass accuracy features ---
    /// Signed mean mass error across matched fragments (ppm or Th)
    pub mass_accuracy_mean: f64,
    /// Absolute mean mass error
    pub abs_mass_accuracy_mean: f64,
    /// Standard deviation of mass errors
    pub mass_accuracy_std: f64,

    // --- RT deviation ---
    /// RT deviation from expected (calibrated) RT (minutes)
    pub rt_deviation: f64,
    /// Absolute RT deviation
    pub abs_rt_deviation: f64,

    // --- MS1 features (HRAM only) ---
    /// Pearson correlation between reference XIC and MS1 precursor XIC
    pub ms1_precursor_coelution: f64,
    /// Cosine similarity of observed vs theoretical isotope envelope
    pub ms1_isotope_cosine: f64,

    // --- Peptide properties ---
    /// Number of modifications on the peptide
    pub modification_count: u8,
    /// Unmodified sequence length
    pub peptide_length: u8,
    /// Number of missed cleavage sites (internal K/R)
    pub missed_cleavages: u8,

    // --- Median polish features (from XIC matrix decomposition) ---
    /// Cosine similarity between median polish row effects and library intensities
    pub median_polish_cosine: f64,
    /// R² of the additive model in linear space
    pub median_polish_rsquared: f64,
    /// Fraction of signal unexplained by additive model
    pub median_polish_residual_ratio: f64,
}

/// Scored entry from coelution-based search.
///
/// Holds all information needed for Mokapot FDR (via PIN file) and
/// blib output (spectral library). Produced by `run_coelution_search()`.
#[derive(Debug, Clone)]
pub struct CoelutionScoredEntry {
    /// Library entry ID
    pub entry_id: u32,
    /// Whether this is a decoy
    pub is_decoy: bool,
    /// Unmodified sequence
    pub sequence: String,
    /// Modified sequence
    pub modified_sequence: String,
    /// Precursor charge
    pub charge: u8,
    /// Precursor m/z
    pub precursor_mz: f64,
    /// Protein accession(s)
    pub protein_ids: Vec<String>,
    /// Scan number at apex
    pub scan_number: u32,
    /// Apex RT (minutes)
    pub apex_rt: f64,
    /// Peak boundaries from DIA-NN-style detection
    pub peak_bounds: XICPeakBounds,
    /// Coelution feature set (~45 features for Mokapot)
    pub features: CoelutionFeatureSet,
    /// Fragment m/z values at apex (for blib output)
    pub fragment_mzs: Vec<f64>,
    /// Fragment intensities at apex (for blib output)
    pub fragment_intensities: Vec<f32>,
    /// Reference XIC values within peak bounds (for blib coefficient series)
    pub reference_xic: Vec<(f64, f64)>,
    /// Source file name (file stem, e.g. "sample1")
    pub file_name: String,
    /// Run-level q-value (from Mokapot on single file, default 1.0)
    pub run_qvalue: f64,
    /// Experiment-level q-value (from Mokapot across files, default 1.0)
    pub experiment_qvalue: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies that IsolationWindow correctly includes/excludes m/z values at boundaries.
    #[test]
    fn test_isolation_window_contains() {
        let window = IsolationWindow::symmetric(500.0, 12.5);
        assert!(window.contains(500.0));
        assert!(window.contains(487.5));
        assert!(window.contains(512.5));
        assert!(!window.contains(487.4));
        assert!(!window.contains(512.6));
    }

    /// Verifies that NeutralLoss::H2O and NH3 return correct monoisotopic masses.
    #[test]
    fn test_neutral_loss_mass() {
        assert!((NeutralLoss::H2O.mass() - 18.010565).abs() < 1e-6);
        assert!((NeutralLoss::NH3.mass() - 17.026549).abs() < 1e-6);
    }

    /// Verifies that IonType::from_char parses b/y/z ions and returns Unknown for invalid chars.
    #[test]
    fn test_ion_type_from_char() {
        assert_eq!(IonType::from_char('b'), IonType::B);
        assert_eq!(IonType::from_char('Y'), IonType::Y);
        assert_eq!(IonType::from_char('z'), IonType::Z);
        assert_eq!(IonType::from_char('?'), IonType::Unknown);
    }

    /// Verifies that MS1 find_peak_ppm returns the most intense peak within ppm tolerance.
    #[test]
    fn test_ms1_spectrum_find_peak() {
        let mut ms1 = MS1Spectrum::new(1, 10.0);
        ms1.mzs = vec![499.995, 500.0, 500.005, 501.0];
        ms1.intensities = vec![50.0, 100.0, 75.0, 25.0];

        // Find peak at 500.0 with 10 ppm tolerance
        let result = ms1.find_peak_ppm(500.0, 10.0);
        assert!(result.is_some());
        let (mz, int) = result.unwrap();
        // Should find the 100.0 intensity peak at 500.0
        assert!((mz - 500.0).abs() < 0.001);
        assert!((int - 100.0).abs() < 0.1);

        // Find peak at 501.0
        let result2 = ms1.find_peak_ppm(501.0, 10.0);
        assert!(result2.is_some());

        // No peak at 600.0
        let result3 = ms1.find_peak_ppm(600.0, 10.0);
        assert!(result3.is_none());
    }

    /// Verifies isotope m/z spacing calculation for charge 2 using neutron mass.
    #[test]
    fn test_isotope_envelope_mz_calculation() {
        // Test isotope m/z calculation for charge 2
        let mzs = IsotopeEnvelope::calculate_isotope_mzs(500.0, 2);
        let gap = IsotopeEnvelope::NEUTRON_MASS / 2.0; // ~0.501434

        assert!((mzs[0] - (500.0 - gap)).abs() < 1e-6); // M-1
        assert!((mzs[1] - 500.0).abs() < 1e-6); // M+0
        assert!((mzs[2] - (500.0 + gap)).abs() < 1e-6); // M+1
        assert!((mzs[3] - (500.0 + 2.0 * gap)).abs() < 1e-6); // M+2
        assert!((mzs[4] - (500.0 + 3.0 * gap)).abs() < 1e-6); // M+3
    }

    /// Verifies isotope envelope extraction recovers intensities and M+0 mass error from MS1.
    #[test]
    fn test_isotope_envelope_extraction() {
        // Create MS1 spectrum with isotope peaks
        let mut ms1 = MS1Spectrum::new(1, 10.0);
        let gap = IsotopeEnvelope::NEUTRON_MASS / 2.0; // charge 2

        // Add isotope peaks with slight mass error on M+0
        let m0_observed = 500.002; // 4 ppm error
        ms1.mzs = vec![
            500.0 - gap,       // M-1
            m0_observed,       // M+0 with +2 mDa error
            500.0 + gap,       // M+1
            500.0 + 2.0 * gap, // M+2
            500.0 + 3.0 * gap, // M+3
        ];
        ms1.intensities = vec![10.0, 100.0, 80.0, 40.0, 15.0];

        // Extract envelope
        let envelope = IsotopeEnvelope::extract(&ms1, 500.0, 2, 10.0);

        // Check intensities
        assert!((envelope.intensities[0] - 10.0).abs() < 0.1); // M-1
        assert!((envelope.intensities[1] - 100.0).abs() < 0.1); // M+0
        assert!((envelope.intensities[2] - 80.0).abs() < 0.1); // M+1
        assert!((envelope.intensities[3] - 40.0).abs() < 0.1); // M+2
        assert!((envelope.intensities[4] - 15.0).abs() < 0.1); // M+3

        // Check M+0 was detected
        assert!(envelope.has_m0());
        assert!(envelope.m0_observed_mz.is_some());
        assert!((envelope.m0_observed_mz.unwrap() - m0_observed).abs() < 1e-6);

        // Check mass error (observed - theoretical = 0.002 Da)
        assert!(envelope.m0_mass_error.is_some());
        assert!((envelope.m0_mass_error.unwrap() - 0.002).abs() < 1e-6);

        // Check ppm error (0.002 / 500.0 * 1e6 = 4 ppm)
        let ppm = envelope.m0_ppm_error(500.0).unwrap();
        assert!((ppm - 4.0).abs() < 0.1);
    }

    /// Verifies that missing isotope peaks are reported as zero while present ones are detected.
    #[test]
    fn test_isotope_envelope_missing_peaks() {
        // MS1 spectrum missing some isotope peaks
        let mut ms1 = MS1Spectrum::new(1, 10.0);
        ms1.mzs = vec![500.0, 500.5]; // Only M+0 and M+1
        ms1.intensities = vec![100.0, 80.0];

        let envelope = IsotopeEnvelope::extract(&ms1, 500.0, 2, 10.0);

        // M-1 and M+2, M+3 should be 0
        assert!(envelope.intensities[0] < 0.1); // M-1 missing
        assert!(envelope.intensities[1] > 90.0); // M+0 present
        assert!(envelope.intensities[2] > 70.0); // M+1 present
        assert!(envelope.intensities[3] < 0.1); // M+2 missing
        assert!(envelope.intensities[4] < 0.1); // M+3 missing

        // M+0 should still be detected
        assert!(envelope.has_m0());
    }
}
