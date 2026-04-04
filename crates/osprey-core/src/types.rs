//! Core data types for Osprey DIA analysis
//!
//! These types represent the fundamental data structures used throughout
//! the analysis pipeline, from spectral library entries to detection results.

use crate::config::FdrLevel;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

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

    /// Check if an m/z value falls within this window (half-open: [lower, upper))
    ///
    /// Half-open convention ensures entries at shared boundaries between adjacent
    /// windows belong to exactly one window, preventing double-counting.
    pub fn contains(&self, mz: f64) -> bool {
        mz >= self.lower_bound() && mz < self.upper_bound()
    }
}

// =============================================================================
// Peak boundary types
// =============================================================================

/// Peak boundaries for blib output and chromatographic peak detection.
///
/// Used by the blib writer to store per-run peak boundaries in the
/// OspreyPeakBoundaries table, and by the chromatographic peak detector.
#[derive(Debug, Clone, Default)]
pub struct PeakBoundaries {
    /// Start retention time (minutes)
    pub start_rt: f64,
    /// End retention time (minutes)
    pub end_rt: f64,
    /// Apex retention time (minutes)
    pub apex_rt: f64,
    /// Coefficient or intensity at apex
    pub apex_coefficient: f64,
    /// Integrated area under peak
    pub integrated_area: f64,
    /// Peak quality metrics
    pub peak_quality: PeakQuality,
}

/// Peak quality metrics
#[derive(Debug, Clone, Default)]
pub struct PeakQuality {
    /// Signal-to-noise ratio
    pub signal_to_noise: f64,
    /// Peak symmetry factor
    pub symmetry: f64,
    /// FWHM in minutes
    pub fwhm: f64,
}

// =============================================================================
// XIC peak types (DIA-NN-style fragment coelution)
// =============================================================================

/// Peak boundaries detected from fragment XICs using DIA-NN-style valley detection.
///
/// Derived from fragment chromatograms with adaptive boundary detection
/// that handles overlapping peaks.
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

/// Feature set for coelution-based scoring.
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

    // --- SG-weighted multi-scan spectral scores ---
    /// XCorr on raw DIA spectra at apex ±2 scans, Savitzky-Golay weighted
    pub sg_weighted_xcorr: f64,
    /// Cosine similarity on raw DIA spectra at apex ±2 scans, SG weighted
    pub sg_weighted_cosine: f64,
    /// Minimum per-fragment R² with the shared elution profile (weakest link)
    pub median_polish_min_fragment_r2: f64,
    /// Mean pairwise correlation of median polish residuals across fragments
    pub median_polish_residual_correlation: f64,
}

/// A CWT candidate peak with boundaries and coelution score.
///
/// Stored in parquet as packed LE bytes for top-N peak storage,
/// enabling inter-replicate peak reconciliation without re-running CWT.
#[derive(Debug, Clone, Copy, Default)]
pub struct CwtCandidate {
    /// Peak apex RT (minutes)
    pub apex_rt: f64,
    /// Start RT of peak (minutes)
    pub start_rt: f64,
    /// End RT of peak (minutes)
    pub end_rt: f64,
    /// Integrated area within boundaries
    pub area: f64,
    /// Signal-to-noise ratio
    pub snr: f64,
    /// Mean pairwise fragment correlation within peak
    pub coelution_score: f64,
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
    /// Run-level precursor q-value (modified_sequence + charge, default 1.0)
    pub run_precursor_qvalue: f64,
    /// Run-level peptide q-value (modified_sequence only, default 1.0)
    pub run_peptide_qvalue: f64,
    /// Run-level protein q-value (picked-protein FDR, default 1.0)
    pub run_protein_qvalue: f64,
    /// Experiment-level precursor q-value (modified_sequence + charge, default 1.0)
    pub experiment_precursor_qvalue: f64,
    /// Experiment-level peptide q-value (modified_sequence only, default 1.0)
    pub experiment_peptide_qvalue: f64,
    /// Experiment-level protein q-value (picked-protein FDR, default 1.0)
    pub experiment_protein_qvalue: f64,
    /// SVM discriminant score (from Percolator/Mokapot, same score used for q-value and PEP)
    pub score: f64,
    /// Posterior error probability (computed on the final SVM score)
    pub pep: f64,
    /// Top-N CWT candidate peaks (sorted by coelution_score descending).
    /// Used for inter-replicate peak reconciliation.
    pub cwt_candidates: Vec<CwtCandidate>,
}

impl CoelutionScoredEntry {
    /// Convert to a lightweight FDR stub, dropping heavy fields.
    ///
    /// The full entry data is already persisted to a `.scores.parquet` cache file
    /// and can be reloaded on demand for Percolator feature extraction,
    /// reconciliation CWT lookup, or blib/report output.
    pub fn to_fdr_entry(&self) -> FdrEntry {
        FdrEntry {
            entry_id: self.entry_id,
            parquet_index: 0, // populated by caller after Parquet write
            is_decoy: self.is_decoy,
            charge: self.charge,
            scan_number: self.scan_number,
            apex_rt: self.apex_rt,
            start_rt: self.peak_bounds.start_rt,
            end_rt: self.peak_bounds.end_rt,
            coelution_sum: self.features.coelution_sum,
            score: self.score,
            run_precursor_qvalue: self.run_precursor_qvalue,
            run_peptide_qvalue: self.run_peptide_qvalue,
            run_protein_qvalue: self.run_protein_qvalue,
            experiment_precursor_qvalue: self.experiment_precursor_qvalue,
            experiment_peptide_qvalue: self.experiment_peptide_qvalue,
            experiment_protein_qvalue: self.experiment_protein_qvalue,
            pep: self.pep,
            modified_sequence: Arc::from(self.modified_sequence.as_str()),
        }
    }
}

/// Lightweight stub for FDR scoring and cross-file reconciliation.
///
/// Keeps only the fields needed for Percolator/Mokapot result mapping,
/// target-decoy competition, consensus RT computation, and reconciliation
/// planning.  Full entry data (features, CWT candidates, fragments) is
/// stored in per-file `.scores.parquet` caches and reloaded on demand.
///
/// ~120 bytes inline per entry, versus ~940 bytes for `CoelutionScoredEntry`.
/// `modified_sequence` uses `Arc<str>` for deduplication: with 240 GPF files
/// the same ~3.5M unique sequences appear across all files, so interning
/// reduces heap from ~6 GB to ~90 MB.  `file_name` is not stored -- entries
/// are already keyed by file in `Vec<(String, Vec<FdrEntry>)>`.
#[derive(Debug, Clone)]
pub struct FdrEntry {
    /// Library entry ID (for psm_id construction and library lookup)
    pub entry_id: u32,
    /// Row index in the per-file Parquet cache (for CWT candidate and feature lookup
    /// after the FDR stub Vec is compacted to remove non-passing entries)
    pub parquet_index: u32,
    /// Whether this is a decoy
    pub is_decoy: bool,
    /// Precursor charge (for psm_id and passing_precursors key)
    pub charge: u8,
    /// Scan number at apex (for psm_id construction)
    pub scan_number: u32,
    /// Apex RT in minutes (for consensus computation and blib)
    pub apex_rt: f64,
    /// Peak start RT in minutes (for consensus and reconciliation overlap)
    pub start_rt: f64,
    /// Peak end RT in minutes (for consensus and reconciliation overlap)
    pub end_rt: f64,
    /// Coelution sum score (for simple FDR sort and consensus weighting)
    pub coelution_sum: f64,
    /// SVM discriminant score (for consensus ranking)
    pub score: f64,
    /// Run-level precursor q-value (modified_sequence + charge, mutated by FDR)
    pub run_precursor_qvalue: f64,
    /// Run-level peptide q-value (modified_sequence only, mutated by FDR)
    pub run_peptide_qvalue: f64,
    /// Run-level protein q-value (picked-protein FDR, default 1.0)
    pub run_protein_qvalue: f64,
    /// Experiment-level precursor q-value (modified_sequence + charge, mutated by FDR)
    pub experiment_precursor_qvalue: f64,
    /// Experiment-level peptide q-value (modified_sequence only, mutated by FDR)
    pub experiment_peptide_qvalue: f64,
    /// Experiment-level protein q-value (picked-protein FDR, default 1.0)
    pub experiment_protein_qvalue: f64,
    /// Posterior error probability (mutated by FDR)
    pub pep: f64,
    /// Modified sequence (interned Arc<str> for psm_id, grouping, passing_precursors)
    pub modified_sequence: Arc<str>,
}

impl FdrEntry {
    /// Effective run-level q-value based on the chosen FDR filtering level.
    ///
    /// - `Precursor`: returns run_precursor_qvalue
    /// - `Peptide`: returns run_peptide_qvalue
    /// - `Both`: returns max(run_precursor_qvalue, run_peptide_qvalue)
    pub fn effective_run_qvalue(&self, level: FdrLevel) -> f64 {
        match level {
            FdrLevel::Precursor => self.run_precursor_qvalue,
            FdrLevel::Peptide => self.run_peptide_qvalue,
            FdrLevel::Both => self.run_precursor_qvalue.max(self.run_peptide_qvalue),
        }
    }

    /// Effective experiment-level q-value based on the chosen FDR filtering level.
    ///
    /// - `Precursor`: returns experiment_precursor_qvalue
    /// - `Peptide`: returns experiment_peptide_qvalue
    /// - `Both`: returns max(experiment_precursor_qvalue, experiment_peptide_qvalue)
    pub fn effective_experiment_qvalue(&self, level: FdrLevel) -> f64 {
        match level {
            FdrLevel::Precursor => self.experiment_precursor_qvalue,
            FdrLevel::Peptide => self.experiment_peptide_qvalue,
            FdrLevel::Both => self
                .experiment_precursor_qvalue
                .max(self.experiment_peptide_qvalue),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies that IsolationWindow uses half-open [lower, upper) convention.
    ///
    /// The upper bound is exclusive so entries at shared boundaries between
    /// adjacent windows belong to exactly one window.
    #[test]
    fn test_isolation_window_contains() {
        let window = IsolationWindow::symmetric(500.0, 12.5);
        // lower = 487.5, upper = 512.5
        assert!(window.contains(500.0)); // interior
        assert!(window.contains(487.5)); // lower bound inclusive
        assert!(!window.contains(512.5)); // upper bound exclusive
        assert!(!window.contains(487.4)); // below lower
        assert!(!window.contains(512.6)); // above upper
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
