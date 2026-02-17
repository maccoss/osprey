//! Batch scoring using BLAS-accelerated matrix operations
//!
//! This module provides high-performance spectral scoring by:
//! 1. Preprocessing library entries and spectra once
//! 2. Storing preprocessed data as matrices
//! 3. Using BLAS matrix multiplication for all-vs-all scoring
//!
//! ## LibCosine Scoring (pyXcorrDIA-compatible)
//!
//! The `LibCosineScorer` implements pyXcorrDIA-style library cosine scoring:
//! - **NO binning** - direct ppm-based fragment matching
//! - sqrt(intensity) preprocessing (down-weights dominant peaks)
//! - L2 normalization and cosine angle calculation
//! - Stores mass errors for matched fragments (for calibration)
//!
//! ## BatchScorer (XCorr)
//!
//! The `BatchScorer` uses binning for BLAS-accelerated XCorr calculations:
//! - Unit resolution: 1.0005079 Da bins, 0.4 offset
//! - HRAM: 0.02 Da bins, 0 offset
//!
//! Performance: 10-20× faster than one-at-a-time scoring for typical DIA data.

// Ensure BLAS is linked
extern crate blas_src;
extern crate openblas_src;

use indicatif::{ProgressBar, ProgressStyle};
use ndarray::{Array1, Array2};
use osprey_core::{
    peptide_isotope_cosine, FragmentToleranceConfig, IsotopeEnvelope, LibraryEntry,
    LibraryFragment, MS1Spectrum, Spectrum,
};
use rayon::prelude::*;
use std::collections::HashMap;

use super::has_topn_fragment_match;

use crate::SpectralScorer;

/// Scaling factor for peptide-centric search (DIA default)
pub const PEPTIDE_CENTRIC_SCALING: f64 = 0.0001;

/// Scaling factor for spectrum-centric search
pub const SPECTRUM_CENTRIC_SCALING: f64 = 0.005;

/// Default number of bins for scoring (200-2000 m/z range at 1 Da resolution)
pub const DEFAULT_NUM_BINS: usize = 2000;

// ============================================================================
// LibCosine Scoring (pyXcorrDIA-compatible, NO binning)
// ============================================================================

/// Result of fragment matching for a single library-spectrum pair
#[derive(Debug, Clone)]
pub struct FragmentMatchResult {
    /// Matched experimental intensities (0 if no match)
    pub matched_exp_intensities: Vec<f64>,
    /// Matched library intensities
    pub matched_lib_intensities: Vec<f64>,
    /// Matched library m/z values
    pub matched_lib_mzs: Vec<f64>,
    /// Mass errors for matched fragments (in configured unit: ppm or Th)
    pub mass_errors: Vec<f64>,
    /// Number of fragments that matched (non-zero experimental intensity)
    pub n_matched: usize,
}

/// LibCosine scorer using pyXcorrDIA-style ppm-based matching (NO binning)
///
/// This scorer matches fragments using a tolerance window, NOT binning.
/// For each library fragment, it finds the best matching experimental peak
/// within the configured tolerance (ppm or Da).
///
/// # Intensity Preprocessing
///
/// Both library and experimental intensities are sqrt-transformed:
/// ```text
/// preprocessed_intensity = sqrt(intensity)
/// ```
/// This down-weights dominant peaks and gives more influence to smaller peaks.
///
/// # Scoring
///
/// The score is the cosine angle between preprocessed vectors:
/// ```text
/// score = (exp · lib) / (|exp| × |lib|)
/// ```
#[derive(Debug, Clone, Default)]
pub struct LibCosineScorer {
    /// Fragment tolerance configuration
    pub tolerance: FragmentToleranceConfig,
}

// Default: 10 ppm fragment tolerance

impl LibCosineScorer {
    /// Create a new LibCosine scorer with default 10 ppm tolerance
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a LibCosine scorer with custom tolerance
    pub fn with_tolerance(tolerance: FragmentToleranceConfig) -> Self {
        Self { tolerance }
    }

    /// Create a LibCosine scorer for HRAM data with ppm tolerance
    pub fn hram(ppm: f64) -> Self {
        Self {
            tolerance: FragmentToleranceConfig::hram(ppm),
        }
    }

    /// Create a LibCosine scorer for unit resolution data with Da tolerance
    pub fn unit_resolution(da: f64) -> Self {
        Self {
            tolerance: FragmentToleranceConfig::unit_resolution(da),
        }
    }

    /// Match fragments from a library entry against an experimental spectrum
    ///
    /// For each library fragment, finds the best matching experimental peak
    /// within the tolerance window. If no match is found, uses 0 intensity.
    ///
    /// Returns matched intensities and mass errors for calibration.
    pub fn match_fragments(
        &self,
        library_entry: &LibraryEntry,
        spectrum: &Spectrum,
    ) -> FragmentMatchResult {
        let mut matched_exp = Vec::with_capacity(library_entry.fragments.len());
        let mut matched_lib = Vec::with_capacity(library_entry.fragments.len());
        let mut matched_mzs = Vec::with_capacity(library_entry.fragments.len());
        let mut mass_errors = Vec::with_capacity(library_entry.fragments.len());
        let mut n_matched = 0;

        for lib_frag in &library_entry.fragments {
            let lib_mz = lib_frag.mz;
            let lib_intensity = lib_frag.relative_intensity as f64;

            // Find closest matching experimental peak within tolerance (pyXcorrDIA approach)
            // Select peak with smallest m/z difference, not highest intensity
            let mut best_intensity = 0.0f64;
            let mut best_mz = None;
            let mut best_mz_diff = f64::INFINITY;

            for (&exp_mz, &exp_intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
                if self.tolerance.within_tolerance(lib_mz, exp_mz) {
                    let mz_diff = (exp_mz - lib_mz).abs();
                    if mz_diff < best_mz_diff {
                        best_mz_diff = mz_diff;
                        best_intensity = exp_intensity as f64;
                        best_mz = Some(exp_mz);
                    }
                }
            }

            // Include this fragment (0 intensity if not matched)
            matched_exp.push(best_intensity);
            matched_lib.push(lib_intensity);
            matched_mzs.push(lib_mz);

            // Calculate mass error in configured unit (ppm or Th)
            if let Some(obs_mz) = best_mz {
                mass_errors.push(self.tolerance.mass_error(lib_mz, obs_mz));
                n_matched += 1;
            }
        }

        FragmentMatchResult {
            matched_exp_intensities: matched_exp,
            matched_lib_intensities: matched_lib,
            matched_lib_mzs: matched_mzs,
            mass_errors,
            n_matched,
        }
    }

    /// Calculate LibCosine score with SMZ preprocessing
    ///
    /// # SMZ Preprocessing
    ///
    /// ```text
    /// preprocessed = sqrt(intensity) × mz²
    /// ```
    ///
    /// # Score Calculation
    ///
    /// ```text
    /// score = cosine_angle(exp_preprocessed, lib_preprocessed)
    /// ```
    pub fn calculate_score(&self, match_result: &FragmentMatchResult) -> f64 {
        if match_result.matched_exp_intensities.is_empty() {
            return 0.0;
        }

        // Apply sqrt preprocessing: sqrt(intensity)
        // This down-weights dominant peaks and gives more influence to smaller peaks
        let exp_preprocessed: Vec<f64> = match_result
            .matched_exp_intensities
            .iter()
            .map(|&int| int.sqrt())
            .collect();

        let lib_preprocessed: Vec<f64> = match_result
            .matched_lib_intensities
            .iter()
            .map(|&int| int.sqrt())
            .collect();

        // Calculate cosine angle
        Self::cosine_angle(&exp_preprocessed, &lib_preprocessed)
    }

    /// Score a library entry against a spectrum
    ///
    /// Returns (score, mass_errors) where mass_errors are for matched fragments only.
    pub fn score_with_errors(
        &self,
        library_entry: &LibraryEntry,
        spectrum: &Spectrum,
    ) -> (f64, Vec<f64>) {
        let match_result = self.match_fragments(library_entry, spectrum);
        let score = self.calculate_score(&match_result);
        (score, match_result.mass_errors)
    }

    /// Score a library entry against a spectrum (simple interface)
    pub fn score(&self, library_entry: &LibraryEntry, spectrum: &Spectrum) -> f64 {
        let match_result = self.match_fragments(library_entry, spectrum);
        self.calculate_score(&match_result)
    }

    /// Calculate cosine angle between two vectors
    fn cosine_angle(a: &[f64], b: &[f64]) -> f64 {
        if a.len() != b.len() || a.is_empty() {
            return 0.0;
        }

        let dot: f64 = a.iter().zip(b.iter()).map(|(x, y)| x * y).sum();
        let norm_a: f64 = a.iter().map(|x| x * x).sum::<f64>().sqrt();
        let norm_b: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt();

        if norm_a < 1e-10 || norm_b < 1e-10 {
            return 0.0;
        }

        let cos_angle = dot / (norm_a * norm_b);
        // Clamp to [0, 1] to handle floating point errors
        cos_angle.clamp(0.0, 1.0)
    }
}

/// Result of LibCosine batch scoring for a single library entry
#[derive(Debug, Clone)]
pub struct LibCosineMatch {
    /// Library entry ID
    pub entry_id: u32,
    /// Whether this is a decoy
    pub is_decoy: bool,
    /// Library retention time
    pub library_rt: f64,
    /// Measured retention time (from best matching spectrum)
    pub measured_rt: f64,
    /// Best LibCosine score
    pub score: f64,
    /// MS2 mass errors for matched fragments (in configured unit)
    pub ms2_mass_errors: Vec<f64>,
    /// Average MS2 mass error
    pub avg_ms2_error: Option<f64>,
    /// Number of matched fragments
    pub n_matched_fragments: usize,
    /// Total number of library fragments
    pub n_library_fragments: usize,
    /// Best matching spectrum index
    pub best_spectrum_idx: usize,
    /// Library precursor m/z
    pub library_precursor_mz: f64,
}

/// Preprocessed library ready for batch scoring
///
/// Stores SMZ-preprocessed and L2-normalized vectors for each library entry
/// as a single matrix for efficient BLAS operations.
/// Uses f32 for memory efficiency - sufficient precision for normalized intensity values.
#[derive(Debug, Clone)]
pub struct PreprocessedLibrary {
    /// Matrix of preprocessed library vectors (n_entries × n_bins)
    /// Each row is an SMZ-preprocessed, L2-normalized spectrum
    pub matrix: Array2<f32>,
    /// Library entry IDs corresponding to each row
    pub entry_ids: Vec<u32>,
    /// Mapping from entry ID to matrix row index
    pub id_to_row: HashMap<u32, usize>,
    /// Number of bins used
    pub num_bins: usize,
    /// Minimum m/z for binning
    pub min_mz: f64,
    /// Bin width
    pub bin_width: f64,
}

impl PreprocessedLibrary {
    /// Create a preprocessed library from library entries
    ///
    /// Applies SMZ preprocessing (sqrt(intensity) × m/z²) and L2 normalization
    /// to each library entry's fragments.
    pub fn from_entries(entries: &[LibraryEntry]) -> Self {
        Self::from_entries_with_config(entries, DEFAULT_NUM_BINS, 200.0, 1.0005)
    }

    /// Create with custom binning configuration
    pub fn from_entries_with_config(
        entries: &[LibraryEntry],
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
    ) -> Self {
        Self::from_entries_with_config_and_context(entries, num_bins, min_mz, bin_width, None)
    }

    /// Create with custom binning configuration and logging context
    pub fn from_entries_with_config_and_context(
        entries: &[LibraryEntry],
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
        context: Option<&str>,
    ) -> Self {
        let n_entries = entries.len();

        // Preprocess all entries in parallel
        let preprocessed: Vec<(u32, Array1<f32>)> = entries
            .par_iter()
            .filter(|e| !e.fragments.is_empty())
            .map(|entry| {
                let binned = Self::preprocess_library_entry(entry, num_bins, min_mz, bin_width);
                (entry.id, binned)
            })
            .collect();

        // Build matrix and mappings
        let n_valid = preprocessed.len();
        let mut matrix: Array2<f32> = Array2::zeros((n_valid, num_bins));
        let mut entry_ids = Vec::with_capacity(n_valid);
        let mut id_to_row = HashMap::with_capacity(n_valid);

        for (row_idx, (entry_id, binned)) in preprocessed.into_iter().enumerate() {
            matrix.row_mut(row_idx).assign(&binned);
            entry_ids.push(entry_id);
            id_to_row.insert(entry_id, row_idx);
        }

        if let Some(ctx) = context {
            log::debug!(
                "PreprocessedLibrary [{}]: {} entries -> {} valid rows ({}×{} matrix)",
                ctx,
                n_entries,
                n_valid,
                n_valid,
                num_bins
            );
        } else {
            log::debug!(
                "PreprocessedLibrary: {} entries -> {} valid rows ({}×{} matrix)",
                n_entries,
                n_valid,
                n_valid,
                num_bins
            );
        }

        Self {
            matrix,
            entry_ids,
            id_to_row,
            num_bins,
            min_mz,
            bin_width,
        }
    }

    /// Preprocess a single library entry (for internal use)
    /// Uses f32 for memory efficiency - source intensities are already f32.
    fn preprocess_library_entry(
        entry: &LibraryEntry,
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
    ) -> Array1<f32> {
        let mut binned: Array1<f32> = Array1::zeros(num_bins);

        // Apply sqrt preprocessing: sqrt(intensity)
        for frag in &entry.fragments {
            if frag.mz >= min_mz && frag.mz < min_mz + num_bins as f64 * bin_width {
                let bin = ((frag.mz - min_mz) / bin_width) as usize;
                if bin < num_bins {
                    let sqrt_int = frag.relative_intensity.sqrt();
                    binned[bin] += sqrt_int;
                }
            }
        }

        // L2 normalize
        let norm: f32 = binned.mapv(|x: f32| x * x).sum();
        let norm = norm.sqrt();
        if norm > 1e-10 {
            binned /= norm;
        }

        binned
    }

    /// Get the number of preprocessed entries
    pub fn len(&self) -> usize {
        self.entry_ids.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.entry_ids.is_empty()
    }

    /// Get a subset of the library for specific entry IDs
    pub fn subset(&self, entry_ids: &[u32]) -> PreprocessedLibrary {
        let rows: Vec<usize> = entry_ids
            .iter()
            .filter_map(|id| self.id_to_row.get(id).copied())
            .collect();

        let n_rows = rows.len();
        let mut matrix: Array2<f32> = Array2::zeros((n_rows, self.num_bins));
        let mut new_entry_ids = Vec::with_capacity(n_rows);
        let mut new_id_to_row = HashMap::with_capacity(n_rows);

        for (new_idx, &old_idx) in rows.iter().enumerate() {
            matrix.row_mut(new_idx).assign(&self.matrix.row(old_idx));
            let id = self.entry_ids[old_idx];
            new_entry_ids.push(id);
            new_id_to_row.insert(id, new_idx);
        }

        PreprocessedLibrary {
            matrix,
            entry_ids: new_entry_ids,
            id_to_row: new_id_to_row,
            num_bins: self.num_bins,
            min_mz: self.min_mz,
            bin_width: self.bin_width,
        }
    }
}

/// Preprocessed spectra ready for batch scoring
///
/// Stores preprocessed spectra as a matrix for efficient BLAS operations.
/// Uses f32 for memory efficiency - sufficient precision for normalized intensity values.
#[derive(Debug, Clone)]
pub struct PreprocessedSpectra {
    /// Matrix of preprocessed spectra (n_spectra × n_bins)
    /// Each row is an SMZ-preprocessed, L2-normalized spectrum
    pub matrix: Array2<f32>,
    /// Spectrum indices corresponding to each row
    pub spectrum_indices: Vec<usize>,
    /// Retention times for each spectrum
    pub retention_times: Vec<f64>,
    /// Mapping from spectrum index to matrix row
    pub idx_to_row: HashMap<usize, usize>,
    /// Number of bins used
    pub num_bins: usize,
    /// Minimum m/z for binning
    pub min_mz: f64,
    /// Bin width
    pub bin_width: f64,
}

impl PreprocessedSpectra {
    /// Create preprocessed spectra from a slice of spectra
    pub fn from_spectra(spectra: &[Spectrum]) -> Self {
        Self::from_spectra_with_config(spectra, DEFAULT_NUM_BINS, 200.0, 1.0005)
    }

    /// Create with custom binning configuration
    pub fn from_spectra_with_config(
        spectra: &[Spectrum],
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
    ) -> Self {
        Self::from_spectra_with_config_and_context(spectra, num_bins, min_mz, bin_width, None)
    }

    /// Create with custom binning configuration and logging context
    pub fn from_spectra_with_config_and_context(
        spectra: &[Spectrum],
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
        context: Option<&str>,
    ) -> Self {
        let n_spectra = spectra.len();

        // Preprocess all spectra in parallel
        let preprocessed: Vec<(usize, f64, Array1<f32>)> = spectra
            .par_iter()
            .enumerate()
            .filter(|(_, s)| !s.mzs.is_empty())
            .map(|(idx, spectrum)| {
                let binned = Self::preprocess_spectrum(spectrum, num_bins, min_mz, bin_width);
                (idx, spectrum.retention_time, binned)
            })
            .collect();

        // Build matrix and mappings
        let n_valid = preprocessed.len();
        let mut matrix: Array2<f32> = Array2::zeros((n_valid, num_bins));
        let mut spectrum_indices = Vec::with_capacity(n_valid);
        let mut retention_times = Vec::with_capacity(n_valid);
        let mut idx_to_row = HashMap::with_capacity(n_valid);

        for (row_idx, (spec_idx, rt, binned)) in preprocessed.into_iter().enumerate() {
            matrix.row_mut(row_idx).assign(&binned);
            spectrum_indices.push(spec_idx);
            retention_times.push(rt);
            idx_to_row.insert(spec_idx, row_idx);
        }

        if let Some(ctx) = context {
            log::debug!(
                "PreprocessedSpectra [{}]: {} spectra -> {} valid rows ({}×{} matrix)",
                ctx,
                n_spectra,
                n_valid,
                n_valid,
                num_bins
            );
        } else {
            log::debug!(
                "PreprocessedSpectra: {} spectra -> {} valid rows ({}×{} matrix)",
                n_spectra,
                n_valid,
                n_valid,
                num_bins
            );
        }

        Self {
            matrix,
            spectrum_indices,
            retention_times,
            idx_to_row,
            num_bins,
            min_mz,
            bin_width,
        }
    }

    /// Preprocess a single spectrum with SMZ transformation and L2 normalization
    /// Uses f32 for memory efficiency - source intensities are already f32.
    fn preprocess_spectrum(
        spectrum: &Spectrum,
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
    ) -> Array1<f32> {
        let mut binned: Array1<f32> = Array1::zeros(num_bins);

        // Apply sqrt preprocessing: sqrt(intensity)
        for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            if mz >= min_mz && mz < min_mz + num_bins as f64 * bin_width {
                let bin = ((mz - min_mz) / bin_width) as usize;
                if bin < num_bins {
                    let sqrt_int = intensity.sqrt();
                    binned[bin] += sqrt_int;
                }
            }
        }

        // L2 normalize
        let norm: f32 = binned.mapv(|x: f32| x * x).sum();
        let norm = norm.sqrt();
        if norm > 1e-10 {
            binned /= norm;
        }

        binned
    }

    /// Get the number of preprocessed spectra
    pub fn len(&self) -> usize {
        self.spectrum_indices.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.spectrum_indices.is_empty()
    }

    /// Get a subset of spectra by original indices
    pub fn subset(&self, indices: &[usize]) -> PreprocessedSpectra {
        let rows: Vec<usize> = indices
            .iter()
            .filter_map(|idx| self.idx_to_row.get(idx).copied())
            .collect();

        let n_rows = rows.len();
        let mut matrix: Array2<f32> = Array2::zeros((n_rows, self.num_bins));
        let mut new_indices = Vec::with_capacity(n_rows);
        let mut new_rts = Vec::with_capacity(n_rows);
        let mut new_idx_to_row = HashMap::with_capacity(n_rows);

        for (new_row, &old_row) in rows.iter().enumerate() {
            matrix.row_mut(new_row).assign(&self.matrix.row(old_row));
            let idx = self.spectrum_indices[old_row];
            new_indices.push(idx);
            new_rts.push(self.retention_times[old_row]);
            new_idx_to_row.insert(idx, new_row);
        }

        PreprocessedSpectra {
            matrix,
            spectrum_indices: new_indices,
            retention_times: new_rts,
            idx_to_row: new_idx_to_row,
            num_bins: self.num_bins,
            min_mz: self.min_mz,
            bin_width: self.bin_width,
        }
    }
}

/// Batch scorer using BLAS-accelerated matrix operations
///
/// Computes LibCosine scores for all library-spectrum pairs in a single
/// matrix multiplication, providing 10-20× speedup over one-at-a-time scoring.
#[derive(Debug, Clone)]
pub struct BatchScorer {
    /// Number of bins for scoring
    num_bins: usize,
    /// Minimum m/z
    min_mz: f64,
    /// Bin width
    bin_width: f64,
}

impl Default for BatchScorer {
    fn default() -> Self {
        Self {
            num_bins: DEFAULT_NUM_BINS,
            min_mz: 200.0,
            bin_width: 1.0005,
        }
    }
}

impl BatchScorer {
    /// Create a new batch scorer with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Create with custom binning configuration
    pub fn with_config(num_bins: usize, min_mz: f64, bin_width: f64) -> Self {
        Self {
            num_bins,
            min_mz,
            bin_width,
        }
    }

    /// Preprocess a library for batch scoring
    pub fn preprocess_library(&self, entries: &[LibraryEntry]) -> PreprocessedLibrary {
        self.preprocess_library_with_context(entries, None)
    }

    /// Preprocess a library for batch scoring with logging context
    pub fn preprocess_library_with_context(
        &self,
        entries: &[LibraryEntry],
        context: Option<&str>,
    ) -> PreprocessedLibrary {
        PreprocessedLibrary::from_entries_with_config_and_context(
            entries,
            self.num_bins,
            self.min_mz,
            self.bin_width,
            context,
        )
    }

    /// Preprocess spectra for batch scoring
    pub fn preprocess_spectra(&self, spectra: &[Spectrum]) -> PreprocessedSpectra {
        self.preprocess_spectra_with_context(spectra, None)
    }

    /// Preprocess spectra for batch scoring with logging context
    pub fn preprocess_spectra_with_context(
        &self,
        spectra: &[Spectrum],
        context: Option<&str>,
    ) -> PreprocessedSpectra {
        PreprocessedSpectra::from_spectra_with_config_and_context(
            spectra,
            self.num_bins,
            self.min_mz,
            self.bin_width,
            context,
        )
    }

    /// Calculate LibCosine scores for all library-spectrum pairs
    ///
    /// Uses BLAS matrix multiplication: scores = library_matrix @ spectra_matrix.T
    ///
    /// Returns a matrix of shape (n_library × n_spectra) where score[i][j] is
    /// the LibCosine score between library entry i and spectrum j.
    ///
    /// # Performance
    /// For 5000 library entries × 500 spectra, this is ~20× faster than
    /// computing each pair individually.
    ///
    /// Note: Internal computation uses f32 for memory efficiency, returns f64 for compatibility.
    pub fn score_all(
        &self,
        library: &PreprocessedLibrary,
        spectra: &PreprocessedSpectra,
    ) -> Array2<f64> {
        // Matrix multiplication: (n_library × n_bins) · (n_bins × n_spectra)
        // Result: (n_library × n_spectra) in f32, then convert to f64
        let scores_f32 = library.matrix.dot(&spectra.matrix.t());
        scores_f32.mapv(|x| x as f64)
    }

    /// Calculate LibCosine scores and return best match for each library entry
    ///
    /// Returns: Vec of (entry_id, best_spectrum_idx, best_score, best_rt)
    pub fn find_best_matches(
        &self,
        library: &PreprocessedLibrary,
        spectra: &PreprocessedSpectra,
    ) -> Vec<(u32, usize, f64, f64)> {
        if library.is_empty() || spectra.is_empty() {
            return Vec::new();
        }

        let scores = self.score_all(library, spectra);

        // For each library entry, find the spectrum with highest score
        library
            .entry_ids
            .par_iter()
            .enumerate()
            .filter_map(|(lib_row, &entry_id)| {
                let row = scores.row(lib_row);
                let (best_col, &best_score) =
                    row.iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1))?;

                if best_score > 0.0 {
                    let spec_idx = spectra.spectrum_indices[best_col];
                    let rt = spectra.retention_times[best_col];
                    Some((entry_id, spec_idx, best_score, rt))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Calculate LibCosine scores with RT filtering
    ///
    /// Only considers library-spectrum pairs where the RT difference is within tolerance.
    /// This is more efficient when most pairs would be filtered out anyway.
    ///
    /// Returns: Vec of (entry_id, best_spectrum_idx, best_score, best_rt)
    pub fn find_best_matches_with_rt_filter(
        &self,
        library: &PreprocessedLibrary,
        spectra: &PreprocessedSpectra,
        library_rts: &[f64],
        rt_tolerance: f64,
    ) -> Vec<(u32, usize, f64, f64)> {
        if library.is_empty() || spectra.is_empty() {
            return Vec::new();
        }

        // Compute full score matrix
        let scores = self.score_all(library, spectra);

        // For each library entry, find best spectrum within RT tolerance
        library
            .entry_ids
            .par_iter()
            .enumerate()
            .filter_map(|(lib_row, &entry_id)| {
                // Find the library RT (indexed by row position in preprocessed library)
                let lib_rt = if lib_row < library_rts.len() {
                    library_rts[lib_row]
                } else {
                    return None;
                };

                let row = scores.row(lib_row);
                let mut best_score = 0.0f64;
                let mut best_col = 0usize;

                for (col, &score) in row.iter().enumerate() {
                    let spec_rt = spectra.retention_times[col];
                    if (spec_rt - lib_rt).abs() <= rt_tolerance && score > best_score {
                        best_score = score;
                        best_col = col;
                    }
                }

                if best_score > 0.0 {
                    let spec_idx = spectra.spectrum_indices[best_col];
                    let rt = spectra.retention_times[best_col];
                    Some((entry_id, spec_idx, best_score, rt))
                } else {
                    None
                }
            })
            .collect()
    }

    /// Score a single library entry against all spectra (for targeted queries)
    ///
    /// Returns: Vec of (spectrum_idx, score, rt) sorted by score descending
    pub fn score_entry_vs_all(
        &self,
        entry: &LibraryEntry,
        spectra: &PreprocessedSpectra,
    ) -> Vec<(usize, f64, f64)> {
        let entry_vec = PreprocessedLibrary::preprocess_library_entry(
            entry,
            self.num_bins,
            self.min_mz,
            self.bin_width,
        );

        // Score against all spectra (f32 matrix multiplication)
        let scores = spectra.matrix.dot(&entry_vec);

        let mut results: Vec<(usize, f64, f64)> = scores
            .iter()
            .enumerate()
            .filter(|(_, &score)| score > 0.0)
            .map(|(row, &score)| {
                let spec_idx = spectra.spectrum_indices[row];
                let rt = spectra.retention_times[row];
                (spec_idx, score as f64, rt) // Cast f32 score to f64 for API compatibility
            })
            .collect();

        // Sort by score descending
        results.sort_by(|a, b| b.1.total_cmp(&a.1));
        results
    }
}

/// Result of batch scoring for calibration
#[derive(Debug, Clone)]
pub struct CalibrationMatch {
    /// Library entry ID
    pub entry_id: u32,
    /// Whether this is a decoy
    pub is_decoy: bool,
    /// Library retention time
    pub library_rt: f64,
    /// Measured retention time (from best matching spectrum)
    pub measured_rt: f64,
    /// Primary score (XCorr, same as xcorr_score)
    pub score: f64,
    /// MS1 (precursor) m/z error (in configured unit: ppm or Th)
    pub ms1_error: Option<f64>,
    /// Library precursor m/z (theoretical)
    pub library_precursor_mz: f64,
    /// Observed precursor m/z (from spectrum isolation window center)
    pub observed_precursor_mz: Option<f64>,
    /// MS2 fragment mass errors (in configured unit: ppm or Th)
    pub ms2_mass_errors: Vec<f64>,
    /// Average MS2 mass error (only for matched fragments)
    pub avg_ms2_error: Option<f64>,
    /// Number of matched fragments
    pub n_matched_fragments: usize,
    /// Total number of library fragments
    pub n_library_fragments: usize,
    /// XCorr score (Comet-style cross-correlation)
    pub xcorr_score: f64,
    /// Comet-style E-value (calculated from XCorr survival function)
    /// Lower is better - represents expected number of random matches at this score
    pub evalue: f64,
    /// Isotope cosine score (MS1 isotope envelope match)
    pub isotope_cosine_score: Option<f64>,
    /// Peptide sequence
    pub sequence: String,
    /// Charge state
    pub charge: u8,
    /// Best matching spectrum scan number
    pub scan_number: u32,
    /// X!Tandem-style hyperscore (ln(n_b!) + ln(n_y!) + Σln(I+1))
    pub hyperscore: f64,
    /// Number of matched b-ions (from hyperscore computation)
    pub n_b_ions: u32,
    /// Number of matched y-ions (from hyperscore computation)
    pub n_y_ions: u32,

    // === NEW: ML Scoring Features ===
    /// Fragment XIC co-elution correlation (sum of positive correlations)
    pub correlation_score: f64,
    /// LibCosine spectral similarity at apex
    pub libcosine_apex: f64,
    /// Number of top-6 library fragments matched at apex (0-6)
    pub top6_matched_apex: u8,
    /// Hyperscore calculated at apex
    pub hyperscore_apex: f64,
    /// Signal-to-noise ratio of reference XIC peak
    /// Signal = apex - background_mean, Noise = background_SD
    pub signal_to_noise: f64,
    /// Peak width in minutes (from detected XIC boundaries, None for non-XIC paths)
    pub peak_width_minutes: Option<f64>,
    /// Linear discriminant score (weighted combination of features)
    pub discriminant_score: f64,
    /// Posterior error probability from KDE
    pub posterior_error: f64,
    /// Q-value from target-decoy competition
    pub q_value: f64,
}

/// Calculate Comet-style E-value from XCorr scores
///
/// Fits a linear regression to the survival function of XCorr scores:
/// log10(cumulative_count) = intercept + slope × xcorr
///
/// Then E-value = 10^(slope × best_xcorr + intercept)
///
/// This matches Comet's approach but operates per-peptide instead of per-spectrum:
/// - Comet: one spectrum scored against many peptides → E-value per peptide
/// - Osprey: one peptide scored against many spectra → E-value per peptide
///
/// # Arguments
/// * `all_xcorr_scores` - All XCorr scores for this peptide against spectra in RT window
/// * `best_xcorr` - The best XCorr score (what we want to evaluate)
///
/// # Returns
/// E-value (lower is better, represents expected random matches at this score level)
pub fn calculate_evalue_from_xcorr_distribution(all_xcorr_scores: &[f64], best_xcorr: f64) -> f64 {
    const BAD_EVALUE: f64 = 999.0; // Comet's sentinel for bad fits
    const HISTO_SIZE: usize = 152; // Same as Comet
    const MIN_POINTS_FOR_FIT: usize = 10; // Need reasonable number for stable fit

    // Need enough scores for a meaningful distribution
    if all_xcorr_scores.len() < MIN_POINTS_FOR_FIT {
        return BAD_EVALUE;
    }

    // Build histogram like Comet: bin_index = (int)(xcorr * 10.0 * 0.005 + 0.5) = (int)(xcorr * 0.05 + 0.5)
    // This gives ~20 XCorr units per bin (bin 1 = XCorr 10-30, bin 2 = XCorr 30-50, etc.)
    // For our typical XCorr range of 0-3, we need finer resolution
    // Scale factor: bin_index = (int)(xcorr * 50.0 + 0.5), so 1 bin = 0.02 XCorr units
    const SCALE: f64 = 50.0;
    let mut histogram = [0i32; HISTO_SIZE];

    for &score in all_xcorr_scores {
        if score > 0.0 {
            let bin = ((score * SCALE) + 0.5) as usize;
            if bin < HISTO_SIZE {
                histogram[bin] += 1;
            } else {
                histogram[HISTO_SIZE - 1] += 1; // Cap at max bin
            }
        }
    }

    // Find maximum correlation score index (highest non-zero bin)
    let mut max_corr = 0usize;
    for i in (0..HISTO_SIZE - 1).rev() {
        if histogram[i] > 0 {
            max_corr = i;
            break;
        }
    }

    // Find the "next" correlation index - where we start seeing gaps in the histogram
    // This helps exclude outlier scores from the fit (Comet's approach)
    let mut next_corr = max_corr;
    let mut found_first_nonzero = false;
    for i in 0..max_corr {
        if histogram[i] > 0 {
            found_first_nonzero = true;
        }
        // Look for consecutive zeros after we've seen data
        if histogram[i] == 0
            && found_first_nonzero
            && i >= 5
            && i + 1 < HISTO_SIZE
            && (histogram[i + 1] == 0 || i + 1 == max_corr)
        {
            next_corr = if i > 0 { i - 1 } else { 0 };
            break;
        }
    }

    // If no gap found, use max_corr - 1 to exclude the outlier
    if next_corr == max_corr && max_corr >= 5 {
        for i in (max_corr.saturating_sub(5)..=max_corr).rev() {
            if histogram[i] == 0 {
                next_corr = i;
                break;
            }
        }
        if next_corr == max_corr {
            next_corr = max_corr.saturating_sub(1);
        }
    }

    // Build cumulative distribution from next_corr down
    let mut cumulative = [0.0f64; HISTO_SIZE];
    cumulative[next_corr] = histogram[next_corr] as f64;
    for i in (0..next_corr).rev() {
        cumulative[i] = cumulative[i + 1] + histogram[i] as f64;
        // Zero out positions where histogram was zero (Comet's approach)
        if histogram[i + 1] == 0 {
            cumulative[i + 1] = 0.0;
        }
    }

    // Convert to log10, handling zeros
    let mut log_cumulative = [0.0f64; HISTO_SIZE];
    for i in (0..=next_corr).rev() {
        if cumulative[i] > 0.0 {
            log_cumulative[i] = cumulative[i].log10();
        } else if i + 1 < HISTO_SIZE && cumulative[i + 1] > 0.0 {
            log_cumulative[i] = cumulative[i + 1].log10();
        }
    }

    // Determine start point for regression (exclude low bins)
    let mut start_corr = next_corr.saturating_sub(5);
    // Count zeros in the fit range and adjust
    let num_zeros: usize = (start_corr..=next_corr)
        .filter(|&i| log_cumulative[i] == 0.0)
        .count();
    start_corr = start_corr.saturating_sub(num_zeros);

    // Linear regression on log10(cumulative) vs bin index
    // Fit: log10(cumulative) = a + b * bin_index
    // We expect b < 0 (higher bin = higher XCorr = fewer matches = lower cumulative)
    let mut sum_x = 0.0f64;
    let mut sum_y = 0.0f64;
    let mut n_points = 0usize;

    for i in start_corr..=next_corr {
        if histogram[i] > 0 {
            sum_x += i as f64;
            sum_y += log_cumulative[i];
            n_points += 1;
        }
    }

    if n_points < 3 {
        return BAD_EVALUE;
    }

    let mean_x = sum_x / n_points as f64;
    let mean_y = sum_y / n_points as f64;

    // Calculate slope and intercept
    let mut sx = 0.0f64; // Sum of squared deviations
    let mut sxy = 0.0f64; // Sum of cross products

    for (i, &log_cum_val) in log_cumulative
        .iter()
        .enumerate()
        .take(next_corr + 1)
        .skip(start_corr)
    {
        if log_cum_val > 0.0 {
            let dx = i as f64 - mean_x;
            let dy = log_cum_val - mean_y;
            sx += dx * dx;
            sxy += dx * dy;
        }
    }

    if sx < 1e-10 {
        return BAD_EVALUE;
    }

    let slope = sxy / sx; // b in log10(cumulative) = a + b * bin_index
    let intercept = mean_y - slope * mean_x; // a

    // Comet checks: if slope >= 0, the fit is bad (we expect negative slope)
    // Higher bins should have fewer matches, so slope should be negative
    if slope >= 0.0 {
        return BAD_EVALUE;
    }

    // Convert best_xcorr to bin scale for the formula
    // E-value = 10^(slope * bin_index + intercept)
    // But we need to account for our SCALE factor
    // bin_index = xcorr * SCALE, so: E-value = 10^(slope * SCALE * xcorr + intercept)
    let log_evalue = slope * SCALE * best_xcorr + intercept;

    // Clamp to reasonable range and compute E-value
    let evalue = 10.0f64.powf(log_evalue.clamp(-10.0, 10.0));

    // Clamp to Comet's range [1e-10, 999]
    evalue.clamp(1e-10, BAD_EVALUE)
}

/// Calculate E-value from a score distribution using survival function.
///
/// Same Comet-style approach as `calculate_evalue_from_xcorr_distribution` but with
/// a caller-specified scale factor instead of a hardcoded XCorr scale.
/// Works for any positive score where higher is better (XCorr, hyperscore, etc.).
///
/// # Arguments
/// * `all_scores` - All scores for this peptide against spectra in RT window
/// * `best_score` - The best score (what we want to evaluate)
/// * `scale` - Fixed scale factor: bin = score × scale. Must be chosen so the expected
///   score range fits within 152 histogram bins. Examples:
///   - XCorr (0-3 range): scale=50.0 → bins 0-150
///   - Hyperscore (0-300 range): scale=0.5 → bins 0-150
///
/// **Important**: scale must be FIXED (not derived from per-peptide data) to preserve
/// discrimination between targets and decoys.
pub fn calculate_evalue_from_score_distribution(
    all_scores: &[f64],
    best_score: f64,
    scale: f64,
) -> f64 {
    const BAD_EVALUE: f64 = 999.0;
    const HISTO_SIZE: usize = 152;
    const MIN_POINTS_FOR_FIT: usize = 10;

    if all_scores.len() < MIN_POINTS_FOR_FIT {
        return BAD_EVALUE;
    }

    let mut histogram = [0i32; HISTO_SIZE];
    for &score in all_scores {
        if score > 0.0 {
            let bin = ((score * scale) + 0.5) as usize;
            if bin < HISTO_SIZE {
                histogram[bin] += 1;
            } else {
                histogram[HISTO_SIZE - 1] += 1;
            }
        }
    }

    // Find highest non-zero bin
    let mut max_corr = 0usize;
    for i in (0..HISTO_SIZE - 1).rev() {
        if histogram[i] > 0 {
            max_corr = i;
            break;
        }
    }

    // Find "next_corr" — exclude outlier tail (Comet approach)
    let mut next_corr = max_corr;
    let mut found_first_nonzero = false;
    for i in 0..max_corr {
        if histogram[i] > 0 {
            found_first_nonzero = true;
        }
        if histogram[i] == 0
            && found_first_nonzero
            && i >= 5
            && i + 1 < HISTO_SIZE
            && (histogram[i + 1] == 0 || i + 1 == max_corr)
        {
            next_corr = if i > 0 { i - 1 } else { 0 };
            break;
        }
    }

    if next_corr == max_corr && max_corr >= 5 {
        for i in (max_corr.saturating_sub(5)..=max_corr).rev() {
            if histogram[i] == 0 {
                next_corr = i;
                break;
            }
        }
        if next_corr == max_corr {
            next_corr = max_corr.saturating_sub(1);
        }
    }

    // Build cumulative distribution (survival function)
    let mut cumulative = [0.0f64; HISTO_SIZE];
    cumulative[next_corr] = histogram[next_corr] as f64;
    for i in (0..next_corr).rev() {
        cumulative[i] = cumulative[i + 1] + histogram[i] as f64;
        if histogram[i + 1] == 0 {
            cumulative[i + 1] = 0.0;
        }
    }

    // Convert to log10
    let mut log_cumulative = [0.0f64; HISTO_SIZE];
    for i in (0..=next_corr).rev() {
        if cumulative[i] > 0.0 {
            log_cumulative[i] = cumulative[i].log10();
        } else if i + 1 < HISTO_SIZE && cumulative[i + 1] > 0.0 {
            log_cumulative[i] = cumulative[i + 1].log10();
        }
    }

    // Determine regression range
    let mut start_corr = next_corr.saturating_sub(5);
    let num_zeros: usize = (start_corr..=next_corr)
        .filter(|&i| log_cumulative[i] == 0.0)
        .count();
    start_corr = start_corr.saturating_sub(num_zeros);

    // Linear regression: log10(cumulative) = a + b * bin_index
    let mut sum_x = 0.0f64;
    let mut sum_y = 0.0f64;
    let mut n_points = 0usize;

    for i in start_corr..=next_corr {
        if histogram[i] > 0 {
            sum_x += i as f64;
            sum_y += log_cumulative[i];
            n_points += 1;
        }
    }

    if n_points < 3 {
        return BAD_EVALUE;
    }

    let mean_x = sum_x / n_points as f64;
    let mean_y = sum_y / n_points as f64;

    let mut sx = 0.0f64;
    let mut sxy = 0.0f64;

    for (i, &log_cum_val) in log_cumulative
        .iter()
        .enumerate()
        .take(next_corr + 1)
        .skip(start_corr)
    {
        if log_cum_val > 0.0 {
            let dx = i as f64 - mean_x;
            let dy = log_cum_val - mean_y;
            sx += dx * dx;
            sxy += dx * dy;
        }
    }

    if sx < 1e-10 {
        return BAD_EVALUE;
    }

    let slope = sxy / sx;
    let intercept = mean_y - slope * mean_x;

    // Expect negative slope (higher bins = fewer matches)
    if slope >= 0.0 {
        return BAD_EVALUE;
    }

    // E-value = 10^(slope * scale * best_score + intercept)
    let log_evalue = slope * scale * best_score + intercept;
    let evalue = 10.0f64.powf(log_evalue.clamp(-10.0, 10.0));
    evalue.clamp(1e-10, BAD_EVALUE)
}

/// Paired target-decoy calibration result for debugging CSV
#[derive(Debug, Clone)]
pub struct PairedCalibrationResult {
    /// Target entry ID
    pub target_entry_id: u32,
    /// Charge state
    pub charge: u8,
    /// Target peptide sequence
    pub target_sequence: String,
    /// Decoy peptide sequence
    pub decoy_sequence: String,
    /// Target primary score (XCorr)
    pub target_score: f64,
    /// Decoy primary score (XCorr)
    pub decoy_score: f64,
    /// Winning primary score (max of target and decoy)
    pub winning_score: f64,
    /// Target XCorr score
    pub target_xcorr: f64,
    /// Decoy XCorr score
    pub decoy_xcorr: f64,
    /// Winning XCorr score (max of target and decoy)
    pub winning_xcorr: f64,
    /// Target E-value (Comet-style, calculated from XCorr)
    pub target_evalue: f64,
    /// Decoy E-value (Comet-style, calculated from XCorr)
    pub decoy_evalue: f64,
    /// Winning E-value (min of target and decoy - lower is better)
    pub winning_evalue: f64,
    /// True if target wins by E-value (target_evalue < decoy_evalue)
    pub target_wins_evalue: bool,
    /// Target isotope cosine score
    pub target_isotope_score: Option<f64>,
    /// Decoy isotope cosine score
    pub decoy_isotope_score: Option<f64>,
    /// Target precursor error (ppm)
    pub target_precursor_error_ppm: Option<f64>,
    /// Decoy precursor error (ppm)
    pub decoy_precursor_error_ppm: Option<f64>,
    /// Target measured RT (minutes)
    pub target_rt: f64,
    /// Decoy measured RT (minutes)
    pub decoy_rt: f64,
    /// Library RT (for target, same for decoy since it's a reversed sequence)
    pub library_rt: f64,
    /// Expected RT (mapped from library RT via calibration, if available)
    pub expected_rt: Option<f64>,
    /// Target delta RT: |expected_rt - target_rt| if calibrated, else |library_rt - target_rt|
    pub target_delta_rt: f64,
    /// Decoy delta RT: |expected_rt - decoy_rt| if calibrated, else |library_rt - decoy_rt|
    pub decoy_delta_rt: f64,
    /// Target number of matched fragments (library fragments that matched experimental peaks)
    pub target_matched_frags: usize,
    /// Decoy number of matched fragments (library fragments that matched experimental peaks)
    pub decoy_matched_frags: usize,
    /// True if target wins the competition (target score > decoy score)
    pub target_wins: bool,
    /// Target hyperscore
    pub target_hyperscore: f64,
    /// Decoy hyperscore
    pub decoy_hyperscore: f64,
    /// Winning hyperscore (max of target and decoy)
    pub winning_hyperscore: f64,
    /// Target matched b-ions
    pub target_n_b: u32,
    /// Target matched y-ions
    pub target_n_y: u32,
    /// Decoy matched b-ions
    pub decoy_n_b: u32,
    /// Decoy matched y-ions
    pub decoy_n_y: u32,
}

/// Pair calibration matches by target-decoy ID relationship
///
/// Decoys are linked to targets via ID bitmask: decoy_id = target_id | 0x80000000
pub fn pair_calibration_matches(
    matches: &[CalibrationMatch],
    expected_rt_fn: Option<&dyn Fn(f64) -> f64>,
) -> Vec<PairedCalibrationResult> {
    use std::collections::HashMap;

    // Build lookup maps by base ID (strip high bit for decoys)
    let mut target_matches: HashMap<u32, &CalibrationMatch> = HashMap::new();
    let mut decoy_matches: HashMap<u32, &CalibrationMatch> = HashMap::new();

    for m in matches {
        let base_id = m.entry_id & 0x7FFFFFFF;
        if m.is_decoy {
            decoy_matches.insert(base_id, m);
        } else {
            target_matches.insert(base_id, m);
        }
    }

    // Create paired results
    let mut paired: Vec<PairedCalibrationResult> = Vec::new();
    for (&target_id, &target_match) in &target_matches {
        if let Some(&decoy_match) = decoy_matches.get(&target_id) {
            // Compute expected RT if calibration function provided
            let expected_rt = expected_rt_fn.map(|f| f(target_match.library_rt));

            // Compute delta RT for target and decoy
            // If calibration available, use |expected - measured|
            // Otherwise, use |library - measured| which is still useful
            let reference_rt = expected_rt.unwrap_or(target_match.library_rt);
            let target_delta_rt = (target_match.measured_rt - reference_rt).abs();
            let decoy_delta_rt = (decoy_match.measured_rt - reference_rt).abs();

            // Compute winning scores (what FDR calculation actually uses)
            let target_wins = target_match.score > decoy_match.score;
            let winning_score = target_match.score.max(decoy_match.score);
            let winning_xcorr = target_match.xcorr_score.max(decoy_match.xcorr_score);

            // Use E-values from CalibrationMatch (calculated from survival function)
            // Lower E-value is better
            let target_evalue = target_match.evalue;
            let decoy_evalue = decoy_match.evalue;
            let winning_evalue = target_evalue.min(decoy_evalue);
            let target_wins_evalue = target_evalue < decoy_evalue;

            paired.push(PairedCalibrationResult {
                target_entry_id: target_id,
                charge: target_match.charge,
                target_sequence: target_match.sequence.clone(),
                decoy_sequence: decoy_match.sequence.clone(),
                target_score: target_match.score,
                decoy_score: decoy_match.score,
                winning_score,
                target_xcorr: target_match.xcorr_score,
                decoy_xcorr: decoy_match.xcorr_score,
                winning_xcorr,
                target_evalue,
                decoy_evalue,
                winning_evalue,
                target_wins_evalue,
                target_isotope_score: target_match.isotope_cosine_score,
                decoy_isotope_score: decoy_match.isotope_cosine_score,
                target_precursor_error_ppm: target_match.ms1_error,
                decoy_precursor_error_ppm: decoy_match.ms1_error,
                target_rt: target_match.measured_rt,
                decoy_rt: decoy_match.measured_rt,
                library_rt: target_match.library_rt,
                expected_rt,
                target_delta_rt,
                decoy_delta_rt,
                target_matched_frags: target_match.n_matched_fragments,
                decoy_matched_frags: decoy_match.n_matched_fragments,
                target_wins,
                target_hyperscore: target_match.hyperscore,
                decoy_hyperscore: decoy_match.hyperscore,
                winning_hyperscore: target_match.hyperscore.max(decoy_match.hyperscore),
                target_n_b: target_match.n_b_ions,
                target_n_y: target_match.n_y_ions,
                decoy_n_b: decoy_match.n_b_ions,
                decoy_n_y: decoy_match.n_y_ions,
            });
        }
    }

    // Sort by WINNING E-value ascending (lower is better) - best matches first
    paired.sort_by(|a, b| a.winning_evalue.total_cmp(&b.winning_evalue));

    paired
}

/// Sample peptides from the library for calibration using 2D stratified sampling.
///
/// For large libraries (hundreds of thousands of entries), scoring all entries for calibration
/// is too slow. Instead, we sample a representative subset using a 2D grid over RT × m/z space.
/// This ensures good coverage across both dimensions, critical for:
/// - LOESS RT calibration (needs points spanning full RT range)
/// - Mass calibration (needs points spanning full m/z range)
///
/// Returns both sampled targets AND their paired decoys.
///
/// # Arguments
/// * `library` - Full library (targets + decoys)
/// * `sample_size` - Number of target peptides to sample (0 = use all)
/// * `seed` - Deterministic offset for reproducibility
pub fn sample_library_for_calibration(
    library: &[LibraryEntry],
    sample_size: usize,
    seed: u64,
) -> Vec<LibraryEntry> {
    use std::collections::HashSet;

    if sample_size == 0 {
        return library.to_vec();
    }

    let targets: Vec<&LibraryEntry> = library.iter().filter(|e| !e.is_decoy).collect();
    let decoys: Vec<&LibraryEntry> = library.iter().filter(|e| e.is_decoy).collect();

    if targets.len() <= sample_size {
        log::info!(
            "Calibration sampling: library has {} targets (<= requested {}), using all {} entries",
            targets.len(),
            sample_size,
            library.len()
        );
        return library.to_vec();
    }

    // Build target_id → decoy map (decoy_id = target_id | 0x80000000)
    let decoy_map: HashMap<u32, &LibraryEntry> =
        decoys.iter().map(|d| (d.id & 0x7FFFFFFF, *d)).collect();

    // 2D stratified sampling: divide RT × m/z space into a grid
    // Use ~sqrt(sample_size)/2 bins per axis for good 2D coverage
    let bins_per_axis = ((sample_size as f64).sqrt() / 2.0).ceil().max(5.0) as usize;

    let (rt_min, rt_max) = targets.iter().fold((f64::MAX, f64::MIN), |(lo, hi), e| {
        (lo.min(e.retention_time), hi.max(e.retention_time))
    });
    let (mz_min, mz_max) = targets.iter().fold((f64::MAX, f64::MIN), |(lo, hi), e| {
        (lo.min(e.precursor_mz), hi.max(e.precursor_mz))
    });

    let rt_range = (rt_max - rt_min).max(1e-6);
    let mz_range = (mz_max - mz_min).max(1e-6);
    let rt_bin_width = rt_range / bins_per_axis as f64;
    let mz_bin_width = mz_range / bins_per_axis as f64;

    // Assign each target to a 2D grid cell
    let mut grid: Vec<Vec<Vec<usize>>> = vec![vec![Vec::new(); bins_per_axis]; bins_per_axis];
    for (i, target) in targets.iter().enumerate() {
        let rt_bin = ((target.retention_time - rt_min) / rt_bin_width).floor() as usize;
        let mz_bin = ((target.precursor_mz - mz_min) / mz_bin_width).floor() as usize;
        let rt_bin = rt_bin.min(bins_per_axis - 1);
        let mz_bin = mz_bin.min(bins_per_axis - 1);
        grid[rt_bin][mz_bin].push(i);
    }

    // Count non-empty cells and compute per-cell quota
    let n_occupied: usize = grid
        .iter()
        .flat_map(|row| row.iter())
        .filter(|cell| !cell.is_empty())
        .count();
    let per_cell = if n_occupied > 0 {
        sample_size / n_occupied
    } else {
        1
    };
    let per_cell = per_cell.max(1);

    // Deterministic stride sampling from each cell
    let offset = seed as usize;
    let mut sampled_ids: HashSet<u32> = HashSet::new();
    let mut sampled: Vec<LibraryEntry> = Vec::with_capacity(sample_size * 2);

    for row in grid.iter().take(bins_per_axis) {
        for cell in row.iter().take(bins_per_axis) {
            if cell.is_empty() {
                continue;
            }

            let n_take = cell.len().min(per_cell);
            let stride = (cell.len() / n_take).max(1);
            let cell_offset = offset % cell.len().max(1);

            for j in 0..n_take {
                let idx = (cell_offset + j * stride) % cell.len();
                let target = targets[cell[idx]];

                if sampled_ids.contains(&target.id) {
                    continue;
                }
                sampled_ids.insert(target.id);
                sampled.push(target.clone());

                if let Some(decoy) = decoy_map.get(&target.id) {
                    sampled.push((*decoy).clone());
                }
            }
        }
    }

    // If we under-sampled (sparse grid), do a second pass with increased per-cell quota
    if sampled_ids.len() < sample_size {
        let remaining = sample_size - sampled_ids.len();
        let extra_per_cell = (remaining / n_occupied.max(1)).max(1);

        'outer: for row in grid.iter().take(bins_per_axis) {
            for cell in row.iter().take(bins_per_axis) {
                if cell.is_empty() {
                    continue;
                }

                let mut added = 0;
                for &target_idx in cell {
                    if added >= extra_per_cell {
                        break;
                    }
                    let target = targets[target_idx];
                    if sampled_ids.contains(&target.id) {
                        continue;
                    }
                    sampled_ids.insert(target.id);
                    sampled.push(target.clone());
                    if let Some(decoy) = decoy_map.get(&target.id) {
                        sampled.push((*decoy).clone());
                    }
                    added += 1;
                }

                if sampled_ids.len() >= sample_size {
                    break 'outer;
                }
            }
        }
    }

    let n_decoys = sampled.len() - sampled_ids.len();
    log::info!(
        "Calibration sampling: {} targets + {} decoys from {} total ({}\u{00d7}{} grid, {} occupied cells)",
        sampled_ids.len(),
        n_decoys,
        library.len(),
        bins_per_axis,
        bins_per_axis,
        n_occupied,
    );

    sampled
}

/// Group spectra by isolation window bounds
///
/// Returns: Vec of ((lower_bound, upper_bound), spectrum_indices)
pub fn group_spectra_by_isolation_window(spectra: &[Spectrum]) -> Vec<((f64, f64), Vec<usize>)> {
    let mut windows: HashMap<(i64, i64), Vec<usize>> = HashMap::new();

    for (idx, spec) in spectra.iter().enumerate() {
        let iso = &spec.isolation_window;
        // Round window bounds to avoid floating point issues (0.1 m/z precision)
        let lower_key = (iso.lower_bound() * 10.0) as i64;
        let upper_key = (iso.upper_bound() * 10.0) as i64;

        windows.entry((lower_key, upper_key)).or_default().push(idx);
    }

    // Convert back to f64 windows with spectrum indices
    windows
        .into_iter()
        .map(|((lower, upper), indices)| ((lower as f64 / 10.0, upper as f64 / 10.0), indices))
        .collect()
}

/// Run windowed batch calibration scoring
///
/// This is the memory-efficient version for large libraries. Instead of
/// preprocessing the entire library at once (which requires ~50GB for 3M entries),
/// it groups spectra by isolation window and only preprocesses library entries
/// that fall within each window.
///
/// For a typical 3mz DIA window with a 500 m/z library range and 3M precursors,
/// each window only contains ~18K precursors (3M × 3/500), reducing memory by ~170×.
///
/// # Arguments
/// * `library` - Full library (targets + decoys)
/// * `spectra` - All MS2 spectra
/// * `fragment_tolerance` - Fragment m/z tolerance (for computing MS1 error in correct unit)
/// * `rt_tolerance` - RT tolerance for matching
///
/// # Returns
/// Calibration matches sorted by score descending
pub fn run_windowed_calibration_scoring(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    fragment_tolerance: FragmentToleranceConfig,
    rt_tolerance: f64,
) -> Vec<CalibrationMatch> {
    // Group spectra by isolation window
    let window_groups = group_spectra_by_isolation_window(spectra);

    if window_groups.is_empty() {
        log::warn!("No spectra with isolation windows found");
        return Vec::new();
    }

    log::info!(
        "Grouped {} spectra into {} isolation windows for windowed scoring",
        spectra.len(),
        window_groups.len()
    );

    // Build library precursor m/z index for fast lookup
    let id_to_entry: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    // Process each window in parallel
    let all_matches: Vec<Vec<CalibrationMatch>> = window_groups
        .par_iter()
        .map(|((lower, upper), spectrum_indices)| {
            // Find library entries whose precursor m/z falls in this window
            let window_entries: Vec<&LibraryEntry> = library
                .iter()
                .filter(|e| e.precursor_mz >= *lower && e.precursor_mz <= *upper)
                .collect();

            if window_entries.is_empty() {
                return Vec::new();
            }

            // Extract spectra for this window
            let window_spectra: Vec<Spectrum> = spectrum_indices
                .iter()
                .map(|&idx| spectra[idx].clone())
                .collect();

            if window_spectra.is_empty() {
                return Vec::new();
            }

            // Preprocess only the entries in this window
            let scorer = BatchScorer::new();
            let cloned_entries: Vec<LibraryEntry> =
                window_entries.iter().map(|e| (*e).clone()).collect();
            let window_context = format!("window {:.1}-{:.1} m/z", lower, upper);
            let preprocessed_library =
                scorer.preprocess_library_with_context(&cloned_entries, Some(&window_context));
            let preprocessed_spectra =
                scorer.preprocess_spectra_with_context(&window_spectra, Some(&window_context));

            // Get library RTs in the same order as preprocessed rows
            let library_rts: Vec<f64> = preprocessed_library
                .entry_ids
                .iter()
                .filter_map(|id| id_to_entry.get(id).map(|e| e.retention_time))
                .collect();

            // Find best matches with RT filtering
            let matches = scorer.find_best_matches_with_rt_filter(
                &preprocessed_library,
                &preprocessed_spectra,
                &library_rts,
                rt_tolerance,
            );

            // Convert to CalibrationMatch with MS1 error (in configured unit)
            matches
                .into_iter()
                .filter_map(|(entry_id, spec_idx, score, measured_rt)| {
                    let entry = id_to_entry.get(&entry_id)?;

                    // Get observed precursor m/z from spectrum (isolation window center)
                    let observed_precursor_mz = if spec_idx < window_spectra.len() {
                        Some(window_spectra[spec_idx].isolation_window.center)
                    } else {
                        None
                    };

                    // Calculate MS1 error in configured unit (ppm for HRAM, Th for unit resolution)
                    let ms1_error = observed_precursor_mz
                        .map(|obs_mz| fragment_tolerance.mass_error(entry.precursor_mz, obs_mz));

                    // Get scan number from spectrum
                    let scan_number = if spec_idx < window_spectra.len() {
                        window_spectra[spec_idx].scan_number
                    } else {
                        0
                    };

                    // For binned BatchScorer, use simple E-value approximation
                    // (we don't have the full XCorr distribution here)
                    let evalue = (-score * 20.0).exp();

                    Some(CalibrationMatch {
                        entry_id,
                        is_decoy: entry.is_decoy,
                        library_rt: entry.retention_time,
                        measured_rt,
                        score,
                        ms1_error,
                        library_precursor_mz: entry.precursor_mz,
                        observed_precursor_mz,
                        // Binned BatchScorer doesn't compute MS2 errors
                        ms2_mass_errors: Vec::new(),
                        avg_ms2_error: None,
                        n_matched_fragments: 0,
                        n_library_fragments: entry.fragments.len(),
                        // For binned BatchScorer, the score IS the XCorr-like score
                        xcorr_score: score,
                        evalue,
                        isotope_cosine_score: None,
                        sequence: entry.sequence.clone(),
                        charge: entry.charge,
                        scan_number,
                        hyperscore: 0.0,
                        n_b_ions: 0,
                        n_y_ions: 0,

                        // ML scoring features (placeholder values)
                        correlation_score: 0.0,
                        libcosine_apex: 0.0,
                        top6_matched_apex: 0,
                        hyperscore_apex: 0.0,
                        signal_to_noise: 0.0,
                        peak_width_minutes: None,
                        discriminant_score: score,
                        posterior_error: 0.0,
                        q_value: 1.0,
                    })
                })
                .collect()
        })
        .collect();

    // Flatten and deduplicate (keep best score per entry)
    let mut best_matches: HashMap<u32, CalibrationMatch> = HashMap::new();
    for matches in all_matches {
        for m in matches {
            best_matches
                .entry(m.entry_id)
                .and_modify(|existing| {
                    if m.score > existing.score {
                        *existing = m.clone();
                    }
                })
                .or_insert(m);
        }
    }

    // Sort by score descending
    let mut results: Vec<CalibrationMatch> = best_matches.into_values().collect();
    results.sort_by(|a, b| b.score.total_cmp(&a.score));

    log::info!(
        "Windowed batch scoring complete: {} unique matches from {} windows",
        results.len(),
        window_groups.len()
    );

    results
}

/// Run batch calibration scoring
///
/// Scores all library entries (targets and decoys) against all spectra
/// within the specified RT tolerance, using BLAS-accelerated matrix operations.
///
/// Returns matches sorted by score descending.
pub fn run_calibration_scoring(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    fragment_tolerance: FragmentToleranceConfig,
    rt_tolerance: f64,
) -> Vec<CalibrationMatch> {
    let scorer = BatchScorer::new();

    // Preprocess library and spectra
    log::info!(
        "Preprocessing {} library entries for batch scoring",
        library.len()
    );
    let preprocessed_library = scorer.preprocess_library(library);

    log::info!("Preprocessing {} spectra for batch scoring", spectra.len());
    let preprocessed_spectra = scorer.preprocess_spectra(spectra);

    // Get library RTs in the same order as preprocessed rows
    let library_rts: Vec<f64> = preprocessed_library
        .entry_ids
        .iter()
        .filter_map(|id| {
            library
                .iter()
                .find(|e| e.id == *id)
                .map(|e| e.retention_time)
        })
        .collect();

    log::info!(
        "Computing LibCosine scores ({} × {} = {} pairs via BLAS)",
        preprocessed_library.len(),
        preprocessed_spectra.len(),
        preprocessed_library.len() * preprocessed_spectra.len()
    );

    // Find best matches with RT filtering
    let matches = scorer.find_best_matches_with_rt_filter(
        &preprocessed_library,
        &preprocessed_spectra,
        &library_rts,
        rt_tolerance,
    );

    // Convert to CalibrationMatch with MS1 PPM error
    let id_to_entry: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    let mut results: Vec<CalibrationMatch> = matches
        .into_iter()
        .filter_map(|(entry_id, spec_idx, score, measured_rt)| {
            let entry = id_to_entry.get(&entry_id)?;

            // Get observed precursor m/z from spectrum (isolation window center)
            let observed_precursor_mz = if spec_idx < spectra.len() {
                Some(spectra[spec_idx].isolation_window.center)
            } else {
                None
            };

            // Calculate MS1 error in configured unit (ppm for HRAM, Th for unit resolution)
            let ms1_error = observed_precursor_mz
                .map(|obs_mz| fragment_tolerance.mass_error(entry.precursor_mz, obs_mz));

            // Get scan number from spectrum
            let scan_number = if spec_idx < spectra.len() {
                spectra[spec_idx].scan_number
            } else {
                0
            };

            // For binned BatchScorer, use simple E-value approximation
            let evalue = (-score * 20.0).exp();

            Some(CalibrationMatch {
                entry_id,
                is_decoy: entry.is_decoy,
                library_rt: entry.retention_time,
                measured_rt,
                score,
                ms1_error,
                library_precursor_mz: entry.precursor_mz,
                observed_precursor_mz,
                // Binned BatchScorer doesn't compute MS2 errors
                ms2_mass_errors: Vec::new(),
                avg_ms2_error: None,
                n_matched_fragments: 0,
                n_library_fragments: entry.fragments.len(),
                // For binned BatchScorer, the score IS the XCorr-like score
                xcorr_score: score,
                evalue,
                isotope_cosine_score: None,
                sequence: entry.sequence.clone(),
                charge: entry.charge,
                scan_number,
                hyperscore: 0.0,
                n_b_ions: 0,
                n_y_ions: 0,

                // ML scoring features (placeholder values)
                correlation_score: 0.0,
                libcosine_apex: 0.0,
                top6_matched_apex: 0,
                hyperscore_apex: 0.0,
                signal_to_noise: 0.0,
                peak_width_minutes: None,
                discriminant_score: score,
                posterior_error: 0.0,
                q_value: 1.0,
            })
        })
        .collect();

    // Sort by score descending
    results.sort_by(|a, b| b.score.total_cmp(&a.score));

    log::info!("Batch scoring complete: {} matches found", results.len());

    results
}

/// MS1 spectrum index for efficient nearest-neighbor lookup
///
/// This is a simplified interface that wraps the MS1Index from osprey-io
/// for use in calibration scoring.
pub trait MS1SpectrumLookup: Sync {
    /// Find the nearest MS1 spectrum to a given retention time
    fn find_nearest(&self, retention_time: f64) -> Option<&MS1Spectrum>;
}

/// Run calibration scoring with XCorr + MS1 isotope extraction
///
/// This is the primary calibration scoring function. For each library entry:
/// 1. Filter by precursor window, RT tolerance, and top-N fragment match
/// 2. Calculate XCorr across all passing spectra to find the best RT
/// 3. Collect MS2 mass errors from top-N fragment matches at best spectrum
/// 4. Calculate Comet-style E-value from the XCorr score distribution
/// 5. Extract M+0 isotope peak from MS1 for accurate mass calibration
///
/// XCorr always uses unit resolution bins (2001 bins, 1.0005 m/z) for speed,
/// regardless of whether the data is unit resolution or HRAM.
///
/// # Arguments
/// * `library` - Full library (targets + decoys)
/// * `spectra` - All MS2 spectra
/// * `ms1_index` - MS1 spectrum index for isotope extraction (optional)
/// * `fragment_tolerance` - Fragment m/z tolerance (ppm or Da)
/// * `precursor_tolerance_ppm` - Precursor tolerance for isotope matching (default: 10 ppm)
/// * `rt_tolerance` - RT tolerance for matching
///
/// # Returns
/// Calibration matches sorted by score descending
pub fn run_xcorr_calibration_scoring<M: MS1SpectrumLookup>(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    ms1_index: Option<&M>,
    fragment_tolerance: FragmentToleranceConfig,
    precursor_tolerance_ppm: f64,
    rt_tolerance: f64,
) -> Vec<CalibrationMatch> {
    // Group spectra by isolation window
    let window_groups = group_spectra_by_isolation_window(spectra);

    if window_groups.is_empty() {
        log::warn!("No spectra with isolation windows found");
        return Vec::new();
    }

    let ms1_available = ms1_index.is_some();
    log::info!(
        "XCorr calibration scoring: {} spectra in {} windows (MS1 available: {})",
        spectra.len(),
        window_groups.len(),
        ms1_available
    );
    log::info!("  - XCorr: unit resolution bins (2001 bins) for all data types");
    log::info!("  - E-value: Comet-style survival function from XCorr distribution");

    // Always use unit resolution bins for calibration XCorr (fast: 2001 bins vs 100K for HRAM)
    let xcorr_scorer = SpectralScorer::new();

    // Set up progress bar for calibration scoring
    let pb = ProgressBar::new(library.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} entries",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Process each window in parallel
    // MS1SpectrumLookup requires Sync, so ms1_index can be shared across threads
    let all_matches: Vec<CalibrationMatch> = if ms1_available {
        // Parallel processing with MS1 extraction
        let ms1 = ms1_index.unwrap();

        window_groups
            .par_iter()
            .flat_map(|((lower, upper), spectrum_indices)| {
                // Find library entries whose precursor m/z falls in this window
                let window_entries: Vec<&LibraryEntry> = library
                    .iter()
                    .filter(|e| e.precursor_mz >= *lower && e.precursor_mz <= *upper)
                    .collect();

                if window_entries.is_empty() {
                    return Vec::new();
                }

                // Extract spectra for this window
                let window_spectra: Vec<&Spectrum> =
                    spectrum_indices.iter().map(|&idx| &spectra[idx]).collect();

                if window_spectra.is_empty() {
                    return Vec::new();
                }

                // === BATCH XCORR PREPROCESSING (pyXcorrDIA approach) ===
                // Preprocess all spectra for XCorr ONCE per window (uses f32 for memory efficiency)
                // Spectra are shared across many library entries, so upfront preprocessing is efficient
                let spec_xcorr_preprocessed: Vec<Vec<f32>> = window_spectra
                    .iter()
                    .map(|spec| xcorr_scorer.preprocess_spectrum_for_xcorr(spec))
                    .collect();

                // NOTE: Library entries are preprocessed lazily (only after passing RT + topN filters)
                // This avoids wasting ~8KB per entry for entries that never get scored

                // Score each library entry against all spectra in this window
                let mut window_matches: Vec<CalibrationMatch> = Vec::new();

                for entry in &window_entries {
                    // === FIRST PASS: find spectra passing RT + topN filters ===
                    let candidate_indices: Vec<usize> = window_spectra
                        .iter()
                        .enumerate()
                        .filter(|(_, spec)| {
                            let rt_diff = (spec.retention_time - entry.retention_time).abs();
                            if rt_diff > rt_tolerance {
                                return false;
                            }
                            has_topn_fragment_match(
                                &entry.fragments,
                                &spec.mzs,
                                fragment_tolerance.tolerance,
                                fragment_tolerance.unit,
                            )
                        })
                        .map(|(idx, _)| idx)
                        .collect();

                    // Skip XCorr preprocessing entirely if no spectra pass filters
                    if candidate_indices.is_empty() {
                        pb.inc(1);
                        continue;
                    }

                    // Preprocess this library entry for XCorr only now that we know it's needed
                    let entry_xcorr_vec = xcorr_scorer.preprocess_library_for_xcorr(entry);

                    // === COLLECT ALL XCORR SCORES FOR E-VALUE CALCULATION ===
                    // Score against spectra that passed RT + topN filters
                    let mut all_xcorr_scores: Vec<f64> = Vec::new();
                    let mut best_xcorr = 0.0f64;
                    let mut best_local_idx = 0usize; // Index into window_spectra
                    let mut best_spec_idx = 0usize; // Index into global spectra
                    let mut best_measured_rt = 0.0f64;

                    for &local_idx in &candidate_indices {
                        let spec = window_spectra[local_idx];

                        // Calculate XCorr from preprocessed vectors
                        let xcorr = SpectralScorer::xcorr_from_preprocessed(
                            &spec_xcorr_preprocessed[local_idx],
                            &entry_xcorr_vec,
                        );

                        // Collect ALL scores for E-value calculation
                        all_xcorr_scores.push(xcorr);

                        if xcorr > best_xcorr {
                            best_xcorr = xcorr;
                            best_local_idx = local_idx;
                            best_spec_idx = spectrum_indices[local_idx];
                            best_measured_rt = spec.retention_time;
                        }
                    }

                    if best_xcorr > 0.0 {
                        let best_spec = window_spectra[best_local_idx];

                        // Collect MS2 mass errors from top-N fragments at best spectrum
                        let (_has_match, ms2_errors) = super::topn_fragment_match_with_errors(
                            &entry.fragments,
                            &best_spec.mzs,
                            fragment_tolerance.tolerance,
                            fragment_tolerance.unit,
                        );
                        let n_matched = ms2_errors.len();

                        // Extract M+0 from MS1 spectrum (pyXcorrDIA approach)
                        let (ms1_error, observed_precursor_mz, isotope_cosine_score) =
                            if let Some(ms1_spec) = ms1.find_nearest(best_measured_rt) {
                                let envelope = IsotopeEnvelope::extract(
                                    ms1_spec,
                                    entry.precursor_mz,
                                    entry.charge,
                                    precursor_tolerance_ppm,
                                );

                                if envelope.has_m0() {
                                    let isotope_score = peptide_isotope_cosine(
                                        &entry.sequence,
                                        &envelope.intensities,
                                    );
                                    let error = envelope.m0_observed_mz.map(|obs_mz| {
                                        fragment_tolerance.mass_error(entry.precursor_mz, obs_mz)
                                    });
                                    (error, envelope.m0_observed_mz, isotope_score)
                                } else {
                                    (None, None, None)
                                }
                            } else {
                                let iso_center = spectra[best_spec_idx].isolation_window.center;
                                let error =
                                    fragment_tolerance.mass_error(entry.precursor_mz, iso_center);
                                (Some(error), Some(iso_center), None)
                            };

                        let avg_ms2_error = if !ms2_errors.is_empty() {
                            Some(ms2_errors.iter().sum::<f64>() / ms2_errors.len() as f64)
                        } else {
                            None
                        };

                        let evalue =
                            calculate_evalue_from_xcorr_distribution(&all_xcorr_scores, best_xcorr);

                        window_matches.push(CalibrationMatch {
                            entry_id: entry.id,
                            is_decoy: entry.is_decoy,
                            library_rt: entry.retention_time,
                            measured_rt: best_measured_rt,
                            score: best_xcorr,
                            ms1_error,
                            library_precursor_mz: entry.precursor_mz,
                            observed_precursor_mz,
                            ms2_mass_errors: ms2_errors,
                            avg_ms2_error,
                            n_matched_fragments: n_matched,
                            n_library_fragments: entry.fragments.len(),
                            xcorr_score: best_xcorr,
                            evalue,
                            isotope_cosine_score,
                            sequence: entry.sequence.clone(),
                            charge: entry.charge,
                            scan_number: spectra[best_spec_idx].scan_number,
                            hyperscore: 0.0,
                            n_b_ions: 0,
                            n_y_ions: 0,

                            // ML scoring features (placeholder values)
                            correlation_score: 0.0,
                            libcosine_apex: 0.0,
                            top6_matched_apex: 0,
                            hyperscore_apex: 0.0,
                            signal_to_noise: 0.0,
                            peak_width_minutes: None,
                            discriminant_score: best_xcorr,
                            posterior_error: 0.0,
                            q_value: 1.0,
                        });
                    }
                    pb.inc(1);
                }

                window_matches
            })
            .collect()
    } else {
        // Parallel processing without MS1 (fall back to isolation window)
        window_groups
            .par_iter()
            .flat_map(|((lower, upper), spectrum_indices)| {
                let window_entries: Vec<&LibraryEntry> = library
                    .iter()
                    .filter(|e| e.precursor_mz >= *lower && e.precursor_mz <= *upper)
                    .collect();

                if window_entries.is_empty() {
                    return Vec::new();
                }

                let window_spectra: Vec<&Spectrum> =
                    spectrum_indices.iter().map(|&idx| &spectra[idx]).collect();

                if window_spectra.is_empty() {
                    return Vec::new();
                }

                // === BATCH XCORR PREPROCESSING (pyXcorrDIA approach) ===
                // Preprocess all spectra for XCorr ONCE per window (uses f32 for memory efficiency)
                // Spectra are shared across many library entries, so upfront preprocessing is efficient
                let spec_xcorr_preprocessed: Vec<Vec<f32>> = window_spectra
                    .iter()
                    .map(|spec| xcorr_scorer.preprocess_spectrum_for_xcorr(spec))
                    .collect();

                // NOTE: Library entries are preprocessed lazily (only after passing RT + topN filters)
                // This avoids wasting ~8KB per entry for entries that never get scored

                let mut window_matches: Vec<CalibrationMatch> = Vec::new();

                for entry in &window_entries {
                    // === FIRST PASS: find spectra passing RT + topN filters ===
                    let candidate_indices: Vec<usize> = window_spectra
                        .iter()
                        .enumerate()
                        .filter(|(_, spec)| {
                            let rt_diff = (spec.retention_time - entry.retention_time).abs();
                            if rt_diff > rt_tolerance {
                                return false;
                            }
                            has_topn_fragment_match(
                                &entry.fragments,
                                &spec.mzs,
                                fragment_tolerance.tolerance,
                                fragment_tolerance.unit,
                            )
                        })
                        .map(|(idx, _)| idx)
                        .collect();

                    // Skip XCorr preprocessing entirely if no spectra pass filters
                    if candidate_indices.is_empty() {
                        pb.inc(1);
                        continue;
                    }

                    // Preprocess this library entry for XCorr only now that we know it's needed
                    let entry_xcorr_vec = xcorr_scorer.preprocess_library_for_xcorr(entry);

                    // === COLLECT ALL XCORR SCORES FOR E-VALUE CALCULATION ===
                    let mut all_xcorr_scores: Vec<f64> = Vec::new();
                    let mut best_xcorr = 0.0f64;
                    let mut best_local_idx = 0usize;
                    let mut best_spec_idx = 0usize;
                    let mut best_measured_rt = 0.0f64;

                    for &local_idx in &candidate_indices {
                        let spec = window_spectra[local_idx];

                        let xcorr = SpectralScorer::xcorr_from_preprocessed(
                            &spec_xcorr_preprocessed[local_idx],
                            &entry_xcorr_vec,
                        );

                        // Collect ALL scores for E-value calculation
                        all_xcorr_scores.push(xcorr);

                        if xcorr > best_xcorr {
                            best_xcorr = xcorr;
                            best_local_idx = local_idx;
                            best_spec_idx = spectrum_indices[local_idx];
                            best_measured_rt = spec.retention_time;
                        }
                    }

                    if best_xcorr > 0.0 {
                        let best_spec = window_spectra[best_local_idx];

                        // Collect MS2 mass errors from top-N fragments at best spectrum
                        let (_has_match, ms2_errors) = super::topn_fragment_match_with_errors(
                            &entry.fragments,
                            &best_spec.mzs,
                            fragment_tolerance.tolerance,
                            fragment_tolerance.unit,
                        );
                        let n_matched = ms2_errors.len();

                        let observed_precursor_mz =
                            Some(spectra[best_spec_idx].isolation_window.center);
                        let ms1_error = observed_precursor_mz.map(|obs_mz| {
                            fragment_tolerance.mass_error(entry.precursor_mz, obs_mz)
                        });

                        let avg_ms2_error = if !ms2_errors.is_empty() {
                            Some(ms2_errors.iter().sum::<f64>() / ms2_errors.len() as f64)
                        } else {
                            None
                        };

                        let evalue =
                            calculate_evalue_from_xcorr_distribution(&all_xcorr_scores, best_xcorr);

                        window_matches.push(CalibrationMatch {
                            entry_id: entry.id,
                            is_decoy: entry.is_decoy,
                            library_rt: entry.retention_time,
                            measured_rt: best_measured_rt,
                            score: best_xcorr,
                            ms1_error,
                            library_precursor_mz: entry.precursor_mz,
                            observed_precursor_mz,
                            ms2_mass_errors: ms2_errors,
                            avg_ms2_error,
                            n_matched_fragments: n_matched,
                            n_library_fragments: entry.fragments.len(),
                            xcorr_score: best_xcorr,
                            evalue,
                            isotope_cosine_score: None,
                            sequence: entry.sequence.clone(),
                            charge: entry.charge,
                            scan_number: spectra[best_spec_idx].scan_number,
                            hyperscore: 0.0,
                            n_b_ions: 0,
                            n_y_ions: 0,

                            // ML scoring features (placeholder values)
                            correlation_score: 0.0,
                            libcosine_apex: 0.0,
                            top6_matched_apex: 0,
                            hyperscore_apex: 0.0,
                            signal_to_noise: 0.0,
                            peak_width_minutes: None,
                            discriminant_score: best_xcorr,
                            posterior_error: 0.0,
                            q_value: 1.0,
                        });
                    }
                    pb.inc(1);
                }

                window_matches
            })
            .collect()
    };

    pb.finish_with_message("Done");

    // Deduplicate (keep best XCorr per entry, since XCorr is used for RT selection)
    let mut best_matches: HashMap<u32, CalibrationMatch> = HashMap::new();
    for m in all_matches {
        best_matches
            .entry(m.entry_id)
            .and_modify(|existing| {
                if m.xcorr_score > existing.xcorr_score {
                    *existing = m.clone();
                }
            })
            .or_insert(m);
    }

    // Sort by XCorr descending (since XCorr is used for RT selection)
    let mut results: Vec<CalibrationMatch> = best_matches.into_values().collect();
    results.sort_by(|a, b| b.xcorr_score.total_cmp(&a.xcorr_score));

    // Log scoring statistics
    let targets = results.iter().filter(|m| !m.is_decoy).count();
    let decoys = results.iter().filter(|m| m.is_decoy).count();

    log::info!(
        "Scoring complete: {} matches ({} targets, {} decoys)",
        results.len(),
        targets,
        decoys
    );

    results
}

/// Minimum number of spectra required for co-elution scoring.
/// Fewer than 3 time points makes Pearson correlation unreliable.
/// Minimum number of spectra required for coelution analysis
pub const MIN_COELUTION_SPECTRA: usize = 3;

/// Count how many of the top-6 library fragments have matches at apex
///
/// Helper function for calibration and coelution feature extraction.
pub fn count_top6_matched_at_apex(
    library_fragments: &[LibraryFragment],
    spectrum_mzs: &[f64],
    tolerance: FragmentToleranceConfig,
) -> u8 {
    // Get top 6 fragments by intensity
    let top6_indices = super::get_top_n_fragment_indices(library_fragments, 6);

    // Count matches using binary search
    let mut count = 0u8;
    for &idx in &top6_indices {
        let lib_mz = library_fragments[idx].mz;
        if super::has_match_within_tolerance(
            lib_mz,
            spectrum_mzs,
            tolerance.tolerance,
            tolerance.unit,
        ) {
            count += 1;
        }
    }
    count
}

/// 4. Correlate each other fragment XIC against the smoothed reference (Pearson)
/// 5. Score = sum of positive correlations
/// 6. Apex RT = RT at the maximum of the reference XIC
/// 7. Collect MS2 mass errors at the apex spectrum for calibration
/// 8. Extract M+0 isotope peak from MS1 for accurate mass calibration
#[allow(clippy::too_many_arguments)]
pub fn run_coelution_calibration_scoring<M: MS1SpectrumLookup>(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    ms1_index: Option<&M>,
    fragment_tolerance: FragmentToleranceConfig,
    precursor_tolerance_ppm: f64,
    rt_tolerance: f64,
    expected_rt_fn: Option<&(dyn Fn(f64) -> f64 + Sync)>,
    xcorr_scorer: Option<&SpectralScorer>,
) -> Vec<CalibrationMatch> {
    // Group spectra by isolation window
    let window_groups = group_spectra_by_isolation_window(spectra);

    if window_groups.is_empty() {
        log::warn!("No spectra with isolation windows found");
        return Vec::new();
    }

    let ms1_available = ms1_index.is_some();
    log::info!(
        "Co-elution calibration scoring: {} spectra in {} windows (MS1 available: {})",
        spectra.len(),
        window_groups.len(),
        ms1_available
    );
    log::info!(
        "  - Fragment co-elution: top-6 XIC correlation ({} {} tolerance)",
        fragment_tolerance.tolerance,
        match fragment_tolerance.unit {
            osprey_core::ToleranceUnit::Ppm => "ppm",
            osprey_core::ToleranceUnit::Mz => "Th",
        }
    );

    let tol_da = match fragment_tolerance.unit {
        osprey_core::ToleranceUnit::Mz => fragment_tolerance.tolerance,
        osprey_core::ToleranceUnit::Ppm => 0.0, // ppm computed per-fragment
    };
    let tol_ppm = match fragment_tolerance.unit {
        osprey_core::ToleranceUnit::Ppm => fragment_tolerance.tolerance,
        osprey_core::ToleranceUnit::Mz => 0.0,
    };

    // Set up progress bar for calibration scoring
    let pb = ProgressBar::new(library.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} entries",
            )
            .unwrap()
            .progress_chars("#>-"),
    );

    // Inner function: score one library entry using fragment co-elution
    // Returns Option<CalibrationMatch>
    let score_entry = |entry: &LibraryEntry,
                       candidate_spectra: &[&Spectrum],
                       candidate_global_indices: &[usize]|
     -> Option<CalibrationMatch> {
        if candidate_spectra.len() < MIN_COELUTION_SPECTRA {
            return None;
        }

        // Extract XICs for top 6 fragments
        let xics =
            super::extract_fragment_xics(&entry.fragments, candidate_spectra, tol_da, tol_ppm, 6);

        if xics.len() < 2 {
            return None; // Need at least 2 fragment XICs for co-elution
        }

        // Find reference XIC: fragment with highest total intensity
        let ref_idx = xics
            .iter()
            .enumerate()
            .max_by(|a, b| {
                let sum_a: f64 = a.1 .1.iter().map(|(_, v)| *v).sum();
                let sum_b: f64 = b.1 .1.iter().map(|(_, v)| *v).sum();
                sum_a.total_cmp(&sum_b)
            })
            .map(|(i, _)| i)?;

        let ref_xic = &xics[ref_idx].1;

        // Pre-filter: count fragments with non-trivial signal
        // (xics only contains fragments with ≥1 non-zero intensity point)
        let n_with_signal = xics
            .iter()
            .filter(|(_, xic)| xic.iter().any(|(_, v)| *v > 0.01))
            .count();
        if n_with_signal < 2 {
            return None; // Need at least 2 fragments with signal
        }

        // Peak-detection-first co-elution scoring:
        // 1. Detect candidate peaks in the reference XIC
        // 2. Score each peak by pairwise fragment correlation within its boundaries
        // 3. Pick the peak with the highest co-elution score
        let candidates = osprey_chromatography::detect_all_xic_peaks(ref_xic, 0.01, 5.0);

        if candidates.is_empty() {
            return None;
        }

        let min_corr_score = 0.5;

        // Score each candidate peak by pairwise fragment correlation
        let best = candidates
            .iter()
            .map(|bp| {
                let si = bp.start_index;
                let ei = bp.end_index;
                let peak_len = ei - si + 1;

                let corr_sum = if peak_len >= 3 {
                    let vals: Vec<Vec<f64>> = xics
                        .iter()
                        .map(|(_, xic)| xic[si..=ei].iter().map(|(_, v)| *v).collect())
                        .collect();
                    let mut sum = 0.0f64;
                    for ii in 0..vals.len() {
                        for jj in (ii + 1)..vals.len() {
                            sum += super::pearson_correlation_raw(&vals[ii], &vals[jj]);
                        }
                    }
                    sum
                } else {
                    0.0
                };

                (bp, corr_sum)
            })
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();

        let coelution_sum = best.1;
        if coelution_sum < min_corr_score {
            return None;
        }

        let ref_start = best.0.start_index;
        let ref_end = best.0.end_index;

        // Apex is the highest-intensity point within the peak boundaries
        let (apex_idx, _apex_val) = ref_xic[ref_start..=ref_end]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
            .map(|(i, &(_, v))| (ref_start + i, v))
            .unwrap_or((best.0.apex_index, 0.0));
        let apex_rt = ref_xic[apex_idx].0;
        let raw_ints: Vec<f64> = ref_xic.iter().map(|(_, v)| *v).collect();
        let signal_to_noise =
            osprey_chromatography::compute_snr(&raw_ints, apex_idx, ref_start, ref_end);

        // Peak width from boundaries (in minutes)
        let peak_width_minutes = ref_xic[ref_end].0 - ref_xic[ref_start].0;

        // Find the spectrum closest to apex for mass errors and scan number
        let apex_spec_local_idx = candidate_spectra
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1.retention_time - apex_rt)
                    .abs()
                    .partial_cmp(&(b.1.retention_time - apex_rt).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)?;

        let apex_spec = candidate_spectra[apex_spec_local_idx];
        let apex_global_idx = candidate_global_indices[apex_spec_local_idx];

        // Collect MS2 mass errors from top-N fragments at apex spectrum
        let (_has_match, ms2_errors) = super::topn_fragment_match_with_errors(
            &entry.fragments,
            &apex_spec.mzs,
            fragment_tolerance.tolerance,
            fragment_tolerance.unit,
        );
        let n_matched = ms2_errors.len();

        let avg_ms2_error = if !ms2_errors.is_empty() {
            Some(ms2_errors.iter().sum::<f64>() / ms2_errors.len() as f64)
        } else {
            None
        };

        // MS1 extraction
        let (ms1_error, observed_precursor_mz, isotope_cosine_score) = if let Some(ms1) = ms1_index
        {
            if let Some(ms1_spec) = ms1.find_nearest(apex_rt) {
                let envelope = IsotopeEnvelope::extract(
                    ms1_spec,
                    entry.precursor_mz,
                    entry.charge,
                    precursor_tolerance_ppm,
                );
                if envelope.has_m0() {
                    let isotope_score =
                        peptide_isotope_cosine(&entry.sequence, &envelope.intensities);
                    let error = envelope
                        .m0_observed_mz
                        .map(|obs_mz| fragment_tolerance.mass_error(entry.precursor_mz, obs_mz));
                    (error, envelope.m0_observed_mz, isotope_score)
                } else {
                    (None, None, None)
                }
            } else {
                let iso_center = spectra[apex_global_idx].isolation_window.center;
                let error = fragment_tolerance.mass_error(entry.precursor_mz, iso_center);
                (Some(error), Some(iso_center), None)
            }
        } else {
            let iso_center = spectra[apex_global_idx].isolation_window.center;
            let error = fragment_tolerance.mass_error(entry.precursor_mz, iso_center);
            (Some(error), Some(iso_center), None)
        };

        // === NEW: Calculate scoring features at apex ===

        // 1. LibCosine score at apex
        let libcosine_scorer = LibCosineScorer::with_tolerance(fragment_tolerance);
        let match_result = libcosine_scorer.match_fragments(entry, apex_spec);
        let libcosine_apex = libcosine_scorer.calculate_score(&match_result);

        // 2. Count top-6 matched ions at apex
        let top6_matched_apex =
            count_top6_matched_at_apex(&entry.fragments, &apex_spec.mzs, fragment_tolerance);

        // 3. Hyperscore at apex
        // Convert intensities from f32 to f64 for hyperscore calculation
        let apex_intensities_f64: Vec<f64> =
            apex_spec.intensities.iter().map(|&x| x as f64).collect();
        let hyperscore_result = super::compute_hyperscore(
            &entry.fragments,
            &apex_spec.mzs,
            &apex_intensities_f64,
            fragment_tolerance.tolerance,
            fragment_tolerance.unit,
        );

        // 4. XCorr at apex (using appropriate HRAM or unit resolution binning)
        let xcorr_at_apex = if let Some(scorer) = xcorr_scorer {
            scorer.xcorr(apex_spec, entry).xcorr
        } else {
            0.0
        };

        Some(CalibrationMatch {
            entry_id: entry.id,
            is_decoy: entry.is_decoy,
            library_rt: entry.retention_time,
            measured_rt: apex_rt,
            score: coelution_sum,
            ms1_error,
            library_precursor_mz: entry.precursor_mz,
            observed_precursor_mz,
            ms2_mass_errors: ms2_errors,
            avg_ms2_error,
            n_matched_fragments: n_matched,
            n_library_fragments: entry.fragments.len(),
            xcorr_score: xcorr_at_apex,
            evalue: 0.0,
            isotope_cosine_score,
            sequence: entry.sequence.clone(),
            charge: entry.charge,
            scan_number: spectra[apex_global_idx].scan_number,
            hyperscore: 0.0,
            n_b_ions: 0,
            n_y_ions: 0,

            // ML scoring features
            correlation_score: coelution_sum,
            libcosine_apex,
            top6_matched_apex,
            hyperscore_apex: hyperscore_result.score,
            signal_to_noise,
            peak_width_minutes: Some(peak_width_minutes),
            discriminant_score: coelution_sum, // Initially use correlation, will be replaced by LDA
            posterior_error: 0.0,              // Will be filled by LDA
            q_value: 1.0,                      // Will be filled by LDA
        })
    };

    // Process each window in parallel
    let all_matches: Vec<CalibrationMatch> = window_groups
        .par_iter()
        .flat_map(|((lower, upper), spectrum_indices)| {
            let window_entries: Vec<&LibraryEntry> = library
                .iter()
                .filter(|e| e.precursor_mz >= *lower && e.precursor_mz <= *upper)
                .collect();

            if window_entries.is_empty() {
                return Vec::new();
            }

            let window_spectra: Vec<&Spectrum> =
                spectrum_indices.iter().map(|&idx| &spectra[idx]).collect();

            if window_spectra.is_empty() {
                return Vec::new();
            }

            let mut window_matches: Vec<CalibrationMatch> = Vec::new();

            for entry in &window_entries {
                // Filter candidates: RT tolerance + top 2-of-6 fragment pre-filter
                let candidate_pairs: Vec<(usize, &Spectrum)> = window_spectra
                    .iter()
                    .enumerate()
                    .filter(|(_, spec)| {
                        let expected_rt = match expected_rt_fn {
                            Some(f) => f(entry.retention_time),
                            None => entry.retention_time,
                        };
                        let rt_diff = (spec.retention_time - expected_rt).abs();
                        if rt_diff > rt_tolerance {
                            return false;
                        }
                        has_topn_fragment_match(
                            &entry.fragments,
                            &spec.mzs,
                            fragment_tolerance.tolerance,
                            fragment_tolerance.unit,
                        )
                    })
                    .map(|(idx, spec)| (idx, *spec))
                    .collect();

                if candidate_pairs.is_empty() {
                    pb.inc(1);
                    continue;
                }

                // Sort by RT for XIC extraction
                let mut sorted_pairs = candidate_pairs;
                sorted_pairs.sort_by(|a, b| {
                    a.1.retention_time
                        .partial_cmp(&b.1.retention_time)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });

                let candidate_spectra: Vec<&Spectrum> =
                    sorted_pairs.iter().map(|(_, s)| *s).collect();
                let candidate_global_indices: Vec<usize> = sorted_pairs
                    .iter()
                    .map(|(local_idx, _)| spectrum_indices[*local_idx])
                    .collect();

                if let Some(m) = score_entry(entry, &candidate_spectra, &candidate_global_indices) {
                    window_matches.push(m);
                }
                pb.inc(1);
            }

            window_matches
        })
        .collect();

    pb.finish_with_message("Done");

    // Deduplicate (keep best co-elution score per entry)
    let mut best_matches: HashMap<u32, CalibrationMatch> = HashMap::new();
    for m in all_matches {
        best_matches
            .entry(m.entry_id)
            .and_modify(|existing| {
                if m.score > existing.score {
                    *existing = m.clone();
                }
            })
            .or_insert(m);
    }

    // Sort by score descending
    let mut results: Vec<CalibrationMatch> = best_matches.into_values().collect();
    results.sort_by(|a, b| b.score.total_cmp(&a.score));

    // Log scoring statistics
    let targets = results.iter().filter(|m| !m.is_decoy).count();
    let decoys = results.iter().filter(|m| m.is_decoy).count();

    log::info!(
        "Scoring complete: {} matches ({} targets, {} decoys)",
        results.len(),
        targets,
        decoys
    );

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{FragmentAnnotation, IsolationWindow, LibraryFragment};

    fn make_test_library_entry(
        id: u32,
        mz: f64,
        rt: f64,
        fragments: Vec<(f64, f32)>,
    ) -> LibraryEntry {
        let mut entry = LibraryEntry::new(id, "PEPTIDE".into(), "PEPTIDE".into(), 2, mz, rt);
        entry.fragments = fragments
            .into_iter()
            .map(|(mz, intensity)| LibraryFragment {
                mz,
                relative_intensity: intensity,
                annotation: FragmentAnnotation::default(),
            })
            .collect();
        entry
    }

    fn make_test_spectrum(idx: u32, rt: f64, peaks: Vec<(f64, f32)>) -> Spectrum {
        let (mzs, intensities): (Vec<f64>, Vec<f32>) = peaks.into_iter().unzip();
        Spectrum {
            scan_number: idx,
            retention_time: rt,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs,
            intensities,
        }
    }

    /// Verifies that PreprocessedLibrary correctly indexes library entries and maps entry IDs to row positions.
    #[test]
    fn test_preprocess_library() {
        let entries = vec![
            make_test_library_entry(1, 500.0, 10.0, vec![(300.0, 100.0), (400.0, 50.0)]),
            make_test_library_entry(2, 600.0, 15.0, vec![(350.0, 80.0), (450.0, 60.0)]),
        ];

        let preprocessed = PreprocessedLibrary::from_entries(&entries);

        assert_eq!(preprocessed.len(), 2);
        assert_eq!(preprocessed.entry_ids, vec![1, 2]);
        assert!(preprocessed.id_to_row.contains_key(&1));
        assert!(preprocessed.id_to_row.contains_key(&2));
    }

    /// Verifies that PreprocessedSpectra correctly stores spectrum indices and retention times.
    #[test]
    fn test_preprocess_spectra() {
        let spectra = vec![
            make_test_spectrum(1, 10.0, vec![(300.0, 100.0), (400.0, 50.0)]),
            make_test_spectrum(2, 15.0, vec![(350.0, 80.0), (450.0, 60.0)]),
        ];

        let preprocessed = PreprocessedSpectra::from_spectra(&spectra);

        assert_eq!(preprocessed.len(), 2);
        assert_eq!(preprocessed.spectrum_indices, vec![0, 1]);
        assert_eq!(preprocessed.retention_times, vec![10.0, 15.0]);
    }

    /// Verifies that batch XCorr scoring returns a score near 1.0 when library and observed spectra have identical peaks.
    #[test]
    fn test_batch_scoring_perfect_match() {
        // Create library entry and spectrum with identical peaks
        let entries = vec![make_test_library_entry(
            1,
            500.0,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)],
        )];

        let spectra = vec![make_test_spectrum(
            1,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)],
        )];

        let scorer = BatchScorer::new();
        let lib = scorer.preprocess_library(&entries);
        let spec = scorer.preprocess_spectra(&spectra);

        let scores = scorer.score_all(&lib, &spec);

        // Perfect match should give score close to 1.0
        assert!(
            scores[[0, 0]] > 0.99,
            "Expected ~1.0, got {}",
            scores[[0, 0]]
        );
    }

    /// Verifies that batch XCorr scoring returns a score near 0.0 when library and observed spectra have no overlapping peaks.
    #[test]
    fn test_batch_scoring_no_match() {
        // Create library entry and spectrum with completely different peaks
        let entries = vec![make_test_library_entry(
            1,
            500.0,
            10.0,
            vec![(300.0, 100.0)],
        )];

        let spectra = vec![make_test_spectrum(
            1,
            10.0,
            vec![(800.0, 100.0)], // Completely different m/z
        )];

        let scorer = BatchScorer::new();
        let lib = scorer.preprocess_library(&entries);
        let spec = scorer.preprocess_spectra(&spectra);

        let scores = scorer.score_all(&lib, &spec);

        // No overlap should give score close to 0.0
        assert!(
            scores[[0, 0]] < 0.01,
            "Expected ~0.0, got {}",
            scores[[0, 0]]
        );
    }

    /// Verifies that find_best_matches correctly pairs each library entry with its highest-scoring observed spectrum.
    #[test]
    fn test_find_best_matches() {
        let entries = vec![
            make_test_library_entry(1, 500.0, 10.0, vec![(300.0, 100.0)]),
            make_test_library_entry(2, 600.0, 15.0, vec![(400.0, 100.0)]),
        ];

        let spectra = vec![
            make_test_spectrum(1, 10.0, vec![(300.0, 100.0)]), // Matches entry 1
            make_test_spectrum(2, 15.0, vec![(400.0, 100.0)]), // Matches entry 2
        ];

        let scorer = BatchScorer::new();
        let lib = scorer.preprocess_library(&entries);
        let spec = scorer.preprocess_spectra(&spectra);

        let matches = scorer.find_best_matches(&lib, &spec);

        assert_eq!(matches.len(), 2);

        // Entry 1 should match spectrum 0 (index in original spectra)
        let entry1_match = matches.iter().find(|(id, _, _, _)| *id == 1).unwrap();
        assert_eq!(entry1_match.1, 0); // spectrum index
        assert!(entry1_match.2 > 0.99); // score

        // Entry 2 should match spectrum 1
        let entry2_match = matches.iter().find(|(id, _, _, _)| *id == 2).unwrap();
        assert_eq!(entry2_match.1, 1);
        assert!(entry2_match.2 > 0.99);
    }

    /// Verifies that subsetting a PreprocessedLibrary by entry IDs retains only the requested entries.
    #[test]
    fn test_library_subset() {
        let entries = vec![
            make_test_library_entry(1, 500.0, 10.0, vec![(300.0, 100.0)]),
            make_test_library_entry(2, 600.0, 15.0, vec![(400.0, 100.0)]),
            make_test_library_entry(3, 700.0, 20.0, vec![(500.0, 100.0)]),
        ];

        let lib = PreprocessedLibrary::from_entries(&entries);
        let subset = lib.subset(&[1, 3]);

        assert_eq!(subset.len(), 2);
        assert!(subset.id_to_row.contains_key(&1));
        assert!(subset.id_to_row.contains_key(&3));
        assert!(!subset.id_to_row.contains_key(&2));
    }

    // =========================================================================
    // LibCosine Scorer Tests (ppm-based matching, NO binning)
    // =========================================================================

    /// Verifies that LibCosineScorer returns a score near 1.0 for an exact peak match using PPM-based tolerance.
    #[test]
    fn test_libcosine_perfect_match() {
        // Library entry with fragments
        let entry = make_test_library_entry(
            1,
            500.0,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)],
        );

        // Spectrum with EXACT same peaks
        let spectrum =
            make_test_spectrum(1, 10.0, vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)]);

        let scorer = LibCosineScorer::hram(10.0); // 10 ppm
        let score = scorer.score(&entry, &spectrum);

        // Perfect match should give score = 1.0
        assert!(score > 0.99, "Expected ~1.0, got {}", score);
    }

    /// Verifies that LibCosineScorer returns a score near 0.0 when no observed peaks fall within PPM tolerance of library fragments.
    #[test]
    fn test_libcosine_no_match() {
        // Library entry with fragments at low m/z
        let entry = make_test_library_entry(1, 500.0, 10.0, vec![(300.0, 100.0)]);

        // Spectrum with peaks at completely different m/z
        let spectrum = make_test_spectrum(1, 10.0, vec![(800.0, 100.0)]);

        let scorer = LibCosineScorer::hram(10.0);
        let score = scorer.score(&entry, &spectrum);

        // No overlap should give score = 0.0
        assert!(score < 0.01, "Expected ~0.0, got {}", score);
    }

    /// Verifies that LibCosineScorer matches peaks within the specified PPM tolerance and rejects those outside it.
    #[test]
    fn test_libcosine_ppm_matching() {
        // Library fragment at 500.0 m/z
        let entry = make_test_library_entry(1, 500.0, 10.0, vec![(500.0, 100.0)]);

        // Spectrum with peak at 500.005 m/z (10 ppm offset from 500.0)
        // 10 ppm of 500.0 = 0.005 m/z
        let spectrum = make_test_spectrum(1, 10.0, vec![(500.005, 100.0)]);

        let scorer = LibCosineScorer::hram(10.0);
        let score = scorer.score(&entry, &spectrum);

        // Should match within 10 ppm tolerance
        assert!(score > 0.99, "Expected match at 10 ppm, got {}", score);

        // With tighter tolerance (5 ppm), should NOT match
        let scorer_tight = LibCosineScorer::hram(5.0);
        let score_tight = scorer_tight.score(&entry, &spectrum);
        assert!(
            score_tight < 0.01,
            "Expected no match at 5 ppm, got {}",
            score_tight
        );
    }

    /// Verifies that score_with_errors returns accurate per-fragment mass errors in PPM units.
    #[test]
    fn test_libcosine_mass_errors() {
        // Library fragment at 500.0 m/z
        let entry = make_test_library_entry(1, 500.0, 10.0, vec![(500.0, 100.0)]);

        // Spectrum with peak at 500.005 m/z (10 ppm offset)
        let spectrum = make_test_spectrum(1, 10.0, vec![(500.005, 100.0)]);

        let scorer = LibCosineScorer::hram(20.0);
        let (score, mass_errors) = scorer.score_with_errors(&entry, &spectrum);

        assert!(score > 0.5);
        assert_eq!(mass_errors.len(), 1);

        // Expected mass error: (500.005 - 500.0) / 500.0 * 1e6 = 10 ppm
        let expected_ppm = 10.0;
        assert!(
            (mass_errors[0] - expected_ppm).abs() < 0.1,
            "Expected ~10 ppm, got {} ppm",
            mass_errors[0]
        );
    }

    /// Verifies that unit-resolution LibCosineScorer matches peaks within Dalton tolerance and rejects those outside it.
    #[test]
    fn test_libcosine_da_tolerance() {
        // Library fragment at 500.0 m/z
        let entry = make_test_library_entry(1, 500.0, 10.0, vec![(500.0, 100.0)]);

        // Spectrum with peak at 500.25 m/z (0.25 Da offset)
        let spectrum = make_test_spectrum(1, 10.0, vec![(500.25, 100.0)]);

        // 0.3 Da tolerance should match
        let scorer_loose = LibCosineScorer::unit_resolution(0.3);
        let score_loose = scorer_loose.score(&entry, &spectrum);
        assert!(
            score_loose > 0.99,
            "Expected match at 0.3 Da, got {}",
            score_loose
        );

        // 0.2 Da tolerance should NOT match
        let scorer_tight = LibCosineScorer::unit_resolution(0.2);
        let score_tight = scorer_tight.score(&entry, &spectrum);
        assert!(
            score_tight < 0.01,
            "Expected no match at 0.2 Da, got {}",
            score_tight
        );
    }

    /// Verifies that match_fragments correctly identifies matched and unmatched library fragments and reports per-fragment mass errors.
    #[test]
    fn test_libcosine_multiple_fragments() {
        // Library with 3 fragments
        let entry = make_test_library_entry(
            1,
            500.0,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)],
        );

        // Spectrum matching only 2 of 3 fragments
        let spectrum = make_test_spectrum(1, 10.0, vec![(300.0, 100.0), (400.0, 50.0)]);

        let scorer = LibCosineScorer::hram(10.0);
        let match_result = scorer.match_fragments(&entry, &spectrum);

        // Should have 3 library fragments, 2 matched
        assert_eq!(match_result.matched_lib_intensities.len(), 3);
        assert_eq!(match_result.n_matched, 2);
        assert_eq!(match_result.mass_errors.len(), 2);

        // Third fragment should have 0 experimental intensity
        assert!(match_result.matched_exp_intensities[2] < 1e-10);
    }

    /// Verifies that fragment matching selects the closest m/z peak rather than the highest intensity peak when multiple candidates fall within tolerance.
    #[test]
    fn test_libcosine_closest_mz_selected() {
        // Library fragment at 500.0 m/z
        let entry = make_test_library_entry(1, 500.0, 10.0, vec![(500.0, 100.0)]);

        // Spectrum with multiple peaks within tolerance, different intensities
        // Should select the CLOSEST m/z, not highest intensity (pyXcorrDIA approach)
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![499.996, 500.000, 500.004], // All within 10 ppm
            intensities: vec![200.0, 50.0, 150.0], // 200 is highest, but 50 is closest
        };

        let scorer = LibCosineScorer::hram(10.0);
        let match_result = scorer.match_fragments(&entry, &spectrum);

        // Should select the peak at 500.000 (intensity 50) because it's CLOSEST to lib_mz
        // NOT the peak at 499.996 (intensity 200) which has higher intensity
        assert_eq!(match_result.n_matched, 1);
        assert!(
            (match_result.matched_exp_intensities[0] - 50.0).abs() < 1e-10,
            "Expected intensity 50.0 (closest m/z), got {}",
            match_result.matched_exp_intensities[0]
        );

        // Mass error should be ~0 ppm since we matched 500.000 to 500.0
        assert!(
            match_result.mass_errors[0].abs() < 0.1,
            "Expected ~0 ppm, got {} ppm",
            match_result.mass_errors[0]
        );
    }

    #[test]
    fn test_count_top6_matched_at_apex() {
        use osprey_core::ToleranceUnit;

        // Create test library fragments with varying intensities
        let fragments = vec![
            LibraryFragment {
                mz: 100.0,
                relative_intensity: 0.5,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 200.0,
                relative_intensity: 1.0, // Top 1
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 0.2, // Not in top 6
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 0.8, // Top 2
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 0.7, // Top 3
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 600.0,
                relative_intensity: 0.6, // Top 4
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 700.0,
                relative_intensity: 0.55, // Top 5
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 800.0,
                relative_intensity: 0.51, // Top 6
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Spectrum with peaks matching some of the top-6 fragments
        // Matches: 200.0, 400.0, 500.0, 600.0 (4 out of top 6)
        let spectrum_mzs = vec![200.0, 350.0, 400.0, 500.0, 600.0];

        let tolerance_config = FragmentToleranceConfig {
            tolerance: 10.0,
            unit: ToleranceUnit::Ppm,
        };

        let count = count_top6_matched_at_apex(&fragments, &spectrum_mzs, tolerance_config);
        assert_eq!(count, 4); // 200, 400, 500, 600 matched

        // Test with no matches
        let empty_spectrum: Vec<f64> = vec![];
        let count = count_top6_matched_at_apex(&fragments, &empty_spectrum, tolerance_config);
        assert_eq!(count, 0);

        // Test with all top-6 matched
        let all_matched_spectrum = vec![200.0, 400.0, 500.0, 600.0, 700.0, 800.0];
        let count = count_top6_matched_at_apex(&fragments, &all_matched_spectrum, tolerance_config);
        assert_eq!(count, 6);

        // Test with fewer than 6 fragments total
        let few_fragments = vec![
            LibraryFragment {
                mz: 100.0,
                relative_intensity: 1.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 200.0,
                relative_intensity: 0.8,
                annotation: FragmentAnnotation::default(),
            },
        ];
        let few_spectrum = vec![100.0, 200.0];
        let count = count_top6_matched_at_apex(&few_fragments, &few_spectrum, tolerance_config);
        assert_eq!(count, 2); // All 2 fragments matched
    }
}
