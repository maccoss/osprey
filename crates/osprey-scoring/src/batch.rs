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
//! - SMZ preprocessing: `sqrt(intensity) × mz²`
//! - L2 normalization and cosine angle calculation
//! - Stores mass errors for matched fragments (for calibration)
//!
//! ## BatchScorer (XCorr/Ridge Regression)
//!
//! The `BatchScorer` uses binning for BLAS-accelerated XCorr calculations:
//! - Unit resolution: 1.0005079 Da bins, 0.4 offset
//! - HRAM: 0.02 Da bins, 0 offset
//!
//! Performance: 10-20× faster than one-at-a-time scoring for typical DIA data.

// Ensure BLAS is linked
extern crate blas_src;
extern crate openblas_src;

use ndarray::{Array1, Array2};
use osprey_core::{FragmentToleranceConfig, IsotopeEnvelope, LibraryEntry, MS1Spectrum, Spectrum, ToleranceUnit, peptide_isotope_cosine};
use rayon::prelude::*;
use std::collections::HashMap;

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
    /// Mass errors for matched fragments (in configured unit: ppm or Da)
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
/// # SMZ Preprocessing
///
/// Both library and experimental intensities are preprocessed:
/// ```text
/// preprocessed_intensity = sqrt(intensity) × mz²
/// ```
///
/// # Scoring
///
/// The score is the cosine angle between preprocessed vectors:
/// ```text
/// score = (exp · lib) / (|exp| × |lib|)
/// ```
#[derive(Debug, Clone)]
pub struct LibCosineScorer {
    /// Fragment tolerance configuration
    pub tolerance: FragmentToleranceConfig,
}

impl Default for LibCosineScorer {
    fn default() -> Self {
        Self {
            tolerance: FragmentToleranceConfig::default(), // 10 ppm
        }
    }
}

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

            // Calculate mass error if matched
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

        // Apply SMZ preprocessing: sqrt(intensity) × mz²
        let exp_preprocessed: Vec<f64> = match_result
            .matched_exp_intensities
            .iter()
            .zip(match_result.matched_lib_mzs.iter())
            .map(|(&int, &mz)| int.sqrt() * mz.powi(2))
            .collect();

        let lib_preprocessed: Vec<f64> = match_result
            .matched_lib_intensities
            .iter()
            .zip(match_result.matched_lib_mzs.iter())
            .map(|(&int, &mz)| int.sqrt() * mz.powi(2))
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
        cos_angle.max(0.0).min(1.0)
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
#[derive(Debug, Clone)]
pub struct PreprocessedLibrary {
    /// Matrix of preprocessed library vectors (n_entries × n_bins)
    /// Each row is an SMZ-preprocessed, L2-normalized spectrum
    pub matrix: Array2<f64>,
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
        let preprocessed: Vec<(u32, Array1<f64>)> = entries
            .par_iter()
            .filter(|e| !e.fragments.is_empty())
            .map(|entry| {
                let binned = Self::preprocess_library_entry(entry, num_bins, min_mz, bin_width);
                (entry.id, binned)
            })
            .collect();

        // Build matrix and mappings
        let n_valid = preprocessed.len();
        let mut matrix = Array2::zeros((n_valid, num_bins));
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
    fn preprocess_library_entry(
        entry: &LibraryEntry,
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
    ) -> Array1<f64> {
        let mut binned = Array1::zeros(num_bins);

        // Apply SMZ preprocessing: sqrt(intensity) × m/z²
        for frag in &entry.fragments {
            if frag.mz >= min_mz && frag.mz < min_mz + num_bins as f64 * bin_width {
                let bin = ((frag.mz - min_mz) / bin_width) as usize;
                if bin < num_bins {
                    let smz = (frag.relative_intensity as f64).sqrt() * frag.mz.powi(2);
                    binned[bin] += smz;
                }
            }
        }

        // L2 normalize
        let norm: f64 = binned.mapv(|x: f64| x * x).sum();
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
        let mut matrix = Array2::zeros((n_rows, self.num_bins));
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
#[derive(Debug, Clone)]
pub struct PreprocessedSpectra {
    /// Matrix of preprocessed spectra (n_spectra × n_bins)
    /// Each row is an SMZ-preprocessed, L2-normalized spectrum
    pub matrix: Array2<f64>,
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
        let preprocessed: Vec<(usize, f64, Array1<f64>)> = spectra
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
        let mut matrix = Array2::zeros((n_valid, num_bins));
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
    fn preprocess_spectrum(
        spectrum: &Spectrum,
        num_bins: usize,
        min_mz: f64,
        bin_width: f64,
    ) -> Array1<f64> {
        let mut binned = Array1::zeros(num_bins);

        // Apply SMZ preprocessing: sqrt(intensity) × m/z²
        for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            if mz >= min_mz && mz < min_mz + num_bins as f64 * bin_width {
                let bin = ((mz - min_mz) / bin_width) as usize;
                if bin < num_bins {
                    let smz = (intensity as f64).sqrt() * mz.powi(2);
                    binned[bin] += smz;
                }
            }
        }

        // L2 normalize
        let norm: f64 = binned.mapv(|x: f64| x * x).sum();
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
        let mut matrix = Array2::zeros((n_rows, self.num_bins));
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
    pub fn preprocess_library_with_context(&self, entries: &[LibraryEntry], context: Option<&str>) -> PreprocessedLibrary {
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
    pub fn preprocess_spectra_with_context(&self, spectra: &[Spectrum], context: Option<&str>) -> PreprocessedSpectra {
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
    pub fn score_all(
        &self,
        library: &PreprocessedLibrary,
        spectra: &PreprocessedSpectra,
    ) -> Array2<f64> {
        // Matrix multiplication: (n_library × n_bins) · (n_bins × n_spectra)
        // Result: (n_library × n_spectra)
        library.matrix.dot(&spectra.matrix.t())
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
                let (best_col, &best_score) = row
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))?;

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

        // Score against all spectra
        let scores = spectra.matrix.dot(&entry_vec);

        let mut results: Vec<(usize, f64, f64)> = scores
            .iter()
            .enumerate()
            .filter(|(_, &score)| score > 0.0)
            .map(|(row, &score)| {
                let spec_idx = spectra.spectrum_indices[row];
                let rt = spectra.retention_times[row];
                (spec_idx, score, rt)
            })
            .collect();

        // Sort by score descending
        results.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
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
    /// LibCosine score
    pub score: f64,
    /// MS1 (precursor) m/z error in PPM (observed - theoretical)
    pub ms1_ppm_error: Option<f64>,
    /// Library precursor m/z (theoretical)
    pub library_precursor_mz: f64,
    /// Observed precursor m/z (from spectrum isolation window center)
    pub observed_precursor_mz: Option<f64>,
    /// MS2 fragment mass errors (in configured unit: ppm or Da)
    pub ms2_mass_errors: Vec<f64>,
    /// Average MS2 mass error (only for matched fragments)
    pub avg_ms2_error: Option<f64>,
    /// Number of matched fragments
    pub n_matched_fragments: usize,
    /// Total number of library fragments
    pub n_library_fragments: usize,
    /// XCorr score (Comet-style cross-correlation)
    pub xcorr_score: f64,
    /// Isotope cosine score (MS1 isotope envelope match)
    pub isotope_cosine_score: Option<f64>,
    /// Peptide sequence
    pub sequence: String,
    /// Charge state
    pub charge: u8,
    /// Best matching spectrum scan number
    pub scan_number: u32,
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
    /// Target LibCosine score
    pub target_libcosine: f64,
    /// Decoy LibCosine score
    pub decoy_libcosine: f64,
    /// Winning LibCosine score (max of target and decoy) - this is what FDR uses
    pub winning_libcosine: f64,
    /// Target XCorr score
    pub target_xcorr: f64,
    /// Decoy XCorr score
    pub decoy_xcorr: f64,
    /// Winning XCorr score (max of target and decoy)
    pub winning_xcorr: f64,
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
    /// True if target wins the competition (target_libcosine > decoy_libcosine)
    pub target_wins: bool,
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
            let winning_libcosine = target_match.score.max(decoy_match.score);
            let winning_xcorr = target_match.xcorr_score.max(decoy_match.xcorr_score);

            paired.push(PairedCalibrationResult {
                target_entry_id: target_id,
                charge: target_match.charge,
                target_sequence: target_match.sequence.clone(),
                decoy_sequence: decoy_match.sequence.clone(),
                target_libcosine: target_match.score,
                decoy_libcosine: decoy_match.score,
                winning_libcosine,
                target_xcorr: target_match.xcorr_score,
                decoy_xcorr: decoy_match.xcorr_score,
                winning_xcorr,
                target_isotope_score: target_match.isotope_cosine_score,
                decoy_isotope_score: decoy_match.isotope_cosine_score,
                target_precursor_error_ppm: target_match.ms1_ppm_error,
                decoy_precursor_error_ppm: decoy_match.ms1_ppm_error,
                target_rt: target_match.measured_rt,
                decoy_rt: decoy_match.measured_rt,
                library_rt: target_match.library_rt,
                expected_rt,
                target_delta_rt,
                decoy_delta_rt,
                target_matched_frags: target_match.n_matched_fragments,
                decoy_matched_frags: decoy_match.n_matched_fragments,
                target_wins,
            });
        }
    }

    // Sort by WINNING LibCosine score descending - this matches FDR calculation order
    paired.sort_by(|a, b| b.winning_libcosine.partial_cmp(&a.winning_libcosine).unwrap_or(std::cmp::Ordering::Equal));

    paired
}

/// Sample peptides from the library for calibration (pyXcorrDIA strategy)
///
/// For large libraries (millions of entries), scoring all entries for calibration
/// is too slow. Instead, we sample a representative subset:
///
/// - Default: 2000 peptides (can be doubled if first attempt fails)
/// - Samples diverse precursor m/z range to ensure good calibration coverage
/// - Returns both sampled targets AND their paired decoys
///
/// # Arguments
/// * `library` - Full library (targets + decoys)
/// * `sample_size` - Number of target peptides to sample (0 = use all)
/// * `seed` - Random seed for reproducibility
///
/// # Returns
/// Sampled subset of library entries (targets + their decoys)
pub fn sample_library_for_calibration(
    library: &[LibraryEntry],
    sample_size: usize,
    seed: u64,
) -> Vec<LibraryEntry> {
    use std::collections::HashSet;

    if sample_size == 0 {
        // Use all entries
        return library.to_vec();
    }

    // Separate targets and decoys
    let targets: Vec<&LibraryEntry> = library.iter().filter(|e| !e.is_decoy).collect();
    let decoys: Vec<&LibraryEntry> = library.iter().filter(|e| e.is_decoy).collect();

    if targets.len() <= sample_size {
        // Library is small enough, use everything
        log::debug!("Library has {} targets, using all (sample_size={})", targets.len(), sample_size);
        return library.to_vec();
    }

    // Build a map from target_id to decoy (decoy_id = target_id | 0x80000000)
    let decoy_map: HashMap<u32, &LibraryEntry> = decoys
        .iter()
        .map(|d| (d.id & 0x7FFFFFFF, *d))
        .collect();

    // Simple deterministic sampling: stride through sorted targets
    // This ensures we get good coverage across m/z range
    let mut sorted_targets = targets.clone();
    sorted_targets.sort_by(|a, b| {
        a.precursor_mz.partial_cmp(&b.precursor_mz).unwrap_or(std::cmp::Ordering::Equal)
    });

    // Use seed to create a starting offset
    let offset = (seed as usize) % sorted_targets.len().max(1);
    let stride = sorted_targets.len() / sample_size;
    let stride = stride.max(1);

    let mut sampled_target_ids: HashSet<u32> = HashSet::new();
    let mut sampled: Vec<LibraryEntry> = Vec::with_capacity(sample_size * 2);

    for i in 0..sample_size {
        let idx = (offset + i * stride) % sorted_targets.len();
        let target = sorted_targets[idx];

        if sampled_target_ids.contains(&target.id) {
            continue; // Already sampled (shouldn't happen with stride)
        }
        sampled_target_ids.insert(target.id);

        // Add target
        sampled.push(target.clone());

        // Add corresponding decoy if exists
        if let Some(decoy) = decoy_map.get(&target.id) {
            sampled.push((*decoy).clone());
        }
    }

    log::info!(
        "Sampled {} targets + {} decoys from {} total entries for calibration",
        sampled_target_ids.len(),
        sampled.len() - sampled_target_ids.len(),
        library.len()
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

        windows.entry((lower_key, upper_key))
            .or_default()
            .push(idx);
    }

    // Convert back to f64 windows with spectrum indices
    windows
        .into_iter()
        .map(|((lower, upper), indices)| {
            ((lower as f64 / 10.0, upper as f64 / 10.0), indices)
        })
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
/// * `rt_tolerance` - RT tolerance for matching
///
/// # Returns
/// Calibration matches sorted by score descending
pub fn run_windowed_calibration_scoring(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
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
            let cloned_entries: Vec<LibraryEntry> = window_entries.iter().map(|e| (*e).clone()).collect();
            let window_context = format!("window {:.1}-{:.1} m/z", lower, upper);
            let preprocessed_library = scorer.preprocess_library_with_context(&cloned_entries, Some(&window_context));
            let preprocessed_spectra = scorer.preprocess_spectra_with_context(&window_spectra, Some(&window_context));

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

            // Convert to CalibrationMatch with MS1 PPM error
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

                    // Calculate MS1 PPM error if we have observed m/z
                    let ms1_ppm_error = observed_precursor_mz.map(|obs_mz| {
                        ((obs_mz - entry.precursor_mz) / entry.precursor_mz) * 1e6
                    });

                    // Get scan number from spectrum
                    let scan_number = if spec_idx < window_spectra.len() {
                        window_spectra[spec_idx].scan_number
                    } else {
                        0
                    };

                    Some(CalibrationMatch {
                        entry_id,
                        is_decoy: entry.is_decoy,
                        library_rt: entry.retention_time,
                        measured_rt,
                        score,
                        ms1_ppm_error,
                        library_precursor_mz: entry.precursor_mz,
                        observed_precursor_mz,
                        // Binned BatchScorer doesn't compute MS2 errors
                        ms2_mass_errors: Vec::new(),
                        avg_ms2_error: None,
                        n_matched_fragments: 0,
                        n_library_fragments: entry.fragments.len(),
                        // For binned BatchScorer, the score IS the XCorr-like score
                        xcorr_score: score,
                        isotope_cosine_score: None,
                        sequence: entry.sequence.clone(),
                        charge: entry.charge,
                        scan_number,
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
    results.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

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
    rt_tolerance: f64,
) -> Vec<CalibrationMatch> {
    let scorer = BatchScorer::new();

    // Preprocess library and spectra
    log::info!("Preprocessing {} library entries for batch scoring", library.len());
    let preprocessed_library = scorer.preprocess_library(library);

    log::info!("Preprocessing {} spectra for batch scoring", spectra.len());
    let preprocessed_spectra = scorer.preprocess_spectra(spectra);

    // Get library RTs in the same order as preprocessed rows
    let library_rts: Vec<f64> = preprocessed_library
        .entry_ids
        .iter()
        .filter_map(|id| library.iter().find(|e| e.id == *id).map(|e| e.retention_time))
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

            // Calculate MS1 PPM error if we have observed m/z
            let ms1_ppm_error = observed_precursor_mz.map(|obs_mz| {
                ((obs_mz - entry.precursor_mz) / entry.precursor_mz) * 1e6
            });

            // Get scan number from spectrum
            let scan_number = if spec_idx < spectra.len() {
                spectra[spec_idx].scan_number
            } else {
                0
            };

            Some(CalibrationMatch {
                entry_id,
                is_decoy: entry.is_decoy,
                library_rt: entry.retention_time,
                measured_rt,
                score,
                ms1_ppm_error,
                library_precursor_mz: entry.precursor_mz,
                observed_precursor_mz,
                // Binned BatchScorer doesn't compute MS2 errors
                ms2_mass_errors: Vec::new(),
                avg_ms2_error: None,
                n_matched_fragments: 0,
                n_library_fragments: entry.fragments.len(),
                // For binned BatchScorer, the score IS the XCorr-like score
                xcorr_score: score,
                isotope_cosine_score: None,
                sequence: entry.sequence.clone(),
                charge: entry.charge,
                scan_number,
            })
        })
        .collect();

    // Sort by score descending
    results.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

    log::info!("Batch scoring complete: {} matches found", results.len());

    results
}

/// Run LibCosine calibration scoring (pyXcorrDIA-compatible)
///
/// This is the proper implementation that matches pyXcorrDIA exactly:
/// - Uses ppm-based fragment matching (NO binning)
/// - Applies SMZ preprocessing: sqrt(intensity) × mz²
/// - Calculates cosine angle between preprocessed vectors
/// - Stores MS2 mass errors for matched fragments (for calibration)
///
/// # Arguments
/// * `library` - Full library (targets + decoys)
/// * `spectra` - All MS2 spectra
/// * `fragment_tolerance` - Fragment m/z tolerance (ppm or Da)
/// * `rt_tolerance` - RT tolerance for matching
///
/// # Returns
/// Calibration matches sorted by score descending
pub fn run_libcosine_calibration_scoring(
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
        "LibCosine scoring: {} spectra in {} isolation windows ({} ppm tolerance)",
        spectra.len(),
        window_groups.len(),
        if fragment_tolerance.unit == ToleranceUnit::Ppm {
            format!("{:.1}", fragment_tolerance.tolerance)
        } else {
            format!("{:.3} Da", fragment_tolerance.tolerance)
        }
    );

    // Create the LibCosine scorer with configured tolerance
    let scorer = LibCosineScorer::with_tolerance(fragment_tolerance);

    // Create XCorr scorer for batch preprocessing
    let xcorr_scorer = SpectralScorer::new();

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
            let window_spectra: Vec<&Spectrum> = spectrum_indices
                .iter()
                .map(|&idx| &spectra[idx])
                .collect();

            if window_spectra.is_empty() {
                return Vec::new();
            }

            // === BATCH XCORR PREPROCESSING (pyXcorrDIA approach) ===
            // Preprocess all spectra for XCorr ONCE per window
            let spec_xcorr_preprocessed: Vec<Vec<f64>> = window_spectra
                .iter()
                .map(|spec| xcorr_scorer.preprocess_spectrum_for_xcorr(spec))
                .collect();

            // Preprocess all library entries for XCorr ONCE per window
            let entry_xcorr_preprocessed: HashMap<u32, Vec<f64>> = window_entries
                .iter()
                .map(|entry| (entry.id, xcorr_scorer.preprocess_library_for_xcorr(entry)))
                .collect();

            // Score each library entry against all spectra in this window
            let mut window_matches: Vec<CalibrationMatch> = Vec::new();

            for entry in &window_entries {
                let mut best_score = 0.0f64;
                let mut best_local_idx = 0usize;  // Index into window_spectra
                let mut best_spec_idx = 0usize;   // Index into global spectra
                let mut best_measured_rt = 0.0f64;
                let mut best_ms2_errors = Vec::new();
                let mut best_n_matched = 0usize;

                // Find best matching spectrum within RT tolerance (by LibCosine)
                for (local_idx, &spec) in window_spectra.iter().enumerate() {
                    let rt_diff = (spec.retention_time - entry.retention_time).abs();
                    if rt_diff > rt_tolerance {
                        continue;
                    }

                    // Score using LibCosine with ppm matching (NO binning)
                    let (score, ms2_errors) = scorer.score_with_errors(entry, spec);

                    if score > best_score {
                        best_score = score;
                        best_local_idx = local_idx;
                        best_spec_idx = spectrum_indices[local_idx];
                        best_measured_rt = spec.retention_time;
                        best_ms2_errors = ms2_errors;
                        best_n_matched = best_ms2_errors.len();
                    }
                }

                if best_score > 0.0 {
                    // Get observed precursor m/z from best spectrum
                    let observed_precursor_mz = Some(spectra[best_spec_idx].isolation_window.center);

                    // Calculate MS1 PPM error
                    let ms1_ppm_error = observed_precursor_mz.map(|obs_mz| {
                        ((obs_mz - entry.precursor_mz) / entry.precursor_mz) * 1e6
                    });

                    // Calculate average MS2 error
                    let avg_ms2_error = if !best_ms2_errors.is_empty() {
                        Some(best_ms2_errors.iter().sum::<f64>() / best_ms2_errors.len() as f64)
                    } else {
                        None
                    };

                    // === LOOKUP XCORR FROM PREPROCESSED VECTORS ===
                    // XCorr was precomputed for ALL spectra, just look up at best LibCosine index
                    let xcorr_score = if let Some(entry_preprocessed) = entry_xcorr_preprocessed.get(&entry.id) {
                        SpectralScorer::xcorr_from_preprocessed(
                            &spec_xcorr_preprocessed[best_local_idx],
                            entry_preprocessed,
                        )
                    } else {
                        0.0
                    };

                    window_matches.push(CalibrationMatch {
                        entry_id: entry.id,
                        is_decoy: entry.is_decoy,
                        library_rt: entry.retention_time,
                        measured_rt: best_measured_rt,
                        score: best_score,
                        ms1_ppm_error,
                        library_precursor_mz: entry.precursor_mz,
                        observed_precursor_mz,
                        ms2_mass_errors: best_ms2_errors,
                        avg_ms2_error,
                        n_matched_fragments: best_n_matched,
                        n_library_fragments: entry.fragments.len(),
                        xcorr_score,
                        isotope_cosine_score: None,
                        sequence: entry.sequence.clone(),
                        charge: entry.charge,
                        scan_number: spectra[best_spec_idx].scan_number,
                    });
                }
            }

            window_matches
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
    results.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(std::cmp::Ordering::Equal));

    // Log scoring statistics
    let total_ms2_errors: usize = results.iter().map(|m| m.ms2_mass_errors.len()).sum();
    let matched_count = results.iter().filter(|m| m.n_matched_fragments > 0).count();
    log::info!(
        "LibCosine scoring complete: {} unique matches ({} with fragment matches, {} total MS2 errors)",
        results.len(),
        matched_count,
        total_ms2_errors
    );

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

/// Run LibCosine calibration scoring with MS1 isotope extraction (pyXcorrDIA-compatible)
///
/// This is the full pyXcorrDIA-compatible implementation that extracts the M+0
/// isotope peak from actual MS1 spectra for accurate MS1 mass calibration.
///
/// ## Difference from `run_libcosine_calibration_scoring`
///
/// - **With MS1**: Extracts M+0 peak from MS1 spectrum → accurate mass error
/// - **Without MS1**: Uses isolation window center → not accurate for calibration
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
pub fn run_libcosine_calibration_scoring_with_ms1<M: MS1SpectrumLookup>(
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
        "LibCosine scoring with MS1 extraction: {} spectra in {} windows (MS1 available: {})",
        spectra.len(),
        window_groups.len(),
        ms1_available
    );

    // Create the LibCosine scorer with configured tolerance
    let scorer = LibCosineScorer::with_tolerance(fragment_tolerance);

    // Create XCorr scorer for batch preprocessing
    let xcorr_scorer = SpectralScorer::new();

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
                let window_spectra: Vec<&Spectrum> = spectrum_indices
                    .iter()
                    .map(|&idx| &spectra[idx])
                    .collect();

                if window_spectra.is_empty() {
                    return Vec::new();
                }

                // === BATCH XCORR PREPROCESSING (pyXcorrDIA approach) ===
                // Preprocess all spectra for XCorr ONCE per window
                let spec_xcorr_preprocessed: Vec<Vec<f64>> = window_spectra
                    .iter()
                    .map(|spec| xcorr_scorer.preprocess_spectrum_for_xcorr(spec))
                    .collect();

                // Preprocess all library entries for XCorr ONCE per window
                let entry_xcorr_preprocessed: HashMap<u32, Vec<f64>> = window_entries
                    .iter()
                    .map(|entry| (entry.id, xcorr_scorer.preprocess_library_for_xcorr(entry)))
                    .collect();

                // Score each library entry against all spectra in this window
                let mut window_matches: Vec<CalibrationMatch> = Vec::new();

                for entry in &window_entries {
                    // Get preprocessed XCorr vector for this entry
                    let entry_xcorr_vec = match entry_xcorr_preprocessed.get(&entry.id) {
                        Some(v) => v,
                        None => continue,
                    };

                    // === FIND BEST RT USING XCORR (not LibCosine) ===
                    // This tests whether XCorr is fundamentally better at RT selection
                    let mut best_xcorr = 0.0f64;
                    let mut best_local_idx = 0usize;  // Index into window_spectra
                    let mut best_spec_idx = 0usize;   // Index into global spectra
                    let mut best_measured_rt = 0.0f64;

                    for (local_idx, &spec) in window_spectra.iter().enumerate() {
                        let rt_diff = (spec.retention_time - entry.retention_time).abs();
                        if rt_diff > rt_tolerance {
                            continue;
                        }

                        // Calculate XCorr from preprocessed vectors
                        let xcorr = SpectralScorer::xcorr_from_preprocessed(
                            &spec_xcorr_preprocessed[local_idx],
                            entry_xcorr_vec,
                        );

                        if xcorr > best_xcorr {
                            best_xcorr = xcorr;
                            best_local_idx = local_idx;
                            best_spec_idx = spectrum_indices[local_idx];
                            best_measured_rt = spec.retention_time;
                        }
                    }

                    if best_xcorr > 0.0 {
                        // Now calculate LibCosine at the RT selected by XCorr
                        let best_spec = window_spectra[best_local_idx];
                        let (libcosine_score, ms2_errors) = scorer.score_with_errors(entry, best_spec);
                        let n_matched = ms2_errors.len();

                        // Extract M+0 from MS1 spectrum (pyXcorrDIA approach)
                        let (ms1_ppm_error, observed_precursor_mz, isotope_cosine_score) = if let Some(ms1_spec) = ms1.find_nearest(best_measured_rt) {
                            // Extract isotope envelope from MS1
                            let envelope = IsotopeEnvelope::extract(
                                ms1_spec,
                                entry.precursor_mz,
                                entry.charge,
                                precursor_tolerance_ppm,
                            );

                            if envelope.has_m0() {
                                // Calculate isotope cosine score using exact elemental composition
                                let isotope_score = peptide_isotope_cosine(&entry.sequence, &envelope.intensities);

                                // Use actual M+0 peak for calibration
                                (
                                    envelope.m0_ppm_error(entry.precursor_mz),
                                    envelope.m0_observed_mz,
                                    isotope_score,
                                )
                            } else {
                                // M+0 not found - no MS1 calibration for this match
                                (None, None, None)
                            }
                        } else {
                            // No MS1 spectrum found - fall back to isolation window (not accurate)
                            let iso_center = spectra[best_spec_idx].isolation_window.center;
                            let ppm = ((iso_center - entry.precursor_mz) / entry.precursor_mz) * 1e6;
                            (Some(ppm), Some(iso_center), None)
                        };

                        // Calculate average MS2 error
                        let avg_ms2_error = if !ms2_errors.is_empty() {
                            Some(ms2_errors.iter().sum::<f64>() / ms2_errors.len() as f64)
                        } else {
                            None
                        };

                        window_matches.push(CalibrationMatch {
                            entry_id: entry.id,
                            is_decoy: entry.is_decoy,
                            library_rt: entry.retention_time,
                            measured_rt: best_measured_rt,
                            score: libcosine_score,  // LibCosine at XCorr-selected RT
                            ms1_ppm_error,
                            library_precursor_mz: entry.precursor_mz,
                            observed_precursor_mz,
                            ms2_mass_errors: ms2_errors,
                            avg_ms2_error,
                            n_matched_fragments: n_matched,
                            n_library_fragments: entry.fragments.len(),
                            xcorr_score: best_xcorr,  // XCorr used for RT selection
                            isotope_cosine_score,
                            sequence: entry.sequence.clone(),
                            charge: entry.charge,
                            scan_number: spectra[best_spec_idx].scan_number,
                        });
                    }
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

                let window_spectra: Vec<&Spectrum> = spectrum_indices
                    .iter()
                    .map(|&idx| &spectra[idx])
                    .collect();

                if window_spectra.is_empty() {
                    return Vec::new();
                }

                // === BATCH XCORR PREPROCESSING (pyXcorrDIA approach) ===
                // Preprocess all spectra for XCorr ONCE per window
                let spec_xcorr_preprocessed: Vec<Vec<f64>> = window_spectra
                    .iter()
                    .map(|spec| xcorr_scorer.preprocess_spectrum_for_xcorr(spec))
                    .collect();

                // Preprocess all library entries for XCorr ONCE per window
                let entry_xcorr_preprocessed: HashMap<u32, Vec<f64>> = window_entries
                    .iter()
                    .map(|entry| (entry.id, xcorr_scorer.preprocess_library_for_xcorr(entry)))
                    .collect();

                let mut window_matches: Vec<CalibrationMatch> = Vec::new();

                for entry in &window_entries {
                    // Get preprocessed XCorr vector for this entry
                    let entry_xcorr_vec = match entry_xcorr_preprocessed.get(&entry.id) {
                        Some(v) => v,
                        None => continue,
                    };

                    // === FIND BEST RT USING XCORR (not LibCosine) ===
                    let mut best_xcorr = 0.0f64;
                    let mut best_local_idx = 0usize;
                    let mut best_spec_idx = 0usize;
                    let mut best_measured_rt = 0.0f64;

                    for (local_idx, &spec) in window_spectra.iter().enumerate() {
                        let rt_diff = (spec.retention_time - entry.retention_time).abs();
                        if rt_diff > rt_tolerance {
                            continue;
                        }

                        let xcorr = SpectralScorer::xcorr_from_preprocessed(
                            &spec_xcorr_preprocessed[local_idx],
                            entry_xcorr_vec,
                        );

                        if xcorr > best_xcorr {
                            best_xcorr = xcorr;
                            best_local_idx = local_idx;
                            best_spec_idx = spectrum_indices[local_idx];
                            best_measured_rt = spec.retention_time;
                        }
                    }

                    if best_xcorr > 0.0 {
                        // Now calculate LibCosine at the RT selected by XCorr
                        let best_spec = window_spectra[best_local_idx];
                        let (libcosine_score, ms2_errors) = scorer.score_with_errors(entry, best_spec);
                        let n_matched = ms2_errors.len();

                        let observed_precursor_mz = Some(spectra[best_spec_idx].isolation_window.center);
                        let ms1_ppm_error = observed_precursor_mz.map(|obs_mz| {
                            ((obs_mz - entry.precursor_mz) / entry.precursor_mz) * 1e6
                        });

                        let avg_ms2_error = if !ms2_errors.is_empty() {
                            Some(ms2_errors.iter().sum::<f64>() / ms2_errors.len() as f64)
                        } else {
                            None
                        };

                        window_matches.push(CalibrationMatch {
                            entry_id: entry.id,
                            is_decoy: entry.is_decoy,
                            library_rt: entry.retention_time,
                            measured_rt: best_measured_rt,
                            score: libcosine_score,  // LibCosine at XCorr-selected RT
                            ms1_ppm_error,
                            library_precursor_mz: entry.precursor_mz,
                            observed_precursor_mz,
                            ms2_mass_errors: ms2_errors,
                            avg_ms2_error,
                            n_matched_fragments: n_matched,
                            n_library_fragments: entry.fragments.len(),
                            xcorr_score: best_xcorr,  // XCorr used for RT selection
                            isotope_cosine_score: None,
                            sequence: entry.sequence.clone(),
                            charge: entry.charge,
                            scan_number: spectra[best_spec_idx].scan_number,
                        });
                    }
                }

                window_matches
            })
            .collect()
    };

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
    results.sort_by(|a, b| b.xcorr_score.partial_cmp(&a.xcorr_score).unwrap_or(std::cmp::Ordering::Equal));

    // Log scoring statistics
    let total_ms2_errors: usize = results.iter().map(|m| m.ms2_mass_errors.len()).sum();
    let ms1_calibrated = results.iter().filter(|m| m.ms1_ppm_error.is_some()).count();
    log::info!(
        "LibCosine scoring complete: {} matches ({} with MS1 calibration, {} MS2 errors)",
        results.len(),
        ms1_calibrated,
        total_ms2_errors
    );

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{FragmentAnnotation, IsolationWindow, LibraryFragment};

    fn make_test_library_entry(id: u32, mz: f64, rt: f64, fragments: Vec<(f64, f32)>) -> LibraryEntry {
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
        assert!(scores[[0, 0]] > 0.99, "Expected ~1.0, got {}", scores[[0, 0]]);
    }

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
        assert!(scores[[0, 0]] < 0.01, "Expected ~0.0, got {}", scores[[0, 0]]);
    }

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
        let spectrum = make_test_spectrum(
            1,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)],
        );

        let scorer = LibCosineScorer::hram(10.0); // 10 ppm
        let score = scorer.score(&entry, &spectrum);

        // Perfect match should give score = 1.0
        assert!(score > 0.99, "Expected ~1.0, got {}", score);
    }

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
        assert!(score_tight < 0.01, "Expected no match at 5 ppm, got {}", score_tight);
    }

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

    #[test]
    fn test_libcosine_da_tolerance() {
        // Library fragment at 500.0 m/z
        let entry = make_test_library_entry(1, 500.0, 10.0, vec![(500.0, 100.0)]);

        // Spectrum with peak at 500.25 m/z (0.25 Da offset)
        let spectrum = make_test_spectrum(1, 10.0, vec![(500.25, 100.0)]);

        // 0.3 Da tolerance should match
        let scorer_loose = LibCosineScorer::unit_resolution(0.3);
        let score_loose = scorer_loose.score(&entry, &spectrum);
        assert!(score_loose > 0.99, "Expected match at 0.3 Da, got {}", score_loose);

        // 0.2 Da tolerance should NOT match
        let scorer_tight = LibCosineScorer::unit_resolution(0.2);
        let score_tight = scorer_tight.score(&entry, &spectrum);
        assert!(score_tight < 0.01, "Expected no match at 0.2 Da, got {}", score_tight);
    }

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
            mzs: vec![499.996, 500.000, 500.004],     // All within 10 ppm
            intensities: vec![200.0, 50.0, 150.0],    // 200 is highest, but 50 is closest
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
}
