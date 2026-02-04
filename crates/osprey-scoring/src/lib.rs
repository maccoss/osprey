//! Osprey Scoring - Feature extraction and peptide scoring
//!
//! This crate provides:
//! - Feature extraction from coefficient time series and spectra
//! - Spectral similarity metrics (LibCosine, XCorr)
//! - Hyperscore calculation
//! - Decoy generation with enzyme-aware sequence reversal
//! - **Batch scoring with BLAS acceleration** (10-20× faster)
//! - **Streaming pipeline** for memory-efficient processing (optional)
//!
//! ## Batch Scoring
//!
//! For high-performance scoring of many library entries against many spectra,
//! use the [`batch`] module:
//!
//! ```ignore
//! use osprey_scoring::batch::{BatchScorer, PreprocessedLibrary, PreprocessedSpectra};
//!
//! let scorer = BatchScorer::new();
//! let lib = scorer.preprocess_library(&library_entries);
//! let spec = scorer.preprocess_spectra(&spectra);
//!
//! // Score all pairs with single BLAS call
//! let scores = scorer.score_all(&lib, &spec);
//! ```
//!
//! ## Streaming Pipeline (requires `streaming` feature)
//!
//! For processing large mzML files with overlapped I/O and preprocessing:
//!
//! ```ignore
//! use osprey_scoring::pipeline::{PipelineCoordinator, run_streaming_pipeline};
//!
//! // From file (async)
//! let matches = run_streaming_pipeline("file.mzML", &library, config).await?;
//!
//! // From pre-loaded spectra
//! let coordinator = PipelineCoordinator::new();
//! let matches = coordinator.run_on_spectra(&spectra, &library);
//! ```

pub mod batch;
pub mod pipeline;

use osprey_core::{
    FeatureSet, FragmentAnnotation, IonType, LibraryEntry, LibraryFragment, Modification, Result,
    Spectrum,
};
use rayon::prelude::*;
use std::collections::HashMap;

/// Spectral similarity scorer implementing LibCosine and XCorr
///
/// Based on pyXcorrDIA scoring methods:
/// - **LibCosine**: sqrt(intensity) × m/z² preprocessing with L2 normalization
/// - **XCorr**: Comet-style windowing normalization with sliding window subtraction
#[derive(Debug, Clone)]
pub struct SpectralScorer {
    /// Mass tolerance for matching fragments (Da)
    tolerance_da: f64,
    /// Mass tolerance for matching fragments (ppm)
    tolerance_ppm: f64,
    /// Number of bins for XCorr (0-2000 m/z range, Comet-style)
    num_bins: usize,
    /// Bin width for XCorr
    bin_width: f64,
    /// Minimum m/z for binning (0.0 to match Comet/pyXcorrDIA)
    #[allow(dead_code)]
    min_mz: f64,
}

impl Default for SpectralScorer {
    fn default() -> Self {
        let bin_width = 1.0005079; // Comet default bin width
        let max_mz = 2000.0;
        // Comet BIN macro: BIN(mass) = (int)(mass / bin_width + (1 - offset))
        // With offset = 0.4: num_bins = (int)(2000 / 1.0005079 + 0.6) + 1
        let num_bins = ((max_mz / bin_width) + 0.6) as usize + 1;
        Self {
            tolerance_da: 0.5, // For unit resolution
            tolerance_ppm: 20.0, // For HRAM
            num_bins,
            bin_width,
            min_mz: 0.0, // Start at 0 m/z (Comet/pyXcorrDIA style)
        }
    }
}

impl SpectralScorer {
    /// Create a new spectral scorer with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Set Da tolerance for unit resolution matching
    pub fn with_tolerance_da(mut self, tolerance: f64) -> Self {
        self.tolerance_da = tolerance;
        self
    }

    /// Set ppm tolerance for HRAM matching
    pub fn with_tolerance_ppm(mut self, tolerance: f64) -> Self {
        self.tolerance_ppm = tolerance;
        self
    }

    /// Compute LibCosine score between observed spectrum and library entry
    ///
    /// LibCosine uses SMZ preprocessing: sqrt(intensity) × m/z²
    /// Then computes cosine similarity between normalized vectors.
    pub fn lib_cosine(&self, observed: &Spectrum, library: &LibraryEntry) -> SpectralScore {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return SpectralScore::default();
        }

        // Match library fragments to observed peaks
        let matches = self.match_fragments(observed, library);

        if matches.is_empty() {
            return SpectralScore::default();
        }

        // Preprocess library fragments: sqrt(intensity) × m/z²
        let mut lib_preprocessed: Vec<f64> = Vec::with_capacity(matches.len());
        let mut obs_preprocessed: Vec<f64> = Vec::with_capacity(matches.len());

        // Keep raw intensities for correlation computation
        let mut lib_intensities: Vec<f64> = Vec::with_capacity(matches.len());
        let mut obs_intensities: Vec<f64> = Vec::with_capacity(matches.len());

        for m in &matches {
            // Library: sqrt(relative_intensity) × m/z²
            let lib_val = (m.lib_intensity as f64).sqrt() * m.lib_mz.powi(2);
            lib_preprocessed.push(lib_val);
            lib_intensities.push(m.lib_intensity as f64);

            // Observed: sqrt(intensity) × m/z²
            let obs_val = (m.obs_intensity as f64).sqrt() * m.obs_mz.powi(2);
            obs_preprocessed.push(obs_val);
            obs_intensities.push(m.obs_intensity as f64);
        }

        // L2 normalize both vectors
        let lib_norm = lib_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();
        let obs_norm = obs_preprocessed.iter().map(|x| x * x).sum::<f64>().sqrt();

        if lib_norm < 1e-10 || obs_norm < 1e-10 {
            return SpectralScore::default();
        }

        for v in lib_preprocessed.iter_mut() {
            *v /= lib_norm;
        }
        for v in obs_preprocessed.iter_mut() {
            *v /= obs_norm;
        }

        // Cosine similarity (dot product of normalized vectors)
        let dot_product: f64 = lib_preprocessed
            .iter()
            .zip(obs_preprocessed.iter())
            .map(|(a, b)| a * b)
            .sum();

        // Compute additional metrics
        let n_matched = matches.len() as u32;
        let n_library = library.fragments.len() as u32;
        let fragment_coverage = n_matched as f64 / n_library as f64;

        // Compute explained intensity (sum of matched intensities / total observed intensity)
        let matched_intensity: f64 = matches.iter().map(|m| m.obs_intensity as f64).sum();
        let total_intensity: f64 = observed.intensities.iter().map(|&i| i as f64).sum();
        let explained_intensity = if total_intensity > 0.0 {
            matched_intensity / total_intensity
        } else {
            0.0
        };

        // Compute Pearson correlation
        let pearson_correlation = Self::pearson_correlation(&lib_intensities, &obs_intensities);

        // Compute Spearman correlation
        let spearman_correlation = Self::spearman_correlation(&lib_intensities, &obs_intensities);

        // Compute consecutive ion series
        let consecutive_ions = self.longest_consecutive_ions(library, &matches);

        // Compute sequence coverage (backbone cleavage sites)
        let sequence_coverage = self.compute_sequence_coverage(library, &matches);

        // Compute base peak rank
        let base_peak_rank = self.compute_base_peak_rank(&matches, &lib_intensities, &obs_intensities);

        // Compute top-3 matches
        let top3_matches = self.compute_top3_matches(library, &matches);

        SpectralScore {
            lib_cosine: dot_product,
            xcorr: 0.0, // Not computed in this method
            dot_product,
            n_matched,
            n_library,
            fragment_coverage,
            explained_intensity,
            pearson_correlation,
            spearman_correlation,
            consecutive_ions,
            sequence_coverage,
            base_peak_rank,
            top3_matches,
        }
    }

    /// Compute Pearson correlation coefficient
    fn pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
        if x.len() < 2 || x.len() != y.len() {
            return 0.0;
        }

        let n = x.len() as f64;
        let mean_x = x.iter().sum::<f64>() / n;
        let mean_y = y.iter().sum::<f64>() / n;

        let mut sum_xy = 0.0;
        let mut sum_x2 = 0.0;
        let mut sum_y2 = 0.0;

        for (xi, yi) in x.iter().zip(y.iter()) {
            let dx = xi - mean_x;
            let dy = yi - mean_y;
            sum_xy += dx * dy;
            sum_x2 += dx * dx;
            sum_y2 += dy * dy;
        }

        let denom = (sum_x2 * sum_y2).sqrt();
        if denom < 1e-10 {
            return 0.0;
        }

        sum_xy / denom
    }

    /// Compute Spearman rank correlation coefficient
    fn spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
        if x.len() < 2 || x.len() != y.len() {
            return 0.0;
        }

        // Get ranks
        let rank_x = Self::get_ranks(x);
        let rank_y = Self::get_ranks(y);

        // Pearson correlation of ranks
        Self::pearson_correlation(&rank_x, &rank_y)
    }

    /// Convert values to ranks (1-based, average for ties)
    fn get_ranks(values: &[f64]) -> Vec<f64> {
        let n = values.len();
        let mut indexed: Vec<(usize, f64)> = values.iter().copied().enumerate().collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let mut ranks = vec![0.0; n];
        let mut i = 0;
        while i < n {
            let mut j = i;
            // Find all tied values
            while j < n && (indexed[j].1 - indexed[i].1).abs() < 1e-10 {
                j += 1;
            }
            // Assign average rank to all tied values
            let avg_rank = (i + j + 1) as f64 / 2.0; // 1-based ranks
            for k in i..j {
                ranks[indexed[k].0] = avg_rank;
            }
            i = j;
        }

        ranks
    }

    /// Find longest consecutive b or y ion series
    fn longest_consecutive_ions(&self, library: &LibraryEntry, matches: &[FragmentMatch]) -> u32 {
        // Build sets of matched b and y ion ordinals
        let mut matched_b: Vec<u8> = Vec::new();
        let mut matched_y: Vec<u8> = Vec::new();

        for m in matches {
            // Find the corresponding library fragment annotation
            for frag in &library.fragments {
                if (frag.mz - m.lib_mz).abs() < 0.001 {
                    match frag.annotation.ion_type {
                        IonType::B => matched_b.push(frag.annotation.ordinal),
                        IonType::Y => matched_y.push(frag.annotation.ordinal),
                        _ => {}
                    }
                    break;
                }
            }
        }

        // Find longest consecutive run in each
        let longest_b = Self::longest_consecutive_run(&matched_b);
        let longest_y = Self::longest_consecutive_run(&matched_y);

        longest_b.max(longest_y)
    }

    /// Find longest consecutive run of ordinals
    fn longest_consecutive_run(ordinals: &[u8]) -> u32 {
        if ordinals.is_empty() {
            return 0;
        }

        let mut sorted: Vec<u8> = ordinals.to_vec();
        sorted.sort();
        sorted.dedup();

        let mut max_run = 1;
        let mut current_run = 1;

        for i in 1..sorted.len() {
            if sorted[i] == sorted[i - 1] + 1 {
                current_run += 1;
                max_run = max_run.max(current_run);
            } else {
                current_run = 1;
            }
        }

        max_run
    }

    /// Compute sequence coverage (fraction of backbone cleavage sites covered)
    fn compute_sequence_coverage(&self, library: &LibraryEntry, matches: &[FragmentMatch]) -> f64 {
        let seq_len = library.sequence.len();
        if seq_len < 2 {
            return 0.0;
        }

        // Number of possible cleavage sites is seq_len - 1
        let n_sites = seq_len - 1;
        let mut covered = vec![false; n_sites];

        for m in matches {
            // Find the corresponding library fragment
            for frag in &library.fragments {
                if (frag.mz - m.lib_mz).abs() < 0.001 {
                    let ordinal = frag.annotation.ordinal as usize;
                    if ordinal > 0 && ordinal <= n_sites {
                        // b-ion of ordinal n covers site n-1 (0-indexed)
                        // y-ion of ordinal n covers site seq_len - n - 1 (0-indexed)
                        match frag.annotation.ion_type {
                            IonType::B => {
                                if ordinal <= n_sites {
                                    covered[ordinal - 1] = true;
                                }
                            }
                            IonType::Y => {
                                let site = seq_len - ordinal - 1;
                                if site < n_sites {
                                    covered[site] = true;
                                }
                            }
                            _ => {}
                        }
                    }
                    break;
                }
            }
        }

        let n_covered = covered.iter().filter(|&&x| x).count();
        n_covered as f64 / n_sites as f64
    }

    /// Compute rank of observed base peak in library ordering
    fn compute_base_peak_rank(
        &self,
        matches: &[FragmentMatch],
        lib_intensities: &[f64],
        obs_intensities: &[f64],
    ) -> u32 {
        if matches.is_empty() || obs_intensities.is_empty() {
            return 0;
        }

        // Find index of max observed intensity
        let max_obs_idx = obs_intensities
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i);

        if let Some(obs_idx) = max_obs_idx {
            // Get library intensity of this match
            let lib_int = lib_intensities[obs_idx];

            // Count how many library fragments have higher intensity
            let rank = lib_intensities.iter().filter(|&&x| x > lib_int).count() + 1;
            rank as u32
        } else {
            0
        }
    }

    /// Count how many of the top-3 library peaks are matched
    fn compute_top3_matches(&self, library: &LibraryEntry, matches: &[FragmentMatch]) -> u32 {
        if library.fragments.is_empty() {
            return 0;
        }

        // Get top 3 library fragments by intensity
        let mut lib_sorted: Vec<(f64, f64)> = library
            .fragments
            .iter()
            .map(|f| (f.mz, f.relative_intensity as f64))
            .collect();
        lib_sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        let top3: Vec<f64> = lib_sorted.iter().take(3).map(|(mz, _)| *mz).collect();

        // Count how many are matched
        let mut count = 0;
        for top_mz in &top3 {
            for m in matches {
                if (m.lib_mz - top_mz).abs() < 0.001 {
                    count += 1;
                    break;
                }
            }
        }

        count
    }

    /// Compute XCorr (cross-correlation) score using Comet-style preprocessing
    ///
    /// XCorr applies:
    /// 1. Windowing normalization (normalize to 50.0 within 10 windows)
    /// 2. Sliding window subtraction (offset=75)
    /// 3. Dot product scoring
    pub fn xcorr(&self, observed: &Spectrum, library: &LibraryEntry) -> SpectralScore {
        if library.fragments.is_empty() || observed.mzs.is_empty() {
            return SpectralScore::default();
        }

        // First get LibCosine for additional metrics
        let lib_cosine_score = self.lib_cosine(observed, library);

        // Bin observed spectrum using Comet-style BIN macro
        // BIN(mass) = (int)(mass / bin_width + (1 - offset)) = (int)(mass / bin_width + 0.6)
        let mut obs_binned = vec![0.0f64; self.num_bins];
        for (&mz, &intensity) in observed.mzs.iter().zip(observed.intensities.iter()) {
            let bin = ((mz / self.bin_width) + 0.6) as usize;
            if bin < self.num_bins {
                // Apply sqrt transformation to experimental spectrum
                obs_binned[bin] += (intensity as f64).sqrt();
            }
        }

        // Apply windowing normalization (Comet's MakeCorrData)
        let windowed = self.apply_windowing_normalization(&obs_binned);

        // Apply sliding window subtraction (fast XCorr preprocessing)
        let xcorr_preprocessed = self.apply_sliding_window(&windowed);

        // Bin library spectrum (theoretical) - use unit intensity like Comet
        // Comet BIN macro: BIN(mass) = (int)(mass / bin_width + (1 - offset))
        // With offset = 0.4: bin = (int)(mz / bin_width + 0.6)
        let mut lib_binned = vec![0.0f64; self.num_bins];
        for frag in &library.fragments {
            // Comet-style binning: bin = (int)(mz / bin_width + 0.6)
            let bin = ((frag.mz / self.bin_width) + 0.6) as usize;
            if bin < self.num_bins {
                // Use unit intensity (1.0) for theoretical spectrum, NOT library intensity
                // NO windowing for theoretical - Comet just looks up preprocessed experimental
                // values at fragment bin positions and sums them
                lib_binned[bin] = 1.0;
            }
        }

        // XCorr = dot product with scaling (NO windowing on theoretical)
        // This matches Comet exactly: score = sum(experimental_preprocessed[frag_bins]) * 0.005
        let xcorr_raw: f64 = lib_binned
            .iter()
            .zip(xcorr_preprocessed.iter())
            .map(|(a, b)| a * b)
            .sum();

        // Scale XCorr (pyXcorrDIA uses 0.005 for spectrum-centric)
        let xcorr_scaled = xcorr_raw * 0.005;

        SpectralScore {
            lib_cosine: lib_cosine_score.lib_cosine,
            xcorr: xcorr_scaled,
            dot_product: lib_cosine_score.dot_product,
            n_matched: lib_cosine_score.n_matched,
            n_library: lib_cosine_score.n_library,
            fragment_coverage: lib_cosine_score.fragment_coverage,
            explained_intensity: lib_cosine_score.explained_intensity,
            pearson_correlation: lib_cosine_score.pearson_correlation,
            spearman_correlation: lib_cosine_score.spearman_correlation,
            consecutive_ions: lib_cosine_score.consecutive_ions,
            sequence_coverage: lib_cosine_score.sequence_coverage,
            base_peak_rank: lib_cosine_score.base_peak_rank,
            top3_matches: lib_cosine_score.top3_matches,
        }
    }

    /// Match library fragments to observed peaks
    fn match_fragments(&self, observed: &Spectrum, library: &LibraryEntry) -> Vec<FragmentMatch> {
        let mut matches = Vec::new();

        for frag in &library.fragments {
            // Find best matching observed peak
            let mut best_match: Option<(usize, f64)> = None;
            let mut best_error = f64::MAX;

            for (i, (&mz, &intensity)) in
                observed.mzs.iter().zip(observed.intensities.iter()).enumerate()
            {
                let error_da = (mz - frag.mz).abs();
                let error_ppm = error_da / frag.mz * 1e6;

                // Check if within tolerance (use both Da and ppm, take minimum)
                let within_tolerance =
                    error_da <= self.tolerance_da || error_ppm <= self.tolerance_ppm;

                if within_tolerance && error_da < best_error {
                    best_error = error_da;
                    best_match = Some((i, intensity as f64));
                }
            }

            if let Some((idx, obs_intensity)) = best_match {
                matches.push(FragmentMatch {
                    lib_mz: frag.mz,
                    lib_intensity: frag.relative_intensity,
                    obs_mz: observed.mzs[idx],
                    obs_intensity: obs_intensity as f32,
                });
            }
        }

        matches
    }

    /// Apply Comet-style windowing normalization
    ///
    /// Divides spectrum into 10 windows and normalizes each to max=50.0
    fn apply_windowing_normalization(&self, spectrum: &[f64]) -> Vec<f64> {
        let mut result = vec![0.0f64; spectrum.len()];
        let num_windows = 10;
        let window_size = (spectrum.len() / num_windows) + 1;

        // Find global max for threshold
        let global_max = spectrum.iter().cloned().fold(0.0f64, f64::max);
        let threshold = global_max * 0.05;

        for window_idx in 0..num_windows {
            let start = window_idx * window_size;
            let end = ((window_idx + 1) * window_size).min(spectrum.len());

            // Find max in this window
            let mut window_max = 0.0f64;
            for i in start..end {
                if spectrum[i] > window_max {
                    window_max = spectrum[i];
                }
            }

            // Normalize this window to 50.0
            if window_max > 0.0 {
                let norm_factor = 50.0 / window_max;
                for i in start..end {
                    if spectrum[i] > threshold {
                        result[i] = spectrum[i] * norm_factor;
                    }
                }
            }
        }

        result
    }

    /// Apply sliding window subtraction for fast XCorr
    ///
    /// Subtracts local average (excluding center) with offset=75
    fn apply_sliding_window(&self, spectrum: &[f64]) -> Vec<f64> {
        let mut result = vec![0.0f64; spectrum.len()];
        let offset: i64 = 75;
        let window_range = 2 * offset + 1;
        let norm_factor = 1.0 / (window_range - 1) as f64;

        let n = spectrum.len() as i64;

        for i in 0..spectrum.len() {
            let i64_i = i as i64;
            let mut sum = 0.0;

            // Sum values in window, excluding center
            for j in (i64_i - offset)..=(i64_i + offset) {
                if j >= 0 && j < n && j != i64_i {
                    sum += spectrum[j as usize];
                }
            }

            // Subtract local average from center value
            result[i] = spectrum[i] - sum * norm_factor;
        }

        result
    }

    // ========================================================================
    // Batch XCorr preprocessing methods (for BLAS-accelerated scoring)
    // ========================================================================

    /// Preprocess a spectrum for XCorr (Comet-style)
    ///
    /// Returns a preprocessed vector ready for dot product with library vectors.
    /// Preprocessing: bin → sqrt → windowing → flanking subtraction
    ///
    /// This allows precomputing once per spectrum and reusing across library entries.
    pub fn preprocess_spectrum_for_xcorr(&self, spectrum: &Spectrum) -> Vec<f64> {
        // Bin observed spectrum with sqrt transformation
        let mut binned = vec![0.0f64; self.num_bins];
        for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            let bin = ((mz / self.bin_width) + 0.6) as usize;
            if bin < self.num_bins {
                binned[bin] += (intensity as f64).sqrt();
            }
        }

        // Apply windowing normalization
        let windowed = self.apply_windowing_normalization(&binned);

        // Apply flanking bin subtraction
        self.apply_sliding_window(&windowed)
    }

    /// Preprocess a library entry for XCorr (Comet-style)
    ///
    /// Returns a preprocessed vector ready for dot product with spectrum vectors.
    /// Preprocessing: bin with unit intensity → windowing (NO flanking subtraction)
    ///
    /// This allows precomputing once per library entry and reusing across spectra.
    pub fn preprocess_library_for_xcorr(&self, entry: &LibraryEntry) -> Vec<f64> {
        // Bin library fragments with unit intensity
        // Comet-style: theoretical spectrum is NOT windowed, just unit intensities at fragment bins
        // The score is simply: sum(experimental_preprocessed[frag_bins]) * 0.005
        let mut binned = vec![0.0f64; self.num_bins];
        for frag in &entry.fragments {
            let bin = ((frag.mz / self.bin_width) + 0.6) as usize;
            if bin < self.num_bins {
                // Use unit intensity (1.0), NOT library intensity (Comet-style)
                // NO windowing applied - this matches Comet exactly
                binned[bin] = 1.0;
            }
        }

        // Return binned vector directly - NO windowing for theoretical spectrum
        // Comet just looks up preprocessed experimental values at fragment bin positions
        binned
    }

    /// Compute XCorr from preprocessed vectors
    ///
    /// This is a simple dot product with scaling, used after preprocessing.
    #[inline]
    pub fn xcorr_from_preprocessed(spectrum_preprocessed: &[f64], library_preprocessed: &[f64]) -> f64 {
        let min_len = spectrum_preprocessed.len().min(library_preprocessed.len());
        let raw: f64 = spectrum_preprocessed[..min_len]
            .iter()
            .zip(library_preprocessed[..min_len].iter())
            .map(|(s, l)| s * l)
            .sum();

        // Scale XCorr (pyXcorrDIA uses 0.005 for spectrum-centric)
        raw * 0.005
    }

    /// Get the number of bins used for XCorr
    pub fn num_bins(&self) -> usize {
        self.num_bins
    }
}

/// Result of spectral similarity scoring
#[derive(Debug, Clone, Default)]
pub struct SpectralScore {
    /// LibCosine score (0-1, cosine similarity with SMZ preprocessing)
    pub lib_cosine: f64,
    /// XCorr score (Comet-style cross-correlation)
    pub xcorr: f64,
    /// Raw dot product
    pub dot_product: f64,
    /// Number of matched fragments
    pub n_matched: u32,
    /// Number of library fragments
    pub n_library: u32,
    /// Fragment coverage (n_matched / n_library)
    pub fragment_coverage: f64,
    /// Fraction of observed intensity explained by matches
    pub explained_intensity: f64,
    /// Pearson correlation between library and observed intensities
    pub pearson_correlation: f64,
    /// Spearman rank correlation between library and observed intensities
    pub spearman_correlation: f64,
    /// Longest consecutive b or y ion series
    pub consecutive_ions: u32,
    /// Sequence coverage (backbone cleavages covered)
    pub sequence_coverage: f64,
    /// Rank of base peak (highest intensity) in library
    pub base_peak_rank: u32,
    /// Number of top-3 library peaks that matched
    pub top3_matches: u32,
}

/// Internal struct for fragment matching
#[derive(Debug)]
struct FragmentMatch {
    lib_mz: f64,
    lib_intensity: f32,
    obs_mz: f64,
    obs_intensity: f32,
}

/// Spectrum aggregator for combining spectra across peak apex region
///
/// FR-5.2.1: Aggregate observed spectrum across peak apex region
#[derive(Debug, Clone)]
pub struct SpectrumAggregator {
    /// Mass tolerance for combining peaks (Da)
    tolerance_da: f64,
    /// Number of scans to include around apex
    apex_window: usize,
}

impl Default for SpectrumAggregator {
    fn default() -> Self {
        Self {
            tolerance_da: 0.5,
            apex_window: 5, // 2 scans on each side of apex
        }
    }
}

impl SpectrumAggregator {
    /// Create a new spectrum aggregator
    pub fn new() -> Self {
        Self::default()
    }

    /// Set mass tolerance for peak combining
    pub fn with_tolerance_da(mut self, tolerance: f64) -> Self {
        self.tolerance_da = tolerance;
        self
    }

    /// Set number of scans to aggregate around apex
    pub fn with_apex_window(mut self, window: usize) -> Self {
        self.apex_window = window;
        self
    }

    /// Aggregate spectra around the peak apex
    ///
    /// Combines peaks from multiple spectra within the apex region
    /// by summing intensities of peaks at similar m/z values.
    pub fn aggregate(&self, spectra: &[Spectrum], apex_rt: f64) -> Spectrum {
        if spectra.is_empty() {
            return Spectrum::new(
                0,
                apex_rt,
                osprey_core::IsolationWindow::symmetric(0.0, 0.0),
            );
        }

        // Find apex scan index
        let apex_idx = spectra
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1.retention_time - apex_rt)
                    .abs()
                    .partial_cmp(&(b.1.retention_time - apex_rt).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Select scans in the window around apex
        let start_idx = apex_idx.saturating_sub(self.apex_window / 2);
        let end_idx = (apex_idx + self.apex_window / 2 + 1).min(spectra.len());
        let window_spectra = &spectra[start_idx..end_idx];

        if window_spectra.is_empty() {
            return spectra[apex_idx].clone();
        }

        // Collect all peaks from window
        let mut all_peaks: Vec<(f64, f64)> = Vec::new();
        for spectrum in window_spectra {
            for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
                all_peaks.push((mz, intensity as f64));
            }
        }

        // Sort by m/z
        all_peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Merge peaks within tolerance
        let mut merged_mzs: Vec<f64> = Vec::new();
        let mut merged_intensities: Vec<f32> = Vec::new();

        let mut i = 0;
        while i < all_peaks.len() {
            let mut sum_mz = all_peaks[i].0;
            let mut sum_intensity = all_peaks[i].1;
            let mut count = 1.0;
            let base_mz = all_peaks[i].0;

            // Merge all peaks within tolerance
            let mut j = i + 1;
            while j < all_peaks.len() && (all_peaks[j].0 - base_mz).abs() <= self.tolerance_da {
                sum_mz += all_peaks[j].0;
                sum_intensity += all_peaks[j].1;
                count += 1.0;
                j += 1;
            }

            // Use intensity-weighted m/z average
            merged_mzs.push(sum_mz / count);
            merged_intensities.push(sum_intensity as f32);

            i = j;
        }

        // Create aggregated spectrum using the apex spectrum's metadata
        let apex_spectrum = &spectra[apex_idx];
        Spectrum {
            scan_number: apex_spectrum.scan_number,
            retention_time: apex_spectrum.retention_time,
            precursor_mz: apex_spectrum.precursor_mz,
            isolation_window: apex_spectrum.isolation_window,
            mzs: merged_mzs,
            intensities: merged_intensities,
        }
    }

    /// Aggregate spectra weighted by coefficient values
    ///
    /// Uses ridge regression coefficients to weight the contribution
    /// of each scan to the aggregated spectrum.
    pub fn aggregate_weighted(
        &self,
        spectra: &[Spectrum],
        coefficients: &[(f64, f64)], // (rt, coefficient) pairs
        apex_rt: f64,
    ) -> Spectrum {
        if spectra.is_empty() || coefficients.is_empty() {
            return Spectrum::new(
                0,
                apex_rt,
                osprey_core::IsolationWindow::symmetric(0.0, 0.0),
            );
        }

        // Find apex scan index
        let apex_idx = spectra
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1.retention_time - apex_rt)
                    .abs()
                    .partial_cmp(&(b.1.retention_time - apex_rt).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0);

        // Select scans in the window around apex
        let start_idx = apex_idx.saturating_sub(self.apex_window / 2);
        let end_idx = (apex_idx + self.apex_window / 2 + 1).min(spectra.len());
        let window_spectra = &spectra[start_idx..end_idx];

        if window_spectra.is_empty() {
            return spectra[apex_idx].clone();
        }

        // Get coefficient weights for each spectrum
        let mut weights: Vec<f64> = Vec::new();
        for spectrum in window_spectra {
            // Find corresponding coefficient
            let weight = coefficients
                .iter()
                .min_by(|a, b| {
                    (a.0 - spectrum.retention_time)
                        .abs()
                        .partial_cmp(&(b.0 - spectrum.retention_time).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(_, c)| c.max(0.0)) // Use max(0, c) for non-negative weights
                .unwrap_or(0.0);
            weights.push(weight);
        }

        // Normalize weights
        let total_weight: f64 = weights.iter().sum();
        if total_weight < 1e-10 {
            // All weights zero, fall back to simple aggregation
            return self.aggregate(spectra, apex_rt);
        }
        for w in weights.iter_mut() {
            *w /= total_weight;
        }

        // Collect all peaks with weights
        let mut all_peaks: Vec<(f64, f64)> = Vec::new();
        for (spectrum, &weight) in window_spectra.iter().zip(weights.iter()) {
            for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
                all_peaks.push((mz, intensity as f64 * weight));
            }
        }

        // Sort by m/z
        all_peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        // Merge peaks within tolerance
        let mut merged_mzs: Vec<f64> = Vec::new();
        let mut merged_intensities: Vec<f32> = Vec::new();

        let mut i = 0;
        while i < all_peaks.len() {
            let mut sum_mz_weighted = all_peaks[i].0 * all_peaks[i].1;
            let mut sum_intensity = all_peaks[i].1;
            let base_mz = all_peaks[i].0;

            // Merge all peaks within tolerance
            let mut j = i + 1;
            while j < all_peaks.len() && (all_peaks[j].0 - base_mz).abs() <= self.tolerance_da {
                sum_mz_weighted += all_peaks[j].0 * all_peaks[j].1;
                sum_intensity += all_peaks[j].1;
                j += 1;
            }

            if sum_intensity > 1e-10 {
                // Use intensity-weighted m/z average
                merged_mzs.push(sum_mz_weighted / sum_intensity);
                merged_intensities.push(sum_intensity as f32);
            }

            i = j;
        }

        // Create aggregated spectrum
        let apex_spectrum = &spectra[apex_idx];
        Spectrum {
            scan_number: apex_spectrum.scan_number,
            retention_time: apex_spectrum.retention_time,
            precursor_mz: apex_spectrum.precursor_mz,
            isolation_window: apex_spectrum.isolation_window,
            mzs: merged_mzs,
            intensities: merged_intensities,
        }
    }
}

/// Feature extractor for peptide candidates
#[derive(Debug)]
pub struct FeatureExtractor {
    /// Spectral scorer for computing similarity metrics
    spectral_scorer: SpectralScorer,
    /// Spectrum aggregator for combining scans
    spectrum_aggregator: SpectrumAggregator,
}

impl Default for FeatureExtractor {
    fn default() -> Self {
        Self {
            spectral_scorer: SpectralScorer::default(),
            spectrum_aggregator: SpectrumAggregator::default(),
        }
    }
}

impl FeatureExtractor {
    /// Create a new feature extractor
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a feature extractor with custom spectral scorer
    pub fn with_scorer(spectral_scorer: SpectralScorer) -> Self {
        Self {
            spectral_scorer,
            spectrum_aggregator: SpectrumAggregator::default(),
        }
    }

    /// Create a feature extractor with custom scorer and aggregator
    pub fn with_scorer_and_aggregator(
        spectral_scorer: SpectralScorer,
        spectrum_aggregator: SpectrumAggregator,
    ) -> Self {
        Self {
            spectral_scorer,
            spectrum_aggregator,
        }
    }

    /// Get the spectrum aggregator
    pub fn aggregator(&self) -> &SpectrumAggregator {
        &self.spectrum_aggregator
    }

    /// Extract features using spectrum aggregation across peak apex region
    ///
    /// FR-5.2.1: This method aggregates spectra around the peak apex before
    /// computing spectral features, providing more robust scoring.
    ///
    /// Parameters:
    /// - `entry`: The library entry being scored
    /// - `coefficient_series`: Vec of (retention_time, coefficient) pairs
    /// - `spectra`: All spectra in the peak region
    /// - `expected_rt`: Optional expected retention time
    pub fn extract_with_aggregation(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        spectra: &[Spectrum],
        expected_rt: Option<f64>,
    ) -> FeatureSet {
        // Find apex RT from coefficient series
        let apex_rt = coefficient_series
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(rt, _)| *rt)
            .unwrap_or(entry.retention_time);

        // Aggregate spectra around apex (FR-5.2.1)
        let aggregated = self.spectrum_aggregator.aggregate_weighted(
            spectra,
            coefficient_series,
            apex_rt,
        );

        // Extract features using the aggregated spectrum
        self.extract_with_expected_rt(entry, coefficient_series, Some(&aggregated), expected_rt)
    }

    /// Extract features for a peptide candidate
    ///
    /// Computes both chromatographic features (from coefficient series) and
    /// spectral features (from apex spectrum comparison to library).
    ///
    /// Parameters:
    /// - `entry`: The library entry being scored
    /// - `coefficient_series`: Vec of (retention_time, coefficient) pairs
    /// - `apex_spectrum`: Optional spectrum at peak apex for spectral scoring
    /// - `expected_rt`: Optional expected retention time for RT deviation calculation
    pub fn extract(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        apex_spectrum: Option<&Spectrum>,
    ) -> FeatureSet {
        self.extract_with_expected_rt(entry, coefficient_series, apex_spectrum, None)
    }

    /// Extract features with expected RT for deviation calculation
    pub fn extract_with_expected_rt(
        &self,
        entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)],
        apex_spectrum: Option<&Spectrum>,
        expected_rt: Option<f64>,
    ) -> FeatureSet {
        let mut features = FeatureSet::default();

        // Chromatographic features from coefficient series (FR-5.1.*)
        if !coefficient_series.is_empty() {
            // FR-5.1.1: Peak apex coefficient (maximum value)
            let (apex_idx, apex_value) = coefficient_series
                .iter()
                .enumerate()
                .max_by(|a, b| a.1 .1.partial_cmp(&b.1 .1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, (_, c))| (i, *c))
                .unwrap_or((0, 0.0));
            features.peak_apex = apex_value;

            let apex_rt = coefficient_series.get(apex_idx).map(|(rt, _)| *rt).unwrap_or(0.0);

            // FR-5.1.2: Integrated peak area (AUC of coefficients)
            features.peak_area = coefficient_series.iter().map(|(_, c)| c).sum();

            // FR-5.1.8: Number of contributing scans
            features.n_contributing_scans = coefficient_series
                .iter()
                .filter(|(_, c)| *c > 0.0)
                .count() as u32;

            // FR-5.1.4: Peak width (FWHM)
            let half_max = features.peak_apex / 2.0;
            let above_half: Vec<_> = coefficient_series
                .iter()
                .filter(|(_, c)| *c >= half_max)
                .collect();
            if above_half.len() >= 2 {
                features.peak_width = above_half.last().unwrap().0 - above_half.first().unwrap().0;
            }

            // FR-5.1.5: Peak symmetry (leading/trailing ratio)
            features.peak_symmetry = self.compute_peak_symmetry(coefficient_series, apex_idx);

            // FR-5.1.6 & FR-5.1.7: RT deviation
            if let Some(expected) = expected_rt.or(Some(entry.retention_time)) {
                features.rt_deviation = apex_rt - expected;
                // Normalized deviation (use FWHM as uncertainty proxy)
                if features.peak_width > 0.0 {
                    features.rt_deviation_normalized = features.rt_deviation / features.peak_width;
                }
            }

            // FR-5.1.9: Coefficient stability (variance near apex)
            features.coefficient_stability = self.compute_coefficient_stability(coefficient_series, apex_idx);

            // FR-5.1.10: Peak boundary sharpness
            features.peak_sharpness = self.compute_peak_sharpness(coefficient_series, apex_idx);

            // FR-5.1.11: Peak prominence (apex / baseline)
            features.peak_prominence = self.compute_peak_prominence(coefficient_series, apex_value);

            // FR-5.1.3: EMG fit quality (simplified - use symmetry as proxy)
            // Full EMG fitting would go here, but for now use a heuristic
            features.emg_fit_quality = self.estimate_emg_quality(coefficient_series, apex_idx);
        }

        // Spectral features from apex spectrum (FR-5.2.*)
        if let Some(spectrum) = apex_spectrum {
            // Compute both LibCosine and XCorr scores
            let spectral_score = self.spectral_scorer.xcorr(spectrum, entry);

            // FR-5.2.4: Dot product (LibCosine is our primary spectral score)
            features.dot_product = spectral_score.lib_cosine;

            // FR-5.2.3: Normalized spectral contrast angle
            features.spectral_contrast_angle =
                (spectral_score.lib_cosine.clamp(-1.0, 1.0)).acos().to_degrees();

            // FR-5.2.7: Fragment coverage (fraction detected)
            features.fragment_coverage = spectral_score.fragment_coverage;

            // FR-5.2.12: Explained intensity fraction
            features.explained_intensity = spectral_score.explained_intensity;

            // FR-5.2.2: Hyperscore (XCorr style)
            features.hyperscore = spectral_score.xcorr;

            // FR-5.2.5: Pearson intensity correlation
            features.pearson_correlation = spectral_score.pearson_correlation;

            // FR-5.2.6: Spearman rank correlation
            features.spearman_correlation = spectral_score.spearman_correlation;

            // FR-5.2.8: Sequence coverage (backbone coverage)
            features.sequence_coverage = spectral_score.sequence_coverage;

            // FR-5.2.9: Consecutive ion count (longest b/y run)
            features.consecutive_ions = spectral_score.consecutive_ions;

            // FR-5.2.10: Base peak match rank
            features.base_peak_rank = spectral_score.base_peak_rank;

            // FR-5.2.11: Top-3 match count
            features.top3_matches = spectral_score.top3_matches;
        }

        // Count modifications
        features.modification_count = entry.modifications.len() as u32;

        features
    }

    /// Compute peak symmetry (leading half area / trailing half area)
    fn compute_peak_symmetry(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        if series.len() < 3 || apex_idx == 0 || apex_idx >= series.len() - 1 {
            return 1.0; // Default to symmetric
        }

        // Sum coefficients before and after apex
        let leading_area: f64 = series[..apex_idx].iter().map(|(_, c)| c).sum();
        let trailing_area: f64 = series[apex_idx + 1..].iter().map(|(_, c)| c).sum();

        if trailing_area > 1e-10 {
            leading_area / trailing_area
        } else if leading_area > 1e-10 {
            f64::MAX // Infinitely asymmetric
        } else {
            1.0
        }
    }

    /// Compute coefficient stability (inverse of variance near apex)
    fn compute_coefficient_stability(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        // Take 5 points centered on apex
        let start = apex_idx.saturating_sub(2);
        let end = (apex_idx + 3).min(series.len());

        if end <= start + 1 {
            return 0.0;
        }

        let window: Vec<f64> = series[start..end].iter().map(|(_, c)| *c).collect();
        let mean: f64 = window.iter().sum::<f64>() / window.len() as f64;
        let variance: f64 = window.iter().map(|c| (c - mean).powi(2)).sum::<f64>() / window.len() as f64;

        // Return inverse of coefficient of variation (higher = more stable)
        if mean > 1e-10 {
            mean / (variance.sqrt() + 1e-10)
        } else {
            0.0
        }
    }

    /// Compute peak sharpness (slope at boundaries)
    fn compute_peak_sharpness(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        if series.len() < 3 {
            return 0.0;
        }

        let apex_coef = series[apex_idx].1;

        // Compute average slope from apex to boundaries
        let mut total_slope = 0.0;
        let mut count = 0;

        // Slope going left from apex
        if apex_idx > 0 {
            let dt = series[apex_idx].0 - series[0].0;
            if dt > 0.0 {
                total_slope += (apex_coef - series[0].1) / dt;
                count += 1;
            }
        }

        // Slope going right from apex
        if apex_idx < series.len() - 1 {
            let dt = series[series.len() - 1].0 - series[apex_idx].0;
            if dt > 0.0 {
                total_slope += (apex_coef - series[series.len() - 1].1) / dt;
                count += 1;
            }
        }

        if count > 0 {
            total_slope / count as f64
        } else {
            0.0
        }
    }

    /// Compute peak prominence (apex / local baseline)
    fn compute_peak_prominence(&self, series: &[(f64, f64)], apex_value: f64) -> f64 {
        if series.is_empty() {
            return 0.0;
        }

        // Estimate baseline from minimum values at edges
        let edge_values: Vec<f64> = series
            .iter()
            .take(3)
            .chain(series.iter().rev().take(3))
            .map(|(_, c)| *c)
            .collect();

        let baseline = edge_values.iter().cloned().fold(f64::MAX, f64::min).max(0.0);

        if baseline > 1e-10 {
            apex_value / baseline
        } else if apex_value > 0.0 {
            apex_value / 1e-10 // Very high prominence
        } else {
            0.0
        }
    }

    /// Estimate EMG fit quality (heuristic based on peak shape)
    fn estimate_emg_quality(&self, series: &[(f64, f64)], apex_idx: usize) -> f64 {
        if series.len() < 5 {
            return 0.0;
        }

        // Simple heuristic: measure how Gaussian-like the peak is
        // by comparing expected and actual half-width points
        let apex_coef = series[apex_idx].1;
        let half_max = apex_coef / 2.0;

        // Find actual half-width points
        let mut left_half_idx = apex_idx;
        let mut right_half_idx = apex_idx;

        for i in (0..apex_idx).rev() {
            if series[i].1 < half_max {
                left_half_idx = i;
                break;
            }
        }

        for i in apex_idx + 1..series.len() {
            if series[i].1 < half_max {
                right_half_idx = i;
                break;
            }
        }

        // For a Gaussian, we expect roughly symmetric FWHM
        let left_half = apex_idx - left_half_idx;
        let right_half = right_half_idx - apex_idx;

        if left_half > 0 && right_half > 0 {
            // Ratio closer to 1.0 = more symmetric = better fit
            let ratio = left_half.min(right_half) as f64 / left_half.max(right_half) as f64;
            ratio // Value between 0 and 1
        } else {
            0.5 // Default if we can't measure
        }
    }
}

/// Decoy generator using enzyme-aware sequence reversal
///
/// Following the pyXcorrDIA approach:
/// - Reverse peptide sequence preserving terminal residue based on enzyme
/// - Recalculate fragment m/z values for the reversed sequence
/// - Precursor m/z stays the same (same amino acid composition)
#[derive(Debug, Clone)]
pub struct DecoyGenerator {
    /// Method for generating decoys
    method: DecoyMethod,
    /// Enzyme used for digestion
    enzyme: Enzyme,
    /// Amino acid masses for fragment calculation
    aa_masses: HashMap<char, f64>,
}

/// Decoy generation method
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum DecoyMethod {
    /// Reverse the peptide sequence (preserving terminal residue)
    #[default]
    Reverse,
    /// Shuffle the peptide sequence
    Shuffle,
}

/// Digestion enzyme (affects terminal preservation during reversal)
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Enzyme {
    /// Trypsin: C-terminal cleavage at K/R (preserves C-terminus)
    #[default]
    Trypsin,
    /// Lys-C: C-terminal cleavage at K (preserves C-terminus)
    LysC,
    /// Lys-N: N-terminal cleavage at K (preserves N-terminus)
    LysN,
    /// Asp-N: N-terminal cleavage at D (preserves N-terminus)
    AspN,
    /// No enzyme specificity
    Unspecific,
}

impl Enzyme {
    /// Returns true if this enzyme cleaves at C-terminus (so C-term should be preserved)
    pub fn preserves_c_terminus(&self) -> bool {
        matches!(self, Enzyme::Trypsin | Enzyme::LysC | Enzyme::Unspecific)
    }
}

impl DecoyGenerator {
    /// Create a new decoy generator with default settings (trypsin, reverse)
    pub fn new(method: DecoyMethod) -> Self {
        Self {
            method,
            enzyme: Enzyme::default(),
            aa_masses: Self::standard_aa_masses(),
        }
    }

    /// Create a decoy generator with specified enzyme
    pub fn with_enzyme(method: DecoyMethod, enzyme: Enzyme) -> Self {
        Self {
            method,
            enzyme,
            aa_masses: Self::standard_aa_masses(),
        }
    }

    /// Standard monoisotopic amino acid residue masses
    fn standard_aa_masses() -> HashMap<char, f64> {
        let mut masses = HashMap::new();
        masses.insert('A', 71.037114);
        masses.insert('R', 156.101111);
        masses.insert('N', 114.042927);
        masses.insert('D', 115.026943);
        masses.insert('C', 103.009185);
        masses.insert('E', 129.042593);
        masses.insert('Q', 128.058578);
        masses.insert('G', 57.021464);
        masses.insert('H', 137.058912);
        masses.insert('I', 113.084064);
        masses.insert('L', 113.084064);
        masses.insert('K', 128.094963);
        masses.insert('M', 131.040485);
        masses.insert('F', 147.068414);
        masses.insert('P', 97.052764);
        masses.insert('S', 87.032028);
        masses.insert('T', 101.047679);
        masses.insert('W', 186.079313);
        masses.insert('Y', 163.063329);
        masses.insert('V', 99.068414);
        masses
    }

    /// Generate a decoy from a target library entry
    pub fn generate(&self, target: &LibraryEntry) -> Result<LibraryEntry> {
        let mut decoy = target.clone();
        decoy.is_decoy = true;

        // Generate decoy sequence
        let (decoy_sequence, position_mapping) = match self.method {
            DecoyMethod::Reverse => self.reverse_sequence(&target.sequence),
            DecoyMethod::Shuffle => self.shuffle_sequence(&target.sequence),
        };

        decoy.sequence = decoy_sequence.clone();

        // Update modified sequence notation
        decoy.modified_sequence = format!("DECOY_{}", target.modified_sequence);

        // Remap modifications to new positions
        decoy.modifications = self.remap_modifications(&target.modifications, &position_mapping);

        // Recalculate fragment m/z values for the reversed sequence
        decoy.fragments = self.recalculate_fragments(target, &position_mapping);

        // Precursor m/z stays the same (same amino acid composition)
        // decoy.precursor_mz = target.precursor_mz;

        // Update ID to indicate decoy (set high bit)
        decoy.id = target.id | 0x80000000;

        // Update protein ID to indicate decoy
        decoy.protein_ids = target
            .protein_ids
            .iter()
            .map(|p| format!("DECOY_{}", p))
            .collect();

        Ok(decoy)
    }

    /// Reverse sequence preserving terminal residue based on enzyme
    ///
    /// For trypsin/LysC (C-terminal cleavage): preserve C-terminus only
    /// For Lys-N/Asp-N (N-terminal cleavage): preserve N-terminus only
    ///
    /// Returns (reversed_sequence, position_mapping)
    /// position_mapping[new_pos] = old_pos
    fn reverse_sequence(&self, sequence: &str) -> (String, Vec<usize>) {
        let chars: Vec<char> = sequence.chars().collect();
        let len = chars.len();

        if len <= 2 {
            // Too short to reverse meaningfully
            return (sequence.to_string(), (0..len).collect());
        }

        let mut reversed: Vec<char> = Vec::with_capacity(len);
        let mut position_mapping: Vec<usize> = Vec::with_capacity(len);

        if self.enzyme.preserves_c_terminus() {
            // Preserve C-terminus only (trypsin, Lys-C)
            // Reverse everything except the C-terminal residue
            // Original: P E P T I D E K  (K is the cleavage site)
            // Reversed: E D I T P E P K
            //
            // Reverse positions 0..len-1, keep position len-1 at the end
            for i in (0..len - 1).rev() {
                reversed.push(chars[i]);
                position_mapping.push(i);
            }

            // Keep C-term
            reversed.push(chars[len - 1]);
            position_mapping.push(len - 1);
        } else {
            // Preserve N-terminus only (Lys-N, Asp-N)
            // Reverse everything except the N-terminal residue
            // Original: K P E P T I D E  (K is the cleavage site)
            // Reversed: K E D I T P E P

            // Keep N-term
            reversed.push(chars[0]);
            position_mapping.push(0);

            // Reverse the rest (positions 1..len)
            for i in (1..len).rev() {
                reversed.push(chars[i]);
                position_mapping.push(i);
            }
        }

        (reversed.into_iter().collect(), position_mapping)
    }

    /// Shuffle sequence preserving terminal residue
    fn shuffle_sequence(&self, sequence: &str) -> (String, Vec<usize>) {
        let chars: Vec<char> = sequence.chars().collect();
        let len = chars.len();

        if len <= 2 {
            return (sequence.to_string(), (0..len).collect());
        }

        let mut shuffled: Vec<char> = Vec::with_capacity(len);
        let mut position_mapping: Vec<usize> = Vec::with_capacity(len);

        if self.enzyme.preserves_c_terminus() {
            // Keep C-term only, shuffle the rest (positions 0..len-1)
            let prefix_len = len - 1;
            for i in 0..prefix_len {
                let j = (i + prefix_len / 2 + 1) % prefix_len;
                shuffled.push(chars[j]);
                position_mapping.push(j);
            }

            shuffled.push(chars[len - 1]);
            position_mapping.push(len - 1);
        } else {
            // Keep N-term only, shuffle the rest (positions 1..len)
            shuffled.push(chars[0]);
            position_mapping.push(0);

            let rest_len = len - 1;
            for i in 0..rest_len {
                let j = 1 + (i + rest_len / 2 + 1) % rest_len;
                shuffled.push(chars[j]);
                position_mapping.push(j);
            }
        }

        (shuffled.into_iter().collect(), position_mapping)
    }

    /// Remap modifications to new positions after sequence reversal
    fn remap_modifications(
        &self,
        modifications: &[Modification],
        position_mapping: &[usize],
    ) -> Vec<Modification> {
        // Create reverse mapping: old_pos -> new_pos
        let mut reverse_map: HashMap<usize, usize> = HashMap::new();
        for (new_pos, &old_pos) in position_mapping.iter().enumerate() {
            reverse_map.insert(old_pos, new_pos);
        }

        modifications
            .iter()
            .filter_map(|m| {
                reverse_map.get(&m.position).map(|&new_pos| Modification {
                    position: new_pos,
                    unimod_id: m.unimod_id,
                    mass_delta: m.mass_delta,
                    name: m.name.clone(),
                })
            })
            .collect()
    }

    /// Recalculate fragment m/z values for the reversed sequence
    ///
    /// When reversing a sequence:
    /// - b-ions become y-ions at complementary positions
    /// - y-ions become b-ions at complementary positions
    ///
    /// For a sequence of length N:
    /// - b{i} -> y{N-i} (reversed)
    /// - y{i} -> b{N-i} (reversed)
    fn recalculate_fragments(
        &self,
        target: &LibraryEntry,
        position_mapping: &[usize],
    ) -> Vec<LibraryFragment> {
        let seq_len = target.sequence.len();
        let decoy_chars: Vec<char> = self.get_decoy_sequence_chars(target, position_mapping);

        // Build modification mass map for decoy (by position)
        let mut mod_masses: HashMap<usize, f64> = HashMap::new();
        for m in &target.modifications {
            // Find new position for this modification
            if let Some(new_pos) = position_mapping
                .iter()
                .position(|&old| old == m.position)
            {
                mod_masses.insert(new_pos, m.mass_delta);
            }
        }

        target
            .fragments
            .iter()
            .filter_map(|frag| {
                self.recalculate_single_fragment(frag, seq_len, &decoy_chars, &mod_masses)
            })
            .collect()
    }

    /// Get decoy sequence characters from position mapping
    fn get_decoy_sequence_chars(
        &self,
        target: &LibraryEntry,
        position_mapping: &[usize],
    ) -> Vec<char> {
        let target_chars: Vec<char> = target.sequence.chars().collect();
        position_mapping
            .iter()
            .map(|&old_pos| target_chars[old_pos])
            .collect()
    }

    /// Recalculate a single fragment for the decoy
    fn recalculate_single_fragment(
        &self,
        frag: &LibraryFragment,
        seq_len: usize,
        decoy_chars: &[char],
        mod_masses: &HashMap<usize, f64>,
    ) -> Option<LibraryFragment> {
        let annotation = &frag.annotation;

        // Handle b and y ions - they swap when sequence is reversed
        let (new_ion_type, new_ordinal) = match annotation.ion_type {
            IonType::B => {
                // b{i} covers residues 0..i (N-terminal)
                // In reversed sequence, this becomes y{seq_len - i}
                (IonType::Y, (seq_len - annotation.ordinal as usize) as u8)
            }
            IonType::Y => {
                // y{i} covers residues (seq_len-i)..seq_len (C-terminal)
                // In reversed sequence, this becomes b{seq_len - i}
                (IonType::B, (seq_len - annotation.ordinal as usize) as u8)
            }
            // For other ion types, keep as-is but don't recalculate
            _ => return Some(frag.clone()),
        };

        if new_ordinal == 0 || new_ordinal as usize > seq_len {
            return None;
        }

        // Calculate new m/z based on the decoy sequence
        let new_mz = self.calculate_fragment_mz(
            new_ion_type,
            new_ordinal as usize,
            annotation.charge,
            decoy_chars,
            mod_masses,
            annotation.neutral_loss.as_ref().map(|nl| nl.mass()),
        )?;

        Some(LibraryFragment {
            mz: new_mz,
            relative_intensity: frag.relative_intensity,
            annotation: FragmentAnnotation {
                ion_type: new_ion_type,
                ordinal: new_ordinal,
                charge: annotation.charge,
                neutral_loss: annotation.neutral_loss.clone(),
            },
        })
    }

    /// Calculate fragment m/z from sequence
    fn calculate_fragment_mz(
        &self,
        ion_type: IonType,
        ordinal: usize,
        charge: u8,
        sequence: &[char],
        mod_masses: &HashMap<usize, f64>,
        neutral_loss: Option<f64>,
    ) -> Option<f64> {
        const PROTON_MASS: f64 = 1.007276;
        const H2O_MASS: f64 = 18.010565;

        let seq_len = sequence.len();

        // Sum amino acid masses
        let (start, end) = match ion_type {
            IonType::B => (0, ordinal),           // b-ions: N-terminal fragment
            IonType::Y => (seq_len - ordinal, seq_len), // y-ions: C-terminal fragment
            _ => return None,
        };

        if end > seq_len {
            return None;
        }

        let mut mass = 0.0;

        // Sum residue masses
        for i in start..end {
            let aa = sequence.get(i)?;
            mass += self.aa_masses.get(aa)?;

            // Add modification mass if present
            if let Some(&mod_mass) = mod_masses.get(&i) {
                mass += mod_mass;
            }
        }

        // Add terminal groups and ion-specific masses
        match ion_type {
            IonType::B => {
                // b-ion: [M + H]+ - loses C-terminal OH
                mass += PROTON_MASS;
            }
            IonType::Y => {
                // y-ion: [M + H]+ + H2O - keeps C-terminal OH, gains H on N-terminal
                mass += H2O_MASS + PROTON_MASS;
            }
            _ => {}
        }

        // Subtract neutral loss if present
        if let Some(loss) = neutral_loss {
            mass -= loss;
        }

        // Calculate m/z for charge state
        let mz = (mass + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64;

        Some(mz)
    }

    /// Generate decoys for all entries in a library (legacy method without collision detection)
    pub fn generate_all(&self, targets: &[LibraryEntry]) -> Result<Vec<LibraryEntry>> {
        targets.iter().map(|t| self.generate(t)).collect()
    }

    /// Generate decoys with collision detection (pyXcorrDIA approach)
    ///
    /// This method:
    /// 1. Builds a set of all target sequences for collision detection
    /// 2. For each target, tries reversal first
    /// 3. If reversed sequence collides with a target, uses cycling fallback
    /// 4. If all methods fail, excludes the target-decoy pair
    ///
    /// Returns (valid_targets, decoys, statistics)
    pub fn generate_all_with_collision_detection(
        &self,
        targets: &[LibraryEntry],
    ) -> (Vec<LibraryEntry>, Vec<LibraryEntry>, DecoyGenerationStats) {
        use std::collections::HashSet;

        // Build set of all target sequences for collision detection
        let target_sequences: HashSet<&str> = targets
            .iter()
            .map(|t| t.sequence.as_str())
            .collect();

        // Result type for each parallel task
        enum DecoyResult {
            Reversed(LibraryEntry, LibraryEntry),      // (target, decoy)
            CyclingFallback(LibraryEntry, LibraryEntry), // (target, decoy)
            SkippedNoFragments,
            FailedFragmentCalculation,
            ExcludedNoUniqueDecoy,
        }

        // Process each target in parallel
        let results: Vec<DecoyResult> = targets
            .par_iter()
            .map(|target| {
                // Skip if target has no fragments (can't generate meaningful decoy)
                if target.fragments.is_empty() {
                    return DecoyResult::SkippedNoFragments;
                }

                // Try reversal first
                let (reversed_seq, position_mapping) = self.reverse_sequence(&target.sequence);

                // Check if reversed sequence is valid:
                // 1. Different from original (not palindromic)
                // 2. Not in target database (no collision)
                if reversed_seq != target.sequence && !target_sequences.contains(reversed_seq.as_str()) {
                    // Reversal succeeded - create decoy
                    match self.create_decoy_from_mapping(target, &reversed_seq, &position_mapping) {
                        Ok(decoy) => {
                            return DecoyResult::Reversed(target.clone(), decoy);
                        }
                        Err(_) => {
                            return DecoyResult::FailedFragmentCalculation;
                        }
                    }
                }

                // Reversal failed - try cycling fallback
                let max_retries = target.sequence.len().min(10);

                for cycle_length in 1..=max_retries {
                    let (cycled_seq, cycle_mapping) = self.cycle_sequence(&target.sequence, cycle_length);

                    // Check if cycled sequence is valid
                    if cycled_seq != target.sequence && !target_sequences.contains(cycled_seq.as_str()) {
                        match self.create_decoy_from_mapping(target, &cycled_seq, &cycle_mapping) {
                            Ok(decoy) => {
                                return DecoyResult::CyclingFallback(target.clone(), decoy);
                            }
                            Err(_) => {
                                continue; // Try next cycle length
                            }
                        }
                    }
                }

                // All methods failed - exclude this target-decoy pair
                DecoyResult::ExcludedNoUniqueDecoy
            })
            .collect();

        // Aggregate results
        let mut valid_targets: Vec<LibraryEntry> = Vec::with_capacity(targets.len());
        let mut decoys: Vec<LibraryEntry> = Vec::with_capacity(targets.len());
        let mut stats = DecoyGenerationStats::default();

        for result in results {
            match result {
                DecoyResult::Reversed(target, decoy) => {
                    valid_targets.push(target);
                    decoys.push(decoy);
                    stats.reversed += 1;
                }
                DecoyResult::CyclingFallback(target, decoy) => {
                    valid_targets.push(target);
                    decoys.push(decoy);
                    stats.cycling_fallback += 1;
                }
                DecoyResult::SkippedNoFragments => {
                    stats.skipped_no_fragments += 1;
                }
                DecoyResult::FailedFragmentCalculation => {
                    stats.failed_fragment_calculation += 1;
                }
                DecoyResult::ExcludedNoUniqueDecoy => {
                    stats.excluded_no_unique_decoy += 1;
                }
            }
        }

        // Log statistics
        log::info!("Decoy generation statistics:");
        log::info!("  Reversed: {} ({:.1}%)", stats.reversed,
            100.0 * stats.reversed as f64 / targets.len() as f64);
        log::info!("  Cycling fallback: {} ({:.1}%)", stats.cycling_fallback,
            100.0 * stats.cycling_fallback as f64 / targets.len() as f64);
        log::info!("  Excluded (no unique decoy): {} ({:.1}%)", stats.excluded_no_unique_decoy,
            100.0 * stats.excluded_no_unique_decoy as f64 / targets.len() as f64);
        if stats.skipped_no_fragments > 0 {
            log::info!("  Skipped (no fragments): {}", stats.skipped_no_fragments);
        }
        if stats.failed_fragment_calculation > 0 {
            log::warn!("  Failed fragment calculation: {}", stats.failed_fragment_calculation);
        }

        (valid_targets, decoys, stats)
    }

    /// Cycle sequence by N positions, preserving terminal residue based on enzyme
    ///
    /// For trypsin (C-term preserved): PEPTIDEK with cycle=2 -> TIDEPEPK
    /// Cycling shifts the middle portion while keeping the cleavage site fixed
    fn cycle_sequence(&self, sequence: &str, cycle_length: usize) -> (String, Vec<usize>) {
        let chars: Vec<char> = sequence.chars().collect();
        let len = chars.len();

        if len <= 2 || cycle_length == 0 {
            return (sequence.to_string(), (0..len).collect());
        }

        let mut cycled: Vec<char> = Vec::with_capacity(len);
        let mut position_mapping: Vec<usize> = Vec::with_capacity(len);

        if self.enzyme.preserves_c_terminus() {
            // Preserve C-terminus only (trypsin, Lys-C)
            // Cycle positions 0..len-1, keep position len-1 at the end
            let middle_len = len - 1;
            let effective_cycle = cycle_length % middle_len;

            for i in 0..middle_len {
                let src_idx = (i + effective_cycle) % middle_len;
                cycled.push(chars[src_idx]);
                position_mapping.push(src_idx);
            }

            // Keep C-term
            cycled.push(chars[len - 1]);
            position_mapping.push(len - 1);
        } else {
            // Preserve N-terminus only (Lys-N, Asp-N)
            // Keep N-term, cycle positions 1..len
            cycled.push(chars[0]);
            position_mapping.push(0);

            let middle_len = len - 1;
            let effective_cycle = cycle_length % middle_len;

            for i in 0..middle_len {
                let src_idx = 1 + ((i + effective_cycle) % middle_len);
                cycled.push(chars[src_idx]);
                position_mapping.push(src_idx);
            }
        }

        (cycled.into_iter().collect(), position_mapping)
    }

    /// Create a decoy entry from a target and position mapping
    fn create_decoy_from_mapping(
        &self,
        target: &LibraryEntry,
        decoy_sequence: &str,
        position_mapping: &[usize],
    ) -> Result<LibraryEntry> {
        let mut decoy = target.clone();
        decoy.is_decoy = true;
        decoy.sequence = decoy_sequence.to_string();
        decoy.modified_sequence = format!("DECOY_{}", target.modified_sequence);

        // Remap modifications to new positions
        decoy.modifications = self.remap_modifications(&target.modifications, position_mapping);

        // Recalculate fragment m/z values
        decoy.fragments = self.recalculate_fragments(target, position_mapping);

        // Update ID to indicate decoy (set high bit)
        decoy.id = target.id | 0x80000000;

        // Update protein ID to indicate decoy
        decoy.protein_ids = target
            .protein_ids
            .iter()
            .map(|p| format!("DECOY_{}", p))
            .collect();

        Ok(decoy)
    }
}

/// Statistics from decoy generation with collision detection
#[derive(Debug, Clone, Default)]
pub struct DecoyGenerationStats {
    /// Number of decoys generated using reversal (primary method)
    pub reversed: usize,
    /// Number of decoys generated using cycling fallback
    pub cycling_fallback: usize,
    /// Number of targets excluded because no unique decoy could be generated
    pub excluded_no_unique_decoy: usize,
    /// Number of targets skipped because they had no fragments
    pub skipped_no_fragments: usize,
    /// Number of targets where fragment calculation failed
    pub failed_fragment_calculation: usize,
}

impl DecoyGenerationStats {
    /// Total number of valid target-decoy pairs created
    pub fn total_pairs(&self) -> usize {
        self.reversed + self.cycling_fallback
    }

    /// Total number of targets that were not paired
    pub fn total_excluded(&self) -> usize {
        self.excluded_no_unique_decoy + self.skipped_no_fragments + self.failed_fragment_calculation
    }
}

impl Default for DecoyGenerator {
    fn default() -> Self {
        Self::new(DecoyMethod::Reverse)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::IsolationWindow;

    #[test]
    fn test_feature_extractor() {
        let extractor = FeatureExtractor::new();
        let entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);

        let series = vec![(9.0, 0.1), (10.0, 1.0), (11.0, 0.5)];
        let features = extractor.extract(&entry, &series, None);

        assert!((features.peak_apex - 1.0).abs() < 1e-6);
        assert!((features.peak_area - 1.6).abs() < 1e-6);
        assert_eq!(features.n_contributing_scans, 3);
    }

    #[test]
    fn test_spectral_scorer_lib_cosine() {
        let scorer = SpectralScorer::new();

        // Create a library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 75.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum that exactly matches library
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 500.0],
            intensities: vec![100.0, 50.0, 75.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        // Perfect match should give cosine close to 1.0
        assert!(score.lib_cosine > 0.99, "LibCosine should be ~1.0 for perfect match, got {}", score.lib_cosine);
        assert_eq!(score.n_matched, 3);
        assert!((score.fragment_coverage - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_spectral_scorer_partial_match() {
        let scorer = SpectralScorer::new();

        // Create a library entry with 4 fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 500.0,
                relative_intensity: 75.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 600.0,
                relative_intensity: 25.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum with only 2 matching peaks
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 700.0], // 700.0 doesn't match
            intensities: vec![100.0, 50.0, 80.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        // Should have 2 matches out of 4
        assert_eq!(score.n_matched, 2);
        assert!((score.fragment_coverage - 0.5).abs() < 1e-6);
        // LibCosine should still be reasonable for the matched fragments
        assert!(score.lib_cosine > 0.9, "LibCosine should be high for matched fragments");
    }

    #[test]
    fn test_spectral_scorer_no_match() {
        let scorer = SpectralScorer::new();

        // Create a library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum with no matching peaks
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![800.0, 900.0], // No matches
            intensities: vec![100.0, 50.0],
        };

        let score = scorer.lib_cosine(&spectrum, &entry);

        assert_eq!(score.n_matched, 0);
        assert!((score.fragment_coverage - 0.0).abs() < 1e-6);
        assert!((score.lib_cosine - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_spectral_scorer_xcorr() {
        let scorer = SpectralScorer::new();

        // Create a library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create an observed spectrum
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0],
            intensities: vec![100.0, 50.0],
        };

        let score = scorer.xcorr(&spectrum, &entry);

        // XCorr should be computed
        // The exact value depends on the preprocessing, but it should be non-zero
        // for matching spectra
        assert!(score.xcorr != 0.0 || score.lib_cosine > 0.9, "Should have some score");
    }

    #[test]
    fn test_feature_extractor_with_spectrum() {
        let extractor = FeatureExtractor::new();

        // Create library entry with fragments
        let mut entry = LibraryEntry::new(1, "PEPTIDE".into(), "PEPTIDE".into(), 2, 500.0, 10.0);
        entry.fragments = vec![
            LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            },
            LibraryFragment {
                mz: 400.0,
                relative_intensity: 50.0,
                annotation: FragmentAnnotation::default(),
            },
        ];

        // Create apex spectrum
        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0],
            intensities: vec![100.0, 50.0],
        };

        let series = vec![(9.0, 0.1), (10.0, 1.0), (11.0, 0.5)];
        let features = extractor.extract(&entry, &series, Some(&spectrum));

        // Chromatographic features should be computed
        assert!((features.peak_apex - 1.0).abs() < 1e-6);

        // Spectral features should be computed
        assert!(features.dot_product > 0.9, "Dot product should be high for matching spectrum");
        assert!((features.fragment_coverage - 1.0).abs() < 1e-6, "All fragments should match");
    }

    #[test]
    fn test_decoy_reverse_trypsin() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);

        // PEPTIDEK has 8 chars: P(0) E(1) P(2) T(3) I(4) D(5) E(6) K(7)
        // For trypsin: keep only C-term (K), reverse positions 0-6
        // PEPTIDE (positions 0-6) -> EDITPEP (reversed)
        // Result: E D I T P E P K = EDITPEPK
        let target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0);
        let decoy = generator.generate(&target).unwrap();

        assert!(decoy.is_decoy);
        assert_eq!(decoy.sequence, "EDITPEPK"); // PEPTIDE->EDITPEP, K stays
        assert!(decoy.modified_sequence.starts_with("DECOY_"));
        assert_eq!(decoy.id, 1 | 0x80000000);
    }

    #[test]
    fn test_decoy_reverse_lysn() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::LysN);

        // KPEPTIDE -> should preserve K at start
        // Rest PEPTIDE reversed: EDITPEP
        // Result: KEDITPEP
        let target = LibraryEntry::new(1, "KPEPTIDE".into(), "KPEPTIDE".into(), 2, 500.0, 10.0);
        let decoy = generator.generate(&target).unwrap();

        assert!(decoy.is_decoy);
        assert_eq!(decoy.sequence, "KEDITPEP"); // K stays, PEPTIDE->EDITPEP
    }

    #[test]
    fn test_decoy_preserves_precursor_mz() {
        let generator = DecoyGenerator::new(DecoyMethod::Reverse);
        let target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.123, 10.0);

        let decoy = generator.generate(&target).unwrap();

        // Same amino acids -> same precursor mass
        assert!((decoy.precursor_mz - 500.123).abs() < 1e-6);
    }

    #[test]
    fn test_decoy_protein_id_prefix() {
        let generator = DecoyGenerator::new(DecoyMethod::Reverse);
        let mut target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0);
        target.protein_ids = vec!["P12345".to_string(), "Q67890".to_string()];

        let decoy = generator.generate(&target).unwrap();

        assert_eq!(decoy.protein_ids.len(), 2);
        assert!(decoy.protein_ids[0].starts_with("DECOY_"));
        assert!(decoy.protein_ids[1].starts_with("DECOY_"));
    }

    #[test]
    fn test_decoy_with_modification() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);
        let mut target = LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0);

        // Add modification at position 3 (0-indexed)
        // PEPTIDEK positions: P(0) E(1) P(2) T(3) I(4) D(5) E(6) K(7)
        target.modifications.push(Modification {
            position: 3, // T
            unimod_id: Some(35),
            mass_delta: 15.994915,
            name: Some("Oxidation".to_string()),
        });

        let decoy = generator.generate(&target).unwrap();

        // Reversed: EDITPEPK
        // Position mapping: [6, 5, 4, 3, 2, 1, 0, 7]
        // new_pos 0 <- old_pos 6 (E)
        // new_pos 1 <- old_pos 5 (D)
        // new_pos 2 <- old_pos 4 (I)
        // new_pos 3 <- old_pos 3 (T)  <- modification stays here
        // new_pos 4 <- old_pos 2 (P)
        // new_pos 5 <- old_pos 1 (E)
        // new_pos 6 <- old_pos 0 (P)
        // new_pos 7 <- old_pos 7 (K)
        assert_eq!(decoy.modifications.len(), 1);
        // The T at original position 3 stays at position 3 (symmetric reversal)
        assert_eq!(decoy.modifications[0].position, 3);
    }

    #[test]
    fn test_short_sequence() {
        let generator = DecoyGenerator::new(DecoyMethod::Reverse);
        let target = LibraryEntry::new(1, "PK".into(), "PK".into(), 2, 200.0, 5.0);

        let decoy = generator.generate(&target).unwrap();

        // Too short to meaningfully reverse
        assert_eq!(decoy.sequence, "PK");
    }

    #[test]
    fn test_enzyme_detection() {
        assert!(Enzyme::Trypsin.preserves_c_terminus());
        assert!(Enzyme::LysC.preserves_c_terminus());
        assert!(!Enzyme::LysN.preserves_c_terminus());
        assert!(!Enzyme::AspN.preserves_c_terminus());
    }

    #[test]
    fn test_spectrum_aggregator() {
        let aggregator = SpectrumAggregator::new().with_tolerance_da(0.5);

        // Create test spectra at different RTs
        let spectra = vec![
            Spectrum {
                scan_number: 1,
                retention_time: 9.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0, 400.0],
                intensities: vec![50.0, 25.0],
            },
            Spectrum {
                scan_number: 2,
                retention_time: 10.0, // Apex
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0, 400.0],
                intensities: vec![100.0, 50.0],
            },
            Spectrum {
                scan_number: 3,
                retention_time: 10.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0, 400.0],
                intensities: vec![75.0, 35.0],
            },
        ];

        let aggregated = aggregator.aggregate(&spectra, 10.0);

        // Should have 2 peaks (300.0 and 400.0)
        assert_eq!(aggregated.mzs.len(), 2);

        // Intensities should be summed
        // 300.0: 50 + 100 + 75 = 225
        // 400.0: 25 + 50 + 35 = 110
        let int_300 = aggregated.intensities[0];
        let int_400 = aggregated.intensities[1];
        assert!((int_300 - 225.0).abs() < 1.0, "Expected ~225, got {}", int_300);
        assert!((int_400 - 110.0).abs() < 1.0, "Expected ~110, got {}", int_400);
    }

    #[test]
    fn test_spectrum_aggregator_weighted() {
        let aggregator = SpectrumAggregator::new().with_tolerance_da(0.5);

        // Create test spectra with corresponding coefficients
        let spectra = vec![
            Spectrum {
                scan_number: 1,
                retention_time: 9.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0],
                intensities: vec![100.0],
            },
            Spectrum {
                scan_number: 2,
                retention_time: 10.0, // Apex
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::symmetric(500.0, 12.5),
                mzs: vec![300.0],
                intensities: vec![100.0],
            },
        ];

        // Coefficients: scan 1 has weight 0.2, scan 2 has weight 0.8
        let coefficients = vec![(9.5, 0.2), (10.0, 0.8)];

        let aggregated = aggregator.aggregate_weighted(&spectra, &coefficients, 10.0);

        // Should have 1 peak
        assert_eq!(aggregated.mzs.len(), 1);

        // Weighted intensity: 100 * 0.2 + 100 * 0.8 = 100 (after normalization)
        // Since we normalize weights to sum to 1.0, result should be 100
        let int = aggregated.intensities[0];
        assert!((int - 100.0).abs() < 1.0, "Expected ~100, got {}", int);
    }

    #[test]
    fn test_decoy_collision_detection() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);

        // Create targets where one reversed sequence collides with another target
        // PEPTIDEK reverses to EDITPEPK
        // Include EDITPEPK as a target so there's a collision
        let mut targets = vec![
            LibraryEntry::new(1, "PEPTIDEK".into(), "PEPTIDEK".into(), 2, 500.0, 10.0),
            LibraryEntry::new(2, "EDITPEPK".into(), "EDITPEPK".into(), 2, 500.0, 11.0), // This IS the reversed form!
            LibraryEntry::new(3, "ANOTHERK".into(), "ANOTHERK".into(), 2, 600.0, 12.0),
        ];

        // Add some fragments so they're not skipped
        for t in &mut targets {
            t.fragments = vec![LibraryFragment {
                mz: 300.0,
                relative_intensity: 100.0,
                annotation: FragmentAnnotation::default(),
            }];
        }

        let (valid_targets, decoys, stats) = generator.generate_all_with_collision_detection(&targets);

        // PEPTIDEK can't use reversal (would collide with EDITPEPK)
        // It should either use cycling or be excluded
        // EDITPEPK reverses to PEPTIDEK (also a collision!)
        // ANOTHERK reverses to REHTONAK (no collision)

        // At minimum, ANOTHERK should get a valid decoy
        assert!(stats.reversed >= 1, "At least ANOTHERK should use reversal");

        // Either collisions are resolved via cycling, or peptides are excluded
        assert!(
            stats.cycling_fallback + stats.excluded_no_unique_decoy >= 2,
            "PEPTIDEK and EDITPEPK should use cycling or be excluded"
        );

        // Total pairs created should match valid_targets and decoys
        assert_eq!(valid_targets.len(), decoys.len());
        assert_eq!(valid_targets.len(), stats.total_pairs());

        // Verify no decoy sequence matches any target sequence
        let target_seqs: std::collections::HashSet<&str> = valid_targets.iter().map(|t| t.sequence.as_str()).collect();
        for decoy in &decoys {
            assert!(
                !target_seqs.contains(decoy.sequence.as_str()),
                "Decoy {} should not match any target",
                decoy.sequence
            );
        }
    }

    #[test]
    fn test_cycle_sequence() {
        let generator = DecoyGenerator::with_enzyme(DecoyMethod::Reverse, Enzyme::Trypsin);

        // PEPTIDEK with cycle=1: EPTIDEPK (shift by 1, keep K)
        let (cycled, _) = generator.cycle_sequence("PEPTIDEK", 1);
        assert_eq!(cycled, "EPTIDEPK");

        // PEPTIDEK with cycle=2: PTIDEPEK
        let (cycled, _) = generator.cycle_sequence("PEPTIDEK", 2);
        assert_eq!(cycled, "PTIDEPEK");
    }
}
