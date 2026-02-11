//! Per-spectrum preprocessing worker for pipelined scoring
//!
//! This module provides preprocessing workers that can process spectra
//! as they stream in, preparing them for XCorr batch scoring.
//!
//! Note: LibCosine scoring uses PPM matching (not binning) and is computed
//! on-demand using LibCosineScorer, not preprocessed vectors.

use osprey_core::{BinConfig, IsolationWindow, Spectrum};

/// Configuration for spectrum preprocessing (XCorr only)
#[derive(Debug, Clone)]
pub struct PreprocessingConfig {
    /// Bin configuration (Comet BIN macro)
    pub bin_config: BinConfig,
}

impl Default for PreprocessingConfig {
    fn default() -> Self {
        Self {
            bin_config: BinConfig::unit_resolution(),
        }
    }
}

impl PreprocessingConfig {
    /// Create config from a BinConfig
    pub fn with_bin_config(bin_config: BinConfig) -> Self {
        Self { bin_config }
    }
}

/// A spectrum that has been preprocessed for XCorr scoring
///
/// Contains XCorr preprocessed vector plus metadata needed for window-based scoring.
/// LibCosine is computed on-demand using LibCosineScorer (PPM matching, no binning).
/// Uses f32 for memory efficiency - sufficient precision for normalized intensity values.
#[derive(Debug, Clone)]
pub struct PreprocessedSpectrum {
    /// Original scan number
    pub scan_number: u32,
    /// Retention time in minutes
    pub retention_time: f64,
    /// Isolation window for grouping
    pub isolation_window: IsolationWindow,
    /// XCorr-preprocessed vector (windowed + flanking subtracted)
    pub xcorr_vector: Vec<f32>,
}

/// Worker that preprocesses spectra for XCorr batch scoring
///
/// Note: This worker only preprocesses for XCorr. LibCosine scoring
/// uses PPM matching via LibCosineScorer and does not use binning.
#[derive(Debug, Clone)]
pub struct PreprocessingWorker {
    pub config: PreprocessingConfig,
}

impl Default for PreprocessingWorker {
    fn default() -> Self {
        Self::new()
    }
}

impl PreprocessingWorker {
    /// Create a new preprocessing worker with default config
    pub fn new() -> Self {
        Self {
            config: PreprocessingConfig::default(),
        }
    }

    /// Create a worker with custom config
    pub fn with_config(config: PreprocessingConfig) -> Self {
        Self { config }
    }

    /// Preprocess a spectrum for XCorr scoring
    ///
    /// This prepares the spectrum for XCorr scoring against library entries.
    /// LibCosine scoring is done separately using LibCosineScorer (PPM matching).
    pub fn process(&self, spectrum: &Spectrum) -> PreprocessedSpectrum {
        let xcorr_vector = self.preprocess_for_xcorr(spectrum);

        PreprocessedSpectrum {
            scan_number: spectrum.scan_number,
            retention_time: spectrum.retention_time,
            isolation_window: spectrum.isolation_window,
            xcorr_vector,
        }
    }

    /// Preprocess spectrum for XCorr scoring (Comet-style)
    /// Uses f32 for memory efficiency - source intensities are already f32.
    fn preprocess_for_xcorr(&self, spectrum: &Spectrum) -> Vec<f32> {
        let n_bins = self.config.bin_config.n_bins;

        // Bin with sqrt transformation using Comet BIN macro
        let mut binned = vec![0.0f32; n_bins];
        for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            if let Some(bin) = self.config.bin_config.mz_to_bin(mz) {
                binned[bin] += intensity.sqrt();
            }
        }

        // Apply windowing normalization
        let windowed = self.apply_windowing_normalization(&binned);

        // Apply flanking bin subtraction
        self.apply_sliding_window(&windowed)
    }

    /// Apply windowing normalization (Comet-style)
    ///
    /// Divides spectrum into 10 windows and normalizes each window's
    /// maximum peak to 50.0, with a 5% threshold.
    fn apply_windowing_normalization(&self, binned: &[f32]) -> Vec<f32> {
        let mut result = binned.to_vec();
        let n_bins = self.config.bin_config.n_bins;
        let num_windows = 10;
        let window_size = n_bins / num_windows;

        for w in 0..num_windows {
            let start = w * window_size;
            let end = ((w + 1) * window_size).min(n_bins);

            // Find max in this window
            let max_val = result[start..end].iter().cloned().fold(0.0f32, f32::max);

            if max_val > 0.0 {
                let threshold = max_val * 0.05;
                let scale = 50.0 / max_val;

                for bin in result[start..end].iter_mut() {
                    if *bin < threshold {
                        *bin = 0.0;
                    } else {
                        *bin *= scale;
                    }
                }
            }
        }

        result
    }

    /// Apply sliding window subtraction for XCorr (Comet-style)
    ///
    /// Uses prefix sum for O(n) performance. Comet divides by (2*offset) = 150.
    fn apply_sliding_window(&self, spectrum: &[f32]) -> Vec<f32> {
        let n = spectrum.len();
        let offset: usize = 75;
        let norm_factor = 1.0f32 / (2 * offset) as f32;

        // Build prefix sum for O(n) window sums
        let mut prefix = vec![0.0f32; n + 1];
        for i in 0..n {
            prefix[i + 1] = prefix[i] + spectrum[i];
        }

        let mut result = vec![0.0f32; n];
        for i in 0..n {
            let left = i.saturating_sub(offset);
            let right = if i + offset < n { i + offset + 1 } else { n };
            let window_sum = prefix[right] - prefix[left];
            let sum_excluding_center = window_sum - spectrum[i];
            result[i] = spectrum[i] - sum_excluding_center * norm_factor;
        }

        result
    }

    /// Get the preprocessing configuration
    pub fn config(&self) -> &PreprocessingConfig {
        &self.config
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_spectrum(scan: u32, rt: f64, peaks: Vec<(f64, f32)>) -> Spectrum {
        Spectrum {
            scan_number: scan,
            retention_time: rt,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: peaks.iter().map(|(mz, _)| *mz).collect(),
            intensities: peaks.iter().map(|(_, int)| *int).collect(),
        }
    }

    /// Verifies that the preprocessing worker preserves scan metadata and produces an XCorr vector with the correct bin count.
    #[test]
    fn test_preprocessing_worker() {
        let worker = PreprocessingWorker::new();
        let spectrum =
            make_test_spectrum(1, 10.0, vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)]);

        let preprocessed = worker.process(&spectrum);

        assert_eq!(preprocessed.scan_number, 1);
        assert!((preprocessed.retention_time - 10.0).abs() < 0.001);
        assert_eq!(
            preprocessed.xcorr_vector.len(),
            worker.config.bin_config.n_bins
        );
    }

    /// Verifies that preprocessing an empty spectrum produces an all-zero XCorr vector without errors.
    #[test]
    fn test_empty_spectrum() {
        let worker = PreprocessingWorker::new();
        let spectrum = make_test_spectrum(1, 10.0, vec![]);

        let preprocessed = worker.process(&spectrum);

        // Should handle empty spectrum gracefully
        assert_eq!(preprocessed.xcorr_vector.iter().sum::<f32>(), 0.0);
    }
}
