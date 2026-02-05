//! Per-spectrum preprocessing worker for pipelined scoring
//!
//! This module provides preprocessing workers that can process spectra
//! as they stream in, preparing them for XCorr batch scoring.
//!
//! Note: LibCosine scoring uses PPM matching (not binning) and is computed
//! on-demand using LibCosineScorer, not preprocessed vectors.

use osprey_core::{IsolationWindow, Spectrum};

/// Configuration for spectrum preprocessing (XCorr only)
#[derive(Debug, Clone)]
pub struct PreprocessingConfig {
    /// Number of bins for XCorr scoring (default: 2000)
    pub num_bins: usize,
    /// Minimum m/z for binning
    pub min_mz: f64,
    /// Bin width for binning
    pub bin_width: f64,
}

impl Default for PreprocessingConfig {
    fn default() -> Self {
        let bin_width = 1.0005079; // Comet default
        let max_mz = 2000.0;
        let num_bins = ((max_mz / bin_width) + 0.6) as usize + 1;
        Self {
            num_bins,
            min_mz: 0.0,
            bin_width,
        }
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
            isolation_window: spectrum.isolation_window.clone(),
            xcorr_vector,
        }
    }

    /// Preprocess spectrum for XCorr scoring (Comet-style)
    /// Uses f32 for memory efficiency - source intensities are already f32.
    fn preprocess_for_xcorr(&self, spectrum: &Spectrum) -> Vec<f32> {
        // Bin with sqrt transformation
        let mut binned = vec![0.0f32; self.config.num_bins];
        for (&mz, &intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            let bin = ((mz / self.config.bin_width) + 0.6) as usize;
            if bin < self.config.num_bins {
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
        let num_windows = 10;
        let window_size = self.config.num_bins / num_windows;

        for w in 0..num_windows {
            let start = w * window_size;
            let end = ((w + 1) * window_size).min(self.config.num_bins);

            // Find max in this window
            let max_val = result[start..end]
                .iter()
                .cloned()
                .fold(0.0f32, f32::max);

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

    /// Apply sliding window subtraction for XCorr
    ///
    /// Subtracts the average of neighboring bins to remove
    /// systematic bias (flanking bin subtraction).
    fn apply_sliding_window(&self, windowed: &[f32]) -> Vec<f32> {
        let offset = 75; // Comet default offset
        let mut result = vec![0.0f32; windowed.len()];

        for i in 0..windowed.len() {
            // Calculate sum of bins in the offset window
            let mut sum = 0.0f32;
            let mut count = 0usize;

            // Left side of window
            if i >= offset {
                for j in (i - offset)..i {
                    sum += windowed[j];
                    count += 1;
                }
            } else {
                for j in 0..i {
                    sum += windowed[j];
                    count += 1;
                }
            }

            // Right side of window
            if i + offset < windowed.len() {
                for j in (i + 1)..=(i + offset) {
                    if j < windowed.len() {
                        sum += windowed[j];
                        count += 1;
                    }
                }
            } else {
                for j in (i + 1)..windowed.len() {
                    sum += windowed[j];
                    count += 1;
                }
            }

            // Subtract average of neighbors from center bin
            let avg = if count > 0 { sum / count as f32 } else { 0.0 };
            result[i] = windowed[i] - avg;
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

    #[test]
    fn test_preprocessing_worker() {
        let worker = PreprocessingWorker::new();
        let spectrum = make_test_spectrum(
            1,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0), (500.0, 75.0)],
        );

        let preprocessed = worker.process(&spectrum);

        assert_eq!(preprocessed.scan_number, 1);
        assert!((preprocessed.retention_time - 10.0).abs() < 0.001);
        assert_eq!(preprocessed.xcorr_vector.len(), worker.config.num_bins);
    }

    #[test]
    fn test_empty_spectrum() {
        let worker = PreprocessingWorker::new();
        let spectrum = make_test_spectrum(1, 10.0, vec![]);

        let preprocessed = worker.process(&spectrum);

        // Should handle empty spectrum gracefully
        assert_eq!(preprocessed.xcorr_vector.iter().sum::<f32>(), 0.0);
    }
}
