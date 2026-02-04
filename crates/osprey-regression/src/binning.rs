//! Spectrum binning for regression
//!
//! Converts continuous m/z spectra into discrete bins for matrix operations.
//! Supports both unit resolution (~1 Th bins) and HRAM (variable ppm bins).

use osprey_core::{BinConfig, BinnedSpectrum, LibraryEntry, Spectrum};

/// Spectrum binner that converts continuous m/z values to discrete bins
#[derive(Debug, Clone)]
pub struct Binner {
    config: BinConfig,
}

impl Binner {
    /// Create a new binner with the given configuration
    pub fn new(config: BinConfig) -> Self {
        Self { config }
    }

    /// Create a unit resolution binner (Comet-style ~1 Th bins)
    pub fn unit_resolution() -> Self {
        Self::new(BinConfig::unit_resolution())
    }

    /// Create an HRAM binner (0.02 Th bins, 0 offset)
    pub fn hram() -> Self {
        Self::new(BinConfig::hram())
    }

    /// Get the bin configuration
    pub fn config(&self) -> &BinConfig {
        &self.config
    }

    /// Get the number of bins
    pub fn n_bins(&self) -> usize {
        self.config.n_bins
    }

    /// Bin an observed spectrum
    pub fn bin_spectrum(&self, spectrum: &Spectrum) -> BinnedSpectrum {
        let mut bin_indices = Vec::new();
        let mut intensities = Vec::new();

        for (mz, intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
            if let Some(bin) = self.config.mz_to_bin(*mz) {
                // Check if we already have this bin (accumulate if so)
                if let Some(pos) = bin_indices.iter().position(|&b| b == bin as u32) {
                    intensities[pos] += intensity;
                } else {
                    bin_indices.push(bin as u32);
                    intensities.push(*intensity);
                }
            }
        }

        BinnedSpectrum::new(bin_indices, intensities)
    }

    /// Bin a library entry's predicted spectrum
    pub fn bin_library_entry(&self, entry: &LibraryEntry) -> BinnedSpectrum {
        let mut bin_indices = Vec::new();
        let mut intensities = Vec::new();

        for fragment in &entry.fragments {
            if let Some(bin) = self.config.mz_to_bin(fragment.mz) {
                // Check if we already have this bin (accumulate if so)
                if let Some(pos) = bin_indices.iter().position(|&b| b == bin as u32) {
                    intensities[pos] += fragment.relative_intensity;
                } else {
                    bin_indices.push(bin as u32);
                    intensities.push(fragment.relative_intensity);
                }
            }
        }

        BinnedSpectrum::new(bin_indices, intensities)
    }

    /// Convert a binned spectrum to a dense vector
    pub fn to_dense(&self, binned: &BinnedSpectrum) -> Vec<f64> {
        let mut dense = vec![0.0; self.config.n_bins];
        for (bin, intensity) in binned.bin_indices.iter().zip(binned.intensities.iter()) {
            dense[*bin as usize] = *intensity as f64;
        }
        dense
    }

    /// Normalize a binned spectrum to unit total intensity
    pub fn normalize(&self, binned: &mut BinnedSpectrum) {
        let total: f32 = binned.intensities.iter().sum();
        if total > 0.0 {
            for intensity in &mut binned.intensities {
                *intensity /= total;
            }
        }
    }

    /// Normalize a dense spectrum vector to unit total intensity
    pub fn normalize_dense(&self, dense: &mut [f64]) {
        let total: f64 = dense.iter().sum();
        if total > 0.0 {
            for v in dense.iter_mut() {
                *v /= total;
            }
        }
    }
}

impl Default for Binner {
    fn default() -> Self {
        Self::unit_resolution()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::IsolationWindow;

    #[test]
    fn test_unit_resolution_binner() {
        let binner = Binner::unit_resolution();
        assert!(binner.n_bins() > 1000);
        assert!(binner.n_bins() < 3000);
    }

    #[test]
    fn test_bin_spectrum() {
        let binner = Binner::unit_resolution();

        let spectrum = Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 500.0, 600.0],
            intensities: vec![100.0, 200.0, 300.0, 400.0],
        };

        let binned = binner.bin_spectrum(&spectrum);
        assert_eq!(binned.len(), 4);
    }

    #[test]
    fn test_normalize() {
        let binner = Binner::unit_resolution();
        let mut binned = BinnedSpectrum::new(vec![0, 1, 2], vec![100.0, 200.0, 300.0]);

        binner.normalize(&mut binned);

        let total: f32 = binned.intensities.iter().sum();
        assert!((total - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_to_dense() {
        let binner = Binner::unit_resolution();
        let binned = BinnedSpectrum::new(vec![0, 10, 100], vec![1.0, 2.0, 3.0]);

        let dense = binner.to_dense(&binned);
        assert_eq!(dense.len(), binner.n_bins());
        assert!((dense[0] - 1.0).abs() < 1e-6);
        assert!((dense[10] - 2.0).abs() < 1e-6);
        assert!((dense[100] - 3.0).abs() < 1e-6);
    }
}
