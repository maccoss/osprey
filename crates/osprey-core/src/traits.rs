//! Traits for Osprey components
//!
//! This module defines the core traits that enable pluggable implementations
//! for library loading, spectrum sources, and regression solving.

use crate::{LibraryEntry, Result, Spectrum};
use std::path::Path;

/// Trait for loading spectral libraries from various formats
pub trait LibraryLoader: Send + Sync {
    /// Load library entries from a file
    fn load(&self, path: &Path) -> Result<Vec<LibraryEntry>>;

    /// Check if this loader supports the given file format
    fn supports_format(&self, path: &Path) -> bool;

    /// Get the name of the library format
    fn format_name(&self) -> &'static str;
}

/// Trait for sources of MS/MS spectra
///
/// Implementations provide an iterator over spectra from a data file.
pub trait SpectrumSource: Iterator<Item = Result<Spectrum>> + Send {
    /// Get the total number of spectra if known
    fn total_spectra(&self) -> Option<usize>;

    /// Reset the source to the beginning
    fn reset(&mut self) -> Result<()>;

    /// Get the file path being read
    fn file_path(&self) -> &Path;
}

/// Trait for regression solvers
///
/// Implementations solve the regularized regression problem to deconvolute
/// mixed spectra into peptide contributions.
pub trait RegressionSolver: Send + Sync {
    /// Solve the regression problem
    ///
    /// Given design matrix A (m bins × k candidates), observed spectrum b (m × 1),
    /// and regularization parameter lambda, returns coefficient vector x (k × 1).
    ///
    /// The problem solved is: minimize ‖Ax - b‖² + λ‖x‖²
    fn solve(
        &self,
        design_matrix: &[f64], // Flattened row-major m×k matrix
        observed: &[f64],      // m-element observed spectrum
        n_bins: usize,         // m
        n_candidates: usize,   // k
        lambda: f64,
    ) -> Result<Vec<f64>>;

    /// Solve with non-negativity constraint
    ///
    /// Same as solve, but with additional constraint x ≥ 0
    fn solve_nonnegative(
        &self,
        design_matrix: &[f64],
        observed: &[f64],
        n_bins: usize,
        n_candidates: usize,
        lambda: f64,
    ) -> Result<Vec<f64>> {
        // Default implementation: just solve without constraint
        // Subclasses should override with proper NNLS
        self.solve(design_matrix, observed, n_bins, n_candidates, lambda)
    }
}

/// Trait for feature extraction
///
/// Implementations compute features from coefficient time series and spectra.
pub trait FeatureExtractor: Send + Sync {
    /// Extract features for a peptide candidate
    fn extract(
        &self,
        library_entry: &LibraryEntry,
        coefficient_series: &[(f64, f64)], // (rt, coefficient) pairs
        apex_spectrum: Option<&Spectrum>,
    ) -> crate::FeatureSet;
}

/// Trait for peak detection
///
/// Implementations find peaks in coefficient time series.
pub trait PeakDetector: Send + Sync {
    /// Detect peaks in the coefficient time series
    ///
    /// Returns a list of (start_rt, apex_rt, end_rt) tuples for each detected peak.
    fn detect_peaks(&self, coefficient_series: &[(f64, f64)]) -> Vec<(f64, f64, f64)>;
}

/// Trait for FDR control
///
/// Implementations compute q-values from target and decoy scores.
pub trait FdrController: Send + Sync {
    /// Compute q-values for target scores
    ///
    /// Given target scores and decoy scores, returns q-values for each target.
    fn compute_qvalues(&self, target_scores: &[f64], decoy_scores: &[f64]) -> Result<Vec<f64>>;
}

#[cfg(test)]
mod tests {
    use super::*;

    struct MockLoader;

    impl LibraryLoader for MockLoader {
        fn load(&self, _path: &Path) -> Result<Vec<LibraryEntry>> {
            Ok(vec![])
        }

        fn supports_format(&self, path: &Path) -> bool {
            path.extension().map_or(false, |e| e == "mock")
        }

        fn format_name(&self) -> &'static str {
            "Mock"
        }
    }

    /// Verifies that LibraryLoader::supports_format matches only the expected file extension.
    #[test]
    fn test_mock_loader() {
        let loader = MockLoader;
        assert!(loader.supports_format(Path::new("test.mock")));
        assert!(!loader.supports_format(Path::new("test.tsv")));
    }
}
