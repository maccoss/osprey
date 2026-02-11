//! Design matrix construction for regression
//!
//! Builds the design matrix A where each column represents a candidate
//! peptide's predicted spectrum in binned form.

use crate::Binner;
use ndarray::{Array1, Array2};
use osprey_core::LibraryEntry;

/// Builder for design matrices
#[derive(Debug)]
pub struct DesignMatrixBuilder {
    binner: Binner,
    normalize_columns: bool,
}

impl DesignMatrixBuilder {
    /// Create a new design matrix builder
    pub fn new(binner: Binner) -> Self {
        Self {
            binner,
            normalize_columns: true,
        }
    }

    /// Set whether to normalize columns to unit sum
    pub fn normalize_columns(mut self, normalize: bool) -> Self {
        self.normalize_columns = normalize;
        self
    }

    /// Get the binner
    pub fn binner(&self) -> &Binner {
        &self.binner
    }

    /// Build a design matrix from candidate library entries
    ///
    /// Returns a matrix of shape (n_bins, n_candidates) where each column
    /// is the binned and optionally normalized predicted spectrum.
    pub fn build(&self, candidates: &[&LibraryEntry]) -> Array2<f64> {
        let n_bins = self.binner.n_bins();
        let n_candidates = candidates.len();

        let mut matrix = Array2::<f64>::zeros((n_bins, n_candidates));

        for (col_idx, entry) in candidates.iter().enumerate() {
            let binned = self.binner.bin_library_entry(entry);
            let mut dense = self.binner.to_dense(&binned);

            if self.normalize_columns {
                self.binner.normalize_dense(&mut dense);
            }

            for (row_idx, value) in dense.iter().enumerate() {
                matrix[[row_idx, col_idx]] = *value;
            }
        }

        matrix
    }

    /// Build a design matrix and return candidate indices
    ///
    /// This is useful when you want to track which library entries
    /// correspond to which columns in the matrix.
    pub fn build_with_indices(&self, candidates: &[&LibraryEntry]) -> (Array2<f64>, Vec<u32>) {
        let indices: Vec<u32> = candidates.iter().map(|e| e.id).collect();
        let matrix = self.build(candidates);
        (matrix, indices)
    }

    /// Bin an observed spectrum into a dense vector
    pub fn bin_observed(&self, spectrum: &osprey_core::Spectrum) -> Array1<f64> {
        let binned = self.binner.bin_spectrum(spectrum);
        let dense = self.binner.to_dense(&binned);
        Array1::from_vec(dense)
    }
}

impl Default for DesignMatrixBuilder {
    fn default() -> Self {
        Self::new(Binner::default())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{FragmentAnnotation, LibraryFragment};

    fn make_test_entry(id: u32, fragments: Vec<(f64, f32)>) -> LibraryEntry {
        let mut entry = LibraryEntry::new(
            id,
            "PEPTIDE".to_string(),
            "PEPTIDE".to_string(),
            2,
            500.0,
            10.0,
        );

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

    /// Verifies that the design matrix has correct dimensions and column-normalized entries summing to 1.0.
    #[test]
    fn test_build_design_matrix() {
        let builder = DesignMatrixBuilder::default();

        let entry1 = make_test_entry(1, vec![(300.0, 100.0), (400.0, 200.0)]);
        let entry2 = make_test_entry(2, vec![(350.0, 150.0), (450.0, 250.0)]);

        let candidates: Vec<&LibraryEntry> = vec![&entry1, &entry2];
        let matrix = builder.build(&candidates);

        assert_eq!(matrix.nrows(), builder.binner().n_bins());
        assert_eq!(matrix.ncols(), 2);

        // Columns should be normalized
        let col1_sum: f64 = matrix.column(0).sum();
        let col2_sum: f64 = matrix.column(1).sum();
        assert!((col1_sum - 1.0).abs() < 1e-6);
        assert!((col2_sum - 1.0).abs() < 1e-6);
    }

    /// Verifies that build_with_indices returns the correct library entry IDs alongside the design matrix.
    #[test]
    fn test_build_with_indices() {
        let builder = DesignMatrixBuilder::default();

        let entry1 = make_test_entry(42, vec![(300.0, 100.0)]);
        let entry2 = make_test_entry(99, vec![(400.0, 200.0)]);

        let candidates: Vec<&LibraryEntry> = vec![&entry1, &entry2];
        let (matrix, indices) = builder.build_with_indices(&candidates);

        assert_eq!(indices, vec![42, 99]);
        assert_eq!(matrix.ncols(), 2);
    }
}
