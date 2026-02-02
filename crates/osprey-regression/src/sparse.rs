//! Sparse matrix operations for HRAM data
//!
//! This module provides sparse matrix operations optimized for high-resolution
//! mass spectrometry data where most bins are empty.
//!
//! For HRAM data, instead of binning into a fixed grid, we use direct peak matching
//! with ppm tolerance. This creates a sparse "match matrix" where each column
//! represents a candidate peptide's contribution to matched peaks.

use ndarray::Array1;
use osprey_core::{LibraryEntry, OspreyError, Result, Spectrum};
use sprs::{CsMat, CsMatI, TriMat};

/// Configuration for HRAM peak matching
#[derive(Debug, Clone, Copy)]
pub struct HramConfig {
    /// Mass tolerance in ppm
    pub tolerance_ppm: f64,
    /// Minimum relative intensity threshold (0-1)
    pub min_intensity: f32,
    /// Maximum number of matches per library peak
    pub max_matches_per_peak: usize,
}

impl Default for HramConfig {
    fn default() -> Self {
        Self {
            tolerance_ppm: 20.0,
            min_intensity: 0.01,
            max_matches_per_peak: 3,
        }
    }
}

impl HramConfig {
    /// Create config with given ppm tolerance
    pub fn with_ppm(tolerance_ppm: f64) -> Self {
        Self {
            tolerance_ppm,
            ..Default::default()
        }
    }

    /// Calculate m/z tolerance at a given m/z
    pub fn tolerance_at_mz(&self, mz: f64) -> f64 {
        mz * self.tolerance_ppm / 1_000_000.0
    }

    /// Check if two m/z values match within tolerance
    pub fn matches(&self, mz1: f64, mz2: f64) -> bool {
        let tolerance = self.tolerance_at_mz((mz1 + mz2) / 2.0);
        (mz1 - mz2).abs() <= tolerance
    }
}

/// Sparse design matrix builder for HRAM data
///
/// Instead of binning, this builder directly matches observed peaks to
/// predicted library peaks within a ppm tolerance.
#[derive(Debug)]
pub struct SparseMatrixBuilder {
    config: HramConfig,
    normalize_columns: bool,
}

impl SparseMatrixBuilder {
    /// Create a new sparse matrix builder
    pub fn new(config: HramConfig) -> Self {
        Self {
            config,
            normalize_columns: true,
        }
    }

    /// Set whether to normalize columns
    pub fn normalize_columns(mut self, normalize: bool) -> Self {
        self.normalize_columns = normalize;
        self
    }

    /// Get the HRAM configuration
    pub fn config(&self) -> &HramConfig {
        &self.config
    }

    /// Build a sparse design matrix from observed spectrum and candidates
    ///
    /// Returns a sparse matrix A where:
    /// - Rows correspond to observed peaks (n_observed)
    /// - Columns correspond to candidate peptides (n_candidates)
    /// - A[i,j] = intensity of library peak j that matches observed peak i
    ///
    /// Also returns the observed intensities vector b (n_observed × 1).
    pub fn build(
        &self,
        spectrum: &Spectrum,
        candidates: &[&LibraryEntry],
    ) -> (CsMat<f64>, Array1<f64>, Vec<u32>) {
        let n_observed = spectrum.mzs.len();
        let n_candidates = candidates.len();

        if n_observed == 0 || n_candidates == 0 {
            return (
                CsMat::empty(sprs::CompressedStorage::CSC, n_observed),
                Array1::zeros(n_observed),
                Vec::new(),
            );
        }

        // Build triplet matrix (row, col, value)
        let mut triplets = TriMat::new((n_observed, n_candidates));

        // Track column sums for normalization
        let mut col_sums = vec![0.0f64; n_candidates];

        // For each candidate, find matching observed peaks
        for (col_idx, candidate) in candidates.iter().enumerate() {
            for fragment in &candidate.fragments {
                // Skip low intensity fragments
                if fragment.relative_intensity < self.config.min_intensity {
                    continue;
                }

                // Find matching observed peaks
                let lib_mz = fragment.mz;
                let tolerance = self.config.tolerance_at_mz(lib_mz);

                // Binary search for potential matches
                let start_idx = spectrum
                    .mzs
                    .partition_point(|&mz| mz < lib_mz - tolerance);
                let end_idx = spectrum.mzs.partition_point(|&mz| mz <= lib_mz + tolerance);

                // Add matches (limited to max_matches_per_peak)
                let mut matches_added = 0;
                for row_idx in start_idx..end_idx {
                    if matches_added >= self.config.max_matches_per_peak {
                        break;
                    }

                    let obs_mz = spectrum.mzs[row_idx];
                    if self.config.matches(lib_mz, obs_mz) {
                        let intensity = fragment.relative_intensity as f64;
                        triplets.add_triplet(row_idx, col_idx, intensity);
                        col_sums[col_idx] += intensity;
                        matches_added += 1;
                    }
                }
            }
        }

        // Convert to CSC format
        let mut matrix: CsMatI<f64, usize> = triplets.to_csc();

        // Normalize columns if requested
        if self.normalize_columns {
            // Get column pointers first
            let col_ptrs: Vec<usize> = matrix.indptr().raw_storage().to_vec();

            // Now we can mutably borrow data
            let data = matrix.data_mut();

            for col in 0..n_candidates {
                if col_sums[col] > 0.0 {
                    let start = col_ptrs[col];
                    let end = col_ptrs[col + 1];
                    for idx in start..end {
                        data[idx] /= col_sums[col];
                    }
                }
            }
        }

        // Build observed intensity vector
        let b = Array1::from_vec(spectrum.intensities.iter().map(|&x| x as f64).collect());

        // Get candidate IDs
        let indices: Vec<u32> = candidates.iter().map(|e| e.id).collect();

        (matrix, b, indices)
    }

    /// Build sparse design matrix and return as triple (matrix, observed, indices)
    pub fn build_with_indices(
        &self,
        spectrum: &Spectrum,
        candidates: &[&LibraryEntry],
    ) -> (CsMat<f64>, Array1<f64>, Vec<u32>) {
        self.build(spectrum, candidates)
    }
}

impl Default for SparseMatrixBuilder {
    fn default() -> Self {
        Self::new(HramConfig::default())
    }
}

/// Sparse ridge regression solver
///
/// Solves the regularized least squares problem for sparse design matrices:
///   minimize ‖Ax - b‖² + λ‖x‖²
#[derive(Debug, Clone)]
pub struct SparseRidgeSolver {
    /// Default regularization parameter
    default_lambda: f64,
    /// Maximum iterations for iterative solver
    max_iter: usize,
    /// Convergence tolerance
    tolerance: f64,
}

impl SparseRidgeSolver {
    /// Create a new sparse ridge solver
    pub fn new(default_lambda: f64) -> Self {
        Self {
            default_lambda,
            max_iter: 100,
            tolerance: 1e-8,
        }
    }

    /// Set maximum iterations
    pub fn with_max_iter(mut self, max_iter: usize) -> Self {
        self.max_iter = max_iter;
        self
    }

    /// Set convergence tolerance
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }

    /// Solve the sparse ridge regression problem
    ///
    /// For sparse matrices, we use conjugate gradient on the normal equations:
    ///   (A'A + λI)x = A'b
    pub fn solve(
        &self,
        a: &CsMat<f64>,
        b: &Array1<f64>,
        lambda: Option<f64>,
    ) -> Result<Array1<f64>> {
        let lambda = lambda.unwrap_or(self.default_lambda);

        let n_rows = a.rows();
        let n_cols = a.cols();

        // Check dimensions
        if n_rows != b.len() {
            return Err(OspreyError::RegressionError(format!(
                "Dimension mismatch: A has {} rows but b has {} elements",
                n_rows,
                b.len()
            )));
        }

        if n_cols == 0 {
            return Ok(Array1::zeros(0));
        }

        // Compute A'b
        let atb = self.sparse_transpose_vec_mult(a, b);

        // For small problems, convert to dense and solve directly
        if n_cols <= 200 {
            return self.solve_dense(a, &atb, lambda, n_cols);
        }

        // For larger problems, use conjugate gradient
        self.solve_cg(a, &atb, lambda, n_cols)
    }

    /// Solve using dense matrix operations (for small problems)
    fn solve_dense(
        &self,
        a: &CsMat<f64>,
        atb: &Array1<f64>,
        lambda: f64,
        n_cols: usize,
    ) -> Result<Array1<f64>> {
        // Compute A'A as dense matrix
        let mut ata = ndarray::Array2::<f64>::zeros((n_cols, n_cols));

        // A'A[i,j] = sum_k A[k,i] * A[k,j]
        // For CSC format, iterate over columns
        for col_i in 0..n_cols {
            let col_i_view = a.outer_view(col_i).unwrap();

            for col_j in col_i..n_cols {
                let col_j_view = a.outer_view(col_j).unwrap();

                // Dot product of sparse columns
                let dot = sparse_dot(&col_i_view, &col_j_view);
                ata[[col_i, col_j]] = dot;
                if col_i != col_j {
                    ata[[col_j, col_i]] = dot;
                }
            }
        }

        // Add regularization
        for i in 0..n_cols {
            ata[[i, i]] += lambda;
        }

        // Solve via Cholesky
        solve_cholesky_dense(&ata, atb)
    }

    /// Solve using conjugate gradient (for larger problems)
    fn solve_cg(
        &self,
        a: &CsMat<f64>,
        atb: &Array1<f64>,
        lambda: f64,
        n_cols: usize,
    ) -> Result<Array1<f64>> {
        // Conjugate gradient for (A'A + λI)x = A'b
        let mut x = Array1::<f64>::zeros(n_cols);
        let mut r = atb.clone(); // r = b - Ax, initially r = b
        let mut p = r.clone();
        let mut rs_old = r.dot(&r);

        for _ in 0..self.max_iter {
            // Compute Ap where the "A" here is (A'A + λI)
            let ap = self.apply_normal_equations(a, &p, lambda);

            let alpha = rs_old / p.dot(&ap);

            // x = x + alpha * p
            x = &x + &(&p * alpha);

            // r = r - alpha * Ap
            r = &r - &(&ap * alpha);

            let rs_new = r.dot(&r);

            if rs_new.sqrt() < self.tolerance {
                break;
            }

            // p = r + (rs_new / rs_old) * p
            let beta = rs_new / rs_old;
            p = &r + &(&p * beta);
            rs_old = rs_new;
        }

        Ok(x)
    }

    /// Apply (A'A + λI) to a vector
    fn apply_normal_equations(&self, a: &CsMat<f64>, x: &Array1<f64>, lambda: f64) -> Array1<f64> {
        // First compute Ax
        let ax = self.sparse_vec_mult(a, x);
        // Then compute A'(Ax)
        let atax = self.sparse_transpose_vec_mult(a, &ax);
        // Add λx
        &atax + &(x * lambda)
    }

    /// Multiply sparse matrix by dense vector: y = A * x
    fn sparse_vec_mult(&self, a: &CsMat<f64>, x: &Array1<f64>) -> Array1<f64> {
        let n_rows = a.rows();
        let mut y = Array1::<f64>::zeros(n_rows);

        // For CSC format: iterate over columns
        for (col_idx, col) in a.outer_iterator().enumerate() {
            let x_val = x[col_idx];
            for (row_idx, &val) in col.iter() {
                y[row_idx] += val * x_val;
            }
        }

        y
    }

    /// Multiply sparse matrix transpose by dense vector: y = A' * x
    fn sparse_transpose_vec_mult(&self, a: &CsMat<f64>, x: &Array1<f64>) -> Array1<f64> {
        let n_cols = a.cols();
        let mut y = Array1::<f64>::zeros(n_cols);

        // For CSC format: iterate over columns
        for (col_idx, col) in a.outer_iterator().enumerate() {
            let mut sum = 0.0;
            for (row_idx, &val) in col.iter() {
                sum += val * x[row_idx];
            }
            y[col_idx] = sum;
        }

        y
    }

    /// Solve with non-negativity constraint
    pub fn solve_nonnegative(
        &self,
        a: &CsMat<f64>,
        b: &Array1<f64>,
        lambda: Option<f64>,
    ) -> Result<Array1<f64>> {
        // First solve unconstrained
        let mut x = self.solve(a, b, lambda)?;

        // Simple projection method
        let max_iter = 50;
        for _ in 0..max_iter {
            let mut changed = false;
            for i in 0..x.len() {
                if x[i] < 0.0 {
                    x[i] = 0.0;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }

        Ok(x)
    }

    /// Compute residual norm ‖Ax - b‖²
    pub fn residual(&self, a: &CsMat<f64>, x: &Array1<f64>, b: &Array1<f64>) -> f64 {
        let ax = self.sparse_vec_mult(a, x);
        let diff = &ax - b;
        diff.dot(&diff)
    }
}

impl Default for SparseRidgeSolver {
    fn default() -> Self {
        Self::new(1.0)
    }
}

/// Compute dot product of two sparse vectors
fn sparse_dot(a: &sprs::CsVecViewI<f64, usize>, b: &sprs::CsVecViewI<f64, usize>) -> f64 {
    let mut sum = 0.0;
    let mut a_iter = a.iter();
    let mut b_iter = b.iter();

    let mut a_next = a_iter.next();
    let mut b_next = b_iter.next();

    while let (Some((a_idx, &a_val)), Some((b_idx, &b_val))) = (a_next, b_next) {
        match a_idx.cmp(&b_idx) {
            std::cmp::Ordering::Less => a_next = a_iter.next(),
            std::cmp::Ordering::Greater => b_next = b_iter.next(),
            std::cmp::Ordering::Equal => {
                sum += a_val * b_val;
                a_next = a_iter.next();
                b_next = b_iter.next();
            }
        }
    }

    sum
}

/// Solve dense system via Cholesky decomposition
fn solve_cholesky_dense(ata_reg: &ndarray::Array2<f64>, atb: &Array1<f64>) -> Result<Array1<f64>> {
    let n = ata_reg.nrows();

    // Cholesky decomposition: A = L * L^T
    let mut l = ndarray::Array2::<f64>::zeros((n, n));

    for i in 0..n {
        for j in 0..=i {
            let mut sum = ata_reg[[i, j]];
            for k in 0..j {
                sum -= l[[i, k]] * l[[j, k]];
            }

            if i == j {
                if sum <= 0.0 {
                    return Err(OspreyError::RegressionError(
                        "Matrix is not positive definite".into(),
                    ));
                }
                l[[i, j]] = sum.sqrt();
            } else {
                l[[i, j]] = sum / l[[j, j]];
            }
        }
    }

    // Forward substitution: L * y = b
    let mut y = Array1::<f64>::zeros(n);
    for i in 0..n {
        let mut sum = atb[i];
        for j in 0..i {
            sum -= l[[i, j]] * y[j];
        }
        y[i] = sum / l[[i, i]];
    }

    // Back substitution: L^T * x = y
    let mut x = Array1::<f64>::zeros(n);
    for i in (0..n).rev() {
        let mut sum = y[i];
        for j in (i + 1)..n {
            sum -= l[[j, i]] * x[j];
        }
        x[i] = sum / l[[i, i]];
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{FragmentAnnotation, IsolationWindow, LibraryFragment};

    fn make_test_spectrum() -> Spectrum {
        Spectrum {
            scan_number: 1,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.1, 400.2, 500.3, 600.4],
            intensities: vec![100.0, 200.0, 300.0, 400.0],
        }
    }

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

    #[test]
    fn test_hram_config_tolerance() {
        let config = HramConfig::with_ppm(20.0);

        // At 500 m/z, 20 ppm = 0.01 Th
        let tol = config.tolerance_at_mz(500.0);
        assert!((tol - 0.01).abs() < 1e-6);

        // At 1000 m/z, 20 ppm = 0.02 Th
        let tol = config.tolerance_at_mz(1000.0);
        assert!((tol - 0.02).abs() < 1e-6);
    }

    #[test]
    fn test_hram_config_matches() {
        let config = HramConfig::with_ppm(20.0);

        // Should match within tolerance
        assert!(config.matches(500.0, 500.005)); // 10 ppm
        assert!(config.matches(500.0, 499.995)); // 10 ppm

        // Should not match outside tolerance
        assert!(!config.matches(500.0, 500.02)); // 40 ppm
    }

    #[test]
    fn test_sparse_matrix_builder() {
        let builder = SparseMatrixBuilder::new(HramConfig::with_ppm(50.0));

        let spectrum = make_test_spectrum();
        // Entry with fragments matching spectrum peaks
        let entry1 = make_test_entry(1, vec![(300.1, 1.0), (400.2, 1.0)]);
        // Entry with fragments NOT matching spectrum peaks
        let entry2 = make_test_entry(2, vec![(350.0, 1.0), (450.0, 1.0)]);

        let candidates: Vec<&LibraryEntry> = vec![&entry1, &entry2];
        let (matrix, b, indices) = builder.build(&spectrum, &candidates);

        // Entry1 should have matches, entry2 should not
        assert_eq!(matrix.cols(), 2);
        assert_eq!(indices, vec![1, 2]);
        assert_eq!(b.len(), 4);

        // Column 0 (entry1) should have non-zeros
        let col0_nnz = matrix.outer_view(0).map(|v| v.nnz()).unwrap_or(0);
        assert!(col0_nnz > 0, "Entry1 should have matches");

        // Column 1 (entry2) should have no matches
        let col1_nnz = matrix.outer_view(1).map(|v| v.nnz()).unwrap_or(0);
        assert_eq!(col1_nnz, 0, "Entry2 should have no matches");
    }

    #[test]
    fn test_sparse_solver_simple() {
        let solver = SparseRidgeSolver::new(0.1);

        // Create a simple sparse matrix
        let mut triplets = TriMat::new((3, 2));
        triplets.add_triplet(0, 0, 1.0);
        triplets.add_triplet(1, 0, 0.5);
        triplets.add_triplet(1, 1, 0.5);
        triplets.add_triplet(2, 1, 1.0);

        let a: CsMat<f64> = triplets.to_csc();
        let b = Array1::from_vec(vec![1.0, 0.5, 1.0]);

        let x = solver.solve(&a, &b, None).unwrap();

        assert_eq!(x.len(), 2);
        // Check solution is reasonable (both positive given the setup)
        assert!(x[0] > 0.0);
        assert!(x[1] > 0.0);
    }

    #[test]
    fn test_sparse_solver_nonnegative() {
        let solver = SparseRidgeSolver::new(0.1);

        // Setup where unconstrained solution would be negative
        let mut triplets = TriMat::new((2, 2));
        triplets.add_triplet(0, 0, 1.0);
        triplets.add_triplet(1, 1, 1.0);

        let a: CsMat<f64> = triplets.to_csc();
        let b = Array1::from_vec(vec![1.0, -1.0]);

        let x = solver.solve_nonnegative(&a, &b, Some(0.0)).unwrap();

        // All coefficients should be non-negative
        for &val in x.iter() {
            assert!(val >= 0.0);
        }
    }
}
