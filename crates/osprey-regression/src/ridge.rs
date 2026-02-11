//! Ridge regression solver for spectrum deconvolution
//!
//! Solves the regularized least squares problem:
//!   minimize ‖Ax - b‖² + λ‖x‖²
//!
//! where A is the design matrix, b is the observed spectrum,
//! x is the coefficient vector, and λ is the regularization parameter.

use crate::cd_nnls::{solve_cd_nnls_f64, CdNnlsParams};
use ndarray::{Array1, Array2};
use osprey_core::{OspreyError, Result};

/// Ridge regression solver
#[derive(Debug, Clone)]
pub struct RidgeSolver {
    /// Default regularization parameter
    default_lambda: f64,
    /// Coordinate descent parameters for NNLS
    cd_params: CdNnlsParams,
}

impl RidgeSolver {
    /// Create a new ridge solver with the given default lambda
    pub fn new(default_lambda: f64) -> Self {
        Self {
            default_lambda,
            cd_params: CdNnlsParams::default(),
        }
    }

    /// Create a new ridge solver with custom CD parameters
    pub fn with_params(default_lambda: f64, cd_params: CdNnlsParams) -> Self {
        Self {
            default_lambda,
            cd_params,
        }
    }

    /// Solve the ridge regression problem
    ///
    /// Given design matrix A (m × k) and observed vector b (m × 1),
    /// returns coefficient vector x (k × 1) that minimizes:
    ///   ‖Ax - b‖² + λ‖x‖²
    pub fn solve(
        &self,
        a: &Array2<f64>,
        b: &Array1<f64>,
        lambda: Option<f64>,
    ) -> Result<Array1<f64>> {
        let lambda = lambda.unwrap_or(self.default_lambda);

        // Check dimensions
        if a.nrows() != b.len() {
            return Err(OspreyError::RegressionError(format!(
                "Dimension mismatch: A has {} rows but b has {} elements",
                a.nrows(),
                b.len()
            )));
        }

        let k = a.ncols();

        // Handle empty case
        if k == 0 {
            return Ok(Array1::zeros(0));
        }

        // Compute A'A
        let ata = a.t().dot(a);

        // Add regularization: A'A + λI
        let mut ata_reg = ata.clone();
        for i in 0..k {
            ata_reg[[i, i]] += lambda;
        }

        // Compute A'b
        let atb = a.t().dot(b);

        // Solve via Cholesky decomposition
        // (A'A + λI) x = A'b
        let x = self.solve_cholesky(&ata_reg, &atb)?;

        Ok(x)
    }

    /// Solve using Cholesky decomposition
    ///
    /// Implements in-place Cholesky factorization and forward/back substitution.
    fn solve_cholesky(&self, ata_reg: &Array2<f64>, atb: &Array1<f64>) -> Result<Array1<f64>> {
        let n = ata_reg.nrows();

        // Cholesky decomposition: A = L * L^T
        // We only store the lower triangular part L
        let mut l = Array2::<f64>::zeros((n, n));

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

    /// Solve with non-negativity constraint using coordinate descent
    ///
    /// Returns coefficient vector x ≥ 0 that minimizes:
    ///   ‖Ax - b‖² + λ‖x‖²
    ///
    /// Uses coordinate descent with active set acceleration for O(m*k) per sweep
    /// instead of O(k³) for Cholesky decomposition.
    pub fn solve_nonnegative(
        &self,
        a: &Array2<f64>,
        b: &Array1<f64>,
        lambda: Option<f64>,
    ) -> Result<Array1<f64>> {
        let lambda = lambda.unwrap_or(self.default_lambda);
        solve_cd_nnls_f64(a, b, lambda, &self.cd_params)
    }

    /// Compute the residual norm ‖Ax - b‖²
    pub fn residual(&self, a: &Array2<f64>, x: &Array1<f64>, b: &Array1<f64>) -> f64 {
        let ax = a.dot(x);
        let diff = &ax - b;
        diff.dot(&diff)
    }
}

impl Default for RidgeSolver {
    fn default() -> Self {
        Self {
            default_lambda: 1.0,
            cd_params: CdNnlsParams::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

    /// Verifies that ridge regression with identity matrix and zero lambda recovers the exact solution.
    #[test]
    fn test_ridge_identity() {
        // Simple test: A = I, b = [1, 2, 3], lambda = 0
        // Solution should be x = [1, 2, 3] / (1 + 0) = [1, 2, 3]
        let solver = RidgeSolver::new(0.0);
        let a = Array2::eye(3);
        let b = array![1.0, 2.0, 3.0];

        let x = solver.solve(&a, &b, None).unwrap();

        assert_abs_diff_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(x[1], 2.0, epsilon = 1e-10);
        assert_abs_diff_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    /// Verifies that ridge regression with lambda=1 shrinks coefficients by a factor of 2 for an identity design matrix.
    #[test]
    fn test_ridge_with_regularization() {
        // With regularization, solution should be shrunk toward zero
        let solver = RidgeSolver::new(1.0);
        let a = Array2::eye(3);
        let b = array![1.0, 2.0, 3.0];

        let x = solver.solve(&a, &b, None).unwrap();

        // With lambda = 1, x = b / (1 + 1) = b / 2
        assert_abs_diff_eq!(x[0], 0.5, epsilon = 1e-10);
        assert_abs_diff_eq!(x[1], 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(x[2], 1.5, epsilon = 1e-10);
    }

    /// Verifies that NNLS projects negative unconstrained solutions to zero.
    #[test]
    fn test_nonnegative_simple() {
        // Test that negative solutions are projected to zero
        let solver = RidgeSolver::new(0.0);
        let a = Array2::eye(3);
        let b = array![1.0, -2.0, 3.0];

        let x = solver.solve_nonnegative(&a, &b, None).unwrap();

        assert!(x[0] >= 0.0);
        assert!(x[1] >= 0.0); // Should be 0, not -2
        assert!(x[2] >= 0.0);
    }

    /// Verifies that NNLS produces all non-negative coefficients with a reasonable residual for correlated columns.
    #[test]
    fn test_nonnegative_non_trivial() {
        // Test with a non-identity matrix where the unconstrained solution
        // has negative values and requires proper NNLS iteration
        let solver = RidgeSolver::new(0.1);

        // Design matrix with correlated columns (like overlapping spectra)
        let a = array![
            [1.0, 0.5, 0.2],
            [0.5, 1.0, 0.3],
            [0.2, 0.3, 1.0],
            [0.1, 0.2, 0.8],
        ];
        // Observed vector that will cause negative coefficients in unconstrained solution
        let b = array![0.8, 0.2, 0.9, 0.7];

        let x = solver.solve_nonnegative(&a, &b, None).unwrap();

        // All coefficients must be non-negative
        for (i, &val) in x.iter().enumerate() {
            assert!(val >= 0.0, "Coefficient {} is negative: {}", i, val);
        }

        // Verify the solution is reasonable by checking residual
        let residual = solver.residual(&a, &x, &b);
        assert!(residual < 1.0, "Residual too high: {}", residual);
    }

    /// Verifies that NNLS achieves a lower or equal residual compared to naive clamping of negative values.
    #[test]
    fn test_nonnegative_better_than_naive_projection() {
        // Verify that NNLS produces a better solution than naive projection
        let solver = RidgeSolver::new(0.1);

        let a = array![[1.0, 0.8, 0.1], [0.8, 1.0, 0.2], [0.1, 0.2, 1.0],];
        let b = array![1.0, -0.5, 0.8];

        // Get proper NNLS solution
        let x_nnls = solver.solve_nonnegative(&a, &b, None).unwrap();

        // Compute naive projection: just clamp negatives to zero
        let x_unconstrained = solver.solve(&a, &b, None).unwrap();
        let x_naive = x_unconstrained.mapv(|v| v.max(0.0));

        // Both should be non-negative
        assert!(x_nnls.iter().all(|&v| v >= 0.0));
        assert!(x_naive.iter().all(|&v| v >= 0.0));

        // NNLS solution should have equal or better residual than naive projection
        let residual_nnls = solver.residual(&a, &x_nnls, &b);
        let residual_naive = solver.residual(&a, &x_naive, &b);

        assert!(
            residual_nnls <= residual_naive + 1e-6,
            "NNLS residual {} should be <= naive residual {}",
            residual_nnls,
            residual_naive
        );
    }

    /// Verifies that NNLS returns the same result as unconstrained ridge when all coefficients are already positive.
    #[test]
    fn test_nonnegative_all_positive_unconstrained() {
        // When unconstrained solution is already all positive,
        // NNLS should return the same result
        let solver = RidgeSolver::new(0.1);
        let a = Array2::eye(3);
        let b = array![1.0, 2.0, 3.0];

        let x_unconstrained = solver.solve(&a, &b, None).unwrap();
        let x_nnls = solver.solve_nonnegative(&a, &b, None).unwrap();

        // Should be the same (or very close)
        for i in 0..3 {
            assert_abs_diff_eq!(x_unconstrained[i], x_nnls[i], epsilon = 1e-6);
        }
    }

    /// Verifies that the residual norm is zero when the solution exactly matches the observed vector.
    #[test]
    fn test_residual() {
        let solver = RidgeSolver::new(0.0);
        let a = Array2::eye(3);
        let b = array![1.0, 2.0, 3.0];
        let x = array![1.0, 2.0, 3.0];

        let res = solver.residual(&a, &x, &b);
        assert_abs_diff_eq!(res, 0.0, epsilon = 1e-10);
    }
}
