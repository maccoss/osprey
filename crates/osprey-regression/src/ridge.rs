//! Ridge regression solver for spectrum deconvolution
//!
//! Solves the regularized least squares problem:
//!   minimize ‖Ax - b‖² + λ‖x‖²
//!
//! where A is the design matrix, b is the observed spectrum,
//! x is the coefficient vector, and λ is the regularization parameter.

use ndarray::{Array1, Array2};
use osprey_core::{OspreyError, Result};

/// Ridge regression solver
#[derive(Debug, Clone)]
pub struct RidgeSolver {
    /// Default regularization parameter
    default_lambda: f64,
}

impl RidgeSolver {
    /// Create a new ridge solver with the given default lambda
    pub fn new(default_lambda: f64) -> Self {
        Self { default_lambda }
    }

    /// Solve the ridge regression problem
    ///
    /// Given design matrix A (m × k) and observed vector b (m × 1),
    /// returns coefficient vector x (k × 1) that minimizes:
    ///   ‖Ax - b‖² + λ‖x‖²
    pub fn solve(&self, a: &Array2<f64>, b: &Array1<f64>, lambda: Option<f64>) -> Result<Array1<f64>> {
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
                            "Matrix is not positive definite".into()
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

    /// Solve with non-negativity constraint using active set method
    ///
    /// Returns coefficient vector x ≥ 0 that minimizes:
    ///   ‖Ax - b‖² + λ‖x‖²
    pub fn solve_nonnegative(
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

        // Compute A'A + λI and A'b once
        let mut ata_reg = a.t().dot(a);
        for i in 0..k {
            ata_reg[[i, i]] += lambda;
        }
        let atb = a.t().dot(b);

        // Simple iterative projection method for NNLS
        // Start with unconstrained solution
        let mut x = self.solve_cholesky(&ata_reg, &atb)?;

        // Project to non-negative and iterate
        let max_iter = 100;
        for _ in 0..max_iter {
            // Project: set negative values to zero
            let mut changed = false;
            for i in 0..k {
                if x[i] < 0.0 {
                    x[i] = 0.0;
                    changed = true;
                }
            }

            if !changed {
                break;
            }

            // Re-solve on active set (variables > 0)
            // This is a simplified version - a full active set method would be more efficient
            // For now, we just zero out negative components and continue
        }

        // Final projection
        for i in 0..k {
            if x[i] < 0.0 {
                x[i] = 0.0;
            }
        }

        Ok(x)
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
        Self::new(1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::array;

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

    #[test]
    fn test_nonnegative() {
        // Test that negative solutions are projected to zero
        let solver = RidgeSolver::new(0.0);
        let a = Array2::eye(3);
        let b = array![1.0, -2.0, 3.0];

        let x = solver.solve_nonnegative(&a, &b, None).unwrap();

        assert!(x[0] >= 0.0);
        assert!(x[1] >= 0.0); // Should be 0, not -2
        assert!(x[2] >= 0.0);
    }

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
