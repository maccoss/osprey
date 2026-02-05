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

    /// Solve with non-negativity constraint using projected gradient descent
    ///
    /// Returns coefficient vector x ≥ 0 that minimizes:
    ///   ‖Ax - b‖² + λ‖x‖²
    ///
    /// Uses projected gradient descent which is well-suited for dense solutions
    /// (many non-zero coefficients), as expected with unit resolution binning
    /// where many precursor→product ion collisions occur.
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

        // Precompute A'A + λI and A'b (reused every iteration)
        let mut ata_reg = a.t().dot(a);
        for i in 0..k {
            ata_reg[[i, i]] += lambda;
        }
        let atb = a.t().dot(b);

        // Initialize with unconstrained solution
        let mut x = self.solve_cholesky(&ata_reg, &atb)?;

        // Check if unconstrained solution is already non-negative
        let all_positive = x.iter().all(|&v| v >= 0.0);
        if all_positive {
            return Ok(x);
        }

        // Project initial solution to non-negative
        x.mapv_inplace(|v| v.max(0.0));

        // Projected gradient descent parameters
        let max_iter = 1000;
        let tol = 1e-8;

        // Step size: use 1/L where L is the Lipschitz constant of the gradient
        // L = largest eigenvalue of A'A + λI
        // We use Frobenius norm as a conservative upper bound
        let lipschitz = self.estimate_lipschitz(&ata_reg);
        let step_size = 1.0 / lipschitz;

        for _iter in 0..max_iter {
            // Gradient: ∇f(x) = (A'A + λI)x - A'b
            let gradient = ata_reg.dot(&x) - &atb;

            // Gradient step: x_new = x - step_size * gradient
            let x_new: Array1<f64> = &x - step_size * &gradient;

            // Project to non-negative orthant: x_new = max(x_new, 0)
            let x_new = x_new.mapv(|v| v.max(0.0));

            // Check convergence: ||x_new - x||₂ < tol
            let diff = &x_new - &x;
            let change = diff.dot(&diff).sqrt();

            x = x_new;

            if change < tol {
                break;
            }
        }

        Ok(x)
    }

    /// Estimate the Lipschitz constant of the gradient
    ///
    /// Uses Frobenius norm as an upper bound for the spectral norm (largest eigenvalue).
    /// This gives a conservative step size that guarantees convergence.
    fn estimate_lipschitz(&self, ata_reg: &Array2<f64>) -> f64 {
        let mut sum = 0.0;
        for val in ata_reg.iter() {
            sum += val * val;
        }
        // Add small epsilon to avoid division by zero
        sum.sqrt().max(1e-10)
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

    #[test]
    fn test_nonnegative_better_than_naive_projection() {
        // Verify that NNLS produces a better solution than naive projection
        let solver = RidgeSolver::new(0.1);

        let a = array![
            [1.0, 0.8, 0.1],
            [0.8, 1.0, 0.2],
            [0.1, 0.2, 1.0],
        ];
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
