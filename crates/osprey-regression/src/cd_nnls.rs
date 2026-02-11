//! Coordinate Descent Non-Negative Least Squares (CD-NNLS) solver
//!
//! Solves: minimize ‖Ax - b‖² + λ‖x‖²  subject to x ≥ 0
//!
//! Uses coordinate descent with active set acceleration for O(m×k) per-sweep
//! cost instead of O(m×k²) for forming the normal equations.
//!
//! ## Algorithm
//!
//! For each variable j, the closed-form coordinate update is:
//!
//! ```text
//! x_j^new = max(0, (a_j · (b - A_{-j} x_{-j})) / (‖a_j‖² + λ))
//!         = max(0, (a_j · residual + x_j · ‖a_j‖²) / (‖a_j‖² + λ))
//! ```
//!
//! where residual = b - Ax is maintained incrementally.
//!
//! ## Active Set Acceleration
//!
//! After warmup sweeps, variables at zero are marked inactive and skipped.
//! Periodic reactivation checks ensure convergence to the true optimum.

use ndarray::{Array1, Array2, Zip};
use osprey_core::{OspreyError, Result};

/// Parameters for the CD-NNLS solver
#[derive(Debug, Clone)]
pub struct CdNnlsParams {
    /// Maximum number of full sweeps
    pub max_sweeps: usize,
    /// Convergence tolerance (max absolute coefficient change)
    pub tolerance: f64,
    /// Relative convergence tolerance
    pub relative_tolerance: f64,
    /// Number of warmup sweeps before activating the active set
    pub warmup_sweeps: usize,
    /// How often to check inactive variables for re-activation
    pub reactivation_interval: usize,
}

impl Default for CdNnlsParams {
    fn default() -> Self {
        Self {
            max_sweeps: 50,
            tolerance: 1e-6,
            relative_tolerance: 1e-4,
            warmup_sweeps: 3,
            reactivation_interval: 5,
        }
    }
}

impl CdNnlsParams {
    /// Create params with f32 tolerance
    pub fn for_f32() -> Self {
        Self {
            tolerance: 1e-6,
            ..Default::default()
        }
    }

    /// Create params with f64 tolerance (tighter)
    pub fn for_f64() -> Self {
        Self {
            tolerance: 1e-8,
            ..Default::default()
        }
    }
}

/// Solve NNLS via coordinate descent (f32 version)
///
/// # Arguments
/// * `a` - Design matrix (m bins × k candidates), column-major access
/// * `b` - Observed spectrum (m × 1)
/// * `lambda` - Regularization parameter (L2 penalty)
/// * `params` - Solver parameters
///
/// # Returns
/// Coefficient vector x (k × 1) with x ≥ 0
pub fn solve_cd_nnls_f32(
    a: &Array2<f32>,
    b: &Array1<f32>,
    lambda: f32,
    params: &CdNnlsParams,
) -> Result<Array1<f32>> {
    let m = a.nrows(); // number of bins
    let k = a.ncols(); // number of candidates

    if k == 0 {
        return Ok(Array1::zeros(0));
    }

    if m != b.len() {
        return Err(OspreyError::RegressionError(format!(
            "Dimension mismatch: A has {} rows but b has {} elements",
            m,
            b.len()
        )));
    }

    // Precompute column norms: ‖a_j‖² + λ
    let mut col_norm_sq = Array1::<f32>::zeros(k);
    for j in 0..k {
        let col = a.column(j);
        col_norm_sq[j] = col.dot(&col) + lambda;
    }

    // Initialize
    let mut x = Array1::<f32>::zeros(k);
    let mut residual = b.to_owned(); // r = b - Ax = b (since x=0)
    let mut active = vec![true; k]; // All active initially

    let tol_f32 = params.tolerance as f32;
    let rel_tol_f32 = params.relative_tolerance as f32;

    for sweep in 0..params.max_sweeps {
        let mut max_change: f32 = 0.0;
        let mut max_coeff: f32 = 0.0;

        // Determine which variables to iterate
        let use_active_set = sweep >= params.warmup_sweeps;

        for j in 0..k {
            // Skip inactive variables (unless still in warmup)
            if use_active_set && !active[j] {
                continue;
            }

            let old_xj = x[j];
            let col_j = a.column(j);

            // Compute correlation: a_j · residual
            // Note: residual = b - Ax, so this gives us -∂f/∂x_j (up to scaling)
            let corr = col_j.dot(&residual);

            // Closed-form coordinate update for ridge NNLS:
            // x_j^* = (a_j · residual + old_xj · ‖a_j‖²) / (‖a_j‖² + λ)
            //       = (a_j · (b - A_{-j} x_{-j})) / (‖a_j‖² + λ)
            let col_norm_only = col_norm_sq[j] - lambda; // ‖a_j‖²
            let unconstrained = (corr + old_xj * col_norm_only) / col_norm_sq[j];
            let new_xj = unconstrained.max(0.0);

            // Update if changed
            let delta = new_xj - old_xj;
            if delta.abs() > 1e-12 {
                // Update residual: r -= delta * a_j (axpy operation)
                Zip::from(&mut residual)
                    .and(&col_j)
                    .for_each(|r, &a_val| *r -= delta * a_val);
                x[j] = new_xj;
            }

            let change = delta.abs();
            max_change = max_change.max(change);
            max_coeff = max_coeff.max(new_xj);
        }

        // Update active set after warmup
        if use_active_set {
            for j in 0..k {
                active[j] = x[j] > 0.0;
            }

            // Periodic re-activation check for inactive variables
            let sweeps_since_warmup = sweep - params.warmup_sweeps;
            if sweeps_since_warmup > 0 && sweeps_since_warmup % params.reactivation_interval == 0 {
                for j in 0..k {
                    if !active[j] {
                        let col_j = a.column(j);
                        // Gradient at x[j]=0 is -a_j · residual
                        // If gradient < -tolerance, the variable wants to become positive
                        let grad = -col_j.dot(&residual);
                        if grad < -tol_f32 {
                            active[j] = true; // re-activate
                        }
                    }
                }
            }
        }

        // Check convergence
        if max_change < tol_f32 {
            break;
        }
        if max_coeff > 1e-10 && max_change / max_coeff < rel_tol_f32 {
            break;
        }
    }

    // Final validation sweep: check all variables once more
    // This ensures we haven't missed anything due to active set pruning
    for j in 0..k {
        if !active[j] {
            let col_j = a.column(j);
            let corr = col_j.dot(&residual);
            let unconstrained = corr / col_norm_sq[j]; // old_xj = 0
            if unconstrained > tol_f32 {
                // This variable should be positive - do one update
                let new_xj = unconstrained.max(0.0);
                let delta = new_xj;
                Zip::from(&mut residual)
                    .and(&col_j)
                    .for_each(|r, &a_val| *r -= delta * a_val);
                x[j] = new_xj;
            }
        }
    }

    Ok(x)
}

/// Solve NNLS via coordinate descent (f64 version)
///
/// # Arguments
/// * `a` - Design matrix (m bins × k candidates), column-major access
/// * `b` - Observed spectrum (m × 1)
/// * `lambda` - Regularization parameter (L2 penalty)
/// * `params` - Solver parameters
///
/// # Returns
/// Coefficient vector x (k × 1) with x ≥ 0
pub fn solve_cd_nnls_f64(
    a: &Array2<f64>,
    b: &Array1<f64>,
    lambda: f64,
    params: &CdNnlsParams,
) -> Result<Array1<f64>> {
    let m = a.nrows();
    let k = a.ncols();

    if k == 0 {
        return Ok(Array1::zeros(0));
    }

    if m != b.len() {
        return Err(OspreyError::RegressionError(format!(
            "Dimension mismatch: A has {} rows but b has {} elements",
            m,
            b.len()
        )));
    }

    // Precompute column norms: ‖a_j‖² + λ
    let mut col_norm_sq = Array1::<f64>::zeros(k);
    for j in 0..k {
        let col = a.column(j);
        col_norm_sq[j] = col.dot(&col) + lambda;
    }

    // Initialize
    let mut x = Array1::<f64>::zeros(k);
    let mut residual = b.to_owned();
    let mut active = vec![true; k];

    for sweep in 0..params.max_sweeps {
        let mut max_change: f64 = 0.0;
        let mut max_coeff: f64 = 0.0;

        let use_active_set = sweep >= params.warmup_sweeps;

        for j in 0..k {
            if use_active_set && !active[j] {
                continue;
            }

            let old_xj = x[j];
            let col_j = a.column(j);

            let corr = col_j.dot(&residual);
            let col_norm_only = col_norm_sq[j] - lambda;
            let unconstrained = (corr + old_xj * col_norm_only) / col_norm_sq[j];
            let new_xj = unconstrained.max(0.0);

            let delta = new_xj - old_xj;
            if delta.abs() > 1e-15 {
                Zip::from(&mut residual)
                    .and(&col_j)
                    .for_each(|r, &a_val| *r -= delta * a_val);
                x[j] = new_xj;
            }

            let change = delta.abs();
            max_change = max_change.max(change);
            max_coeff = max_coeff.max(new_xj);
        }

        if use_active_set {
            for j in 0..k {
                active[j] = x[j] > 0.0;
            }

            let sweeps_since_warmup = sweep - params.warmup_sweeps;
            if sweeps_since_warmup > 0 && sweeps_since_warmup % params.reactivation_interval == 0 {
                for j in 0..k {
                    if !active[j] {
                        let col_j = a.column(j);
                        let grad = -col_j.dot(&residual);
                        if grad < -params.tolerance {
                            active[j] = true;
                        }
                    }
                }
            }
        }

        if max_change < params.tolerance {
            break;
        }
        if max_coeff > 1e-12 && max_change / max_coeff < params.relative_tolerance {
            break;
        }
    }

    // Final validation sweep
    for j in 0..k {
        if !active[j] {
            let col_j = a.column(j);
            let corr = col_j.dot(&residual);
            let unconstrained = corr / col_norm_sq[j];
            if unconstrained > params.tolerance {
                let new_xj = unconstrained.max(0.0);
                Zip::from(&mut residual)
                    .and(&col_j)
                    .for_each(|r, &a_val| *r -= new_xj * a_val);
                x[j] = new_xj;
            }
        }
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr1;

    /// Verifies that CD-NNLS with identity matrix and zero lambda recovers the exact solution.
    #[test]
    fn test_cd_identity_matrix() {
        // A = I, b = [1, 2, 3], lambda = 0 -> x = [1, 2, 3]
        let a = Array2::<f32>::eye(3);
        let b = arr1(&[1.0f32, 2.0, 3.0]);
        let params = CdNnlsParams::for_f32();

        let x = solve_cd_nnls_f32(&a, &b, 0.0, &params).unwrap();

        assert!((x[0] - 1.0).abs() < 1e-5, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-5, "x[1] = {}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-5, "x[2] = {}", x[2]);
    }

    /// Verifies that CD-NNLS with lambda=1 shrinks identity matrix coefficients by a factor of 2.
    #[test]
    fn test_cd_with_regularization() {
        // A = I, b = [1, 2, 3], lambda = 1 -> x = [0.5, 1.0, 1.5]
        // Solution: x_j = b_j / (1 + lambda) = b_j / 2
        let a = Array2::<f32>::eye(3);
        let b = arr1(&[1.0f32, 2.0, 3.0]);
        let params = CdNnlsParams::for_f32();

        let x = solve_cd_nnls_f32(&a, &b, 1.0, &params).unwrap();

        assert!((x[0] - 0.5).abs() < 1e-5, "x[0] = {}", x[0]);
        assert!((x[1] - 1.0).abs() < 1e-5, "x[1] = {}", x[1]);
        assert!((x[2] - 1.5).abs() < 1e-5, "x[2] = {}", x[2]);
    }

    /// Verifies that CD-NNLS clamps negative unconstrained solutions to zero while preserving positive ones.
    #[test]
    fn test_cd_nnls_projection() {
        // A = I, b = [1, -2, 3], lambda = 0 -> x = [1, 0, 3]
        // The unconstrained solution would be [1, -2, 3], but x[1] is clamped to 0
        let a = Array2::<f32>::eye(3);
        let b = arr1(&[1.0f32, -2.0, 3.0]);
        let params = CdNnlsParams::for_f32();

        let x = solve_cd_nnls_f32(&a, &b, 0.0, &params).unwrap();

        assert!((x[0] - 1.0).abs() < 1e-5, "x[0] = {}", x[0]);
        assert!(x[1] < 1e-6, "x[1] should be 0, got {}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-5, "x[2] = {}", x[2]);
    }

    /// Verifies that CD-NNLS handles correlated (collinear) columns with non-negative coefficients and reasonable residual.
    #[test]
    fn test_cd_correlated_columns() {
        // Two columns with similar patterns
        // This tests that CD handles collinearity correctly
        // Column-major: first 4 values are col 0, next 4 are col 1
        use ndarray::ShapeBuilder;
        let a = Array2::from_shape_vec(
            (4, 2).f(), // Column-major (Fortran) order
            vec![1.0_f32, 1.0, 0.0, 1.0, 0.9, 0.9, 0.1, 0.1],
        )
        .unwrap();
        // a is: col 0 = [1, 1, 0, 1], col 1 = [0.9, 0.9, 0.1, 0.1]

        let b = arr1(&[1.0f32, 1.0, 0.0, 1.0]);
        let params = CdNnlsParams::for_f32();

        let x = solve_cd_nnls_f32(&a, &b, 0.1, &params).unwrap();

        // Both coefficients should be non-negative
        assert!(x[0] >= 0.0, "x[0] = {}", x[0]);
        assert!(x[1] >= 0.0, "x[1] = {}", x[1]);

        // Check residual is reasonable (correlated columns make this harder)
        let pred = a.dot(&x);
        let residual: f32 = (&b - &pred).mapv(|v| v * v).sum();
        assert!(residual < 2.0, "Residual too large: {}", residual);
    }

    /// Verifies that CD-NNLS with 100 candidates correctly identifies 5 true peptides via active set acceleration.
    #[test]
    fn test_cd_large_k_sparse_solution() {
        // k = 100 candidates, only 5 have signal
        // This tests active set efficiency
        let m = 200;
        let k = 100;

        let mut a = Array2::<f32>::zeros((m, k));
        let mut b = Array1::<f32>::zeros(m);

        // Create 5 "true" peptides with distinct patterns
        let true_indices = [3, 17, 42, 68, 91];
        for (i, &j) in true_indices.iter().enumerate() {
            let start = i * 40;
            let end = start + 30;
            for row in start..end.min(m) {
                a[[row, j]] = 1.0;
            }
            // Add to observed spectrum with coefficient 0.5
            for row in start..end.min(m) {
                b[row] += 0.5;
            }
        }

        // Add noise columns (random patterns, but lower intensity)
        for j in 0..k {
            if !true_indices.contains(&j) {
                let offset = (j * 7) % m;
                for row in offset..(offset + 10).min(m) {
                    a[[row, j]] = 0.3;
                }
            }
        }

        let params = CdNnlsParams::for_f32();
        let x = solve_cd_nnls_f32(&a, &b, 0.1, &params).unwrap();

        // Count non-zero coefficients (with a meaningful threshold)
        let nnz = x.iter().filter(|&&v| v > 0.01).count();

        // Should find the true peptides plus maybe some noise
        // With correlated noise columns, some spurious coefficients are expected
        assert!(
            nnz >= 3,
            "Too few non-zero coefficients: {} (expected at least 3 of the 5 true peptides)",
            nnz
        );

        // The true peptides should have higher coefficients than noise
        let mut sorted_x: Vec<f32> = x.to_vec();
        sorted_x.sort_by(|a, b| b.partial_cmp(a).unwrap());

        // Top 5 coefficients should be significant (true peptides)
        assert!(
            sorted_x[0] > 0.1,
            "Largest coefficient should be significant"
        );
        assert!(sorted_x[4] > 0.05, "5th largest should be significant too");

        // The true peptides should have non-zero coefficients
        for &j in &true_indices {
            assert!(
                x[j] > 0.1,
                "True peptide {} has low coefficient: {}",
                j,
                x[j]
            );
        }
    }

    /// Verifies that the f64 version of CD-NNLS produces the same results as the f32 version with tighter tolerance.
    #[test]
    fn test_cd_f64_version() {
        // Same as f32 test but with f64
        let a = Array2::<f64>::eye(3);
        let b = arr1(&[1.0f64, 2.0, 3.0]);
        let params = CdNnlsParams::for_f64();

        let x = solve_cd_nnls_f64(&a, &b, 0.0, &params).unwrap();

        assert!((x[0] - 1.0).abs() < 1e-10, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-10, "x[1] = {}", x[1]);
        assert!((x[2] - 3.0).abs() < 1e-10, "x[2] = {}", x[2]);
    }

    /// Verifies that CD-NNLS returns an empty coefficient vector when given zero candidates.
    #[test]
    fn test_cd_empty_input() {
        let a = Array2::<f32>::zeros((10, 0));
        let b = Array1::<f32>::zeros(10);
        let params = CdNnlsParams::for_f32();

        let x = solve_cd_nnls_f32(&a, &b, 0.1, &params).unwrap();
        assert_eq!(x.len(), 0);
    }

    /// Verifies that CD-NNLS returns an error when design matrix rows and observation vector length differ.
    #[test]
    fn test_cd_dimension_mismatch() {
        let a = Array2::<f32>::zeros((10, 5));
        let b = Array1::<f32>::zeros(8); // Wrong size
        let params = CdNnlsParams::for_f32();

        let result = solve_cd_nnls_f32(&a, &b, 0.1, &params);
        assert!(result.is_err());
    }
}
