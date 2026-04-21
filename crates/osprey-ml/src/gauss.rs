//! Gauss-Jordan elimination for solution of systems of linear equations
//!
//! LDA requires solving the generalized eigenvalue problem for scatter matrices
//! Sb and Sw. We can actually solve this as the standard eigenvalue problem for
//! the matrix inv(Sw).dot(Sb) - or we solve the linear sysem Sw.dot(x) = Sb,
//! then calculate the eigenvalue for x. This is the approach we take.

//! Originally from Sage (https://github.com/lazear/sage)
//! Copyright (c) 2022 Michael Lazear
//! Licensed under the MIT License

use super::matrix::Matrix;

#[derive(Debug)]
pub struct Gauss {
    pub left: Matrix,
    pub right: Matrix,
}

impl Matrix {
    fn swap_rows(&mut self, i: usize, j: usize) {
        for k in 0..self.cols {
            let tmp = self[(i, k)];
            self[(i, k)] = self[(j, k)];
            self[(j, k)] = tmp;
        }
    }
}

impl Gauss {
    pub fn solve_inner(left: Matrix, right: Matrix, eps: f64) -> Option<Matrix> {
        let mut g = Gauss { left, right };
        g.fill_zero(eps);
        g.echelon();
        g.reduce();
        g.backfill();

        // If `left` is the identity matrix, then `right` contains
        // the solution to the system of equations
        match g.left_solved() {
            true => Some(g.right),
            false => None,
        }
    }

    pub fn solve(left: Matrix, right: Matrix) -> Option<Matrix> {
        let mut eps = 1E-8;
        while eps <= 1.0 {
            if let Some(mat) = Gauss::solve_inner(left.clone(), right.clone(), eps) {
                return Some(mat);
            }
            eps *= 10.0;
        }
        None
    }
    /// This SO answer details how to handle covariance matrices with zeros on
    /// diagonals, which can ruin solving
    /// https://stackoverflow.com/a/35958102
    /// "thus instead of using Sigma = Cov(X) you do Sigma = Cov(X) + eps * I,
    ///  where eps is prefedefined small constant, and I is identity matrix.
    ///  Consequently you never have a zero values on the diagonal,
    ///  and it is easy to prove that for reasonable epsilon, this will be inversible"
    fn fill_zero(&mut self, eps: f64) {
        for i in 0..self.left.cols {
            self.left[(i, i)] += eps;
        }
    }

    // Is `left` an identity matrix, or else contains rows of all zeros?
    fn left_solved(&self) -> bool {
        const TOL: f64 = 1E-8;
        let n = self.left.cols;
        for i in 0..n {
            for j in 0..n {
                let x = self.left[(i, j)];
                if i == j {
                    if (x - 1.0).abs() > TOL && x.abs() > TOL {
                        log::debug!("Finding solution to linear system failed: left side of matrix [{},{}] = {}", i, j, x);
                        return false;
                    }
                } else if x.abs() > TOL {
                    log::debug!("Finding solution to linear system failed: left side of matrix [{},{}] = {}", i, j, x);
                    return false;
                }
            }
        }
        true
    }

    fn echelon(&mut self) {
        const EPS: f64 = 1E-12;
        let (m, n) = self.left.shape();
        let mut h = 0;
        let mut k = 0;

        while h < m && k < n {
            // Partial pivoting: pick the row with the largest-magnitude value
            // in the current pivot column. Using the signed value skips
            // large-magnitude negatives in favor of small positives, which
            // hurts numerical stability.
            let mut max = (h, 0.0f64);
            for i in h..m {
                let abs_val = self.left[(i, k)].abs();
                if abs_val > max.1 {
                    max = (i, abs_val);
                }
            }
            let i = max.0;
            if max.1 < EPS {
                k += 1;
                continue;
            }

            // Swap rows (partial pivoting)
            if h != max.0 {
                self.left.swap_rows(h, i);
                self.right.swap_rows(h, i);
            }

            // Clear rows below pivot row
            for i in h + 1..m {
                let factor = self.left[(i, k)] / self.left[(h, k)];
                self.left[(i, k)] = 0.0;
                for j in k + 1..n {
                    self.left[(i, j)] -= self.left[(h, j)] * factor;
                }
                for j in 0..self.right.cols {
                    self.right[(i, j)] -= self.right[(h, j)] * factor;
                }
            }
            h += 1;
            k += 1;
        }
    }

    // Reduce left matrix to reduced echelon form - diagonal is all ones
    fn reduce(&mut self) {
        for i in (0..self.left.rows).rev() {
            for j in 0..self.left.cols {
                let x = self.left[(i, j)];
                if x == 0.0 {
                    continue;
                }
                for k in j..self.left.cols {
                    self.left[(i, k)] /= x;
                }
                for k in 0..self.right.cols {
                    self.right[(i, k)] /= x;
                }
                break;
            }
        }
    }

    // Solve the upper triangular matrix
    fn backfill(&mut self) {
        for i in (0..self.left.rows).rev() {
            for j in 0..self.left.cols {
                if self.left[(i, j)] == 0.0 {
                    continue;
                }
                for k in 0..i {
                    let factor = self.left[(k, j)] / self.left[(i, j)];
                    for h in 0..self.left.cols {
                        self.left[(k, h)] -= self.left[(i, h)] * factor;
                    }
                    for h in 0..self.right.cols {
                        self.right[(k, h)] -= self.right[(i, h)] * factor;
                    }
                }
                break;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn all_close(a: &[f64], b: &[f64], eps: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < eps)
    }

    /// Regression test for abs-value pivot selection.
    ///
    /// Column 0 contains a large-magnitude negative (-4) and a zero. Signed-
    /// max pivoting picks the zero (because 0 > -4), skips the column, and
    /// produces a permuted-identity result that `left_solved` rejects through
    /// every eps in the ladder. Abs-max pivoting picks -4 and solves exactly.
    #[test]
    fn gauss_negative_pivot() {
        let left = Matrix::new([-4.0, 2.0, 0.0, 3.0], 2, 2);
        let right = Matrix::new([1.0, 6.0], 2, 1);
        // Use solve_inner with eps=0 for an exact answer; the outer `solve`
        // adds a 1e-8 diagonal perturbation that prevents bit-exact matches.
        let solution =
            Gauss::solve_inner(left, right, 0.0).expect("solve_inner should succeed");
        assert_eq!(solution.shape(), (2, 1));
        let scores = [solution[(0, 0)], solution[(1, 0)]];
        assert!(all_close(&scores, &[0.75, 2.0], 1e-12), "got {:?}", scores);
    }

    /// Regression test for rank-deficient input.
    ///
    /// `left` has rank 1 (row 2 = 2 * row 1). With `eps = 0` (no diagonal
    /// perturbation), `solve_inner` must return `None` because echelon
    /// produces a zero row on `left`, and the surviving non-zero off-diagonal
    /// entry fails the `left_solved` tolerance check.
    #[test]
    fn gauss_near_singular_returns_none() {
        let left = Matrix::new([1.0, 2.0, 2.0, 4.0], 2, 2);
        let right = Matrix::new([1.0, 2.0], 2, 1);
        assert!(Gauss::solve_inner(left, right, 0.0).is_none());
    }
}
