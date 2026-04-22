//! Linear SVM for binary classification
//!
//! Implements a linear Support Vector Machine using dual coordinate descent
//! for L2-regularized L2-loss SVM (squared hinge loss).
//!
//! Reference: Hsieh et al. (2008) "A Dual Coordinate Descent Method for
//! Large-scale Linear SVM", ICML. This is the algorithm used by Liblinear
//! and sklearn's LinearSVC.
//!
//! Used in Percolator-style semi-supervised scoring for proteomics FDR control.

use super::matrix::Matrix;
use rayon::prelude::*;

/// Linear SVM classifier trained via dual coordinate descent (Liblinear algorithm)
#[derive(Debug, Clone)]
pub struct LinearSvm {
    /// Feature weights (hyperplane normal vector)
    weights: Vec<f64>,
    /// Bias term (intercept)
    bias: f64,
}

impl LinearSvm {
    /// Train a linear SVM using dual coordinate descent for L2-regularized
    /// L2-loss SVM (squared hinge loss).
    ///
    /// This is the algorithm from Hsieh et al. (2008), used by Liblinear and
    /// sklearn's LinearSVC. It converges at O(log(1/ε)) — exponentially faster
    /// than Pegasos SGD's O(1/ε).
    ///
    /// # Arguments
    /// * `features` - Feature matrix (rows = samples, cols = features)
    /// * `labels` - Labels: `true` = decoy (y=-1), `false` = target (y=+1)
    /// * `c` - Cost parameter (higher = less regularization, harder margin)
    /// * `seed` - Random seed for reproducible shuffling
    ///
    /// # Returns
    /// Trained LinearSvm model
    pub fn fit(features: &Matrix, labels: &[bool], c: f64, seed: u64) -> Self {
        assert_eq!(features.rows, labels.len());
        let n = features.rows;
        let p = features.cols;

        if n == 0 || p == 0 {
            return LinearSvm {
                weights: vec![0.0; p],
                bias: 0.0,
            };
        }

        // Convert labels: target (false) → +1, decoy (true) → -1
        let y: Vec<f64> = labels.iter().map(|&d| if d { -1.0 } else { 1.0 }).collect();

        // Dual coordinate descent for L2-regularized L2-loss SVM
        // Dual problem: min_α (1/2) α^T Q̄ α - e^T α, s.t. α_i ≥ 0
        // where Q̄_ij = y_i y_j (x_i · x_j) + δ_ij/(2C)
        // Primal-dual: w = Σ α_i y_i x_i

        let inv_2c = 1.0 / (2.0 * c);

        // Precompute diagonal: D_ii = ||x_i||² + 1.0 (bias feature) + 1/(2C)
        let mut diag: Vec<f64> = Vec::with_capacity(n);
        for i in 0..n {
            let row = features.row_slice(i);
            let norm_sq: f64 = row.iter().map(|v| v * v).sum();
            diag.push(norm_sq + 1.0 + inv_2c);
        }

        // Initialize dual variables and primal weight vector
        // w has p+1 elements: w[0..p] = feature weights, w[p] = bias
        let mut alpha = vec![0.0f64; n];
        let mut w = vec![0.0f64; p + 1];

        // RNG for index permutation
        let mut rng = Xorshift64::new(seed);
        let mut indices: Vec<usize> = (0..n).collect();

        // Convergence: stop when max projected gradient < eps * initial max PG
        let eps = 0.01; // Liblinear default tolerance
        let max_iter = 200; // Cap iterations (linear convergence typically needs 20-50)
        let mut initial_max_pg = f64::NEG_INFINITY;
        let mut last_max_pg = 0.0f64;
        let mut converged = false;

        for iter in 0..max_iter {
            fisher_yates_shuffle(&mut indices, &mut rng);

            let mut max_pg_violation = 0.0f64;

            for &i in &indices {
                let row = features.row_slice(i);

                // w · x_i (augmented: includes bias w[p] * 1.0)
                let wx: f64 = w[..p].iter().zip(row).map(|(wi, xi)| wi * xi).sum::<f64>() + w[p];

                // Gradient: g = y_i * (w · x_i) - 1 + α_i/(2C)
                let g = y[i] * wx - 1.0 + alpha[i] * inv_2c;

                // Projected gradient for convergence check
                let pg = if alpha[i] == 0.0 { g.min(0.0) } else { g };

                max_pg_violation = max_pg_violation.max(pg.abs());

                if pg.abs() > 1e-12 {
                    let alpha_old = alpha[i];
                    alpha[i] = (alpha[i] - g / diag[i]).max(0.0);
                    let d = (alpha[i] - alpha_old) * y[i];

                    // Update w += d * x_i (augmented)
                    for (wj, &xj) in w[..p].iter_mut().zip(row.iter()) {
                        *wj += d * xj;
                    }
                    w[p] += d; // bias feature = 1.0
                }
            }

            // Set initial max PG on first iteration for relative convergence
            if iter == 0 {
                initial_max_pg = max_pg_violation;
                if initial_max_pg <= 0.0 {
                    converged = true;
                    break;
                }
            }

            last_max_pg = max_pg_violation;
            if max_pg_violation < eps * initial_max_pg {
                log::debug!("SVM converged after {} iterations (n={})", iter + 1, n);
                converged = true;
                break;
            }
        }

        if !converged {
            log::debug!(
                "SVM hit max iter ({}) (n={}, C={:.4}, pg={:.2}/{:.2})",
                max_iter,
                n,
                c,
                last_max_pg,
                eps * initial_max_pg
            );
        }

        // Split augmented weight vector into weights and bias
        let bias = w[p];
        w.truncate(p);

        LinearSvm { weights: w, bias }
    }

    /// Compute decision function values: w · x + b for each sample
    ///
    /// Higher values indicate more target-like (positive class)
    pub fn decision_function(&self, features: &Matrix) -> Vec<f64> {
        let scores = features.dotv(&self.weights);
        scores.iter().map(|&s| s + self.bias).collect()
    }

    /// Score a single feature vector: w · x + b
    #[inline]
    pub fn score_single(&self, features: &[f64]) -> f64 {
        let mut score = self.bias;
        for (w, x) in self.weights.iter().zip(features.iter()) {
            score += w * x;
        }
        score
    }

    /// Get the learned feature weights
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Get the bias term
    pub fn bias(&self) -> f64 {
        self.bias
    }
}

/// Select the best C value via internal cross-validation
///
/// For each candidate C, trains SVM on train folds and evaluates on held-out fold.
/// Returns the C value that yields the most targets passing the FDR threshold.
///
/// # Arguments
/// * `features` - Feature matrix for the training set
/// * `labels` - Labels for the training set (true = decoy)
/// * `entry_ids` - Entry IDs for target-decoy pairing
/// * `c_values` - Candidate C values to try
/// * `fold_assignments` - Pre-computed fold assignments for each sample
/// * `n_folds` - Number of folds
/// * `seed` - Random seed
/// * `fdr_threshold` - FDR threshold for counting passing targets
#[allow(clippy::too_many_arguments)]
pub fn grid_search_c(
    features: &Matrix,
    labels: &[bool],
    entry_ids: &[u32],
    c_values: &[f64],
    fold_assignments: &[usize],
    n_folds: usize,
    seed: u64,
    fdr_threshold: f64,
) -> (f64, Vec<usize>) {
    // Evaluate each C value in parallel
    let results: Vec<(f64, usize)> = c_values
        .par_iter()
        .map(|&c| {
            let mut total_passing = 0usize;

            for fold in 0..n_folds {
                // Split into train/test
                let train_indices: Vec<usize> = (0..features.rows)
                    .filter(|&i| fold_assignments[i] != fold)
                    .collect();
                let test_indices: Vec<usize> = (0..features.rows)
                    .filter(|&i| fold_assignments[i] == fold)
                    .collect();

                if train_indices.is_empty() || test_indices.is_empty() {
                    continue;
                }

                // Extract train features and labels
                let train_features = extract_rows(features, &train_indices);
                let train_labels: Vec<bool> = train_indices.iter().map(|&i| labels[i]).collect();

                // Train SVM
                let model = LinearSvm::fit(&train_features, &train_labels, c, seed);

                // Score test set
                let test_features = extract_rows(features, &test_indices);
                let test_scores = model.decision_function(&test_features);
                let test_labels: Vec<bool> = test_indices.iter().map(|&i| labels[i]).collect();
                let test_entry_ids: Vec<u32> = test_indices.iter().map(|&i| entry_ids[i]).collect();

                // Count passing targets via target-decoy competition
                let n_pass = count_passing_targets_svm(
                    &test_scores,
                    &test_labels,
                    &test_entry_ids,
                    fdr_threshold,
                );
                total_passing += n_pass;
            }

            log::debug!(
                "  Grid search C={:.4}: {} passing targets",
                c,
                total_passing
            );

            (c, total_passing)
        })
        .collect();

    // Find best C (highest passing targets, first C as tiebreaker).
    // Iterator::max_by_key returns the LAST element on a tie (per stdlib
    // docs), so we scan manually with a strict `>` to get first-tied,
    // matching the comment above and OspreySharp's GridSearchC. This is
    // a parity fix: a tie between e.g. C=1 and C=10 previously yielded
    // C=10 in Rust but C=1 in C#, and the former's higher-complexity
    // model drove per-fold SVM weight drift between the two
    // implementations on Stellar single-file.
    let per_c_counts: Vec<usize> = results.iter().map(|&(_, count)| count).collect();
    let mut best_c = c_values[0];
    let mut best_passing = 0usize;
    for (c, count) in results {
        if count > best_passing {
            best_passing = count;
            best_c = c;
        }
    }

    log::debug!("  Best C={:.4} ({} passing targets)", best_c, best_passing);
    (best_c, per_c_counts)
}

/// Count targets passing FDR threshold using paired target-decoy competition
fn count_passing_targets_svm(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    fdr_threshold: f64,
) -> usize {
    use std::collections::HashMap;

    // Group by base_id, keep best score per target/decoy
    let mut targets: HashMap<u32, (usize, f64)> = HashMap::new();
    let mut decoys: HashMap<u32, f64> = HashMap::new();

    for (i, (&is_decoy, &eid)) in labels.iter().zip(entry_ids).enumerate() {
        let base_id = eid & 0x7FFFFFFF;
        if is_decoy {
            decoys
                .entry(base_id)
                .and_modify(|s| {
                    if scores[i] > *s {
                        *s = scores[i];
                    }
                })
                .or_insert(scores[i]);
        } else {
            targets
                .entry(base_id)
                .and_modify(|(_, s)| {
                    if scores[i] > *s {
                        *s = scores[i];
                    }
                })
                .or_insert((i, scores[i]));
        }
    }

    // Compete each pair: winner enters ranked list
    // Include base_id for deterministic tiebreaking in sort
    let mut winners: Vec<(f64, bool, u32)> = Vec::with_capacity(targets.len());
    for (&base_id, &(_idx, target_score)) in &targets {
        let decoy_score = decoys.get(&base_id).copied().unwrap_or(f64::NEG_INFINITY);
        if target_score > decoy_score {
            winners.push((target_score, false, base_id)); // target wins
        } else {
            winners.push((decoy_score, true, base_id)); // decoy wins (including ties)
        }
    }
    // Add unpaired decoys
    for (&base_id, &decoy_score) in &decoys {
        if !targets.contains_key(&base_id) {
            winners.push((decoy_score, true, base_id));
        }
    }

    // Sort by score descending, then base_id ascending for deterministic tiebreaking
    winners.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.2.cmp(&b.2)));

    // Walk and compute FDR = (n_decoy + 1) / n_target
    let mut n_target = 0usize;
    let mut n_decoy = 0usize;
    let mut max_passing = 0usize;

    for &(_, is_decoy, _) in &winners {
        if is_decoy {
            n_decoy += 1;
        } else {
            n_target += 1;
        }

        if n_target > 0 {
            let fdr = (n_decoy + 1) as f64 / n_target as f64;
            if fdr <= fdr_threshold {
                max_passing = n_target;
            }
        }
    }

    max_passing
}

/// Extract specific rows from a matrix into a new matrix
fn extract_rows(matrix: &Matrix, row_indices: &[usize]) -> Matrix {
    let n_cols = matrix.cols;
    let data: Vec<f64> = row_indices
        .iter()
        .flat_map(|&row| matrix.row_slice(row).iter().copied())
        .collect();

    Matrix::new(data, row_indices.len(), n_cols)
}

/// Standardize features to zero mean and unit variance
#[derive(Debug, Clone)]
pub struct FeatureStandardizer {
    means: Vec<f64>,
    stds: Vec<f64>,
}

impl FeatureStandardizer {
    /// Per-feature means used for standardization.
    pub fn means(&self) -> &[f64] {
        &self.means
    }

    /// Per-feature standard deviations used for standardization.
    /// Values below 1e-12 are replaced by 1.0 during fit to guard
    /// against zero-variance features.
    pub fn stds(&self) -> &[f64] {
        &self.stds
    }

    /// Compute mean and std for each feature column
    pub fn fit(features: &Matrix) -> Self {
        let n = features.rows as f64;
        let p = features.cols;

        let mut means = vec![0.0; p];
        let mut stds = vec![0.0; p];

        // Compute means
        for row in 0..features.rows {
            for col in 0..p {
                means[col] += features[(row, col)];
            }
        }
        for m in means.iter_mut() {
            *m /= n;
        }

        // Compute standard deviations
        for row in 0..features.rows {
            for col in 0..p {
                let diff = features[(row, col)] - means[col];
                stds[col] += diff * diff;
            }
        }
        for s in stds.iter_mut() {
            *s = (*s / n).sqrt();
            // Avoid division by zero for zero-variance features
            if *s < 1e-12 {
                *s = 1.0;
            }
        }

        FeatureStandardizer { means, stds }
    }

    /// Transform features using pre-computed mean/std: (x - mean) / std
    pub fn transform(&self, features: &Matrix) -> Matrix {
        let mut data = features.data.clone();
        let p = features.cols;

        for row in 0..features.rows {
            for col in 0..p {
                let idx = row * p + col;
                data[idx] = (data[idx] - self.means[col]) / self.stds[col];
            }
        }

        Matrix::new(data, features.rows, p)
    }

    /// Fit and transform in one step
    pub fn fit_transform(features: &Matrix) -> (Self, Matrix) {
        let standardizer = Self::fit(features);
        let transformed = standardizer.transform(features);
        (standardizer, transformed)
    }

    /// Fit from a collection of feature slices (avoids building a full Matrix).
    ///
    /// Each entry in `rows` must have exactly `n_features` elements.
    pub fn fit_from_slices(rows: &[Vec<f64>], n_features: usize) -> Self {
        let n = rows.len() as f64;
        let mut means = vec![0.0; n_features];
        let mut stds = vec![0.0; n_features];

        for row in rows {
            for (col, &val) in row.iter().enumerate() {
                means[col] += val;
            }
        }
        for m in means.iter_mut() {
            *m /= n;
        }

        for row in rows {
            for (col, &val) in row.iter().enumerate() {
                let diff = val - means[col];
                stds[col] += diff * diff;
            }
        }
        for s in stds.iter_mut() {
            *s = (*s / n).sqrt();
            if *s < 1e-12 {
                *s = 1.0;
            }
        }

        FeatureStandardizer { means, stds }
    }

    /// Transform a single feature vector in-place.
    pub fn transform_slice(&self, features: &mut [f64]) {
        for (col, val) in features.iter_mut().enumerate() {
            *val = (*val - self.means[col]) / self.stds[col];
        }
    }

    /// Get the number of features this standardizer was fit on.
    pub fn n_features(&self) -> usize {
        self.means.len()
    }
}

/// Simple xorshift64 PRNG for deterministic shuffling
struct Xorshift64 {
    state: u64,
}

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Ensure non-zero state
        Xorshift64 {
            state: if seed == 0 { 1 } else { seed },
        }
    }

    fn next(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }
}

/// Fisher-Yates shuffle using our PRNG
fn fisher_yates_shuffle(slice: &mut [usize], rng: &mut Xorshift64) {
    for i in (1..slice.len()).rev() {
        let j = (rng.next() as usize) % (i + 1);
        slice.swap(i, j);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linearly_separable() {
        // Targets: high feature values, Decoys: low feature values
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            5.0, 5.0,  // target
            4.0, 6.0,  // target
            6.0, 4.0,  // target
            5.5, 5.5,  // target
            1.0, 1.0,  // decoy
            0.0, 2.0,  // decoy
            2.0, 0.0,  // decoy
            1.5, 1.5,  // decoy
        ], 8, 2);

        let labels = vec![false, false, false, false, true, true, true, true];
        let model = LinearSvm::fit(&features, &labels, 1.0, 42);
        let scores = model.decision_function(&features);

        // All targets should score higher than all decoys
        let target_min = scores[..4].iter().cloned().fold(f64::INFINITY, f64::min);
        let decoy_max = scores[4..]
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            target_min > decoy_max,
            "target_min={}, decoy_max={}",
            target_min,
            decoy_max
        );
    }

    #[test]
    fn test_overlapping_classes() {
        // Some overlap between classes
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            5.0, 5.0,  // target
            4.0, 4.0,  // target
            3.0, 3.0,  // target (near boundary)
            1.0, 1.0,  // decoy
            2.0, 2.0,  // decoy
            2.5, 2.5,  // decoy (near boundary)
        ], 6, 2);

        let labels = vec![false, false, false, true, true, true];
        let model = LinearSvm::fit(&features, &labels, 1.0, 42);
        let scores = model.decision_function(&features);

        // Average target score should be higher than average decoy score
        let avg_target: f64 = scores[..3].iter().sum::<f64>() / 3.0;
        let avg_decoy: f64 = scores[3..].iter().sum::<f64>() / 3.0;
        assert!(
            avg_target > avg_decoy,
            "avg_target={}, avg_decoy={}",
            avg_target,
            avg_decoy
        );
    }

    #[test]
    fn test_weights_direction() {
        // Feature 0 is discriminative, feature 1 is noise
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            10.0, 0.5,  // target
            9.0,  0.3,  // target
            8.0,  0.7,  // target
            1.0,  0.4,  // decoy
            2.0,  0.6,  // decoy
            3.0,  0.2,  // decoy
        ], 6, 2);

        let labels = vec![false, false, false, true, true, true];
        let model = LinearSvm::fit(&features, &labels, 1.0, 42);

        // Weight for feature 0 should be positive (targets have higher values)
        assert!(model.weights()[0] > 0.0, "weights: {:?}", model.weights());
        // Feature 0 (discriminative) should have larger absolute weight than
        // feature 1 (noise), though with only 6 samples the margin may be modest
        assert!(
            model.weights()[0].abs() > model.weights()[1].abs(),
            "weights: {:?}",
            model.weights()
        );
    }

    #[test]
    fn test_empty_input() {
        let features = Matrix::zeros(0, 3);
        let labels: Vec<bool> = vec![];
        let model = LinearSvm::fit(&features, &labels, 1.0, 42);
        assert_eq!(model.weights().len(), 3);
        assert_eq!(model.bias(), 0.0);
    }

    #[test]
    fn test_deterministic_with_seed() {
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            5.0, 5.0,
            1.0, 1.0,
            4.0, 6.0,
            2.0, 0.0,
        ], 4, 2);
        let labels = vec![false, true, false, true];

        let model1 = LinearSvm::fit(&features, &labels, 1.0, 42);
        let model2 = LinearSvm::fit(&features, &labels, 1.0, 42);

        assert_eq!(model1.weights(), model2.weights());
        assert_eq!(model1.bias(), model2.bias());
    }

    #[test]
    fn test_decision_function() {
        let model = LinearSvm {
            weights: vec![1.0, 2.0],
            bias: -3.0,
        };
        let features = Matrix::new(vec![1.0, 1.0, 2.0, 2.0], 2, 2);
        let scores = model.decision_function(&features);

        // score[0] = 1*1 + 2*1 + (-3) = 0
        // score[1] = 1*2 + 2*2 + (-3) = 3
        assert!((scores[0] - 0.0).abs() < 1e-10);
        assert!((scores[1] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_feature_standardizer() {
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            10.0, 100.0,
            20.0, 200.0,
            30.0, 300.0,
        ], 3, 2);

        let (standardizer, transformed) = FeatureStandardizer::fit_transform(&features);

        // Means should be 20.0 and 200.0
        assert!((standardizer.means[0] - 20.0).abs() < 1e-10);
        assert!((standardizer.means[1] - 200.0).abs() < 1e-10);

        // Transformed mean should be ~0
        let col0_mean: f64 = (0..3).map(|r| transformed[(r, 0)]).sum::<f64>() / 3.0;
        let col1_mean: f64 = (0..3).map(|r| transformed[(r, 1)]).sum::<f64>() / 3.0;
        assert!(col0_mean.abs() < 1e-10, "col0_mean = {}", col0_mean);
        assert!(col1_mean.abs() < 1e-10, "col1_mean = {}", col1_mean);

        // Transformed std should be ~1
        let col0_std: f64 = ((0..3)
            .map(|r| (transformed[(r, 0)] - col0_mean).powi(2))
            .sum::<f64>()
            / 3.0)
            .sqrt();
        assert!((col0_std - 1.0).abs() < 1e-10, "col0_std = {}", col0_std);
    }

    #[test]
    fn test_standardizer_zero_variance() {
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            5.0, 1.0,
            5.0, 2.0,
            5.0, 3.0,
        ], 3, 2);

        let standardizer = FeatureStandardizer::fit(&features);

        // Feature 0 has zero variance — std should be set to 1.0
        assert!((standardizer.stds[0] - 1.0).abs() < 1e-10);

        // Should not produce NaN
        let transformed = standardizer.transform(&features);
        for row in 0..3 {
            assert!(
                !transformed[(row, 0)].is_nan(),
                "NaN in zero-variance feature"
            );
        }
    }

    #[test]
    fn test_grid_search_c() {
        // Create well-separated data with entry IDs for competition
        #[rustfmt::skip]
        let features = Matrix::new(vec![
            5.0, 5.0,  // target 1
            4.5, 5.5,  // target 2
            5.5, 4.5,  // target 3
            1.0, 1.0,  // decoy 1
            0.5, 1.5,  // decoy 2
            1.5, 0.5,  // decoy 3
        ], 6, 2);

        let labels = vec![false, false, false, true, true, true];
        let entry_ids = vec![1, 2, 3, 1 | 0x80000000, 2 | 0x80000000, 3 | 0x80000000];
        let fold_assignments = vec![0, 1, 2, 0, 1, 2];
        let c_values = vec![0.01, 0.1, 1.0, 10.0];

        let (best_c, _per_c_counts) = grid_search_c(
            &features,
            &labels,
            &entry_ids,
            &c_values,
            &fold_assignments,
            3,
            42,
            0.10, // 10% FDR for small test
        );

        // Should pick a reasonable C (any of them should work for well-separated data)
        assert!(
            c_values.contains(&best_c),
            "best_c={} not in candidates",
            best_c
        );
    }

    #[test]
    fn test_count_passing_targets() {
        // 3 targets beat their decoys, 1 decoy beats its target
        let scores = vec![0.9, 0.8, 0.7, 0.3, 0.2, 0.1, 0.85];
        let labels = vec![false, false, false, true, true, true, true];
        let entry_ids = vec![
            1,
            2,
            3,
            1 | 0x80000000,
            2 | 0x80000000,
            3 | 0x80000000,
            4 | 0x80000000,
        ];

        let n_pass = count_passing_targets_svm(&scores, &labels, &entry_ids, 0.10);

        // Pairs: (1: 0.9 vs 0.3 → target), (2: 0.8 vs 0.2 → target), (3: 0.7 vs 0.1 → target)
        // Unpaired decoy 4 with score 0.85
        // Winners sorted: target(0.9), decoy_4(0.85), target(0.8), target(0.7)
        // Walk: 1T/0D=0%, 1T/1D=100%, 2T/1D=50%, 3T/1D=33% → none pass 10%
        // Actually with +1: FDR = (0+1)/1=100%, (1+1)/1=200%, (1+1)/2=100%, (1+1)/3=66%
        // Hmm, the +1 makes it harder. For 3 targets with 0 decoys among pairs:
        // Wait, unpaired decoy enters the winners list too.
        // Let me recalculate...
        // Winners: (0.9, target_1), (0.85, decoy_4), (0.8, target_2), (0.7, target_3)
        // Walk: 1T,0D → FDR=(0+1)/1=100%; 1T,1D → FDR=(1+1)/1=200%; 2T,1D → (1+1)/2=100%; 3T,1D → (1+1)/3=66%
        // None pass 10% FDR
        assert_eq!(n_pass, 0);
    }

    #[test]
    fn test_count_passing_targets_clean() {
        // All targets beat their decoys, no unpaired decoys
        let scores = vec![0.9, 0.8, 0.7, 0.6, 0.5, 0.1, 0.05, 0.02, 0.01, 0.005];
        let labels = vec![
            false, false, false, false, false, true, true, true, true, true,
        ];
        let entry_ids: Vec<u32> = vec![
            1,
            2,
            3,
            4,
            5,
            1 | 0x80000000,
            2 | 0x80000000,
            3 | 0x80000000,
            4 | 0x80000000,
            5 | 0x80000000,
        ];

        let n_pass = count_passing_targets_svm(&scores, &labels, &entry_ids, 0.50);

        // All targets win their competitions
        // Winners: 5 targets, 0 decoys
        // FDR = (0+1)/1=100%, (0+1)/2=50%, (0+1)/3=33%, (0+1)/4=25%, (0+1)/5=20%
        // At 50% FDR: max_passing where FDR ≤ 50% → positions 2-5 all pass → 5 targets
        assert_eq!(n_pass, 5);
    }

    #[test]
    fn test_xorshift_deterministic() {
        let mut rng1 = Xorshift64::new(42);
        let mut rng2 = Xorshift64::new(42);

        for _ in 0..100 {
            assert_eq!(rng1.next(), rng2.next());
        }
    }

    #[test]
    fn test_shuffle_deterministic() {
        let mut a = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let mut b = a.clone();

        let mut rng1 = Xorshift64::new(123);
        let mut rng2 = Xorshift64::new(123);

        fisher_yates_shuffle(&mut a, &mut rng1);
        fisher_yates_shuffle(&mut b, &mut rng2);

        assert_eq!(a, b);
        // Should be shuffled (extremely unlikely to stay sorted)
        assert_ne!(a, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    }

    #[test]
    fn test_dual_cd_convergence() {
        // Moderately sized problem: decision boundary at x0 + x1 = 10
        let mut data = Vec::new();
        let mut labels = Vec::new();
        let mut rng = Xorshift64::new(12345);
        let n = 200;
        for _ in 0..n {
            let x0 = (rng.next() % 1000) as f64 / 100.0; // 0-10
            let x1 = (rng.next() % 1000) as f64 / 100.0;
            let is_target = x0 + x1 > 10.0;
            data.push(x0);
            data.push(x1);
            labels.push(!is_target); // false=target, true=decoy
        }
        let features = Matrix::new(data, n, 2);
        let model = LinearSvm::fit(&features, &labels, 1.0, 42);

        // Both weights should be positive (boundary is x0 + x1 = const)
        assert!(model.weights()[0] > 0.0, "w0={}", model.weights()[0]);
        assert!(model.weights()[1] > 0.0, "w1={}", model.weights()[1]);

        // Weights should be similar magnitude (both features equally important)
        let ratio = model.weights()[0] / model.weights()[1];
        assert!(ratio > 0.3 && ratio < 3.0, "weight ratio={:.3}", ratio);
    }

    /// Verifies count_passing_targets_svm() returns deterministic results with tied scores.
    ///
    /// Previously, HashMap iteration order caused different sort orderings for entries
    /// with tied scores, producing different FDR walk results. The fix adds base_id
    /// as a secondary sort key. This test creates tied scores and verifies consistency.
    #[test]
    fn test_count_passing_targets_svm_deterministic() {
        // Create 10 targets and 10 decoys where some targets have tied scores.
        // Targets 1-5 score 0.9 (tied), targets 6-10 score 0.6.
        // Decoys all score 0.1 (all targets win competition).
        let mut scores = Vec::new();
        let mut labels = Vec::new();
        let mut entry_ids = Vec::new();

        // Targets with tied scores
        for id in 1..=5u32 {
            scores.push(0.9);
            labels.push(false);
            entry_ids.push(id);
        }
        for id in 6..=10u32 {
            scores.push(0.6);
            labels.push(false);
            entry_ids.push(id);
        }
        // Decoys
        for id in 1..=10u32 {
            scores.push(0.1);
            labels.push(true);
            entry_ids.push(id | 0x80000000);
        }

        let first_result = count_passing_targets_svm(&scores, &labels, &entry_ids, 0.50);

        // Run multiple times — must always produce the same count
        for _ in 0..20 {
            let result = count_passing_targets_svm(&scores, &labels, &entry_ids, 0.50);
            assert_eq!(
                result, first_result,
                "count_passing_targets_svm must be deterministic with tied scores"
            );
        }

        // With 10 targets winning, 0 decoys in winners, FDR = (0+1)/n_target:
        // After 1 target: (0+1)/1=100%, after 2: 50%, after 3: 33%, ..., after 10: 10%
        // At 50% FDR, max_passing where FDR ≤ 50% → n_target=2 gives 50%, so passes
        assert!(
            first_result >= 2,
            "Expected at least 2 passing, got {}",
            first_result
        );
    }
}
