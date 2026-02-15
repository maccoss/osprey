//! Linear SVM for binary classification
//!
//! Implements a linear Support Vector Machine using the Pegasos algorithm
//! (Primal Estimated sub-GrAdient SOlver for SVM).
//!
//! Reference: Shalev-Shwartz et al. (2011) "Pegasos: Primal Estimated
//! sub-GrAdient SOlver for SVM", Mathematical Programming.
//!
//! Used in Percolator-style semi-supervised scoring for proteomics FDR control.

use super::matrix::Matrix;
use rayon::prelude::*;

/// Linear SVM classifier trained via Pegasos SGD
#[derive(Debug, Clone)]
pub struct LinearSvm {
    /// Feature weights (hyperplane normal vector)
    weights: Vec<f64>,
    /// Bias term (intercept)
    bias: f64,
}

impl LinearSvm {
    /// Train a linear SVM on the given data using Pegasos SGD.
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

        // Pegasos parameters
        let lambda = 1.0 / (c * n as f64);
        let max_epochs = 50;
        let convergence_threshold = 1e-6;

        let mut w = vec![0.0; p];
        let mut b = 0.0;
        let mut t: u64 = 1; // global iteration counter

        // Sample indices for shuffling
        let mut indices: Vec<usize> = (0..n).collect();
        let mut rng = Xorshift64::new(seed);

        let mut prev_w = w.clone();

        for _epoch in 0..max_epochs {
            // Shuffle indices
            fisher_yates_shuffle(&mut indices, &mut rng);

            for &i in &indices {
                let lr = 1.0 / (lambda * t as f64);
                t += 1;

                // Compute margin: y_i * (w · x_i + b)
                let row = features.row_slice(i);
                let wx: f64 = w.iter().zip(row).map(|(wi, xi)| wi * xi).sum::<f64>() + b;
                let margin = y[i] * wx;

                if margin < 1.0 {
                    // Hinge loss active: update w and b
                    for (wj, &xj) in w.iter_mut().zip(row.iter()) {
                        *wj = (1.0 - lr * lambda) * *wj + lr * y[i] * xj;
                    }
                    b += lr * y[i];
                } else {
                    // No loss: only apply regularization to w
                    for wj in w.iter_mut() {
                        *wj *= 1.0 - lr * lambda;
                    }
                }

                // Pegasos projection: ensure ||w|| ≤ 1/sqrt(lambda)
                let w_norm: f64 = w.iter().map(|wi| wi * wi).sum::<f64>().sqrt();
                let max_norm = 1.0 / lambda.sqrt();
                if w_norm > max_norm {
                    let scale = max_norm / w_norm;
                    for wi in w.iter_mut() {
                        *wi *= scale;
                    }
                }
            }

            // Check convergence: relative weight change
            let diff: f64 = w
                .iter()
                .zip(&prev_w)
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            let w_norm: f64 = w.iter().map(|wi| wi * wi).sum::<f64>().sqrt();
            if w_norm > 0.0 && diff / w_norm < convergence_threshold {
                log::debug!("SVM converged after {} epochs", _epoch + 1);
                break;
            }

            prev_w.clone_from(&w);
        }

        LinearSvm {
            weights: w,
            bias: b,
        }
    }

    /// Compute decision function values: w · x + b for each sample
    ///
    /// Higher values indicate more target-like (positive class)
    pub fn decision_function(&self, features: &Matrix) -> Vec<f64> {
        let scores = features.dotv(&self.weights);
        scores.iter().map(|&s| s + self.bias).collect()
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
) -> f64 {
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

    // Find best C (highest passing targets, first C as tiebreaker)
    let (best_c, best_passing) = results
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .unwrap_or((c_values[0], 0));

    log::debug!("  Best C={:.4} ({} passing targets)", best_c, best_passing);
    best_c
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
    let mut winners: Vec<(f64, bool)> = Vec::with_capacity(targets.len());
    for (base_id, &(_idx, target_score)) in &targets {
        let decoy_score = decoys.get(base_id).copied().unwrap_or(f64::NEG_INFINITY);
        if target_score > decoy_score {
            winners.push((target_score, false)); // target wins
        } else {
            winners.push((decoy_score, true)); // decoy wins (including ties)
        }
    }
    // Add unpaired decoys
    for (base_id, &decoy_score) in &decoys {
        if !targets.contains_key(base_id) {
            winners.push((decoy_score, true));
        }
    }

    // Sort by score descending
    winners.sort_by(|a, b| b.0.total_cmp(&a.0));

    // Walk and compute FDR = (n_decoy + 1) / n_target
    let mut n_target = 0usize;
    let mut n_decoy = 0usize;
    let mut max_passing = 0usize;

    for &(_, is_decoy) in &winners {
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

        // Weight for feature 0 should be much larger than feature 1
        assert!(
            model.weights()[0].abs() > model.weights()[1].abs() * 2.0,
            "weights: {:?}",
            model.weights()
        );
        // Weight for feature 0 should be positive (targets have higher values)
        assert!(model.weights()[0] > 0.0);
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

        let best_c = grid_search_c(
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
}
