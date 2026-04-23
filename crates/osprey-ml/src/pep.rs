//! Posterior Error Probability estimation using KDE + isotonic regression
//!
//! Implements the Percolator/qvality-style PEP calculation:
//! 1. Target-decoy competition (done externally, before calling this module)
//! 2. KDE density estimation for target and decoy score distributions
//! 3. Bayes' rule for posterior probability of being incorrect
//! 4. Isotonic regression (PAVA) for monotonicity enforcement
//!
//! References:
//! - Käll et al. (2008) "Posterior error probabilities and false discovery rates"
//! - The qvality algorithm as implemented in triqler/mokapot

use super::*;

/// Posterior error probability estimator using KDE + isotonic regression
///
/// Computes PEP for competition winners (not all entries).
/// PEP(s) = π₀ · f_decoy(s) / [π₀ · f_decoy(s) + (1 - π₀) · f_target(s)]
pub struct PepEstimator {
    /// PEP values at evenly-spaced score bins
    bins: Vec<f64>,
    /// Minimum score in the binning range
    min_score: f64,
    /// Score step between bins
    score_step: f64,
}

impl PepEstimator {
    /// Fit PEP model on competition winners
    ///
    /// # Arguments
    /// * `scores` - Scores of competition winners (not required to be sorted)
    /// * `is_decoy` - Whether each winner is a decoy
    /// * `n_bins` - Number of bins for score discretization (default: 1000)
    pub fn fit(scores: &[f64], is_decoy: &[bool], n_bins: usize) -> Self {
        assert_eq!(scores.len(), is_decoy.len());

        if scores.is_empty() {
            return PepEstimator {
                bins: vec![1.0],
                min_score: 0.0,
                score_step: 1.0,
            };
        }

        // Separate target and decoy scores
        let decoy_scores: Vec<f64> = scores
            .iter()
            .zip(is_decoy)
            .filter(|&(_, &d)| d)
            .map(|(&s, _)| s)
            .collect();

        let target_scores: Vec<f64> = scores
            .iter()
            .zip(is_decoy)
            .filter(|&(_, &d)| !d)
            .map(|(&s, _)| s)
            .collect();

        if target_scores.is_empty() || decoy_scores.is_empty() {
            // Can't estimate PEP without both classes
            return PepEstimator {
                bins: vec![1.0; n_bins],
                min_score: scores.iter().cloned().fold(f64::INFINITY, f64::min),
                score_step: 1.0,
            };
        }

        // π₀: prior probability that a target is incorrect
        // Estimated as n_decoy / n_target, clamped for stability
        let pi0 = (decoy_scores.len() as f64 / target_scores.len() as f64).clamp(0.01, 0.99);

        // KDE for each distribution
        let decoy_kde = Kde::new(&decoy_scores);
        let target_kde = Kde::new(&target_scores);

        // Score range for binning
        let min_score = scores.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_score = scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let score_step = if n_bins > 1 {
            (max_score - min_score) / (n_bins - 1) as f64
        } else {
            1.0
        };

        // Compute raw PEP at each bin
        let mut bins: Vec<f64> = (0..n_bins)
            .map(|i| {
                let score = min_score + i as f64 * score_step;
                let f_decoy = decoy_kde.pdf(score) * pi0;
                let f_target = target_kde.pdf(score) * (1.0 - pi0);
                let denom = f_decoy + f_target;
                if denom > 0.0 {
                    (f_decoy / denom).clamp(0.0, 1.0)
                } else {
                    1.0
                }
            })
            .collect();

        // Enforce monotonicity: PEP must be non-increasing as score increases
        // (higher scores = more confident = lower PEP)
        isotonic_regression_decreasing(&mut bins);

        PepEstimator {
            bins,
            min_score,
            score_step,
        }
    }

    /// Fit with default 1000 bins
    pub fn fit_default(scores: &[f64], is_decoy: &[bool]) -> Self {
        Self::fit(scores, is_decoy, 1000)
    }

    /// Look up PEP for a given score using linear interpolation
    pub fn posterior_error(&self, score: f64) -> f64 {
        if self.bins.is_empty() {
            return 1.0;
        }

        let bin_lo = self.bins.len().saturating_sub(1).min(
            ((score - self.min_score) / self.score_step)
                .max(0.0)
                .floor() as usize,
        );
        let bin_hi = self.bins.len().saturating_sub(1).min(bin_lo + 1);

        let lower = self.bins[bin_lo];
        let upper = self.bins[bin_hi];

        // Linear interpolation
        let bin_lo_score = bin_lo as f64 * self.score_step + self.min_score;
        let frac = ((score - bin_lo_score) / self.score_step).clamp(0.0, 1.0);

        let pep = lower + (upper - lower) * frac;
        pep.clamp(0.0, 1.0)
    }
}

/// Gaussian KDE for density estimation
struct Kde {
    sample: Vec<f64>,
    bandwidth: f64,
    constant: f64,
}

impl Kde {
    /// Create KDE with Silverman's rule bandwidth
    fn new(sample: &[f64]) -> Self {
        let sigma = std(sample);
        let n = sample.len() as f64;
        let factor = 4.0 / 3.0;
        let exponent = 1.0 / 5.0;
        let bandwidth = sigma * (factor / n).powf(exponent);
        let constant = (2.0 * std::f64::consts::PI).sqrt() * bandwidth * n;

        Kde {
            sample: sample.to_vec(),
            bandwidth: if bandwidth > 0.0 { bandwidth } else { 1.0 },
            constant: if constant > 0.0 { constant } else { 1.0 },
        }
    }

    /// Evaluate PDF at a point.
    ///
    /// Summation is serial (not `par_iter`) because Rayon's work-stealing
    /// reduction tree is non-deterministic: scheduling differences lead to
    /// different tree shapes, and floating-point `+` is non-associative,
    /// so parallel reduction produces pep values that differ by ~1 ULP
    /// across runs. Stage 5 calls this once for each of ~1000 bins over ~241k
    /// samples (~sub-second total), so the perf cost of serial reduction
    /// is negligible compared to the Percolator SVM pass that precedes it.
    fn pdf(&self, x: f64) -> f64 {
        let h = self.bandwidth;
        let sum: f64 = self
            .sample
            .iter()
            .fold(0.0, |acc, &xi| acc + (-0.5 * ((x - xi) / h).powi(2)).exp());
        sum / self.constant
    }
}

/// Pool Adjacent Violators Algorithm (PAVA) for isotonic regression
///
/// Enforces that the output is monotonically non-increasing
/// (higher index = higher score = lower PEP).
///
/// This is more principled than a simple cumulative-max sweep:
/// when violations are found, PAVA replaces the violating block
/// with the average of its values, producing smoother estimates.
pub fn isotonic_regression_decreasing(values: &mut [f64]) {
    if values.len() <= 1 {
        return;
    }

    let n = values.len();

    // We want non-increasing: values[0] >= values[1] >= ... >= values[n-1]
    // PAVA: process left to right, merge blocks when violation found
    // A "block" is a contiguous range with a single averaged value.
    //
    // Store blocks as (start_index, end_index_exclusive, value)
    let mut blocks: Vec<(usize, usize, f64)> = Vec::with_capacity(n);

    for (i, &val) in values.iter().enumerate() {
        // Start a new block with this single element
        blocks.push((i, i + 1, val));

        // Merge with previous block while we have a violation
        // (current block value > previous block value means non-increasing is violated)
        while blocks.len() >= 2 {
            let len = blocks.len();
            let prev = &blocks[len - 2];
            let curr = &blocks[len - 1];

            if curr.2 > prev.2 {
                // Violation: current block has higher value than previous
                // Merge: compute weighted average
                let prev_count = (prev.1 - prev.0) as f64;
                let curr_count = (curr.1 - curr.0) as f64;
                let avg = (prev.2 * prev_count + curr.2 * curr_count) / (prev_count + curr_count);

                let new_start = prev.0;
                let new_end = curr.1;

                blocks.pop();
                blocks.pop();
                blocks.push((new_start, new_end, avg));
            } else {
                break;
            }
        }
    }

    // Write merged block values back
    for (start, end, value) in blocks {
        for v in &mut values[start..end] {
            *v = value;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isotonic_regression_already_decreasing() {
        let mut values = vec![1.0, 0.8, 0.6, 0.4, 0.2];
        isotonic_regression_decreasing(&mut values);
        assert_eq!(values, vec![1.0, 0.8, 0.6, 0.4, 0.2]);
    }

    #[test]
    fn test_isotonic_regression_single_violation() {
        // values[2] > values[1] is a violation
        let mut values = vec![1.0, 0.3, 0.5, 0.2, 0.1];
        isotonic_regression_decreasing(&mut values);

        // After PAVA: indices 1,2 should be averaged → (0.3+0.5)/2 = 0.4
        assert!((values[1] - 0.4).abs() < 1e-10);
        assert!((values[2] - 0.4).abs() < 1e-10);

        // Check monotonicity
        for i in 1..values.len() {
            assert!(
                values[i] <= values[i - 1] + 1e-10,
                "values[{}]={} > values[{}]={}",
                i,
                values[i],
                i - 1,
                values[i - 1]
            );
        }
    }

    #[test]
    fn test_isotonic_regression_all_increasing() {
        let mut values = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        isotonic_regression_decreasing(&mut values);

        // All should become the average: (0.1+0.2+0.3+0.4+0.5)/5 = 0.3
        for &v in &values {
            assert!((v - 0.3).abs() < 1e-10, "v = {}", v);
        }
    }

    #[test]
    fn test_isotonic_regression_empty_and_single() {
        let mut empty: Vec<f64> = vec![];
        isotonic_regression_decreasing(&mut empty);

        let mut single = vec![0.5];
        isotonic_regression_decreasing(&mut single);
        assert_eq!(single, vec![0.5]);
    }

    #[test]
    fn test_pep_well_separated() {
        // Targets have high scores, decoys have low scores
        let mut scores = Vec::new();
        let mut is_decoy = Vec::new();

        // 100 targets with scores around 5.0
        for i in 0..100 {
            scores.push(4.0 + (i as f64) * 0.02);
            is_decoy.push(false);
        }
        // 100 decoys with scores around 1.0
        for i in 0..100 {
            scores.push(0.0 + (i as f64) * 0.02);
            is_decoy.push(true);
        }

        let estimator = PepEstimator::fit_default(&scores, &is_decoy);

        // High-scoring target should have low PEP
        let pep_high = estimator.posterior_error(5.0);
        assert!(pep_high < 0.1, "PEP at high score: {}", pep_high);

        // Low-scoring decoy should have high PEP
        let pep_low = estimator.posterior_error(0.5);
        assert!(pep_low > 0.5, "PEP at low score: {}", pep_low);
    }

    #[test]
    fn test_pep_monotonicity() {
        // Create scores with some overlap
        let mut scores = Vec::new();
        let mut is_decoy = Vec::new();

        for i in 0..50 {
            scores.push(3.0 + (i as f64) * 0.1);
            is_decoy.push(false);
        }
        for i in 0..50 {
            scores.push(1.0 + (i as f64) * 0.1);
            is_decoy.push(true);
        }

        let estimator = PepEstimator::fit_default(&scores, &is_decoy);

        // PEP should be monotonically non-increasing with score
        let test_scores: Vec<f64> = (0..20).map(|i| 0.5 + i as f64 * 0.5).collect();
        let peps: Vec<f64> = test_scores
            .iter()
            .map(|&s| estimator.posterior_error(s))
            .collect();

        for i in 1..peps.len() {
            assert!(
                peps[i] <= peps[i - 1] + 1e-6,
                "PEP not monotonic: score={}, PEP={} > prev PEP={}",
                test_scores[i],
                peps[i],
                peps[i - 1]
            );
        }
    }

    #[test]
    fn test_pep_empty_input() {
        let estimator = PepEstimator::fit_default(&[], &[]);
        assert!((estimator.posterior_error(1.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pep_all_targets() {
        let scores = vec![5.0, 4.0, 3.0];
        let is_decoy = vec![false, false, false];
        let estimator = PepEstimator::fit_default(&scores, &is_decoy);
        // With no decoys, PEP should be 1.0 (can't estimate)
        assert!((estimator.posterior_error(5.0) - 1.0).abs() < 1e-10);
    }

    /// Regression for Copilot PR #16 review: `PepEstimator::fit_default`
    /// must produce bit-identical bins across consecutive calls with
    /// the same input. This is the exact invariant the serial-KDE-sum
    /// fix provides: `par_iter().sum()` would yield different
    /// reduction-tree shapes across calls depending on Rayon's
    /// work-stealing state, drifting bins at ~1 ULP per accumulated
    /// term. Serial `iter().fold()` is fully deterministic on a fixed
    /// input, so two calls must give identical output.
    ///
    /// This test intentionally does NOT permute the input: serial
    /// summation is still order-dependent at the last ULP for floating
    /// point, so a permutation-invariance test would be too strict. The
    /// fix guarantees same-input reproducibility, not multiset
    /// invariance — same-input reproducibility is what the downstream
    /// pipeline relies on for cross-run PEP stability.
    #[test]
    fn fit_default_is_deterministic() {
        let n = 200;
        let mut scores: Vec<f64> = Vec::with_capacity(n);
        let mut is_decoy: Vec<bool> = Vec::with_capacity(n);
        for i in 0..n {
            let t = (i as f64) * 0.037 + 1.5;
            let d = (i as f64) * 0.037 + 0.5;
            scores.push(t);
            is_decoy.push(false);
            scores.push(d);
            is_decoy.push(true);
        }

        let a = PepEstimator::fit_default(&scores, &is_decoy);
        let b = PepEstimator::fit_default(&scores, &is_decoy);

        assert_eq!(
            a.bins, b.bins,
            "PepEstimator::fit_default is non-deterministic -- \
             KDE summation likely has a parallel reduction again"
        );
        assert_eq!(a.min_score, b.min_score);
        assert_eq!(a.score_step, b.score_step);
    }
}
