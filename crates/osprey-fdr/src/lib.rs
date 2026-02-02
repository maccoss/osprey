//! Osprey FDR - False Discovery Rate control
//!
//! This crate provides:
//! - Target-decoy competition
//! - Q-value calculation
//! - Run-level and experiment-level FDR
//! - Integration with mokapot for semi-supervised learning

pub mod mokapot;

pub use mokapot::{MokapotResult, MokapotRunner, PsmFeatures};

use osprey_core::Result;

/// FDR controller using target-decoy competition
#[derive(Debug, Default)]
pub struct FdrController {
    /// FDR threshold
    fdr_threshold: f64,
}

impl FdrController {
    /// Create a new FDR controller
    pub fn new(fdr_threshold: f64) -> Self {
        Self { fdr_threshold }
    }

    /// Get the FDR threshold
    pub fn threshold(&self) -> f64 {
        self.fdr_threshold
    }

    /// Compute q-values from target and decoy scores
    ///
    /// Uses the simple formula: q = (# decoys with score >= s) / (# targets with score >= s)
    ///
    /// Returns q-values for each target (in same order as input)
    pub fn compute_qvalues(
        &self,
        target_scores: &[f64],
        decoy_scores: &[f64],
    ) -> Result<Vec<f64>> {
        if target_scores.is_empty() {
            return Ok(Vec::new());
        }

        // Combine and sort all scores
        let mut all_scores: Vec<(f64, bool)> = target_scores
            .iter()
            .map(|&s| (s, true))
            .chain(decoy_scores.iter().map(|&s| (s, false)))
            .collect();

        // Sort descending by score (higher is better)
        all_scores.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));

        // Compute FDR at each score threshold
        let mut target_count = 0usize;
        let mut decoy_count = 0usize;
        let mut fdr_at_score: Vec<(f64, f64)> = Vec::with_capacity(all_scores.len());

        for (score, is_target) in &all_scores {
            if *is_target {
                target_count += 1;
            } else {
                decoy_count += 1;
            }

            // FDR = decoys / targets (or 1.0 if no targets)
            let fdr = if target_count > 0 {
                decoy_count as f64 / target_count as f64
            } else {
                1.0
            };

            fdr_at_score.push((*score, fdr));
        }

        // Convert FDR to q-value (minimum FDR at this score or higher)
        // Work backwards to ensure monotonicity
        let mut min_fdr = 1.0f64;
        let mut qvalue_at_score: Vec<(f64, f64)> = Vec::with_capacity(fdr_at_score.len());

        for (score, fdr) in fdr_at_score.iter().rev() {
            min_fdr = min_fdr.min(*fdr);
            qvalue_at_score.push((*score, min_fdr));
        }
        qvalue_at_score.reverse();

        // Map back to original target order
        let qvalues: Vec<f64> = target_scores
            .iter()
            .map(|&score| {
                // Find q-value for this score (binary search would be faster)
                qvalue_at_score
                    .iter()
                    .find(|(s, _)| (*s - score).abs() < 1e-10)
                    .map(|(_, q)| *q)
                    .unwrap_or(1.0)
            })
            .collect();

        Ok(qvalues)
    }

    /// Filter targets by q-value threshold
    pub fn filter_by_qvalue<T: Clone>(
        &self,
        items: &[T],
        qvalues: &[f64],
    ) -> Vec<T> {
        items
            .iter()
            .zip(qvalues.iter())
            .filter(|(_, &q)| q <= self.fdr_threshold)
            .map(|(item, _)| item.clone())
            .collect()
    }

    /// Count detections at various FDR thresholds
    pub fn count_at_thresholds(&self, qvalues: &[f64]) -> FdrCounts {
        FdrCounts {
            at_001: qvalues.iter().filter(|&&q| q <= 0.001).count(),
            at_01: qvalues.iter().filter(|&&q| q <= 0.01).count(),
            at_05: qvalues.iter().filter(|&&q| q <= 0.05).count(),
            at_10: qvalues.iter().filter(|&&q| q <= 0.10).count(),
            total: qvalues.len(),
        }
    }
}

/// Detection counts at various FDR thresholds
#[derive(Debug, Clone, Default)]
pub struct FdrCounts {
    /// Count at 0.1% FDR
    pub at_001: usize,
    /// Count at 1% FDR
    pub at_01: usize,
    /// Count at 5% FDR
    pub at_05: usize,
    /// Count at 10% FDR
    pub at_10: usize,
    /// Total count
    pub total: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qvalue_computation() {
        let controller = FdrController::new(0.01);

        // Target scores: high values
        let targets = vec![10.0, 8.0, 6.0, 4.0, 2.0];
        // Decoy scores: lower values (one high decoy)
        let decoys = vec![9.0, 3.0, 1.0];

        let qvalues = controller.compute_qvalues(&targets, &decoys).unwrap();

        // Highest target (10.0) should have very low q-value
        assert!(qvalues[0] < 0.5);

        // Q-values should be monotonically related to scores
        // (not strictly monotonic due to decoy competition)
    }

    #[test]
    fn test_filter_by_qvalue() {
        let controller = FdrController::new(0.05);
        let items = vec!["a", "b", "c", "d"];
        let qvalues = vec![0.01, 0.03, 0.10, 0.20];

        let filtered = controller.filter_by_qvalue(&items, &qvalues);
        assert_eq!(filtered, vec!["a", "b"]);
    }

    #[test]
    fn test_count_at_thresholds() {
        let controller = FdrController::new(0.01);
        let qvalues = vec![0.001, 0.005, 0.02, 0.03, 0.08, 0.15];

        let counts = controller.count_at_thresholds(&qvalues);
        assert_eq!(counts.at_001, 1);
        assert_eq!(counts.at_01, 2);
        assert_eq!(counts.at_05, 4);
        assert_eq!(counts.at_10, 5);
        assert_eq!(counts.total, 6);
    }
}
