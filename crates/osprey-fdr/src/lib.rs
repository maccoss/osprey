//! Osprey FDR - False Discovery Rate control
//!
//! This crate provides:
//! - Target-decoy competition
//! - Q-value calculation
//! - Run-level and experiment-level FDR
//! - Integration with mokapot for semi-supervised learning

pub mod mokapot;

pub use mokapot::{
    get_pin_feature_names, pin_feature_value, MokapotResult, MokapotRunner, PsmFeatures,
    NUM_PIN_FEATURES,
};

use osprey_core::Result;
use std::collections::HashMap;

/// FDR controller using target-decoy competition
#[derive(Debug, Default)]
pub struct FdrController {
    /// FDR threshold
    fdr_threshold: f64,
}

/// Result of target-decoy competition
#[derive(Debug, Clone)]
pub struct CompetitionResult<T> {
    /// Winning items that pass FDR threshold (targets only)
    pub passing_targets: Vec<T>,
    /// Number of target winners
    pub n_target_wins: usize,
    /// Number of decoy winners
    pub n_decoy_wins: usize,
    /// FDR at the threshold
    pub fdr_at_threshold: f64,
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

    /// Run target-decoy competition and return targets passing FDR threshold
    ///
    /// This is the proper competition approach (matching pyXcorrDIA):
    /// 1. Each target competes with its paired decoy - higher score wins
    /// 2. Winners are sorted by score descending
    /// 3. Walk down list computing FDR = decoy_wins / target_wins at each position
    /// 4. Find the MAXIMUM cumulative_targets at any position where FDR <= threshold
    /// 5. Return all targets up to that count
    ///
    /// # Arguments
    /// * `matches` - Iterator of (item, score, is_decoy, target_id) tuples
    ///   - target_id links targets to decoys (decoy has same target_id as its target)
    ///
    /// # Returns
    /// CompetitionResult with passing targets and statistics
    pub fn compete_and_filter<T, I>(&self, matches: I) -> CompetitionResult<T>
    where
        T: Clone,
        I: IntoIterator<Item = (T, f64, bool, u32)>,
    {
        // Group by target_id (mask off high bit to get base ID)
        let mut target_scores: HashMap<u32, (T, f64)> = HashMap::new();
        let mut decoy_scores: HashMap<u32, f64> = HashMap::new();

        for (item, score, is_decoy, entry_id) in matches {
            // Get base target ID (clear the decoy flag bit)
            let target_id = entry_id & 0x7FFFFFFF;

            if is_decoy {
                decoy_scores
                    .entry(target_id)
                    .and_modify(|s| {
                        if score > *s {
                            *s = score;
                        }
                    })
                    .or_insert(score);
            } else {
                target_scores
                    .entry(target_id)
                    .and_modify(|(existing_item, s)| {
                        if score > *s {
                            *existing_item = item.clone();
                            *s = score;
                        }
                    })
                    .or_insert((item, score));
            }
        }

        log::debug!(
            "FDR input: {} unique targets, {} unique decoys (by base ID)",
            target_scores.len(),
            decoy_scores.len()
        );

        // Competition: for each target, compare with its decoy
        // Winner advances with (score, is_target_winner, Option<item>)
        // Note: Ties go to decoy (conservative for FDR estimation)
        let mut winners: Vec<(f64, bool, Option<T>)> = Vec::with_capacity(target_scores.len());

        // Debug counters
        let mut target_wins_count = 0usize;
        let mut decoy_wins_count = 0usize;
        let mut missing_decoy_count = 0usize;

        for (target_id, (item, target_score)) in target_scores {
            let decoy_score = decoy_scores
                .get(&target_id)
                .copied()
                .unwrap_or(f64::NEG_INFINITY);

            if !decoy_scores.contains_key(&target_id) {
                missing_decoy_count += 1;
            }

            if target_score > decoy_score {
                // Target wins (strict greater than)
                winners.push((target_score, true, Some(item)));
                target_wins_count += 1;
            } else {
                // Decoy wins (including ties - conservative)
                winners.push((decoy_score, false, None));
                decoy_wins_count += 1;
            }
        }

        log::debug!(
            "Competition results: {} targets won, {} decoys won, {} missing decoys",
            target_wins_count,
            decoy_wins_count,
            missing_decoy_count
        );

        // Sort winners by score descending (highest scores first)
        winners.sort_by(|a, b| b.0.total_cmp(&a.0));

        // Log score ranges and first decoy position for diagnostics
        if !winners.is_empty() {
            let first_decoy_rank = winners.iter().position(|(_, is_target, _)| !*is_target);
            let target_scores: Vec<f64> = winners
                .iter()
                .filter(|(_, is_t, _)| *is_t)
                .map(|(s, _, _)| *s)
                .collect();
            let decoy_scores: Vec<f64> = winners
                .iter()
                .filter(|(_, is_t, _)| !*is_t)
                .map(|(s, _, _)| *s)
                .collect();

            log::debug!(
                "FDR walk: {} winners, first decoy at rank {} | target scores: [{:.2}, {:.2}] | decoy scores: [{:.2}, {:.2}]",
                winners.len(),
                first_decoy_rank.map_or("none".to_string(), |r| (r + 1).to_string()),
                target_scores.first().copied().unwrap_or(0.0),
                target_scores.last().copied().unwrap_or(0.0),
                decoy_scores.first().copied().unwrap_or(0.0),
                decoy_scores.last().copied().unwrap_or(0.0),
            );
        }

        // First pass: walk down and find MAX cumulative_targets at any position where FDR <= threshold
        // This matches pyXcorrDIA's approach: valid['cumulative_targets'].max()
        let mut n_target_wins = 0usize;
        let mut n_decoy_wins = 0usize;
        let mut max_targets_at_valid_fdr = 0usize;
        let mut fdr_at_threshold = 0.0;

        for (_, is_target_winner, _) in winners.iter() {
            if *is_target_winner {
                n_target_wins += 1;
            } else {
                n_decoy_wins += 1;
            }

            // FDR at this position
            let fdr = if n_target_wins > 0 {
                n_decoy_wins as f64 / n_target_wins as f64
            } else {
                1.0
            };

            // Track the maximum targets at any position where FDR is valid
            if fdr <= self.fdr_threshold {
                max_targets_at_valid_fdr = n_target_wins;
                fdr_at_threshold = fdr;
            }
        }

        // Second pass: collect the first max_targets_at_valid_fdr target winners
        let mut passing_targets: Vec<T> = Vec::with_capacity(max_targets_at_valid_fdr);
        let mut targets_collected = 0usize;

        for (_, is_target_winner, item_opt) in winners {
            if is_target_winner {
                targets_collected += 1;
                if targets_collected <= max_targets_at_valid_fdr {
                    if let Some(item) = item_opt {
                        passing_targets.push(item);
                    }
                } else {
                    break; // We've collected enough targets
                }
            }
        }

        CompetitionResult {
            passing_targets,
            n_target_wins,
            n_decoy_wins,
            fdr_at_threshold,
        }
    }

    /// Compute q-values from target and decoy scores (legacy method)
    ///
    /// Note: For proper target-decoy competition where targets and decoys are paired,
    /// use `compete_and_filter` instead. This method treats targets and decoys as
    /// independent pools.
    ///
    /// Returns q-values for each target (in same order as input)
    pub fn compute_qvalues(&self, target_scores: &[f64], decoy_scores: &[f64]) -> Result<Vec<f64>> {
        if target_scores.is_empty() {
            return Ok(Vec::new());
        }

        // Track original indices for targets so we can write q-values directly
        // Format: (score, is_target, original_index_if_target)
        let mut all_scores: Vec<(f64, bool, Option<usize>)> = target_scores
            .iter()
            .enumerate()
            .map(|(idx, &s)| (s, true, Some(idx)))
            .chain(decoy_scores.iter().map(|&s| (s, false, None)))
            .collect();

        // Sort descending by score (higher is better)
        all_scores.sort_by(|a, b| b.0.total_cmp(&a.0));

        // Compute FDR at each position and track target indices
        let mut target_count = 0usize;
        let mut decoy_count = 0usize;
        let mut fdr_and_indices: Vec<(f64, Option<usize>)> = Vec::with_capacity(all_scores.len());

        for (_, is_target, orig_idx) in &all_scores {
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

            fdr_and_indices.push((fdr, *orig_idx));
        }

        // Convert FDR to q-value (minimum FDR at this score or lower rank)
        // Work backwards to ensure monotonicity, and directly assign to output
        let mut qvalues = vec![1.0f64; target_scores.len()];
        let mut min_fdr = 1.0f64;

        for (fdr, orig_idx) in fdr_and_indices.iter().rev() {
            min_fdr = min_fdr.min(*fdr);
            if let Some(idx) = orig_idx {
                qvalues[*idx] = min_fdr;
            }
        }

        Ok(qvalues)
    }

    /// Filter targets by q-value threshold
    pub fn filter_by_qvalue<T: Clone>(&self, items: &[T], qvalues: &[f64]) -> Vec<T> {
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

    /// Verifies that q-values are computed correctly from target and decoy score distributions.
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

    /// Verifies that items are correctly filtered to only those at or below the q-value threshold.
    #[test]
    fn test_filter_by_qvalue() {
        let controller = FdrController::new(0.05);
        let items = vec!["a", "b", "c", "d"];
        let qvalues = vec![0.01, 0.03, 0.10, 0.20];

        let filtered = controller.filter_by_qvalue(&items, &qvalues);
        assert_eq!(filtered, vec!["a", "b"]);
    }

    /// Verifies that detection counts at standard FDR thresholds (0.1%, 1%, 5%, 10%) are tallied correctly.
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

    // ============================================
    // Tests for compete_and_filter (pyXcorrDIA style)
    // ============================================

    /// Verifies that target-decoy competition selects the target when its score exceeds the decoy.
    #[test]
    fn test_competition_target_wins() {
        // Target has higher score than its decoy -> target wins
        let controller = FdrController::new(0.10); // 10% FDR

        // Format: (item, score, is_decoy, entry_id)
        // Target ID = 1, Decoy ID = 1 | 0x80000000
        let matches = vec![
            ("target_1", 0.9, false, 1),            // Target, score 0.9
            ("decoy_1", 0.7, true, 1 | 0x80000000), // Decoy, score 0.7
        ];

        let result = controller.compete_and_filter(matches);

        // Target wins (0.9 > 0.7)
        assert_eq!(result.n_target_wins, 1);
        assert_eq!(result.n_decoy_wins, 0);
        assert_eq!(result.passing_targets.len(), 1);
        assert_eq!(result.passing_targets[0], "target_1");
    }

    /// Verifies that the decoy wins competition when its score is higher than the target.
    #[test]
    fn test_competition_decoy_wins() {
        // Decoy has higher score than target -> decoy wins
        let controller = FdrController::new(0.10);

        let matches = vec![
            ("target_1", 0.5, false, 1),
            ("decoy_1", 0.8, true, 1 | 0x80000000),
        ];

        let result = controller.compete_and_filter(matches);

        // Decoy wins (0.8 > 0.5)
        assert_eq!(result.n_target_wins, 0);
        assert_eq!(result.n_decoy_wins, 1);
        assert_eq!(result.passing_targets.len(), 0); // No targets pass
    }

    /// Verifies that tied scores are conservatively awarded to the decoy for FDR estimation.
    #[test]
    fn test_competition_tie_goes_to_decoy() {
        // Equal scores -> decoy wins (conservative for FDR)
        let controller = FdrController::new(0.10);

        let matches = vec![
            ("target_1", 0.75, false, 1),
            ("decoy_1", 0.75, true, 1 | 0x80000000), // Same score
        ];

        let result = controller.compete_and_filter(matches);

        // Tie goes to decoy
        assert_eq!(result.n_target_wins, 0);
        assert_eq!(result.n_decoy_wins, 1);
        assert_eq!(result.passing_targets.len(), 0);
    }

    /// Verifies that cumulative FDR is computed correctly as decoy_wins/target_wins and only targets below the threshold pass.
    #[test]
    fn test_competition_fdr_calculation() {
        // Test FDR = decoy_wins / target_wins at each threshold
        let controller = FdrController::new(0.10); // 10% FDR

        // 5 pairs: 4 targets win, 1 decoy wins
        // Scores designed so FDR stays below 10% for top targets
        let matches = vec![
            // Pair 1: target wins with score 0.95
            ("target_1", 0.95, false, 1),
            ("decoy_1", 0.50, true, 1 | 0x80000000),
            // Pair 2: target wins with score 0.90
            ("target_2", 0.90, false, 2),
            ("decoy_2", 0.40, true, 2 | 0x80000000),
            // Pair 3: target wins with score 0.85
            ("target_3", 0.85, false, 3),
            ("decoy_3", 0.30, true, 3 | 0x80000000),
            // Pair 4: target wins with score 0.80
            ("target_4", 0.80, false, 4),
            ("decoy_4", 0.20, true, 4 | 0x80000000),
            // Pair 5: DECOY wins with score 0.82 (higher than target's 0.75)
            ("target_5", 0.75, false, 5),
            ("decoy_5", 0.82, true, 5 | 0x80000000),
        ];

        let result = controller.compete_and_filter(matches);

        // Winners sorted by score: 0.95(T), 0.90(T), 0.85(T), 0.82(D), 0.80(T)
        // Walk down:
        //   0.95(T): 1T, 0D -> FDR = 0/1 = 0%
        //   0.90(T): 2T, 0D -> FDR = 0/2 = 0%
        //   0.85(T): 3T, 0D -> FDR = 0/3 = 0%
        //   0.82(D): 3T, 1D -> FDR = 1/3 = 33% > 10%, stop adding targets
        //   0.80(T): 4T, 1D -> FDR = 1/4 = 25% > 10%

        assert_eq!(result.n_target_wins, 4);
        assert_eq!(result.n_decoy_wins, 1);
        // Only 3 targets pass at 10% FDR (before decoy win pushed FDR > 10%)
        assert_eq!(result.passing_targets.len(), 3);
    }

    /// Verifies that all targets pass when every target beats its paired decoy, yielding 0% FDR.
    #[test]
    fn test_competition_multiple_pairs_all_targets_win() {
        // All targets beat their decoys
        let controller = FdrController::new(0.01); // 1% FDR

        let matches = vec![
            ("t1", 0.9, false, 1),
            ("d1", 0.1, true, 1 | 0x80000000),
            ("t2", 0.8, false, 2),
            ("d2", 0.2, true, 2 | 0x80000000),
            ("t3", 0.7, false, 3),
            ("d3", 0.3, true, 3 | 0x80000000),
        ];

        let result = controller.compete_and_filter(matches);

        // All targets win, 0 decoys win -> FDR = 0/3 = 0%
        assert_eq!(result.n_target_wins, 3);
        assert_eq!(result.n_decoy_wins, 0);
        assert_eq!(result.passing_targets.len(), 3);
        assert!(result.fdr_at_threshold < 0.001); // FDR is 0
    }

    /// Verifies that no targets pass when all decoys outscore their paired targets.
    #[test]
    fn test_competition_all_decoys_win() {
        // All decoys beat their targets
        let controller = FdrController::new(0.01);

        let matches = vec![
            ("t1", 0.1, false, 1),
            ("d1", 0.9, true, 1 | 0x80000000),
            ("t2", 0.2, false, 2),
            ("d2", 0.8, true, 2 | 0x80000000),
            ("t3", 0.3, false, 3),
            ("d3", 0.7, true, 3 | 0x80000000),
        ];

        let result = controller.compete_and_filter(matches);

        // All decoys win -> no targets pass
        assert_eq!(result.n_target_wins, 0);
        assert_eq!(result.n_decoy_wins, 3);
        assert_eq!(result.passing_targets.len(), 0);
    }

    /// Verifies that when a peptide is scored multiple times, only the best score per target ID is retained for competition.
    #[test]
    fn test_competition_keeps_best_score_per_peptide() {
        // Same peptide scored multiple times, keep best
        let controller = FdrController::new(0.10);

        let matches = vec![
            // Target 1 scored twice - keep best (0.9)
            ("t1_v1", 0.7, false, 1),
            ("t1_v2", 0.9, false, 1), // Better score
            // Decoy 1 scored twice - keep best (0.5)
            ("d1_v1", 0.3, true, 1 | 0x80000000),
            ("d1_v2", 0.5, true, 1 | 0x80000000), // Better score
        ];

        let result = controller.compete_and_filter(matches);

        // Target's best (0.9) vs Decoy's best (0.5) -> Target wins
        assert_eq!(result.n_target_wins, 1);
        assert_eq!(result.n_decoy_wins, 0);
        assert_eq!(result.passing_targets.len(), 1);
        assert_eq!(result.passing_targets[0], "t1_v2"); // The better-scoring version
    }

    /// Verifies that competition with no input produces zero wins and an empty result set.
    #[test]
    fn test_competition_empty_input() {
        let controller = FdrController::new(0.01);
        let matches: Vec<(&str, f64, bool, u32)> = vec![];

        let result = controller.compete_and_filter(matches);

        assert_eq!(result.n_target_wins, 0);
        assert_eq!(result.n_decoy_wins, 0);
        assert_eq!(result.passing_targets.len(), 0);
    }

    /// Verifies that a target with no corresponding decoy wins by default against negative infinity.
    #[test]
    fn test_competition_target_without_decoy() {
        // Target has no corresponding decoy (decoy score defaults to 0)
        let controller = FdrController::new(0.10);

        let matches = vec![
            ("t1", 0.5, false, 1), // Target only, no decoy for ID 1
        ];

        let result = controller.compete_and_filter(matches);

        // Target (0.5) vs missing decoy (0.0) -> Target wins
        assert_eq!(result.n_target_wins, 1);
        assert_eq!(result.n_decoy_wins, 0);
        assert_eq!(result.passing_targets.len(), 1);
    }

    /// Verifies that FDR recovery after a decoy spike allows later targets to pass, matching pyXcorrDIA's max cumulative targets behavior.
    #[test]
    fn test_competition_fdr_recovers_after_spike() {
        // Key test: FDR spikes above threshold early, but recovers as more targets accumulate
        // This tests the pyXcorrDIA behavior: max(cumulative_targets) where FDR <= threshold
        let controller = FdrController::new(0.10); // 10% FDR

        // Create a scenario where:
        // - 2 targets win at top scores
        // - 1 decoy wins (pushing FDR to 50%)
        // - Then 10 more targets win (FDR drops to 1/12 = 8.3%)
        let mut matches = vec![
            // Pair 1: target wins with score 0.99
            ("target_1", 0.99, false, 1),
            ("decoy_1", 0.10, true, 1 | 0x80000000),
            // Pair 2: target wins with score 0.98
            ("target_2", 0.98, false, 2),
            ("decoy_2", 0.10, true, 2 | 0x80000000),
            // Pair 3: DECOY wins with score 0.97 (higher than target's 0.50)
            ("target_3", 0.50, false, 3),
            ("decoy_3", 0.97, true, 3 | 0x80000000),
        ];

        // Add 10 more target winners with scores 0.90 down to 0.81
        for i in 4..=13 {
            let score = 0.90 - (i - 4) as f64 * 0.01;
            matches.push((
                Box::leak(format!("target_{}", i).into_boxed_str()),
                score,
                false,
                i as u32,
            ));
            matches.push((
                Box::leak(format!("decoy_{}", i).into_boxed_str()),
                0.05,
                true,
                (i as u32) | 0x80000000,
            ));
        }

        let result = controller.compete_and_filter(matches);

        // Walk through:
        // - 0.99(T): 1T, 0D -> FDR = 0% ✓
        // - 0.98(T): 2T, 0D -> FDR = 0% ✓
        // - 0.97(D): 2T, 1D -> FDR = 50% ✗ (spike!)
        // - 0.90(T): 3T, 1D -> FDR = 33% ✗
        // - 0.89(T): 4T, 1D -> FDR = 25% ✗
        // - 0.88(T): 5T, 1D -> FDR = 20% ✗
        // - 0.87(T): 6T, 1D -> FDR = 16.7% ✗
        // - 0.86(T): 7T, 1D -> FDR = 14.3% ✗
        // - 0.85(T): 8T, 1D -> FDR = 12.5% ✗
        // - 0.84(T): 9T, 1D -> FDR = 11.1% ✗
        // - 0.83(T): 10T, 1D -> FDR = 10% ✓ (recovered!)
        // - 0.82(T): 11T, 1D -> FDR = 9.1% ✓
        // - 0.81(T): 12T, 1D -> FDR = 8.3% ✓

        assert_eq!(result.n_target_wins, 12); // All 12 targets won
        assert_eq!(result.n_decoy_wins, 1); // 1 decoy won

        // OLD BUGGY BEHAVIOR: Only 2 targets would pass (before decoy spike)
        // NEW CORRECT BEHAVIOR: 12 targets pass (max cumulative at FDR <= 10%)
        assert_eq!(result.passing_targets.len(), 12);
        assert!(result.fdr_at_threshold <= 0.10);
    }
}
