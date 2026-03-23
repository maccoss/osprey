//! Inter-replicate peak reconciliation for consistent peak integration across replicates.
//!
//! After initial per-run FDR, this module:
//! 1. Computes consensus library RTs for peptides detected across multiple runs
//! 2. Refits per-run LOESS calibration using consensus peptides
//! 3. Reconciles peak boundaries across runs (CWT peak overlap or forced integration)
//! 4. Runs a second FDR pass on the reconciled consensus set

use osprey_chromatography::RTCalibration;
use osprey_core::{CwtCandidate, FdrEntry};
use std::collections::HashMap;

/// Consensus RT for a peptide across all runs.
///
/// Computed independently for targets and decoys to maintain fair
/// target-decoy competition in the second FDR pass.
#[derive(Debug, Clone)]
pub struct PeptideConsensusRT {
    /// Modified sequence (peptide-level grouping)
    pub modified_sequence: String,
    /// Whether this is a decoy consensus
    pub is_decoy: bool,
    /// Consensus library RT (weighted median from all runs)
    pub consensus_library_rt: f64,
    /// Median peak width across runs where detected (minutes, measured RT space)
    pub median_peak_width: f64,
    /// Number of runs where this peptide was detected
    pub n_runs_detected: usize,
}

/// Compute consensus library RTs for peptides detected across runs.
///
/// For each target peptide passing `consensus_fdr` in any run, and its
/// paired decoy, computes a consensus library RT using the weighted median
/// of per-run detections mapped back to library RT space.
///
/// Targets and decoys get independent consensus RTs to maintain fair
/// competition — both sides benefit equally from inter-replicate consensus.
///
/// # Arguments
/// * `per_file_entries` - Per-file scored entries (after first-pass FDR)
/// * `per_file_calibrations` - Per-file RT calibrations for inverse prediction
/// * `consensus_fdr` - FDR threshold for selecting consensus peptides (typically 0.01)
pub fn compute_consensus_rts(
    per_file_entries: &[(String, Vec<FdrEntry>)],
    per_file_calibrations: &HashMap<String, RTCalibration>,
    consensus_fdr: f64,
) -> Vec<PeptideConsensusRT> {
    // 1. Collect target peptides passing FDR in any run (run-level OR experiment-level)
    let mut target_peptides: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if !entry.is_decoy && entry.run_qvalue.min(entry.experiment_qvalue) <= consensus_fdr {
                target_peptides.insert(entry.modified_sequence.clone());
            }
        }
    }

    if target_peptides.is_empty() {
        return Vec::new();
    }

    log::info!(
        "Inter-replicate reconciliation: {} target peptides pass {:.1}% FDR in at least one run",
        target_peptides.len(),
        consensus_fdr * 100.0
    );

    // 2. Collect detections for targets and their paired decoys
    //    Group by (modified_sequence, is_decoy)
    //    Each detection: (file_name, apex_rt, coelution_sum, peak_width)
    type DetectionVec = Vec<(String, f64, f64, f64)>;
    let mut detections: HashMap<(String, bool), DetectionVec> = HashMap::new();

    // Identify decoy modified_sequences paired to target peptides
    // Decoys have DECOY_ prefix on modified_sequence
    let mut decoy_peptides: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if entry.is_decoy {
                // Check if the paired target is in our consensus set
                // Decoy modified_sequence has DECOY_ prefix
                let target_seq = entry
                    .modified_sequence
                    .strip_prefix("DECOY_")
                    .unwrap_or(&entry.modified_sequence);
                if target_peptides.contains(target_seq) {
                    decoy_peptides.insert(entry.modified_sequence.clone());
                }
            }
        }
    }

    log::info!(
        "Inter-replicate reconciliation: {} paired decoy peptides included",
        decoy_peptides.len()
    );

    // Collect all detections for consensus peptides (targets and decoys)
    for (file_name, entries) in per_file_entries {
        for entry in entries {
            let include = if entry.is_decoy {
                decoy_peptides.contains(&entry.modified_sequence)
            } else {
                target_peptides.contains(&entry.modified_sequence)
            };

            if include {
                let peak_width = entry.end_rt - entry.start_rt;
                detections
                    .entry((entry.modified_sequence.clone(), entry.is_decoy))
                    .or_default()
                    .push((
                        file_name.clone(),
                        entry.apex_rt,
                        entry.coelution_sum,
                        peak_width,
                    ));
            }
        }
    }

    // 3. Compute consensus RT for each peptide
    let mut consensus: Vec<PeptideConsensusRT> = Vec::new();

    for ((modified_sequence, is_decoy), dets) in &detections {
        if dets.is_empty() {
            continue;
        }

        // Map each detection's apex_rt to library RT space
        let mut library_rt_weights: Vec<(f64, f64)> = Vec::new();
        let mut peak_widths: Vec<f64> = Vec::new();

        for (file_name, apex_rt, coelution_sum, peak_width) in dets {
            if let Some(cal) = per_file_calibrations.get(file_name) {
                let library_rt = cal.inverse_predict(*apex_rt);
                if library_rt.is_finite() && *coelution_sum > 0.0 {
                    library_rt_weights.push((library_rt, *coelution_sum));
                    peak_widths.push(*peak_width);
                }
            }
        }

        if library_rt_weights.is_empty() {
            continue;
        }

        let consensus_library_rt = weighted_median(&library_rt_weights);
        let median_peak_width = simple_median(&peak_widths);
        let n_runs_detected = library_rt_weights.len();

        consensus.push(PeptideConsensusRT {
            modified_sequence: modified_sequence.clone(),
            is_decoy: *is_decoy,
            consensus_library_rt,
            median_peak_width,
            n_runs_detected,
        });
    }

    // Sort for deterministic output
    consensus.sort_by(|a, b| {
        a.is_decoy
            .cmp(&b.is_decoy)
            .then(a.modified_sequence.cmp(&b.modified_sequence))
    });

    let n_targets = consensus.iter().filter(|c| !c.is_decoy).count();
    let n_decoys = consensus.iter().filter(|c| c.is_decoy).count();
    let multi_run = consensus.iter().filter(|c| c.n_runs_detected > 1).count();
    log::info!(
        "Inter-replicate consensus: {} targets, {} decoys ({} detected in multiple runs)",
        n_targets,
        n_decoys,
        multi_run
    );

    consensus
}

/// Refit RT calibration using consensus peptides for a given run.
///
/// Uses (consensus_library_rt → run's detected apex_rt) pairs for target peptides
/// detected in this run. Produces a tighter calibration curve based on high-quality
/// inter-replicate consensus points.
///
/// Only uses target peptides (not decoys) for the calibration fit, since targets
/// have biologically meaningful RTs while decoy RTs are arbitrary.
///
/// # Returns
/// A refined RTCalibration, or None if too few consensus points for this run.
pub fn refit_calibration_with_consensus(
    consensus: &[PeptideConsensusRT],
    entries: &[FdrEntry],
    consensus_fdr: f64,
) -> Option<RTCalibration> {
    use osprey_chromatography::RTCalibrator;

    // Build (consensus_library_rt → measured_apex_rt) pairs for target peptides in this run
    let consensus_map: HashMap<&str, f64> = consensus
        .iter()
        .filter(|c| !c.is_decoy)
        .map(|c| (c.modified_sequence.as_str(), c.consensus_library_rt))
        .collect();

    let mut library_rts = Vec::new();
    let mut measured_rts = Vec::new();

    for entry in entries {
        if entry.is_decoy || entry.run_qvalue.min(entry.experiment_qvalue) > consensus_fdr {
            continue;
        }
        if let Some(&consensus_lib_rt) = consensus_map.get(entry.modified_sequence.as_str()) {
            library_rts.push(consensus_lib_rt);
            measured_rts.push(entry.apex_rt);
        }
    }

    if library_rts.len() < 20 {
        log::warn!(
            "Too few consensus points ({}) for second-pass calibration, skipping",
            library_rts.len()
        );
        return None;
    }

    // Disable outlier removal (retention=1.0) — these are FDR-controlled detections,
    // not noisy initial matches. Robustness iterations still downweight any bad points.
    let calibrator = RTCalibrator::new()
        .with_bandwidth(0.3)
        .with_outlier_retention(1.0);
    match calibrator.fit(&library_rts, &measured_rts) {
        Ok(cal) => {
            let stats = cal.stats();
            log::info!(
                "Second-pass RT calibration: {} points, R²={:.4}, residual_SD={:.3} min",
                stats.n_points,
                stats.r_squared,
                stats.residual_std
            );
            Some(cal)
        }
        Err(e) => {
            log::warn!("Second-pass RT calibration failed: {}", e);
            None
        }
    }
}

/// Result of inter-replicate peak reconciliation for a single entry.
#[derive(Debug, Clone)]
pub enum ReconcileAction {
    /// Keep existing peak (already matches consensus RT)
    Keep,
    /// Use a different stored CWT candidate that overlaps consensus RT
    UseCwtPeak {
        /// Index into the entry's cwt_candidates
        candidate_idx: usize,
        /// CWT peak start RT
        start_rt: f64,
        /// CWT peak apex RT
        apex_rt: f64,
        /// CWT peak end RT
        end_rt: f64,
    },
    /// No CWT peak overlaps — forced integration at consensus RT
    ForcedIntegration {
        /// Expected measured RT from refined calibration
        expected_rt: f64,
        /// Integration half-width (half of median_peak_width)
        half_width: f64,
    },
}

/// Determine what reconciliation action to take for an entry given its consensus RT.
///
/// Checks if the entry's current peak or any stored CWT candidate overlaps
/// the expected measured RT. If not, returns ForcedIntegration.
pub fn determine_reconcile_action(
    start_rt: f64,
    end_rt: f64,
    cwt_candidates: &[CwtCandidate],
    expected_measured_rt: f64,
    median_peak_width: f64,
) -> ReconcileAction {
    // Check if the current peak already contains the expected RT
    if start_rt <= expected_measured_rt && expected_measured_rt <= end_rt {
        return ReconcileAction::Keep;
    }

    // Check stored CWT candidates for overlap
    for (idx, cand) in cwt_candidates.iter().enumerate() {
        if cand.start_rt <= expected_measured_rt && expected_measured_rt <= cand.end_rt {
            return ReconcileAction::UseCwtPeak {
                candidate_idx: idx,
                start_rt: cand.start_rt,
                apex_rt: cand.apex_rt,
                end_rt: cand.end_rt,
            };
        }
    }

    // No overlap — forced integration
    ReconcileAction::ForcedIntegration {
        expected_rt: expected_measured_rt,
        half_width: median_peak_width / 2.0,
    }
}

/// Plan inter-replicate peak reconciliation for all entries across all runs.
///
/// Returns a map of (file_name, entry_index) → ReconcileAction for entries
/// that need modification. Entries not in the map should keep their current peaks.
///
/// This function only plans the reconciliation — actual re-scoring happens
/// in the pipeline where spectra are available.
///
/// `per_file_cwt_candidates` provides the CWT candidates per entry for each file,
/// loaded from parquet caches on demand. If a file has no CWT data, entries default
/// to forced integration when not matching the current peak.
pub fn plan_reconciliation(
    consensus: &[PeptideConsensusRT],
    per_file_entries: &[(String, Vec<FdrEntry>)],
    per_file_cwt_candidates: &HashMap<String, Vec<Vec<CwtCandidate>>>,
    per_file_refined_cal: &HashMap<String, RTCalibration>,
    per_file_original_cal: &HashMap<String, RTCalibration>,
) -> HashMap<(String, usize), ReconcileAction> {
    let mut actions: HashMap<(String, usize), ReconcileAction> = HashMap::new();

    // Build consensus lookup by (modified_sequence, is_decoy)
    let consensus_map: HashMap<(&str, bool), &PeptideConsensusRT> = consensus
        .iter()
        .map(|c| ((c.modified_sequence.as_str(), c.is_decoy), c))
        .collect();

    let mut stats_keep = 0usize;
    let mut stats_cwt = 0usize;
    let mut stats_forced = 0usize;
    let mut stats_no_consensus = 0usize;

    let empty_cwt: Vec<Vec<CwtCandidate>> = Vec::new();

    for (file_name, entries) in per_file_entries {
        // Use refined calibration if available, fall back to original
        let cal = per_file_refined_cal
            .get(file_name)
            .or_else(|| per_file_original_cal.get(file_name));

        let cal = match cal {
            Some(c) => c,
            None => continue,
        };

        let file_cwt = per_file_cwt_candidates.get(file_name).unwrap_or(&empty_cwt);

        for (entry_idx, entry) in entries.iter().enumerate() {
            let key = (entry.modified_sequence.as_str(), entry.is_decoy);
            let consensus_entry = match consensus_map.get(&key) {
                Some(c) => c,
                None => {
                    stats_no_consensus += 1;
                    continue;
                }
            };

            let expected_rt = cal.predict(consensus_entry.consensus_library_rt);
            let cwt = file_cwt.get(entry_idx).map(|v| v.as_slice()).unwrap_or(&[]);
            let action = determine_reconcile_action(
                entry.start_rt,
                entry.end_rt,
                cwt,
                expected_rt,
                consensus_entry.median_peak_width,
            );

            match &action {
                ReconcileAction::Keep => stats_keep += 1,
                ReconcileAction::UseCwtPeak { .. } => stats_cwt += 1,
                ReconcileAction::ForcedIntegration { .. } => stats_forced += 1,
            }

            // Only store non-Keep actions
            if !matches!(action, ReconcileAction::Keep) {
                actions.insert((file_name.clone(), entry_idx), action);
            }
        }
    }

    log::info!(
        "Inter-replicate reconciliation plan: {} keep, {} use alternate CWT peak, {} forced integration, {} not in consensus",
        stats_keep,
        stats_cwt,
        stats_forced,
        stats_no_consensus
    );

    actions
}

/// Compute weighted median of (value, weight) pairs.
///
/// Sorts by value, accumulates weights until cumulative >= total/2,
/// then returns the value at that point.
fn weighted_median(pairs: &[(f64, f64)]) -> f64 {
    if pairs.is_empty() {
        return 0.0;
    }
    if pairs.len() == 1 {
        return pairs[0].0;
    }

    let mut sorted: Vec<(f64, f64)> = pairs.to_vec();
    sorted.sort_by(|a, b| a.0.total_cmp(&b.0));

    let total_weight: f64 = sorted.iter().map(|(_, w)| w).sum();
    let half = total_weight / 2.0;

    let mut cumulative = 0.0;
    for &(value, weight) in &sorted {
        cumulative += weight;
        if cumulative >= half {
            return value;
        }
    }

    // Fallback (shouldn't reach here)
    sorted.last().unwrap().0
}

/// Simple median of a slice.
fn simple_median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weighted_median_single() {
        assert!((weighted_median(&[(5.0, 1.0)]) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_weighted_median_equal_weights() {
        // Equal weights → regular median
        let pairs = vec![(1.0, 1.0), (2.0, 1.0), (3.0, 1.0)];
        let result = weighted_median(&pairs);
        assert!(
            (result - 2.0).abs() < 1e-10,
            "Weighted median with equal weights should be regular median, got {}",
            result
        );
    }

    #[test]
    fn test_weighted_median_skewed_weights() {
        // Heavy weight on 10.0 should pull median toward 10.0
        let pairs = vec![(5.0, 1.0), (10.0, 100.0), (15.0, 1.0)];
        let result = weighted_median(&pairs);
        assert!(
            (result - 10.0).abs() < 1e-10,
            "Weighted median should be ~10.0 with heavy weight there, got {}",
            result
        );
    }

    #[test]
    fn test_simple_median_odd() {
        assert!((simple_median(&[3.0, 1.0, 2.0]) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_simple_median_even() {
        assert!((simple_median(&[1.0, 2.0, 3.0, 4.0]) - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_weighted_median_empty() {
        assert!((weighted_median(&[]) - 0.0).abs() < 1e-10);
    }
}
