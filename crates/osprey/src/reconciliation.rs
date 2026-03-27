//! Inter-replicate peak reconciliation for consistent peak integration across replicates.
//!
//! After initial per-run FDR, this module:
//! 1. Computes consensus library RTs for peptides detected across multiple runs
//! 2. Refits per-run LOESS calibration using consensus peptides
//! 3. Reconciles peak boundaries across runs (CWT peak overlap or forced integration)
//! 4. Runs a second FDR pass on the reconciled consensus set

use osprey_chromatography::RTCalibration;
use osprey_core::{CwtCandidate, FdrEntry};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

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
/// For each target peptide passing `consensus_fdr` at experiment level, and its
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
    // 1. Collect target peptides passing run-level FDR in at least one replicate
    let mut target_peptides: std::collections::HashSet<Arc<str>> = std::collections::HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if !entry.is_decoy && entry.run_qvalue <= consensus_fdr {
                target_peptides.insert(entry.modified_sequence.clone());
            }
        }
    }

    if target_peptides.is_empty() {
        return Vec::new();
    }

    log::info!(
        "Inter-replicate reconciliation: {} target peptides pass {:.1}% run-level FDR in at least one replicate",
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
    let mut decoy_peptides: std::collections::HashSet<Arc<str>> = std::collections::HashSet::new();
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
                decoy_peptides.contains(&*entry.modified_sequence)
            } else {
                target_peptides.contains(&*entry.modified_sequence)
            };

            if include {
                let peak_width = entry.end_rt - entry.start_rt;
                detections
                    .entry((entry.modified_sequence.to_string(), entry.is_decoy))
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
        if entry.is_decoy || entry.experiment_qvalue > consensus_fdr {
            continue;
        }
        if let Some(&consensus_lib_rt) = consensus_map.get(&*entry.modified_sequence) {
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
    experiment_fdr: f64,
) -> HashMap<(String, usize), ReconcileAction> {
    // Build consensus lookup by (modified_sequence, is_decoy)
    let consensus_map: HashMap<(&str, bool), &PeptideConsensusRT> = consensus
        .iter()
        .map(|c| ((c.modified_sequence.as_str(), c.is_decoy), c))
        .collect();

    // Build set of precursors (base_sequence, charge) that passed run-level FDR in any replicate.
    // Only these precursors (and their paired decoys) will be reconciled.
    // "base_sequence" strips the DECOY_ prefix so targets and decoys share the same key.
    let mut passing_precursors: std::collections::HashSet<(&str, u8)> =
        std::collections::HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if !entry.is_decoy && entry.run_qvalue <= experiment_fdr {
                passing_precursors.insert((&entry.modified_sequence, entry.charge));
            }
        }
    }

    let stats_keep_target = AtomicUsize::new(0);
    let stats_keep_decoy = AtomicUsize::new(0);
    let stats_cwt_target = AtomicUsize::new(0);
    let stats_cwt_decoy = AtomicUsize::new(0);
    let stats_forced_target = AtomicUsize::new(0);
    let stats_forced_decoy = AtomicUsize::new(0);

    let empty_cwt: Vec<Vec<CwtCandidate>> = Vec::new();

    // Process files in parallel — each file's entries are independent
    let actions: HashMap<(String, usize), ReconcileAction> = per_file_entries
        .par_iter()
        .flat_map(|(file_name, entries)| {
            // Use refined calibration if available, fall back to original
            let cal = per_file_refined_cal
                .get(file_name)
                .or_else(|| per_file_original_cal.get(file_name));

            let cal = match cal {
                Some(c) => c,
                None => return Vec::new(),
            };

            let file_cwt = per_file_cwt_candidates.get(file_name).unwrap_or(&empty_cwt);

            let mut file_actions = Vec::new();
            for (entry_idx, entry) in entries.iter().enumerate() {
                // Check if peptide is in consensus
                let key = (&*entry.modified_sequence, entry.is_decoy);
                let consensus_entry = match consensus_map.get(&key) {
                    Some(c) => c,
                    None => continue,
                };

                // Only reconcile precursors (peptide+charge) that passed experiment FDR,
                // or their paired decoys. This avoids wasting time on charge states that
                // weren't detected and won't appear in output.
                let base_seq = entry
                    .modified_sequence
                    .strip_prefix("DECOY_")
                    .unwrap_or(&entry.modified_sequence);
                if !passing_precursors.contains(&(base_seq, entry.charge)) {
                    continue;
                }

                let expected_rt = cal.predict(consensus_entry.consensus_library_rt);
                let cwt = file_cwt.get(entry_idx).map(|v| v.as_slice()).unwrap_or(&[]);
                let action = determine_reconcile_action(
                    entry.start_rt,
                    entry.end_rt,
                    cwt,
                    expected_rt,
                    consensus_entry.median_peak_width,
                );

                if entry.is_decoy {
                    match &action {
                        ReconcileAction::Keep => stats_keep_decoy.fetch_add(1, Ordering::Relaxed),
                        ReconcileAction::UseCwtPeak { .. } => {
                            stats_cwt_decoy.fetch_add(1, Ordering::Relaxed)
                        }
                        ReconcileAction::ForcedIntegration { .. } => {
                            stats_forced_decoy.fetch_add(1, Ordering::Relaxed)
                        }
                    };
                } else {
                    match &action {
                        ReconcileAction::Keep => stats_keep_target.fetch_add(1, Ordering::Relaxed),
                        ReconcileAction::UseCwtPeak { .. } => {
                            stats_cwt_target.fetch_add(1, Ordering::Relaxed)
                        }
                        ReconcileAction::ForcedIntegration { .. } => {
                            stats_forced_target.fetch_add(1, Ordering::Relaxed)
                        }
                    };
                }

                // Only store non-Keep actions
                if !matches!(action, ReconcileAction::Keep) {
                    file_actions.push(((file_name.clone(), entry_idx), action));
                }
            }
            file_actions
        })
        .collect();

    let n_keep_t = stats_keep_target.load(Ordering::Relaxed);
    let n_keep_d = stats_keep_decoy.load(Ordering::Relaxed);
    let n_cwt_t = stats_cwt_target.load(Ordering::Relaxed);
    let n_cwt_d = stats_cwt_decoy.load(Ordering::Relaxed);
    let n_forced_t = stats_forced_target.load(Ordering::Relaxed);
    let n_forced_d = stats_forced_decoy.load(Ordering::Relaxed);
    let n_keep = n_keep_t + n_keep_d;
    let n_cwt = n_cwt_t + n_cwt_d;
    let n_forced = n_forced_t + n_forced_d;
    let n_evaluated = n_keep + n_cwt + n_forced;
    let n_rescore = n_cwt + n_forced;
    let n_precursors = passing_precursors.len();
    let n_files = per_file_entries.len();
    log::info!(
        "Reconciliation: {} passing precursors × {} files = {} entries evaluated (across {} consensus peptides)",
        n_precursors,
        n_files,
        n_evaluated,
        consensus_map.len() / 2, // /2 because map has both target and decoy keys
    );
    log::info!(
        "Reconciliation plan: {} already correct ({} targets, {} decoys), {} need re-scoring ({} targets, {} decoys)",
        n_keep,
        n_keep_t,
        n_keep_d,
        n_rescore,
        n_cwt_t + n_forced_t,
        n_cwt_d + n_forced_d,
    );
    log::info!(
        "  Re-scoring breakdown: {} use alternate CWT peak, {} forced integration",
        n_cwt,
        n_forced,
    );

    actions
}

/// A precursor that needs gap-filling in a specific file.
///
/// Gap-fill targets are precursors that passed run-level FDR in at least one
/// replicate but were not detected in the initial search for a particular file. Both the target
/// and its paired decoy are included to maintain fair target-decoy competition.
#[derive(Debug, Clone)]
pub struct GapFillTarget {
    /// Library entry ID for the target
    pub target_entry_id: u32,
    /// Library entry ID for the paired decoy (target_entry_id | 0x80000000)
    pub decoy_entry_id: u32,
    /// Expected measured RT in this file (consensus library RT mapped through calibration)
    pub expected_rt: f64,
    /// Integration half-width (median_peak_width / 2)
    pub half_width: f64,
    /// Modified sequence (for logging/debugging)
    pub modified_sequence: Arc<str>,
    /// Charge state
    pub charge: u8,
}

/// Identify precursors that passed run-level FDR in any replicate but are missing from specific files.
///
/// For each file, finds passing precursors (target + paired decoy) that have no
/// entry in that file's `per_file_entries`. Returns gap-fill targets grouped by file,
/// with expected RTs computed from the consensus library RT mapped through each
/// file's calibration.
///
/// # Arguments
/// * `consensus` - Consensus RTs from `compute_consensus_rts()`
/// * `per_file_entries` - Per-file FdrEntry stubs (after first-pass FDR)
/// * `per_file_refined_cal` - Refined per-file calibrations (from consensus refit)
/// * `per_file_original_cal` - Original per-file calibrations
/// * `experiment_fdr` - Experiment-level FDR threshold
/// * `lib_lookup` - (modified_sequence, charge) → (target_entry_id, decoy_entry_id)
pub fn identify_gap_fill_targets(
    consensus: &[PeptideConsensusRT],
    per_file_entries: &[(String, Vec<FdrEntry>)],
    per_file_refined_cal: &HashMap<String, RTCalibration>,
    per_file_original_cal: &HashMap<String, RTCalibration>,
    experiment_fdr: f64,
    lib_lookup: &HashMap<(Arc<str>, u8), (u32, u32)>,
) -> HashMap<String, Vec<GapFillTarget>> {
    // 1. Build passing precursors: (modified_sequence, charge) for targets passing run-level
    //    FDR in at least one replicate. This ensures gap-fill covers all precursors that
    //    will be eligible for the second-pass experiment-level FDR.
    let mut passing_precursors: HashSet<(Arc<str>, u8)> = HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if !entry.is_decoy && entry.run_qvalue <= experiment_fdr {
                passing_precursors.insert((entry.modified_sequence.clone(), entry.charge));
            }
        }
    }

    if passing_precursors.is_empty() {
        return HashMap::new();
    }

    // 2. Build consensus lookup by (modified_sequence, is_decoy)
    let consensus_map: HashMap<(&str, bool), &PeptideConsensusRT> = consensus
        .iter()
        .map(|c| ((c.modified_sequence.as_str(), c.is_decoy), c))
        .collect();

    // 3. For each file, find missing precursors and create gap-fill targets
    let result: HashMap<String, Vec<GapFillTarget>> = per_file_entries
        .par_iter()
        .filter_map(|(file_name, entries)| {
            // Get calibration for this file
            let cal = per_file_refined_cal
                .get(file_name)
                .or_else(|| per_file_original_cal.get(file_name));
            let cal = match cal {
                Some(c) => c,
                None => return None,
            };

            // Build set of (modified_sequence, charge) present in this file (targets only)
            let present: HashSet<(Arc<str>, u8)> = entries
                .iter()
                .filter(|e| !e.is_decoy)
                .map(|e| (e.modified_sequence.clone(), e.charge))
                .collect();

            // Find missing precursors
            let mut targets: Vec<GapFillTarget> = Vec::new();
            for (mod_seq, charge) in &passing_precursors {
                if present.contains(&(mod_seq.clone(), *charge)) {
                    continue;
                }

                // Look up library entry IDs
                let (target_id, decoy_id) = match lib_lookup.get(&(mod_seq.clone(), *charge)) {
                    Some(ids) => *ids,
                    None => continue,
                };

                // Look up consensus RT for this target peptide
                let consensus_entry = match consensus_map.get(&(mod_seq.as_ref(), false)) {
                    Some(c) => c,
                    None => continue,
                };

                // Map consensus library RT to expected measured RT in this file
                let expected_rt = cal.predict(consensus_entry.consensus_library_rt);
                let half_width = consensus_entry.median_peak_width / 2.0;

                targets.push(GapFillTarget {
                    target_entry_id: target_id,
                    decoy_entry_id: decoy_id,
                    expected_rt,
                    half_width,
                    modified_sequence: mod_seq.clone(),
                    charge: *charge,
                });
            }

            if targets.is_empty() {
                None
            } else {
                Some((file_name.clone(), targets))
            }
        })
        .collect();

    // Log summary
    let total_targets: usize = result.values().map(|v| v.len()).sum();
    let n_files_with_gaps = result.len();
    if total_targets > 0 {
        log::info!(
            "Gap-fill: {} precursor-file pairs to score across {} files ({} passing precursors, {} files total)",
            total_targets,
            n_files_with_gaps,
            passing_precursors.len(),
            per_file_entries.len(),
        );
    }

    result
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
    use osprey_chromatography::RTCalibration;
    use osprey_core::CwtCandidate;

    // ---- Helper: create an identity RTCalibration (library RT = measured RT) ----
    fn make_identity_calibration() -> RTCalibration {
        let params = osprey_chromatography::RTModelParams {
            library_rts: vec![0.0, 10.0, 20.0, 30.0, 40.0],
            fitted_rts: vec![0.0, 10.0, 20.0, 30.0, 40.0],
            abs_residuals: vec![0.1, 0.1, 0.1, 0.1, 0.1],
        };
        RTCalibration::from_model_params(&params, 0.1).unwrap()
    }

    // ---- Helper: create an FdrEntry with key fields ----
    #[allow(clippy::too_many_arguments)]
    fn make_fdr_entry(
        entry_id: u32,
        modified_sequence: &str,
        charge: u8,
        is_decoy: bool,
        apex_rt: f64,
        start_rt: f64,
        end_rt: f64,
        coelution_sum: f64,
        experiment_qvalue: f64,
    ) -> FdrEntry {
        FdrEntry {
            entry_id,
            is_decoy,
            charge,
            scan_number: 1,
            apex_rt,
            start_rt,
            end_rt,
            coelution_sum,
            score: 0.0,
            run_qvalue: experiment_qvalue,
            experiment_qvalue,
            pep: 1.0,
            modified_sequence: Arc::from(modified_sequence),
        }
    }

    // ---- weighted_median / simple_median tests ----

    #[test]
    fn test_weighted_median_single() {
        assert!((weighted_median(&[(5.0, 1.0)]) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_weighted_median_equal_weights() {
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

    // ---- determine_reconcile_action tests ----

    #[test]
    fn test_reconcile_action_keep_when_peak_contains_expected_rt() {
        let action = determine_reconcile_action(
            9.0,  // start_rt
            11.0, // end_rt
            &[],  // no CWT candidates
            10.0, // expected RT is within [9, 11]
            1.0,  // median_peak_width
        );
        assert!(matches!(action, ReconcileAction::Keep));
    }

    #[test]
    fn test_reconcile_action_keep_at_boundary() {
        // Expected RT exactly at peak start
        let action = determine_reconcile_action(10.0, 12.0, &[], 10.0, 1.0);
        assert!(matches!(action, ReconcileAction::Keep));
        // Expected RT exactly at peak end
        let action = determine_reconcile_action(10.0, 12.0, &[], 12.0, 1.0);
        assert!(matches!(action, ReconcileAction::Keep));
    }

    #[test]
    fn test_reconcile_action_uses_cwt_candidate() {
        let cwt = vec![
            CwtCandidate {
                apex_rt: 15.0,
                start_rt: 14.0,
                end_rt: 16.0,
                area: 1000.0,
                snr: 5.0,
                coelution_score: 0.9,
            },
            CwtCandidate {
                apex_rt: 20.0,
                start_rt: 19.0,
                end_rt: 21.0,
                area: 500.0,
                snr: 3.0,
                coelution_score: 0.8,
            },
        ];
        // Current peak at 10 ± 1, expected RT at 15.0 → CWT candidate 0 overlaps
        let action = determine_reconcile_action(9.0, 11.0, &cwt, 15.0, 1.0);
        match action {
            ReconcileAction::UseCwtPeak {
                candidate_idx,
                apex_rt,
                ..
            } => {
                assert_eq!(candidate_idx, 0);
                assert!((apex_rt - 15.0).abs() < 1e-10);
            }
            other => panic!("Expected UseCwtPeak, got {:?}", other),
        }
    }

    #[test]
    fn test_reconcile_action_uses_second_cwt_candidate() {
        let cwt = vec![
            CwtCandidate {
                apex_rt: 5.0,
                start_rt: 4.0,
                end_rt: 6.0,
                area: 1000.0,
                snr: 5.0,
                coelution_score: 0.9,
            },
            CwtCandidate {
                apex_rt: 20.0,
                start_rt: 19.0,
                end_rt: 21.0,
                area: 500.0,
                snr: 3.0,
                coelution_score: 0.8,
            },
        ];
        // Expected RT at 20.0 → first CWT doesn't overlap, second does
        let action = determine_reconcile_action(9.0, 11.0, &cwt, 20.0, 1.0);
        match action {
            ReconcileAction::UseCwtPeak {
                candidate_idx,
                apex_rt,
                ..
            } => {
                assert_eq!(candidate_idx, 1);
                assert!((apex_rt - 20.0).abs() < 1e-10);
            }
            other => panic!("Expected UseCwtPeak(1), got {:?}", other),
        }
    }

    #[test]
    fn test_reconcile_action_forced_integration() {
        // No CWT candidates overlap expected RT
        let cwt = vec![CwtCandidate {
            apex_rt: 5.0,
            start_rt: 4.0,
            end_rt: 6.0,
            area: 1000.0,
            snr: 5.0,
            coelution_score: 0.9,
        }];
        let action = determine_reconcile_action(9.0, 11.0, &cwt, 20.0, 2.0);
        match action {
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                assert!((expected_rt - 20.0).abs() < 1e-10);
                assert!((half_width - 1.0).abs() < 1e-10); // median_peak_width/2 = 2.0/2
            }
            other => panic!("Expected ForcedIntegration, got {:?}", other),
        }
    }

    #[test]
    fn test_reconcile_action_forced_when_no_cwt() {
        let action = determine_reconcile_action(9.0, 11.0, &[], 20.0, 3.0);
        match action {
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                assert!((expected_rt - 20.0).abs() < 1e-10);
                assert!((half_width - 1.5).abs() < 1e-10); // 3.0/2
            }
            other => panic!("Expected ForcedIntegration, got {:?}", other),
        }
    }

    // ---- compute_consensus_rts tests ----

    #[test]
    fn test_consensus_rts_basic() {
        // Two files, both detecting peptide PEPTIDEK at slightly different RTs
        let cal = make_identity_calibration();
        let per_file_calibrations: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal),
        ]
        .into_iter()
        .collect();

        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![
                    make_fdr_entry(1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005),
                    make_fdr_entry(
                        1 | 0x80000000,
                        "DECOY_PEPTIDEK",
                        2,
                        true,
                        15.2,
                        14.7,
                        15.7,
                        3.0,
                        1.0,
                    ),
                ],
            ),
            (
                "file2".to_string(),
                vec![
                    make_fdr_entry(1, "PEPTIDEK", 2, false, 15.1, 14.6, 15.6, 7.5, 0.008),
                    make_fdr_entry(
                        1 | 0x80000000,
                        "DECOY_PEPTIDEK",
                        2,
                        true,
                        15.3,
                        14.8,
                        15.8,
                        2.5,
                        1.0,
                    ),
                ],
            ),
        ];

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01);

        // Should have consensus for target and decoy
        assert_eq!(consensus.len(), 2);

        let target = consensus.iter().find(|c| !c.is_decoy).unwrap();
        assert_eq!(target.modified_sequence, "PEPTIDEK");
        assert_eq!(target.n_runs_detected, 2);
        // With identity calibration, consensus library RT ≈ weighted median of apex RTs
        assert!(
            (target.consensus_library_rt - 15.0).abs() < 0.5,
            "Consensus RT {} should be near 15.0",
            target.consensus_library_rt
        );
    }

    #[test]
    fn test_consensus_rts_excludes_non_passing() {
        let cal = make_identity_calibration();
        let per_file_calibrations: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let per_file_entries = vec![(
            "file1".to_string(),
            vec![
                // This entry has experiment_qvalue > 0.01, should not contribute
                make_fdr_entry(1, "FAILPEPTIDE", 2, false, 10.0, 9.5, 10.5, 5.0, 0.05),
            ],
        )];

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01);
        assert!(
            consensus.is_empty(),
            "Non-passing peptides should be excluded"
        );
    }

    #[test]
    fn test_consensus_rts_outlier_robustness() {
        // 3 files: 2 agree on RT~15.0, 1 outlier at RT~30.0
        // Weighted median should be robust to the outlier
        let cal = make_identity_calibration();
        let per_file_calibrations: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
            ("file3".to_string(), cal),
        ]
        .into_iter()
        .collect();

        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005,
                )],
            ),
            (
                "file2".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.1, 14.6, 15.6, 7.5, 0.008,
                )],
            ),
            (
                "file3".to_string(),
                vec![
                    // Outlier: wrong peak at RT 30.0 with low coelution_sum
                    make_fdr_entry(1, "PEPTIDEK", 2, false, 30.0, 29.5, 30.5, 3.0, 0.009),
                ],
            ),
        ];

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01);
        let target = consensus.iter().find(|c| !c.is_decoy).unwrap();
        // Weighted median should be near 15.0, not pulled to 30.0
        assert!(
            (target.consensus_library_rt - 15.0).abs() < 1.0,
            "Consensus RT {} should be near 15.0, not pulled by outlier",
            target.consensus_library_rt
        );
    }

    // ---- plan_reconciliation tests ----

    #[test]
    fn test_plan_reconciliation_keep_correct_peaks() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 1.0,
            n_runs_detected: 2,
        }];

        // Entry already at correct RT
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005,
            )],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> =
            vec![("file1".to_string(), vec![vec![]])]
                .into_iter()
                .collect();

        let actions = plan_reconciliation(
            &consensus,
            &per_file_entries,
            &cwt,
            &refined_cal,
            &original_cal,
            0.01,
        );

        // Entry at correct peak → Keep → not in the map
        assert!(
            actions.is_empty(),
            "Entry at correct peak should not appear in actions"
        );
    }

    #[test]
    fn test_plan_reconciliation_forced_integration_wrong_peak() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 20.0, // consensus at RT 20
            median_peak_width: 2.0,
            n_runs_detected: 2,
        }];

        // Entry at wrong RT (10.0 instead of 20.0), no CWT candidates
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 10.0, 9.5, 10.5, 5.0, 0.005,
            )],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> =
            vec![("file1".to_string(), vec![vec![]])]
                .into_iter()
                .collect();

        let actions = plan_reconciliation(
            &consensus,
            &per_file_entries,
            &cwt,
            &refined_cal,
            &original_cal,
            0.01,
        );

        let action = actions.get(&("file1".to_string(), 0)).unwrap();
        match action {
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                assert!(
                    (*expected_rt - 20.0).abs() < 0.5,
                    "Expected RT {} should be near 20.0",
                    expected_rt
                );
                assert!(
                    (*half_width - 1.0).abs() < 1e-10,
                    "Half-width should be 1.0 (median_peak_width/2)"
                );
            }
            other => panic!("Expected ForcedIntegration, got {:?}", other),
        }
    }

    #[test]
    fn test_plan_reconciliation_uses_cwt_candidate() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 20.0,
            median_peak_width: 2.0,
            n_runs_detected: 2,
        }];

        // Entry at wrong RT, but CWT candidate at correct RT
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 10.0, 9.5, 10.5, 5.0, 0.005,
            )],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> = vec![(
            "file1".to_string(),
            vec![vec![CwtCandidate {
                apex_rt: 20.0,
                start_rt: 19.0,
                end_rt: 21.0,
                area: 800.0,
                snr: 4.0,
                coelution_score: 0.85,
            }]],
        )]
        .into_iter()
        .collect();

        let actions = plan_reconciliation(
            &consensus,
            &per_file_entries,
            &cwt,
            &refined_cal,
            &original_cal,
            0.01,
        );

        let action = actions.get(&("file1".to_string(), 0)).unwrap();
        match action {
            ReconcileAction::UseCwtPeak { apex_rt, .. } => {
                assert!(
                    (*apex_rt - 20.0).abs() < 1e-10,
                    "Should use CWT candidate at RT 20.0"
                );
            }
            other => panic!("Expected UseCwtPeak, got {:?}", other),
        }
    }

    #[test]
    fn test_plan_reconciliation_skips_non_passing_precursors() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "FAILPEPTIDE".to_string(),
            is_decoy: false,
            consensus_library_rt: 20.0,
            median_peak_width: 2.0,
            n_runs_detected: 2,
        }];

        // Entry with experiment_qvalue > 0.01 (doesn't pass)
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1,
                "FAILPEPTIDE",
                2,
                false,
                10.0,
                9.5,
                10.5,
                5.0,
                0.05, // fails FDR
            )],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> =
            vec![("file1".to_string(), vec![vec![]])]
                .into_iter()
                .collect();

        let actions = plan_reconciliation(
            &consensus,
            &per_file_entries,
            &cwt,
            &refined_cal,
            &original_cal,
            0.01,
        );

        assert!(
            actions.is_empty(),
            "Non-passing precursors should not be reconciled"
        );
    }

    // ---- identify_gap_fill_targets tests ----

    #[test]
    fn test_gap_fill_identifies_missing_precursors() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![
            PeptideConsensusRT {
                modified_sequence: "PEPTIDEK".to_string(),
                is_decoy: false,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 2,
            },
            PeptideConsensusRT {
                modified_sequence: "DECOY_PEPTIDEK".to_string(),
                is_decoy: true,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 2,
            },
        ];

        // file1 has PEPTIDEK, file2 does NOT
        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![
                    make_fdr_entry(1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005),
                    make_fdr_entry(
                        1 | 0x80000000,
                        "DECOY_PEPTIDEK",
                        2,
                        true,
                        15.2,
                        14.7,
                        15.7,
                        3.0,
                        1.0,
                    ),
                ],
            ),
            (
                "file2".to_string(),
                vec![
                    // file2 has a DIFFERENT peptide, not PEPTIDEK
                    make_fdr_entry(2, "OTHERPEPTIDE", 2, false, 20.0, 19.5, 20.5, 6.0, 0.003),
                ],
            ),
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        // file1 has PEPTIDEK → no gap-fill needed
        assert!(
            !gap_fill.contains_key("file1"),
            "file1 already has PEPTIDEK, no gap-fill needed"
        );

        // file2 is missing PEPTIDEK → should have a gap-fill target
        let file2_targets = gap_fill
            .get("file2")
            .expect("file2 should have gap-fill targets");
        assert_eq!(file2_targets.len(), 1);
        assert_eq!(file2_targets[0].target_entry_id, 1);
        assert_eq!(file2_targets[0].decoy_entry_id, 1 | 0x80000000);
        assert_eq!(&*file2_targets[0].modified_sequence, "PEPTIDEK");
        assert_eq!(file2_targets[0].charge, 2);
        // With identity calibration, expected_rt should be near consensus library RT
        assert!(
            (file2_targets[0].expected_rt - 15.0).abs() < 0.5,
            "Expected RT {} should be near consensus 15.0",
            file2_targets[0].expected_rt
        );
        assert!(
            (file2_targets[0].half_width - 0.5).abs() < 1e-10,
            "Half-width should be 0.5 (1.0/2)"
        );
    }

    #[test]
    fn test_gap_fill_no_targets_when_all_present() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 1.0,
            n_runs_detected: 2,
        }];

        // Both files have PEPTIDEK
        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005,
                )],
            ),
            (
                "file2".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.1, 14.6, 15.6, 7.5, 0.008,
                )],
            ),
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        assert!(
            gap_fill.is_empty(),
            "No gap-fill needed when all files have the precursor"
        );
    }

    #[test]
    fn test_gap_fill_excludes_non_passing_precursors() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 1.0,
            n_runs_detected: 1,
        }];

        // file1 has PEPTIDEK but it FAILS experiment-level FDR
        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.05, // fails
                )],
            ),
            ("file2".to_string(), vec![]),
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        assert!(
            gap_fill.is_empty(),
            "Non-passing precursors should not generate gap-fill targets"
        );
    }

    #[test]
    fn test_gap_fill_skips_missing_library_entry() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 1.0,
            n_runs_detected: 2,
        }];

        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005,
                )],
            ),
            ("file2".to_string(), vec![]),
        ];

        // Empty lib_lookup — PEPTIDEK not in library
        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> = HashMap::new();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        assert!(
            gap_fill.is_empty(),
            "Should skip precursors not in library lookup"
        );
    }

    #[test]
    fn test_gap_fill_multiple_precursors_multiple_files() {
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
            ("file3".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![
            PeptideConsensusRT {
                modified_sequence: "PEPTIDEK".to_string(),
                is_decoy: false,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 2,
            },
            PeptideConsensusRT {
                modified_sequence: "ANOTHERPEP".to_string(),
                is_decoy: false,
                consensus_library_rt: 25.0,
                median_peak_width: 1.5,
                n_runs_detected: 2,
            },
        ];

        // file1: has both peptides
        // file2: has PEPTIDEK only (missing ANOTHERPEP)
        // file3: has ANOTHERPEP only (missing PEPTIDEK)
        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![
                    make_fdr_entry(1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005),
                    make_fdr_entry(2, "ANOTHERPEP", 2, false, 25.0, 24.5, 25.5, 7.0, 0.006),
                ],
            ),
            (
                "file2".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.1, 14.6, 15.6, 7.5, 0.008,
                )],
            ),
            (
                "file3".to_string(),
                vec![make_fdr_entry(
                    2,
                    "ANOTHERPEP",
                    2,
                    false,
                    25.1,
                    24.6,
                    25.6,
                    6.5,
                    0.007,
                )],
            ),
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> = vec![
            ((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000)),
            ((Arc::from("ANOTHERPEP"), 2u8), (2u32, 2 | 0x80000000)),
        ]
        .into_iter()
        .collect();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        // file1: has both → no gap-fill
        assert!(!gap_fill.contains_key("file1"));

        // file2: missing ANOTHERPEP → 1 gap-fill target
        let f2 = gap_fill.get("file2").expect("file2 should have gap-fill");
        assert_eq!(f2.len(), 1);
        assert_eq!(&*f2[0].modified_sequence, "ANOTHERPEP");
        assert_eq!(f2[0].target_entry_id, 2);

        // file3: missing PEPTIDEK → 1 gap-fill target
        let f3 = gap_fill.get("file3").expect("file3 should have gap-fill");
        assert_eq!(f3.len(), 1);
        assert_eq!(&*f3[0].modified_sequence, "PEPTIDEK");
        assert_eq!(f3[0].target_entry_id, 1);
    }

    #[test]
    fn test_gap_fill_skips_file_without_calibration() {
        let cal = make_identity_calibration();
        // Only file1 has calibration, file2 does not
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 1.0,
            n_runs_detected: 2,
        }];

        let per_file_entries = vec![
            (
                "file1".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005,
                )],
            ),
            ("file2".to_string(), vec![]), // missing PEPTIDEK AND no calibration
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        // file2 has no calibration → should be skipped entirely
        assert!(
            !gap_fill.contains_key("file2"),
            "Files without calibration should be skipped"
        );
    }

    #[test]
    fn test_gap_fill_single_file_no_gaps() {
        // Single file experiment → no gap-fill possible
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> = refined_cal.clone();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 1.0,
            n_runs_detected: 1,
        }];

        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.0, 14.5, 15.5, 8.0, 0.005,
            )],
        )];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        let gap_fill = identify_gap_fill_targets(
            &consensus,
            &per_file_entries,
            &refined_cal,
            &original_cal,
            0.01,
            &lib_lookup,
        );

        assert!(
            gap_fill.is_empty(),
            "Single file with the precursor → no gaps"
        );
    }
}
