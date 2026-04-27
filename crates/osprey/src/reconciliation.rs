//! Inter-replicate peak reconciliation for consistent peak integration across replicates.
//!
//! After initial per-run FDR, this module:
//! 1. Computes consensus library RTs for peptides detected across multiple runs
//! 2. Refits per-run LOESS calibration using consensus peptides
//! 3. Reconciles peak boundaries across runs (CWT peak overlap or forced integration)
//! 4. Runs a second FDR pass on the reconciled consensus set

use osprey_chromatography::RTCalibration;
use osprey_core::{CwtCandidate, FdrEntry, FdrLevel};
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
    /// MAD of this peptide's apex RTs across runs, in library RT space (minutes).
    /// `None` when fewer than 3 detections contribute (insufficient to estimate).
    /// Captures within-peptide RT reproducibility, typically 3-5x tighter than
    /// the cross-peptide calibration MAD. Used by `plan_reconciliation` as a
    /// peptide-specific RT tolerance.
    pub apex_library_rt_mad: Option<f64>,
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
/// * `protein_fdr_threshold` - If > 0, rescue borderline peptides whose first-pass
///   protein q-value is <= this threshold. Lets peptides from strong proteins
///   contribute to consensus RT computation even if their own peptide q-value is
///   borderline. Typically set to `config.protein_fdr`. Pass 0.0 to disable.
pub fn compute_consensus_rts(
    per_file_entries: &[(String, Vec<FdrEntry>)],
    per_file_calibrations: &HashMap<String, RTCalibration>,
    consensus_fdr: f64,
    protein_fdr_threshold: f64,
) -> Vec<PeptideConsensusRT> {
    // A target detection qualifies for consensus RT computation if:
    //   (1) Its per-entry precursor q-value is within consensus_fdr — this is
    //       a HARD gate. The detection's RT is only trustworthy if the
    //       detection itself has decent direct evidence. Protein-level FDR
    //       cannot rescue poor precursor-level evidence, because the
    //       consensus RT is driven by each entry's own apex_rt.
    //   (2) AND either the peptide-level q-value also passes, OR the first-
    //       pass protein group passes (the borderline-peptide rescue).
    //
    // The earlier looser gate allowed protein FDR to rescue entries with
    // poor precursor-level evidence. That was observed to poison the
    // consensus when a strong protein had a wrong-peak detection in some
    // charge state: the low-quality wrong-peak apex got pulled into the
    // weighted median and shifted the consensus to the interferer. The
    // tighter gate rejects such entries and lets the remaining good
    // detections anchor the consensus.
    let qualifies = |entry: &FdrEntry| -> bool {
        if entry.is_decoy {
            return false;
        }
        if entry.run_precursor_qvalue > consensus_fdr {
            return false;
        }
        entry.run_peptide_qvalue <= consensus_fdr
            || (protein_fdr_threshold > 0.0 && entry.run_protein_qvalue <= protein_fdr_threshold)
    };

    // 1. Collect target peptides passing run-level FDR (or rescued by protein FDR)
    let mut target_peptides: std::collections::HashSet<Arc<str>> = std::collections::HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if qualifies(entry) {
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
    //    Group by (modified_sequence, is_decoy).
    //    Each detection: (file_name, apex_rt, score, peak_width, coelution_sum).
    //    `score` is the SVM discriminant (signed) — used for the consensus
    //    weight via `sigmoid(score)`. `coelution_sum` is retained only as a
    //    sanity filter (> 0.0 → fragments correlate positively).
    type DetectionVec = Vec<(String, f64, f64, f64, f64)>;
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

    // Collect detections for consensus peptides.
    // For targets: only include detections that qualify for consensus (pass peptide
    // FDR or rescued by protein FDR). This ensures the consensus RT and peak width
    // are computed from actual good detections, not noise peaks in replicates where
    // the peptide wasn't confidently found.
    // For decoys: include all (paired to passing targets for reconciliation planning).
    for (file_name, entries) in per_file_entries {
        for entry in entries {
            let include = if entry.is_decoy {
                decoy_peptides.contains(&*entry.modified_sequence)
            } else {
                target_peptides.contains(&*entry.modified_sequence) && qualifies(entry)
            };

            if include {
                let peak_width = entry.end_rt - entry.start_rt;
                detections
                    .entry((entry.modified_sequence.to_string(), entry.is_decoy))
                    .or_default()
                    .push((
                        file_name.clone(),
                        entry.apex_rt,
                        entry.score,
                        peak_width,
                        entry.coelution_sum,
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
        let mut peak_width_weights: Vec<(f64, f64)> = Vec::new();

        for (file_name, apex_rt, score, peak_width, coelution_sum) in dets {
            if let Some(cal) = per_file_calibrations.get(file_name) {
                let library_rt = cal.inverse_predict(*apex_rt);
                // Sanity filter: require at least minimally positive fragment
                // co-elution. This rejects anti-correlated noise integrations
                // (e.g., forced integrations at empty RT regions).
                if library_rt.is_finite() && *coelution_sum > 0.0 {
                    // Weight by sigmoid(SVM score) — the SVM discriminant is
                    // our strongest quality signal, and sigmoid maps it to
                    // (0, 1) monotonically. A wrong-peak detection with
                    // negative score gets near-zero weight and cannot poison
                    // the weighted median; a strong detection with positive
                    // score dominates. Floor at 1e-6 so every detection keeps
                    // a non-zero weight (avoids degenerate zero-total-weight
                    // edge cases when all scores are very negative).
                    let weight = (1.0 / (1.0 + (-*score).exp())).max(1e-6);
                    library_rt_weights.push((library_rt, weight));
                    peak_width_weights.push((*peak_width, weight));
                }
            }
        }

        if library_rt_weights.is_empty() {
            continue;
        }

        let consensus_library_rt = weighted_median(&library_rt_weights);
        let median_peak_width = weighted_median(&peak_width_weights);
        let n_runs_detected = library_rt_weights.len();

        // Within-peptide RT MAD in library space. Requires >= 3 detections
        // for a stable estimate (with only 2 points, MAD is just half the
        // range and not robust). Captures how tightly this peptide reproduces
        // across runs — typically much tighter than the cross-peptide
        // calibration MAD, so this gives a better RT tolerance downstream.
        let apex_library_rt_mad = if n_runs_detected >= 3 {
            let mut abs_devs: Vec<f64> = library_rt_weights
                .iter()
                .map(|(lib_rt, _)| (lib_rt - consensus_library_rt).abs())
                .collect();
            abs_devs.sort_by(|a, b| a.total_cmp(b));
            let mid = abs_devs.len() / 2;
            let mad = if abs_devs.len() % 2 == 0 {
                0.5 * (abs_devs[mid - 1] + abs_devs[mid])
            } else {
                abs_devs[mid]
            };
            Some(mad)
        } else {
            None
        };

        if crate::trace::is_traced(modified_sequence) {
            let kind = if *is_decoy { "DECOY" } else { "TARGET" };
            log::info!(
                "[trace] consensus {} {} — detections ({} runs contributing):",
                modified_sequence,
                kind,
                n_runs_detected,
            );
            for (file_name, apex_rt, score, peak_width, coelution_sum) in dets {
                let lib_rt_str = per_file_calibrations
                    .get(file_name)
                    .map(|cal| format!("lib_rt={:.4}", cal.inverse_predict(*apex_rt)))
                    .unwrap_or_else(|| "lib_rt=<no-cal>".into());
                let weight = (1.0 / (1.0 + (-*score).exp())).max(1e-6);
                log::info!(
                    "[trace]   file={} apex_rt={:.4} {} score={:.3} weight={:.4} coelution_sum={:.3} peak_width={:.3}",
                    file_name,
                    apex_rt,
                    lib_rt_str,
                    score,
                    weight,
                    coelution_sum,
                    peak_width,
                );
            }
            log::info!(
                "[trace] consensus {} {} → consensus_library_rt={:.4}, median_peak_width={:.3}",
                modified_sequence,
                kind,
                consensus_library_rt,
                median_peak_width,
            );
        }

        consensus.push(PeptideConsensusRT {
            modified_sequence: modified_sequence.clone(),
            is_decoy: *is_decoy,
            consensus_library_rt,
            median_peak_width,
            n_runs_detected,
            apex_library_rt_mad,
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
        if entry.is_decoy || entry.effective_experiment_qvalue(FdrLevel::Both) > consensus_fdr {
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
    //
    // classical_robust_iterations must mirror the value used at Stage 4 on
    // both Rust and C# sides; otherwise the refit fitted_values diverge
    // cross-impl. Stage 4 reads OSPREY_LOESS_CLASSICAL_ROBUST (default = on,
    // see pipeline.rs); read it the same way here.
    let classical_robust = std::env::var("OSPREY_LOESS_CLASSICAL_ROBUST").as_deref() != Ok("0");
    let calibrator = RTCalibrator::with_config(osprey_chromatography::RTCalibratorConfig {
        bandwidth: 0.3,
        outlier_retention: 1.0,
        classical_robust_iterations: classical_robust,
        ..osprey_chromatography::RTCalibratorConfig::default()
    });
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
/// Uses **apex proximity** to decide whether the current peak is at the right RT.
/// The tolerance comes from the refined per-file calibration's local residual,
/// which reflects the actual run-to-run RT variation at that gradient position
/// (derived from thousands of consensus peptides, not just this one).
///
/// Previous versions checked boundary containment (`start_rt <= expected <= end_rt`),
/// which incorrectly classified peaks with wrong apex but wide tails as "Keep".
pub fn determine_reconcile_action(
    apex_rt: f64,
    cwt_candidates: &[CwtCandidate],
    expected_measured_rt: f64,
    rt_tolerance: f64,
    half_width: f64,
) -> ReconcileAction {
    // Check if the current peak's apex is within tolerance of the expected RT
    if (apex_rt - expected_measured_rt).abs() <= rt_tolerance {
        return ReconcileAction::Keep;
    }

    // Check CWT candidates by apex proximity (not boundary overlap).
    // Pick the candidate whose apex is closest to expected_measured_rt, if within tolerance.
    let best_cwt = cwt_candidates
        .iter()
        .enumerate()
        .filter(|(_, cand)| (cand.apex_rt - expected_measured_rt).abs() <= rt_tolerance)
        .min_by(|(_, a), (_, b)| {
            let da = (a.apex_rt - expected_measured_rt).abs();
            let db = (b.apex_rt - expected_measured_rt).abs();
            da.total_cmp(&db)
        });

    if let Some((idx, cand)) = best_cwt {
        return ReconcileAction::UseCwtPeak {
            candidate_idx: idx,
            start_rt: cand.start_rt,
            apex_rt: cand.apex_rt,
            end_rt: cand.end_rt,
        };
    }

    // No peak with apex near expected RT — forced integration
    ReconcileAction::ForcedIntegration {
        expected_rt: expected_measured_rt,
        half_width,
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

    // Global within-peptide RT MAD (library RT space) — median of per-peptide
    // apex MADs across all target peptides that have >= 3 detections.
    //
    // After cross-run RT alignment (LOESS refit on consensus peptides), the
    // *within-peptide* variance is what matters for reconciliation tolerance.
    // We expect it to be roughly peptide-independent: once aligned, the scatter
    // of a well-behaved peptide's apex RT across replicates reflects the
    // chromatographic reproducibility of the instrument, not anything unique
    // to that peptide. The median across peptides is a far more stable
    // estimator than any single peptide's MAD (which is very noisy with only
    // 3-5 replicates).
    //
    // This replaces the earlier cross-peptide calibration MAD, which conflated
    // LOESS fit residuals (peptide-to-peptide variation) with true run-to-run
    // RT reproducibility. Typical global within-peptide MAD is 3-5x tighter
    // than calibration MAD, giving a correspondingly tighter rt_tolerance.
    //
    // If no peptides have >= 3 detections (e.g., 2-replicate experiment), fall
    // back to a conservative 0.05 min.
    let global_within_peptide_mad_lib = {
        let mut peptide_mads: Vec<f64> = consensus
            .iter()
            .filter(|c| !c.is_decoy)
            .filter_map(|c| c.apex_library_rt_mad)
            .collect();
        if peptide_mads.is_empty() {
            log::warn!(
                "Reconciliation: no peptides have >= 3 detections; using 0.05 min fallback MAD"
            );
            0.05
        } else {
            peptide_mads.sort_by(|a, b| a.total_cmp(b));
            let mid = peptide_mads.len() / 2;
            let med = if peptide_mads.len() % 2 == 0 {
                0.5 * (peptide_mads[mid - 1] + peptide_mads[mid])
            } else {
                peptide_mads[mid]
            };
            log::info!(
                "Reconciliation: global within-peptide RT MAD = {:.4} min ({} peptides with >= 3 detections contributed)",
                med,
                peptide_mads.len()
            );
            med
        }
    };

    // Build set of precursors (base_sequence, charge) that will appear in
    // the blib output. Every such precursor needs reconciliation in ALL
    // files so its per-file boundaries agree with the consensus. We use the
    // MINIMUM of (run_precursor, run_peptide, experiment_precursor,
    // experiment_peptide) q-values — most permissive — because the blib
    // admits a precursor if ANY of those levels passes at the configured
    // FDR gate. Using the strict max (e.g., `FdrLevel::Both`) would leave
    // out precursors that clear peptide-level FDR but just barely fail
    // precursor-level; those precursors still ship in the blib and must be
    // reconciled so sibling charges don't end up at inconsistent RTs within
    // the same run. Decoys are included through the paired-decoy logic
    // below, not this gate. "base_sequence" strips the DECOY_ prefix so
    // targets and decoys share the same key.
    let mut passing_precursors: std::collections::HashSet<(&str, u8)> =
        std::collections::HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if entry.is_decoy {
                continue;
            }
            let best_q = entry
                .run_precursor_qvalue
                .min(entry.run_peptide_qvalue)
                .min(entry.experiment_precursor_qvalue)
                .min(entry.experiment_peptide_qvalue);
            if best_q <= experiment_fdr {
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

            // Compute a per-file safety ceiling for RT tolerance from the
            // refined calibration's residuals. This is a cross-peptide measure
            // and is usually LOOSER than any individual peptide's within-peptide
            // reproducibility — we use it only as an upper bound so that a
            // peptide with pathologically low apparent MAD (e.g., n=3 with all
            // detections at exactly the same RT) cannot produce a sub-floor
            // tolerance that rejects every tiny drift.
            //
            // The per-entry tolerance itself is derived from within-peptide MAD
            // (see the inner loop below).
            //
            // Sigma-clipped MAD: the raw MAD is inflated when wrong-peak
            // detections feed their bad apex_rt into the refit calibration. We
            // clip at 3-sigma from the raw MAD and recompute MAD from the
            // survivors. Capped at the original (first-pass) calibration
            // tolerance so each successive calibration pass can only tighten.
            let file_cal_tolerance_ceiling = {
                let raw_mad = cal.stats().mad;
                let clip_threshold = raw_mad * 1.4826 * 3.0;
                let clipped_mad = sigma_clipped_mad(cal.abs_residuals(), clip_threshold);
                let refined_tolerance = (clipped_mad * 1.4826 * 3.0).max(0.1);

                let original_cal_for_cap = per_file_original_cal.get(file_name);
                let cap = original_cal_for_cap
                    .map(|oc| (oc.stats().mad * 1.4826 * 3.0).max(0.1))
                    .unwrap_or(refined_tolerance);

                let final_tol = refined_tolerance.min(cap);
                log::debug!(
                    "Reconciliation RT tolerance ceiling for {}: {:.3} min (raw_MAD={:.3}, clipped_MAD={:.3}, original_cap={:.3})",
                    file_name,
                    final_tol,
                    raw_mad,
                    clipped_mad,
                    cap,
                );
                final_tol
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
                let cwt = file_cwt
                    .get(entry.parquet_index as usize)
                    .map(|v| v.as_slice())
                    .unwrap_or(&[]);

                // RT tolerance from the global within-peptide MAD (library RT
                // space, treated as approximately equal to measured space given
                // local calibration slope ~ 1). We use the *global* median
                // rather than this peptide's individual MAD because:
                //
                //   (a) After alignment, within-peptide scatter is roughly
                //       peptide-independent — it reflects instrument/LC
                //       reproducibility, not anything unique to the peptide.
                //   (b) An individual peptide's MAD from 3-5 replicates is a
                //       very noisy estimator; the cross-peptide median is far
                //       more stable.
                //
                // The 3 × 1.4826 factor converts MAD to a ~3-sigma tolerance.
                // Floor of 0.1 min accommodates scan-resolution rounding.
                // Ceiling is the file's cross-peptide calibration tolerance as
                // a safety net against pathological MAD estimates.
                let peptide_tolerance = (global_within_peptide_mad_lib * 1.4826 * 3.0)
                    .max(0.1)
                    .min(file_cal_tolerance_ceiling);

                let action = determine_reconcile_action(
                    entry.apex_rt,
                    cwt,
                    expected_rt,
                    peptide_tolerance,
                    consensus_entry.median_peak_width / 2.0,
                );

                if crate::trace::is_traced(&entry.modified_sequence) {
                    let kind = if entry.is_decoy { "DECOY" } else { "TARGET" };
                    let action_str = match &action {
                        ReconcileAction::Keep => "Keep".to_string(),
                        ReconcileAction::UseCwtPeak {
                            candidate_idx,
                            apex_rt,
                            ..
                        } => format!("UseCwtPeak(idx={}, apex={:.4})", candidate_idx, apex_rt),
                        ReconcileAction::ForcedIntegration {
                            expected_rt,
                            half_width,
                        } => format!(
                            "ForcedIntegration(center={:.4}, half_width={:.3})",
                            expected_rt, half_width
                        ),
                    };
                    let own_mad = consensus_entry
                        .apex_library_rt_mad
                        .map(|m| format!("{:.4}", m))
                        .unwrap_or_else(|| "n/a".into());
                    log::info!(
                        "[trace] plan {} z={} {} file={} current_apex={:.4} expected_rt={:.4} residual={:+.3} tolerance={:.3} (global_MAD={:.4}, this_peptide_MAD={} from n={}, ceiling={:.3}) n_cwt={} → {}",
                        entry.modified_sequence,
                        entry.charge,
                        kind,
                        file_name,
                        entry.apex_rt,
                        expected_rt,
                        entry.apex_rt - expected_rt,
                        peptide_tolerance,
                        global_within_peptide_mad_lib,
                        own_mad,
                        consensus_entry.n_runs_detected,
                        file_cal_tolerance_ceiling,
                        cwt.len(),
                        action_str,
                    );
                    for (i, c) in cwt.iter().enumerate() {
                        log::info!(
                            "[trace]   cwt[{}]: apex={:.4} (res={:+.3}) coelution={:.4} area={:.1}",
                            i,
                            c.apex_rt,
                            c.apex_rt - expected_rt,
                            c.coelution_score,
                            c.area,
                        );
                    }
                }

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
#[allow(clippy::too_many_arguments)]
pub fn identify_gap_fill_targets(
    consensus: &[PeptideConsensusRT],
    per_file_entries: &[(String, Vec<FdrEntry>)],
    per_file_refined_cal: &HashMap<String, RTCalibration>,
    per_file_original_cal: &HashMap<String, RTCalibration>,
    experiment_fdr: f64,
    lib_lookup: &HashMap<(Arc<str>, u8), (u32, u32)>,
    lib_precursor_mz: &HashMap<u32, f64>,
    per_file_isolation_mz: &HashMap<String, Vec<(f64, f64)>>,
) -> HashMap<String, Vec<GapFillTarget>> {
    // 1. Build passing precursors: (modified_sequence, charge) for targets
    //    that will end up in the blib output. See the matching rationale in
    //    `plan_reconciliation` — we use the minimum of (run, experiment) ×
    //    (precursor, peptide) q-values so a peptide that clears peptide-
    //    level FDR (even with middling precursor-level q) still gets its
    //    missing charge states gap-filled to the correct consensus RT.
    let mut passing_precursors: HashSet<(Arc<str>, u8)> = HashSet::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            if entry.is_decoy {
                continue;
            }
            let best_q = entry
                .run_precursor_qvalue
                .min(entry.run_peptide_qvalue)
                .min(entry.experiment_precursor_qvalue)
                .min(entry.experiment_peptide_qvalue);
            if best_q <= experiment_fdr {
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
    let filtered_by_mz = AtomicUsize::new(0);
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

            // This file's isolation window m/z coverage (list of [lower, upper] intervals).
            // When available, gap-fill candidates are filtered to only those whose
            // library precursor m/z falls inside at least one window — otherwise
            // forcing an integration in this file would land on noise (the precursor's
            // fragments aren't being selected for MS2 in this m/z range). This is the
            // general fix for GPF datasets where each file covers a disjoint m/z range,
            // and also helps normal multi-replicate acquisitions with partial overlap.
            let iso_windows = per_file_isolation_mz.get(file_name);

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

                // m/z range filter: skip precursors whose m/z is not covered by any
                // isolation window in this file. Without this, GPF runs would force
                // integrations at precursors that this file physically cannot observe.
                if let Some(windows) = iso_windows {
                    let precursor_mz = match lib_precursor_mz.get(&target_id) {
                        Some(&mz) => mz,
                        None => continue,
                    };
                    let in_range = windows
                        .iter()
                        .any(|&(lo, hi)| precursor_mz >= lo && precursor_mz < hi);
                    if !in_range {
                        filtered_by_mz.fetch_add(1, Ordering::Relaxed);
                        continue;
                    }
                }

                // Look up consensus RT for this target peptide
                let consensus_entry = match consensus_map.get(&(mod_seq.as_ref(), false)) {
                    Some(c) => c,
                    None => continue,
                };

                // Map consensus library RT to expected measured RT in this file
                let expected_rt = cal.predict(consensus_entry.consensus_library_rt);
                let half_width = consensus_entry.median_peak_width / 2.0;

                if crate::trace::is_traced(mod_seq) {
                    log::info!(
                        "[trace] gap-fill: {} z={} file={} expected_rt={:.4} half_width={:.3} (missing from this file, passed in another)",
                        mod_seq,
                        charge,
                        file_name,
                        expected_rt,
                        half_width,
                    );
                }

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
    let n_filtered = filtered_by_mz.load(Ordering::Relaxed);
    if total_targets > 0 || n_filtered > 0 {
        log::info!(
            "Gap-fill: {} precursor-file pairs to score across {} files ({} passing precursors, {} files total)",
            total_targets,
            n_files_with_gaps,
            passing_precursors.len(),
            per_file_entries.len(),
        );
    }
    if n_filtered > 0 {
        log::info!(
            "Gap-fill: {} candidates skipped because precursor m/z is outside the file's isolation windows (GPF or disjoint m/z ranges)",
            n_filtered,
        );
    }

    result
}

/// Compute MAD from absolute residuals after sigma-clipping outliers.
///
/// Removes residuals above `clip_threshold`, then returns the median of the
/// survivors. If fewer than 20 points survive the clip, returns the raw
/// median of all residuals as a fallback (avoids unstable estimates from
/// tiny samples).
fn sigma_clipped_mad(abs_residuals: &[f64], clip_threshold: f64) -> f64 {
    if abs_residuals.is_empty() {
        return 0.0;
    }
    let mut clipped: Vec<f64> = abs_residuals
        .iter()
        .copied()
        .filter(|&r| r <= clip_threshold)
        .collect();
    if clipped.len() < 20 {
        // Too few survivors; fall back to raw median
        let mut all: Vec<f64> = abs_residuals.to_vec();
        all.sort_by(|a, b| a.total_cmp(b));
        return all[all.len() / 2];
    }
    clipped.sort_by(|a, b| a.total_cmp(b));
    clipped[clipped.len() / 2]
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

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_chromatography::RTCalibration;
    use osprey_core::CwtCandidate;

    // ---- sigma_clipped_mad tests ----

    #[test]
    fn test_sigma_clipped_mad_no_outliers() {
        // All residuals are small and tightly distributed
        let residuals = vec![
            0.05, 0.08, 0.10, 0.03, 0.07, 0.12, 0.06, 0.09, 0.04, 0.11, 0.05, 0.08, 0.10, 0.03,
            0.07, 0.12, 0.06, 0.09, 0.04, 0.11,
        ];
        let raw_mad = {
            let mut s = residuals.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            s[s.len() / 2]
        };
        let clip_threshold = raw_mad * 1.4826 * 3.0;
        let clipped = sigma_clipped_mad(&residuals, clip_threshold);
        // No outliers so clipped MAD should equal raw MAD
        assert!(
            (clipped - raw_mad).abs() < 1e-10,
            "clipped={}, raw={}",
            clipped,
            raw_mad
        );
    }

    #[test]
    fn test_sigma_clipped_mad_with_outliers() {
        // 80 inlier residuals at ~0.05 min, 20 outlier residuals at ~2.0 min
        let mut residuals: Vec<f64> = (0..80).map(|i| 0.03 + (i as f64) * 0.001).collect();
        residuals.extend(vec![2.0; 20]);
        let raw_mad = {
            let mut s = residuals.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            s[s.len() / 2]
        };
        // Raw MAD is pulled toward outliers (P50 of 100 values, ~0.06)
        let clip_threshold = raw_mad * 1.4826 * 3.0;
        let clipped = sigma_clipped_mad(&residuals, clip_threshold);
        // Clipped MAD should be close to the inlier distribution (~0.05)
        assert!(
            clipped < 0.08,
            "clipped MAD should be near inlier level (~0.05), got {}",
            clipped
        );
        // Outliers at 2.0 should not affect the clipped MAD
        assert!(
            clipped < raw_mad,
            "clipped MAD ({}) should be <= raw MAD ({})",
            clipped,
            raw_mad
        );
    }

    #[test]
    fn test_sigma_clipped_mad_too_few_survivors() {
        // >80% outliers: 5 inliers at 0.05, 95 outliers at 2.0
        let mut residuals: Vec<f64> = vec![0.05; 5];
        residuals.extend(vec![2.0; 95]);
        let raw_mad = {
            let mut s = residuals.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            s[s.len() / 2]
        };
        let clip_threshold = raw_mad * 1.4826 * 3.0;
        let clipped = sigma_clipped_mad(&residuals, clip_threshold);
        // Fewer than 20 inliers survive the clip so falls back to raw MAD
        assert!(
            (clipped - raw_mad).abs() < 1e-10,
            "should fall back to raw MAD ({}) when too few survivors, got {}",
            raw_mad,
            clipped
        );
    }

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
            parquet_index: 0,
            is_decoy,
            charge,
            scan_number: 1,
            apex_rt,
            start_rt,
            end_rt,
            coelution_sum,
            score: 0.0,
            run_precursor_qvalue: experiment_qvalue,
            run_peptide_qvalue: experiment_qvalue,
            run_protein_qvalue: 1.0,
            experiment_precursor_qvalue: experiment_qvalue,
            experiment_peptide_qvalue: experiment_qvalue,
            experiment_protein_qvalue: 1.0,
            pep: 1.0,
            modified_sequence: Arc::from(modified_sequence),
        }
    }

    // ---- Helper: detailed FdrEntry builder for consensus-gating / weighting tests ----
    #[allow(clippy::too_many_arguments)]
    fn make_consensus_entry(
        entry_id: u32,
        modified_sequence: &str,
        charge: u8,
        is_decoy: bool,
        apex_rt: f64,
        peak_width: f64,
        coelution_sum: f64,
        score: f64,
        run_precursor_qvalue: f64,
        run_peptide_qvalue: f64,
        run_protein_qvalue: f64,
    ) -> FdrEntry {
        FdrEntry {
            entry_id,
            parquet_index: 0,
            is_decoy,
            charge,
            scan_number: 1,
            apex_rt,
            start_rt: apex_rt - peak_width / 2.0,
            end_rt: apex_rt + peak_width / 2.0,
            coelution_sum,
            score,
            run_precursor_qvalue,
            run_peptide_qvalue,
            run_protein_qvalue,
            experiment_precursor_qvalue: run_precursor_qvalue,
            experiment_peptide_qvalue: run_peptide_qvalue,
            experiment_protein_qvalue: run_protein_qvalue,
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
    fn test_weighted_median_empty() {
        assert!((weighted_median(&[]) - 0.0).abs() < 1e-10);
    }

    // ---- compute_consensus_rts regression tests ----

    #[test]
    fn test_consensus_rejects_low_precursor_q_despite_protein_rescue() {
        // Regression: on Stellar, DAQVVGMTTTGAAK had three z=3 wrong-peak
        // detections with precursor_q = 0.099/0.091/0.016 that got rescued
        // into consensus via protein FDR because protein_q was low. Their
        // wrong-peak apex (8.46) out-weighted two correct z=2 detections
        // (8.67), pulling the consensus to the interferer. The tightened
        // qualification gate now rejects entries with precursor_q above
        // consensus_fdr regardless of protein rescue, so only the two good
        // z=2 detections anchor the consensus.
        let cal = make_identity_calibration();
        let cals: HashMap<String, RTCalibration> = ["file_20", "file_21", "file_22"]
            .iter()
            .map(|f| (f.to_string(), cal.clone()))
            .collect();

        let per_file_entries: Vec<(String, Vec<FdrEntry>)> = vec![
            (
                "file_20".to_string(),
                vec![
                    // z=3 wrong peak at 8.46, protein rescues (was poisoning)
                    make_consensus_entry(
                        1, "PEPTIDEK", 3, false, 8.46, 0.12, 5.48, -4.1, 0.099, 0.126, 0.005,
                    ),
                ],
            ),
            (
                "file_21".to_string(),
                vec![
                    // z=2 correct peak at 8.67 with strong precursor evidence
                    make_consensus_entry(
                        2, "PEPTIDEK", 2, false, 8.67, 0.24, 10.25, 0.32, 0.009, 0.010, 0.005,
                    ),
                    // z=3 wrong peak, rescued by protein, should now be rejected
                    make_consensus_entry(
                        3, "PEPTIDEK", 3, false, 8.44, 0.12, 7.11, -3.6, 0.091, 0.099, 0.005,
                    ),
                ],
            ),
            (
                "file_22".to_string(),
                vec![
                    // z=2 correct peak at 8.69
                    make_consensus_entry(
                        4, "PEPTIDEK", 2, false, 8.69, 0.18, 7.39, 1.61, 0.002, 0.002, 0.005,
                    ),
                    // z=3 wrong peak, rescued by protein, should now be rejected
                    make_consensus_entry(
                        5, "PEPTIDEK", 3, false, 8.46, 0.15, 6.03, -0.26, 0.016, 0.016, 0.005,
                    ),
                ],
            ),
        ];

        // With protein_fdr_threshold=0.01, the OLD gate would have rescued all
        // three z=3 entries. The NEW gate requires precursor_q <= consensus_fdr
        // (0.01) as a hard precondition, so they are rejected.
        let consensus = compute_consensus_rts(&per_file_entries, &cals, 0.01, 0.01);

        let peptide = consensus
            .iter()
            .find(|c| c.modified_sequence == "PEPTIDEK" && !c.is_decoy)
            .expect("PEPTIDEK should be in consensus (2 z=2 detections qualify)");

        assert_eq!(
            peptide.n_runs_detected, 2,
            "Only the two z=2 detections should anchor consensus; the three z=3 wrong-peak detections should be rejected by the precursor-q gate"
        );
        assert!(
            (peptide.consensus_library_rt - 8.68).abs() < 0.05,
            "Consensus should be on the correct peak at ~8.68, got {:.4}",
            peptide.consensus_library_rt
        );
    }

    #[test]
    fn test_consensus_weighting_downweights_negative_score_detections() {
        // Even if a low-quality wrong-peak detection somehow slips through the
        // qualification gate (e.g., its precursor-q is borderline but passing),
        // the sigmoid(score) weighting should crush its influence on the
        // weighted median. Here two detections have identical precursor-q but
        // one has score -4.0 (sigmoid ~ 0.018) and the other +1.5 (sigmoid ~
        // 0.818), a ~45x weight ratio.
        let cal = make_identity_calibration();
        let cals: HashMap<String, RTCalibration> = ["file_A", "file_B"]
            .iter()
            .map(|f| (f.to_string(), cal.clone()))
            .collect();

        let per_file_entries: Vec<(String, Vec<FdrEntry>)> = vec![
            (
                "file_A".to_string(),
                vec![make_consensus_entry(
                    1, "PEPTIDEK", 2, false, 8.40, 0.10, 3.0, -4.0, 0.005, 0.005, 1.0,
                )],
            ),
            (
                "file_B".to_string(),
                vec![make_consensus_entry(
                    2, "PEPTIDEK", 2, false, 8.70, 0.10, 3.0, 1.5, 0.005, 0.005, 1.0,
                )],
            ),
        ];

        let consensus = compute_consensus_rts(&per_file_entries, &cals, 0.01, 0.0);
        let peptide = consensus
            .iter()
            .find(|c| c.modified_sequence == "PEPTIDEK" && !c.is_decoy)
            .expect("PEPTIDEK should be in consensus");

        // The high-score detection's weighted contribution should dominate, so
        // the weighted median collapses onto 8.70 (not 8.40 and not the midpoint).
        assert!(
            (peptide.consensus_library_rt - 8.70).abs() < 1e-6,
            "Consensus should track the high-score detection at 8.70, got {:.4}",
            peptide.consensus_library_rt
        );
    }

    // ---- determine_reconcile_action tests ----
    //
    // These tests use the apex-proximity logic: Keep if |apex - expected| <= tolerance,
    // UseCwtPeak if a CWT candidate's apex is within tolerance, else ForcedIntegration.

    #[test]
    fn test_reconcile_action_keep_when_apex_near_expected() {
        // Apex at 10.0, expected at 10.0, tolerance 0.5 -> Keep
        let action = determine_reconcile_action(10.0, &[], 10.0, 0.5, 0.1);
        assert!(matches!(action, ReconcileAction::Keep));
    }

    #[test]
    fn test_reconcile_action_keep_at_tolerance_boundary() {
        // Apex exactly at tolerance distance -> Keep
        let action = determine_reconcile_action(10.5, &[], 10.0, 0.5, 0.1);
        assert!(matches!(action, ReconcileAction::Keep));
        // Just inside tolerance from the other side
        let action = determine_reconcile_action(9.5, &[], 10.0, 0.5, 0.1);
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
        // Current apex at 10.0, expected at 15.0, tolerance 0.5
        // Apex is 5.0 away (outside tolerance), CWT candidate 0 has apex 15.0 (within tolerance)
        let action = determine_reconcile_action(10.0, &cwt, 15.0, 0.5, 0.1);
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
        // Expected RT at 20.0, tolerance 0.5. Apex at 10.0 is far away.
        // CWT 0 (apex 5.0) is 15 min off. CWT 1 (apex 20.0) is within tolerance.
        let action = determine_reconcile_action(10.0, &cwt, 20.0, 0.5, 0.1);
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
        // CWT candidate apex at 5.0, expected at 20.0, tolerance 0.5
        // Neither current apex (10.0) nor CWT (5.0) is within tolerance of 20.0
        let cwt = vec![CwtCandidate {
            apex_rt: 5.0,
            start_rt: 4.0,
            end_rt: 6.0,
            area: 1000.0,
            snr: 5.0,
            coelution_score: 0.9,
        }];
        let action = determine_reconcile_action(10.0, &cwt, 20.0, 0.5, 1.0);
        match action {
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                assert!((expected_rt - 20.0).abs() < 1e-10);
                assert!((half_width - 1.0).abs() < 1e-10);
            }
            other => panic!("Expected ForcedIntegration, got {:?}", other),
        }
    }

    #[test]
    fn test_reconcile_action_forced_when_no_cwt() {
        // Apex at 10.0, expected at 20.0, tolerance 0.5, no CWT -> ForcedIntegration
        let action = determine_reconcile_action(10.0, &[], 20.0, 0.5, 1.5);
        match action {
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                assert!((expected_rt - 20.0).abs() < 1e-10);
                assert!((half_width - 1.5).abs() < 1e-10);
            }
            other => panic!("Expected ForcedIntegration, got {:?}", other),
        }
    }

    #[test]
    fn test_reconcile_action_wrong_apex_wide_boundaries() {
        // Regression test for the Seer dataset bug: peak with apex at 6.72 but
        // boundaries [6.49, 7.97] that happened to span the expected RT (7.98).
        // Old boundary-containment logic returned Keep; new apex-proximity logic
        // correctly returns ForcedIntegration since apex is 1.26 min off.
        let action = determine_reconcile_action(
            6.72,  // apex at wrong RT
            &[],   // no CWT candidates
            7.98,  // expected RT from consensus
            0.3,   // tolerance from refined calibration
            0.075, // half_width (median_peak_width / 2)
        );
        match action {
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                assert!((expected_rt - 7.98).abs() < 1e-10);
                assert!((half_width - 0.075).abs() < 1e-10);
            }
            other => panic!(
                "Expected ForcedIntegration for wrong-apex peak, got {:?}",
                other
            ),
        }
    }

    #[test]
    fn test_reconcile_action_cwt_picks_closest_apex() {
        // Two CWT candidates both within tolerance, picks the one with closest apex
        let cwt = vec![
            CwtCandidate {
                apex_rt: 10.3,
                start_rt: 10.0,
                end_rt: 10.6,
                area: 1000.0,
                snr: 5.0,
                coelution_score: 0.9,
            },
            CwtCandidate {
                apex_rt: 10.05,
                start_rt: 9.9,
                end_rt: 10.2,
                area: 500.0,
                snr: 3.0,
                coelution_score: 0.8,
            },
        ];
        // Current apex at 5.0 (far away), expected at 10.1, tolerance 0.5
        // CWT 0 apex 10.3 is 0.2 away (within tolerance)
        // CWT 1 apex 10.05 is 0.05 away (closer, also within tolerance)
        // Should pick CWT 1 (closest apex)
        let action = determine_reconcile_action(5.0, &cwt, 10.1, 0.5, 0.1);
        match action {
            ReconcileAction::UseCwtPeak {
                candidate_idx,
                apex_rt,
                ..
            } => {
                assert_eq!(candidate_idx, 1, "Should pick CWT with closest apex");
                assert!((apex_rt - 10.05).abs() < 1e-10);
            }
            other => panic!("Expected UseCwtPeak(1), got {:?}", other),
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

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01, 0.0);

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

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01, 0.0);
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

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01, 0.0);
        let target = consensus.iter().find(|c| !c.is_decoy).unwrap();
        // Weighted median should be near 15.0, not pulled to 30.0
        assert!(
            (target.consensus_library_rt - 15.0).abs() < 1.0,
            "Consensus RT {} should be near 15.0, not pulled by outlier",
            target.consensus_library_rt
        );
    }

    #[test]
    fn test_consensus_peak_width_excludes_non_passing_detections() {
        // 2 good detections (pass FDR, narrow 0.2 min peaks)
        // + 1 bad detection (fails FDR with run_qvalue=0.05, wide 2.0 min peak)
        // The bad detection should be excluded from consensus entirely,
        // so median_peak_width should be ~0.2, not influenced by the 2.0 min peak.
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
                    1, "PEPTIDEK", 2, false, 15.0, 14.9, 15.1, 100.0, 0.005,
                )], // width 0.2, passes FDR
            ),
            (
                "file2".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.1, 15.0, 15.2, 90.0, 0.008,
                )], // width 0.2, passes FDR
            ),
            (
                "file3".to_string(),
                vec![make_fdr_entry(
                    1, "PEPTIDEK", 2, false, 15.0, 14.0, 16.0, 2.0, 0.05,
                )], // width 2.0, FAILS FDR (run_qvalue=0.05 > consensus_fdr=0.01)
            ),
        ];

        let consensus = compute_consensus_rts(&per_file_entries, &per_file_calibrations, 0.01, 0.0);
        let target = consensus.iter().find(|c| !c.is_decoy).unwrap();
        // Bad detection excluded — only 0.2 min peaks contribute
        assert!(
            target.median_peak_width < 0.3,
            "Peak width {} should be ~0.2 min (bad detection excluded by FDR filter)",
            target.median_peak_width
        );
        // Only 2 runs should contribute to the consensus
        assert_eq!(target.n_runs_detected, 2, "Only 2 runs pass FDR, not 3");
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
            apex_library_rt_mad: None,
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
            apex_library_rt_mad: None,
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
            apex_library_rt_mad: None,
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
            apex_library_rt_mad: None,
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

    #[test]
    fn test_plan_reconciliation_reconciles_blib_eligible_precursors() {
        // Regression (Stellar TPTLQPTPEVHNGLR z=2): a precursor whose run-
        // level precursor and peptide q-values and experiment-level precursor
        // q-value all exceed 0.01, but whose experiment-level PEPTIDE q-value
        // is 0.0008 (well below 0.01), still ends up in the blib — every
        // observation gets admitted via peptide-level experiment FDR.
        // Reconciliation must cover these precursors in every file so their
        // per-file apexes align with the consensus; otherwise a wrong-peak
        // first-pass detection in one file stays at that wrong RT and
        // breaks the same-boundaries-across-charges invariant.
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 10.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: Some(0.02),
        }];

        // Entry mirrors sample 21 z=2 for TPTLQPTPEVHNGLR: poor run and
        // experiment precursor q-values, but strong experiment-peptide q
        // driven by other observations of the same peptide in other files.
        let mut entry = make_consensus_entry(
            1, "PEPTIDEK", 2, false, 11.0, 0.2, 5.0, -2.0, 0.12, 0.13, 1.0,
        );
        entry.experiment_precursor_qvalue = 0.04;
        entry.experiment_peptide_qvalue = 0.0008;

        let per_file_entries = vec![("file1".to_string(), vec![entry])];
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
            !actions.is_empty(),
            "Precursor passing only at experiment-peptide level must still \
             be reconciled; blib admits it via peptide-level FDR, so \
             reconciliation must keep its boundaries consistent across files."
        );
    }

    // ---- Global MAD tolerance tests ----
    //
    // These tests verify that `plan_reconciliation` uses a global MAD-derived
    // tolerance (robust to outliers) instead of per-point local interpolation.
    //
    // The bug these guard against: previously we used `cal.local_tolerance(query_rt)`
    // which interpolated abs_residuals at the query RT. A peptide with a wrong apex
    // contributed a large residual to the calibration training set, inflating the
    // local tolerance at that RT, which then allowed the wrong-RT detection to pass
    // the apex proximity check. Using a global MAD eliminates this feedback loop.

    /// Build an RTCalibration with specific abs_residuals for testing.
    /// The stats().mad = median of abs_residuals.
    fn make_calibration_with_residuals(abs_residuals: Vec<f64>) -> RTCalibration {
        let n = abs_residuals.len();
        let library_rts: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let fitted_rts: Vec<f64> = library_rts.clone();
        let params = osprey_chromatography::RTModelParams {
            library_rts,
            fitted_rts,
            abs_residuals,
        };
        RTCalibration::from_model_params(&params, 0.1).unwrap()
    }

    #[test]
    fn test_plan_reconciliation_tolerance_from_global_mad() {
        // Calibration with tight residuals except for one huge outlier at position 15.
        // Global MAD (median) should be near 0.05, giving tolerance ≈ 3 × 0.05 × 1.4826 = 0.22 min.
        let mut residuals = vec![0.05_f64; 30];
        residuals[15] = 1.50; // One huge outlier at library_rt = 15.0
        let cal = make_calibration_with_residuals(residuals);

        // Sanity check: MAD is robust to the single outlier
        let mad = cal.stats().mad;
        assert!(
            (mad - 0.05).abs() < 1e-6,
            "MAD should be 0.05 (median), not affected by single outlier, got {}",
            mad
        );

        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        // Consensus at library RT = 15.0 (exactly where the outlier residual sits)
        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: None,
        }];

        // Entry with apex 0.5 min off from consensus (should FAIL the tight tolerance).
        // With old local_tolerance: the outlier residual 1.5 at library RT 15 would
        // be interpolated -> tolerance ~4.5 min -> Keep (wrong).
        // With global MAD tolerance: ~0.22 min -> 0.5 > 0.22 -> re-score (correct).
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.5, 15.4, 15.6, 8.0, 0.005,
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

        // Should be re-scored (not in "Keep" state) because 0.5 min deviation > 0.22 min tolerance
        assert!(
            actions.contains_key(&("file1".to_string(), 0)),
            "Entry 0.5 min off consensus should be flagged for re-scoring with global MAD tolerance"
        );
    }

    #[test]
    fn test_plan_reconciliation_tolerance_minimum_floor() {
        // Calibration with essentially zero residuals (very tight fit).
        // Without a floor, tolerance would be ~0, failing legitimate small deviations.
        // The 0.1 min minimum floor prevents this.
        let cal = make_calibration_with_residuals(vec![0.001_f64; 30]);
        assert!(cal.stats().mad < 0.01, "MAD should be near zero");

        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: None,
        }];

        // Entry with apex 0.05 min off (within the 0.1 min minimum floor).
        // Should Keep because 0.05 <= floor(0.1).
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.05, 14.95, 15.15, 8.0, 0.005,
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

        // Entry within 0.1 min floor → Keep → not in actions map
        assert!(
            actions.is_empty(),
            "Entry within 0.1 min floor should be kept (actions should be empty)"
        );
    }

    #[test]
    fn test_plan_reconciliation_tolerance_matches_expected_mad_formula() {
        // Peptide MAD of 0.1 min -> tolerance = 3 × 0.1 × 1.4826 ≈ 0.445 min.
        // An entry 0.5 min off should fail; an entry 0.4 min off should pass.
        // The file calibration MAD is set wide enough that it does NOT clip
        // the peptide-derived tolerance (it acts only as a ceiling).
        let cal = make_calibration_with_residuals(vec![0.2_f64; 30]);
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let expected_tolerance = 3.0 * 0.1 * 1.4826; // ≈ 0.445
        assert!(expected_tolerance > 0.4 && expected_tolerance < 0.5);

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: Some(0.1),
        }];

        // Entry 1: 0.4 min off → within tolerance (0.445) → Keep
        // Entry 2: 0.5 min off → exceeds tolerance → re-score
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![
                make_fdr_entry(1, "PEPTIDEK", 2, false, 15.4, 15.3, 15.5, 8.0, 0.005),
                make_fdr_entry(2, "PEPTIDEK", 3, false, 15.5, 15.4, 15.6, 8.0, 0.005),
            ],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> =
            vec![("file1".to_string(), vec![vec![], vec![]])]
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

        // Entry 0 (0.4 min off): should Keep → not in actions
        assert!(
            !actions.contains_key(&("file1".to_string(), 0)),
            "Entry 0.4 min off should Keep (peptide-MAD tolerance 0.445)"
        );
        // Entry 1 (0.5 min off): should be re-scored
        assert!(
            actions.contains_key(&("file1".to_string(), 1)),
            "Entry 0.5 min off should be flagged (exceeds peptide-MAD tolerance 0.445)"
        );
    }

    #[test]
    fn test_plan_reconciliation_uses_global_within_peptide_mad_tighter_than_calibration() {
        // Single consensus peptide with within-peptide MAD = 0.02 min → global
        // median MAD = 0.02 → tolerance = 3 × 0.02 × 1.4826 = 0.089 floored to 0.1.
        // File cal MAD = 0.1 → ceiling = 0.445. Global-MAD tolerance (0.1 after floor)
        // wins because it's below the ceiling. An entry 0.15 min off — which the
        // old cross-peptide calibration formula (0.445) would Keep — is now re-scored.
        let cal = make_calibration_with_residuals(vec![0.1_f64; 30]);
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: Some(0.02),
        }];

        // Entry 0.15 min off from expected: |0.15| > 0.1 → re-score
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.15, 15.05, 15.25, 8.0, 0.005,
            )],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> = vec![(
            "file1".to_string(),
            vec![vec![CwtCandidate {
                apex_rt: 15.02,
                start_rt: 14.92,
                end_rt: 15.12,
                area: 1000.0,
                snr: 5.0,
                coelution_score: 0.9,
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

        assert!(
            actions.contains_key(&("file1".to_string(), 0)),
            "Entry 0.15 min off should be flagged when peptide MAD is 0.02 (tolerance ~0.1 after floor), even though file cal MAD would allow 0.445"
        );
    }

    #[test]
    fn test_plan_reconciliation_global_mad_is_median_across_peptides() {
        // Five peptides with MADs 0.01, 0.03, 0.05, 0.20, 0.50 → median = 0.05.
        // Tolerance = max(0.1, 3 × 0.05 × 1.4826 = 0.222) = 0.222 min.
        // An outlier peptide (MAD 0.50) does not drag the tolerance up.
        // An entry 0.3 min off for ANY peptide should be re-scored.
        let cal = make_calibration_with_residuals(vec![0.2_f64; 30]);
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus: Vec<PeptideConsensusRT> = ["PEP1", "PEP2", "PEP3", "PEP4", "PEP5"]
            .iter()
            .zip([0.01, 0.03, 0.05, 0.20, 0.50])
            .map(|(seq, mad)| PeptideConsensusRT {
                modified_sequence: seq.to_string(),
                is_decoy: false,
                consensus_library_rt: 15.0,
                median_peak_width: 0.2,
                n_runs_detected: 3,
                apex_library_rt_mad: Some(mad),
            })
            .collect();

        // For PEP3 (MAD 0.05, the median), an entry 0.3 min off should exceed
        // the global tolerance (0.222). An entry 0.15 min off should Keep.
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![
                make_fdr_entry(1, "PEP3", 2, false, 15.3, 15.2, 15.4, 8.0, 0.005),
                make_fdr_entry(2, "PEP3", 3, false, 15.15, 15.05, 15.25, 8.0, 0.005),
            ],
        )];

        let cwt: HashMap<String, Vec<Vec<CwtCandidate>>> =
            vec![("file1".to_string(), vec![vec![], vec![]])]
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
            actions.contains_key(&("file1".to_string(), 0)),
            "Entry 0.3 min off should be flagged (exceeds global-MAD tolerance 0.222)"
        );
        assert!(
            !actions.contains_key(&("file1".to_string(), 1)),
            "Entry 0.15 min off should Keep (within global-MAD tolerance 0.222)"
        );
    }

    #[test]
    fn test_plan_reconciliation_outliers_in_training_do_not_inflate_tolerance() {
        // Regression test for the "self-fulfilling tolerance" bug.
        //
        // Scenario: the refined calibration training set includes a peptide whose
        // detection was at the wrong RT (apex 1.2 min off the truth). That point
        // contributes a large abs_residual (~1.2) to the calibration. Under the
        // OLD local_tolerance interpolation, the query RT near this bad point
        // would pick up tolerance ~3.6 min, allowing the 1.2 min deviation to
        // pass the apex proximity check ("Keep").
        //
        // With the NEW global MAD, a single large residual among 30 small ones
        // barely moves the median, so tolerance stays tight (~0.22 min).

        // 30 tight residuals + 1 huge outlier at position 15
        let mut residuals = vec![0.05_f64; 30];
        residuals[15] = 1.20; // simulates a wrong-RT peptide in the training set

        let cal = make_calibration_with_residuals(residuals);
        let mad = cal.stats().mad;
        assert!(
            mad < 0.1,
            "MAD should be robust to the outlier, got {}",
            mad
        );

        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        // Consensus at library_rt = 15.0 (right at the outlier location)
        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: None,
        }];

        // This entry simulates the bad peptide itself: apex at 16.2 when consensus expects 15.0.
        // Boundaries are wide enough that old boundary-containment logic would keep it too.
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 16.2, 14.8, 16.3, 8.0, 0.005,
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

        // With global MAD tolerance, this 1.2 min deviation must be caught.
        assert!(
            actions.contains_key(&("file1".to_string(), 0)),
            "Wrong-RT peptide must be flagged for re-scoring — the outlier in the \
             calibration training set should not inflate the tolerance"
        );
    }

    #[test]
    fn test_plan_reconciliation_tolerance_capped_by_original_calibration() {
        // Regression test: the reconciliation tolerance must never exceed the
        // original (first-pass) calibration tolerance. If the refined calibration
        // is contaminated by wrong-peak detections, its MAD is inflated. The cap
        // at the original tolerance prevents the contamination from disabling
        // reconciliation entirely.
        //
        // Refined calibration: huge residuals → raw MAD ~1.0, sigma-clipped MAD
        //   still large. Without the cap, tolerance would be ~4.5 min.
        // Original calibration: tight residuals → MAD ~0.05, tolerance ~0.22 min.
        // The cap restricts the tolerance to 0.22 min.

        let refined_residuals = vec![1.0_f64; 30]; // all large (contaminated)
        let refined_cal = make_calibration_with_residuals(refined_residuals);
        let original_cal = make_calibration_with_residuals(vec![0.05_f64; 30]);

        let refined_cals: HashMap<String, RTCalibration> = vec![("file1".to_string(), refined_cal)]
            .into_iter()
            .collect();
        let original_cals: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), original_cal)]
                .into_iter()
                .collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.2,
            n_runs_detected: 3,
            apex_library_rt_mad: None,
        }];

        // Entry at 15.5 — 0.5 min from expected. With the uncapped refined
        // tolerance (~4.5 min) this would be Keep. With the original-calibration
        // cap (~0.22 min) it must be flagged for re-scoring.
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.5, 15.0, 16.0, 8.0, 0.005,
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
            &refined_cals,
            &original_cals,
            0.01,
        );

        assert!(
            actions.contains_key(&("file1".to_string(), 0)),
            "Entry 0.5 min from expected must be flagged when original calibration \
             tolerance is ~0.22 min, even though refined tolerance would allow it"
        );
    }

    #[test]
    fn test_plan_reconciliation_wrong_peak_1min_off_is_caught() {
        // Regression test: a peak 1.0 min from the consensus-predicted RT must
        // be flagged for re-scoring, not classified as Keep. This guards against
        // the tolerance inflation bug that allowed 36% of files on the Seer
        // dataset to have tolerances > 0.4 min (some up to 2.3 min).

        let cal = make_calibration_with_residuals(vec![0.05_f64; 30]);
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.5,
            n_runs_detected: 3,
            apex_library_rt_mad: None,
        }];

        // Peak at 16.0 — 1.0 min from expected (15.0). Clearly wrong.
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 16.0, 15.5, 16.5, 8.0, 0.005,
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
            actions.contains_key(&("file1".to_string(), 0)),
            "Peak 1.0 min from expected must be flagged for re-scoring — \
             this was the core reconciliation failure on multi-replicate data"
        );
    }

    #[test]
    fn test_plan_reconciliation_correct_peak_near_expected_is_kept() {
        // Counterpart to the wrong-peak test: a peak 0.05 min from expected
        // should be classified as Keep, not re-scored. Guards against the
        // tolerance being too tight and re-scoring correct peaks.

        let cal = make_calibration_with_residuals(vec![0.05_f64; 30]);
        let refined_cal: HashMap<String, RTCalibration> = vec![("file1".to_string(), cal.clone())]
            .into_iter()
            .collect();
        let original_cal: HashMap<String, RTCalibration> =
            vec![("file1".to_string(), cal)].into_iter().collect();

        let consensus = vec![PeptideConsensusRT {
            modified_sequence: "PEPTIDEK".to_string(),
            is_decoy: false,
            consensus_library_rt: 15.0,
            median_peak_width: 0.5,
            n_runs_detected: 3,
            apex_library_rt_mad: None,
        }];

        // Peak at 15.05 — 0.05 min from expected. Should be kept.
        let per_file_entries = vec![(
            "file1".to_string(),
            vec![make_fdr_entry(
                1, "PEPTIDEK", 2, false, 15.05, 14.8, 15.3, 8.0, 0.005,
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
            !actions.contains_key(&("file1".to_string(), 0)),
            "Peak 0.05 min from expected should be Keep (within tolerance)"
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
                apex_library_rt_mad: None,
            },
            PeptideConsensusRT {
                modified_sequence: "DECOY_PEPTIDEK".to_string(),
                is_decoy: true,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 2,
                apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
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
            apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
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
            apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
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
            apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
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
                apex_library_rt_mad: None,
            },
            PeptideConsensusRT {
                modified_sequence: "ANOTHERPEP".to_string(),
                is_decoy: false,
                consensus_library_rt: 25.0,
                median_peak_width: 1.5,
                n_runs_detected: 2,
                apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
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
            apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
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
            apex_library_rt_mad: None,
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
            &HashMap::new(),
            &HashMap::new(),
        );

        assert!(
            gap_fill.is_empty(),
            "Single file with the precursor → no gaps"
        );
    }

    #[test]
    fn test_gap_fill_filters_precursors_outside_isolation_windows() {
        // GPF-style scenario: 2 files cover disjoint m/z ranges.
        // file1 covers 400-500, file2 covers 500-600. A precursor detected in
        // file1 at m/z 450 should NOT generate a gap-fill target in file2,
        // because file2 cannot observe m/z 450.
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal = refined_cal.clone();

        let consensus = vec![
            PeptideConsensusRT {
                modified_sequence: "PEPTIDEK".to_string(),
                is_decoy: false,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 1,
                apex_library_rt_mad: None,
            },
            PeptideConsensusRT {
                modified_sequence: "DECOY_PEPTIDEK".to_string(),
                is_decoy: true,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 1,
                apex_library_rt_mad: None,
            },
        ];

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
                    2,
                    "OTHERPEPTIDE",
                    2,
                    false,
                    20.0,
                    19.5,
                    20.5,
                    6.0,
                    0.003,
                )],
            ),
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        // PEPTIDEK's precursor m/z = 450 (lives in file1's range)
        let lib_precursor_mz: HashMap<u32, f64> = vec![(1u32, 450.0f64)].into_iter().collect();

        // Disjoint GPF windows: file1 = [400,500), file2 = [500,600)
        let per_file_isolation_mz: HashMap<String, Vec<(f64, f64)>> = vec![
            ("file1".to_string(), vec![(400.0, 500.0)]),
            ("file2".to_string(), vec![(500.0, 600.0)]),
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
            &lib_precursor_mz,
            &per_file_isolation_mz,
        );

        // file2 should have NO gap-fill target for PEPTIDEK (m/z 450 not in [500,600))
        assert!(
            !gap_fill.contains_key("file2"),
            "file2 should not gap-fill a precursor whose m/z is outside its isolation windows"
        );
    }

    #[test]
    fn test_gap_fill_allows_precursors_inside_isolation_windows() {
        // Counterpart to the GPF test: precursor m/z IS inside file2's
        // isolation range, so the gap-fill target should be generated.
        let cal = make_identity_calibration();
        let refined_cal: HashMap<String, RTCalibration> = vec![
            ("file1".to_string(), cal.clone()),
            ("file2".to_string(), cal.clone()),
        ]
        .into_iter()
        .collect();
        let original_cal = refined_cal.clone();

        let consensus = vec![
            PeptideConsensusRT {
                modified_sequence: "PEPTIDEK".to_string(),
                is_decoy: false,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 1,
                apex_library_rt_mad: None,
            },
            PeptideConsensusRT {
                modified_sequence: "DECOY_PEPTIDEK".to_string(),
                is_decoy: true,
                consensus_library_rt: 15.0,
                median_peak_width: 1.0,
                n_runs_detected: 1,
                apex_library_rt_mad: None,
            },
        ];

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
                    2,
                    "OTHERPEPTIDE",
                    2,
                    false,
                    20.0,
                    19.5,
                    20.5,
                    6.0,
                    0.003,
                )],
            ),
        ];

        let lib_lookup: HashMap<(Arc<str>, u8), (u32, u32)> =
            vec![((Arc::from("PEPTIDEK"), 2u8), (1u32, 1 | 0x80000000))]
                .into_iter()
                .collect();

        // PEPTIDEK's precursor m/z = 450
        let lib_precursor_mz: HashMap<u32, f64> = vec![(1u32, 450.0f64)].into_iter().collect();

        // Overlapping replicate-style ranges: both files see 400-500
        let per_file_isolation_mz: HashMap<String, Vec<(f64, f64)>> = vec![
            ("file1".to_string(), vec![(400.0, 500.0)]),
            ("file2".to_string(), vec![(400.0, 500.0)]),
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
            &lib_precursor_mz,
            &per_file_isolation_mz,
        );

        // file2 SHOULD gap-fill PEPTIDEK because m/z 450 is inside its window
        let file2_targets = gap_fill
            .get("file2")
            .expect("file2 should gap-fill PEPTIDEK (m/z in range)");
        assert_eq!(file2_targets.len(), 1);
        assert_eq!(file2_targets[0].target_entry_id, 1);
    }
}
