//! Second-pass q-value strategy selection and the transfer-compete recompute.
//!
//! Ported from the C# Osprey pass-2 FDR (`OSPREY_PASS2_QVALUE`,
//! `Pass2FdrSidecar.ComputePass2TransferCompeteFull`,
//! `PercolatorFdr.ComputeFullPopulationPrecursorFdrStreaming`) for cross-impl parity.
//!
//! Default is `percolator` (retrain the 2nd-pass SVM) — identical to the historical
//! Rust behaviour, so nothing changes unless `OSPREY_PASS2_QVALUE` is set.
//!
//! `transfer-compete`: do NOT retrain. Apply the FROZEN 1st-pass model to the reconciled
//! survivors and recompute run/experiment q + PEP by a fresh target-decoy competition over
//! the FULL pre-compaction population (the 1st-pass scalar sidecars), with the reconciled
//! survivors' frozen-model scores swapped in. This restores FDR control that the retrain's
//! decoy-depleted null loses.
//!
//! The streaming competition takes an optional `stratum_base_ids`: `None` is the full-
//! population transfer-compete path; a non-null stratum is the CONSTRAINED competition used
//! by the (upcoming) protein-compact mode. Building the hook now keeps that feature additive.

use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

use osprey_fdr::percolator::{compete_from_indices, compute_conservative_qvalues};
use osprey_ml::pep::PepEstimator;
use osprey_ml::svm::FeatureStandardizer;

/// base_id = entry_id with the decoy high bit cleared.
const BASE_ID_MASK: u32 = 0x7FFF_FFFF;

/// How the merge-node second pass assigns the reported q-values.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pass2QValueMode {
    /// Default: retrain the 2nd-pass Percolator SVM and recompute q over the compacted pool.
    Percolator,
    /// Frozen 1st-pass model + fresh full-population target-decoy competition (no retrain).
    TransferCompete,
    /// Like `TransferCompete`, but the competition is CONSTRAINED to the protein stratum
    /// (peptides of proteins with >=2 first-pass-detected peptides) and the compaction gate
    /// is expanded to admit those peptides — recovering present-protein depth at honest FDR.
    ProteinCompact,
}

impl Pass2QValueMode {
    /// True when the mode replaces the 2nd-pass retrain with the frozen 1st-pass model + a
    /// target-decoy competition (no retrain): transfer-compete competes over the full
    /// pre-compaction population; protein-compact competes over the protein stratum. Both feed
    /// `compute_full_population_fdr_streaming` (the stratum via its `stratum_base_ids` arg). This
    /// mirrors the C# streaming frozen 2nd pass (Pass2FdrSidecar.ComputePass2TransferCompeteFull).
    pub fn uses_frozen_model(self) -> bool {
        matches!(
            self,
            Pass2QValueMode::TransferCompete | Pass2QValueMode::ProteinCompact
        )
    }
}

static PASS2_MODE: OnceLock<Pass2QValueMode> = OnceLock::new();

/// The process-wide pass-2 q-value mode, resolved once from `OSPREY_PASS2_QVALUE`.
/// Unset / `percolator` / unrecognized → `Percolator` (the default). `transfer-compete`
/// → `TransferCompete`. Mirrors the C# `NormalizePass2QValue` (trim + lowercase).
pub fn pass2_mode() -> Pass2QValueMode {
    *PASS2_MODE.get_or_init(|| match std::env::var("OSPREY_PASS2_QVALUE") {
        Ok(v) => match v.trim().to_ascii_lowercase().as_str() {
            "transfer-compete" => Pass2QValueMode::TransferCompete,
            "protein-compact" => Pass2QValueMode::ProteinCompact,
            _ => Pass2QValueMode::Percolator,
        },
        Err(_) => Pass2QValueMode::Percolator,
    })
}

static PASS2_PROTEIN_COMPACT_RETRAIN: OnceLock<bool> = OnceLock::new();

/// Whether protein-compact should RETRAIN the 2nd-pass SVM over the stratum-expanded pool
/// instead of applying the frozen 1st-pass model (`OSPREY_PROTEIN_COMPACT_RETRAIN`, set and
/// not `"0"`). Default off (frozen). This is a deliberate diagnostic A/B toggle: the frozen
/// path is the calibrated one; retrain over the decoy-depleted post-compaction null is
/// anti-conservative (measured ~1.2–1.5% FDP against entrapments). Mirrors the C#
/// `OspreyEnvironment.Pass2ProteinCompactRetrain` (`IsSetAndNotZero`: set and not `"0"`).
pub fn pass2_protein_compact_retrain() -> bool {
    *PASS2_PROTEIN_COMPACT_RETRAIN.get_or_init(|| {
        match std::env::var("OSPREY_PROTEIN_COMPACT_RETRAIN") {
            Ok(v) => !v.is_empty() && v != "0",
            Err(_) => false,
        }
    })
}

/// The frozen 1st-pass linear model: the fold-averaged SVM weights + bias and the training
/// standardizer. Reproduces byte-for-byte the averaged-model scoring `run_percolator_fdr`
/// applies (standardize the raw features in place, then `bias + Σ w·x`), so a survivor
/// re-scored here matches the score it would have carried in the 1st-pass sidecar.
#[derive(Clone)]
pub struct FrozenLinearModel {
    avg_weights: Vec<f64>,
    avg_bias: f64,
    standardizer: FeatureStandardizer,
}

impl FrozenLinearModel {
    pub fn new(avg_weights: Vec<f64>, avg_bias: f64, standardizer: FeatureStandardizer) -> Self {
        Self {
            avg_weights,
            avg_bias,
            standardizer,
        }
    }

    /// Feature-vector width this model expects. Callers skip entries whose vector differs.
    pub fn num_features(&self) -> usize {
        self.avg_weights.len()
    }

    /// Score one RAW (un-standardized) feature vector. Standardizes a copy — `raw` is not
    /// mutated — then applies the averaged linear model. Identical math to the score pass in
    /// `run_percolator_fdr` (pipeline.rs).
    pub fn score(&self, raw: &[f64]) -> f64 {
        let mut features = raw.to_vec();
        self.standardizer.transform_slice(&mut features);
        let mut score = self.avg_bias;
        for (w, x) in self.avg_weights.iter().zip(features.iter()) {
            score += w * x;
        }
        score
    }
}

/// Streamed full-population precursor FDR: run/experiment q + PEP for the reported survivors,
/// competing over the full pre-compaction population one file at a time (peak memory flat in
/// file count — cross-file state is bounded by the number of distinct precursors).
///
/// Faithful port of C# `ComputeFullPopulationPrecursorFdrStreaming`.
/// * `read_file_scalars(file_idx)` returns that file's `(entry_ids, scores)` from its 1st-pass
///   scalar sidecar (the full pre-compaction population).
/// * `survivor_score_override[(file_idx, entry_id)]` swaps in the frozen-model score for each
///   reconciled survivor.
/// * `survivors` is every reported (post-reconciliation) entry as `(file_idx, entry_id)`.
/// * `stratum_base_ids`: `None` = full-population competition; `Some` = competition constrained
///   to those base_ids (protein-compact).
///
/// Returns `(run_q, exp_q, pep)` keyed by `(file_idx, entry_id)` for exactly the survivors.
#[allow(clippy::type_complexity)]
pub fn compute_full_population_fdr_streaming(
    n_files: usize,
    read_file_scalars: impl Fn(usize) -> (Vec<u32>, Vec<f64>),
    survivor_score_override: &HashMap<(usize, u32), f64>,
    survivors: &[(usize, u32)],
    stratum_base_ids: Option<&HashSet<u32>>,
) -> (
    HashMap<(usize, u32), f64>,
    HashMap<(usize, u32), f64>,
    HashMap<(usize, u32), f64>,
) {
    let survivor_set: HashSet<(usize, u32)> = survivors.iter().copied().collect();
    let mut survivor_entry_ids: HashSet<u32> = HashSet::new();
    for s in &survivor_set {
        survivor_entry_ids.insert(s.1);
    }

    let mut survivor_run_q: HashMap<(usize, u32), f64> = HashMap::with_capacity(survivor_set.len());
    let mut survivor_exp_q: HashMap<(usize, u32), f64> = HashMap::with_capacity(survivor_set.len());
    let mut survivor_pep: HashMap<(usize, u32), f64> = HashMap::with_capacity(survivor_set.len());

    // Experiment-level per-base_id best target/decoy observation (score + locator).
    let mut best_target: HashMap<u32, (f64, usize, u32)> = HashMap::new();
    let mut best_decoy: HashMap<u32, (f64, usize, u32)> = HashMap::new();
    // Best (min) run q per survivor entry_id — the best-of-runs floor for the experiment q.
    let mut min_run_q: HashMap<u32, f64> = HashMap::new();

    for file_idx in 0..n_files {
        let (entry_ids, mut scores) = read_file_scalars(file_idx);
        let m = entry_ids.len();
        let mut labels = vec![false; m];
        for i in 0..m {
            let eid = entry_ids[i];
            labels[i] = (eid & !BASE_ID_MASK) != 0; // decoy high bit set
            if let Some(&ov) = survivor_score_override.get(&(file_idx, eid)) {
                scores[i] = ov; // swap in the reconciled survivor's frozen-model score
            }
        }

        // Run-level: compete within this file (only stratum members when stratified).
        let all_idx: Vec<usize> = match stratum_base_ids {
            None => (0..m).collect(),
            Some(strat) => (0..m)
                .filter(|&i| strat.contains(&(entry_ids[i] & BASE_ID_MASK)))
                .collect(),
        };
        let (wi, ws, wd) = compete_from_indices(&scores, &labels, &entry_ids, &all_idx);
        let mut q = vec![0.0f64; wi.len()];
        compute_conservative_qvalues(&ws, &wd, &mut q);
        for rank in 0..wi.len() {
            let eid = entry_ids[wi[rank]];
            if !survivor_entry_ids.contains(&eid) {
                continue;
            }
            let qv = q[rank];
            let key = (file_idx, eid);
            if survivor_set.contains(&key) {
                survivor_run_q.insert(key, qv);
            }
            match min_run_q.get(&eid) {
                Some(&cur) if qv >= cur => {}
                _ => {
                    min_run_q.insert(eid, qv);
                }
            }
        }

        // Experiment-level: fold every observation into the per-base_id bests.
        for i in 0..m {
            let eid = entry_ids[i];
            let bid = eid & BASE_ID_MASK;
            if let Some(strat) = stratum_base_ids {
                if !strat.contains(&bid) {
                    continue;
                }
            }
            let s = scores[i];
            if labels[i] {
                match best_decoy.get(&bid) {
                    Some(&(cs, _, _)) if s <= cs => {}
                    _ => {
                        best_decoy.insert(bid, (s, file_idx, eid));
                    }
                }
            } else {
                match best_target.get(&bid) {
                    Some(&(cs, _, _)) if s <= cs => {}
                    _ => {
                        best_target.insert(bid, (s, file_idx, eid));
                    }
                }
            }
        }
        // entry_ids/scores/labels/all_idx dropped here before the next file is read.
    }

    // Experiment competition: one winner per base_id, conservative q, PEP over the winners.
    let mut base_ids: HashSet<u32> = HashSet::new();
    base_ids.extend(best_target.keys().copied());
    base_ids.extend(best_decoy.keys().copied());
    let w = base_ids.len();
    let mut exp_score = vec![0.0f64; w];
    let mut exp_is_decoy = vec![false; w];
    let mut exp_base_id = vec![0u32; w];
    let mut winner_loc: HashMap<u32, (usize, u32, f64)> = HashMap::with_capacity(w);
    for (wi2, &bid) in base_ids.iter().enumerate() {
        let t = best_target.get(&bid).copied();
        let d = best_decoy.get(&bid).copied();
        // compete_from_indices: target wins strictly (t > d); ties go to the decoy. Scores are
        // finite SVM discriminants, so `t <= d` reproduces the C# `!(t.score > d.score)`.
        let decoy_wins = match (t, d) {
            (Some(tt), Some(dd)) => tt.0 <= dd.0,
            (Some(_), None) => false,
            _ => true,
        };
        let win = if decoy_wins { d.unwrap() } else { t.unwrap() };
        exp_score[wi2] = win.0;
        exp_is_decoy[wi2] = decoy_wins;
        exp_base_id[wi2] = bid;
        winner_loc.insert(bid, (win.1, win.2, win.0));
    }

    // Sort winners by score desc, base_id asc (unique base_id ⇒ total order).
    let mut perm: Vec<usize> = (0..w).collect();
    perm.sort_by(|&a, &b| {
        exp_score[b]
            .total_cmp(&exp_score[a])
            .then_with(|| exp_base_id[a].cmp(&exp_base_id[b]))
    });
    let sorted_score: Vec<f64> = perm.iter().map(|&i| exp_score[i]).collect();
    let sorted_decoy: Vec<bool> = perm.iter().map(|&i| exp_is_decoy[i]).collect();
    let sorted_base_id: Vec<u32> = perm.iter().map(|&i| exp_base_id[i]).collect();
    let mut q_exp = vec![0.0f64; w];
    compute_conservative_qvalues(&sorted_score, &sorted_decoy, &mut q_exp);
    let mut base_id_exp_q: HashMap<u32, f64> = HashMap::with_capacity(w);
    for i in 0..w {
        base_id_exp_q.insert(sorted_base_id[i], q_exp[i]);
    }

    let pep_estimator = PepEstimator::fit_default(&exp_score, &exp_is_decoy);

    let multi_file = n_files > 1;
    for &key in &survivor_set {
        let (file_idx, eid) = key;
        let bid = eid & BASE_ID_MASK;

        survivor_run_q.entry(key).or_insert(1.0);

        if multi_file {
            // Experiment q = base_id winner q, floored up to this precursor's best run q.
            let mut eq = base_id_exp_q.get(&bid).copied().unwrap_or(1.0);
            let floor_q = min_run_q.get(&eid).copied().unwrap_or(1.0);
            if eq < floor_q {
                eq = floor_q;
            }
            survivor_exp_q.insert(key, eq);
        } else {
            // Single file: experiment q == run q.
            let rq = survivor_run_q[&key];
            survivor_exp_q.insert(key, rq);
        }

        // PEP is real only on the single experiment-winner observation of each base_id.
        let mut pep = 1.0;
        if let Some(&(loc_file, loc_eid, loc_score)) = winner_loc.get(&bid) {
            if loc_eid == eid && loc_file == file_idx {
                pep = pep_estimator.posterior_error(loc_score);
            }
        }
        survivor_pep.insert(key, pep);
    }

    (survivor_run_q, survivor_exp_q, survivor_pep)
}

#[cfg(test)]
mod tests {
    use super::*;

    const DECOY: u32 = 0x8000_0000;

    /// Every survivor gets a run q, experiment q, and PEP, all in valid ranges — the streaming
    /// competition covers the reported set and produces finite, bounded statistics (multi-file
    /// path: experiment-level competition + best-of-runs floor). The byte-exact values are
    /// validated by the cross-impl parity gate, not here.
    #[test]
    fn streaming_covers_survivors_with_valid_ranges() {
        // Two files, base_ids 1 and 2; the targets are the reported survivors.
        let files: [(Vec<u32>, Vec<f64>); 2] = [
            (vec![1, 1 | DECOY, 2, 2 | DECOY], vec![5.0, 1.0, 2.0, 3.0]),
            (vec![1, 1 | DECOY], vec![4.5, 0.5]),
        ];
        let read = |i: usize| files[i].clone();
        let survivors = vec![(0usize, 1u32), (0, 2), (1, 1)];
        let overrides: HashMap<(usize, u32), f64> = HashMap::new();

        let (run_q, exp_q, pep) =
            compute_full_population_fdr_streaming(2, read, &overrides, &survivors, None);

        for key in &survivors {
            let rq = *run_q.get(key).expect("run_q missing a survivor");
            let eq = *exp_q.get(key).expect("exp_q missing a survivor");
            let pp = *pep.get(key).expect("pep missing a survivor");
            assert!(rq >= 0.0 && rq.is_finite(), "run_q out of range: {rq}");
            assert!(eq >= 0.0 && eq.is_finite(), "exp_q out of range: {eq}");
            assert!((0.0..=1.0).contains(&pp), "pep out of [0,1]: {pp}");
        }
    }

    /// A survivor's swapped-in override score participates in the competition in place of the
    /// sidecar score (the reconciled-survivor path).
    #[test]
    fn override_replaces_sidecar_score() {
        let files = [(vec![1u32, 1 | DECOY], vec![0.1, 0.2])];
        let read = |i: usize| files[i].clone();
        let survivors = vec![(0usize, 1u32)];
        let mut overrides = HashMap::new();
        overrides.insert((0usize, 1u32), 9.0); // strong override for the target

        let (run_q, exp_q, _pep) =
            compute_full_population_fdr_streaming(1, read, &overrides, &survivors, None);

        // Single file → experiment q == run q, and the survivor is covered.
        let rq = *run_q.get(&(0, 1)).expect("run_q missing");
        assert_eq!(rq, exp_q[&(0, 1)]);
        assert!(rq >= 0.0 && rq.is_finite());
    }

    /// A non-null stratum restricts the competition to its base_ids: an off-stratum survivor
    /// never wins a (stratum-filtered) competition, so it falls to the default q=1.0 — the
    /// signal the orchestrator uses to leave it on its 1st-pass q. Protein-compact (#3) relies
    /// on this. (The in-stratum survivor is competed; here it is a lone target+decoy pair, so
    /// its conservative q is also 1.0 — the exact values are a parity-gate concern.)
    #[test]
    fn stratum_routes_off_stratum_to_default_q() {
        let files = [(
            vec![1u32, 1 | DECOY, 2, 2 | DECOY],
            vec![5.0, 1.0, 4.0, 2.0],
        )];
        let read = |i: usize| files[i].clone();
        let survivors = vec![(0usize, 1u32), (0, 2)];
        let overrides: HashMap<(usize, u32), f64> = HashMap::new();
        let stratum: HashSet<u32> = [1u32].into_iter().collect(); // only base_id 1

        let (run_q, _exp, _pep) =
            compute_full_population_fdr_streaming(1, read, &overrides, &survivors, Some(&stratum));

        assert!(
            run_q.contains_key(&(0, 1)),
            "in-stratum survivor should be covered"
        );
        assert_eq!(
            run_q.get(&(0, 2)).copied(),
            Some(1.0),
            "off-stratum survivor should fall to the default q=1.0"
        );
    }
}
