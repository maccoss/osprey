//! Second-pass q-value strategy selection and the transfer-compete recompute.
//!
//! Ported from the C# Osprey pass-2 FDR (`OSPREY_PASS2_QVALUE`,
//! `Pass2FdrSidecar.ComputePass2TransferCompeteFull`,
//! `PercolatorFdr.ComputeFullPopulationPrecursorFdrStreaming`) for cross-impl parity.
//!
//! Default is `protein-compact`. The former `percolator` default — retrain the 2nd-pass SVM
//! and recompute q over the COMPACTED pool — has been removed: compaction depletes the decoys,
//! so the retrained null is thin and the reported q is anti-conservative (measured ~9% true FDP
//! at a nominal 1% on an 82-file SEA-AD entrapment set). The linear model trained by the 1st-pass
//! SVM is now the model for pass 2 in every mode; only the diagnostic
//! `OSPREY_PROTEIN_COMPACT_RETRAIN` A/B still retrains.
//!
//! `transfer-compete`: do NOT retrain. Apply the FROZEN 1st-pass model to the reconciled
//! survivors and recompute run/experiment q + PEP by a fresh target-decoy competition over
//! the FULL pre-compaction population (the 1st-pass scalar sidecars), with the reconciled
//! survivors' frozen-model scores swapped in. This restores FDR control that the retrain's
//! decoy-depleted null loses.
//!
//! The streaming competition takes an optional `stratum_base_ids`: `None` is the full-
//! population transfer-compete path; a non-null stratum is the STRATIFIED competition used by
//! protein-compact. The stratum and its complement compete separately, each against its own
//! paired decoys, at run scope and at experiment scope (C# issue #4581): the stratum is chosen
//! by target evidence, so it takes most of the high-scoring targets and only the base rate of
//! the high-scoring decoys, and the decoys it leaves behind must still count against the
//! complement's targets. Port of C# `StreamingFdr.ComputeComplementExperimentQ` and the split
//! `CompeteOneFile`.

use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;

use osprey_fdr::percolator::{compete_from_indices, compute_conservative_qvalues};
use osprey_ml::pep::PepEstimator;
use osprey_ml::svm::FeatureStandardizer;

/// base_id = entry_id with the decoy high bit cleared.
const BASE_ID_MASK: u32 = 0x7FFF_FFFF;

/// How the merge-node second pass assigns the reported q-values.
/// Every mode applies the FROZEN 1st-pass linear model; they differ only in the population the
/// target-decoy competition runs over. Nothing here retrains.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Pass2QValueMode {
    /// Frozen 1st-pass model + fresh full-population target-decoy competition (no retrain).
    TransferCompete,
    /// Default. Like `TransferCompete`, but the competition is CONSTRAINED to the protein
    /// stratum (peptides of proteins with >=2 first-pass-detected peptides) and the compaction
    /// gate is expanded to admit those peptides — recovering present-protein depth.
    ProteinCompact,
}

static PASS2_MODE: OnceLock<Pass2QValueMode> = OnceLock::new();

/// The process-wide pass-2 q-value mode, resolved once from `OSPREY_PASS2_QVALUE`.
/// Unset → `ProteinCompact` (the default). Mirrors the C# `NormalizePass2QValue`
/// (trim + lowercase).
///
/// An unrecognized value PANICS rather than falling back. Silently substituting the default
/// would report numbers the caller did not ask for under an output name that says otherwise —
/// and the removed `percolator` token in particular is one that existing sweep scripts still
/// pass. These environment variables are developer-side instrumentation, not a supported
/// interface, so a removed name is expected to fail loudly rather than be aliased forever.
pub fn pass2_mode() -> Pass2QValueMode {
    *PASS2_MODE.get_or_init(|| match std::env::var("OSPREY_PASS2_QVALUE") {
        Ok(v) if !v.trim().is_empty() => match v.trim().to_ascii_lowercase().as_str() {
            "transfer-compete" => Pass2QValueMode::TransferCompete,
            "protein-compact" => Pass2QValueMode::ProteinCompact,
            other => panic!(
                "OSPREY_PASS2_QVALUE='{other}' is not a recognized mode. Recognized: \
                 'transfer-compete', 'protein-compact'. Unset it for the default \
                 (protein-compact). The 'percolator' mode was REMOVED: it retrained the \
                 2nd-pass SVM on a compaction-depleted decoy pool, which reports \
                 anti-conservative q-values."
            ),
        },
        _ => Pass2QValueMode::ProteinCompact,
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
/// * `stratum_base_ids`: `None` = full-population competition; `Some` = STRATIFIED competition
///   (protein-compact): these base_ids and their complement compete separately, each against
///   its own paired decoys. Within each file the run competition splits in two by base_id; the
///   experiment fold below is the stratum only.
/// * `complement_experiment_q`: protein-compact only — the experiment q of the stratum's
///   complement, from [`compute_complement_experiment_q`], merged into the winner-q map so a
///   complement survivor's experiment q is its own competition's (still floored to its best run
///   q). Its base_ids are disjoint from the stratum's, so the two never share a key. `None`
///   leaves complement entries at experiment q 1.0.
///
/// Returns `(run_q, exp_q, pep, exp_aggregate)`. The first three are keyed by
/// `(file_idx, entry_id)` for exactly the survivors; `exp_aggregate` is keyed by full
/// `entry_id` (one value per entry, not per observation) and is the score this competition
/// ranked that entry on — the caller writes it beside the experiment q it produced.
///
/// An `entry_id` ABSENT from `exp_aggregate` took no part in the experiment fold, which under
/// protein-compact means an OFF-STRATUM entry; the caller carries its pass-1 aggregate, because
/// the complement's experiment competition ranked it on that aggregate. Absence (an `Option`
/// from the lookup), never an in-band sentinel: the score is a signed discriminant, so 0.0 is
/// an ordinary mid-distribution value and a consumer building a score-space acceptance
/// boundary would take a minimum over fabricated zeros. NaN is rejected for the same reason it
/// is on the C# side — it propagates silently, and the sidecar comparators test
/// `(a - b).abs() <= tolerance`, which is FALSE for NaN against NaN, turning byte-identical
/// files into a red gate.
#[allow(clippy::type_complexity)]
pub fn compute_full_population_fdr_streaming(
    n_files: usize,
    read_file_scalars: impl Fn(usize) -> (Vec<u32>, Vec<f64>),
    survivor_score_override: &HashMap<(usize, u32), f64>,
    survivors: &[(usize, u32)],
    stratum_base_ids: Option<&HashSet<u32>>,
    complement_experiment_q: Option<&HashMap<u32, f64>>,
) -> (
    HashMap<(usize, u32), f64>,
    HashMap<(usize, u32), f64>,
    HashMap<(usize, u32), f64>,
    HashMap<u32, f64>,
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

        // Run-level: compete within this file, conservative q on the winners. Stratified
        // (protein-compact), the file splits in TWO by base_id - the stratum and its complement -
        // and each half competes against its OWN paired decoys. Both halves are reported, so
        // both need a calibrated q.
        //
        // The complement is not optional. The stratum is chosen by TARGET evidence (proteins with
        // 2 or more detected target peptides), so it takes most of the high-scoring targets and
        // only the base rate of the high-scoring decoys; the decoys it leaves behind are
        // concentrated in the complement. When the complement's targets kept the q they earned
        // against the whole population, those decoys counted against no one. Competing the
        // complement against its own decoys gives every decoy the same protein evidence as its
        // paired target - each counts in the half its target is in.
        //
        // A peak Stage 6 changed (reconciliation moved it, or gap-fill created it) competes on its
        // new score in whichever half its base_id belongs to. Membership is protein evidence
        // alone, decided before Stage 6, and base_id keeps a target and its paired decoy in the
        // same half always.
        let mut compete_half = |indices: &[usize]| {
            let (wi, ws, wd) = compete_from_indices(&scores, &labels, &entry_ids, indices);
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
        };
        match stratum_base_ids {
            None => {
                let all_idx: Vec<usize> = (0..m).collect();
                compete_half(&all_idx);
            }
            Some(strat) => {
                let (in_idx, off_idx): (Vec<usize>, Vec<usize>) =
                    (0..m).partition(|&i| strat.contains(&(entry_ids[i] & BASE_ID_MASK)));
                compete_half(&in_idx);
                compete_half(&off_idx);
            }
        }

        // Experiment-level: fold every observation into the per-base_id bests. When stratified
        // this is the STRATUM ONLY. The complement's experiment q is competed from the
        // analysis-wide 1st-pass experiment records instead (`compute_complement_experiment_q`),
        // which hold every base_id the complement needs - its non-survivors included - without a
        // per-file best for each of them.
        //
        // A pass-1 aggregate is the right score for the complement because reconciliation
        // anchors on the best-scoring peak and corrects the others toward it, so a changed peak
        // is never the one that supplied the maximum.
        for i in 0..m {
            let eid = entry_ids[i];
            let bid = eid & BASE_ID_MASK;
            let in_stratum = match stratum_base_ids {
                None => true,
                Some(strat) => strat.contains(&bid),
            };
            if !in_stratum {
                continue;
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
    // The experiment aggregate score PER FULL ENTRY_ID (sidecar v4) — the score this
    // competition ranked each entry on. Keyed by entry_id and not by base_id, exactly as the
    // 1st-pass producer (`compute_fdr_from_stubs`) keys it: a target and its decoy are
    // distinct entries with distinct aggregates even though the competition below pairs them.
    // `winner_loc` cannot serve this purpose — it holds only the WINNER of each pair, so
    // reading it for a decoy would hand back the target's score whenever the target won.
    //
    // `best_target` / `best_decoy` already carry exactly what is needed: each is the max over
    // that entry's observations across every file, the same max-over-rows reduction the
    // 1st-pass producer performs. Rust has no mean-best-N aggregation mode, so this max IS the
    // aggregate — no branch needed. Mirrors C# `StreamingFdr`'s `aggByEntryId`.
    let mut agg_by_entry_id: HashMap<u32, f64> = HashMap::with_capacity(w * 2);
    for (wi2, &bid) in base_ids.iter().enumerate() {
        let t = best_target.get(&bid).copied();
        let d = best_decoy.get(&bid).copied();
        if let Some(tt) = t {
            agg_by_entry_id.insert(tt.2, tt.0);
        }
        if let Some(dd) = d {
            agg_by_entry_id.insert(dd.2, dd.0);
        }
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

    let mut exp_q_by_winner_id = winner_q_by_entry_id(&exp_score, &exp_is_decoy, &exp_base_id);
    // protein-compact: the stratum's complement competed separately. Its base_ids are disjoint
    // from the stratum's, so the two maps never share a key.
    if let Some(complement) = complement_experiment_q {
        for (&eid, &q) in complement {
            exp_q_by_winner_id.insert(eid, q);
        }
    }

    let pep_estimator = PepEstimator::fit_default(&exp_score, &exp_is_decoy);

    let multi_file = n_files > 1;
    for &key in &survivor_set {
        let (file_idx, eid) = key;
        let bid = eid & BASE_ID_MASK;

        survivor_run_q.entry(key).or_insert(1.0);

        if multi_file {
            // Experiment q = the winner's q, floored up to this precursor's best run q. Full
            // entry_id: an entry takes the experiment q only when it was the side that WON its
            // own target/decoy competition. The loser keeps 1.0.
            let mut eq = exp_q_by_winner_id.get(&eid).copied().unwrap_or(1.0);
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

    (
        survivor_run_q,
        survivor_exp_q,
        survivor_pep,
        agg_by_entry_id,
    )
}

/// The experiment q of the stratum's COMPLEMENT under protein-compact: one target/decoy
/// competition per base_id over every base_id NOT in `stratum_base_ids`, on the per-entry
/// experiment aggregate score the 1st pass ranked it on, with the same tie rule, sort and
/// conservative q as the stratum's own experiment competition. Keyed by the WINNER's full
/// entry_id; a loser is absent and reads as 1.0.
///
/// The stratum holds the base_ids of proteins with 2 or more detected TARGET peptides, so it
/// takes most of the high-scoring targets and only the base rate of the high-scoring decoys.
/// The complement is therefore decoy-rich, and the pass-1 experiment q its targets used to keep
/// was earned against a null diluted by the stratum's true targets. Competing the complement
/// against its own decoys makes the two halves each calibrated, and so the reported union.
///
/// `aggregate_scores` is every analysis-wide 1st-pass experiment record as
/// `(entry_id, aggregate score)`, targets and decoys alike. Order does not matter: the
/// competition is a per-base_id reduction followed by a total sort. Port of C#
/// `StreamingFdr.ComputeComplementExperimentQ`.
pub fn compute_complement_experiment_q(
    aggregate_scores: impl IntoIterator<Item = (u32, f64)>,
    stratum_base_ids: &HashSet<u32>,
) -> HashMap<u32, f64> {
    let mut best_target: HashMap<u32, f64> = HashMap::new();
    let mut best_decoy: HashMap<u32, f64> = HashMap::new();
    for (eid, score) in aggregate_scores {
        let bid = eid & BASE_ID_MASK;
        if stratum_base_ids.contains(&bid) {
            continue;
        }
        let best = if (eid & !BASE_ID_MASK) != 0 {
            &mut best_decoy
        } else {
            &mut best_target
        };
        match best.get(&bid) {
            Some(&cur) if score <= cur => {}
            _ => {
                best.insert(bid, score);
            }
        }
    }
    let mut base_ids: HashSet<u32> = HashSet::new();
    base_ids.extend(best_target.keys().copied());
    base_ids.extend(best_decoy.keys().copied());
    let w = base_ids.len();
    let mut exp_score = Vec::with_capacity(w);
    let mut exp_is_decoy = Vec::with_capacity(w);
    let mut exp_base_id = Vec::with_capacity(w);
    for bid in base_ids {
        let t = best_target.get(&bid).copied();
        let d = best_decoy.get(&bid).copied();
        // The stratum's rule: target wins strictly (t > d); ties go to the decoy.
        let (score, decoy_wins) = match (t, d) {
            (Some(tt), Some(dd)) => {
                if tt > dd {
                    (tt, false)
                } else {
                    (dd, true)
                }
            }
            (Some(tt), None) => (tt, false),
            (None, Some(dd)) => (dd, true),
            (None, None) => unreachable!("base_id came from one of the two maps"),
        };
        exp_score.push(score);
        exp_is_decoy.push(decoy_wins);
        exp_base_id.push(bid);
    }
    winner_q_by_entry_id(&exp_score, &exp_is_decoy, &exp_base_id)
}

/// Conservative experiment q for one winner per base_id, keyed by the WINNER's full entry_id.
/// Shared by the stratum's competition and its complement's so the two halves rank, break ties
/// and estimate q identically.
fn winner_q_by_entry_id(
    exp_score: &[f64],
    exp_is_decoy: &[bool],
    exp_base_id: &[u32],
) -> HashMap<u32, f64> {
    let w = exp_score.len();
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
    // Keyed by the WINNER's full entry_id, decoy bit intact — never the shared base_id.
    // `sorted_decoy` carries which side won, so the winner's entry_id is reconstructible here
    // without a second array, leaving the competition's sort and tie-break untouched. See
    // `compute_experiment_precursor_qvalues` (osprey-fdr) for the defect this closes: on
    // base_id, a target whose DECOY won inherited the winner's q.
    let mut exp_q_by_winner_id: HashMap<u32, f64> = HashMap::with_capacity(w);
    for i in 0..w {
        let winner_id = if sorted_decoy[i] {
            sorted_base_id[i] | !BASE_ID_MASK
        } else {
            sorted_base_id[i]
        };
        exp_q_by_winner_id.insert(winner_id, q_exp[i]);
    }

    exp_q_by_winner_id
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

        let (run_q, exp_q, pep, agg) =
            compute_full_population_fdr_streaming(2, read, &overrides, &survivors, None, None);

        // Every entry that took part in the (unstratified) experiment fold has an aggregate:
        // the max of its score over the files it was observed in. base_id 1's target is seen
        // at 5.0 and 4.5, so its aggregate is the larger.
        assert_eq!(agg.get(&1u32).copied(), Some(5.0));
        assert_eq!(agg.get(&(1 | DECOY)).copied(), Some(1.0));
        assert_eq!(agg.get(&2u32).copied(), Some(2.0));
        assert_eq!(agg.get(&(2 | DECOY)).copied(), Some(3.0));

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

        let (run_q, exp_q, _pep, agg) =
            compute_full_population_fdr_streaming(1, read, &overrides, &survivors, None, None);

        // Single file → experiment q == run q, and the survivor is covered.
        let rq = *run_q.get(&(0, 1)).expect("run_q missing");
        assert_eq!(rq, exp_q[&(0, 1)]);
        assert!(rq >= 0.0 && rq.is_finite());
        // The aggregate is the score the competition ranked the entry on — the override, not
        // the sidecar score it replaced.
        assert_eq!(agg.get(&1u32).copied(), Some(9.0));
    }

    /// Stratified (protein-compact), the run competition splits by base_id into the stratum and
    /// its complement, and each half competes against its own paired decoys. The stratum holds
    /// four strong targets; the complement holds three weaker targets and the one decoy that
    /// outscores them. In one population that decoy is diluted by the stratum's targets; in its
    /// own half it counts against the complement's targets alone.
    #[test]
    fn stratum_and_complement_compete_separately_at_run_scope() {
        let mut ids = Vec::new();
        let mut scores = Vec::new();
        for (bid, t) in [
            (1u32, 10.0),
            (2, 9.5),
            (3, 9.0),
            (4, 8.5),
            (5, 3.0),
            (6, 2.9),
            (7, 2.8),
        ] {
            ids.extend([bid, bid | DECOY]);
            scores.extend([t, 0.1]);
        }
        ids.push(9 | DECOY); // unpaired complement decoy that outscores the complement targets
        scores.push(5.0);
        let files = [(ids, scores)];
        let survivors = vec![(0usize, 1u32), (0, 5), (0, 7)];
        let overrides: HashMap<(usize, u32), f64> = HashMap::new();
        let stratum: HashSet<u32> = [1u32, 2, 3, 4].into_iter().collect();

        let read = |i: usize| files[i].clone();
        let (run_q, exp_q, _pep, agg) = compute_full_population_fdr_streaming(
            1,
            read,
            &overrides,
            &survivors,
            Some(&stratum),
            None,
        );
        // Stratum half: four target winners -> (0+1)/4.
        assert_eq!(run_q[&(0, 1)], 0.25);
        // Complement half: D(5.0) T(3.0) T(2.9) T(2.8) -> forward 1, 2/1, 2/2, 2/3; the
        // backward sweep gives every complement target 2/3.
        assert_eq!(run_q[&(0, 5)], 2.0 / 3.0);
        assert_eq!(run_q[&(0, 7)], 2.0 / 3.0);
        // Single file: experiment q == run q in both halves.
        assert_eq!(exp_q[&(0, 5)], 2.0 / 3.0);
        // The experiment fold is the stratum only, so an off-stratum entry has no aggregate
        // here; the caller carries its pass-1 aggregate.
        assert_eq!(agg.get(&1u32).copied(), Some(10.0));
        assert_eq!(agg.get(&5u32).copied(), None);

        // Unstratified, the same file is one competition and the decoy sits behind seven
        // targets' worth of null: the complement targets get 2/7.
        let read = |i: usize| files[i].clone();
        let (run_q_all, _, _, _) =
            compute_full_population_fdr_streaming(1, read, &overrides, &survivors, None, None);
        assert_eq!(run_q_all[&(0, 5)], 2.0 / 7.0);
        assert_eq!(run_q_all[&(0, 1)], 0.25);
    }

    /// Multi-file, a complement survivor takes the complement's experiment q (merged into the
    /// winner map), floored to its best run q; without the complement map it reads 1.0.
    #[test]
    fn complement_experiment_q_reaches_off_stratum_survivors() {
        // Stratum: base 1. Complement: bases 3, 4, 5, three strong targets per file, so each
        // complement target's run q is (0+1)/3 in each file.
        let file = |shift: f64| -> (Vec<u32>, Vec<f64>) {
            let mut ids = Vec::new();
            let mut scores = Vec::new();
            for (bid, t) in [(1u32, 9.0), (3, 6.0), (4, 5.0), (5, 4.0)] {
                ids.extend([bid, bid | DECOY]);
                scores.extend([t - shift, 0.1]);
            }
            (ids, scores)
        };
        let files = [file(0.0), file(0.5)];
        let survivors = vec![(0usize, 1u32), (0, 3), (1, 3)];
        let overrides: HashMap<(usize, u32), f64> = HashMap::new();
        let stratum: HashSet<u32> = [1u32].into_iter().collect();

        let run = |complement: Option<&HashMap<u32, f64>>| {
            let read = |i: usize| files[i].clone();
            compute_full_population_fdr_streaming(
                2,
                read,
                &overrides,
                &survivors,
                Some(&stratum),
                complement,
            )
        };

        let below_floor: HashMap<u32, f64> = [(3u32, 0.25)].into_iter().collect();
        let (run_q, exp_q, _, _) = run(Some(&below_floor));
        assert_eq!(run_q[&(0, 3)], 1.0 / 3.0);
        assert_eq!(exp_q[&(0, 3)], 1.0 / 3.0, "floored up to the best run q");

        let above_floor: HashMap<u32, f64> = [(3u32, 0.5)].into_iter().collect();
        let (_, exp_q, _, _) = run(Some(&above_floor));
        assert_eq!(
            exp_q[&(1, 3)],
            0.5,
            "the complement's own q when above the floor"
        );

        // No complement map: the loser default 1.0.
        let (_, exp_q, _, _) = run(None);
        assert_eq!(exp_q[&(0, 3)], 1.0);
    }

    /// The complement's experiment competition excludes stratum base_ids, takes the best
    /// aggregate per side per base_id, gives ties to the decoy, and keys q by the winner's full
    /// entry_id (the loser is absent).
    #[test]
    fn complement_experiment_q_competes_off_stratum_only() {
        let stratum: HashSet<u32> = [1u32].into_iter().collect();
        let records = vec![
            (1u32, 100.0),     // stratum: excluded
            (1 | DECOY, 50.0), // stratum: excluded
            (2, 9.0),
            (2, 3.0), // a lower aggregate for the same entry does not replace the best
            (2 | DECOY, 1.0),
            (6, 8.5),
            (7, 8.4),
            (8, 8.3),
            (3, 8.0),
            (3 | DECOY, 8.0), // tie -> decoy wins
            (4, 7.0),
            (5 | DECOY, 6.0), // unpaired decoy wins its base_id
        ];
        let q = compute_complement_experiment_q(records, &stratum);

        assert!(!q.contains_key(&1u32) && !q.contains_key(&(1 | DECOY)));
        // Sorted winners: T2 9.0, T6, T7, T8, D3 8.0, T4 7.0, D5 6.0.
        // Forward (d+1)/t: 1, 1/2, 1/3, 1/4, 2/4, 2/5, 3/5 -> backward min:
        // .25 .25 .25 .25 .4 .4 .6
        assert_eq!(q.get(&2u32).copied(), Some(0.25));
        assert_eq!(q.get(&8u32).copied(), Some(0.25));
        assert_eq!(q.get(&(3 | DECOY)).copied(), Some(0.4));
        assert_eq!(q.get(&3u32), None, "the target lost its tie");
        assert_eq!(q.get(&4u32).copied(), Some(0.4));
        assert_eq!(q.get(&(5 | DECOY)).copied(), Some(0.6));
        assert_eq!(q.get(&(2 | DECOY)), None);
    }
}
