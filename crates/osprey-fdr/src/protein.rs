//! Protein parsimony and picked-protein FDR control
//!
//! Implements native Rust protein inference:
//! - Protein grouping: proteins with identical peptide sets are merged
//! - Subset elimination: groups whose peptides are a strict subset of another are removed
//! - Shared peptide handling: All (default), Razor, or Unique modes
//! - Picked-protein FDR: target-decoy competition at the protein group level
//!
//! Protein information comes from the spectral library (`LibraryEntry.protein_ids`),
//! not from a FASTA file. Decoy proteins are identified by the `DECOY_` prefix
//! added by `DecoyGenerator`.

use osprey_core::config::SharedPeptideMode;
use osprey_core::types::{FdrEntry, LibraryEntry};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::sync::Arc;

/// Index type for protein groups
pub type ProteinGroupId = u32;

/// A protein group: proteins sharing identical peptide sets after parsimony
#[derive(Debug, Clone)]
pub struct ProteinGroup {
    /// Unique ID for this group
    pub id: ProteinGroupId,
    /// Protein accessions in this group (identical peptide sets)
    pub accessions: Vec<String>,
    /// Peptide sequences unique to this group (modified_sequence)
    pub unique_peptides: HashSet<String>,
    /// Peptide sequences shared with other groups (modified_sequence)
    pub shared_peptides: HashSet<String>,
}

/// Result of protein parsimony analysis
#[derive(Debug)]
pub struct ProteinParsimonyResult {
    /// Protein groups after parsimony
    pub groups: Vec<ProteinGroup>,
    /// Peptide (modified_sequence) -> protein group ID(s) mapping
    pub peptide_to_groups: HashMap<String, Vec<ProteinGroupId>>,
}

/// Result of picked-protein FDR computation.
///
/// Only target winners appear in `group_qvalues` and `group_scores`. Decoy winners
/// are statistical machinery for the cumulative FDR computation and are not
/// exposed. Protein-level posterior error probability (PEP) is intentionally
/// not computed — use peptide-level PEP instead for downstream confidence.
#[derive(Debug)]
pub struct ProteinFdrResult {
    /// Protein group ID -> q-value (target winners only)
    pub group_qvalues: HashMap<ProteinGroupId, f64>,
    /// Protein group ID -> best peptide SVM score (target side, winners only)
    pub group_scores: HashMap<ProteinGroupId, f64>,
    /// Peptide (modified_sequence) -> best protein q-value among its groups
    pub peptide_qvalues: HashMap<String, f64>,
}

/// Decoy protein prefix used by DecoyGenerator
const DECOY_PREFIX: &str = "DECOY_";

/// Build protein parsimony from the spectral library.
///
/// Groups proteins with identical peptide sets, eliminates subsets, and
/// classifies peptides as unique or shared. Only target entries are used.
///
/// If `detected_peptides` is provided, only library entries whose
/// modified_sequence is in the detected set are included. This ensures
/// protein groups reflect what was actually observed rather than the
/// full theoretical library, since two proteins may be distinguishable
/// in theory but indistinguishable based on the peptides that were detected.
pub fn build_protein_parsimony(
    library: &[LibraryEntry],
    mode: SharedPeptideMode,
    detected_peptides: Option<&HashSet<String>>,
) -> ProteinParsimonyResult {
    // Step 1: Build bipartite graph from target entries only
    let mut peptide_to_proteins: HashMap<String, HashSet<String>> = HashMap::new();
    let mut protein_to_peptides: HashMap<String, HashSet<String>> = HashMap::new();

    for entry in library {
        if entry.is_decoy {
            continue;
        }
        // If a detected set is provided, skip undetected peptides
        if let Some(detected) = detected_peptides {
            if !detected.contains(&entry.modified_sequence) {
                continue;
            }
        }
        for protein in &entry.protein_ids {
            peptide_to_proteins
                .entry(entry.modified_sequence.clone())
                .or_default()
                .insert(protein.clone());
            protein_to_peptides
                .entry(protein.clone())
                .or_default()
                .insert(entry.modified_sequence.clone());
        }
    }

    log::info!(
        "Protein parsimony: {} proteins, {} peptides in bipartite graph (peptides from precursor-level FDR)",
        protein_to_peptides.len(),
        peptide_to_proteins.len()
    );

    // Step 2: Group proteins with identical peptide sets
    let mut peptide_set_to_accessions: HashMap<BTreeSet<String>, Vec<String>> = HashMap::new();
    for (protein, peptides) in &protein_to_peptides {
        let key: BTreeSet<String> = peptides.iter().cloned().collect();
        peptide_set_to_accessions
            .entry(key)
            .or_default()
            .push(protein.clone());
    }

    let mut groups: Vec<(BTreeSet<String>, Vec<String>)> =
        peptide_set_to_accessions.into_iter().collect();

    // Sort by peptide count descending for subset elimination
    groups.sort_by_key(|g| std::cmp::Reverse(g.0.len()));

    log::info!(
        "Protein parsimony: {} protein groups after identical-set merging",
        groups.len()
    );

    // Step 3: Subset elimination
    let mut retained: Vec<(BTreeSet<String>, Vec<String>)> = Vec::new();
    for (peptide_set, accessions) in groups {
        let is_subset = retained
            .iter()
            .any(|(larger_set, _)| peptide_set.is_subset(larger_set) && peptide_set != *larger_set);
        if !is_subset {
            retained.push((peptide_set, accessions));
        }
    }

    log::info!(
        "Protein parsimony: {} protein groups after subset elimination",
        retained.len()
    );

    // Step 4: Assign group IDs and build peptide -> group mapping
    let mut result_groups: Vec<ProteinGroup> = Vec::with_capacity(retained.len());
    let mut peptide_to_groups: HashMap<String, Vec<ProteinGroupId>> = HashMap::new();

    for (idx, (peptide_set, accessions)) in retained.into_iter().enumerate() {
        let gid = idx as ProteinGroupId;
        for peptide in &peptide_set {
            peptide_to_groups
                .entry(peptide.clone())
                .or_default()
                .push(gid);
        }
        result_groups.push(ProteinGroup {
            id: gid,
            accessions,
            unique_peptides: HashSet::new(),
            shared_peptides: HashSet::new(),
        });
    }

    // Classify peptides as unique or shared
    for (peptide, group_ids) in &peptide_to_groups {
        if group_ids.len() == 1 {
            result_groups[group_ids[0] as usize]
                .unique_peptides
                .insert(peptide.clone());
        } else {
            for &gid in group_ids {
                result_groups[gid as usize]
                    .shared_peptides
                    .insert(peptide.clone());
            }
        }
    }

    let n_unique: usize = result_groups.iter().map(|g| g.unique_peptides.len()).sum();
    let n_shared: usize = peptide_to_groups
        .values()
        .filter(|gids| gids.len() > 1)
        .count();
    log::info!(
        "Protein parsimony: {} unique peptides, {} shared peptides",
        n_unique,
        n_shared
    );

    // Step 5: Apply shared peptide mode
    match mode {
        SharedPeptideMode::All => {
            // No reassignment needed; shared peptides remain mapped to all groups
        }
        SharedPeptideMode::Razor => {
            // Proper iterative greedy set cover:
            //   1. Among all remaining shared peptides, find the group with the
            //      most unique peptides (across all groups that still have any
            //      shared peptides to assign).
            //   2. Assign ALL shared peptides belonging to that group to it.
            //   3. Remove those peptides from the shared sets of all other groups.
            //   4. Repeat until no shared peptides remain.
            //
            // This is order-independent: the "winner" is determined globally at
            // each step, not by iteration order over a HashMap. Tiebreak between
            // groups with equal unique counts is lowest group ID (deterministic).

            // Collect all shared peptides up front (deterministic order by sorting)
            let mut shared_peptides: Vec<String> = peptide_to_groups
                .iter()
                .filter(|(_, gids)| gids.len() > 1)
                .map(|(pep, _)| pep.clone())
                .collect();
            shared_peptides.sort();

            // Track which peptides still need assignment (all shared at start)
            let mut unassigned: HashSet<String> = shared_peptides.into_iter().collect();

            while !unassigned.is_empty() {
                // Find the group with the most unique peptides that still has
                // at least one unassigned shared peptide.
                // Tiebreak: lowest group ID.
                let best_gid: Option<ProteinGroupId> = result_groups
                    .iter()
                    .filter(|g| {
                        g.shared_peptides
                            .iter()
                            .any(|p| unassigned.contains(p.as_str()))
                    })
                    .max_by_key(|g| (g.unique_peptides.len(), std::cmp::Reverse(g.id)))
                    .map(|g| g.id);

                let best_gid = match best_gid {
                    Some(id) => id,
                    None => break, // no groups with unassigned shared peptides
                };

                // Collect the peptides this group will claim in this round.
                // Sort for deterministic processing order.
                let mut claimed: Vec<String> = result_groups[best_gid as usize]
                    .shared_peptides
                    .iter()
                    .filter(|p| unassigned.contains(p.as_str()))
                    .cloned()
                    .collect();
                claimed.sort();

                for peptide in &claimed {
                    // Remove from all groups' shared sets
                    if let Some(group_ids) = peptide_to_groups.get(peptide).cloned() {
                        for gid in &group_ids {
                            result_groups[*gid as usize].shared_peptides.remove(peptide);
                        }
                    }
                    // Add to the winning group's unique set
                    result_groups[best_gid as usize]
                        .unique_peptides
                        .insert(peptide.clone());
                    // Update peptide_to_groups to point only to the winning group
                    peptide_to_groups.insert(peptide.clone(), vec![best_gid]);
                    // Mark as assigned
                    unassigned.remove(peptide);
                }
            }
        }
        SharedPeptideMode::Unique => {
            // Drop shared peptides from the mapping entirely
            let shared_peptides: Vec<String> = peptide_to_groups
                .iter()
                .filter(|(_, gids)| gids.len() > 1)
                .map(|(pep, _)| pep.clone())
                .collect();

            for peptide in &shared_peptides {
                let group_ids = peptide_to_groups.remove(peptide).unwrap_or_default();
                for gid in group_ids {
                    result_groups[gid as usize]
                        .shared_peptides
                        .remove(peptide.as_str());
                }
            }
        }
    }

    ProteinParsimonyResult {
        groups: result_groups,
        peptide_to_groups,
    }
}

/// Compute picked-protein FDR (Savitski et al. 2015).
///
/// For each protein group, this algorithm:
///
/// 1. **Score by best peptide**: Computes a target-side score (max SVM score over
///    target peptides passing the gate) and a decoy-side score (max SVM score over
///    DECOY_-prefixed peptides passing the gate). Uses the **single best peptide**
///    per side — not a sum — because sum aggregation is length-biased (Savitski
///    explicitly tested and rejected it).
///
/// 2. **Pairwise picking**: Each protein group produces exactly one winner. If
///    `target_score >= decoy_score`, the target wins; otherwise the decoy wins.
///    Groups with only a target side win as target; groups with only a decoy side
///    win as decoy. Groups with no passing peptides are skipped.
///
/// 3. **Classical cumulative FDR on winners**: Sorted by score descending (tiebreak
///    group_id ascending), `q = cum_decoys / max(1, cum_targets)` at each position.
///    A backward sweep enforces monotonicity (lower score → non-decreasing q-value).
///
/// 4. **Target winners only**: `group_qvalues` contains only target winners. Decoy
///    winners are statistical machinery for the FDR computation and are not
///    exposed to downstream code.
///
/// 5. **Peptide propagation**: Each peptide's q-value is the best (lowest) q-value
///    among the protein groups it belongs to (important for shared peptides in
///    `SharedPeptideMode::All`).
///
/// ## Why Picked-Protein?
///
/// Classical protein-level TDC suffers from decoy over-representation: as genuine
/// targets dominate, random decoy matches accumulate disproportionately in the
/// low-scoring region. Pairwise picking eliminates this through **structural
/// symmetry** — each target-decoy pair produces exactly one winner, so the pool
/// of winners is balanced by construction and classical cumulative FDR on the
/// winner list is well-calibrated.
///
/// ## Parameters
///
/// - `parsimony`: protein parsimony result (peptide-to-group mapping).
/// - `best_scores`: per-peptide best SVM score and best peptide-level q-value
///   across all files (see [`collect_best_peptide_scores`]).
/// - `qvalue_gate`: only **target** peptides with `best_qvalue <= qvalue_gate`
///   contribute (Savitski uses 0.01 = `run_fdr`). Decoy peptides are **not**
///   gated — they form the null distribution and must be unbiased.
///
/// ## Scoring: SVM discriminant (after trying q-value and PEP)
///
/// Protein scores are the **maximum peptide SVM discriminant** across the
/// group's peptides. We arrived at SVM via two failed alternatives:
///
/// - **Peptide q-value** (attempted): ranking by `min(run_peptide_qvalue)`
///   collapsed the decoy null distribution. About 99% of decoy peptides have
///   `q = 1.0` because they lost peptide-level TDC, so almost no decoy
///   protein could compete. On the Stellar 3-file HeLa dataset this produced
///   only 2 decoy winners out of 6102 — a ~0.03% decoy rate, i.e. the
///   cumulative FDR pinned at near-zero and every target trivially passed
///   1% FDR. Under-calibrated and useless as a gate.
///
/// - **Peptide PEP** (attempted): ranking by `min(pep)` had the same
///   failure mode for a subtler reason. Osprey's [`PepEstimator`] is fit on
///   TDC winners only — its bins cover `[min_winner_score, max_winner_score]`
///   and `posterior_error()` clamps any out-of-range query to the nearest
///   bin. Losing decoys have SVM scores below the winner range, so they all
///   clamp to `bins[0] ≈ 1.0`. Fitting a second estimator over all entries
///   was considered but adds a separate PEP model just for protein FDR and
///   doesn't meaningfully differ from using raw SVM scores (PEP is a
///   monotone function of score).
///
/// - **SVM discriminant** (kept): every entry, winner or loser, target or
///   decoy, has a well-defined SVM score on the same scale. On Stellar we
///   see 348 decoy winners out of 5988 (~5.8% decoy rate) — a real overlap
///   between target and decoy score distributions, correctly calibrated
///   cumulative FDR, and 5634 target groups at 1% FDR. The earlier concern
///   that "SVM scores separate too cleanly at the extreme tails" turned out
///   to be unfounded on real data: the tail overlap is enough to produce a
///   meaningful null.
///
/// The target-side gate is still based on peptide q-value
/// (`best_qvalue <= qvalue_gate`) so the set of "reportable" target proteins
/// matches Savitski's convention. Only the ranking within that gate uses
/// SVM scores.
///
/// ## Why the decoy side is not gated
///
/// A naive implementation would apply the same q-value gate to both sides for
/// symmetry. That is wrong. In target-decoy competition, decoy peptides typically
/// have high q-values because they lose the competition. Gating decoys at
/// `q <= 0.01` leaves only the handful of decoys that happened to outscore their
/// target in peptide-level TDC — a severe survivorship bias that eliminates
/// almost the entire decoy null distribution.
///
/// The target gate (Savitski's convention) restricts analysis to "proteins we'd
/// actually report". The decoy side uses the full decoy peptide pool as the
/// null distribution, which is what makes picked-protein FDR well-calibrated.
///
/// ## Reference
/// Savitski MM, Wilhelm M, Hahne H, Kuster B, Bantscheff M. "A Scalable Approach
/// for Protein False Discovery Rate Estimation in Large Proteomic Data Sets."
/// Mol Cell Proteomics. 2015;14(9):2394-2404.
pub fn compute_protein_fdr(
    parsimony: &ProteinParsimonyResult,
    best_scores: &HashMap<Arc<str>, PeptideScore>,
    qvalue_gate: f64,
) -> ProteinFdrResult {
    // Step 1: For each protein group, compute target and decoy best-peptide
    // SVM scores (maximum, i.e. strongest). Ranking is by SVM discriminant;
    // higher = better.
    // - Target side: only target peptides with best_qvalue <= qvalue_gate
    //   contribute (restricts analysis to "proteins we'd actually report");
    //   the gate uses peptide q-value per Savitski's convention, but ranking
    //   among eligible targets uses the raw SVM score.
    // - Decoy side: ALL decoy peptides contribute (forms the null distribution).
    //   Gating decoys would create survivorship bias — see doc comment above.
    #[derive(Default, Clone, Copy)]
    struct GroupScore {
        target_score: Option<f64>,
        decoy_score: Option<f64>,
    }
    let mut group_scores_map: HashMap<ProteinGroupId, GroupScore> = HashMap::new();

    for (peptide, group_ids) in &parsimony.peptide_to_groups {
        // Target side: peptide sequence directly, GATED on peptide q-value,
        // RANKED by SVM score.
        if let Some(ps) = best_scores.get(peptide.as_str()) {
            if !ps.is_decoy && ps.best_qvalue <= qvalue_gate {
                for &gid in group_ids {
                    let gs = group_scores_map.entry(gid).or_default();
                    gs.target_score = Some(match gs.target_score {
                        Some(s) => s.max(ps.score),
                        None => ps.score,
                    });
                }
            }
        }
        // Decoy side: DECOY_-prefixed peptide sequence, NOT GATED.
        // All decoy peptides contribute their SVM score to the null distribution.
        let decoy_key = format!("{}{}", DECOY_PREFIX, peptide);
        if let Some(ps) = best_scores.get(decoy_key.as_str()) {
            if ps.is_decoy {
                for &gid in group_ids {
                    let gs = group_scores_map.entry(gid).or_default();
                    gs.decoy_score = Some(match gs.decoy_score {
                        Some(s) => s.max(ps.score),
                        None => ps.score,
                    });
                }
            }
        }
    }

    // Step 2: Pair picking. Each protein group produces exactly one winner.
    // Target wins if target_score >= decoy_score (ties go to target).
    #[derive(Clone, Copy)]
    struct Winner {
        group_id: ProteinGroupId,
        score: f64, // best peptide SVM score for this winner
        is_decoy: bool,
    }
    let mut winners: Vec<Winner> = Vec::new();
    // Iterate in deterministic group order to avoid HashMap iteration leak
    for group in &parsimony.groups {
        let gs = match group_scores_map.get(&group.id) {
            Some(x) => *x,
            None => continue, // no peptides passed the gate for this group
        };
        match (gs.target_score, gs.decoy_score) {
            (Some(t), Some(d)) => {
                if t >= d {
                    winners.push(Winner {
                        group_id: group.id,
                        score: t,
                        is_decoy: false,
                    });
                } else {
                    winners.push(Winner {
                        group_id: group.id,
                        score: d,
                        is_decoy: true,
                    });
                }
            }
            (Some(t), None) => winners.push(Winner {
                group_id: group.id,
                score: t,
                is_decoy: false,
            }),
            (None, Some(d)) => winners.push(Winner {
                group_id: group.id,
                score: d,
                is_decoy: true,
            }),
            (None, None) => {}
        }
    }

    // Step 3: Classical cumulative FDR on winners, sorted by SVM score
    // DESCENDING (highest = strongest = first). Tiebreak by group_id
    // ascending for determinism.
    winners.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then(a.group_id.cmp(&b.group_id))
    });

    let mut raw_qvalues: Vec<f64> = Vec::with_capacity(winners.len());
    let mut cum_targets = 0usize;
    let mut cum_decoys = 0usize;
    for w in &winners {
        if w.is_decoy {
            cum_decoys += 1;
        } else {
            cum_targets += 1;
        }
        let q = if cum_targets > 0 {
            (cum_decoys as f64 / cum_targets as f64).min(1.0)
        } else {
            1.0
        };
        raw_qvalues.push(q);
    }

    // Step 4: Backward sweep for monotonicity. Build group_qvalues with
    // target winners only.
    let mut group_qvalues: HashMap<ProteinGroupId, f64> = HashMap::new();
    let mut group_scores: HashMap<ProteinGroupId, f64> = HashMap::new();
    let mut min_q = 1.0f64;
    for i in (0..winners.len()).rev() {
        min_q = min_q.min(raw_qvalues[i]);
        let w = &winners[i];
        if !w.is_decoy {
            group_qvalues.insert(w.group_id, min_q);
            group_scores.insert(w.group_id, w.score);
        }
    }

    // Step 5: Propagate to peptides. Each peptide's protein q-value is the best
    // (lowest) across protein groups it belongs to. Peptides belonging only to
    // groups that lost the pair (decoy won, or had no winner) get q = 1.0.
    let mut peptide_qvalues: HashMap<String, f64> = HashMap::new();
    for (peptide, group_ids) in &parsimony.peptide_to_groups {
        let best_q = group_ids
            .iter()
            .filter_map(|gid| group_qvalues.get(gid))
            .copied()
            .fold(1.0f64, f64::min);
        peptide_qvalues.insert(peptide.clone(), best_q);
    }

    let n_target_winners = winners.iter().filter(|w| !w.is_decoy).count();
    let n_decoy_winners = winners.iter().filter(|w| w.is_decoy).count();
    let n_at_1pct = group_qvalues.values().filter(|&&q| q <= 0.01).count();
    log::info!(
        "Picked-protein FDR: {} winners ({} target, {} decoy); {} target groups at 1% FDR",
        winners.len(),
        n_target_winners,
        n_decoy_winners,
        n_at_1pct,
    );

    ProteinFdrResult {
        group_qvalues,
        group_scores,
        peptide_qvalues,
    }
}

/// Write a protein-level CSV report.
///
/// Each row is a protein group that passes the FDR threshold. Columns include
/// protein accessions, gene names, q-value, PEP, best peptide score, number of
/// peptides, and the contributing peptide sequences.
///
/// Gene names are looked up from the library entries. If multiple accessions in a
/// group have different gene names, they are joined with semicolons.
pub fn write_protein_report(
    path: &std::path::Path,
    parsimony: &ProteinParsimonyResult,
    fdr_result: &ProteinFdrResult,
    protein_fdr_threshold: f64,
    library: &[LibraryEntry],
) -> std::io::Result<()> {
    use std::io::Write;

    // Build accession -> gene name lookup from library
    let mut accession_to_genes: HashMap<&str, HashSet<&str>> = HashMap::new();
    for entry in library {
        if entry.is_decoy {
            continue;
        }
        for (i, acc) in entry.protein_ids.iter().enumerate() {
            if let Some(gene) = entry.gene_names.get(i) {
                if !gene.is_empty() {
                    accession_to_genes
                        .entry(acc.as_str())
                        .or_default()
                        .insert(gene.as_str());
                }
            }
        }
    }

    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "protein_group,protein_accessions,gene_names,protein_qvalue,best_peptide_score,n_unique_peptides,n_shared_peptides,unique_peptides,shared_peptides"
    )?;

    // Collect passing groups, sorted by q-value
    let mut passing_groups: Vec<&ProteinGroup> = parsimony
        .groups
        .iter()
        .filter(|g| {
            fdr_result
                .group_qvalues
                .get(&g.id)
                .is_some_and(|&q| q <= protein_fdr_threshold)
        })
        .collect();
    passing_groups.sort_by(|a, b| {
        let qa = fdr_result.group_qvalues.get(&a.id).copied().unwrap_or(1.0);
        let qb = fdr_result.group_qvalues.get(&b.id).copied().unwrap_or(1.0);
        qa.total_cmp(&qb)
    });

    for group in &passing_groups {
        let q = fdr_result
            .group_qvalues
            .get(&group.id)
            .copied()
            .unwrap_or(1.0);
        let score = fdr_result
            .group_scores
            .get(&group.id)
            .copied()
            .unwrap_or(0.0);
        let accessions = group.accessions.join(";");

        // Collect gene names for all accessions in this group
        let mut genes: Vec<&str> = group
            .accessions
            .iter()
            .flat_map(|acc| {
                accession_to_genes
                    .get(acc.as_str())
                    .into_iter()
                    .flatten()
                    .copied()
            })
            .collect();
        genes.sort();
        genes.dedup();
        let gene_str = genes.join(";");

        let unique_peps: Vec<&String> = {
            let mut v: Vec<&String> = group.unique_peptides.iter().collect();
            v.sort();
            v
        };
        let shared_peps: Vec<&String> = {
            let mut v: Vec<&String> = group.shared_peptides.iter().collect();
            v.sort();
            v
        };
        let unique_str = unique_peps
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join(";");
        let shared_str = shared_peps
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join(";");

        writeln!(
            file,
            "{},{},{},{:.6},{:.4},{},{},{},{}",
            group.id,
            accessions,
            gene_str,
            q,
            score,
            group.unique_peptides.len(),
            group.shared_peptides.len(),
            unique_str,
            shared_str,
        )?;
    }

    log::info!(
        "Wrote protein report: {} protein groups to {}",
        passing_groups.len(),
        path.display()
    );

    Ok(())
}

/// Peptide-level data for protein FDR scoring.
#[derive(Debug, Clone)]
pub struct PeptideScore {
    /// Best (maximum) SVM discriminant score for this peptide across all files.
    /// Used as the **ranking** score in picked-protein FDR (max across group members).
    pub score: f64,
    /// Whether this peptide is a decoy
    pub is_decoy: bool,
    /// Best (lowest) run-level **peptide** q-value for this peptide across all files.
    /// Used as the target-side **gate** in picked-protein FDR (Savitski 2015 convention):
    /// only targets with `best_qvalue <= gate` are eligible to be winners.
    pub best_qvalue: f64,
}

/// Collect peptide-level scores for protein FDR.
///
/// For each unique modified_sequence, keeps the best SVM score and best peptide-level
/// q-value across all files. Can be called at two points in the pipeline:
///
/// 1. **Before compaction** (first-pass protein FDR): uses first-pass SVM scores and
///    first-pass peptide q-values. The full target + decoy pool is present, so
///    picked-protein has symmetric competition.
///
/// 2. **After second-pass FDR** (authoritative protein FDR): uses reconciliation-
///    corrected second-pass SVM scores and second-pass peptide q-values. The pool
///    is restricted to compacted entries (first-pass passing peptides + their
///    paired decoys + any rescued by protein-aware compaction).
///
/// Uses `entry.run_peptide_qvalue` (not precursor q-value) because picked-protein
/// gates on peptide-level FDR per Savitski 2015.
///
/// Returns a map from modified_sequence -> PeptideScore.
pub fn collect_best_peptide_scores(
    per_file_entries: &[(String, Vec<FdrEntry>)],
) -> HashMap<Arc<str>, PeptideScore> {
    let mut best: HashMap<Arc<str>, PeptideScore> = HashMap::new();
    for (_, entries) in per_file_entries {
        for entry in entries {
            best.entry(entry.modified_sequence.clone())
                .and_modify(|ps| {
                    if entry.score > ps.score {
                        ps.score = entry.score;
                    }
                    if entry.run_peptide_qvalue < ps.best_qvalue {
                        ps.best_qvalue = entry.run_peptide_qvalue;
                    }
                })
                .or_insert(PeptideScore {
                    score: entry.score,
                    is_decoy: entry.is_decoy,
                    best_qvalue: entry.run_peptide_qvalue,
                });
        }
    }

    let n_target = best.values().filter(|ps| !ps.is_decoy).count();
    let n_decoy = best.values().filter(|ps| ps.is_decoy).count();
    log::info!(
        "Peptide scores for protein FDR: {} target peptides, {} decoy peptides",
        n_target,
        n_decoy,
    );

    best
}

/// Propagate protein q-values to FdrEntry stubs.
///
/// Sets `run_protein_qvalue` and/or `experiment_protein_qvalue` on each entry
/// based on the protein FDR results.
pub fn propagate_protein_qvalues(
    per_file_entries: &mut [(String, Vec<FdrEntry>)],
    protein_fdr: &ProteinFdrResult,
    set_run: bool,
    set_experiment: bool,
) {
    for (_, entries) in per_file_entries.iter_mut() {
        for entry in entries.iter_mut() {
            let q = protein_fdr
                .peptide_qvalues
                .get(&*entry.modified_sequence)
                .copied()
                .unwrap_or(1.0);
            if set_run {
                entry.run_protein_qvalue = q;
            }
            if set_experiment {
                entry.experiment_protein_qvalue = q;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_lib_entry(
        id: u32,
        modified_sequence: &str,
        protein_ids: Vec<&str>,
        is_decoy: bool,
    ) -> LibraryEntry {
        LibraryEntry {
            id,
            sequence: modified_sequence.to_string(),
            modified_sequence: modified_sequence.to_string(),
            modifications: vec![],
            charge: 2,
            precursor_mz: 500.0,
            retention_time: 30.0,
            rt_calibrated: false,
            fragments: vec![],
            protein_ids: protein_ids.into_iter().map(|s| s.to_string()).collect(),
            gene_names: vec![],
            is_decoy,
        }
    }

    #[test]
    fn test_basic_parsimony_grouping() {
        // Three proteins: P1 has {A, B, C}, P2 has {A, B, C} (identical), P3 has {D, E}
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1", "P2"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1", "P2"], false),
            make_lib_entry(3, "PEPTIDEC", vec!["P1", "P2"], false),
            make_lib_entry(4, "PEPTIDED", vec!["P3"], false),
            make_lib_entry(5, "PEPTIDEE", vec!["P3"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        // P1 and P2 should be merged into one group
        assert_eq!(result.groups.len(), 2);

        // One group has 3 peptides, other has 2
        let mut sizes: Vec<usize> = result
            .groups
            .iter()
            .map(|g| g.unique_peptides.len() + g.shared_peptides.len())
            .collect();
        sizes.sort();
        assert_eq!(sizes, vec![2, 3]);
    }

    #[test]
    fn test_subset_elimination() {
        // P1 has {A, B, C}, P2 has {A, B} (subset of P1) -> P2 eliminated
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1", "P2"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1", "P2"], false),
            make_lib_entry(3, "PEPTIDEC", vec!["P1"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        // Only P1 should remain (P2's peptide set is a strict subset)
        assert_eq!(result.groups.len(), 1);
        assert!(result.groups[0].accessions.contains(&"P1".to_string()));
    }

    #[test]
    fn test_shared_peptides_all_mode() {
        // P1 has {A, B, shared}, P2 has {C, D, shared}
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1"], false),
            make_lib_entry(3, "SHARED", vec!["P1", "P2"], false),
            make_lib_entry(4, "PEPTIDEC", vec!["P2"], false),
            make_lib_entry(5, "PEPTIDED", vec!["P2"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        assert_eq!(result.groups.len(), 2);
        // SHARED should map to both groups
        assert_eq!(result.peptide_to_groups["SHARED"].len(), 2);
    }

    #[test]
    fn test_shared_peptides_razor_mode() {
        // P1 has {A, B, C, shared}, P2 has {D, shared}
        // P1 has more unique peptides, so SHARED should be assigned to P1
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1"], false),
            make_lib_entry(3, "PEPTIDEC", vec!["P1"], false),
            make_lib_entry(4, "SHARED", vec!["P1", "P2"], false),
            make_lib_entry(5, "PEPTIDED", vec!["P2"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::Razor, None);

        assert_eq!(result.groups.len(), 2);
        // SHARED should map to only one group (the one with more unique peptides)
        assert_eq!(result.peptide_to_groups["SHARED"].len(), 1);

        // The group with P1 should have SHARED as unique now
        let p1_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(p1_group.unique_peptides.contains("SHARED"));
        assert!(p1_group.shared_peptides.is_empty());
    }

    #[test]
    fn test_shared_peptides_razor_iterative_greedy() {
        // Test that the greedy set cover is deterministic and path-independent.
        //
        // Setup:
        //   P1 has unique peptides {A, B} and shares {X} with P2, {Y} with P3
        //   P2 has unique peptide {C} and shares {X} with P1
        //   P3 has unique peptide {D} and shares {Y} with P1
        //
        // Initial unique counts: P1=2, P2=1, P3=1
        //
        // Iteration 1: P1 wins (2 unique). Claims both X and Y.
        // After: P1={A, B, X, Y}, P2={C}, P3={D}, no shared remaining.
        //
        // The key property: even though X would have been claimable by P2 and
        // Y by P3 in the old single-pass code (which could pick either order),
        // the iterative greedy deterministically assigns both to P1 because it
        // picks the globally-best group first.
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1"], false),
            make_lib_entry(3, "PEPTIDEX", vec!["P1", "P2"], false),
            make_lib_entry(4, "PEPTIDEY", vec!["P1", "P3"], false),
            make_lib_entry(5, "PEPTIDEC", vec!["P2"], false),
            make_lib_entry(6, "PEPTIDED", vec!["P3"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::Razor, None);

        assert_eq!(result.groups.len(), 3);

        // Both shared peptides should have been razored to P1 (the group with
        // more unique peptides).
        assert_eq!(result.peptide_to_groups["PEPTIDEX"].len(), 1);
        assert_eq!(result.peptide_to_groups["PEPTIDEY"].len(), 1);

        let p1_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(p1_group.unique_peptides.contains("PEPTIDEX"));
        assert!(p1_group.unique_peptides.contains("PEPTIDEY"));

        let p2_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P2".to_string()))
            .unwrap();
        assert!(!p2_group.unique_peptides.contains("PEPTIDEX"));
        assert!(p2_group.shared_peptides.is_empty());

        let p3_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P3".to_string()))
            .unwrap();
        assert!(!p3_group.unique_peptides.contains("PEPTIDEY"));
        assert!(p3_group.shared_peptides.is_empty());
    }

    #[test]
    fn test_shared_peptides_razor_cascading_assignment() {
        // Test a case where the iterative greedy matters for determinism:
        //
        // P1 has unique {A, B, C} and shares {X, Y} with P2 and P3 respectively
        // P2 has unique {D} and shares {X} with P1, {Z} with P3
        // P3 has unique {E} and shares {Y} with P1, {Z} with P2
        //
        // Initial unique counts: P1=3, P2=1, P3=1
        //
        // Iteration 1: P1 wins (3 unique). Claims X and Y.
        // Now: P1 has 5 unique, P2 has 1 unique + shared {Z}, P3 has 1 unique + shared {Z}
        //
        // Iteration 2: P2 and P3 tie at 1 unique. P2 wins by lower group ID (P2 was added before P3 in library).
        // P2 claims Z. Unassigned is now empty.
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1"], false),
            make_lib_entry(3, "PEPTIDEC", vec!["P1"], false),
            make_lib_entry(4, "PEPTIDEX", vec!["P1", "P2"], false),
            make_lib_entry(5, "PEPTIDEY", vec!["P1", "P3"], false),
            make_lib_entry(6, "PEPTIDED", vec!["P2"], false),
            make_lib_entry(7, "PEPTIDEZ", vec!["P2", "P3"], false),
            make_lib_entry(8, "PEPTIDEE", vec!["P3"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::Razor, None);

        // All peptides should be assigned to exactly one group
        for p in ["PEPTIDEX", "PEPTIDEY", "PEPTIDEZ"] {
            assert_eq!(
                result.peptide_to_groups[p].len(),
                1,
                "{} should map to exactly one group",
                p
            );
        }

        // P1 should have claimed X and Y (iteration 1)
        let p1_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(p1_group.unique_peptides.contains("PEPTIDEX"));
        assert!(p1_group.unique_peptides.contains("PEPTIDEY"));
        assert!(p1_group.shared_peptides.is_empty());

        // P2 or P3 should have claimed Z (iteration 2)
        // Whichever one wins, the other should have its shared set empty.
        let p2_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P2".to_string()))
            .unwrap();
        let p3_group = result
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P3".to_string()))
            .unwrap();
        assert!(p2_group.shared_peptides.is_empty());
        assert!(p3_group.shared_peptides.is_empty());

        // Exactly one of P2/P3 has Z as unique
        let p2_has_z = p2_group.unique_peptides.contains("PEPTIDEZ");
        let p3_has_z = p3_group.unique_peptides.contains("PEPTIDEZ");
        assert!(
            p2_has_z ^ p3_has_z,
            "Exactly one of P2/P3 should have PEPTIDEZ"
        );
    }

    #[test]
    fn test_shared_peptides_razor_deterministic() {
        // Run the same parsimony 10 times with the same input and verify
        // the output is byte-identical (no HashMap iteration order leak).
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P1"], false),
            make_lib_entry(3, "SHARED1", vec!["P1", "P2"], false),
            make_lib_entry(4, "SHARED2", vec!["P1", "P2", "P3"], false),
            make_lib_entry(5, "PEPTIDEC", vec!["P2"], false),
            make_lib_entry(6, "PEPTIDED", vec!["P3"], false),
        ];

        let first = build_protein_parsimony(&library, SharedPeptideMode::Razor, None);
        for _ in 0..10 {
            let next = build_protein_parsimony(&library, SharedPeptideMode::Razor, None);
            // Same number of groups
            assert_eq!(first.groups.len(), next.groups.len());
            // Same assignments for all shared peptides
            for (pep, gids) in &first.peptide_to_groups {
                assert_eq!(
                    &next.peptide_to_groups[pep], gids,
                    "Razor assignment must be deterministic for peptide {}",
                    pep
                );
            }
        }
    }

    #[test]
    fn test_shared_peptides_unique_mode() {
        // Same setup, but Unique mode drops SHARED entirely
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "SHARED", vec!["P1", "P2"], false),
            make_lib_entry(3, "PEPTIDEC", vec!["P2"], false),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::Unique, None);

        // SHARED should not be in the mapping
        assert!(!result.peptide_to_groups.contains_key("SHARED"));
    }

    #[test]
    fn test_decoy_entries_excluded_from_parsimony() {
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "DECOY_PEPTIDEA", vec!["DECOY_P1"], true),
        ];

        let result = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        // Only target entries used; decoys excluded
        assert_eq!(result.groups.len(), 1);
        assert_eq!(result.groups[0].accessions, vec!["P1".to_string()]);
    }

    #[test]
    fn test_picked_protein_fdr() {
        // Two protein groups: P1 (target wins) and P2 (decoy wins)
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P2"], false),
        ];

        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        // Best scores: P1's target peptide scores 5.0, P2's target scores 1.0
        // Decoy versions: DECOY_PEPTIDEA scores 2.0, DECOY_PEPTIDEB scores 3.0
        let mut best_scores: HashMap<Arc<str>, PeptideScore> = HashMap::new();
        best_scores.insert(
            Arc::from("PEPTIDEA"),
            PeptideScore {
                score: 5.0,
                is_decoy: false,
                best_qvalue: 0.001,
            },
        );
        best_scores.insert(
            Arc::from("DECOY_PEPTIDEA"),
            PeptideScore {
                score: 2.0,
                is_decoy: true,
                best_qvalue: 0.5,
            },
        );
        best_scores.insert(
            Arc::from("PEPTIDEB"),
            PeptideScore {
                score: 1.0,
                is_decoy: false,
                best_qvalue: 0.005,
            },
        );
        best_scores.insert(
            Arc::from("DECOY_PEPTIDEB"),
            PeptideScore {
                score: 3.0,
                is_decoy: true,
                best_qvalue: 0.3,
            },
        );

        let fdr_result = compute_protein_fdr(&parsimony, &best_scores, 1.0);

        // With composite scoring, both groups have peptides contributing.
        // The exact q-values depend on the composite + best_quality TDC.
        let p1_group = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(fdr_result.group_qvalues.contains_key(&p1_group.id));

        // PEPTIDEA should have a protein q-value
        assert!(fdr_result.peptide_qvalues.contains_key("PEPTIDEA"));
    }

    #[test]
    fn test_peptide_inherits_best_protein_qvalue() {
        // A peptide in two protein groups gets the best (lowest) q-value
        let library = vec![
            make_lib_entry(1, "UNIQUEA", vec!["P1"], false),
            make_lib_entry(2, "SHARED", vec!["P1", "P2"], false),
            make_lib_entry(3, "UNIQUEB", vec!["P2"], false),
        ];

        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        // Both proteins are targets that win; P1 higher scoring
        let mut best_scores: HashMap<Arc<str>, PeptideScore> = HashMap::new();
        best_scores.insert(
            Arc::from("UNIQUEA"),
            PeptideScore {
                score: 10.0,
                is_decoy: false,
                best_qvalue: 0.001,
            },
        );
        best_scores.insert(
            Arc::from("SHARED"),
            PeptideScore {
                score: 8.0,
                is_decoy: false,
                best_qvalue: 0.001,
            },
        );
        best_scores.insert(
            Arc::from("UNIQUEB"),
            PeptideScore {
                score: 3.0,
                is_decoy: false,
                best_qvalue: 0.005,
            },
        );
        best_scores.insert(
            Arc::from("DECOY_UNIQUEA"),
            PeptideScore {
                score: 1.0,
                is_decoy: true,
                best_qvalue: 0.5,
            },
        );
        best_scores.insert(
            Arc::from("DECOY_SHARED"),
            PeptideScore {
                score: 1.0,
                is_decoy: true,
                best_qvalue: 0.5,
            },
        );
        best_scores.insert(
            Arc::from("DECOY_UNIQUEB"),
            PeptideScore {
                score: 1.0,
                is_decoy: true,
                best_qvalue: 0.5,
            },
        );

        let fdr_result = compute_protein_fdr(&parsimony, &best_scores, 1.0);

        // SHARED maps to both groups; should get the better q-value
        let shared_q = fdr_result.peptide_qvalues["SHARED"];
        let p1_group = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        let p1_q = fdr_result.group_qvalues[&p1_group.id];
        assert!(shared_q <= p1_q + f64::EPSILON);
    }

    #[test]
    fn test_protein_with_no_peptides_gets_default_qvalue() {
        // Protein group with no observed peptides gets q-value 1.0
        let library = vec![make_lib_entry(1, "PEPTIDEA", vec!["P1"], false)];

        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        // No scores at all
        let best_scores: HashMap<Arc<str>, PeptideScore> = HashMap::new();
        let fdr_result = compute_protein_fdr(&parsimony, &best_scores, 1.0);

        // PEPTIDEA should get q-value 1.0 (no scores)
        let q = fdr_result
            .peptide_qvalues
            .get("PEPTIDEA")
            .copied()
            .unwrap_or(1.0);
        assert!((q - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn test_collect_best_peptide_scores_uses_svm_score() {
        // collect_best_peptide_scores uses raw SVM score (higher = better)
        let entries = vec![(
            "file1".to_string(),
            vec![
                FdrEntry {
                    entry_id: 1,
                    parquet_index: 0,
                    is_decoy: false,
                    charge: 2,
                    scan_number: 100,
                    apex_rt: 30.0,
                    start_rt: 29.0,
                    end_rt: 31.0,
                    coelution_sum: 5.0,
                    score: 3.0,
                    run_precursor_qvalue: 0.001,
                    run_peptide_qvalue: 0.002,
                    run_protein_qvalue: 1.0,
                    experiment_precursor_qvalue: 1.0,
                    experiment_peptide_qvalue: 1.0,
                    experiment_protein_qvalue: 1.0,
                    pep: 0.001, // good PEP
                    modified_sequence: Arc::from("PEPTIDEA"),
                },
                FdrEntry {
                    entry_id: 1,
                    parquet_index: 0,
                    is_decoy: false,
                    charge: 2,
                    scan_number: 200,
                    apex_rt: 31.0,
                    start_rt: 30.0,
                    end_rt: 32.0,
                    coelution_sum: 6.0,
                    score: 5.0,
                    run_precursor_qvalue: 0.05,
                    run_peptide_qvalue: 0.06,
                    run_protein_qvalue: 1.0,
                    experiment_precursor_qvalue: 1.0,
                    experiment_peptide_qvalue: 1.0,
                    experiment_protein_qvalue: 1.0,
                    pep: 0.20, // worse PEP
                    modified_sequence: Arc::from("PEPTIDEA"),
                },
            ],
        )];

        let best = collect_best_peptide_scores(&entries);
        // Best is the one with higher SVM score (5.0)
        let ps = &best[&Arc::from("PEPTIDEA") as &Arc<str>];
        assert!(
            (ps.score - 5.0).abs() < 1e-10,
            "Expected 5.0, got {}",
            ps.score
        );
        assert!(!ps.is_decoy);
    }

    // ============================================================
    // Picked-Protein FDR (Savitski 2015) regression tests
    // ============================================================

    /// Helper to build a PeptideScore with all the boilerplate.
    fn ps(score: f64, is_decoy: bool, best_qvalue: f64) -> PeptideScore {
        PeptideScore {
            score,
            is_decoy,
            best_qvalue,
        }
    }

    #[test]
    fn test_picked_protein_fdr_target_wins_pair() {
        // Group with target q=0.001 and decoy q=0.01 → target wins
        // (scoring is by peptide q-value; lower = better)
        let library = vec![make_lib_entry(1, "PEPTIDEA", vec!["P1"], false)];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PEPTIDEA"), ps(2.5, false, 0.001));
        scores.insert(Arc::from("DECOY_PEPTIDEA"), ps(1.8, true, 0.01));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let p1 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(
            result.group_qvalues.contains_key(&p1.id),
            "target should win and appear in group_qvalues"
        );
        assert!(
            result.peptide_qvalues["PEPTIDEA"] <= 1.0,
            "peptide should have a q-value (got {})",
            result.peptide_qvalues["PEPTIDEA"]
        );
    }

    #[test]
    fn test_picked_protein_fdr_decoy_wins_pair() {
        // P1: target q=0.008, decoy q=0.001 → decoy wins (decoy had better peptide q)
        // P2: target q=0.001, decoy q=0.005 → target wins
        // Target P1 is NOT in group_qvalues (decoy winners are excluded from output)
        let library = vec![
            make_lib_entry(1, "PEPTIDEA", vec!["P1"], false),
            make_lib_entry(2, "PEPTIDEB", vec!["P2"], false),
        ];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PEPTIDEA"), ps(0.5, false, 0.008));
        scores.insert(Arc::from("DECOY_PEPTIDEA"), ps(1.2, true, 0.001));
        scores.insert(Arc::from("PEPTIDEB"), ps(5.0, false, 0.001));
        scores.insert(Arc::from("DECOY_PEPTIDEB"), ps(0.8, true, 0.005));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let p1 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        let p2 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P2".to_string()))
            .unwrap();

        assert!(
            !result.group_qvalues.contains_key(&p1.id),
            "P1 should not be in group_qvalues (decoy won its pair)"
        );
        assert!(
            result.group_qvalues.contains_key(&p2.id),
            "P2 should be in group_qvalues (target won its pair)"
        );
    }

    #[test]
    fn test_picked_protein_fdr_target_only_wins() {
        // Group with only a target peptide (no decoy match) → target wins
        let library = vec![make_lib_entry(1, "PEPTIDEA", vec!["P1"], false)];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PEPTIDEA"), ps(2.0, false, 0.005));
        // No DECOY_PEPTIDEA in the score map

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let p1 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(
            result.group_qvalues.contains_key(&p1.id),
            "target-only group should win"
        );
    }

    #[test]
    fn test_picked_protein_fdr_target_fails_gate_decoy_wins() {
        // When the target peptide fails the gate but the decoy peptide exists
        // (regardless of its q-value), the decoy side is still populated and
        // the group produces a decoy winner — the target is NOT in group_qvalues.
        //
        // This verifies the asymmetric gating: targets are filtered by q-value
        // but decoys are not, so picking can still run and the null distribution
        // is preserved.
        let library = vec![make_lib_entry(1, "PEPTIDEA", vec!["P1"], false)];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        // Target peptide FAILS the 0.01 gate (q=0.05 > 0.01)
        scores.insert(Arc::from("PEPTIDEA"), ps(5.0, false, 0.05));
        // Decoy peptide has a reasonable q-value and is NOT gated
        scores.insert(Arc::from("DECOY_PEPTIDEA"), ps(1.0, true, 0.02));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let p1 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        assert!(
            !result.group_qvalues.contains_key(&p1.id),
            "target failed gate → group should not appear in group_qvalues (decoy won)"
        );
    }

    #[test]
    fn test_picked_protein_fdr_decoys_not_gated() {
        // Regression test: decoy peptides must NOT be filtered by the q-value gate.
        // Gating decoys creates a survivorship bias that eliminates almost the
        // entire decoy null distribution, producing near-zero decoy winners and
        // a badly calibrated FDR.
        //
        // 10 target proteins with SVM scores [1.0, 1.1, ... 1.9], all passing
        // the q-value gate (q=0.001). Their decoys have raw SVM scores drawn
        // from a mixed distribution — most lose peptide-level TDC (q=1.0) but
        // are still included in the protein FDR null. A couple of decoy SVM
        // scores are larger than the matched target's SVM score, so they
        // should beat the target under SVM-based picked-protein competition.
        // If the decoy side were (incorrectly) gated by q-value, those decoys
        // would be filtered out and every target would win.
        let library: Vec<LibraryEntry> = (0..10)
            .map(|i| {
                make_lib_entry(
                    i + 1,
                    &format!("PEP{}", i),
                    vec![format!("P{}", i).as_str()],
                    false,
                )
            })
            .collect();
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        // All targets: passing gate (q=0.001), SVM scores [1.0, 1.1, ... 1.9]
        for i in 0..10 {
            scores.insert(
                Arc::from(format!("PEP{}", i).as_str()),
                ps(1.0 + (i as f64) * 0.1, false, 0.001),
            );
        }
        // Decoys: mostly high q (lost peptide TDC) — would be filtered by a
        // buggy gate. A few have SVM scores above their matched target:
        //   Decoy P2 SVM=2.5 vs target 1.2 → decoy wins
        //   Decoy P6 SVM=2.0 vs target 1.6 → decoy wins
        //   Decoy P4 SVM=1.15 vs target 1.4 → target still wins
        let decoy_svm_scores = [0.5, 0.7, 2.5, 0.9, 1.15, 0.3, 2.0, 1.05, 0.2, 0.8];
        for (i, &decoy_score) in decoy_svm_scores.iter().enumerate() {
            scores.insert(
                Arc::from(format!("DECOY_PEP{}", i).as_str()),
                ps(decoy_score, true, 1.0),
            );
        }

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let target_winners = result.group_qvalues.len();

        // P2 and P6 lose to their decoys → 8 target winners expected.
        assert!(
            target_winners < 10,
            "at least one decoy should have won its pair (got {} target winners out of 10)",
            target_winners
        );
        assert!(
            target_winners >= 5,
            "most targets should still win (got {} target winners)",
            target_winners
        );
    }

    #[test]
    fn test_picked_protein_fdr_best_not_sum() {
        // A protein with 5 target peptides at SVM scores [0.5, 1.0, 1.5, 2.0, 2.5]
        // uses the MAX (2.5) as its score, not the sum (7.5). This is the key
        // Savitski property: no length bias — a protein with many weak peptides
        // does not accumulate a better score than a protein with one strong
        // peptide.
        let library = vec![
            make_lib_entry(1, "PEP1", vec!["P1"], false),
            make_lib_entry(2, "PEP2", vec!["P1"], false),
            make_lib_entry(3, "PEP3", vec!["P1"], false),
            make_lib_entry(4, "PEP4", vec!["P1"], false),
            make_lib_entry(5, "PEP5", vec!["P1"], false),
        ];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PEP1"), ps(0.5, false, 0.005));
        scores.insert(Arc::from("PEP2"), ps(1.0, false, 0.005));
        scores.insert(Arc::from("PEP3"), ps(1.5, false, 0.005));
        scores.insert(Arc::from("PEP4"), ps(2.0, false, 0.005));
        scores.insert(Arc::from("PEP5"), ps(2.5, false, 0.005));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let p1 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        let score = result.group_scores[&p1.id];
        assert!(
            (score - 2.5).abs() < 1e-10,
            "score should be 2.5 (max SVM), not 7.5 (sum), got {}",
            score
        );
    }

    #[test]
    fn test_picked_protein_fdr_monotonic() {
        // Q-values must be non-decreasing as protein score (peptide q-value) increases
        // — i.e. the strongest protein (lowest q) has the smallest protein q-value.
        let library = vec![
            make_lib_entry(1, "PA", vec!["P1"], false),
            make_lib_entry(2, "PB", vec!["P2"], false),
            make_lib_entry(3, "PC", vec!["P3"], false),
            make_lib_entry(4, "PD", vec!["P4"], false),
        ];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        // 4 targets with ascending peptide q-values (PA is strongest)
        scores.insert(Arc::from("PA"), ps(5.0, false, 0.001));
        scores.insert(Arc::from("PB"), ps(4.0, false, 0.002));
        scores.insert(Arc::from("PC"), ps(3.0, false, 0.003));
        scores.insert(Arc::from("PD"), ps(2.0, false, 0.004));
        // 4 decoys with high q-values — they all lose to their targets
        scores.insert(Arc::from("DECOY_PA"), ps(1.0, true, 0.5));
        scores.insert(Arc::from("DECOY_PB"), ps(1.5, true, 0.5));
        scores.insert(Arc::from("DECOY_PC"), ps(2.5, true, 0.5));
        scores.insert(Arc::from("DECOY_PD"), ps(0.5, true, 0.5));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        // PA (highest score) should have q <= PB q <= PC q <= PD q
        let get_q = |name: &str| -> f64 {
            let g = parsimony
                .groups
                .iter()
                .find(|g| g.accessions.contains(&name.to_string()))
                .unwrap();
            result.group_qvalues.get(&g.id).copied().unwrap_or(1.0)
        };

        let q_pa = get_q("P1");
        let q_pb = get_q("P2");
        // Only check monotonicity where the target was a winner
        assert!(
            q_pa <= q_pb + f64::EPSILON,
            "q-value at higher score (PA={}) should be <= q-value at lower score (PB={})",
            q_pa,
            q_pb
        );
    }

    #[test]
    fn test_picked_protein_fdr_cumulative() {
        // Controlled q-value setup:
        //   P1: target q=0.001, decoy q=1.0   → target wins at q_score=0.001
        //   P2: target q=0.002, decoy q=1.0   → target wins at q_score=0.002
        //   P3: target q=0.005, decoy q=0.003 → decoy wins at q_score=0.003
        //   P4: target q=0.004, decoy q=1.0   → target wins at q_score=0.004
        //
        // Sorted winners ascending by q (lowest = best):
        //   T@0.001, T@0.002, D@0.003, T@0.004
        // Cumulative:
        //   T@0.001: cum_t=1, cum_d=0, raw_q = 0/1 = 0
        //   T@0.002: cum_t=2, cum_d=0, raw_q = 0/2 = 0
        //   D@0.003: cum_t=2, cum_d=1, raw_q = 1/2 = 0.5
        //   T@0.004: cum_t=3, cum_d=1, raw_q = 1/3 ≈ 0.333
        // Backward sweep: T@0.004=0.333, D@0.003=0.333, T@0.002=0, T@0.001=0
        // Target winners in group_qvalues: P1=0, P2=0, P4=0.333, P3 excluded.
        let library = vec![
            make_lib_entry(1, "PA", vec!["P1"], false),
            make_lib_entry(2, "PB", vec!["P2"], false),
            make_lib_entry(3, "PC", vec!["P3"], false),
            make_lib_entry(4, "PD", vec!["P4"], false),
        ];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PA"), ps(10.0, false, 0.001));
        scores.insert(Arc::from("DECOY_PA"), ps(2.0, true, 1.0));
        scores.insert(Arc::from("PB"), ps(8.0, false, 0.002));
        scores.insert(Arc::from("DECOY_PB"), ps(3.0, true, 1.0));
        // P3 target loses: decoy q (0.003) beats target q (0.005)
        scores.insert(Arc::from("PC"), ps(5.0, false, 0.005));
        scores.insert(Arc::from("DECOY_PC"), ps(7.0, true, 0.003));
        scores.insert(Arc::from("PD"), ps(6.0, false, 0.004));
        scores.insert(Arc::from("DECOY_PD"), ps(1.0, true, 1.0));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let get_q = |name: &str| -> f64 {
            let g = parsimony
                .groups
                .iter()
                .find(|g| g.accessions.contains(&name.to_string()))
                .unwrap();
            result.group_qvalues.get(&g.id).copied().unwrap_or(1.0)
        };

        // P1 (T@10) and P2 (T@8) should have q=0 (no decoys above them)
        let q_p1 = get_q("P1");
        let q_p2 = get_q("P2");
        assert!(q_p1 < 0.01, "P1 q should be ~0, got {}", q_p1);
        assert!(q_p2 < 0.01, "P2 q should be ~0, got {}", q_p2);

        // P4 (T@6) should have q = 1/3 = 0.333 (1 decoy above, 3 targets at-or-above)
        let q_p4 = get_q("P4");
        assert!(
            (q_p4 - 1.0 / 3.0).abs() < 0.01,
            "P4 q should be ~0.333, got {}",
            q_p4
        );

        // P3 (T@5, lost to decoy) should NOT be in group_qvalues
        let p3 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P3".to_string()))
            .unwrap();
        assert!(
            !result.group_qvalues.contains_key(&p3.id),
            "P3 should not be in group_qvalues (decoy won its pair)"
        );
    }

    #[test]
    fn test_picked_protein_fdr_deterministic() {
        // Running picked-protein 10 times on the same input produces byte-identical q-values
        let library = vec![
            make_lib_entry(1, "PA", vec!["P1"], false),
            make_lib_entry(2, "PB", vec!["P2"], false),
            make_lib_entry(3, "PC", vec!["P3"], false),
        ];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PA"), ps(3.0, false, 0.001));
        scores.insert(Arc::from("PB"), ps(2.0, false, 0.002));
        scores.insert(Arc::from("PC"), ps(1.0, false, 0.003));
        scores.insert(Arc::from("DECOY_PA"), ps(0.5, true, 0.5));
        scores.insert(Arc::from("DECOY_PB"), ps(0.5, true, 0.5));
        scores.insert(Arc::from("DECOY_PC"), ps(0.5, true, 0.5));

        let first = compute_protein_fdr(&parsimony, &scores, 0.01);
        for _ in 0..10 {
            let next = compute_protein_fdr(&parsimony, &scores, 0.01);
            for (gid, q) in &first.group_qvalues {
                let next_q = next.group_qvalues[gid];
                assert_eq!(
                    q.to_bits(),
                    next_q.to_bits(),
                    "q-value mismatch for group {}: {} vs {}",
                    gid,
                    q,
                    next_q
                );
            }
        }
    }

    #[test]
    fn test_picked_protein_fdr_shared_peptide_contributes_to_all() {
        // Shared peptide in All mode contributes to every protein it belongs to.
        // P1 has unique PA (SVM 1.0) + shared SH (SVM 5.0)
        // P2 has unique PB (SVM 1.0) + shared SH (SVM 5.0)
        // Both proteins should get SH's SVM score of 5.0 as their best (max).
        let library = vec![
            make_lib_entry(1, "PA", vec!["P1"], false),
            make_lib_entry(2, "SH", vec!["P1", "P2"], false),
            make_lib_entry(3, "PB", vec!["P2"], false),
        ];
        let parsimony = build_protein_parsimony(&library, SharedPeptideMode::All, None);

        let mut scores = HashMap::new();
        scores.insert(Arc::from("PA"), ps(1.0, false, 0.005));
        scores.insert(Arc::from("PB"), ps(1.0, false, 0.005));
        scores.insert(Arc::from("SH"), ps(5.0, false, 0.005));
        scores.insert(Arc::from("DECOY_PA"), ps(0.5, true, 0.5));
        scores.insert(Arc::from("DECOY_PB"), ps(0.5, true, 0.5));
        scores.insert(Arc::from("DECOY_SH"), ps(0.5, true, 0.5));

        let result = compute_protein_fdr(&parsimony, &scores, 0.01);

        let p1 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P1".to_string()))
            .unwrap();
        let p2 = parsimony
            .groups
            .iter()
            .find(|g| g.accessions.contains(&"P2".to_string()))
            .unwrap();

        // Both groups should have score 5.0 (max SVM from shared peptide)
        let p1_score = result.group_scores[&p1.id];
        let p2_score = result.group_scores[&p2.id];
        assert!(
            (p1_score - 5.0).abs() < 1e-10,
            "P1 best-peptide SVM should be 5.0 from shared peptide, got {}",
            p1_score
        );
        assert!(
            (p2_score - 5.0).abs() < 1e-10,
            "P2 best-peptide SVM should be 5.0 from shared peptide, got {}",
            p2_score
        );
    }
}
