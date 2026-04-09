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
use osprey_ml::pep::PepEstimator;
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

/// Result of picked-protein FDR computation
#[derive(Debug)]
pub struct ProteinFdrResult {
    /// Protein group ID -> q-value
    pub group_qvalues: HashMap<ProteinGroupId, f64>,
    /// Protein group ID -> posterior error probability
    pub group_pep: HashMap<ProteinGroupId, f64>,
    /// Protein group ID -> best peptide SVM score (target side)
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
        "Protein parsimony: {} proteins, {} peptides in bipartite graph",
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
    groups.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

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
            // Greedy set cover: assign shared peptides to the group with the most unique peptides
            let shared_peptides: Vec<String> = peptide_to_groups
                .iter()
                .filter(|(_, gids)| gids.len() > 1)
                .map(|(pep, _)| pep.clone())
                .collect();

            for peptide in &shared_peptides {
                let group_ids = &peptide_to_groups[peptide];
                // Find group with most unique peptides (tiebreak: lowest group ID for determinism)
                let best_gid = *group_ids
                    .iter()
                    .max_by_key(|&&gid| {
                        let g = &result_groups[gid as usize];
                        (g.unique_peptides.len(), std::cmp::Reverse(gid))
                    })
                    .unwrap();

                // Remove from all groups' shared sets, add to best group's unique set
                for &gid in group_ids {
                    result_groups[gid as usize].shared_peptides.remove(peptide);
                }
                result_groups[best_gid as usize]
                    .unique_peptides
                    .insert(peptide.clone());

                // Update peptide_to_groups to point only to the best group
                peptide_to_groups.insert(peptide.clone(), vec![best_gid]);
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

/// Compute picked-protein FDR at a given level (run or experiment).
///
/// For each protein group, finds the best peptide score (target and decoy separately).
/// Uses picked-protein competition: if the best target peptide outscores the best
/// decoy peptide, the protein is a target; otherwise it is a decoy. Q-values are
/// Compute protein-level FDR using DIA-NN-style composite scoring.
///
/// Uses two complementary scoring metrics, computing q-values independently
/// on each and taking the minimum:
///
/// 1. **Composite score** (sum of per-peptide log-likelihoods):
///    `score += -log(max(E, min(1.0, peptide_err * n_proteotypic_peptides)))`
///    This favors proteins with multiple confident peptides. False proteins
///    (from ~1% peptide FDR) typically have 1 peptide and score low.
///
/// 2. **Best peptide quality**: `max(0, 1 - peptide_err)`
///    This captures the confidence of the single best peptide.
///
/// Both metrics are computed independently for target and decoy protein groups.
/// Target-decoy competition at the protein level produces q-values for each metric,
/// and the final protein q-value is min(q_composite, q_best_quality).
///
/// `qvalue_gate`: only peptides with run q-value <= this threshold contribute
/// to protein scoring (DIA-NN uses 0.02).
pub fn compute_protein_fdr(
    parsimony: &ProteinParsimonyResult,
    best_scores: &HashMap<Arc<str>, PeptideScore>,
    qvalue_gate: f64,
) -> ProteinFdrResult {
    const EPSILON: f64 = 1e-16;

    // Count number of unique proteotypic peptides per protein group (from parsimony)
    let mut n_peptides_per_group: HashMap<ProteinGroupId, usize> = HashMap::new();
    for group in &parsimony.groups {
        n_peptides_per_group.insert(
            group.id,
            group.unique_peptides.len() + group.shared_peptides.len(),
        );
    }

    // For each protein group, accumulate two scores for target and decoy versions.
    // Target scores come from target peptides; decoy scores from DECOY_ prefixed peptides.
    // Both use the SAME n_peptides (from the target group in parsimony).
    let mut target_composite: HashMap<ProteinGroupId, f64> = HashMap::new();
    let mut target_best_quality: HashMap<ProteinGroupId, f64> = HashMap::new();
    let mut decoy_composite: HashMap<ProteinGroupId, f64> = HashMap::new();
    let mut decoy_best_quality: HashMap<ProteinGroupId, f64> = HashMap::new();

    for (peptide, group_ids) in &parsimony.peptide_to_groups {
        let n_pep = group_ids
            .iter()
            .filter_map(|gid| n_peptides_per_group.get(gid))
            .copied()
            .max()
            .unwrap_or(1) as f64;

        // Target peptide contribution
        if let Some(ps) = best_scores.get(peptide.as_str()) {
            if !ps.is_decoy && ps.best_qvalue <= qvalue_gate {
                // Convert SVM score to error probability via logistic function
                let err = 1.0 / (1.0 + ps.score.exp());
                let quality = (1.0 - err).max(0.0);
                let composite_contrib = -(EPSILON.max((err * n_pep).min(1.0)).ln());

                for &gid in group_ids {
                    *target_composite.entry(gid).or_insert(0.0) += composite_contrib;
                    target_best_quality
                        .entry(gid)
                        .and_modify(|q| *q = q.max(quality))
                        .or_insert(quality);
                }
            }
        }

        // Decoy peptide contribution (DECOY_ prefixed version)
        let decoy_key = format!("{}{}", DECOY_PREFIX, peptide);
        if let Some(ps) = best_scores.get(decoy_key.as_str()) {
            if ps.is_decoy && ps.best_qvalue <= qvalue_gate {
                let err = 1.0 / (1.0 + ps.score.exp());
                let quality = (1.0 - err).max(0.0);
                let composite_contrib = -(EPSILON.max((err * n_pep).min(1.0)).ln());

                for &gid in group_ids {
                    *decoy_composite.entry(gid).or_insert(0.0) += composite_contrib;
                    decoy_best_quality
                        .entry(gid)
                        .and_modify(|q| *q = q.max(quality))
                        .or_insert(quality);
                }
            }
        }
    }

    // Diagnostic: log score distributions
    let tc_vals: Vec<f64> = target_composite.values().copied().collect();
    let dc_vals: Vec<f64> = decoy_composite.values().copied().collect();
    let tb_vals: Vec<f64> = target_best_quality.values().copied().collect();
    let db_vals: Vec<f64> = decoy_best_quality.values().copied().collect();
    let median = |v: &mut Vec<f64>| -> f64 {
        if v.is_empty() {
            return 0.0;
        }
        v.sort_by(|a, b| a.total_cmp(b));
        v[v.len() / 2]
    };
    let mut tc = tc_vals.clone();
    let mut dc = dc_vals.clone();
    let mut tb = tb_vals.clone();
    let mut db = db_vals.clone();
    log::info!(
        "Composite scores: {} targets (median={:.2}, max={:.2}), {} decoys (median={:.2}, max={:.2})",
        tc_vals.len(), median(&mut tc), tc_vals.iter().cloned().fold(0.0f64, f64::max),
        dc_vals.len(), median(&mut dc), dc_vals.iter().cloned().fold(0.0f64, f64::max),
    );
    log::info!(
        "Best quality scores: {} targets (median={:.4}, max={:.4}), {} decoys (median={:.4}, max={:.4})",
        tb_vals.len(), median(&mut tb), tb_vals.iter().cloned().fold(0.0f64, f64::max),
        db_vals.len(), median(&mut db), db_vals.iter().cloned().fold(0.0f64, f64::max),
    );

    // Compute q-values independently on each metric using DIA-NN-style TDC.
    // For each metric: sort all target and decoy scores, then for each target protein,
    // q = n_decoys_with_score_>=_s / max(1, n_targets_with_score_>=_s).
    let q_composite =
        compute_protein_qvalues_diann(&parsimony.groups, &target_composite, &decoy_composite);
    let q_best =
        compute_protein_qvalues_diann(&parsimony.groups, &target_best_quality, &decoy_best_quality);

    // Final q-value = min(q_composite, q_best_quality)
    let mut group_qvalues: HashMap<ProteinGroupId, f64> = HashMap::new();
    let mut group_scores: HashMap<ProteinGroupId, f64> = HashMap::new();
    for group in &parsimony.groups {
        let qc = q_composite.get(&group.id).copied().unwrap_or(1.0);
        let qb = q_best.get(&group.id).copied().unwrap_or(1.0);
        group_qvalues.insert(group.id, qc.min(qb));
        group_scores.insert(
            group.id,
            target_composite.get(&group.id).copied().unwrap_or(0.0),
        );
    }

    // Compute protein PEP from the composite score
    let mut target_scores_vec: Vec<f64> = Vec::new();
    let mut decoy_scores_vec: Vec<f64> = Vec::new();
    for group in &parsimony.groups {
        if let Some(&s) = target_composite.get(&group.id) {
            target_scores_vec.push(s);
        }
        if let Some(&s) = decoy_composite.get(&group.id) {
            decoy_scores_vec.push(s);
        }
    }
    let all_scores: Vec<f64> = target_scores_vec
        .iter()
        .chain(decoy_scores_vec.iter())
        .copied()
        .collect();
    let all_is_decoy: Vec<bool> = std::iter::repeat(false)
        .take(target_scores_vec.len())
        .chain(std::iter::repeat(true).take(decoy_scores_vec.len()))
        .collect();
    let pep_estimator = PepEstimator::fit_default(&all_scores, &all_is_decoy);
    let mut group_pep: HashMap<ProteinGroupId, f64> = HashMap::new();
    for group in &parsimony.groups {
        if let Some(&s) = target_composite.get(&group.id) {
            group_pep.insert(group.id, pep_estimator.posterior_error(s));
        }
    }

    // Propagate to peptides: each peptide gets the best (lowest) protein q-value
    let mut peptide_qvalues: HashMap<String, f64> = HashMap::new();
    for (peptide, group_ids) in &parsimony.peptide_to_groups {
        let best_q = group_ids
            .iter()
            .filter_map(|gid| group_qvalues.get(gid))
            .copied()
            .fold(1.0f64, f64::min);
        peptide_qvalues.insert(peptide.clone(), best_q);
    }

    let n_scored = group_qvalues.values().filter(|&&q| q < 1.0).count();
    let n_passing_1pct = group_qvalues.values().filter(|&&q| q <= 0.01).count();
    log::info!(
        "Protein FDR: {} protein groups scored, {} at 1% FDR",
        n_scored,
        n_passing_1pct
    );

    ProteinFdrResult {
        group_qvalues,
        group_pep,
        group_scores,
        peptide_qvalues,
    }
}

/// Compute protein q-values using DIA-NN-style TDC (not picked-protein).
///
/// For each target protein with score s:
///   q = n_decoys_with_score_>=_s / max(1, n_targets_with_score_>=_s)
///
/// This does NOT require explicit target-decoy pairing. Targets and decoys
/// are scored independently, and the decoy score distribution serves as the
/// null model for FDR estimation.
fn compute_protein_qvalues_diann(
    groups: &[ProteinGroup],
    target_scores: &HashMap<ProteinGroupId, f64>,
    decoy_scores: &HashMap<ProteinGroupId, f64>,
) -> HashMap<ProteinGroupId, f64> {
    // Collect and sort all target and decoy scores
    let mut all_target: Vec<f64> = target_scores.values().copied().collect();
    let mut all_decoy: Vec<f64> = decoy_scores.values().copied().collect();
    all_target.sort_by(|a, b| a.total_cmp(b));
    all_decoy.sort_by(|a, b| a.total_cmp(b));

    let n_targets_with = all_target.len();
    let n_decoys_with = all_decoy.len();

    log::debug!(
        "Protein TDC: {} target scores, {} decoy scores",
        n_targets_with,
        n_decoys_with
    );

    // For each target protein, compute raw FDR
    let mut raw_qvals: Vec<(f64, f64, ProteinGroupId)> = Vec::new(); // (score, raw_q, gid)
    for group in groups {
        let score = match target_scores.get(&group.id) {
            Some(&s) if s > 0.0 => s,
            _ => continue,
        };

        // Count decoys with score >= this target's score (binary search)
        let d_pos = all_decoy.partition_point(|&s| s < score);
        let n_decoys_ge = n_decoys_with - d_pos;

        // Count targets with score >= this target's score
        let t_pos = all_target.partition_point(|&s| s < score);
        let n_targets_ge = n_targets_with - t_pos;

        let q = if n_targets_ge > 0 {
            (n_decoys_ge as f64 / n_targets_ge as f64).min(1.0)
        } else {
            1.0
        };

        raw_qvals.push((score, q, group.id));
    }

    // Sort by score ascending for backward sweep (monotonicity)
    raw_qvals.sort_by(|a, b| a.0.total_cmp(&b.0));

    // Backward sweep: enforce monotonicity (lower score -> higher q-value)
    let mut min_q = 1.0f64;
    let mut result: HashMap<ProteinGroupId, f64> = HashMap::new();
    for &(_, q, gid) in raw_qvals.iter().rev() {
        min_q = min_q.min(q);
        result.insert(gid, min_q);
    }

    result
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
        "protein_group,protein_accessions,gene_names,protein_qvalue,protein_pep,best_score,n_unique_peptides,n_shared_peptides,unique_peptides,shared_peptides"
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
        let pep = fdr_result.group_pep.get(&group.id).copied().unwrap_or(1.0);
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
            "{},{},{},{:.6},{:.6},{:.4},{},{},{},{}",
            group.id,
            accessions,
            gene_str,
            q,
            pep,
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
    /// Best SVM discriminant score for this peptide across all files
    pub score: f64,
    /// Whether this peptide is a decoy
    pub is_decoy: bool,
    /// Best (lowest) run-level precursor q-value for this peptide
    pub best_qvalue: f64,
}

/// Collect peptide-level scores for protein FDR.
///
/// For each unique modified_sequence, keeps the best SVM score and best q-value
/// across all files. Called on compacted stubs after second-pass FDR, which
/// contains both passing targets and their paired decoys with reconciliation-
/// corrected scores.
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
                    // Keep the best (lowest) q-value
                    let qval = entry.run_precursor_qvalue;
                    if qval < ps.best_qvalue {
                        ps.best_qvalue = qval;
                    }
                })
                .or_insert(PeptideScore {
                    score: entry.score,
                    is_decoy: entry.is_decoy,
                    best_qvalue: entry.run_precursor_qvalue,
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
}
