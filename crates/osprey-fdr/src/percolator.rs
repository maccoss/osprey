//! Native Percolator implementation for semi-supervised FDR control
//!
//! Implements the Percolator algorithm (Käll et al. 2007) as refined by
//! mokapot (Fondrie & Noble, 2021):
//! - 3-fold cross-validation with peptide-grouped fold assignment
//! - Iterative linear SVM training on high-confidence targets vs all decoys
//! - Grid search for SVM cost parameter C
//! - Per-run and experiment-level FDR with conservative (n_decoy+1)/n_target formula
//! - Posterior error probability via KDE + isotonic regression
//!
//! Osprey is peptide-centric: the units are precursors (peptide + charge)
//! and peptides (unique modified sequences), not PSMs.

use osprey_ml::matrix::Matrix;
use osprey_ml::pep::PepEstimator;
use osprey_ml::svm::{self, FeatureStandardizer, LinearSvm};
use rayon::prelude::*;
use std::collections::HashMap;

/// Configuration for the Percolator runner
#[derive(Debug, Clone)]
pub struct PercolatorConfig {
    /// FDR threshold for selecting positive training set (default: 0.01)
    pub train_fdr: f64,
    /// FDR threshold for final results (default: 0.01)
    pub test_fdr: f64,
    /// Maximum SVM training iterations per fold (default: 10)
    pub max_iterations: usize,
    /// Number of cross-validation folds (default: 3)
    pub n_folds: usize,
    /// Random seed for reproducibility (default: 42)
    pub seed: u64,
    /// Grid search C values for SVM cost parameter
    pub c_values: Vec<f64>,
}

impl Default for PercolatorConfig {
    fn default() -> Self {
        PercolatorConfig {
            train_fdr: 0.01,
            test_fdr: 0.01,
            max_iterations: 10,
            n_folds: 3,
            seed: 42,
            c_values: vec![0.01, 0.1, 1.0, 10.0],
        }
    }
}

/// Input entry for Percolator scoring
#[derive(Debug, Clone)]
pub struct PercolatorEntry {
    /// Unique precursor ID (e.g., "filename_libidx")
    pub id: String,
    /// Source file name (for per-run FDR)
    pub file_name: String,
    /// Modified peptide sequence (for fold grouping and peptide-level FDR)
    pub peptide: String,
    /// Precursor charge state
    pub charge: u8,
    /// Whether this is a decoy
    pub is_decoy: bool,
    /// Entry ID for target-decoy pairing (high bit = decoy)
    pub entry_id: u32,
    /// Raw feature values (37 or 45 depending on search mode)
    pub features: Vec<f64>,
}

/// Result for a single entry after Percolator scoring
#[derive(Debug, Clone)]
pub struct PercolatorResult {
    /// Unique precursor ID (matches PercolatorEntry.id)
    pub id: String,
    /// SVM decision function score
    pub score: f64,
    /// Per-run precursor-level q-value
    pub run_precursor_qvalue: f64,
    /// Per-run peptide-level q-value (best precursor per peptide)
    pub run_peptide_qvalue: f64,
    /// Experiment-wide precursor-level q-value
    pub experiment_precursor_qvalue: f64,
    /// Experiment-wide peptide-level q-value
    pub experiment_peptide_qvalue: f64,
    /// Posterior error probability
    pub pep: f64,
}

/// Full results from Percolator analysis
#[derive(Debug)]
pub struct PercolatorResults {
    /// Per-entry results
    pub entries: Vec<PercolatorResult>,
    /// Feature weights from best model per fold
    pub fold_weights: Vec<Vec<f64>>,
    /// Number of iterations used per fold
    pub iterations_per_fold: Vec<usize>,
}

/// Run the Percolator algorithm on a collection of entries
///
/// # Arguments
/// * `entries` - All precursor entries across all files (targets and decoys)
/// * `config` - Percolator configuration
///
/// # Returns
/// PercolatorResults with scores, q-values, and PEPs for all entries
pub fn run_percolator(
    entries: &[PercolatorEntry],
    config: &PercolatorConfig,
) -> Result<PercolatorResults, String> {
    if entries.is_empty() {
        return Ok(PercolatorResults {
            entries: Vec::new(),
            fold_weights: Vec::new(),
            iterations_per_fold: Vec::new(),
        });
    }

    let n = entries.len();
    let n_features = entries[0].features.len();
    let n_targets = entries.iter().filter(|e| !e.is_decoy).count();
    let n_decoys = entries.iter().filter(|e| e.is_decoy).count();

    log::info!(
        "Percolator: {} entries ({} targets, {} decoys), {} features, {}-fold CV",
        n,
        n_targets,
        n_decoys,
        n_features,
        config.n_folds
    );

    // 1. Build feature matrix
    let feature_data: Vec<f64> = entries
        .iter()
        .flat_map(|e| e.features.iter().copied())
        .collect();
    let features = Matrix::new(feature_data, n, n_features);
    let labels: Vec<bool> = entries.iter().map(|e| e.is_decoy).collect();
    let entry_ids: Vec<u32> = entries.iter().map(|e| e.entry_id).collect();
    let peptides: Vec<String> = entries.iter().map(|e| e.peptide.clone()).collect();

    // 2. Standardize features
    let (_standardizer, std_features) = FeatureStandardizer::fit_transform(&features);
    log::debug!("  Features standardized to zero mean, unit variance");

    // 3. Assign folds (group by peptide sequence to prevent leakage)
    let fold_assignments = create_stratified_folds_by_peptide(&labels, &peptides, config.n_folds);

    // 4. Find best initial feature
    let (best_feat_idx, best_feat_passing) =
        find_best_initial_feature(&std_features, &labels, &entry_ids, config.train_fdr);
    let initial_scores: Vec<f64> = (0..n).map(|i| std_features[(i, best_feat_idx)]).collect();
    log::info!(
        "  Best initial feature: index {} ({} targets at {:.0}% FDR)",
        best_feat_idx,
        best_feat_passing,
        config.train_fdr * 100.0
    );

    // 5. Train per-fold models via cross-validation (parallel)
    let mut final_scores = vec![0.0; n];
    let mut fold_weights: Vec<Vec<f64>> = Vec::new();
    let mut iterations_per_fold: Vec<usize> = Vec::new();

    // Train all folds in parallel — each fold is independent
    let fold_results: Vec<(usize, LinearSvm, usize)> = (0..config.n_folds)
        .into_par_iter()
        .map(|fold| {
            let train_indices: Vec<usize> =
                (0..n).filter(|&i| fold_assignments[i] != fold).collect();
            let test_indices: Vec<usize> =
                (0..n).filter(|&i| fold_assignments[i] == fold).collect();

            log::info!(
                "  Fold {}/{}: {} train, {} test",
                fold + 1,
                config.n_folds,
                train_indices.len(),
                test_indices.len()
            );

            let (best_model, n_iterations) = train_fold(
                &std_features,
                &labels,
                &entry_ids,
                &peptides,
                &train_indices,
                &initial_scores,
                config,
            );

            (fold, best_model, n_iterations)
        })
        .collect();

    // Merge results sequentially (deterministic order)
    for (fold, best_model, n_iterations) in &fold_results {
        let test_indices: Vec<usize> = (0..n).filter(|&i| fold_assignments[i] == *fold).collect();
        let test_features = extract_rows(&std_features, &test_indices);
        let test_scores = best_model.decision_function(&test_features);
        for (i, &idx) in test_indices.iter().enumerate() {
            final_scores[idx] = test_scores[i];
        }

        fold_weights.push(best_model.weights().to_vec());
        iterations_per_fold.push(*n_iterations);
        log::info!(
            "    Fold {} best model weights: [{}]",
            fold + 1,
            best_model
                .weights()
                .iter()
                .map(|w| format!("{:.4}", w))
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    // 5b. Calibrate scores between folds (Granholm et al. 2012)
    calibrate_scores_between_folds(
        &mut final_scores,
        &fold_assignments,
        &labels,
        &entry_ids,
        config.n_folds,
        config.train_fdr,
    );

    // 6. Compute PEP on competition winners
    let (winner_indices, winner_scores, winner_is_decoy) =
        compete_all(&final_scores, &labels, &entry_ids);

    let pep_estimator = PepEstimator::fit_default(&winner_scores, &winner_is_decoy);
    let mut peps = vec![1.0; n];
    for &idx in &winner_indices {
        peps[idx] = pep_estimator.posterior_error(final_scores[idx]);
    }

    // 7. Compute q-values at precursor and peptide levels
    let file_names: Vec<String> = entries.iter().map(|e| e.file_name.clone()).collect();
    let unique_files: Vec<String> = {
        let mut files: Vec<String> = file_names.to_vec();
        files.sort();
        files.dedup();
        files
    };
    let is_single_file = unique_files.len() <= 1;

    // Per-run precursor q-values
    let run_precursor_qvalues =
        compute_per_run_precursor_qvalues(&final_scores, &labels, &entry_ids, &file_names);

    // Per-run peptide q-values
    let run_peptide_qvalues =
        compute_per_run_peptide_qvalues(&final_scores, &labels, &entry_ids, &file_names, &peptides);

    // Experiment-level q-values
    let (exp_precursor_qvalues, exp_peptide_qvalues) = if is_single_file {
        // Single file: copy run-level
        (run_precursor_qvalues.clone(), run_peptide_qvalues.clone())
    } else {
        let exp_prec = compute_experiment_precursor_qvalues(&final_scores, &labels, &entry_ids);
        let exp_pept =
            compute_experiment_peptide_qvalues(&final_scores, &labels, &entry_ids, &peptides);
        (exp_prec, exp_pept)
    };

    // 8. Build results
    let results: Vec<PercolatorResult> = entries
        .iter()
        .enumerate()
        .map(|(i, entry)| PercolatorResult {
            id: entry.id.clone(),
            score: final_scores[i],
            run_precursor_qvalue: run_precursor_qvalues[i],
            run_peptide_qvalue: run_peptide_qvalues[i],
            experiment_precursor_qvalue: exp_precursor_qvalues[i],
            experiment_peptide_qvalue: exp_peptide_qvalues[i],
            pep: peps[i],
        })
        .collect();

    // Log summary statistics
    let n_prec_pass = results
        .iter()
        .zip(entries)
        .filter(|(r, e)| !e.is_decoy && r.experiment_precursor_qvalue <= config.test_fdr)
        .count();
    let n_pept_pass = results
        .iter()
        .zip(entries)
        .filter(|(r, e)| !e.is_decoy && r.experiment_peptide_qvalue <= config.test_fdr)
        .count();
    log::info!(
        "Percolator results: {} precursors, {} peptides at {:.0}% FDR",
        n_prec_pass,
        n_pept_pass,
        config.test_fdr * 100.0
    );

    Ok(PercolatorResults {
        entries: results,
        fold_weights,
        iterations_per_fold,
    })
}

/// Train SVM iteratively for one CV fold
#[allow(clippy::too_many_arguments)]
fn train_fold(
    std_features: &Matrix,
    labels: &[bool],
    entry_ids: &[u32],
    peptides: &[String],
    train_indices: &[usize],
    initial_scores: &[f64],
    config: &PercolatorConfig,
) -> (LinearSvm, usize) {
    let n_features = std_features.cols;
    let mut current_scores = initial_scores.to_vec();

    // Track best model across iterations
    let mut best_model = LinearSvm::fit(&Matrix::zeros(0, n_features), &[], 1.0, config.seed);
    let mut best_iteration = 0usize;
    let mut consecutive_no_improve = 0usize;

    // Start with no trained model (zero weights), so any SVM that passes targets will be better
    let train_labels: Vec<bool> = train_indices.iter().map(|&i| labels[i]).collect();
    let train_entry_ids: Vec<u32> = train_indices.iter().map(|&i| entry_ids[i]).collect();
    let mut best_passing = 0usize;

    for iteration in 0..config.max_iterations {
        // i. Select positive training set
        let train_current_scores: Vec<f64> =
            train_indices.iter().map(|&i| current_scores[i]).collect();
        let selected_target_indices = select_positive_training_set(
            &train_current_scores,
            &train_labels,
            &train_entry_ids,
            config.train_fdr,
            50, // min_positive
        );

        if selected_target_indices.is_empty() {
            log::warn!(
                "    Iteration {}: no targets selected, stopping",
                iteration + 1
            );
            break;
        }

        // Collect all decoy indices in training set
        let decoy_indices: Vec<usize> = (0..train_indices.len())
            .filter(|&i| train_labels[i])
            .collect();

        // Build SVM training set: selected targets + all decoys
        let mut svm_indices: Vec<usize> = selected_target_indices.clone();
        svm_indices.extend_from_slice(&decoy_indices);

        // Map to global indices for feature extraction
        let svm_global_indices: Vec<usize> =
            svm_indices.iter().map(|&i| train_indices[i]).collect();
        let svm_features = extract_rows(std_features, &svm_global_indices);
        let svm_labels: Vec<bool> = svm_indices.iter().map(|&i| train_labels[i]).collect();
        let svm_entry_ids: Vec<u32> = svm_indices.iter().map(|&i| train_entry_ids[i]).collect();

        // ii. Grid search for best C
        let svm_fold_assignments = create_stratified_folds_by_peptide(
            &svm_labels,
            &svm_indices
                .iter()
                .map(|&i| peptides[train_indices[i]].clone())
                .collect::<Vec<_>>(),
            config.n_folds,
        );

        let best_c = svm::grid_search_c(
            &svm_features,
            &svm_labels,
            &svm_entry_ids,
            &config.c_values,
            &svm_fold_assignments,
            config.n_folds,
            config.seed,
            config.train_fdr,
        );

        // iii. Train SVM with best C
        let model = LinearSvm::fit(&svm_features, &svm_labels, best_c, config.seed);

        // iv. Score ALL training set entries with new model
        let train_features = extract_rows(std_features, train_indices);
        let new_train_scores = model.decision_function(&train_features);

        // Update current scores for training set
        for (i, &idx) in train_indices.iter().enumerate() {
            current_scores[idx] = new_train_scores[i];
        }

        // v. Count passing targets (non-conservative for iteration tracking)
        let n_passing = count_passing(
            &new_train_scores,
            &train_labels,
            &train_entry_ids,
            config.train_fdr,
        );

        log::info!(
            "    Iteration {}: C={:.4}, {} selected targets, {} pass {:.0}% FDR",
            iteration + 1,
            best_c,
            selected_target_indices.len(),
            n_passing,
            config.train_fdr * 100.0
        );

        if n_passing > best_passing {
            best_model = model;
            best_passing = n_passing;
            best_iteration = iteration + 1;
            consecutive_no_improve = 0;
            log::info!("      -> New best: {} passing", best_passing);
        } else {
            consecutive_no_improve += 1;
            log::info!(
                "      -> No improvement (best={} from iter {}, {} consecutive)",
                best_passing,
                best_iteration,
                consecutive_no_improve
            );
        }

        if consecutive_no_improve >= 2 {
            log::info!("    Stopping early: 2 consecutive non-improvements");
            break;
        }
    }

    let final_iteration = best_iteration.max(1);
    log::info!(
        "    Fold complete: best iteration {}, {} passing",
        best_iteration,
        best_passing
    );

    (best_model, final_iteration)
}

/// Select high-confidence targets for positive training set
fn select_positive_training_set(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    fdr_threshold: f64,
    min_targets: usize,
) -> Vec<usize> {
    // Target-decoy competition within the provided data
    let (winner_local_indices, winner_scores, winner_is_decoy) =
        compete_subset(scores, labels, entry_ids);

    // Compute q-values on winners (non-conservative for training selection)
    let mut q_values = vec![1.0; winner_local_indices.len()];
    if !winner_local_indices.is_empty() {
        compute_qvalues(&winner_scores, &winner_is_decoy, &mut q_values);
    }

    // Select target winners passing threshold
    let select_at_threshold = |threshold: f64| -> Vec<usize> {
        winner_local_indices
            .iter()
            .enumerate()
            .filter(|&(rank, &idx)| !labels[idx] && q_values[rank] <= threshold)
            .map(|(_, &idx)| idx)
            .collect()
    };

    let mut selected = select_at_threshold(fdr_threshold);

    // Relax threshold if too few
    if selected.len() < min_targets {
        for &threshold in &[0.05, 0.10, 0.25, 0.50] {
            selected = select_at_threshold(threshold);
            if selected.len() >= min_targets {
                log::debug!(
                    "      Relaxed FDR to {:.0}% to get {} targets",
                    threshold * 100.0,
                    selected.len()
                );
                break;
            }
        }
    }

    selected
}

/// Find the best single feature for initial scoring
fn find_best_initial_feature(
    features: &Matrix,
    labels: &[bool],
    entry_ids: &[u32],
    fdr_threshold: f64,
) -> (usize, usize) {
    let n = features.rows;
    let p = features.cols;
    let mut best_idx = 0;
    let mut best_passing = 0;

    for feat in 0..p {
        let scores: Vec<f64> = (0..n).map(|i| features[(i, feat)]).collect();
        let n_pass = count_passing(&scores, labels, entry_ids, fdr_threshold);
        if n_pass > best_passing {
            best_passing = n_pass;
            best_idx = feat;
        }
    }

    (best_idx, best_passing)
}

// ============================================================
// Target-decoy competition and q-value computation
// ============================================================

/// Target-decoy competition on all entries, returning winners sorted by score descending
fn compete_all(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
) -> (Vec<usize>, Vec<f64>, Vec<bool>) {
    let all_indices: Vec<usize> = (0..scores.len()).collect();
    compete_from_indices(scores, labels, entry_ids, &all_indices)
}

/// Target-decoy competition on a subset, returning winners
fn compete_subset(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
) -> (Vec<usize>, Vec<f64>, Vec<bool>) {
    let all_indices: Vec<usize> = (0..scores.len()).collect();
    compete_from_indices(scores, labels, entry_ids, &all_indices)
}

/// Core competition logic: group by base_id, compete, return winners sorted by score desc
fn compete_from_indices(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    indices: &[usize],
) -> (Vec<usize>, Vec<f64>, Vec<bool>) {
    // Group by base_id = entry_id & 0x7FFFFFFF
    let mut targets: HashMap<u32, (usize, f64)> = HashMap::new();
    let mut decoys: HashMap<u32, (usize, f64)> = HashMap::new();

    for &idx in indices {
        let base_id = entry_ids[idx] & 0x7FFFFFFF;
        if labels[idx] {
            decoys
                .entry(base_id)
                .and_modify(|(best_idx, best_score)| {
                    if scores[idx] > *best_score {
                        *best_idx = idx;
                        *best_score = scores[idx];
                    }
                })
                .or_insert((idx, scores[idx]));
        } else {
            targets
                .entry(base_id)
                .and_modify(|(best_idx, best_score)| {
                    if scores[idx] > *best_score {
                        *best_idx = idx;
                        *best_score = scores[idx];
                    }
                })
                .or_insert((idx, scores[idx]));
        }
    }

    // Compete pairs: higher score wins, ties → decoy
    let mut winners: Vec<(usize, f64, bool)> = Vec::with_capacity(targets.len());
    for (base_id, &(t_idx, t_score)) in &targets {
        if let Some(&(d_idx, d_score)) = decoys.get(base_id) {
            if t_score > d_score {
                winners.push((t_idx, t_score, false));
            } else {
                winners.push((d_idx, d_score, true)); // tie → decoy
            }
        } else {
            winners.push((t_idx, t_score, false)); // unpaired target wins
        }
    }
    // Unpaired decoys
    for (base_id, &(d_idx, d_score)) in &decoys {
        if !targets.contains_key(base_id) {
            winners.push((d_idx, d_score, true));
        }
    }

    // Sort by score descending
    winners.sort_by(|a, b| b.1.total_cmp(&a.1));

    let winner_indices: Vec<usize> = winners.iter().map(|w| w.0).collect();
    let winner_scores: Vec<f64> = winners.iter().map(|w| w.1).collect();
    let winner_is_decoy: Vec<bool> = winners.iter().map(|w| w.2).collect();

    (winner_indices, winner_scores, winner_is_decoy)
}

/// Compute conservative q-values: FDR = (n_decoy + 1) / n_target
///
/// Input must be sorted by score descending (winners from competition).
fn compute_conservative_qvalues(_scores: &[f64], is_decoy: &[bool], q_values: &mut [f64]) {
    assert_eq!(is_decoy.len(), q_values.len());
    let n = is_decoy.len();

    let mut n_target = 0usize;
    let mut n_decoy = 0usize;

    // Forward pass: compute FDR at each position
    for i in 0..n {
        if is_decoy[i] {
            n_decoy += 1;
        } else {
            n_target += 1;
        }
        q_values[i] = if n_target > 0 {
            (n_decoy + 1) as f64 / n_target as f64
        } else {
            1.0
        };
    }

    // Backward pass: make monotonically non-increasing (q-value property)
    let mut q_min = 1.0f64;
    for q in q_values.iter_mut().rev() {
        q_min = q_min.min(*q);
        *q = q_min;
    }
}

/// Count targets passing FDR threshold using conservative formula
#[allow(dead_code)]
fn count_passing_conservative(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    fdr_threshold: f64,
) -> usize {
    let (winner_indices, _, winner_is_decoy) = compete_from_indices(
        scores,
        labels,
        entry_ids,
        &(0..scores.len()).collect::<Vec<_>>(),
    );

    let mut q_values = vec![1.0; winner_indices.len()];
    let winner_scores: Vec<f64> = winner_indices.iter().map(|&i| scores[i]).collect();
    compute_conservative_qvalues(&winner_scores, &winner_is_decoy, &mut q_values);

    winner_indices
        .iter()
        .enumerate()
        .filter(|&(rank, &idx)| !labels[idx] && q_values[rank] <= fdr_threshold)
        .count()
}

/// Compute non-conservative q-values: FDR = n_decoy / n_target
///
/// Used internally for iteration tracking and positive training set selection.
/// Less biased for small sample sizes than conservative formula.
/// Input must be sorted by score descending.
fn compute_qvalues(_scores: &[f64], is_decoy: &[bool], q_values: &mut [f64]) {
    assert_eq!(is_decoy.len(), q_values.len());
    let n = is_decoy.len();

    let mut n_target = 0usize;
    let mut n_decoy = 0usize;

    for i in 0..n {
        if is_decoy[i] {
            n_decoy += 1;
        } else {
            n_target += 1;
        }
        q_values[i] = if n_target > 0 {
            n_decoy as f64 / n_target as f64
        } else {
            1.0
        };
    }

    // Backward pass: monotonically non-increasing
    let mut q_min = 1.0f64;
    for q in q_values.iter_mut().rev() {
        q_min = q_min.min(*q);
        *q = q_min;
    }
}

/// Count targets passing FDR threshold using non-conservative formula
///
/// Used for internal iteration tracking and best-feature selection during training.
/// The non-conservative n_decoy/n_target formula avoids the +1 bias that makes
/// it impossible to reach low FDR thresholds with small sample sizes.
fn count_passing(scores: &[f64], labels: &[bool], entry_ids: &[u32], fdr_threshold: f64) -> usize {
    let (winner_indices, _, winner_is_decoy) = compete_from_indices(
        scores,
        labels,
        entry_ids,
        &(0..scores.len()).collect::<Vec<_>>(),
    );

    let mut q_values = vec![1.0; winner_indices.len()];
    let winner_scores: Vec<f64> = winner_indices.iter().map(|&i| scores[i]).collect();
    compute_qvalues(&winner_scores, &winner_is_decoy, &mut q_values);

    winner_indices
        .iter()
        .enumerate()
        .filter(|&(rank, &idx)| !labels[idx] && q_values[rank] <= fdr_threshold)
        .count()
}

// ============================================================
// Score calibration between CV folds (Granholm et al. 2012)
// ============================================================

/// Calibrate scores between CV folds using the Granholm et al. (2012) method.
///
/// Each fold's SVM produces scores on a different scale. This linear transform
/// normalizes them so that:
/// - The score at the FDR threshold maps to 0
/// - The median decoy score maps to -1
///
/// This is the same method used by C++ Percolator (`normalizeScores`) and
/// mokapot (`calibrate_scores`).
///
/// Reference: Granholm V, Noble WS, Käll L. BMC Bioinformatics. 2012;13 Suppl 16:S3.
fn calibrate_scores_between_folds(
    final_scores: &mut [f64],
    fold_assignments: &[usize],
    labels: &[bool],
    entry_ids: &[u32],
    n_folds: usize,
    fdr_threshold: f64,
) {
    for fold in 0..n_folds {
        let test_indices: Vec<usize> = (0..final_scores.len())
            .filter(|&i| fold_assignments[i] == fold)
            .collect();

        if test_indices.is_empty() {
            continue;
        }

        // Find the score at the FDR threshold via target-decoy competition
        let test_scores: Vec<f64> = test_indices.iter().map(|&i| final_scores[i]).collect();
        let test_labels: Vec<bool> = test_indices.iter().map(|&i| labels[i]).collect();
        let test_eids: Vec<u32> = test_indices.iter().map(|&i| entry_ids[i]).collect();

        let threshold_score =
            match find_score_at_fdr(&test_scores, &test_labels, &test_eids, fdr_threshold) {
                Some(s) => s,
                None => {
                    log::warn!(
                        "  Fold {}: no targets pass {:.0}% FDR, skipping calibration",
                        fold + 1,
                        fdr_threshold * 100.0
                    );
                    continue;
                }
            };

        // Find median decoy score in this fold
        let mut decoy_scores: Vec<f64> = test_indices
            .iter()
            .filter(|&&i| labels[i])
            .map(|&i| final_scores[i])
            .collect();
        decoy_scores.sort_by(|a, b| a.total_cmp(b));

        let median_decoy = if decoy_scores.is_empty() {
            threshold_score - 1.0 // fallback
        } else {
            decoy_scores[decoy_scores.len() / 2]
        };

        // Linear transform: threshold_score → 0, median_decoy → -1
        let denom = threshold_score - median_decoy;
        if denom <= 0.0 {
            log::warn!(
                "  Fold {}: median decoy score ({:.4}) >= FDR threshold score ({:.4}), skipping calibration",
                fold + 1,
                median_decoy,
                threshold_score
            );
            continue;
        }

        for &idx in &test_indices {
            final_scores[idx] = (final_scores[idx] - threshold_score) / denom;
        }

        log::debug!(
            "  Fold {}: calibrated scores (FDR threshold={:.4}, median decoy={:.4}, scale={:.4})",
            fold + 1,
            threshold_score,
            median_decoy,
            denom
        );
    }
}

/// Find the minimum target score that passes the FDR threshold after competition.
///
/// Returns None if no targets pass the threshold.
fn find_score_at_fdr(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    fdr_threshold: f64,
) -> Option<f64> {
    let (winner_indices, _, winner_is_decoy) = compete_from_indices(
        scores,
        labels,
        entry_ids,
        &(0..scores.len()).collect::<Vec<_>>(),
    );

    if winner_indices.is_empty() {
        return None;
    }

    // Compute q-values on winners (non-conservative for calibration)
    let mut q_values = vec![1.0; winner_indices.len()];
    let winner_scores: Vec<f64> = winner_indices.iter().map(|&i| scores[i]).collect();
    compute_qvalues(&winner_scores, &winner_is_decoy, &mut q_values);

    // Find the minimum score among targets passing the FDR threshold
    let mut min_passing_score: Option<f64> = None;
    for (rank, &idx) in winner_indices.iter().enumerate() {
        if !labels[idx] && q_values[rank] <= fdr_threshold {
            match min_passing_score {
                None => min_passing_score = Some(scores[idx]),
                Some(current_min) => {
                    if scores[idx] < current_min {
                        min_passing_score = Some(scores[idx]);
                    }
                }
            }
        }
    }

    min_passing_score
}

// ============================================================
// Per-run and experiment-level q-value computation
// ============================================================

/// Compute per-run precursor-level q-values
fn compute_per_run_precursor_qvalues(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    file_names: &[String],
) -> Vec<f64> {
    let n = scores.len();
    let mut qvalues = vec![1.0; n];

    // Group entries by file
    let mut file_groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, fname) in file_names.iter().enumerate() {
        file_groups.entry(fname.as_str()).or_default().push(i);
    }

    for indices in file_groups.values() {
        let file_scores: Vec<f64> = indices.iter().map(|&i| scores[i]).collect();
        let file_labels: Vec<bool> = indices.iter().map(|&i| labels[i]).collect();
        let file_entry_ids: Vec<u32> = indices.iter().map(|&i| entry_ids[i]).collect();

        let (winner_local, _, winner_is_decoy) = compete_from_indices(
            &file_scores,
            &file_labels,
            &file_entry_ids,
            &(0..indices.len()).collect::<Vec<_>>(),
        );

        let mut q = vec![1.0; winner_local.len()];
        let ws: Vec<f64> = winner_local.iter().map(|&i| file_scores[i]).collect();
        compute_conservative_qvalues(&ws, &winner_is_decoy, &mut q);

        // Map back: winners get their q-value
        for (rank, &local_idx) in winner_local.iter().enumerate() {
            let global_idx = indices[local_idx];
            qvalues[global_idx] = q[rank];
        }
    }

    qvalues
}

/// Compute per-run peptide-level q-values
///
/// Keep best-scoring precursor per peptide within each file, then compete
fn compute_per_run_peptide_qvalues(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    file_names: &[String],
    peptides: &[String],
) -> Vec<f64> {
    let n = scores.len();
    let mut qvalues = vec![1.0; n];

    let mut file_groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, fname) in file_names.iter().enumerate() {
        file_groups.entry(fname.as_str()).or_default().push(i);
    }

    for indices in file_groups.values() {
        // Find best precursor per peptide
        let best_per_peptide = best_precursor_per_peptide(indices, scores, labels, peptides);

        // Competition on peptide representatives
        let pept_scores: Vec<f64> = best_per_peptide.iter().map(|&i| scores[i]).collect();
        let pept_labels: Vec<bool> = best_per_peptide.iter().map(|&i| labels[i]).collect();
        let pept_entry_ids: Vec<u32> = best_per_peptide.iter().map(|&i| entry_ids[i]).collect();

        let (winner_local, _, winner_is_decoy) = compete_from_indices(
            &pept_scores,
            &pept_labels,
            &pept_entry_ids,
            &(0..best_per_peptide.len()).collect::<Vec<_>>(),
        );

        let mut q = vec![1.0; winner_local.len()];
        let ws: Vec<f64> = winner_local.iter().map(|&i| pept_scores[i]).collect();
        compute_conservative_qvalues(&ws, &winner_is_decoy, &mut q);

        // Map q-values back to all precursors of each peptide
        // Build peptide → q-value map from winners
        let mut peptide_qvalue: HashMap<&str, f64> = HashMap::new();
        for (rank, &local_idx) in winner_local.iter().enumerate() {
            let global_idx = best_per_peptide[local_idx];
            let pept = &peptides[global_idx];
            peptide_qvalue.insert(pept.as_str(), q[rank]);
        }

        // Assign to all precursors of that peptide in this file
        for &idx in indices {
            if let Some(&q) = peptide_qvalue.get(peptides[idx].as_str()) {
                qvalues[idx] = q;
            }
        }
    }

    qvalues
}

/// Compute experiment-level precursor q-values (across all files)
fn compute_experiment_precursor_qvalues(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
) -> Vec<f64> {
    let n = scores.len();
    let mut qvalues = vec![1.0; n];

    let (winner_indices, _, winner_is_decoy) = compete_all(scores, labels, entry_ids);

    let mut q = vec![1.0; winner_indices.len()];
    let ws: Vec<f64> = winner_indices.iter().map(|&i| scores[i]).collect();
    compute_conservative_qvalues(&ws, &winner_is_decoy, &mut q);

    for (rank, &idx) in winner_indices.iter().enumerate() {
        qvalues[idx] = q[rank];
    }

    qvalues
}

/// Compute experiment-level peptide q-values
fn compute_experiment_peptide_qvalues(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    peptides: &[String],
) -> Vec<f64> {
    let n = scores.len();
    let mut qvalues = vec![1.0; n];

    let all_indices: Vec<usize> = (0..n).collect();
    let best_per_peptide = best_precursor_per_peptide(&all_indices, scores, labels, peptides);

    let pept_scores: Vec<f64> = best_per_peptide.iter().map(|&i| scores[i]).collect();
    let pept_labels: Vec<bool> = best_per_peptide.iter().map(|&i| labels[i]).collect();
    let pept_entry_ids: Vec<u32> = best_per_peptide.iter().map(|&i| entry_ids[i]).collect();

    let (winner_local, _, winner_is_decoy) = compete_from_indices(
        &pept_scores,
        &pept_labels,
        &pept_entry_ids,
        &(0..best_per_peptide.len()).collect::<Vec<_>>(),
    );

    let mut q = vec![1.0; winner_local.len()];
    let ws: Vec<f64> = winner_local.iter().map(|&i| pept_scores[i]).collect();
    compute_conservative_qvalues(&ws, &winner_is_decoy, &mut q);

    let mut peptide_qvalue: HashMap<&str, f64> = HashMap::new();
    for (rank, &local_idx) in winner_local.iter().enumerate() {
        let global_idx = best_per_peptide[local_idx];
        peptide_qvalue.insert(peptides[global_idx].as_str(), q[rank]);
    }

    for i in 0..n {
        if let Some(&qv) = peptide_qvalue.get(peptides[i].as_str()) {
            qvalues[i] = qv;
        }
    }

    qvalues
}

/// Find the best-scoring precursor per peptide from a set of indices
fn best_precursor_per_peptide(
    indices: &[usize],
    scores: &[f64],
    _labels: &[bool],
    peptides: &[String],
) -> Vec<usize> {
    let mut best: HashMap<&str, (usize, f64)> = HashMap::new();

    for &idx in indices {
        let pept = peptides[idx].as_str();
        best.entry(pept)
            .and_modify(|(best_idx, best_score)| {
                if scores[idx] > *best_score {
                    *best_idx = idx;
                    *best_score = scores[idx];
                }
            })
            .or_insert((idx, scores[idx]));
    }

    let mut result: Vec<usize> = best.values().map(|&(idx, _)| idx).collect();
    result.sort(); // deterministic order
    result
}

// ============================================================
// Fold assignment
// ============================================================

/// Create stratified fold assignments grouped by peptide sequence
///
/// All charge states of the same peptide go into the same fold to prevent leakage.
/// Targets and decoys are assigned independently for class balance.
fn create_stratified_folds_by_peptide(
    labels: &[bool],
    peptides: &[String],
    n_folds: usize,
) -> Vec<usize> {
    let mut fold_assignments = vec![0; labels.len()];

    // Group indices by peptide sequence and target/decoy
    let mut target_peptides: HashMap<&str, Vec<usize>> = HashMap::new();
    let mut decoy_peptides: HashMap<&str, Vec<usize>> = HashMap::new();

    for (i, (pept, &is_decoy)) in peptides.iter().zip(labels).enumerate() {
        if is_decoy {
            decoy_peptides.entry(pept.as_str()).or_default().push(i);
        } else {
            target_peptides.entry(pept.as_str()).or_default().push(i);
        }
    }

    // Sort for deterministic assignment
    let mut target_groups: Vec<(&str, &Vec<usize>)> =
        target_peptides.iter().map(|(&k, v)| (k, v)).collect();
    let mut decoy_groups: Vec<(&str, &Vec<usize>)> =
        decoy_peptides.iter().map(|(&k, v)| (k, v)).collect();
    target_groups.sort_by_key(|&(k, _)| k);
    decoy_groups.sort_by_key(|&(k, _)| k);

    for (i, (_, indices)) in target_groups.iter().enumerate() {
        let fold = i % n_folds;
        for &idx in *indices {
            fold_assignments[idx] = fold;
        }
    }

    for (i, (_, indices)) in decoy_groups.iter().enumerate() {
        let fold = i % n_folds;
        for &idx in *indices {
            fold_assignments[idx] = fold;
        }
    }

    fold_assignments
}

/// Extract specific rows from a matrix
fn extract_rows(matrix: &Matrix, row_indices: &[usize]) -> Matrix {
    let n_cols = matrix.cols;
    let data: Vec<f64> = row_indices
        .iter()
        .flat_map(|&row| matrix.row_slice(row).iter().copied())
        .collect();
    Matrix::new(data, row_indices.len(), n_cols)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry(
        id: &str,
        file: &str,
        peptide: &str,
        charge: u8,
        is_decoy: bool,
        entry_id: u32,
        features: Vec<f64>,
    ) -> PercolatorEntry {
        PercolatorEntry {
            id: id.to_string(),
            file_name: file.to_string(),
            peptide: peptide.to_string(),
            charge,
            is_decoy,
            entry_id,
            features,
        }
    }

    #[test]
    fn test_fold_assignment_peptide_grouping() {
        let labels = vec![false, false, false, true, true, true];
        let peptides = vec![
            "PEPTIDEK".to_string(),
            "PEPTIDEK".to_string(), // same peptide, different charge
            "ANOTHERONE".to_string(),
            "PEPTIDEK".to_string(), // decoy
            "PEPTIDEK".to_string(), // decoy, same peptide
            "ANOTHERONE".to_string(),
        ];

        let folds = create_stratified_folds_by_peptide(&labels, &peptides, 3);

        // All charge states of PEPTIDEK (targets) should be in same fold
        assert_eq!(folds[0], folds[1], "Same peptide targets should share fold");

        // All charge states of PEPTIDEK (decoys) should be in same fold
        assert_eq!(folds[3], folds[4], "Same peptide decoys should share fold");
    }

    #[test]
    fn test_conservative_qvalues() {
        // 3 targets, 1 decoy (already sorted by score desc)
        let scores = vec![10.0, 9.0, 8.0, 7.0];
        let is_decoy = vec![false, false, true, false];
        let mut q = vec![0.0; 4];

        compute_conservative_qvalues(&scores, &is_decoy, &mut q);

        // FDR walk with +1:
        // pos 0: 1T, 0D → (0+1)/1 = 1.0
        // pos 1: 2T, 0D → (0+1)/2 = 0.5
        // pos 2: 2T, 1D → (1+1)/2 = 1.0
        // pos 3: 3T, 1D → (1+1)/3 = 0.667
        // Backward pass (min): [0.5, 0.5, 0.667, 0.667]
        assert!((q[0] - 0.5).abs() < 1e-10, "q[0]={}", q[0]);
        assert!((q[1] - 0.5).abs() < 1e-10, "q[1]={}", q[1]);
        assert!((q[2] - 2.0 / 3.0).abs() < 1e-10, "q[2]={}", q[2]);
        assert!((q[3] - 2.0 / 3.0).abs() < 1e-10, "q[3]={}", q[3]);
    }

    #[test]
    fn test_compete_and_count() {
        // 5 targets easily beat 5 decoys
        let scores = vec![10.0, 9.0, 8.0, 7.0, 6.0, 1.0, 0.5, 0.2, 0.1, 0.05];
        let labels = vec![
            false, false, false, false, false, true, true, true, true, true,
        ];
        let entry_ids: Vec<u32> = vec![
            1,
            2,
            3,
            4,
            5,
            1 | 0x80000000,
            2 | 0x80000000,
            3 | 0x80000000,
            4 | 0x80000000,
            5 | 0x80000000,
        ];

        let n = count_passing_conservative(&scores, &labels, &entry_ids, 0.50);
        // 5 targets win, 0 decoys → FDR=(0+1)/n
        // All pass at 50% when n≥2
        assert_eq!(n, 5);
    }

    #[test]
    fn test_percolator_basic() {
        // Build synthetic well-separated data
        let mut entries = Vec::new();

        // 20 targets with high feature values
        for i in 0..20 {
            entries.push(make_entry(
                &format!("file1_{}", i),
                "file1",
                &format!("PEPTIDE{}", i),
                2,
                false,
                (i + 1) as u32,
                vec![4.0 + (i as f64) * 0.1, 5.0 - (i as f64) * 0.05],
            ));
        }

        // 20 paired decoys with low feature values
        for i in 0..20 {
            entries.push(make_entry(
                &format!("file1_d{}", i),
                "file1",
                &format!("DECOY{}", i),
                2,
                true,
                ((i + 1) as u32) | 0x80000000,
                vec![0.5 + (i as f64) * 0.05, 1.0 + (i as f64) * 0.02],
            ));
        }

        let config = PercolatorConfig {
            max_iterations: 3,
            ..Default::default()
        };

        let results = run_percolator(&entries, &config).unwrap();

        // Should have results for all entries
        assert_eq!(results.entries.len(), 40);

        // Targets should generally have higher scores than decoys
        let target_scores: Vec<f64> = results
            .entries
            .iter()
            .zip(&entries)
            .filter(|(_, e)| !e.is_decoy)
            .map(|(r, _)| r.score)
            .collect();
        let decoy_scores: Vec<f64> = results
            .entries
            .iter()
            .zip(&entries)
            .filter(|(_, e)| e.is_decoy)
            .map(|(r, _)| r.score)
            .collect();

        let avg_target = target_scores.iter().sum::<f64>() / target_scores.len() as f64;
        let avg_decoy = decoy_scores.iter().sum::<f64>() / decoy_scores.len() as f64;
        assert!(
            avg_target > avg_decoy,
            "avg_target={} should be > avg_decoy={}",
            avg_target,
            avg_decoy
        );

        // Should have fold weights
        assert_eq!(results.fold_weights.len(), 3);
    }

    #[test]
    fn test_best_precursor_per_peptide() {
        let indices = vec![0, 1, 2, 3, 4];
        let scores = vec![5.0, 3.0, 7.0, 2.0, 6.0];
        let labels = vec![false, false, false, true, true];
        let peptides = vec![
            "PEPK".to_string(),
            "PEPK".to_string(), // same peptide, lower score
            "OTHER".to_string(),
            "PEPK".to_string(),  // decoy
            "OTHER".to_string(), // decoy
        ];

        let best = best_precursor_per_peptide(&indices, &scores, &labels, &peptides);

        // Should have one entry per unique peptide: PEPK (best=0, score=5.0) and OTHER (best=4, score=6.0)
        // Note: decoy PEPK has score 2.0 but target PEPK has score 5.0 - both are in the map
        // since we group by peptide string (not by target/decoy)
        assert_eq!(best.len(), 2);
        // Best PEPK = index 0 (score 5.0, but wait index 2 is OTHER with score 7.0)
        // Actually: PEPK has indices 0(5.0), 1(3.0), 3(2.0) → best=0
        // OTHER has indices 2(7.0), 4(6.0) → best=2
        assert!(best.contains(&0));
        assert!(best.contains(&2));
    }

    #[test]
    fn test_percolator_empty() {
        let config = PercolatorConfig::default();
        let results = run_percolator(&[], &config).unwrap();
        assert!(results.entries.is_empty());
    }
}
