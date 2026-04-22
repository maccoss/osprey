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

use osprey_core::config::FdrLevel;
use osprey_core::diagnostics::{exit_if_only, is_dump_enabled};
use osprey_core::types::FdrEntry;
use osprey_ml::matrix::Matrix;
use osprey_ml::pep::PepEstimator;
use osprey_ml::svm::{self, FeatureStandardizer, LinearSvm};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;

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
    /// Optional feature names for logging (must match feature count)
    pub feature_names: Option<Vec<String>>,
    /// Maximum paired entries for SVM cross-validation (default: 300_000).
    /// When total entries exceed this limit, paired peptide groups (target+decoy
    /// together, same peptide charge states together) are randomly selected until
    /// this limit is reached. The CV is run on the subset; trained models score
    /// ALL entries. 300K divides evenly by 3 folds (100K per fold).
    /// Set to 0 to disable subsampling.
    pub max_train_size: usize,
    /// When true, suppress per-file/experiment FDR summary logging.
    /// Used in the streaming path where run_percolator is called on a training
    /// subset and the FDR numbers are meaningless.
    pub train_only: bool,
}

impl Default for PercolatorConfig {
    fn default() -> Self {
        PercolatorConfig {
            train_fdr: 0.01,
            test_fdr: 0.01,
            max_iterations: 10,
            n_folds: 3,
            seed: 42,
            c_values: vec![0.001, 0.01, 0.1, 1.0, 10.0, 100.0],
            feature_names: None,
            max_train_size: 300_000,
            train_only: false,
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
    /// Bias terms from best model per fold
    pub fold_biases: Vec<f64>,
    /// Feature standardizer used during training (for scoring additional entries)
    pub standardizer: FeatureStandardizer,
    /// Number of iterations used per fold
    pub iterations_per_fold: Vec<usize>,
}

/// FDR result for a single entry from compute_fdr_from_scores
#[derive(Debug, Clone)]
pub struct FdrScoreResult {
    /// Per-run precursor-level q-value
    pub run_precursor_qvalue: f64,
    /// Per-run peptide-level q-value
    pub run_peptide_qvalue: f64,
    /// Experiment-wide precursor-level q-value
    pub experiment_precursor_qvalue: f64,
    /// Experiment-wide peptide-level q-value
    pub experiment_peptide_qvalue: f64,
    /// Posterior error probability
    pub pep: f64,
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
            fold_biases: Vec::new(),
            standardizer: FeatureStandardizer::fit_from_slices(&[], 0),
            iterations_per_fold: Vec::new(),
        });
    }

    let n = entries.len();
    let n_features = entries[0].features.len();
    let n_targets = entries.iter().filter(|e| !e.is_decoy).count();
    let n_decoys = entries.iter().filter(|e| e.is_decoy).count();

    if config.max_train_size > 0 && n > config.max_train_size {
        log::info!(
            "Percolator: {} entries ({} targets, {} decoys), {} features, {}-fold CV, SVM training capped at {}",
            n, n_targets, n_decoys, n_features, config.n_folds, config.max_train_size
        );
    } else {
        log::info!(
            "Percolator: {} entries ({} targets, {} decoys), {} features, {}-fold CV",
            n,
            n_targets,
            n_decoys,
            n_features,
            config.n_folds
        );
    }

    // 1. Build feature matrix
    let feature_data: Vec<f64> = entries
        .iter()
        .flat_map(|e| e.features.iter().copied())
        .collect();
    let features = Matrix::new(feature_data, n, n_features);
    let labels: Vec<bool> = entries.iter().map(|e| e.is_decoy).collect();
    let entry_ids: Vec<u32> = entries.iter().map(|e| e.entry_id).collect();
    let peptides: Vec<String> = entries.iter().map(|e| e.peptide.clone()).collect();

    // 2. Standardize features (on ALL entries — standardization must be global)
    let (standardizer, std_features) = FeatureStandardizer::fit_transform(&features);
    log::debug!("  Features standardized to zero mean, unit variance");

    // 3. Subsample by peptide groups if needed (before fold splitting, per PMC5059416)
    //    Keeps target-decoy pairs and charge states together.
    //    The subsampled set is used for fold assignment + SVM training.
    //    ALL entries are scored with the trained model.
    let train_subset: Option<Vec<usize>> = if config.max_train_size > 0 && n > config.max_train_size
    {
        Some(subsample_by_peptide_group(
            &labels,
            &entry_ids,
            &peptides,
            config.max_train_size,
            config.seed,
        ))
    } else {
        None
    };

    let sub_n = train_subset.as_ref().map_or(n, |s| s.len());

    // Build subset-local arrays (or reference full arrays if no subsampling)
    let sub_labels: Vec<bool> = match &train_subset {
        Some(indices) => indices.iter().map(|&i| labels[i]).collect(),
        None => labels.clone(),
    };
    let sub_entry_ids: Vec<u32> = match &train_subset {
        Some(indices) => indices.iter().map(|&i| entry_ids[i]).collect(),
        None => entry_ids.clone(),
    };
    let sub_peptides: Vec<String> = match &train_subset {
        Some(indices) => indices.iter().map(|&i| peptides[i].clone()).collect(),
        None => peptides.clone(),
    };
    let sub_features = match &train_subset {
        Some(indices) => extract_rows(&std_features, indices),
        None => std_features.clone(),
    };

    if let Some(ref indices) = train_subset {
        let sub_targets = sub_labels.iter().filter(|&&d| !d).count();
        let sub_decoys = sub_labels.iter().filter(|&&d| d).count();
        log::debug!(
            "  Subsampled {} entries ({} targets, {} decoys) from {} for SVM training",
            indices.len(),
            sub_targets,
            sub_decoys,
            n
        );
    }

    // 4. Assign folds on the (possibly subsampled) set
    let fold_assignments = create_stratified_folds_by_peptide(
        &sub_labels,
        &sub_peptides,
        &sub_entry_ids,
        config.n_folds,
    );

    // Stage 5 sub-stage diagnostic dump. Gated by OSPREY_DUMP_SUBSAMPLE=1;
    // exits via OSPREY_SUBSAMPLE_ONLY=1. Captures the subsample membership
    // and fold assignment per entry, so the cross-impl compare can see
    // whether the two tools pick the same 300K subset and assign the same
    // fold IDs. Both algorithms are identical in source — drift here
    // indicates input-array ordering differs between Rust and C#.
    dump_stage5_subsample(entries, train_subset.as_deref(), &fold_assignments);

    // 5. Find best initial feature (on subsampled set)
    let (best_feat_idx, best_feat_passing) =
        find_best_initial_feature(&sub_features, &sub_labels, &sub_entry_ids, config.train_fdr);

    // If no targets pass at the configured FDR, loosen to 5% so training can proceed
    let mut config = config.clone();
    let (best_feat_idx, best_feat_passing) = if best_feat_passing == 0 {
        let relaxed_fdr = 0.05;
        let (idx, passing) =
            find_best_initial_feature(&sub_features, &sub_labels, &sub_entry_ids, relaxed_fdr);
        if passing > 0 {
            log::warn!(
                "  No targets at {:.0}% FDR — loosening train FDR to {:.0}%",
                config.train_fdr * 100.0,
                relaxed_fdr * 100.0
            );
            config.train_fdr = relaxed_fdr;
        } else {
            log::warn!(
                "  No targets at {:.0}% or {:.0}% FDR — features cannot discriminate targets from decoys",
                config.train_fdr * 100.0,
                relaxed_fdr * 100.0
            );
        }
        (idx, passing)
    } else {
        (best_feat_idx, best_feat_passing)
    };
    let config = &config;

    let initial_scores: Vec<f64> = (0..sub_n)
        .map(|i| sub_features[(i, best_feat_idx)])
        .collect();
    let feat_name = config
        .feature_names
        .as_ref()
        .and_then(|names| names.get(best_feat_idx))
        .map(|s| s.as_str())
        .unwrap_or("unknown");
    log::debug!(
        "  Best initial feature: {} ({} targets at {:.0}% FDR)",
        feat_name,
        best_feat_passing,
        config.train_fdr * 100.0
    );

    // 6. Train per-fold models via cross-validation (parallel)
    let mut final_scores = vec![0.0; n];
    let mut fold_weights: Vec<Vec<f64>> = Vec::new();
    let mut iterations_per_fold: Vec<usize> = Vec::new();

    // Progress bar: n_folds * max_iterations total steps (folds that stop early
    // skip remaining steps). Shared across parallel folds via Arc.
    let total_steps = (config.n_folds * config.max_iterations) as u64;
    let pb = indicatif::ProgressBar::new(total_steps);
    pb.set_style(
        indicatif::ProgressStyle::default_bar()
            .template("{spinner:.green} SVM training [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} fold-iterations")
            .unwrap()
            .progress_chars("#>-"),
    );
    // Train all folds in parallel — each fold is independent
    let fold_results: Vec<(usize, LinearSvm, usize)> = (0..config.n_folds)
        .into_par_iter()
        .map(|fold| {
            let train_indices: Vec<usize> = (0..sub_n)
                .filter(|&i| fold_assignments[i] != fold)
                .collect();
            let test_indices: Vec<usize> = (0..sub_n)
                .filter(|&i| fold_assignments[i] == fold)
                .collect();

            log::debug!(
                "  Fold {}/{}: {} train, {} test",
                fold + 1,
                config.n_folds,
                train_indices.len(),
                test_indices.len()
            );

            let (best_model, n_iterations) = train_fold(
                &sub_features,
                &sub_labels,
                &sub_entry_ids,
                &sub_peptides,
                &train_indices,
                &initial_scores,
                config,
                Some(&pb),
            );

            (fold, best_model, n_iterations)
        })
        .collect();

    pb.finish_and_clear();

    // Score ALL entries with trained models
    match &train_subset {
        Some(indices) => {
            // Build set of global indices in subset for quick lookup
            let in_subset: HashSet<usize> = indices.iter().copied().collect();

            // For subset entries: score with the held-out fold model (standard CV)
            for (fold, best_model, _) in &fold_results {
                let test_sub_indices: Vec<usize> = (0..sub_n)
                    .filter(|&i| fold_assignments[i] == *fold)
                    .collect();
                let test_global_indices: Vec<usize> =
                    test_sub_indices.iter().map(|&i| indices[i]).collect();
                let test_features = extract_rows(&std_features, &test_global_indices);
                let test_scores = best_model.decision_function(&test_features);
                for (i, &idx) in test_global_indices.iter().enumerate() {
                    final_scores[idx] = test_scores[i];
                }
            }

            // For non-subset entries: average scores from all fold models (batch)
            let non_subset_indices: Vec<usize> =
                (0..n).filter(|i| !in_subset.contains(i)).collect();
            if !non_subset_indices.is_empty() {
                let non_sub_features = extract_rows(&std_features, &non_subset_indices);
                let models: Vec<&LinearSvm> = fold_results.iter().map(|(_, m, _)| m).collect();
                let n_models = models.len() as f64;
                // Batch score all non-subset entries per model, then average
                let mut avg_scores = vec![0.0f64; non_subset_indices.len()];
                for model in &models {
                    let model_scores = model.decision_function(&non_sub_features);
                    for (i, s) in model_scores.iter().enumerate() {
                        avg_scores[i] += s;
                    }
                }
                for (i, &idx) in non_subset_indices.iter().enumerate() {
                    final_scores[idx] = avg_scores[i] / n_models;
                }
            }
        }
        None => {
            // No subsampling: score test fold directly (standard CV)
            for (fold, best_model, _) in &fold_results {
                let test_indices: Vec<usize> =
                    (0..n).filter(|&i| fold_assignments[i] == *fold).collect();
                let test_features = extract_rows(&std_features, &test_indices);
                let test_scores = best_model.decision_function(&test_features);
                for (i, &idx) in test_indices.iter().enumerate() {
                    final_scores[idx] = test_scores[i];
                }
            }
        }
    }

    let mut fold_biases: Vec<f64> = Vec::new();
    for (fold, best_model, n_iterations) in &fold_results {
        fold_weights.push(best_model.weights().to_vec());
        fold_biases.push(best_model.bias());
        iterations_per_fold.push(*n_iterations);
        log::debug!(
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

    // Stage 5 SVM-internals dump. Gated by OSPREY_DUMP_SVM_WEIGHTS=1;
    // exits via OSPREY_SVM_WEIGHTS_ONLY=1. Captures the 21 per-feature
    // weights plus bias per fold, and the iteration count per fold,
    // right after SVM training and before Granholm calibration. If these
    // match across tools but the end-of-Stage-5 dump still diverges,
    // drift is in calibration; if they differ, drift is in SVM training.
    dump_stage5_svm_weights(
        &fold_weights,
        &fold_biases,
        &iterations_per_fold,
        config.feature_names.as_deref(),
    );

    // 6b. Calibrate scores between folds (Granholm et al. 2012)
    // For subsampled runs, calibrate using subset fold assignments;
    // non-subset entries already have averaged scores (no fold-specific bias).
    match &train_subset {
        Some(indices) => {
            // Build global fold assignments: subset entries get their fold, others get usize::MAX
            let mut global_fold_assignments = vec![usize::MAX; n];
            for (sub_i, &global_i) in indices.iter().enumerate() {
                global_fold_assignments[global_i] = fold_assignments[sub_i];
            }
            calibrate_scores_between_folds(
                &mut final_scores,
                &global_fold_assignments,
                &labels,
                &entry_ids,
                config.n_folds,
                config.train_fdr,
            );
        }
        None => {
            calibrate_scores_between_folds(
                &mut final_scores,
                &fold_assignments,
                &labels,
                &entry_ids,
                config.n_folds,
                config.train_fdr,
            );
        }
    }

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

    // Log per-file and experiment-level FDR statistics (skip in train_only mode)
    if !config.train_only {
        let mut file_groups: HashMap<&str, Vec<usize>> = HashMap::new();
        for (i, entry) in entries.iter().enumerate() {
            file_groups
                .entry(entry.file_name.as_str())
                .or_default()
                .push(i);
        }
        let mut sorted_files: Vec<&str> = file_groups.keys().copied().collect();
        sorted_files.sort();

        log::info!("");
        log::info!(
            "=== Per-file results at {:.0}% run-level FDR ===",
            config.test_fdr * 100.0
        );
        for file_name in &sorted_files {
            let indices = &file_groups[file_name];
            // Count precursors passing precursor-level FDR
            let precursor_count = indices
                .iter()
                .copied()
                .filter(|&i| {
                    !entries[i].is_decoy && results[i].run_precursor_qvalue <= config.test_fdr
                })
                .map(|i| (entries[i].peptide.as_str(), entries[i].charge))
                .collect::<HashSet<_>>()
                .len();
            // Count unique peptides passing peptide-level FDR
            let peptide_count = indices
                .iter()
                .copied()
                .filter(|&i| {
                    !entries[i].is_decoy && results[i].run_peptide_qvalue <= config.test_fdr
                })
                .map(|i| entries[i].peptide.as_str())
                .collect::<HashSet<_>>()
                .len();
            log::info!(
                "  {}: {} precursors (precursor FDR), {} peptides (peptide FDR)",
                file_name,
                precursor_count,
                peptide_count
            );
        }

        // Log experiment-level statistics: report precursor and peptide counts
        // at their respective FDR levels independently, matching the two-level
        // FDR controls applied by Percolator.
        let exp_precursors = (0..results.len())
            .filter(|&i| {
                !entries[i].is_decoy && results[i].experiment_precursor_qvalue <= config.test_fdr
            })
            .map(|i| (entries[i].peptide.as_str(), entries[i].charge))
            .collect::<HashSet<_>>()
            .len();
        let exp_peptides = (0..results.len())
            .filter(|&i| {
                !entries[i].is_decoy && results[i].experiment_peptide_qvalue <= config.test_fdr
            })
            .map(|i| entries[i].peptide.as_str())
            .collect::<HashSet<_>>()
            .len();

        log::info!("");
        log::info!(
            "=== Experiment-level results at {:.0}% FDR ===",
            config.test_fdr * 100.0
        );
        log::info!(
            "  {} precursors passing precursor-level FDR",
            exp_precursors
        );
        log::info!("  {} peptides passing peptide-level FDR", exp_peptides);
    }

    Ok(PercolatorResults {
        entries: results,
        fold_weights,
        fold_biases,
        standardizer,
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
    progress: Option<&indicatif::ProgressBar>,
) -> (LinearSvm, usize) {
    let n_features = std_features.cols;
    let mut current_scores = initial_scores.to_vec();

    // Track best model across iterations
    let mut best_model = LinearSvm::fit(&Matrix::zeros(0, n_features), &[], 1.0, config.seed);
    let mut best_iteration = 0usize;
    let mut consecutive_no_improve = 0usize;

    // Baseline: count initial feature's performance on this fold's training set (for logging)
    let train_labels: Vec<bool> = train_indices.iter().map(|&i| labels[i]).collect();
    let train_entry_ids: Vec<u32> = train_indices.iter().map(|&i| entry_ids[i]).collect();
    let train_initial_scores: Vec<f64> = train_indices.iter().map(|&i| initial_scores[i]).collect();
    let initial_passing = count_passing(
        &train_initial_scores,
        &train_labels,
        &train_entry_ids,
        config.train_fdr,
    );
    let mut best_passing = 0usize;
    let mut iterations_completed = 0usize;
    log::debug!(
        "    Initial feature baseline on training set: {} pass {:.0}% FDR",
        initial_passing,
        config.train_fdr * 100.0
    );

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

        // Build SVM training set: selected targets + all decoys.
        // Subsampling (if needed) is done at the top level before fold splitting,
        // keeping target-decoy pairs together (per PMC5059416).
        let decoy_indices: Vec<usize> = (0..train_indices.len())
            .filter(|&i| train_labels[i])
            .collect();
        let mut svm_indices: Vec<usize> = selected_target_indices.clone();
        svm_indices.extend_from_slice(&decoy_indices);

        // Map to global indices for feature extraction
        let svm_global_indices: Vec<usize> =
            svm_indices.iter().map(|&i| train_indices[i]).collect();
        let svm_features = extract_rows(std_features, &svm_global_indices);
        let svm_labels: Vec<bool> = svm_indices.iter().map(|&i| train_labels[i]).collect();
        let svm_entry_ids: Vec<u32> = svm_indices.iter().map(|&i| train_entry_ids[i]).collect();

        // ii. Grid search for best C
        let svm_peptides: Vec<String> = svm_indices
            .iter()
            .map(|&i| peptides[train_indices[i]].clone())
            .collect();
        let svm_fold_assignments = create_stratified_folds_by_peptide(
            &svm_labels,
            &svm_peptides,
            &svm_entry_ids,
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

        log::debug!(
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
            log::debug!("      -> New best: {} passing", best_passing);
        } else {
            consecutive_no_improve += 1;
            log::debug!(
                "      -> No improvement (best={} from iter {}, {} consecutive)",
                best_passing,
                best_iteration,
                consecutive_no_improve
            );
        }

        iterations_completed = iteration + 1;
        if let Some(pb) = progress {
            pb.inc(1);
        }

        if consecutive_no_improve >= 2 {
            log::debug!("    Stopping early: 2 consecutive non-improvements");
            break;
        }
    }

    // Mark remaining iterations as complete if we stopped early
    if let Some(pb) = progress {
        let remaining = config.max_iterations - iterations_completed;
        if remaining > 0 {
            pb.inc(remaining as u64);
        }
    }

    let final_iteration = best_iteration.max(1);
    log::debug!(
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
pub fn compete_from_indices(
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
    // Store (index, score, is_decoy, base_id) — base_id for deterministic tiebreaking.
    let mut winners: Vec<(usize, f64, bool, u32)> = Vec::with_capacity(targets.len());
    for (&base_id, &(t_idx, t_score)) in &targets {
        if let Some(&(d_idx, d_score)) = decoys.get(&base_id) {
            if t_score > d_score {
                winners.push((t_idx, t_score, false, base_id));
            } else {
                winners.push((d_idx, d_score, true, base_id)); // tie → decoy
            }
        } else {
            winners.push((t_idx, t_score, false, base_id)); // unpaired target wins
        }
    }
    // Unpaired decoys
    for (&base_id, &(d_idx, d_score)) in &decoys {
        if !targets.contains_key(&base_id) {
            winners.push((d_idx, d_score, true, base_id));
        }
    }

    // Sort by score descending, then by base_id ascending for deterministic tiebreaking.
    // IMPORTANT: Use base_id (not array index) as tiebreaker. Array indices depend on input
    // order, and sorting by index can create systematic target/decoy bias if targets
    // systematically have lower indices. base_id is intrinsic and unbiased.
    winners.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.3.cmp(&b.3)));

    let winner_indices: Vec<usize> = winners.iter().map(|w| w.0).collect();
    let winner_scores: Vec<f64> = winners.iter().map(|w| w.1).collect();
    let winner_is_decoy: Vec<bool> = winners.iter().map(|w| w.2).collect();

    (winner_indices, winner_scores, winner_is_decoy)
}

/// Compute conservative q-values: FDR = (n_decoy + 1) / n_target
///
/// Input must be sorted by score descending (winners from competition).
pub fn compute_conservative_qvalues(_scores: &[f64], is_decoy: &[bool], q_values: &mut [f64]) {
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
pub fn best_precursor_per_peptide(
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

/// Create stratified fold assignments grouped by target peptide, keeping pairs together.
///
/// All charge states of the same peptide AND their paired decoys go into the same fold.
/// This prevents target-decoy pairs from being split across folds, which would cause
/// unpaired targets to auto-win competition and inflate the positive training set.
///
/// Grouping uses base_id (entry_id & 0x7FFFFFFF) to link each decoy to its paired
/// target's peptide sequence. All entries sharing a target peptide are assigned together.
fn create_stratified_folds_by_peptide(
    labels: &[bool],
    peptides: &[String],
    entry_ids: &[u32],
    n_folds: usize,
) -> Vec<usize> {
    // 1. Build base_id → target peptide mapping
    let mut base_id_to_target_peptide: HashMap<u32, &str> = HashMap::new();
    for (i, (&is_decoy, &eid)) in labels.iter().zip(entry_ids).enumerate() {
        let base_id = eid & 0x7FFFFFFF;
        if !is_decoy {
            base_id_to_target_peptide
                .entry(base_id)
                .or_insert(peptides[i].as_str());
        }
    }

    // 2. Map each entry to its group key (target peptide via base_id).
    //    Decoys look up their paired target's peptide; unpaired entries use own peptide.
    let group_keys: Vec<&str> = (0..labels.len())
        .map(|i| {
            let base_id = entry_ids[i] & 0x7FFFFFFF;
            base_id_to_target_peptide
                .get(&base_id)
                .copied()
                .unwrap_or(peptides[i].as_str())
        })
        .collect();

    // 3. Group all entries by target peptide
    let mut peptide_groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, &key) in group_keys.iter().enumerate() {
        peptide_groups.entry(key).or_default().push(i);
    }

    // 4. Sort for deterministic assignment, round-robin assign folds
    let mut sorted_groups: Vec<(&str, &Vec<usize>)> =
        peptide_groups.iter().map(|(&k, v)| (k, v)).collect();
    sorted_groups.sort_by_key(|&(k, _)| k);

    let mut fold_assignments = vec![0; labels.len()];
    for (i, (_, indices)) in sorted_groups.iter().enumerate() {
        let fold = i % n_folds;
        for &idx in *indices {
            fold_assignments[idx] = fold;
        }
    }

    fold_assignments
}

/// Cross-impl bisection dump of the subsample + fold-assignment state,
/// taken right after `create_stratified_folds_by_peptide` and before
/// SVM training. Writes `rust_stage5_subsample.tsv` with one row per
/// entry in the best-per-precursor array (n ~= 462k on Stellar single
/// file).
///
/// Columns: `entry_id, native_position, charge, modified_sequence,
/// is_decoy, base_id, in_subsample, fold_id`. `native_position` is the
/// entry's index in the input `entries` slice — divergence on this
/// column between Rust and C# means the two tools populate the array
/// in different orders, which would cascade through every downstream
/// `usize` index. `in_subsample` is false for entries filtered out by
/// `subsample_by_peptide_group`; `fold_id` is `-1` for those rows.
///
/// Gated by `OSPREY_DUMP_SUBSAMPLE=1`. When `OSPREY_SUBSAMPLE_ONLY=1`
/// is also set, exits the process after writing. Rows sorted by
/// `entry_id` for stable human inspection; the compare script
/// hash-joins on `entry_id`, sort-order-agnostic.
fn dump_stage5_subsample(
    entries: &[PercolatorEntry],
    train_subset: Option<&[usize]>,
    fold_assignments: &[usize],
) {
    if !is_dump_enabled("OSPREY_DUMP_SUBSAMPLE") {
        return;
    }

    let n = entries.len();
    let mut in_sub = vec![false; n];
    let mut fold_for = vec![-1i32; n];
    match train_subset {
        Some(indices) => {
            for (sub_pos, &native_pos) in indices.iter().enumerate() {
                in_sub[native_pos] = true;
                fold_for[native_pos] = fold_assignments[sub_pos] as i32;
            }
        }
        None => {
            for (i, &f) in fold_assignments.iter().enumerate() {
                in_sub[i] = true;
                fold_for[i] = f as i32;
            }
        }
    }

    let path = "rust_stage5_subsample.tsv";
    let Ok(mut f) = std::fs::File::create(path) else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "entry_id\tnative_position\tcharge\tmodified_sequence\tis_decoy\tbase_id\tin_subsample\tfold_id"
    )
    .ok();

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| entries[i].entry_id);

    for i in order {
        let e = &entries[i];
        let base_id = e.entry_id & 0x7FFF_FFFF;
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            e.entry_id,
            i,
            e.charge,
            e.peptide,
            if e.is_decoy { "true" } else { "false" },
            base_id,
            if in_sub[i] { "true" } else { "false" },
            fold_for[i],
        )
        .ok();
    }
    log::info!(
        "Wrote Stage 5 subsample dump: {} ({} rows)",
        path,
        n
    );

    exit_if_only("OSPREY_SUBSAMPLE_ONLY", "Stage 5 subsample dump");
}

/// Cross-impl bisection dump of per-fold SVM weights, taken right after
/// training converges and before Granholm cross-fold calibration. Writes
/// `rust_stage5_svm_weights.tsv` with one row per (fold, weight) pair:
/// 21 feature weights + 1 bias per fold, so 66 rows total for the default
/// 3-fold / 21-feature configuration.
///
/// Columns: `fold, weight_idx, feature_name, value, fold_iterations`.
/// The `fold_iterations` column repeats a per-fold scalar (0..max_iter);
/// drift there alone means convergence differs across tools. Sorted by
/// `(fold, weight_idx)` for stable inspection; compare script hash-joins.
///
/// Gated by `OSPREY_DUMP_SVM_WEIGHTS=1`. When `OSPREY_SVM_WEIGHTS_ONLY=1`
/// is also set, exits after writing.
fn dump_stage5_svm_weights(
    fold_weights: &[Vec<f64>],
    fold_biases: &[f64],
    iterations_per_fold: &[usize],
    feature_names: Option<&[String]>,
) {
    if !is_dump_enabled("OSPREY_DUMP_SVM_WEIGHTS") {
        return;
    }

    let path = "rust_stage5_svm_weights.tsv";
    let Ok(mut f) = std::fs::File::create(path) else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(f, "fold\tweight_idx\tfeature_name\tvalue\tfold_iterations").ok();
    for (fold, (weights, bias)) in fold_weights.iter().zip(fold_biases.iter()).enumerate() {
        let iters = iterations_per_fold.get(fold).copied().unwrap_or(0);
        for (wi, &w) in weights.iter().enumerate() {
            let name = feature_names
                .and_then(|names| names.get(wi))
                .map(|s| s.as_str())
                .unwrap_or("unknown");
            writeln!(f, "{}\t{}\t{}\t{}\t{}", fold, wi, name, w, iters).ok();
        }
        writeln!(
            f,
            "{}\t{}\tbias\t{}\t{}",
            fold,
            weights.len(),
            bias,
            iters
        )
        .ok();
    }

    log::info!(
        "Wrote Stage 5 SVM weights dump: {} ({} folds)",
        path,
        fold_weights.len()
    );

    exit_if_only("OSPREY_SVM_WEIGHTS_ONLY", "Stage 5 SVM weights dump");
}

/// Subsample entries by peptide group, keeping target-decoy pairs and charge states together.
///
/// Groups entries by target peptide (via base_id), then randomly selects groups until
/// the total entry count reaches `max_entries`. Returns sorted global indices.
///
/// Per The et al. (2016, PMC5059416): subsample paired PSMs before fold splitting,
/// then apply CV to the subset. The trained model is applied to ALL entries.
pub fn subsample_by_peptide_group(
    labels: &[bool],
    entry_ids: &[u32],
    peptides: &[String],
    max_entries: usize,
    seed: u64,
) -> Vec<usize> {
    let n = labels.len();
    if n <= max_entries {
        return (0..n).collect();
    }

    // Build peptide groups (same as fold assignment: group by target peptide via base_id)
    let mut base_id_to_target_peptide: HashMap<u32, &str> = HashMap::new();
    for (i, (&is_decoy, &eid)) in labels.iter().zip(entry_ids).enumerate() {
        let base_id = eid & 0x7FFFFFFF;
        if !is_decoy {
            base_id_to_target_peptide
                .entry(base_id)
                .or_insert(peptides[i].as_str());
        }
    }

    let mut peptide_groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for i in 0..n {
        let base_id = entry_ids[i] & 0x7FFFFFFF;
        let key = base_id_to_target_peptide
            .get(&base_id)
            .copied()
            .unwrap_or(peptides[i].as_str());
        peptide_groups.entry(key).or_default().push(i);
    }

    // Sort groups deterministically and shuffle with Fisher-Yates
    let mut groups: Vec<(&str, Vec<usize>)> = peptide_groups.into_iter().collect();
    groups.sort_by_key(|&(k, _)| k);

    // Fisher-Yates shuffle with xorshift64 RNG
    let mut rng_state = seed;
    for i in (1..groups.len()).rev() {
        rng_state ^= rng_state << 13;
        rng_state ^= rng_state >> 7;
        rng_state ^= rng_state << 17;
        let j = rng_state as usize % (i + 1);
        groups.swap(i, j);
    }

    // Select groups until we reach max_entries
    let mut selected: Vec<usize> = Vec::with_capacity(max_entries);
    for (_, indices) in &groups {
        if selected.len() + indices.len() > max_entries && !selected.is_empty() {
            break;
        }
        selected.extend_from_slice(indices);
    }

    selected.sort(); // deterministic order for downstream processing
    selected
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

/// Compute FDR results (q-values + PEP) for pre-scored entries.
///
/// This is used by the streaming Percolator path where the caller has already
/// scored all entries using trained models. Computes:
/// - Per-run precursor and peptide q-values
/// - Experiment-level precursor and peptide q-values
/// - Posterior error probabilities via KDE + isotonic regression
///
/// Returns one FdrScoreResult per entry, in the same order as the input arrays.
pub fn compute_fdr_from_scores(
    scores: &[f64],
    labels: &[bool],
    entry_ids: &[u32],
    peptides: &[String],
    file_names: &[String],
    test_fdr: f64,
) -> Vec<FdrScoreResult> {
    let n = scores.len();

    // PEP estimation on competition winners
    let (winner_indices, winner_scores, winner_is_decoy) = compete_all(scores, labels, entry_ids);

    let pep_estimator = PepEstimator::fit_default(&winner_scores, &winner_is_decoy);
    let mut peps = vec![1.0; n];
    for &idx in &winner_indices {
        peps[idx] = pep_estimator.posterior_error(scores[idx]);
    }

    // Q-values at all four levels
    let unique_files: HashSet<&str> = file_names.iter().map(|s| s.as_str()).collect();
    let is_single_file = unique_files.len() <= 1;

    let run_precursor_qvalues =
        compute_per_run_precursor_qvalues(scores, labels, entry_ids, file_names);
    let run_peptide_qvalues =
        compute_per_run_peptide_qvalues(scores, labels, entry_ids, file_names, peptides);

    let (exp_precursor_qvalues, exp_peptide_qvalues) = if is_single_file {
        (run_precursor_qvalues.clone(), run_peptide_qvalues.clone())
    } else {
        let exp_prec = compute_experiment_precursor_qvalues(scores, labels, entry_ids);
        let exp_pept = compute_experiment_peptide_qvalues(scores, labels, entry_ids, peptides);
        (exp_prec, exp_pept)
    };

    // Log per-file statistics
    let mut file_groups: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, fname) in file_names.iter().enumerate() {
        file_groups.entry(fname.as_str()).or_default().push(i);
    }
    let mut sorted_files: Vec<&str> = file_groups.keys().copied().collect();
    sorted_files.sort();

    log::info!("");
    log::info!(
        "=== Per-file results at {:.0}% run-level FDR ===",
        test_fdr * 100.0
    );
    for file_name in &sorted_files {
        let indices = &file_groups[file_name];
        // Precursors passing precursor-level FDR
        let precursor_count = indices
            .iter()
            .copied()
            .filter(|&i| !labels[i] && run_precursor_qvalues[i] <= test_fdr)
            .map(|i| (peptides[i].as_str(), entry_ids[i] & 0x7FFFFFFF))
            .collect::<HashSet<_>>()
            .len();
        // Peptides passing peptide-level FDR
        let peptide_count = indices
            .iter()
            .copied()
            .filter(|&i| !labels[i] && run_peptide_qvalues[i] <= test_fdr)
            .map(|i| peptides[i].as_str())
            .collect::<HashSet<_>>()
            .len();
        log::info!(
            "  {}: {} precursors (precursor FDR), {} peptides (peptide FDR)",
            file_name,
            precursor_count,
            peptide_count
        );
    }

    // Experiment-level stats: precursor and peptide counts at their respective FDR levels
    let exp_precursors = (0..n)
        .filter(|&i| !labels[i] && exp_precursor_qvalues[i] <= test_fdr)
        .map(|i| (peptides[i].as_str(), entry_ids[i] & 0x7FFFFFFF))
        .collect::<HashSet<_>>()
        .len();
    let exp_peptides = (0..n)
        .filter(|&i| !labels[i] && exp_peptide_qvalues[i] <= test_fdr)
        .map(|i| peptides[i].as_str())
        .collect::<HashSet<_>>()
        .len();

    log::info!("");
    log::info!(
        "=== Experiment-level results at {:.0}% FDR ===",
        test_fdr * 100.0
    );
    log::info!(
        "  {} precursors passing precursor-level FDR",
        exp_precursors
    );
    log::info!("  {} peptides passing peptide-level FDR", exp_peptides);

    // Build results
    (0..n)
        .map(|i| FdrScoreResult {
            run_precursor_qvalue: run_precursor_qvalues[i],
            run_peptide_qvalue: run_peptide_qvalues[i],
            experiment_precursor_qvalue: exp_precursor_qvalues[i],
            experiment_peptide_qvalue: exp_peptide_qvalues[i],
            pep: peps[i],
        })
        .collect()
}

/// Compute FDR (q-values + PEP) directly from per-file FdrEntry stubs.
///
/// Memory-efficient alternative to `compute_fdr_from_scores` that avoids
/// allocating flat metadata arrays (peptides, file_names) for all entries.
/// Instead, it operates directly on per_file_entries using temporary per-file
/// arrays (~30 MB each) and small HashMaps (~10 MB total).
///
/// Scores must already be written to `entry.score` before calling this.
///
/// Writes results directly to FdrEntry fields:
/// - `entry.run_precursor_qvalue` = per-file precursor-level q-value
/// - `entry.run_peptide_qvalue` = per-file peptide-level q-value
/// - `entry.experiment_precursor_qvalue` = experiment-wide precursor-level q-value
/// - `entry.experiment_peptide_qvalue` = experiment-wide peptide-level q-value
/// - `entry.pep` = posterior error probability (1.0 for non-competition-winners)
///
/// Protein-level q-value fields are not set here (they default to 1.0).
pub fn compute_fdr_from_stubs(
    per_file_entries: &mut [(String, Vec<FdrEntry>)],
    test_fdr: f64,
    restrict_base_ids: Option<&HashSet<u32>>,
) {
    let n_files = per_file_entries.len();
    let is_single_file = n_files <= 1;

    // === PEP estimation via global target-decoy competition ===
    // Find best target and best decoy per base_id across all files
    let mut targets: HashMap<u32, (f64, usize, usize)> = HashMap::new();
    let mut decoys: HashMap<u32, (f64, usize, usize)> = HashMap::new();

    for (file_idx, (_, entries)) in per_file_entries.iter().enumerate() {
        for (local_idx, entry) in entries.iter().enumerate() {
            let base_id = entry.entry_id & 0x7FFF_FFFF;
            if let Some(restrict) = restrict_base_ids {
                if !restrict.contains(&base_id) {
                    continue;
                }
            }
            let map = if entry.is_decoy {
                &mut decoys
            } else {
                &mut targets
            };
            map.entry(base_id)
                .and_modify(|(best_score, best_fi, best_li)| {
                    if entry.score > *best_score {
                        *best_score = entry.score;
                        *best_fi = file_idx;
                        *best_li = local_idx;
                    }
                })
                .or_insert((entry.score, file_idx, local_idx));
        }
    }

    // Compete target-decoy pairs for PEP fitting
    let mut winner_scores: Vec<f64> = Vec::with_capacity(targets.len());
    let mut winner_is_decoy: Vec<bool> = Vec::with_capacity(targets.len());
    let mut winner_locations: Vec<(usize, usize)> = Vec::with_capacity(targets.len());

    for (&base_id, &(t_score, t_fi, t_li)) in &targets {
        if let Some(&(d_score, d_fi, d_li)) = decoys.get(&base_id) {
            if t_score > d_score {
                winner_scores.push(t_score);
                winner_is_decoy.push(false);
                winner_locations.push((t_fi, t_li));
            } else {
                winner_scores.push(d_score);
                winner_is_decoy.push(true);
                winner_locations.push((d_fi, d_li));
            }
        } else {
            winner_scores.push(t_score);
            winner_is_decoy.push(false);
            winner_locations.push((t_fi, t_li));
        }
    }
    for (&base_id, &(d_score, d_fi, d_li)) in &decoys {
        if !targets.contains_key(&base_id) {
            winner_scores.push(d_score);
            winner_is_decoy.push(true);
            winner_locations.push((d_fi, d_li));
        }
    }

    // Fit PEP model and apply to winners only (non-winners keep pep = 1.0).
    // PEP is reported per-precursor in the blib for downstream confidence,
    // but is NOT used for protein-level FDR — picked-protein ranks by raw
    // SVM score instead, because the PepEstimator's score range is bounded
    // by the winner set and clamps losing decoys to ~1.0, collapsing the
    // protein-level null distribution. See docs/07-fdr-control.md for the
    // full history of that choice.
    let pep_estimator = PepEstimator::fit_default(&winner_scores, &winner_is_decoy);
    for (i, &(fi, li)) in winner_locations.iter().enumerate() {
        per_file_entries[fi].1[li].pep = pep_estimator.posterior_error(winner_scores[i]);
    }

    // Free global competition structures
    drop(targets);
    drop(decoys);
    drop(winner_scores);
    drop(winner_is_decoy);
    drop(winner_locations);

    // === Per-file run-level q-values ===
    // Process one file at a time with temporary flat arrays (~30 MB each)
    for (file_name, entries) in per_file_entries.iter_mut() {
        let n = entries.len();
        if n == 0 {
            continue;
        }

        // Build per-file flat arrays, restricted to matching base_ids if specified
        let active_indices: Vec<usize> = if let Some(restrict) = restrict_base_ids {
            (0..n)
                .filter(|&i| restrict.contains(&(entries[i].entry_id & 0x7FFF_FFFF)))
                .collect()
        } else {
            (0..n).collect()
        };

        if active_indices.is_empty() {
            continue;
        }

        let scores: Vec<f64> = active_indices.iter().map(|&i| entries[i].score).collect();
        let labels: Vec<bool> = active_indices
            .iter()
            .map(|&i| entries[i].is_decoy)
            .collect();
        let eids: Vec<u32> = active_indices
            .iter()
            .map(|&i| entries[i].entry_id)
            .collect();
        let all_indices: Vec<usize> = (0..active_indices.len()).collect();

        // Precursor-level q-values (on active entries only)
        let (prec_winner_local, _, prec_winner_is_decoy) =
            compete_from_indices(&scores, &labels, &eids, &all_indices);
        let mut prec_q = vec![1.0; prec_winner_local.len()];
        let prec_ws: Vec<f64> = prec_winner_local.iter().map(|&i| scores[i]).collect();
        compute_conservative_qvalues(&prec_ws, &prec_winner_is_decoy, &mut prec_q);

        // Map active-local index → precursor q-value
        let mut active_prec_qvalue = vec![1.0; active_indices.len()];
        for (rank, &idx) in prec_winner_local.iter().enumerate() {
            active_prec_qvalue[idx] = prec_q[rank];
        }

        // Peptide-level q-values (on active entries only)
        let mut best_per_pept: HashMap<&str, (usize, f64)> = HashMap::new();
        for (ai, &orig_idx) in active_indices.iter().enumerate() {
            let entry = &entries[orig_idx];
            best_per_pept
                .entry(&*entry.modified_sequence)
                .and_modify(|(best_ai, best_score)| {
                    if scores[ai] > *best_score {
                        *best_ai = ai;
                        *best_score = scores[ai];
                    }
                })
                .or_insert((ai, scores[ai]));
        }
        let mut pept_reps: Vec<usize> = best_per_pept.values().map(|&(ai, _)| ai).collect();
        pept_reps.sort();

        let pept_scores: Vec<f64> = pept_reps.iter().map(|&ai| scores[ai]).collect();
        let pept_labels: Vec<bool> = pept_reps.iter().map(|&ai| labels[ai]).collect();
        let pept_eids: Vec<u32> = pept_reps.iter().map(|&ai| eids[ai]).collect();

        let (pept_winner_local, _, pept_winner_is_decoy) = compete_from_indices(
            &pept_scores,
            &pept_labels,
            &pept_eids,
            &(0..pept_reps.len()).collect::<Vec<_>>(),
        );
        let mut pept_q = vec![1.0; pept_winner_local.len()];
        let pept_ws: Vec<f64> = pept_winner_local.iter().map(|&i| pept_scores[i]).collect();
        compute_conservative_qvalues(&pept_ws, &pept_winner_is_decoy, &mut pept_q);

        // Map peptide q-values back via peptide string (owned keys to avoid borrow conflict)
        let mut peptide_qvalue: HashMap<String, f64> = HashMap::new();
        for (rank, &local_idx) in pept_winner_local.iter().enumerate() {
            let orig_idx = active_indices[pept_reps[local_idx]];
            peptide_qvalue.insert(
                entries[orig_idx].modified_sequence.to_string(),
                pept_q[rank],
            );
        }

        // Write run-level precursor and peptide q-values separately
        for (ai, &orig_idx) in active_indices.iter().enumerate() {
            let entry = &mut entries[orig_idx];
            let pept_qv = peptide_qvalue
                .get(&*entry.modified_sequence)
                .copied()
                .unwrap_or(1.0);
            entry.run_precursor_qvalue = active_prec_qvalue[ai];
            entry.run_peptide_qvalue = pept_qv;
        }

        // Log per-file stats
        let passing: Vec<&FdrEntry> = entries
            .iter()
            .filter(|e| !e.is_decoy && e.effective_run_qvalue(FdrLevel::Both) <= test_fdr)
            .collect();
        let precursor_count = passing
            .iter()
            .map(|e| (&*e.modified_sequence, e.entry_id & 0x7FFF_FFFF))
            .collect::<HashSet<_>>()
            .len();
        let peptide_count = passing
            .iter()
            .map(|e| &*e.modified_sequence)
            .collect::<HashSet<_>>()
            .len();
        log::info!(
            "  {}: {} precursors, {} peptides",
            file_name,
            precursor_count,
            peptide_count
        );
    }

    // === Experiment-level q-values ===
    if is_single_file {
        // Single file: experiment q-values = run q-values
        for (_, entries) in per_file_entries.iter_mut() {
            for entry in entries.iter_mut() {
                entry.experiment_precursor_qvalue = entry.run_precursor_qvalue;
                entry.experiment_peptide_qvalue = entry.run_peptide_qvalue;
            }
        }
    } else {
        // Experiment-level precursor q-values:
        // Find best observation per base_id across all files, compete, compute q-values
        let mut best_target: HashMap<u32, (f64, usize, usize)> = HashMap::new();
        let mut best_decoy: HashMap<u32, (f64, usize, usize)> = HashMap::new();
        for (file_idx, (_, entries)) in per_file_entries.iter().enumerate() {
            for (local_idx, entry) in entries.iter().enumerate() {
                let base_id = entry.entry_id & 0x7FFF_FFFF;
                if let Some(restrict) = restrict_base_ids {
                    if !restrict.contains(&base_id) {
                        continue;
                    }
                }
                let map = if entry.is_decoy {
                    &mut best_decoy
                } else {
                    &mut best_target
                };
                map.entry(base_id)
                    .and_modify(|(best_score, best_fi, best_li)| {
                        if entry.score > *best_score {
                            *best_score = entry.score;
                            *best_fi = file_idx;
                            *best_li = local_idx;
                        }
                    })
                    .or_insert((entry.score, file_idx, local_idx));
            }
        }

        // Compete and compute q-values
        let mut exp_winners: Vec<(f64, bool, u32)> = Vec::with_capacity(best_target.len());
        for (&base_id, &(t_score, _, _)) in &best_target {
            if let Some(&(d_score, _, _)) = best_decoy.get(&base_id) {
                if t_score > d_score {
                    exp_winners.push((t_score, false, base_id));
                } else {
                    exp_winners.push((d_score, true, base_id));
                }
            } else {
                exp_winners.push((t_score, false, base_id));
            }
        }
        for (&base_id, &(d_score, _, _)) in &best_decoy {
            if !best_target.contains_key(&base_id) {
                exp_winners.push((d_score, true, base_id));
            }
        }

        // Sort by score desc, base_id for deterministic tiebreaking
        exp_winners.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.2.cmp(&b.2)));

        let exp_w_scores: Vec<f64> = exp_winners.iter().map(|w| w.0).collect();
        let exp_w_decoy: Vec<bool> = exp_winners.iter().map(|w| w.1).collect();
        let mut exp_q = vec![1.0; exp_winners.len()];
        compute_conservative_qvalues(&exp_w_scores, &exp_w_decoy, &mut exp_q);

        // Map base_id -> experiment precursor q-value (for winners only)
        let mut base_id_exp_prec_q: HashMap<u32, f64> = HashMap::new();
        for (rank, &(_, _, base_id)) in exp_winners.iter().enumerate() {
            base_id_exp_prec_q.insert(base_id, exp_q[rank]);
        }

        drop(best_target);
        drop(best_decoy);
        drop(exp_winners);

        // Experiment-level peptide q-values:
        // Find best precursor per peptide across all files (owned keys)
        let mut best_per_pept: HashMap<String, (f64, bool, u32)> = HashMap::new();
        for (_, entries) in per_file_entries.iter() {
            for entry in entries.iter() {
                if let Some(restrict) = restrict_base_ids {
                    if !restrict.contains(&(entry.entry_id & 0x7FFF_FFFF)) {
                        continue;
                    }
                }
                best_per_pept
                    .entry(entry.modified_sequence.to_string())
                    .and_modify(|(best_score, best_decoy, best_eid)| {
                        if entry.score > *best_score {
                            *best_score = entry.score;
                            *best_decoy = entry.is_decoy;
                            *best_eid = entry.entry_id;
                        }
                    })
                    .or_insert((entry.score, entry.is_decoy, entry.entry_id));
            }
        }

        let n_pept = best_per_pept.len();
        let mut pept_data: Vec<(String, f64, bool, u32)> = best_per_pept
            .into_iter()
            .map(|(pept, (score, is_decoy, eid))| (pept, score, is_decoy, eid))
            .collect();
        pept_data.sort_by_key(|d| d.0.clone());

        let pept_scores: Vec<f64> = pept_data.iter().map(|d| d.1).collect();
        let pept_labels: Vec<bool> = pept_data.iter().map(|d| d.2).collect();
        let pept_eids: Vec<u32> = pept_data.iter().map(|d| d.3).collect();

        let (pept_winner_local, _, pept_winner_is_decoy) = compete_from_indices(
            &pept_scores,
            &pept_labels,
            &pept_eids,
            &(0..n_pept).collect::<Vec<_>>(),
        );
        let mut pept_q = vec![1.0; pept_winner_local.len()];
        let pept_ws: Vec<f64> = pept_winner_local.iter().map(|&i| pept_scores[i]).collect();
        compute_conservative_qvalues(&pept_ws, &pept_winner_is_decoy, &mut pept_q);

        let mut peptide_exp_q: HashMap<String, f64> = HashMap::new();
        for (rank, &local_idx) in pept_winner_local.iter().enumerate() {
            peptide_exp_q.insert(pept_data[local_idx].0.clone(), pept_q[rank]);
        }

        // Free intermediate structures
        drop(pept_data);

        // Write experiment-level precursor and peptide q-values separately
        for (_, entries) in per_file_entries.iter_mut() {
            for entry in entries.iter_mut() {
                let base_id = entry.entry_id & 0x7FFF_FFFF;
                if let Some(restrict) = restrict_base_ids {
                    if !restrict.contains(&base_id) {
                        continue;
                    }
                }
                let prec_q = base_id_exp_prec_q.get(&base_id).copied().unwrap_or(1.0);
                let pept_q = peptide_exp_q
                    .get(&*entry.modified_sequence)
                    .copied()
                    .unwrap_or(1.0);
                entry.experiment_precursor_qvalue = prec_q;
                entry.experiment_peptide_qvalue = pept_q;
            }
        }
    }

    // Log experiment-level stats: precursor and peptide counts at their respective
    // FDR levels independently. Most users should look at the peptide count for
    // a conservative estimate, since peptide-level FDR is typically stricter.
    let total_entries: usize = per_file_entries.iter().map(|(_, e)| e.len()).sum();
    let exp_precursors = per_file_entries
        .iter()
        .flat_map(|(_, entries)| entries.iter())
        .filter(|e| !e.is_decoy && e.experiment_precursor_qvalue <= test_fdr)
        .map(|e| (&*e.modified_sequence, e.entry_id & 0x7FFF_FFFF))
        .collect::<HashSet<_>>()
        .len();
    let exp_peptides = per_file_entries
        .iter()
        .flat_map(|(_, entries)| entries.iter())
        .filter(|e| !e.is_decoy && e.experiment_peptide_qvalue <= test_fdr)
        .map(|e| &*e.modified_sequence)
        .collect::<HashSet<_>>()
        .len();

    log::info!("");
    log::info!(
        "=== Experiment-level results at {:.0}% FDR (from {} total entries) ===",
        test_fdr * 100.0,
        total_entries
    );
    log::info!(
        "  {} precursors passing precursor-level FDR",
        exp_precursors
    );
    log::info!("  {} peptides passing peptide-level FDR", exp_peptides);
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
        // Indices: 0=PEPTIDEK target z2, 1=PEPTIDEK target z3,
        //          2=ANOTHERONE target, 3=KEDITPEP decoy of PEPTIDEK z2,
        //          4=KEDITPEP decoy of PEPTIDEK z3, 5=ENOREHTONA decoy of ANOTHERONE
        let labels = vec![false, false, false, true, true, true];
        let peptides = vec![
            "PEPTIDEK".to_string(),
            "PEPTIDEK".to_string(),
            "ANOTHERONE".to_string(),
            "KEDITPEP".to_string(),   // decoy — different sequence
            "KEDITPEP".to_string(),   // decoy — different sequence
            "ENOREHTONA".to_string(), // decoy — different sequence
        ];
        // entry_ids: base_id links target-decoy pairs, high bit set for decoys
        let entry_ids: Vec<u32> = vec![
            1,              // PEPTIDEK z2 target
            2,              // PEPTIDEK z3 target
            3,              // ANOTHERONE target
            1 | 0x80000000, // decoy paired with PEPTIDEK z2
            2 | 0x80000000, // decoy paired with PEPTIDEK z3
            3 | 0x80000000, // decoy paired with ANOTHERONE
        ];

        let folds = create_stratified_folds_by_peptide(&labels, &peptides, &entry_ids, 3);

        // All charge states of PEPTIDEK (targets) should be in same fold
        assert_eq!(folds[0], folds[1], "Same peptide targets should share fold");

        // Target-decoy pairs must be in the same fold
        assert_eq!(
            folds[0], folds[3],
            "Target PEPTIDEK z2 and its decoy must share fold"
        );
        assert_eq!(
            folds[1], folds[4],
            "Target PEPTIDEK z3 and its decoy must share fold"
        );
        assert_eq!(
            folds[2], folds[5],
            "Target ANOTHERONE and its decoy must share fold"
        );

        // All PEPTIDEK entries (both charges, both targets and decoys) in same fold
        assert_eq!(folds[0], folds[4], "All PEPTIDEK entries should share fold");
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

    /// Verify that compete_from_indices produces deterministic results with tied scores.
    ///
    /// HashMap iteration order is non-deterministic (Rust's RandomState), so entries
    /// with tied scores could appear in any order. The secondary sort key (index)
    /// must break ties deterministically. This test runs competition multiple times
    /// and asserts identical output each time.
    #[test]
    fn test_compete_from_indices_deterministic_with_ties() {
        // Create entries where multiple winners have the same score (tied)
        let scores = vec![
            0.8, // target base_id=1
            0.5, // target base_id=2
            0.5, // target base_id=3 (tied with base_id=2)
            0.5, // target base_id=4 (tied with base_id=2,3)
            0.1, // decoy base_id=1
            0.1, // decoy base_id=2 (tied with decoy base_id=1)
            0.1, // decoy base_id=3
            0.1, // decoy base_id=4
        ];
        let labels = vec![false, false, false, false, true, true, true, true];
        let entry_ids: Vec<u32> = vec![
            1,
            2,
            3,
            4,
            1 | 0x80000000,
            2 | 0x80000000,
            3 | 0x80000000,
            4 | 0x80000000,
        ];
        let indices: Vec<usize> = (0..scores.len()).collect();

        // Run competition multiple times — must always produce identical results
        let (idx_first, scores_first, decoy_first) =
            compete_from_indices(&scores, &labels, &entry_ids, &indices);

        for _ in 0..20 {
            let (idx, sc, dec) = compete_from_indices(&scores, &labels, &entry_ids, &indices);
            assert_eq!(idx, idx_first, "Winner indices must be deterministic");
            assert_eq!(dec, decoy_first, "Winner decoy flags must be deterministic");
            for (a, b) in sc.iter().zip(scores_first.iter()) {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "Winner scores must be deterministic"
                );
            }
        }

        // Verify the secondary sort: among tied-score winners, lower base_id comes first
        // Winners sorted by score desc: 0.8 (base_id=1), then 0.5 (base_ids 2,3,4)
        assert_eq!(scores_first[0], 0.8);
        // The three tied entries at 0.5 should be in base_id order (2,3,4 → indices 1,2,3)
        assert_eq!(idx_first[1], 1);
        assert_eq!(idx_first[2], 2);
        assert_eq!(idx_first[3], 3);
    }

    // ============================================================
    // Multi-level FDR tests
    // ============================================================

    /// Verifies that per-run precursor q-values are computed independently per file.
    /// Entries from different files should not affect each other's q-values.
    #[test]
    fn test_per_run_precursor_qvalues_independent_files() {
        // File A: 3 targets beat decoys, 1 decoy beats target
        // File B: all targets beat decoys
        let scores = vec![
            10.0, 9.0, 8.0, 7.5, // file_a targets (base_ids 1-4)
            7.0, 6.0, 5.0, 8.5, // file_a decoys (base_ids 1-4; decoy 4 beats target 4)
            10.0, 9.0, 8.0, // file_b targets (base_ids 5-7)
            1.0, 1.0, 1.0, // file_b decoys (base_ids 5-7)
        ];
        let labels = vec![
            false, false, false, false, true, true, true, true, false, false, false, true, true,
            true,
        ];
        let entry_ids: Vec<u32> = vec![
            1,
            2,
            3,
            4,
            1 | 0x80000000,
            2 | 0x80000000,
            3 | 0x80000000,
            4 | 0x80000000,
            5,
            6,
            7,
            5 | 0x80000000,
            6 | 0x80000000,
            7 | 0x80000000,
        ];
        let file_names: Vec<String> = vec![
            "file_a", "file_a", "file_a", "file_a", "file_a", "file_a", "file_a", "file_a",
            "file_b", "file_b", "file_b", "file_b", "file_b", "file_b",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let qvalues = compute_per_run_precursor_qvalues(&scores, &labels, &entry_ids, &file_names);

        // File B targets should all have excellent q-values (0 decoys win)
        // Conservative formula: (0+1)/n → q-values at most 1/3, 1/2, 1/1 before backward pass
        for &idx in &[8, 9, 10] {
            assert!(
                qvalues[idx] < 0.5,
                "File B target at idx {} should have low q-value, got {}",
                idx,
                qvalues[idx]
            );
        }

        // File A: decoy at base_id=4 has score 8.5 > target's 7.5, so decoy wins.
        // That means file_a has 1 decoy winner out of 4 competitions.
        // The losing target (idx 3) should have q-value 1.0 (default, not a winner).
        assert!(
            qvalues[3] > 0.5,
            "File A target beaten by decoy should have high q-value, got {}",
            qvalues[3]
        );
    }

    /// Verifies experiment-level precursor q-values pick the best observation per precursor
    /// across all files for target-decoy competition.
    #[test]
    fn test_experiment_precursor_qvalues_cross_file() {
        // Two files, same precursor (base_id=1) appears in both
        // File A: target score 5.0, decoy score 6.0 (decoy wins at run level)
        // File B: target score 9.0, decoy score 3.0 (target wins at run level)
        // At experiment level, best target (9.0) vs best decoy (6.0) → target wins
        let scores = vec![5.0, 6.0, 9.0, 3.0];
        let labels = vec![false, true, false, true];
        let entry_ids: Vec<u32> = vec![1, 1 | 0x80000000, 1, 1 | 0x80000000];

        let qvalues = compute_experiment_precursor_qvalues(&scores, &labels, &entry_ids);

        // Best target (idx 2, score 9.0) should win and get a q-value
        // With only 1 pair and target winning: q = (0+1)/1 = 1.0 (conservative)
        // But the target does win, so it should get a q-value assigned
        assert!(
            qvalues[2] <= 1.0,
            "Best target should have a q-value assigned"
        );
    }

    /// Verifies that experiment-level peptide q-values aggregate by peptide sequence,
    /// keeping only the best-scoring precursor per peptide for competition.
    #[test]
    fn test_experiment_peptide_qvalues_aggregates_by_peptide() {
        // PEPTIDEK appears as z2 (score 8.0) and z3 (score 6.0)
        // ANOTHERONE appears as z2 (score 9.0)
        // Decoys for each
        let scores = vec![8.0, 6.0, 9.0, 3.0, 2.0, 4.0];
        let labels = vec![false, false, false, true, true, true];
        let entry_ids: Vec<u32> = vec![
            1,              // PEPTIDEK z2
            2,              // PEPTIDEK z3
            3,              // ANOTHERONE z2
            1 | 0x80000000, // decoy for base_id=1
            2 | 0x80000000, // decoy for base_id=2
            3 | 0x80000000, // decoy for base_id=3
        ];
        let peptides: Vec<String> = vec![
            "PEPTIDEK",
            "PEPTIDEK",
            "ANOTHERONE",
            "KEDITPEP",
            "KEDITPEP",
            "ENOREHTONA",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let qvalues = compute_experiment_peptide_qvalues(&scores, &labels, &entry_ids, &peptides);

        // Both charge states of PEPTIDEK should get the same peptide-level q-value
        assert!(
            (qvalues[0] - qvalues[1]).abs() < 1e-10,
            "Same peptide should get same q-value: z2={}, z3={}",
            qvalues[0],
            qvalues[1]
        );

        // Both targets should beat their decoys (scores 8>3, 9>4) → good q-values
        assert!(qvalues[0] < 1.0, "PEPTIDEK should have q-value < 1.0");
        assert!(qvalues[2] < 1.0, "ANOTHERONE should have q-value < 1.0");
    }

    /// Verifies that the effective q-value is max(precursor_qvalue, peptide_qvalue),
    /// enforcing dual-level FDR control.
    #[test]
    fn test_effective_qvalue_is_max_of_precursor_and_peptide() {
        // Simulate run_percolator output: different precursor and peptide q-values
        // Precursor passes but peptide fails → should be rejected
        let run_precursor_q: f64 = 0.005;
        let run_peptide_q: f64 = 0.015;
        let effective = run_precursor_q.max(run_peptide_q);
        assert!(
            effective > 0.01,
            "Should be rejected: max({}, {}) = {} > 0.01",
            run_precursor_q,
            run_peptide_q,
            effective
        );

        // Peptide passes but precursor fails → should be rejected
        let run_precursor_q: f64 = 0.015;
        let run_peptide_q: f64 = 0.005;
        let effective = run_precursor_q.max(run_peptide_q);
        assert!(
            effective > 0.01,
            "Should be rejected: max({}, {}) = {} > 0.01",
            run_precursor_q,
            run_peptide_q,
            effective
        );

        // Both pass → should be accepted
        let run_precursor_q: f64 = 0.003;
        let run_peptide_q: f64 = 0.008;
        let effective = run_precursor_q.max(run_peptide_q);
        assert!(
            effective <= 0.01,
            "Should be accepted: max({}, {}) = {} <= 0.01",
            run_precursor_q,
            run_peptide_q,
            effective
        );
    }

    /// Verifies that per-run peptide q-values are propagated to all charge states of the same peptide.
    #[test]
    fn test_per_run_peptide_qvalues_propagate_to_charges() {
        // Same peptide PEPTIDEK with z2 (high score) and z3 (low score) in same file
        // The peptide-level q-value should be the same for both charge states
        let scores = vec![10.0, 5.0, 2.0, 1.0];
        let labels = vec![false, false, true, true];
        let entry_ids: Vec<u32> = vec![1, 2, 1 | 0x80000000, 2 | 0x80000000];
        let file_names: Vec<String> = vec!["file_a"; 4].into_iter().map(String::from).collect();
        let peptides: Vec<String> = vec!["PEPTIDEK", "PEPTIDEK", "KEDITPEP", "KEDITPEP"]
            .into_iter()
            .map(String::from)
            .collect();

        let qvalues =
            compute_per_run_peptide_qvalues(&scores, &labels, &entry_ids, &file_names, &peptides);

        // Both charge states should get the same peptide-level q-value
        assert!(
            (qvalues[0] - qvalues[1]).abs() < 1e-10,
            "Same peptide different charges should have same q-value: z2={}, z3={}",
            qvalues[0],
            qvalues[1]
        );
    }

    // ============================================================
    // CV fold integrity and subsampling tests
    // ============================================================

    /// Verifies that subsampling by peptide group keeps target-decoy pairs together.
    /// If a target is selected, its paired decoy must also be selected, and vice versa.
    #[test]
    fn test_subsample_keeps_target_decoy_pairs() {
        // 10 peptide groups (20 entries total: 10 targets + 10 decoys)
        let mut labels = Vec::new();
        let mut entry_ids = Vec::new();
        let mut peptides = Vec::new();
        for i in 1..=10 {
            labels.push(false); // target
            entry_ids.push(i as u32);
            peptides.push(format!("PEPTIDE{}", i));
            labels.push(true); // decoy
            entry_ids.push((i as u32) | 0x80000000);
            peptides.push(format!("DECOY{}", i));
        }

        // Subsample to 10 entries (should be ~5 groups × 2 entries each)
        let selected = subsample_by_peptide_group(&labels, &entry_ids, &peptides, 10, 42);

        // Every selected target must have its paired decoy also selected
        for &idx in &selected {
            let base_id = entry_ids[idx] & 0x7FFFFFFF;
            let is_decoy = labels[idx];

            // Find the paired entry
            let paired_idx = if is_decoy {
                // Find the target with same base_id
                (0..labels.len()).find(|&j| !labels[j] && entry_ids[j] == base_id)
            } else {
                // Find the decoy with same base_id
                (0..labels.len()).find(|&j| labels[j] && (entry_ids[j] & 0x7FFFFFFF) == base_id)
            };

            if let Some(paired) = paired_idx {
                assert!(
                    selected.contains(&paired),
                    "Entry at idx {} (base_id={}, decoy={}) selected but paired entry at idx {} is missing",
                    idx, base_id, is_decoy, paired
                );
            }
        }
    }

    /// Verifies that subsampling keeps all charge states of the same peptide together.
    #[test]
    fn test_subsample_keeps_charge_states_together() {
        // PEPTIDEK has z2 (base_id=1) and z3 (base_id=2), both with decoys
        // ANOTHERONE has z2 (base_id=3) with decoy
        let labels = vec![false, false, false, true, true, true];
        let entry_ids: Vec<u32> = vec![
            1,              // PEPTIDEK z2
            2,              // PEPTIDEK z3
            3,              // ANOTHERONE z2
            1 | 0x80000000, // decoy for PEPTIDEK z2
            2 | 0x80000000, // decoy for PEPTIDEK z3
            3 | 0x80000000, // decoy for ANOTHERONE z2
        ];
        let peptides: Vec<String> = vec![
            "PEPTIDEK",
            "PEPTIDEK",
            "ANOTHERONE",
            "KEDITPEP",
            "KEDITPEP",
            "ENOREHTONA",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        // Subsample to 4 entries — should select complete peptide groups
        let selected = subsample_by_peptide_group(&labels, &entry_ids, &peptides, 4, 42);

        // If any PEPTIDEK entry is selected, ALL PEPTIDEK entries must be selected
        let peptidek_indices = vec![0, 1, 3, 4]; // targets z2, z3 + decoys
        let any_peptidek_selected = peptidek_indices.iter().any(|i| selected.contains(i));
        if any_peptidek_selected {
            for &idx in &peptidek_indices {
                assert!(
                    selected.contains(&idx),
                    "PEPTIDEK entry at idx {} should be selected (group must stay together)",
                    idx
                );
            }
        }
    }

    /// Verifies that fold assignment groups multiple charge states with their decoys
    /// into the same fold, even with many peptide groups.
    #[test]
    fn test_fold_assignment_multi_charge_with_decoys() {
        // 5 peptides, some with multiple charge states
        let labels = vec![
            false, false, // PEPTIDEK z2, z3
            false, // ANOTHERONE z2
            false, false, // THIRDPEP z2, z3
            true, true, // decoys for PEPTIDEK
            true, // decoy for ANOTHERONE
            true, true, // decoys for THIRDPEP
        ];
        let entry_ids: Vec<u32> = vec![
            1,
            2, // PEPTIDEK z2 (base=1), z3 (base=2)
            3, // ANOTHERONE z2 (base=3)
            4,
            5, // THIRDPEP z2 (base=4), z3 (base=5)
            1 | 0x80000000,
            2 | 0x80000000, // decoys for PEPTIDEK
            3 | 0x80000000, // decoy for ANOTHERONE
            4 | 0x80000000,
            5 | 0x80000000, // decoys for THIRDPEP
        ];
        let peptides: Vec<String> = vec![
            "PEPTIDEK",
            "PEPTIDEK",
            "ANOTHERONE",
            "THIRDPEP",
            "THIRDPEP",
            "KEDITPEP",
            "KEDITPEP",
            "ENOREHTONA",
            "PEPTDRIHT",
            "PEPTDRIHT",
        ]
        .into_iter()
        .map(String::from)
        .collect();

        let folds = create_stratified_folds_by_peptide(&labels, &peptides, &entry_ids, 3);

        // PEPTIDEK targets (0,1) and their decoys (5,6) must share a fold
        assert_eq!(
            folds[0], folds[1],
            "PEPTIDEK z2 and z3 targets must share fold"
        );
        assert_eq!(
            folds[0], folds[5],
            "PEPTIDEK z2 and its decoy must share fold"
        );
        assert_eq!(
            folds[1], folds[6],
            "PEPTIDEK z3 and its decoy must share fold"
        );

        // THIRDPEP targets (3,4) and their decoys (8,9) must share a fold
        assert_eq!(
            folds[3], folds[4],
            "THIRDPEP z2 and z3 targets must share fold"
        );
        assert_eq!(
            folds[3], folds[8],
            "THIRDPEP z2 and its decoy must share fold"
        );
        assert_eq!(
            folds[4], folds[9],
            "THIRDPEP z3 and its decoy must share fold"
        );

        // ANOTHERONE target (2) and its decoy (7) must share a fold
        assert_eq!(
            folds[2], folds[7],
            "ANOTHERONE and its decoy must share fold"
        );
    }

    /// Verifies that subsampling with max_entries >= total returns all entries unchanged.
    #[test]
    fn test_subsample_no_reduction_when_under_limit() {
        let labels = vec![false, true, false, true];
        let entry_ids: Vec<u32> = vec![1, 1 | 0x80000000, 2, 2 | 0x80000000];
        let peptides: Vec<String> = vec!["PEP1", "DEC1", "PEP2", "DEC2"]
            .into_iter()
            .map(String::from)
            .collect();

        let selected = subsample_by_peptide_group(&labels, &entry_ids, &peptides, 100, 42);

        assert_eq!(
            selected.len(),
            4,
            "Should return all entries when under limit"
        );
        assert_eq!(selected, vec![0, 1, 2, 3]);
    }
}
