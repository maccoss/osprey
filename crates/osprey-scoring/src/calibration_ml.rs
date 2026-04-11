//! Machine learning-based scoring for calibration matches
//!
//! Uses Linear Discriminant Analysis (LDA) to learn optimal feature weights
//! from target-decoy competition, replacing simple correlation-based scoring.

use crate::batch::CalibrationMatch;
use osprey_core::OspreyError;
use osprey_ml::{kde, linear_discriminant::LinearDiscriminantAnalysis, matrix::Matrix, qvalue};
use std::collections::HashMap;

/// Train LDA on calibration matches and score them using 3-fold cross-validation
///
/// # Arguments
/// * `matches` - Calibration matches with extracted features
/// * `use_isotope_feature` - Whether to include isotope_cosine as 5th feature (HRAM mode)
///
/// # Returns
/// Number of matches passing 1% FDR threshold
pub fn train_and_score_calibration(
    matches: &mut [CalibrationMatch],
    use_isotope_feature: bool,
) -> Result<usize, OspreyError> {
    if matches.is_empty() {
        return Ok(0);
    }

    let n_features = if use_isotope_feature { 5 } else { 4 };
    log::debug!(
        "Training LDA with 3-fold cross-validation on {} calibration matches ({} features)",
        matches.len(),
        n_features
    );

    // 1. Extract feature matrix
    let features = extract_feature_matrix(matches, use_isotope_feature);

    // 2. Build target/decoy labels (true = decoy)
    let decoy_labels: Vec<bool> = matches.iter().map(|m| m.is_decoy).collect();

    // 3. Extract entry IDs for paired target-decoy competition
    let entry_ids: Vec<u32> = matches.iter().map(|m| m.entry_id).collect();

    // 4. Extract sequences for fold grouping (keep charge states together)
    let sequences: Vec<String> = matches.iter().map(|m| m.sequence.clone()).collect();

    // 5. Train LDA with cross-validation and non-negative constraints
    let discriminants = train_lda_with_nonnegative_cv(
        &features,
        &decoy_labels,
        &entry_ids,
        &sequences,
        use_isotope_feature,
    )?;

    // 5. KDE for posterior error probabilities
    let kde_model = kde::Builder::default().build(&discriminants, &decoy_labels);

    // 6. Assign discriminant scores and posterior errors
    for (i, m) in matches.iter_mut().enumerate() {
        m.discriminant_score = discriminants[i];
        m.posterior_error = kde_model.posterior_error(discriminants[i]);
    }

    // 7. Calculate q-values using paired target-decoy competition
    //    IMPORTANT: Must happen BEFORE sorting matches, since entry_ids and
    //    decoy_labels are in the original order (aligned with discriminants).
    let all_indices: Vec<usize> = (0..matches.len()).collect();
    let winner_indices =
        compete_calibration_pairs(&discriminants, &entry_ids, &decoy_labels, &all_indices);

    let winner_is_decoy: Vec<bool> = winner_indices.iter().map(|&i| decoy_labels[i]).collect();
    let mut winner_q = vec![0.0; winner_indices.len()];
    let n_passing = qvalue::calculate_q_values(&winner_is_decoy, &mut winner_q);

    // Map q-values back: winners get their q-value, losers get 1.0
    let mut q_values = vec![1.0; matches.len()];
    for (rank, &idx) in winner_indices.iter().enumerate() {
        q_values[idx] = winner_q[rank];
    }

    // 8. Assign q-values back to matches
    for (i, m) in matches.iter_mut().enumerate() {
        m.q_value = q_values[i];
    }

    // 9. Sort by discriminant score descending (best matches first)
    matches.sort_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));

    // DEBUG: Log discriminant score distributions
    let target_scores: Vec<f64> = matches
        .iter()
        .filter(|m| !m.is_decoy)
        .map(|m| m.discriminant_score)
        .collect();
    let decoy_scores: Vec<f64> = matches
        .iter()
        .filter(|m| m.is_decoy)
        .map(|m| m.discriminant_score)
        .collect();

    if !target_scores.is_empty() && !decoy_scores.is_empty() {
        let target_mean = target_scores.iter().sum::<f64>() / target_scores.len() as f64;
        let decoy_mean = decoy_scores.iter().sum::<f64>() / decoy_scores.len() as f64;
        let target_max = target_scores
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let target_min = target_scores.iter().cloned().fold(f64::INFINITY, f64::min);
        let decoy_max = decoy_scores
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let decoy_min = decoy_scores.iter().cloned().fold(f64::INFINITY, f64::min);

        log::debug!(
            "  Discriminant scores: target_mean={:.3}, decoy_mean={:.3}",
            target_mean,
            decoy_mean
        );
        log::debug!(
            "  Discriminant ranges: targets=[{:.3}, {:.3}], decoys=[{:.3}, {:.3}]",
            target_min,
            target_max,
            decoy_min,
            decoy_max
        );
    }

    log::info!("LDA scoring complete: {} matches passing 1% FDR", n_passing);

    // Log median feature values for diagnostics
    log_median_features(matches, use_isotope_feature);

    Ok(n_passing)
}

/// Train single LDA model on all data with non-negative weight constraints
///
/// DEPRECATED: Use train_lda_with_nonnegative_cv instead for proper cross-validation
///
/// # Arguments
/// * `features` - Full feature matrix
/// * `decoy_labels` - Target/decoy labels (true = decoy)
/// * `use_isotope_feature` - Whether isotope feature is included
///
/// # Returns
/// Discriminant scores for all samples
#[allow(dead_code)]
fn train_lda_single_model(
    features: &Matrix,
    decoy_labels: &[bool],
    use_isotope_feature: bool,
) -> Result<Vec<f64>, OspreyError> {
    // Train LDA on all data
    let lda = LinearDiscriminantAnalysis::fit(features, decoy_labels)
        .ok_or_else(|| OspreyError::config("LDA training failed"))?;

    // Get eigenvector and clip negative weights to zero
    let mut weights = lda.eigenvector().to_vec();

    // Count negative weights before clipping
    let n_negative = weights.iter().filter(|&w| *w < 0.0).count();
    if n_negative > 0 {
        log::warn!("  Clipping {} negative feature weights to zero", n_negative);
    }

    // Clip to non-negative
    for w in weights.iter_mut() {
        if *w < 0.0 {
            *w = 0.0;
        }
    }

    // Renormalize to unit length
    let norm: f64 = weights.iter().map(|w| w * w).sum::<f64>().sqrt();
    if norm > 1e-10 {
        for w in weights.iter_mut() {
            *w /= norm;
        }
    } else {
        // If all weights were negative and got clipped to zero, fall back to equal weights
        log::warn!("  All weights were negative - using equal weights");
        let equal_weight = 1.0 / (weights.len() as f64).sqrt();
        for w in weights.iter_mut() {
            *w = equal_weight;
        }
    }

    // Log final weights
    let feature_names = if use_isotope_feature {
        vec!["correlation", "libcosine", "top6", "xcorr", "isotope"]
    } else {
        vec!["correlation", "libcosine", "top6", "xcorr"]
    };

    log::debug!(
        "  LDA feature weights (non-negative): {}",
        feature_names
            .iter()
            .enumerate()
            .map(|(i, name)| format!("{}={:.3}", name, weights[i]))
            .collect::<Vec<_>>()
            .join(", ")
    );

    // Create new LDA with non-negative weights
    let lda_nonneg =
        LinearDiscriminantAnalysis::from_weights(weights).map_err(OspreyError::config)?;

    // Score all samples
    let discriminants = lda_nonneg.predict(features);

    Ok(discriminants)
}

/// Train LDA with target selection and 3-fold cross-validation
///
/// Uses Percolator-style methodology with best-iteration tracking:
/// 1. Evaluate all single features, pick the best as baseline
/// 2. Select high-confidence targets (q < train_fdr) as positive training set
/// 3. Train LDA on selected targets vs ALL decoys using 3-fold CV
/// 4. Score all data with consensus weights
/// 5. Keep the best result across iterations (revert if LDA degrades)
/// 6. Stop early if 2 consecutive iterations fail to improve
///
/// Key insight: iterating with updated target selection causes selection bias
/// (the LDA drifts weight away from the initial best feature). We cap iterations
/// and always revert to the best result seen.
///
/// # Arguments
/// * `features` - Full feature matrix (normalized)
/// * `decoy_labels` - Target/decoy labels (true = decoy)
/// * `sequences` - Peptide sequences (for grouping charge states in folds)
/// * `_use_isotope_feature` - Reserved for future use
///
/// # Returns
/// Discriminant scores for all samples (from the best iteration)
fn train_lda_with_nonnegative_cv(
    features: &Matrix,
    decoy_labels: &[bool],
    entry_ids: &[u32],
    sequences: &[String],
    _use_isotope_feature: bool,
) -> Result<Vec<f64>, OspreyError> {
    const N_FOLDS: usize = 3;
    const MAX_ITERATIONS: usize = 3;
    let mut train_fdr: f64 = 0.01;
    const MIN_POSITIVE_EXAMPLES: usize = 50;

    let n_samples = features.rows;
    let n_features = features.cols;
    let feature_names = ["correlation", "libcosine", "top6", "xcorr"];

    // Create fold assignments (stable across iterations - Percolator-RESET grouping)
    let fold_assignments = create_stratified_folds_by_peptide(decoy_labels, sequences, N_FOLDS);

    // Find the best single feature to use as baseline (like Percolator)
    let mut best_feat_idx = 0;
    let mut best_feat_passing = 0;
    for feat_idx in 0..n_features {
        let feat_scores: Vec<f64> = (0..n_samples).map(|i| features[(i, feat_idx)]).collect();
        let n_pass = count_passing_targets(&feat_scores, decoy_labels, entry_ids, train_fdr);
        log::debug!(
            "  Initial feature '{}': {} pass {:.0}% FDR",
            feature_names[feat_idx],
            n_pass,
            train_fdr * 100.0
        );
        if n_pass > best_feat_passing {
            best_feat_passing = n_pass;
            best_feat_idx = feat_idx;
        }
    }

    // If no targets pass at configured FDR, loosen to 5% so training can proceed
    if best_feat_passing == 0 {
        let relaxed_fdr = 0.05;
        let mut relaxed_best_idx = 0;
        let mut relaxed_best_passing = 0;
        for feat_idx in 0..n_features {
            let feat_scores: Vec<f64> = (0..n_samples).map(|i| features[(i, feat_idx)]).collect();
            let n_pass = count_passing_targets(&feat_scores, decoy_labels, entry_ids, relaxed_fdr);
            if n_pass > relaxed_best_passing {
                relaxed_best_passing = n_pass;
                relaxed_best_idx = feat_idx;
            }
        }
        if relaxed_best_passing > 0 {
            log::warn!(
                "  No targets at {:.0}% FDR — loosening train FDR to {:.0}%",
                train_fdr * 100.0,
                relaxed_fdr * 100.0
            );
            train_fdr = relaxed_fdr;
            best_feat_idx = relaxed_best_idx;
            best_feat_passing = relaxed_best_passing;
        } else {
            log::warn!(
                "  No targets at {:.0}% or {:.0}% FDR — features cannot discriminate targets from decoys",
                train_fdr * 100.0,
                relaxed_fdr * 100.0
            );
        }
    }

    // Initialize with best single feature as both current and best-so-far
    let baseline_scores: Vec<f64> = (0..n_samples)
        .map(|i| features[(i, best_feat_idx)])
        .collect();

    log::debug!(
        "  Baseline: '{}' = {} pass {:.0}% FDR",
        feature_names[best_feat_idx],
        best_feat_passing,
        train_fdr * 100.0
    );

    // Best-so-far tracking: always revert to the iteration with the most passing targets
    let mut best_scores = baseline_scores.clone();
    let mut best_passing = best_feat_passing;
    let mut best_iteration = 0; // 0 = baseline

    // Current scores used for target selection in next iteration
    let mut current_scores = baseline_scores;

    // Track history for logging
    let mut iteration_passing: Vec<usize> = vec![best_feat_passing];
    let mut consecutive_no_improve = 0;

    for iteration in 0..MAX_ITERATIONS {
        let mut fold_weights: Vec<Vec<f64>> = Vec::new();
        let mut total_selected_targets = 0;

        // Phase 1: Use CV to estimate stable weights
        for fold_idx in 0..N_FOLDS {
            let train_indices: Vec<usize> = (0..n_samples)
                .filter(|&i| fold_assignments[i] != fold_idx)
                .collect();

            // Select high-confidence targets from TRAIN set using current scores
            let selected_target_indices = select_positive_training_set(
                &current_scores,
                decoy_labels,
                entry_ids,
                &train_indices,
                train_fdr,
                MIN_POSITIVE_EXAMPLES,
            );

            total_selected_targets += selected_target_indices.len();

            // Collect ALL decoy indices from train set
            let decoy_indices: Vec<usize> = train_indices
                .iter()
                .filter(|&&i| decoy_labels[i])
                .copied()
                .collect();

            // Build training set: selected targets + all decoys
            let mut training_indices = selected_target_indices.clone();
            training_indices.extend_from_slice(&decoy_indices);

            let train_features = extract_rows(features, &training_indices);
            let train_labels: Vec<bool> =
                training_indices.iter().map(|&i| decoy_labels[i]).collect();

            log::debug!(
                "    Fold {}/{}: {} selected targets + {} decoys = {} training samples",
                fold_idx + 1,
                N_FOLDS,
                selected_target_indices.len(),
                decoy_indices.len(),
                training_indices.len()
            );

            // Train LDA on clean training set
            // If LDA fails (singular scatter matrix), skip this fold
            match LinearDiscriminantAnalysis::fit(&train_features, &train_labels) {
                Some(lda) => {
                    let weights = lda.eigenvector().to_vec();

                    log::debug!(
                        "    Fold {}/{} weights: {}",
                        fold_idx + 1,
                        N_FOLDS,
                        feature_names
                            .iter()
                            .enumerate()
                            .map(|(i, name)| format!("{}={:.3}", name, weights[i]))
                            .collect::<Vec<_>>()
                            .join(", ")
                    );

                    fold_weights.push(weights);
                }
                None => {
                    log::warn!(
                        "    Fold {}/{}: LDA fit failed (singular matrix, {} targets, {} decoys), skipping",
                        fold_idx + 1,
                        N_FOLDS,
                        selected_target_indices.len(),
                        decoy_indices.len()
                    );
                }
            }
        }

        // Phase 2: Average weights across folds → consensus weights
        // If all folds failed, skip this iteration entirely
        if fold_weights.is_empty() {
            log::warn!(
                "  Iteration {}: All folds failed, using best single feature as baseline",
                iteration + 1
            );
            break;
        }
        if fold_weights.len() < N_FOLDS {
            log::warn!(
                "  Iteration {}: {}/{} folds succeeded, averaging available weights",
                iteration + 1,
                fold_weights.len(),
                N_FOLDS
            );
        }
        let mut consensus_weights = average_weights(&fold_weights);

        // Clip negative weights to zero (non-negative constraint)
        let n_negative = consensus_weights.iter().filter(|&&w| w < 0.0).count();
        if n_negative > 0 {
            log::debug!(
                "    Clipping {} negative consensus weights to zero",
                n_negative
            );
            for w in consensus_weights.iter_mut() {
                if *w < 0.0 {
                    *w = 0.0;
                }
            }
            // Renormalize
            let norm: f64 = consensus_weights.iter().map(|w| w * w).sum::<f64>().sqrt();
            if norm > 1e-10 {
                for w in consensus_weights.iter_mut() {
                    *w /= norm;
                }
            }
        }

        // Phase 3: Score ALL data with consensus weights
        let lda_consensus = LinearDiscriminantAnalysis::from_weights(consensus_weights.clone())
            .map_err(OspreyError::config)?;
        let new_scores = lda_consensus.predict(features);

        let n_passing = count_passing_targets(&new_scores, decoy_labels, entry_ids, train_fdr);

        log::debug!(
            "  Iteration {}: weights=[{}], selected {} train targets, {} pass {:.0}% FDR",
            iteration + 1,
            feature_names
                .iter()
                .enumerate()
                .map(|(i, name)| format!("{}={:.3}", name, consensus_weights[i]))
                .collect::<Vec<_>>()
                .join(", "),
            total_selected_targets / N_FOLDS,
            n_passing,
            train_fdr * 100.0,
        );

        iteration_passing.push(n_passing);

        // Track best iteration: keep scores from whichever iteration gave the most passing
        if n_passing > best_passing {
            best_scores = new_scores.clone();
            best_passing = n_passing;
            best_iteration = iteration + 1;
            consecutive_no_improve = 0;
            log::debug!(
                "    -> New best: {} pass {:.0}% FDR (iteration {})",
                best_passing,
                train_fdr * 100.0,
                best_iteration
            );
        } else {
            consecutive_no_improve += 1;
            log::debug!(
                "    -> No improvement (best remains {} from iteration {}, {} consecutive non-improvements)",
                best_passing, best_iteration, consecutive_no_improve
            );
        }

        // Update current scores for next iteration's target selection
        current_scores = new_scores;

        // Stop early if 2 consecutive iterations didn't improve
        if consecutive_no_improve >= 2 {
            log::debug!(
                "  Stopping early: {} consecutive iterations without improvement",
                consecutive_no_improve
            );
            break;
        }
    }

    // Log final result
    log::debug!(
        "  Iteration history (pass {:.0}% FDR): {:?}",
        train_fdr * 100.0,
        iteration_passing
    );
    log::debug!(
        "  Using iteration {} ({} pass {:.0}% FDR)",
        best_iteration,
        best_passing,
        train_fdr * 100.0
    );

    Ok(best_scores)
}

/// Select high-confidence targets from a subset of data for positive training set
///
/// This is the key Percolator/Mokapot innovation: don't train on all targets (noisy),
/// train only on targets that pass an FDR threshold (clean positive examples).
/// Uses paired target-decoy competition within the subset.
///
/// # Arguments
/// * `scores` - Current discriminant scores for ALL samples
/// * `decoy_labels` - Target/decoy labels for ALL samples
/// * `entry_ids` - Entry IDs for ALL samples (for pairing via `id & 0x7FFFFFFF`)
/// * `subset_indices` - Indices of samples to consider (train set)
/// * `fdr_threshold` - Q-value threshold for selecting targets (e.g., 0.01)
/// * `min_targets` - Minimum number of targets to select
///
/// # Returns
/// Indices (into the original arrays) of selected high-confidence targets
fn select_positive_training_set(
    scores: &[f64],
    decoy_labels: &[bool],
    entry_ids: &[u32],
    subset_indices: &[usize],
    fdr_threshold: f64,
    min_targets: usize,
) -> Vec<usize> {
    // Paired competition within subset: only winners enter FDR walk
    let winner_indices = compete_calibration_pairs(scores, entry_ids, decoy_labels, subset_indices);

    // Compute q-values on winners (already sorted by score descending)
    let winner_is_decoy: Vec<bool> = winner_indices.iter().map(|&i| decoy_labels[i]).collect();
    let mut q_values = vec![1.0; winner_indices.len()];
    if !winner_indices.is_empty() {
        qvalue::calculate_q_values(&winner_is_decoy, &mut q_values);
    }

    // Select targets with q ≤ threshold
    let select_at_threshold = |threshold: f64| -> Vec<usize> {
        winner_indices
            .iter()
            .enumerate()
            .filter(|&(rank, &idx)| !decoy_labels[idx] && q_values[rank] <= threshold)
            .map(|(_, &idx)| idx)
            .collect()
    };

    let mut selected = select_at_threshold(fdr_threshold);

    // If too few, relax threshold progressively
    if selected.len() < min_targets {
        let relaxed_thresholds = [0.05, 0.10, 0.25, 0.50];
        for &threshold in &relaxed_thresholds {
            selected = select_at_threshold(threshold);
            if selected.len() >= min_targets {
                log::debug!(
                    "      Relaxed training FDR to {:.0}% to get {} targets",
                    threshold * 100.0,
                    selected.len()
                );
                break;
            }
        }
    }

    selected
}

/// Paired target-decoy competition for calibration matches.
///
/// Groups matches by `entry_id & 0x7FFFFFFF`, competes each pair on the given
/// score, and returns only winner indices (sorted by score descending).
/// Singletons auto-win. Ties go to decoy (conservative for FDR estimation).
fn compete_calibration_pairs(
    scores: &[f64],
    entry_ids: &[u32],
    is_decoy: &[bool],
    subset_indices: &[usize],
) -> Vec<usize> {
    // Group indices by base_id = entry_id & 0x7FFFFFFF
    let mut groups: HashMap<u32, (Option<usize>, Option<usize>)> = HashMap::new();
    for &idx in subset_indices {
        let base_id = entry_ids[idx] & 0x7FFFFFFF;
        let slot = groups.entry(base_id).or_insert((None, None));
        if is_decoy[idx] {
            // Keep the decoy with the higher score if multiple per base_id
            if slot.1.map_or(true, |prev| scores[idx] > scores[prev]) {
                slot.1 = Some(idx);
            }
        } else if slot.0.map_or(true, |prev| scores[idx] > scores[prev]) {
            slot.0 = Some(idx);
        }
    }

    // Compete each pair
    let mut winners: Vec<usize> = Vec::with_capacity(groups.len());
    for (target_opt, decoy_opt) in groups.values() {
        match (*target_opt, *decoy_opt) {
            (Some(t), None) => winners.push(t),
            (None, Some(d)) => winners.push(d),
            (Some(t), Some(d)) => {
                if scores[t] > scores[d] {
                    winners.push(t); // strict >: ties go to decoy
                } else {
                    winners.push(d);
                }
            }
            (None, None) => {}
        }
    }

    // Sort winners by score descending, then base_id ascending for deterministic tiebreaking.
    // IMPORTANT: Use base_id (not array index) as tiebreaker. Array indices depend on input
    // order (HashMap iteration), and if the input is sorted by entry_id, all targets get low
    // indices → systematic target-before-decoy bias → artificially inflated FDR estimates.
    // base_id is intrinsic to the entry and doesn't correlate with target/decoy status.
    winners.sort_by(|&a, &b| {
        scores[b]
            .total_cmp(&scores[a])
            .then((entry_ids[a] & 0x7FFFFFFF).cmp(&(entry_ids[b] & 0x7FFFFFFF)))
    });
    winners
}

/// Count targets passing a given FDR threshold using paired competition
fn count_passing_targets(
    scores: &[f64],
    decoy_labels: &[bool],
    entry_ids: &[u32],
    fdr_threshold: f64,
) -> usize {
    let all_indices: Vec<usize> = (0..scores.len()).collect();

    // Paired competition: only winners enter FDR walk
    let winner_indices = compete_calibration_pairs(scores, entry_ids, decoy_labels, &all_indices);

    // Compute q-values on winners (already sorted by score descending)
    let winner_is_decoy: Vec<bool> = winner_indices.iter().map(|&i| decoy_labels[i]).collect();
    let mut winner_q = vec![0.0; winner_indices.len()];
    qvalue::calculate_q_values(&winner_is_decoy, &mut winner_q);

    // Count passing targets
    let mut n_passing = 0;
    for (rank, &idx) in winner_indices.iter().enumerate() {
        if !decoy_labels[idx] && winner_q[rank] <= fdr_threshold {
            n_passing += 1;
        }
    }

    n_passing
}

/// Average feature weights across multiple models
fn average_weights(weights_list: &[Vec<f64>]) -> Vec<f64> {
    if weights_list.is_empty() {
        return Vec::new();
    }

    let n_features = weights_list[0].len();
    let n_models = weights_list.len();

    let mut avg = vec![0.0; n_features];
    for weights in weights_list {
        for (i, &w) in weights.iter().enumerate() {
            avg[i] += w;
        }
    }

    for w in avg.iter_mut() {
        *w /= n_models as f64;
    }

    // Renormalize to unit length
    let norm: f64 = avg.iter().map(|w| w * w).sum::<f64>().sqrt();
    if norm > 1e-10 {
        for w in avg.iter_mut() {
            *w /= norm;
        }
    }

    avg
}

/// Create stratified fold assignments for cross-validation
///
/// Implements Percolator-RESET approach:
/// - Groups PSMs by peptide sequence (keeps charge states together)
/// - Stratified by target/decoy to maintain class balance
/// - Deterministic assignment for stability across iterations
fn create_stratified_folds_by_peptide(
    labels: &[bool],
    sequences: &[String],
    n_folds: usize,
) -> Vec<usize> {
    use std::collections::HashMap;

    let mut fold_assignments = vec![0; labels.len()];

    // Group indices by peptide sequence and target/decoy status
    let mut target_peptides: HashMap<String, Vec<usize>> = HashMap::new();
    let mut decoy_peptides: HashMap<String, Vec<usize>> = HashMap::new();

    for (i, (seq, &is_decoy)) in sequences.iter().zip(labels).enumerate() {
        if is_decoy {
            decoy_peptides.entry(seq.clone()).or_default().push(i);
        } else {
            target_peptides.entry(seq.clone()).or_default().push(i);
        }
    }

    // Convert to sorted vectors for consistent ordering
    let mut target_groups: Vec<(String, Vec<usize>)> = target_peptides.into_iter().collect();
    let mut decoy_groups: Vec<(String, Vec<usize>)> = decoy_peptides.into_iter().collect();

    target_groups.sort_by(|a, b| a.0.cmp(&b.0));
    decoy_groups.sort_by(|a, b| a.0.cmp(&b.0));

    // Assign peptide groups to folds round-robin
    for (i, (_seq, indices)) in target_groups.iter().enumerate() {
        let fold = i % n_folds;
        for &idx in indices {
            fold_assignments[idx] = fold;
        }
    }

    for (i, (_seq, indices)) in decoy_groups.iter().enumerate() {
        let fold = i % n_folds;
        for &idx in indices {
            fold_assignments[idx] = fold;
        }
    }

    log::debug!(
        "  Created {} folds: {} target peptides ({} PSMs), {} decoy peptides ({} PSMs)",
        n_folds,
        target_groups.len(),
        target_groups.iter().map(|(_, v)| v.len()).sum::<usize>(),
        decoy_groups.len(),
        decoy_groups.iter().map(|(_, v)| v.len()).sum::<usize>()
    );

    fold_assignments
}

/// Extract specific rows from a matrix
fn extract_rows(matrix: &Matrix, row_indices: &[usize]) -> Matrix {
    let n_cols = matrix.cols;
    let data: Vec<f64> = row_indices
        .iter()
        .flat_map(|&row| (0..n_cols).map(move |col| matrix[(row, col)]))
        .collect();

    Matrix::new(data, row_indices.len(), n_cols)
}

/// Extract feature matrix from calibration matches
///
/// Features are normalized to similar ranges for fair weighting:
/// - correlation: 0-1 (typical range 0-6, divided by 6)
/// - libcosine_apex: 0-1 (already normalized)
/// - top6_matched_apex: 0-1 (count 0-6, divided by 6)
/// - xcorr: ~0-1 (typical range 0-3, divided by 3)
fn extract_feature_matrix(matches: &[CalibrationMatch], _use_isotope_feature: bool) -> Matrix {
    let n_features = 4;

    let features: Vec<f64> = matches
        .iter()
        .flat_map(|m| {
            vec![
                // Normalize correlation score (typical range 0-6)
                (m.correlation_score / 6.0).clamp(0.0, 1.0),
                // LibCosine already 0-1
                m.libcosine_apex.clamp(0.0, 1.0),
                // Top-6 count normalized
                (m.top6_matched_apex as f64 / 6.0).clamp(0.0, 1.0),
                // XCorr normalized (typical range 0-3 for unit resolution, 0-1 for HRAM)
                (m.xcorr_score / 3.0).clamp(0.0, 1.0),
            ]
        })
        .collect();

    Matrix::new(features, matches.len(), n_features)
}

/// Log median feature values for diagnostics
fn log_median_features(matches: &[CalibrationMatch], use_isotope: bool) {
    if matches.is_empty() {
        return;
    }

    let mut corrs: Vec<f64> = matches.iter().map(|m| m.correlation_score).collect();
    let mut libcos: Vec<f64> = matches.iter().map(|m| m.libcosine_apex).collect();
    let mut top6s: Vec<f64> = matches.iter().map(|m| m.top6_matched_apex as f64).collect();
    let mut xcorrs: Vec<f64> = matches.iter().map(|m| m.xcorr_score).collect();
    let mut discrims: Vec<f64> = matches.iter().map(|m| m.discriminant_score).collect();

    corrs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    libcos.sort_by(|a, b| a.partial_cmp(b).unwrap());
    top6s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    xcorrs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    discrims.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mid = matches.len() / 2;

    log::debug!(
        "  Median features: corr={:.3}, libcos={:.3}, top6={:.1}, xcorr={:.4}, discrim={:.3}",
        corrs[mid],
        libcos[mid],
        top6s[mid],
        xcorrs[mid],
        discrims[mid]
    );

    if use_isotope {
        let mut isos: Vec<f64> = matches
            .iter()
            .filter_map(|m| m.isotope_cosine_score)
            .collect();
        if !isos.is_empty() {
            isos.sort_by(|a, b| a.partial_cmp(b).unwrap());
            log::debug!("  Median isotope_cosine: {:.3}", isos[isos.len() / 2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_match(
        entry_id: u32,
        correlation: f64,
        libcosine: f64,
        top6: u8,
        hyperscore: f64,
        isotope: Option<f64>,
        is_decoy: bool,
    ) -> CalibrationMatch {
        CalibrationMatch {
            entry_id,
            is_decoy,
            library_rt: 10.0,
            measured_rt: 10.0,
            score: correlation,
            ms1_error: None,
            library_precursor_mz: 400.0,
            observed_precursor_mz: Some(400.0),
            ms2_mass_errors: vec![],
            avg_ms2_error: None,
            n_matched_fragments: top6 as usize,
            n_library_fragments: 6,
            xcorr_score: 0.0,
            evalue: 1.0,
            isotope_cosine_score: isotope,
            sequence: "PEPTIDE".to_string(),
            charge: 2,
            scan_number: 100,
            hyperscore: 0.0,
            n_b_ions: 0,
            n_y_ions: 0,
            correlation_score: correlation,
            libcosine_apex: libcosine,
            top6_matched_apex: top6,
            hyperscore_apex: hyperscore,
            signal_to_noise: 10.0,          // Default value for tests
            peak_width_minutes: Some(0.25), // Default 15 sec for tests
            discriminant_score: 0.0,
            posterior_error: 0.0,
            q_value: 1.0,
        }
    }

    #[test]
    fn test_feature_matrix_extraction() {
        let matches = vec![
            make_test_match(1, 3.0, 0.8, 5, 10.0, Some(0.9), false),
            make_test_match(1 | 0x80000000, 1.5, 0.5, 3, 5.0, Some(0.7), true),
        ];

        // Test 4 features: correlation, libcosine, top6, snr
        let matrix = extract_feature_matrix(&matches, false);
        assert_eq!(matrix.rows, 2);
        assert_eq!(matrix.cols, 4);

        // Even with isotope flag, still 4 features (isotope removed from LDA)
        let matrix = extract_feature_matrix(&matches, true);
        assert_eq!(matrix.rows, 2);
        assert_eq!(matrix.cols, 4);

        // Check normalization (correlation 3.0 -> 0.5)
        let corr_norm = matrix[(0, 0)];
        assert!((corr_norm - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_lda_training() {
        // Create synthetic data: targets have higher scores
        // Use distinct sequences so fold assignment distributes across folds
        // Each target is paired with a decoy via entry_id (high bit = decoy)
        let mut matches = vec![];

        // 10 good targets with distinct sequences and unique entry_ids
        for i in 0..10 {
            let entry_id = (i + 1) as u32; // IDs 1-10
            let mut m = make_test_match(entry_id, 4.0, 0.9, 6, 15.0, Some(0.95), false);
            m.sequence = format!("TARGET{}", i);
            matches.push(m);
        }

        // 10 poor decoys with distinct sequences, paired with targets via entry_id
        for i in 0..10 {
            let entry_id = ((i + 1) as u32) | 0x80000000; // Paired with target i+1
            let mut m = make_test_match(entry_id, 1.0, 0.3, 2, 3.0, Some(0.5), true);
            m.sequence = format!("DECOY{}", i);
            matches.push(m);
        }

        let result = train_and_score_calibration(&mut matches, true);
        assert!(result.is_ok());

        // Targets should score higher than decoys after LDA
        let target_scores: Vec<f64> = matches
            .iter()
            .filter(|m| !m.is_decoy)
            .map(|m| m.discriminant_score)
            .collect();
        let decoy_scores: Vec<f64> = matches
            .iter()
            .filter(|m| m.is_decoy)
            .map(|m| m.discriminant_score)
            .collect();

        let avg_target = target_scores.iter().sum::<f64>() / target_scores.len() as f64;
        let avg_decoy = decoy_scores.iter().sum::<f64>() / decoy_scores.len() as f64;

        assert!(
            avg_target > avg_decoy,
            "Targets should score higher than decoys"
        );
    }

    #[test]
    fn test_compete_calibration_pairs_target_wins() {
        // Target has higher score → target wins
        let scores = vec![0.8, 0.3];
        let entry_ids = vec![1_u32, 1 | 0x80000000];
        let is_decoy = vec![false, true];
        let indices: Vec<usize> = vec![0, 1];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert_eq!(winners.len(), 1);
        assert_eq!(winners[0], 0); // target wins
    }

    #[test]
    fn test_compete_calibration_pairs_decoy_wins() {
        // Decoy has higher score → decoy wins
        let scores = vec![0.3, 0.8];
        let entry_ids = vec![1_u32, 1 | 0x80000000];
        let is_decoy = vec![false, true];
        let indices: Vec<usize> = vec![0, 1];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert_eq!(winners.len(), 1);
        assert_eq!(winners[0], 1); // decoy wins
    }

    #[test]
    fn test_compete_calibration_pairs_tie_goes_to_decoy() {
        // Equal scores → decoy wins (conservative)
        let scores = vec![0.5, 0.5];
        let entry_ids = vec![1_u32, 1 | 0x80000000];
        let is_decoy = vec![false, true];
        let indices: Vec<usize> = vec![0, 1];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert_eq!(winners.len(), 1);
        assert_eq!(winners[0], 1); // decoy wins on tie
    }

    #[test]
    fn test_compete_calibration_pairs_singleton_auto_wins() {
        // Unpaired entries auto-win
        let scores = vec![0.8, 0.6];
        let entry_ids = vec![1_u32, 2_u32]; // Different base IDs, no pairing
        let is_decoy = vec![false, false];
        let indices: Vec<usize> = vec![0, 1];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert_eq!(winners.len(), 2); // Both auto-win
    }

    #[test]
    fn test_compete_calibration_pairs_multiple_pairs() {
        // Two pairs: pair 1 target wins, pair 2 decoy wins
        let scores = vec![0.9, 0.2, 0.3, 0.7];
        let entry_ids = vec![1_u32, 1 | 0x80000000, 2_u32, 2 | 0x80000000];
        let is_decoy = vec![false, true, false, true];
        let indices: Vec<usize> = vec![0, 1, 2, 3];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert_eq!(winners.len(), 2);
        // Winners sorted by score descending: target(0.9), decoy(0.7)
        assert_eq!(winners[0], 0); // target from pair 1 (score 0.9)
        assert_eq!(winners[1], 3); // decoy from pair 2 (score 0.7)
    }

    #[test]
    fn test_compete_calibration_pairs_empty_input() {
        let scores: Vec<f64> = vec![];
        let entry_ids: Vec<u32> = vec![];
        let is_decoy: Vec<bool> = vec![];
        let indices: Vec<usize> = vec![];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert!(winners.is_empty());
    }

    #[test]
    fn test_compete_calibration_pairs_subset() {
        // Only a subset of indices are considered
        let scores = vec![0.9, 0.2, 0.3, 0.7];
        let entry_ids = vec![1_u32, 1 | 0x80000000, 2_u32, 2 | 0x80000000];
        let is_decoy = vec![false, true, false, true];
        // Only consider pair 2 (indices 2 and 3)
        let indices: Vec<usize> = vec![2, 3];

        let winners = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
        assert_eq!(winners.len(), 1);
        assert_eq!(winners[0], 3); // decoy wins (0.7 > 0.3)
    }

    /// Verify that compete_calibration_pairs produces deterministic results with tied scores.
    ///
    /// Winners are built from HashMap iteration (non-deterministic order), then sorted
    /// by score with index as tiebreaker. This test verifies tied-score entries always
    /// appear in the same order across multiple calls.
    #[test]
    fn test_compete_calibration_pairs_deterministic_with_ties() {
        // 4 pairs, targets all score 0.5 (tied), decoys all score 0.1
        let scores = vec![
            0.5, 0.5, 0.5, 0.5, // targets (indices 0-3)
            0.1, 0.1, 0.1, 0.1, // decoys (indices 4-7)
        ];
        let entry_ids: Vec<u32> = vec![
            10,
            20,
            30,
            40, // targets
            10 | 0x80000000,
            20 | 0x80000000,
            30 | 0x80000000,
            40 | 0x80000000, // decoys
        ];
        let is_decoy = vec![false, false, false, false, true, true, true, true];
        let indices: Vec<usize> = (0..8).collect();

        let first_result = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);

        // Run multiple times — HashMap iteration order varies, but output must be stable
        for _ in 0..20 {
            let result = compete_calibration_pairs(&scores, &entry_ids, &is_decoy, &indices);
            assert_eq!(
                result, first_result,
                "compete_calibration_pairs must be deterministic with tied scores"
            );
        }

        // All targets win (0.5 > 0.1), so winners are target indices
        assert_eq!(first_result.len(), 4);
        // With tied scores, winners must be sorted by base_id ascending
        assert_eq!(first_result, vec![0, 1, 2, 3]);
    }
}
