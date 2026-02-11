//! Machine learning-based scoring for calibration matches
//!
//! Uses Linear Discriminant Analysis (LDA) to learn optimal feature weights
//! from target-decoy competition, replacing simple correlation-based scoring.

use crate::batch::CalibrationMatch;
use osprey_core::OspreyError;
use osprey_ml::{kde, linear_discriminant::LinearDiscriminantAnalysis, matrix::Matrix, qvalue};

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
    log::info!(
        "Training LDA with 3-fold cross-validation on {} calibration matches ({} features)",
        matches.len(),
        n_features
    );

    // 1. Extract feature matrix
    let features = extract_feature_matrix(matches, use_isotope_feature);

    // 2. Build target/decoy labels (true = decoy)
    let decoy_labels: Vec<bool> = matches.iter().map(|m| m.is_decoy).collect();

    // 3. Extract sequences for fold grouping (keep charge states together)
    let sequences: Vec<String> = matches.iter().map(|m| m.sequence.clone()).collect();

    // 4. Train LDA with cross-validation and non-negative constraints
    let discriminants = train_lda_with_nonnegative_cv(
        &features,
        &decoy_labels,
        &sequences,
        use_isotope_feature
    )?;

    // 5. KDE for posterior error probabilities
    let kde_model = kde::Builder::default().build(&discriminants, &decoy_labels);

    // 6. Assign discriminant scores and posterior errors
    for (i, m) in matches.iter_mut().enumerate() {
        m.discriminant_score = discriminants[i];
        m.posterior_error = kde_model.posterior_error(discriminants[i]);
    }

    // 7. Sort by discriminant score descending (best matches first)
    matches.sort_by(|a, b| b.discriminant_score.total_cmp(&a.discriminant_score));

    // DEBUG: Log discriminant score distributions
    let target_scores: Vec<f64> = matches.iter()
        .filter(|m| !m.is_decoy)
        .map(|m| m.discriminant_score)
        .collect();
    let decoy_scores: Vec<f64> = matches.iter()
        .filter(|m| m.is_decoy)
        .map(|m| m.discriminant_score)
        .collect();

    if !target_scores.is_empty() && !decoy_scores.is_empty() {
        let target_mean = target_scores.iter().sum::<f64>() / target_scores.len() as f64;
        let decoy_mean = decoy_scores.iter().sum::<f64>() / decoy_scores.len() as f64;
        let target_max = target_scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let target_min = target_scores.iter().cloned().fold(f64::INFINITY, f64::min);
        let decoy_max = decoy_scores.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let decoy_min = decoy_scores.iter().cloned().fold(f64::INFINITY, f64::min);

        log::info!("  Discriminant scores: target_mean={:.3}, decoy_mean={:.3}", target_mean, decoy_mean);
        log::info!("  Discriminant ranges: targets=[{:.3}, {:.3}], decoys=[{:.3}, {:.3}]",
            target_min, target_max, decoy_min, decoy_max);
    }

    // 8. Calculate q-values using target-decoy competition
    let is_decoy: Vec<bool> = matches.iter().map(|m| m.is_decoy).collect();
    let mut q_values = vec![0.0; matches.len()];
    let n_passing = qvalue::calculate_q_values(&is_decoy, &mut q_values);

    // 9. Assign q-values back to matches
    for (i, m) in matches.iter_mut().enumerate() {
        m.q_value = q_values[i];
    }

    // DEBUG: Log q-value statistics
    let min_q = q_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let first_10_q: Vec<f64> = q_values.iter().take(10).cloned().collect();
    log::info!("  Q-value stats: min_q={:.4}, first_10={:?}", min_q, first_10_q);

    log::info!(
        "LDA scoring complete: {} matches passing 1% FDR",
        n_passing
    );

    // DIAGNOSTIC: Compare LDA performance to correlation-only
    // Note: The LDA training already tracks the best single feature internally
    // and guarantees the result is at least as good as the best single feature.
    let n_passing_corr_only = compare_to_correlation_only(matches);
    log::info!(
        "  Comparison: Correlation-only would yield {} matches passing 1% FDR (LDA: {})",
        n_passing_corr_only, n_passing
    );

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
        vec!["correlation", "libcosine", "top6", "snr", "isotope"]
    } else {
        vec!["correlation", "libcosine", "top6", "snr"]
    };

    log::info!("  LDA feature weights (non-negative): {}",
        feature_names.iter().enumerate()
            .map(|(i, name)| format!("{}={:.3}", name, weights[i]))
            .collect::<Vec<_>>()
            .join(", ")
    );

    // Create new LDA with non-negative weights
    let lda_nonneg = LinearDiscriminantAnalysis::from_weights(weights)
        .map_err(|e| OspreyError::config(e))?;

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
    sequences: &[String],
    _use_isotope_feature: bool,
) -> Result<Vec<f64>, OspreyError> {
    const N_FOLDS: usize = 3;
    const MAX_ITERATIONS: usize = 3;
    const TRAIN_FDR: f64 = 0.01;
    const MIN_POSITIVE_EXAMPLES: usize = 50;

    let n_samples = features.rows;
    let n_features = features.cols;
    let feature_names = ["correlation", "libcosine", "top6", "snr"];

    // Create fold assignments (stable across iterations - Percolator-RESET grouping)
    let fold_assignments = create_stratified_folds_by_peptide(
        decoy_labels,
        sequences,
        N_FOLDS,
    );

    // Find the best single feature to use as baseline (like Percolator)
    let mut best_feat_idx = 0;
    let mut best_feat_passing = 0;
    for feat_idx in 0..n_features {
        let feat_scores: Vec<f64> = (0..n_samples)
            .map(|i| features[(i, feat_idx)])
            .collect();
        let n_pass = count_passing_targets(&feat_scores, decoy_labels, 0.01);
        log::info!("  Initial feature '{}': {} pass 1% FDR", feature_names[feat_idx], n_pass);
        if n_pass > best_feat_passing {
            best_feat_passing = n_pass;
            best_feat_idx = feat_idx;
        }
    }

    // Initialize with best single feature as both current and best-so-far
    let baseline_scores: Vec<f64> = (0..n_samples)
        .map(|i| features[(i, best_feat_idx)])
        .collect();

    log::info!(
        "  Baseline: '{}' = {} pass 1% FDR",
        feature_names[best_feat_idx], best_feat_passing
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
                &train_indices,
                TRAIN_FDR,
                MIN_POSITIVE_EXAMPLES,
            );

            total_selected_targets += selected_target_indices.len();

            // Collect ALL decoy indices from train set
            let decoy_indices: Vec<usize> = train_indices.iter()
                .filter(|&&i| decoy_labels[i])
                .copied()
                .collect();

            // Build training set: selected targets + all decoys
            let mut training_indices = selected_target_indices.clone();
            training_indices.extend_from_slice(&decoy_indices);

            let train_features = extract_rows(features, &training_indices);
            let train_labels: Vec<bool> = training_indices.iter()
                .map(|&i| decoy_labels[i])
                .collect();

            log::debug!(
                "    Fold {}/{}: {} selected targets + {} decoys = {} training samples",
                fold_idx + 1, N_FOLDS,
                selected_target_indices.len(),
                decoy_indices.len(),
                training_indices.len()
            );

            // Train LDA on clean training set
            let lda = LinearDiscriminantAnalysis::fit(&train_features, &train_labels)
                .ok_or_else(|| OspreyError::config(
                    format!("LDA training failed on fold {}/{}", fold_idx + 1, N_FOLDS)
                ))?;

            let weights = lda.eigenvector().to_vec();

            log::debug!("    Fold {}/{} weights: {}",
                fold_idx + 1, N_FOLDS,
                feature_names.iter().enumerate()
                    .map(|(i, name)| format!("{}={:.3}", name, weights[i]))
                    .collect::<Vec<_>>()
                    .join(", ")
            );

            fold_weights.push(weights);
        }

        // Phase 2: Average weights across folds → consensus weights
        let mut consensus_weights = average_weights(&fold_weights);

        // Clip negative weights to zero (non-negative constraint)
        let n_negative = consensus_weights.iter().filter(|&&w| w < 0.0).count();
        if n_negative > 0 {
            log::debug!("    Clipping {} negative consensus weights to zero", n_negative);
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
            .map_err(|e| OspreyError::config(e))?;
        let new_scores = lda_consensus.predict(features);

        let n_passing = count_passing_targets(&new_scores, decoy_labels, 0.01);

        log::info!(
            "  Iteration {}: weights=[{}], selected {} train targets, {} pass 1% FDR",
            iteration + 1,
            feature_names.iter().enumerate()
                .map(|(i, name)| format!("{}={:.3}", name, consensus_weights[i]))
                .collect::<Vec<_>>()
                .join(", "),
            total_selected_targets / N_FOLDS,
            n_passing,
        );

        iteration_passing.push(n_passing);

        // Track best iteration: keep scores from whichever iteration gave the most passing
        if n_passing > best_passing {
            best_scores = new_scores.clone();
            best_passing = n_passing;
            best_iteration = iteration + 1;
            consecutive_no_improve = 0;
            log::info!("    -> New best: {} pass 1% FDR (iteration {})", best_passing, best_iteration);
        } else {
            consecutive_no_improve += 1;
            log::info!(
                "    -> No improvement (best remains {} from iteration {}, {} consecutive non-improvements)",
                best_passing, best_iteration, consecutive_no_improve
            );
        }

        // Update current scores for next iteration's target selection
        current_scores = new_scores;

        // Stop early if 2 consecutive iterations didn't improve
        if consecutive_no_improve >= 2 {
            log::info!(
                "  Stopping early: {} consecutive iterations without improvement",
                consecutive_no_improve
            );
            break;
        }
    }

    // Log final result
    log::info!("  Iteration history (pass 1% FDR): {:?}", iteration_passing);
    log::info!(
        "  Using iteration {} ({} pass 1% FDR)",
        best_iteration, best_passing
    );

    Ok(best_scores)
}

/// Select high-confidence targets from a subset of data for positive training set
///
/// This is the key Percolator/Mokapot innovation: don't train on all targets (noisy),
/// train only on targets that pass an FDR threshold (clean positive examples).
///
/// # Arguments
/// * `scores` - Current discriminant scores for ALL samples
/// * `decoy_labels` - Target/decoy labels for ALL samples
/// * `subset_indices` - Indices of samples to consider (train set)
/// * `fdr_threshold` - Q-value threshold for selecting targets (e.g., 0.01)
/// * `min_targets` - Minimum number of targets to select
///
/// # Returns
/// Indices (into the original arrays) of selected high-confidence targets
fn select_positive_training_set(
    scores: &[f64],
    decoy_labels: &[bool],
    subset_indices: &[usize],
    fdr_threshold: f64,
    min_targets: usize,
) -> Vec<usize> {
    // Sort subset by score descending
    let mut sorted_indices: Vec<usize> = subset_indices.to_vec();
    sorted_indices.sort_by(|&a, &b| scores[b].partial_cmp(&scores[a]).unwrap());

    // Calculate q-values using target-decoy competition within subset
    let mut targets_so_far = 0usize;
    let mut decoys_so_far = 0usize;
    let mut q_values = vec![1.0; sorted_indices.len()];

    for (rank, &idx) in sorted_indices.iter().enumerate() {
        if decoy_labels[idx] {
            decoys_so_far += 1;
        } else {
            targets_so_far += 1;
        }

        if targets_so_far > 0 {
            q_values[rank] = decoys_so_far as f64 / targets_so_far as f64;
        }
    }

    // Make q-values monotonically non-decreasing from the bottom
    let mut min_q = f64::INFINITY;
    for q in q_values.iter_mut().rev() {
        if *q < min_q {
            min_q = *q;
        } else {
            *q = min_q;
        }
    }

    // Select targets with q < threshold
    let mut selected: Vec<usize> = Vec::new();
    for (rank, &idx) in sorted_indices.iter().enumerate() {
        if !decoy_labels[idx] && q_values[rank] <= fdr_threshold {
            selected.push(idx);
        }
    }

    // If too few, relax threshold progressively
    if selected.len() < min_targets {
        let relaxed_thresholds = [0.05, 0.10, 0.25, 0.50];
        for &threshold in &relaxed_thresholds {
            selected.clear();
            for (rank, &idx) in sorted_indices.iter().enumerate() {
                if !decoy_labels[idx] && q_values[rank] <= threshold {
                    selected.push(idx);
                }
            }
            if selected.len() >= min_targets {
                log::debug!("      Relaxed training FDR to {:.0}% to get {} targets", threshold * 100.0, selected.len());
                break;
            }
        }
    }

    selected
}

/// Count targets passing a given FDR threshold
fn count_passing_targets(scores: &[f64], decoy_labels: &[bool], fdr_threshold: f64) -> usize {
    // Sort by score descending
    let mut indices: Vec<usize> = (0..scores.len()).collect();
    indices.sort_by(|&a, &b| scores[b].partial_cmp(&scores[a]).unwrap());

    let mut targets = 0usize;
    let mut decoys = 0usize;
    let mut n_passing = 0usize;
    let mut min_q = f64::INFINITY;

    // Forward pass to compute q-values and count passing
    let mut q_values = vec![1.0; indices.len()];
    for (rank, &idx) in indices.iter().enumerate() {
        if decoy_labels[idx] {
            decoys += 1;
        } else {
            targets += 1;
        }
        if targets > 0 {
            q_values[rank] = decoys as f64 / targets as f64;
        }
    }

    // Make monotonic from bottom
    for q in q_values.iter_mut().rev() {
        if *q < min_q {
            min_q = *q;
        } else {
            *q = min_q;
        }
    }

    // Count targets passing
    for (rank, &idx) in indices.iter().enumerate() {
        if !decoy_labels[idx] && q_values[rank] <= fdr_threshold {
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
        .flat_map(|&row| {
            (0..n_cols).map(move |col| matrix[(row, col)])
        })
        .collect();

    Matrix::new(data, row_indices.len(), n_cols)
}

/// Extract feature matrix from calibration matches
///
/// Features are normalized to similar ranges for fair weighting:
/// - correlation: 0-1 (typical range 0-6, divided by 6)
/// - libcosine_apex: 0-1 (already normalized)
/// - top6_matched_apex: 0-1 (count 0-6, divided by 6)
/// - signal_to_noise: ~0-1 (ln(x+1) transform, typical max ~100 → ln~5)
///
/// Note: hyperscore and isotope are NOT used in calibration LDA:
/// - hyperscore dominates target-decoy discrimination but doesn't correlate with RT quality
/// - isotope shows negative correlation and doesn't help RT calibration
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
                // Signal-to-noise with log transform (typical range 1-100+)
                (m.signal_to_noise.ln_1p() / 5.0).clamp(0.0, 1.0),
            ]
        })
        .collect();

    Matrix::new(features, matches.len(), n_features)
}

/// Compare LDA performance to correlation-only scoring
///
/// Returns the number of matches that would pass 1% FDR using correlation score alone
fn compare_to_correlation_only(matches: &[CalibrationMatch]) -> usize {
    // Create a sorted copy by correlation_score
    let mut indices: Vec<usize> = (0..matches.len()).collect();
    indices.sort_by(|&a, &b| {
        matches[b].correlation_score.total_cmp(&matches[a].correlation_score)
    });

    // Build is_decoy array in correlation-sorted order
    let is_decoy_corr: Vec<bool> = indices.iter()
        .map(|&i| matches[i].is_decoy)
        .collect();

    // Calculate q-values for correlation-only ranking
    let mut q_values_corr = vec![0.0; matches.len()];
    let n_passing_corr = qvalue::calculate_q_values(&is_decoy_corr, &mut q_values_corr);

    // DEBUG: Log correlation-only statistics
    let min_q_corr = q_values_corr.iter().cloned().fold(f64::INFINITY, f64::min);
    log::debug!(
        "  Correlation-only: min_q={:.4}, n_passing={}",
        min_q_corr, n_passing_corr
    );

    n_passing_corr
}

/// Log median feature values for diagnostics
fn log_median_features(matches: &[CalibrationMatch], use_isotope: bool) {
    if matches.is_empty() {
        return;
    }

    let mut corrs: Vec<f64> = matches.iter().map(|m| m.correlation_score).collect();
    let mut libcos: Vec<f64> = matches.iter().map(|m| m.libcosine_apex).collect();
    let mut top6s: Vec<f64> = matches.iter().map(|m| m.top6_matched_apex as f64).collect();
    let mut snrs: Vec<f64> = matches.iter().map(|m| m.signal_to_noise).collect();
    let mut discrims: Vec<f64> = matches.iter().map(|m| m.discriminant_score).collect();

    corrs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    libcos.sort_by(|a, b| a.partial_cmp(b).unwrap());
    top6s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    snrs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    discrims.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mid = matches.len() / 2;

    log::info!(
        "  Median features: corr={:.3}, libcos={:.3}, top6={:.1}, snr={:.1}, discrim={:.3}",
        corrs[mid],
        libcos[mid],
        top6s[mid],
        snrs[mid],
        discrims[mid]
    );

    if use_isotope {
        let mut isos: Vec<f64> = matches
            .iter()
            .filter_map(|m| m.isotope_cosine_score)
            .collect();
        if !isos.is_empty() {
            isos.sort_by(|a, b| a.partial_cmp(b).unwrap());
            log::info!("  Median isotope_cosine: {:.3}", isos[isos.len() / 2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_match(
        correlation: f64,
        libcosine: f64,
        top6: u8,
        hyperscore: f64,
        isotope: Option<f64>,
        is_decoy: bool,
    ) -> CalibrationMatch {
        CalibrationMatch {
            entry_id: 1,
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
            signal_to_noise: 10.0,  // Default value for tests
            discriminant_score: 0.0,
            posterior_error: 0.0,
            q_value: 1.0,
        }
    }

    #[test]
    fn test_feature_matrix_extraction() {
        let matches = vec![
            make_test_match(3.0, 0.8, 5, 10.0, Some(0.9), false),
            make_test_match(1.5, 0.5, 3, 5.0, Some(0.7), true),
        ];

        // Test 4 features: correlation, libcosine, top6, snr (isotope not used in LDA)
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
        let mut matches = vec![];

        // 10 good targets with distinct sequences
        for i in 0..10 {
            let mut m = make_test_match(4.0, 0.9, 6, 15.0, Some(0.95), false);
            m.sequence = format!("TARGET{}", i);
            matches.push(m);
        }

        // 10 poor decoys with distinct sequences
        for i in 0..10 {
            let mut m = make_test_match(1.0, 0.3, 2, 3.0, Some(0.5), true);
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
}
