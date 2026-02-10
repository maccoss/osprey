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
    let n_passing_corr_only = compare_to_correlation_only(matches);
    log::info!(
        "  Comparison: Correlation-only would yield {} matches passing 1% FDR (LDA: {})",
        n_passing_corr_only, n_passing
    );

    // SAFETY CHECK: If LDA performs significantly worse than correlation-only, fall back
    if n_passing_corr_only > 100 && n_passing < n_passing_corr_only / 10 {
        log::warn!(
            "LDA yielded {}x fewer matches than correlation-only ({} vs {}) - falling back to correlation scoring",
            n_passing_corr_only / n_passing.max(1),
            n_passing,
            n_passing_corr_only
        );
        log::warn!("This suggests LDA learned incorrect feature weights for this dataset");

        // Re-score using correlation only
        return fallback_to_correlation_only(matches);
    }

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

/// Train LDA with 3-fold cross-validation and non-negative weight constraints
///
/// Follows Mokapot/Percolator methodology:
/// - 3-fold stratified CV for unbiased FDR estimates
/// - Non-negative weights applied to each fold (clip and renormalize)
/// - Peptides grouped by sequence (Percolator-RESET: keeps charge states together)
/// - Returns out-of-fold predictions for all samples
///
/// # Arguments
/// * `features` - Full feature matrix
/// * `decoy_labels` - Target/decoy labels (true = decoy)
/// * `sequences` - Peptide sequences (for grouping charge states together)
/// * `use_isotope_feature` - Whether isotope feature is included (for logging only, not used in features)
///
/// # Returns
/// Discriminant scores for all samples from cross-validation
fn train_lda_with_nonnegative_cv(
    features: &Matrix,
    decoy_labels: &[bool],
    sequences: &[String],
    _use_isotope_feature: bool,
) -> Result<Vec<f64>, OspreyError> {
    const N_FOLDS: usize = 3;
    const MAX_ITERATIONS: usize = 10;
    const CONVERGENCE_THRESHOLD: f64 = 0.2; // CV < 0.2 for all features = converged

    let n_samples = features.rows;
    let feature_names = vec!["correlation", "libcosine", "top6", "snr"];

    // Storage for weights across iterations to check convergence
    let mut iteration_weights: Vec<Vec<f64>> = Vec::new();
    let mut best_discriminants = vec![0.0; n_samples];

    // Iterate until convergence or max iterations
    for iteration in 0..MAX_ITERATIONS {
        log::debug!("  === LDA Iteration {}/{} ===", iteration + 1, MAX_ITERATIONS);

        // Create fold indices (stratified by target/decoy, grouped by peptide sequence)
        // Re-shuffle folds each iteration for robustness
        let fold_assignments = create_stratified_folds_by_peptide_shuffled(
            decoy_labels,
            sequences,
            N_FOLDS,
            iteration,
        );

        // Storage for cross-validated predictions
        let mut cv_discriminants = vec![0.0; n_samples];

        // Storage for fold weights within this iteration
        let mut fold_weights: Vec<Vec<f64>> = Vec::new();

        // Train and predict for each fold
        for fold_idx in 0..N_FOLDS {
            // Split into train and test sets
            let (train_features, test_features, train_labels, test_indices) =
                split_fold(features, decoy_labels, &fold_assignments, fold_idx);

            // Train LDA on training set
            let lda = LinearDiscriminantAnalysis::fit(&train_features, &train_labels)
                .ok_or_else(|| OspreyError::config(
                    format!("LDA training failed on fold {}/{}", fold_idx + 1, N_FOLDS)
                ))?;

            // Get weights and apply non-negative constraint
            let mut weights = lda.eigenvector().to_vec();

            // Count negative weights before clipping
            let n_negative = weights.iter().filter(|&w| *w < 0.0).count();
            if n_negative > 0 {
                log::debug!("    Fold {}/{}: Clipping {} negative weights", fold_idx + 1, N_FOLDS, n_negative);
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
                // If all weights were negative, use equal weights
                log::warn!("    Fold {}/{}: All weights negative - using equal weights", fold_idx + 1, N_FOLDS);
                let equal_weight = 1.0 / (weights.len() as f64).sqrt();
                for w in weights.iter_mut() {
                    *w = equal_weight;
                }
            }

            // Log fold weights
            log::debug!("    Fold {}/{} LDA weights: {}",
                fold_idx + 1, N_FOLDS,
                feature_names.iter().enumerate()
                    .map(|(i, name)| format!("{}={:.3}", name, weights[i]))
                    .collect::<Vec<_>>()
                    .join(", ")
            );

            // Store weights for convergence check within iteration
            fold_weights.push(weights.clone());

            // Create LDA with non-negative weights
            let lda_nonneg = LinearDiscriminantAnalysis::from_weights(weights)
                .map_err(|e| OspreyError::config(e))?;

            // Score test set
            let fold_predictions = lda_nonneg.predict(&test_features);

            // Store predictions in original order
            for (i, &original_idx) in test_indices.iter().enumerate() {
                cv_discriminants[original_idx] = fold_predictions[i];
            }
        }

        // Average weights across folds for this iteration
        let avg_weights = average_weights(&fold_weights);
        iteration_weights.push(avg_weights.clone());

        // Log iteration summary
        log::debug!("  Iteration {} avg weights: {}",
            iteration + 1,
            feature_names.iter().enumerate()
                .map(|(i, name)| format!("{}={:.3}", name, avg_weights[i]))
                .collect::<Vec<_>>()
                .join(", ")
        );

        // Store best discriminants from this iteration
        best_discriminants = cv_discriminants;

        // Check convergence across iterations (need at least 2 iterations)
        if iteration >= 1 {
            let converged = check_iteration_convergence(&iteration_weights, &feature_names, CONVERGENCE_THRESHOLD);
            if converged {
                log::info!("  LDA converged after {} iterations", iteration + 1);
                break;
            }
        }
    }

    // Final check: warn if didn't converge
    if iteration_weights.len() >= MAX_ITERATIONS {
        log::warn!("  LDA did not converge after {} iterations - using final iteration weights", MAX_ITERATIONS);
    }

    // Log final consensus weights (average across all iterations)
    let consensus_weights = average_weights(&iteration_weights);
    log::info!("  Cross-validated LDA complete ({}-fold, {} iterations, non-negative constraints)",
        N_FOLDS, iteration_weights.len()
    );
    log::info!("  Final model feature weights: {}",
        feature_names.iter().enumerate()
            .map(|(i, name)| format!("{}={:.3}", name, consensus_weights[i]))
            .collect::<Vec<_>>()
            .join(", ")
    );

    Ok(best_discriminants)
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

/// Check convergence of feature weights across iterations
///
/// Returns true if all features have CV < threshold across iterations
fn check_iteration_convergence(
    iteration_weights: &[Vec<f64>],
    feature_names: &[&str],
    threshold: f64,
) -> bool {
    if iteration_weights.len() < 2 {
        return false;
    }

    let n_features = iteration_weights[0].len();
    let n_iterations = iteration_weights.len();
    let mut all_converged = true;

    // Calculate coefficient of variation for each feature across iterations
    for feat_idx in 0..n_features {
        let weights: Vec<f64> = iteration_weights.iter().map(|w| w[feat_idx]).collect();
        let mean = weights.iter().sum::<f64>() / n_iterations as f64;
        let variance = weights.iter().map(|w| (w - mean).powi(2)).sum::<f64>() / n_iterations as f64;
        let std_dev = variance.sqrt();
        let cv = if mean > 1e-10 { std_dev / mean } else { 0.0 };

        if cv > threshold && mean > 0.05 {
            all_converged = false;
            log::debug!(
                "    Feature '{}' not yet converged: CV={:.3} (threshold={:.2})",
                feature_names[feat_idx], cv, threshold
            );
        }
    }

    all_converged
}

/// Create stratified fold assignments for cross-validation with shuffling
///
/// Implements Percolator-RESET approach:
/// - Groups PSMs by peptide sequence (keeps charge states together)
/// - Stratified by target/decoy to maintain class balance
/// - Shuffles fold assignments based on iteration for robustness
/// - Minimizes leakage between folds
fn create_stratified_folds_by_peptide_shuffled(
    labels: &[bool],
    sequences: &[String],
    n_folds: usize,
    iteration: usize,
) -> Vec<usize> {
    use std::collections::HashMap;

    let mut fold_assignments = vec![0; labels.len()];

    // Group indices by peptide sequence and target/decoy status
    let mut target_peptides: HashMap<String, Vec<usize>> = HashMap::new();
    let mut decoy_peptides: HashMap<String, Vec<usize>> = HashMap::new();

    for (i, (seq, &is_decoy)) in sequences.iter().zip(labels).enumerate() {
        if is_decoy {
            decoy_peptides.entry(seq.clone()).or_insert_with(Vec::new).push(i);
        } else {
            target_peptides.entry(seq.clone()).or_insert_with(Vec::new).push(i);
        }
    }

    // Convert to sorted vectors for consistent ordering
    let mut target_groups: Vec<(String, Vec<usize>)> = target_peptides.into_iter().collect();
    let mut decoy_groups: Vec<(String, Vec<usize>)> = decoy_peptides.into_iter().collect();

    // Sort by sequence for deterministic ordering
    target_groups.sort_by(|a, b| a.0.cmp(&b.0));
    decoy_groups.sort_by(|a, b| a.0.cmp(&b.0));

    // Assign target peptide groups to folds with rotation based on iteration
    // This changes fold assignments each iteration while keeping peptides together
    for (i, (_seq, indices)) in target_groups.iter().enumerate() {
        let fold = (i + iteration) % n_folds;
        for &idx in indices {
            fold_assignments[idx] = fold;
        }
    }

    // Assign decoy peptide groups to folds with rotation
    for (i, (_seq, indices)) in decoy_groups.iter().enumerate() {
        let fold = (i + iteration) % n_folds;
        for &idx in indices {
            fold_assignments[idx] = fold;
        }
    }

    log::debug!(
        "    Created {} folds (iteration {}): {} target peptides ({} PSMs), {} decoy peptides ({} PSMs)",
        n_folds,
        iteration + 1,
        target_groups.len(),
        target_groups.iter().map(|(_, v)| v.len()).sum::<usize>(),
        decoy_groups.len(),
        decoy_groups.iter().map(|(_, v)| v.len()).sum::<usize>()
    );

    fold_assignments
}

/// Split data into train and test sets for a specific fold
///
/// # Returns
/// (train_features, test_features, train_labels, test_indices)
fn split_fold(
    features: &Matrix,
    labels: &[bool],
    fold_assignments: &[usize],
    test_fold: usize,
) -> (Matrix, Matrix, Vec<bool>, Vec<usize>) {
    let mut train_rows = Vec::new();
    let mut test_rows = Vec::new();
    let mut train_labels = Vec::new();
    let mut test_indices = Vec::new();

    for (i, &fold) in fold_assignments.iter().enumerate() {
        if fold == test_fold {
            test_rows.push(i);
            test_indices.push(i);
        } else {
            train_rows.push(i);
            train_labels.push(labels[i]);
        }
    }

    // Extract rows for train and test sets
    let train_features = extract_rows(features, &train_rows);
    let test_features = extract_rows(features, &test_rows);

    log::debug!(
        "  Fold {}: {} train, {} test",
        test_fold + 1, train_rows.len(), test_rows.len()
    );

    (train_features, test_features, train_labels, test_indices)
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

/// Fall back to correlation-only scoring when LDA fails
///
/// Sorts by correlation, calculates q-values, and updates match scores
fn fallback_to_correlation_only(matches: &mut [CalibrationMatch]) -> Result<usize, OspreyError> {
    // Sort by correlation_score descending
    matches.sort_by(|a, b| b.correlation_score.total_cmp(&a.correlation_score));

    // Build is_decoy array
    let is_decoy: Vec<bool> = matches.iter().map(|m| m.is_decoy).collect();

    // Calculate q-values
    let mut q_values = vec![0.0; matches.len()];
    let n_passing = qvalue::calculate_q_values(&is_decoy, &mut q_values);

    // Update matches with correlation-based scores and q-values
    for (i, m) in matches.iter_mut().enumerate() {
        m.discriminant_score = m.correlation_score; // Use correlation as discriminant
        m.q_value = q_values[i];
    }

    log::info!(
        "Correlation-only scoring: {} matches passing 1% FDR",
        n_passing
    );

    Ok(n_passing)
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
        let mut matches = vec![];

        // 10 good targets
        for _ in 0..10 {
            matches.push(make_test_match(4.0, 0.9, 6, 15.0, Some(0.95), false));
        }

        // 10 poor decoys
        for _ in 0..10 {
            matches.push(make_test_match(1.0, 0.3, 2, 3.0, Some(0.5), true));
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
