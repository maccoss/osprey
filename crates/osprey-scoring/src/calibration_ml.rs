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

    log::info!(
        "Training LDA with 3-fold cross-validation on {} calibration matches (4 features)",
        matches.len()
    );

    // 1. Extract feature matrix
    let features = extract_feature_matrix(matches, use_isotope_feature);

    // 2. Build target/decoy labels (true = decoy)
    let decoy_labels: Vec<bool> = matches.iter().map(|m| m.is_decoy).collect();

    // 3. Perform 3-fold cross-validation
    let discriminants = train_lda_with_cross_validation(&features, &decoy_labels, use_isotope_feature)?;

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

    // Log median feature values for diagnostics
    log_median_features(matches, use_isotope_feature);

    Ok(n_passing)
}

/// Perform 3-fold cross-validation LDA training
///
/// # Arguments
/// * `features` - Full feature matrix
/// * `decoy_labels` - Target/decoy labels (true = decoy)
/// * `use_isotope_feature` - Whether isotope feature is included
///
/// # Returns
/// Discriminant scores for all samples from cross-validation
fn train_lda_with_cross_validation(
    features: &Matrix,
    decoy_labels: &[bool],
    _use_isotope_feature: bool,
) -> Result<Vec<f64>, OspreyError> {
    const N_FOLDS: usize = 3;
    let n_samples = features.rows;

    // Create fold indices (stratified by target/decoy to maintain class balance)
    let fold_assignments = create_stratified_folds(decoy_labels, N_FOLDS);

    // Storage for cross-validated predictions
    let mut cv_discriminants = vec![0.0; n_samples];

    let feature_names = vec!["correlation", "libcosine", "top6", "snr"];

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

        // Log fold weights for diagnostics
        log::debug!("  Fold {}/{} LDA weights: {}",
            fold_idx + 1, N_FOLDS,
            feature_names.iter().enumerate()
                .map(|(i, name)| format!("{}={:.3}", name, lda.eigenvector()[i]))
                .collect::<Vec<_>>()
                .join(", ")
        );

        // Score test set
        let fold_predictions = lda.predict(&test_features);

        // Store predictions in original order
        for (i, &original_idx) in test_indices.iter().enumerate() {
            cv_discriminants[original_idx] = fold_predictions[i];
        }
    }

    // Train final model on all data for logging feature weights
    let final_lda = LinearDiscriminantAnalysis::fit(features, decoy_labels)
        .ok_or_else(|| OspreyError::config("Final LDA training failed"))?;

    log::info!("  Cross-validated LDA complete");
    log::info!("  Final model feature weights: {}",
        feature_names.iter().enumerate()
            .map(|(i, name)| format!("{}={:.3}", name, final_lda.eigenvector()[i]))
            .collect::<Vec<_>>()
            .join(", ")
    );

    Ok(cv_discriminants)
}

/// Create stratified fold assignments for cross-validation
///
/// Ensures each fold has balanced target/decoy ratios
fn create_stratified_folds(labels: &[bool], n_folds: usize) -> Vec<usize> {
    let n_samples = labels.len();
    let mut fold_assignments = vec![0; n_samples];

    // Separate target and decoy indices
    let target_indices: Vec<usize> = labels.iter()
        .enumerate()
        .filter(|(_, &is_decoy)| !is_decoy)
        .map(|(i, _)| i)
        .collect();

    let decoy_indices: Vec<usize> = labels.iter()
        .enumerate()
        .filter(|(_, &is_decoy)| is_decoy)
        .map(|(i, _)| i)
        .collect();

    // Assign targets to folds round-robin
    for (i, &idx) in target_indices.iter().enumerate() {
        fold_assignments[idx] = i % n_folds;
    }

    // Assign decoys to folds round-robin
    for (i, &idx) in decoy_indices.iter().enumerate() {
        fold_assignments[idx] = i % n_folds;
    }

    log::debug!(
        "  Created {} folds: {} targets, {} decoys",
        n_folds, target_indices.len(), decoy_indices.len()
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
