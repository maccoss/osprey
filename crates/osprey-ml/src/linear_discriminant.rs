//! Originally from Sage (https://github.com/lazear/sage)
//! Copyright (c) 2022 Michael Lazear
//! Licensed under the MIT License
//!
//! Linear Discriminant Analysis for binary classification
//!
//! Simplified version for Osprey - core LDA algorithm only

use super::gauss::Gauss;
use super::matrix::Matrix;

/// Linear Discriminant Analysis classifier
pub struct LinearDiscriminantAnalysis {
    eigenvector: Vec<f64>,
}

impl LinearDiscriminantAnalysis {
    /// Create LDA from pre-computed weights (eigenvector)
    ///
    /// # Arguments
    /// * `weights` - Feature weights (eigenvector)
    ///
    /// # Returns
    /// LDA model with given weights
    pub fn from_weights(weights: Vec<f64>) -> Result<LinearDiscriminantAnalysis, String> {
        if weights.is_empty() {
            return Err("Weights vector cannot be empty".to_string());
        }
        Ok(LinearDiscriminantAnalysis {
            eigenvector: weights,
        })
    }

    /// Fit LDA model to data
    ///
    /// # Arguments
    /// * `features` - Feature matrix (rows = samples, cols = features)
    /// * `decoy` - Boolean labels (true = decoy, false = target)
    ///
    /// # Returns
    /// Trained LDA model, or None if fitting failed
    pub fn fit(features: &Matrix, decoy: &[bool]) -> Option<LinearDiscriminantAnalysis> {
        assert_eq!(features.rows, decoy.len());

        // Calculate class means, and overall mean
        let x_bar = features.mean();
        let mut scatter_within = Matrix::zeros(features.cols, features.cols);
        let mut scatter_between = Matrix::zeros(features.cols, features.cols);

        let mut class_means = Vec::new();

        for class in [true, false] {
            let count = decoy.iter().filter(|&label| *label == class).count();

            let class_data = (0..features.rows)
                .zip(decoy)
                .filter(|&(_, label)| *label == class)
                .flat_map(|(row, _)| features.row(row))
                .collect::<Vec<_>>();

            let mut class_data = Matrix::new(class_data, count, features.cols);
            let class_mean = class_data.mean();

            for row in 0..class_data.rows {
                for col in 0..class_data.cols {
                    class_data[(row, col)] -= class_mean[col];
                }
            }

            let cov = class_data.transpose().dot(&class_data) / class_data.rows as f64;
            scatter_within += cov;

            let diff = Matrix::col_vector(
                class_mean
                    .iter()
                    .zip(x_bar.iter())
                    .map(|(x, y)| x - y)
                    .collect::<Vec<_>>(),
            );

            scatter_between += diff.dot(&diff.transpose());
            class_means.extend(class_mean);
        }

        // Use overall mean as the initial vector for power method
        let mut evec =
            Gauss::solve(scatter_within, scatter_between).map(|mat| mat.power_method(&x_bar))?;

        // Ensure target class scores higher than decoy class
        // (flip eigenvector sign if needed)
        let class_means = Matrix::new(class_means, 2, features.cols);
        let coef = class_means.dotv(&evec);
        if coef[1] < coef[0] {
            evec.iter_mut().for_each(|c| *c *= -1.0);
        }

        log::debug!("LDA eigenvector: {:?}", &evec[..evec.len().min(10)]);

        Some(LinearDiscriminantAnalysis { eigenvector: evec })
    }

    /// Score samples using the learned discriminant
    ///
    /// # Arguments
    /// * `features` - Feature matrix to score
    ///
    /// # Returns
    /// Discriminant scores (higher = more likely target)
    pub fn predict(&self, features: &Matrix) -> Vec<f64> {
        features.dotv(&self.eigenvector)
    }

    /// Get the learned eigenvector (feature weights)
    ///
    /// # Returns
    /// Reference to the eigenvector
    pub fn eigenvector(&self) -> &[f64] {
        &self.eigenvector
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::norm;

    fn all_close(a: &[f64], b: &[f64], eps: f64) -> bool {
        a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() < eps)
    }

    #[test]
    fn linear_discriminant() {
        let a = Matrix::new([1., 2., 3., 4.], 2, 2);
        let eigenvector = [0.4159736, 0.90937671];
        assert!(all_close(
            &a.power_method(&[0.54, 0.34]),
            &eigenvector,
            1E-5
        ));

        #[rustfmt::skip]
        let feats = Matrix::new(
            [
                5., 4., 3., 2.,
                4., 5., 4., 3.,
                6., 3., 4., 5.,
                1., 0., 2., 9.,
                5., 4., 4., 3.,
                2., 1., 1., 9.5,
                1., 0., 2., 8.,
                3., 2., -2., 10.,
            ],
            8,
            4,
        );

        let lda = LinearDiscriminantAnalysis::fit(
            &feats,
            &[false, false, false, true, false, true, true, true],
        )
        .expect("error training LDA");

        let mut scores = lda.predict(&feats);
        let norm_val = norm(&scores);
        scores = scores.into_iter().map(|s| s / norm_val).collect();

        let expected = [
            0.49706043,
            0.48920177,
            0.48920177,
            -0.07209359,
            0.51204672,
            -0.02849527,
            -0.04924864,
            -0.06055943,
        ];

        assert!(
            all_close(&scores, &expected, 1E-8),
            "{:?} {:?}",
            scores,
            expected
        );
    }
}
