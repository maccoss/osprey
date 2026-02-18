//! Machine Learning for Osprey
//!
//! This crate provides machine learning algorithms for FDR control:
//! - Linear Support Vector Machine (SVM) for Percolator-style scoring
//! - Posterior Error Probability (PEP) estimation
//! - Linear Discriminant Analysis (LDA)
//! - Kernel Density Estimation (KDE)
//! - Q-value calculation
//!
//! ## Attribution
//!
//! This code is adapted from Sage (https://github.com/lazear/sage)
//! Copyright (c) 2022 Michael Lazear
//! Licensed under the MIT License

pub mod gauss;
pub mod kde;
pub mod linear_discriminant;
pub mod matrix;
pub mod pep;
pub mod qvalue;
pub mod svm;

#[allow(dead_code)]
fn all_close(lhs: &[f64], rhs: &[f64], eps: f64) -> bool {
    lhs.iter()
        .zip(rhs.iter())
        .all(|(l, r)| (l - r).abs() <= eps)
}

pub fn norm(slice: &[f64]) -> f64 {
    slice.iter().fold(0.0, |acc, x| acc + x.powi(2)).sqrt()
}

pub fn mean(slice: &[f64]) -> f64 {
    slice.iter().sum::<f64>() / slice.len() as f64
}

pub fn std(slice: &[f64]) -> f64 {
    let mean = mean(slice);
    let x = slice.iter().fold(0.0, |acc, x| acc + (x - mean).powi(2));
    (x / slice.len() as f64).sqrt()
}
