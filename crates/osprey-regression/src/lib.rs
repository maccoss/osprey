//! Osprey Regression - Spectrum deconvolution via ridge regression
//!
//! This crate implements the core regression algorithm for deconvoluting
//! mixed DIA spectra into individual peptide contributions.
//!
//! The main components are:
//! - **Binning**: Convert continuous spectra to discrete bins (unit resolution)
//! - **Matrix**: Build design matrices from library spectra
//! - **Ridge**: Solve regularized least squares problems
//! - **Sparse**: HRAM support with sparse matrices and ppm-based matching

pub mod binning;
pub mod matrix;
pub mod ridge;
pub mod sparse;

pub use binning::Binner;
pub use matrix::DesignMatrixBuilder;
pub use ridge::RidgeSolver;
pub use sparse::{HramConfig, SparseMatrixBuilder, SparseRidgeSolver};
