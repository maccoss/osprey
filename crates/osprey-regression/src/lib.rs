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
//! - **Optimized**: High-performance regression with pre-binned library and f32 ops

// Ensure BLAS is linked for ndarray operations
extern crate blas_src;
extern crate openblas_src;

pub mod binning;
pub mod matrix;
pub mod optimized;
pub mod ridge;
pub mod sparse;

pub use binning::Binner;
pub use matrix::DesignMatrixBuilder;
pub use optimized::{BinnedLibrary, BinnedSpectraCache, OptimizedProcessor, OptimizedRegressionResult, OptimizedSolver};
pub use ridge::RidgeSolver;
pub use sparse::{HramConfig, SparseMatrixBuilder, SparseRidgeSolver};
