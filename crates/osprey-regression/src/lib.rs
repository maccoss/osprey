//! Osprey Regression - Spectrum deconvolution via ridge regression
//!
//! This crate implements the core regression algorithm for deconvoluting
//! mixed DIA spectra into individual peptide contributions.
//!
//! The main components are:
//! - **Binning**: Convert continuous spectra to discrete bins (unit resolution)
//! - **Matrix**: Build design matrices from library spectra
//! - **Ridge**: Solve regularized least squares problems
//! - **Solver**: High-performance regression with pre-binned library and f32 ops
//! - **CD-NNLS**: Coordinate descent solver for non-negative least squares

// Ensure BLAS is linked for ndarray operations
extern crate blas_src;
extern crate openblas_src;

pub mod binning;
pub mod cd_nnls;
pub mod matrix;
pub mod ridge;
pub mod solver;
pub use binning::Binner;
pub use cd_nnls::CdNnlsParams;
pub use matrix::DesignMatrixBuilder;
pub use solver::OptimizedSolver;
pub use ridge::RidgeSolver;
