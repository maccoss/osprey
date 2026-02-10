//! Osprey - Peptide-centric DIA analysis tool
//!
//! Osprey is an open-source tool for peptide detection and quantification
//! in data-independent acquisition (DIA) mass spectrometry data.
//!
//! # Overview
//!
//! Osprey uses ridge regression to deconvolute mixed MS/MS spectra,
//! aggregates evidence across the chromatographic dimension, and uses
//! machine learning to score peptide detections with rigorous FDR control.
//!
//! # Quick Start
//!
//! ```no_run
//! use osprey::{OspreyConfig, LibrarySource, run_analysis};
//! use std::path::PathBuf;
//!
//! let config = OspreyConfig {
//!     input_files: vec![PathBuf::from("sample.mzML")],
//!     library_source: LibrarySource::DiannTsv(PathBuf::from("library.tsv")),
//!     output_blib: PathBuf::from("results.blib"),
//!     ..Default::default()
//! };
//!
//! let results = run_analysis(config)?;
//! println!("Detected {} peptides", results.len());
//! # Ok::<(), osprey::OspreyError>(())
//! ```

// Re-export core types
pub use osprey_core::*;

// Re-export I/O
pub use osprey_io::{load_library, BlibLoader, DiannTsvLoader, ElibLoader, MzmlReader, ReportWriter};

// Re-export regression
pub use osprey_regression::{
    Binner, DesignMatrixBuilder, RidgeSolver,
};

// Re-export chromatography
pub use osprey_chromatography::{
    EmgFitter, EmgParameters, PeakDetector, RTCalibration, RTCalibrationStats, RTCalibrator,
    RTCalibratorConfig, RTStratifiedSampler,
};

// Re-export scoring
pub use osprey_scoring::{DecoyGenerator, DecoyMethod, Enzyme, FeatureExtractor};

// Re-export FDR
pub use osprey_fdr::{FdrController, FdrCounts, MokapotResult, MokapotRunner, PsmFeatures};

// Pipeline
mod pipeline;
pub use pipeline::run_analysis;
