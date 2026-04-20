//! Osprey - Peptide-centric DIA analysis tool
//!
//! Osprey is an open-source tool for peptide detection and quantification
//! in data-independent acquisition (DIA) mass spectrometry data.
//!
//! # Overview
//!
//! Osprey uses fragment XIC co-elution analysis to detect peptides in DIA data,
//! with machine learning scoring and rigorous FDR control.
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
//! run_analysis(config)?;
//! # Ok::<(), osprey::OspreyError>(())
//! ```

// Re-export core types
pub use osprey_core::*;

// Re-export I/O
pub use osprey_io::{load_library, BlibLoader, DiannTsvLoader, ElibLoader, MzmlReader};

// Re-export chromatography
pub use osprey_chromatography::{
    EmgFitter, EmgParameters, PeakDetector, RTCalibration, RTCalibrationStats, RTCalibrator,
    RTCalibratorConfig, RTStratifiedSampler,
};

// Re-export scoring
pub use osprey_scoring::{DecoyGenerator, DecoyMethod, Enzyme};

// Re-export FDR
pub use osprey_fdr::{FdrController, FdrCounts, MokapotResult, MokapotRunner};

// Two-tier logging (clean terminal + verbose log file)
pub mod logging;

// Cross-implementation bisection diagnostic dumps
pub mod diagnostics;

// Cross-run peak reconciliation
pub mod reconciliation;

// Runtime peptide trace (env-var gated)
pub mod trace;

// Pipeline
mod pipeline;
pub use pipeline::run_analysis;
