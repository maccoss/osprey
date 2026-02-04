//! Osprey I/O - File parsing and output for Osprey
//!
//! This crate provides parsers for input files (mzML, spectral libraries)
//! and writers for output files (blib, reports).
//!
//! ## MS1 Spectra for Calibration
//!
//! For MS1 mass calibration (pyXcorrDIA-compatible), use:
//! - `load_ms1_spectra()` - Load all MS1 spectra from an mzML file
//! - `MS1Index` - Index for efficient nearest-neighbor lookup by RT

pub mod library;
pub mod mzml;
pub mod output;

pub use library::{load_library, BlibLoader, DiannTsvLoader, ElibLoader};
pub use mzml::{load_ms1_spectra, MS1Index, MzmlReader};
pub use output::{BlibWriter, ReportWriter};
