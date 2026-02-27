//! Osprey I/O - File parsing and output for Osprey
//!
//! This crate provides parsers for input files (mzML, spectral libraries)
//! and writers for output files (blib, reports).
//!
//! ## Loading mzML Files
//!
//! For best performance, use `load_all_spectra()` to read both MS1 and MS2 in a single pass:
//! - `load_all_spectra()` - Load both MS1 and MS2 spectra efficiently (recommended)
//! - `load_ms1_spectra()` - Load only MS1 spectra
//! - `MS1Index` - Index for efficient nearest-neighbor lookup by RT

pub mod library;
pub mod mzml;
pub mod output;

pub use library::{load_library, BlibLoader, DiannTsvLoader, ElibLoader};
pub use mzml::{load_all_spectra, load_ms1_spectra, MS1Index, MzmlReader};
pub use output::BlibWriter;
