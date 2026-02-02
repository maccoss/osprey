//! Osprey I/O - File parsing and output for Osprey
//!
//! This crate provides parsers for input files (mzML, spectral libraries)
//! and writers for output files (blib, reports).

pub mod library;
pub mod mzml;
pub mod output;

pub use library::{load_library, BlibLoader, DiannTsvLoader, ElibLoader};
pub use mzml::MzmlReader;
pub use output::{BlibWriter, ReportWriter};
