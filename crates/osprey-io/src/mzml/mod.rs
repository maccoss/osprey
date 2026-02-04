//! mzML file parsing
//!
//! This module provides parsing of mzML files using the mzdata crate,
//! converting to Osprey's internal spectrum representation.
//!
//! ## MS1 Spectra for Calibration
//!
//! For MS1 mass calibration (pyXcorrDIA-compatible), use:
//! - `load_ms1_spectra()` - Load all MS1 spectra from an mzML file
//! - `MS1Index` - Index for efficient nearest-neighbor lookup by RT
//!
//! ## Streaming Parser (optional, requires `streaming` feature)
//!
//! For async streaming parsing with pipelined preprocessing:
//! - `StreamingMzmlReader` - Async event-based parser
//! - Sends spectra through channels as they are parsed
//! - Memory-efficient for large files

mod parser;

#[cfg(feature = "streaming")]
pub mod streaming;

pub use parser::{load_ms1_spectra, MS1Index, MzmlReader};

#[cfg(feature = "streaming")]
pub use streaming::{ParseStats, StreamingMzmlError, StreamingMzmlReader};
