//! mzML file parsing
//!
//! This module provides parsing of mzML files using the mzdata crate,
//! converting to Osprey's internal spectrum representation.
//!
//! ## Single-Pass Loading (Recommended)
//!
//! For best performance, use `load_all_spectra()` to read both MS1 and MS2
//! spectra in a single pass:
//! - `load_all_spectra()` - Load both MS1 and MS2 spectra efficiently
//!
//! ## MS1 Spectra for Calibration
//!
//! For MS1 mass calibration (pyXcorrDIA-compatible):
//! - `load_ms1_spectra()` - Load only MS1 spectra (parses file again if MS2 already loaded)
//! - `MS1Index` - Index for efficient nearest-neighbor lookup by RT

mod parser;
mod spectra_cache;

pub use parser::{load_all_spectra, load_ms1_spectra, MS1Index, MzmlReader};
pub use spectra_cache::{load_spectra_cache, save_spectra_cache, spectra_cache_path};
