//! mzML file parsing
//!
//! This module provides parsing of mzML files using the mzdata crate,
//! converting to Osprey's internal spectrum representation.

mod parser;

pub use parser::MzmlReader;
