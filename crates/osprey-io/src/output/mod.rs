//! Output writers for Osprey results
//!
//! This module provides writers for various output formats:
//! - blib format for Skyline integration
//! - FDRBench peptide/precursor input TSV

mod blib;
pub mod fdrbench;

pub use blib::{unimod_id_to_mass, BlibWriter};
pub use fdrbench::write_fdrbench_peptide_input;
