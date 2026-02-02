//! Output writers for Osprey results
//!
//! This module provides writers for various output formats:
//! - TSV reports for human-readable results
//! - blib format for Skyline integration

mod blib;
mod report;

pub use blib::BlibWriter;
pub use report::{write_coefficient_summary, ReportWriter};
