//! Output writers for Osprey results
//!
//! This module provides writers for various output formats:
//! - blib format for Skyline integration

mod blib;

pub use blib::{unimod_id_to_mass, BlibWriter};
