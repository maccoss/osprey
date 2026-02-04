//! Osprey Core - Core types and traits for DIA analysis
//!
//! This crate provides the foundational types, configuration structures,
//! error handling, and traits used throughout the Osprey analysis pipeline.

pub mod config;
pub mod error;
pub mod isotope;
pub mod traits;
pub mod types;

pub use config::*;
pub use error::{OspreyError, Result};
pub use isotope::*;
pub use traits::*;
pub use types::*;
