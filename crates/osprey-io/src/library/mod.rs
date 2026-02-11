//! Spectral library loading
//!
//! This module provides parsers for various spectral library formats:
//! - DIA-NN TSV format
//! - BiblioSpec blib format
//! - EncyclopeDIA elib format

mod blib;
mod diann;
mod elib;

pub use blib::BlibLoader;
pub use diann::DiannTsvLoader;
pub use elib::ElibLoader;

use osprey_core::{LibraryEntry, LibrarySource, OspreyError, Result};

/// Load a library from any supported source
pub fn load_library(source: &LibrarySource) -> Result<Vec<LibraryEntry>> {
    match source {
        LibrarySource::DiannTsv(path) => {
            let loader = DiannTsvLoader::new();
            loader.load(path)
        }
        LibrarySource::Blib(path) => {
            let loader = BlibLoader::new();
            loader.load(path)
        }
        LibrarySource::Elib(path) => {
            let loader = ElibLoader::new();
            loader.load(path)
        }
        LibrarySource::SkylineDocument(path) => {
            // TODO: Implement Skyline document extraction
            Err(OspreyError::LibraryLoadError(format!(
                "Skyline document format not yet implemented: {}",
                path.display()
            )))
        }
    }
}
