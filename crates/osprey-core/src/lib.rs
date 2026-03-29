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

/// Copy a local temp file to its final destination, verify, and delete the local copy.
///
/// Writes to network filesystems (NAS/CIFS) should first write to a local temp
/// directory, then call this to copy to the final path. This avoids corrupt files
/// from network write issues. Uses `std::fs::copy` and verifies file size matches.
pub fn copy_and_verify(src: &std::path::Path, dst: &std::path::Path) -> Result<()> {
    let src_size = std::fs::metadata(src)
        .map(|m| m.len())
        .map_err(|e| OspreyError::OutputError(format!("Cannot stat '{}': {}", src.display(), e)))?;

    if src_size == 0 {
        return Err(OspreyError::OutputError(format!(
            "Source file '{}' is 0 bytes — refusing to copy",
            src.display()
        )));
    }

    std::fs::copy(src, dst).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to copy '{}' ({} bytes) to '{}': {}",
            src.display(),
            src_size,
            dst.display(),
            e
        ))
    })?;

    // Verify destination size matches source
    let dst_size = std::fs::metadata(dst).map(|m| m.len()).unwrap_or(0);
    if dst_size != src_size {
        let _ = std::fs::remove_file(dst);
        return Err(OspreyError::OutputError(format!(
            "Copy verification failed: source '{}' is {} bytes but destination '{}' is {} bytes",
            src.display(),
            src_size,
            dst.display(),
            dst_size
        )));
    }

    // Clean up local temp file
    let _ = std::fs::remove_file(src);
    Ok(())
}
