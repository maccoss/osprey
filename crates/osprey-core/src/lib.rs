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
/// First tries `std::fs::copy` (kernel-optimized). If the copied byte count doesn't
/// match the source size, falls back to manual buffered `std::io::copy`. This handles
/// CIFS/NFS filesystems where `copy_file_range` may silently produce truncated files.
pub fn copy_and_verify(src: &std::path::Path, dst: &std::path::Path) -> Result<()> {
    let src_size = std::fs::metadata(src)
        .map(|m| m.len())
        .map_err(|e| OspreyError::OutputError(format!("Cannot stat '{}': {}", src.display(), e)))?;

    if src_size == 0 {
        // Source is empty — this shouldn't happen for valid cache files
        let _ = std::fs::remove_file(src);
        return Err(OspreyError::OutputError(format!(
            "Source file '{}' is 0 bytes — refusing to copy",
            src.display()
        )));
    }

    // Try std::fs::copy first (kernel-optimized, fast on local filesystems)
    let bytes_copied = std::fs::copy(src, dst).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to copy '{}' ({} bytes) to '{}': {}",
            src.display(),
            src_size,
            dst.display(),
            e
        ))
    })?;

    if bytes_copied == src_size {
        // Success — clean up source
        let _ = std::fs::remove_file(src);
        return Ok(());
    }

    // std::fs::copy returned wrong byte count — fall back to manual buffered copy.
    // This handles CIFS where copy_file_range may silently truncate large files.
    eprintln!(
        "WARNING: std::fs::copy returned {} bytes but source is {} bytes for '{}' — retrying with buffered copy",
        bytes_copied, src_size, dst.display()
    );

    let mut reader = std::io::BufReader::with_capacity(
        8 * 1024 * 1024, // 8 MB buffer for large files
        std::fs::File::open(src).map_err(|e| {
            OspreyError::OutputError(format!("Failed to reopen '{}': {}", src.display(), e))
        })?,
    );
    let mut writer = std::io::BufWriter::with_capacity(
        8 * 1024 * 1024,
        std::fs::File::create(dst).map_err(|e| {
            OspreyError::OutputError(format!("Failed to recreate '{}': {}", dst.display(), e))
        })?,
    );
    let manual_bytes = std::io::copy(&mut reader, &mut writer).map_err(|e| {
        OspreyError::OutputError(format!(
            "Buffered copy failed '{}' to '{}': {}",
            src.display(),
            dst.display(),
            e
        ))
    })?;
    use std::io::Write;
    writer.flush().map_err(|e| {
        OspreyError::OutputError(format!("Failed to flush '{}': {}", dst.display(), e))
    })?;
    drop(writer);

    if manual_bytes != src_size {
        let _ = std::fs::remove_file(dst);
        return Err(OspreyError::OutputError(format!(
            "Buffered copy also failed: wrote {} of {} bytes to '{}'",
            manual_bytes,
            src_size,
            dst.display()
        )));
    }

    // Success — clean up source
    let _ = std::fs::remove_file(src);
    Ok(())
}
