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

/// Move a file from `src` to `dst` safely across filesystems.
///
/// Tries `rename` first (instant on same filesystem). Falls back to a manual
/// buffered copy using `std::io::copy` (not `std::fs::copy`) because the
/// kernel-level `copy_file_range`/`sendfile` used by `std::fs::copy` can
/// silently produce 0-byte files on CIFS/NFS network mounts.
///
/// After a successful copy, the source file is deleted.
pub fn move_file_safe(src: &std::path::Path, dst: &std::path::Path) -> Result<()> {
    // Check if src and dst are on the same filesystem before attempting rename.
    // On CIFS/NFS, rename from local to network can silently produce 0-byte files.
    let same_fs = src
        .parent()
        .zip(dst.parent())
        .map(|(sp, dp)| {
            // Simple heuristic: same parent directory means same filesystem
            sp == dp
        })
        .unwrap_or(false);

    if same_fs && std::fs::rename(src, dst).is_ok() {
        return Ok(());
    }

    // Cross-filesystem: manual buffered copy with verification
    let mut reader = std::io::BufReader::new(std::fs::File::open(src).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to open '{}' for copy: {}",
            src.display(),
            e
        ))
    })?);
    let mut writer = std::io::BufWriter::new(std::fs::File::create(dst).map_err(|e| {
        OspreyError::OutputError(format!("Failed to create '{}': {}", dst.display(), e))
    })?);
    std::io::copy(&mut reader, &mut writer).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to copy '{}' to '{}': {}",
            src.display(),
            dst.display(),
            e
        ))
    })?;
    // Flush to ensure all data is written to the network filesystem
    use std::io::Write;
    writer.flush().map_err(|e| {
        OspreyError::OutputError(format!("Failed to flush '{}': {}", dst.display(), e))
    })?;
    drop(writer);

    // Verify the copy succeeded (file size matches)
    let src_size = std::fs::metadata(src).map(|m| m.len()).unwrap_or(0);
    let dst_size = std::fs::metadata(dst).map(|m| m.len()).unwrap_or(0);
    if dst_size != src_size {
        let _ = std::fs::remove_file(dst);
        return Err(OspreyError::OutputError(format!(
            "Copy verification failed: '{}' is {} bytes but '{}' is {} bytes",
            src.display(),
            src_size,
            dst.display(),
            dst_size
        )));
    }

    // Delete source
    let _ = std::fs::remove_file(src);
    Ok(())
}
