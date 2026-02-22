//! Binary cache for parsed mzML spectra
//!
//! Saves parsed spectra (MS2 + MS1) to a bincode file alongside the mzML,
//! so the second load (post-FDR re-scoring) can skip full mzML parsing.

use std::fs;
use std::io::{BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use osprey_core::types::Spectrum;

use crate::mzml::MS1Index;

/// Magic bytes + version for cache format validation
const CACHE_MAGIC: &[u8; 4] = b"OSPC";
const CACHE_VERSION: u32 = 1;

/// Header stored at the beginning of the cache file
#[derive(serde::Serialize, serde::Deserialize)]
struct CacheHeader {
    magic: [u8; 4],
    version: u32,
    /// mzML file size in bytes (for invalidation)
    source_size: u64,
    /// mzML file modification time as duration since UNIX_EPOCH
    source_modified_secs: u64,
}

/// Returns the cache file path for a given mzML file (sibling `.spectra.bin`)
pub fn cache_path_for(mzml_path: &Path) -> PathBuf {
    mzml_path.with_extension("spectra.bin")
}

/// Save parsed spectra to a binary cache file.
///
/// The cache is written alongside the mzML file with a `.spectra.bin` extension.
/// Returns Ok(()) on success, or an error (which callers should log and ignore).
pub fn save_spectra_cache(
    mzml_path: &Path,
    spectra: &[Spectrum],
    ms1_index: &MS1Index,
) -> Result<(), String> {
    let cache_path = cache_path_for(mzml_path);

    // Get source file metadata for invalidation
    let metadata =
        fs::metadata(mzml_path).map_err(|e| format!("Failed to read mzML metadata: {}", e))?;
    let source_size = metadata.len();
    let source_modified_secs = metadata
        .modified()
        .unwrap_or(SystemTime::UNIX_EPOCH)
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let header = CacheHeader {
        magic: *CACHE_MAGIC,
        version: CACHE_VERSION,
        source_size,
        source_modified_secs,
    };

    let file = fs::File::create(&cache_path).map_err(|e| {
        format!(
            "Failed to create cache file '{}': {}",
            cache_path.display(),
            e
        )
    })?;
    let mut writer = BufWriter::new(file);

    bincode::serialize_into(&mut writer, &header)
        .map_err(|e| format!("Failed to write cache header: {}", e))?;
    bincode::serialize_into(&mut writer, &spectra)
        .map_err(|e| format!("Failed to write MS2 spectra to cache: {}", e))?;
    bincode::serialize_into(&mut writer, &ms1_index)
        .map_err(|e| format!("Failed to write MS1 index to cache: {}", e))?;

    Ok(())
}

/// Try to load spectra from cache. Returns None if cache is missing, stale, or corrupt.
pub fn load_spectra_cache(mzml_path: &Path) -> Option<(Vec<Spectrum>, MS1Index)> {
    let cache_path = cache_path_for(mzml_path);

    // Check cache file exists
    let cache_file = fs::File::open(&cache_path).ok()?;
    let mut reader = BufReader::new(cache_file);

    // Read and validate header
    let header: CacheHeader = bincode::deserialize_from(&mut reader).ok()?;

    if &header.magic != CACHE_MAGIC || header.version != CACHE_VERSION {
        log::debug!("Cache magic/version mismatch, ignoring cache");
        return None;
    }

    // Validate against source file
    let metadata = fs::metadata(mzml_path).ok()?;
    let source_size = metadata.len();
    let source_modified_secs = metadata
        .modified()
        .unwrap_or(SystemTime::UNIX_EPOCH)
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    if header.source_size != source_size || header.source_modified_secs != source_modified_secs {
        log::debug!(
            "Cache stale (size: {} vs {}, mtime: {} vs {}), ignoring",
            header.source_size,
            source_size,
            header.source_modified_secs,
            source_modified_secs
        );
        return None;
    }

    // Deserialize spectra
    let spectra: Vec<Spectrum> = bincode::deserialize_from(&mut reader).ok()?;
    let ms1_index: MS1Index = bincode::deserialize_from(&mut reader).ok()?;

    Some((spectra, ms1_index))
}
