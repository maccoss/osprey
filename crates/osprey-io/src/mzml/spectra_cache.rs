//! Binary spectra cache for fast re-loading.
//!
//! After parsing mzML (which takes ~5 minutes for large Astral files),
//! this module writes a compact `.spectra.bin` file that can be reloaded
//! in seconds. Used by the post-FDR re-scoring phase to avoid re-parsing mzML.
//!
//! Format (little-endian):
//! ```text
//! [magic: 8 bytes "OSPRSPC\0"]
//! [version: u32]
//! [source_size: u64]    source file length, or 0 when unknown
//! [source_mtime: i64]   source last-write time, Unix milliseconds UTC,
//!                       or 0 when unknown
//! [n_ms2: u32]
//! [n_ms1: u32]
//! For each MS2 spectrum:
//!   [scan_number: u32] [retention_time: f64] [precursor_mz: f64]
//!   [iso_center: f64] [iso_lower: f64] [iso_upper: f64]
//!   [n_peaks: u32] [mzs: f64 × n_peaks] [intensities: f32 × n_peaks]
//! For each MS1 spectrum:
//!   [scan_number: u32] [retention_time: f64]
//!   [n_peaks: u32] [mzs: f64 × n_peaks] [intensities: f32 × n_peaks]
//! ```
//!
//! The `source_size` / `source_mtime` fingerprint lets a cache that lives
//! apart from its data file (`--cache-dir` / `--output-dir`) be rejected
//! when the source mzML changes. The fingerprint byte layout is shared with
//! OspreySharp's `SpectraCache` (Track A): the m-time is Unix milliseconds
//! UTC specifically so both tools derive identical values for the same file
//! and can share one cache across the cross-impl parity gate.

use osprey_core::{IsolationWindow, MS1Spectrum, OspreyError, Result, Spectrum};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::time::UNIX_EPOCH;

use super::MS1Index;

const MAGIC: &[u8; 8] = b"OSPRSPC\0";
// VERSION 2 (2026-05-09): mzML load now sorts non-monotonic centroids before
// caching, so caches written by VERSION 1 may contain unsorted peaks that
// produce undefined-behavior divergence in fragment matching. Old caches are
// invalidated on this version bump so the fresh load path (which sorts)
// re-populates them.
// VERSION 3 (2026-06-09): header now carries the source file's size +
// last-write-time (Unix ms) so a cache that lives apart from its data file
// (--cache-dir / --output-dir) is rejected when the source changes. Old
// caches re-populate on this bump. Layout matches OspreySharp's VERSION 3.
const VERSION: u32 = 3;

/// Compute the source-file fingerprint written into the cache header:
/// `(size, mtime_ms)` where `mtime_ms` is the last-write time as Unix
/// milliseconds UTC. Returns `(0, 0)` when the source is `None`, missing,
/// or cannot be stat'd — the load path treats a zero size as "no
/// fingerprint recorded" and trusts the cache.
fn source_fingerprint(source_path: Option<&Path>) -> (u64, i64) {
    let Some(path) = source_path else {
        return (0, 0);
    };
    let Ok(meta) = std::fs::metadata(path) else {
        return (0, 0);
    };
    let size = meta.len();
    let mtime_ms = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| d.as_millis() as i64)
        .unwrap_or(0);
    (size, mtime_ms)
}

/// Save spectra to a binary cache file for fast reload.
///
/// `source_path`, when supplied, is stat'd to record a size + mtime
/// fingerprint in the header so [`load_spectra_cache`] can reject the cache
/// if the source mzML later changes. Pass `None` (or a path that can't be
/// stat'd) to write a `(0, 0)` fingerprint, which the loader treats as
/// "trust the cache".
pub fn save_spectra_cache(
    path: &Path,
    spectra: &[Spectrum],
    ms1_index: &MS1Index,
    source_path: Option<&Path>,
) -> Result<()> {
    // Write to local temp dir first, then copy to final destination.
    let tmp_path = std::env::temp_dir().join(format!(
        "osprey_{}_{}",
        std::process::id(),
        path.file_name()
            .unwrap_or(std::ffi::OsStr::new("spectra.bin"))
            .to_string_lossy()
    ));
    let file = File::create(&tmp_path).map_err(|e| {
        OspreyError::IoError(std::io::Error::other(format!(
            "Failed to create spectra cache '{}': {}",
            tmp_path.display(),
            e
        )))
    })?;
    let mut w = BufWriter::new(file);

    let ms1_spectra = ms1_index.spectra();

    let (source_size, source_mtime_ms) = source_fingerprint(source_path);

    // Header
    w.write_all(MAGIC).map_err(write_err)?;
    w.write_all(&VERSION.to_le_bytes()).map_err(write_err)?;
    w.write_all(&source_size.to_le_bytes()).map_err(write_err)?;
    w.write_all(&source_mtime_ms.to_le_bytes())
        .map_err(write_err)?;
    w.write_all(&(spectra.len() as u32).to_le_bytes())
        .map_err(write_err)?;
    w.write_all(&(ms1_spectra.len() as u32).to_le_bytes())
        .map_err(write_err)?;

    // MS2 spectra
    for s in spectra {
        w.write_all(&s.scan_number.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.retention_time.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.precursor_mz.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.isolation_window.center.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.isolation_window.lower_offset.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.isolation_window.upper_offset.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&(s.mzs.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        // Write mzs as raw bytes
        let mz_bytes: &[u8] =
            unsafe { std::slice::from_raw_parts(s.mzs.as_ptr() as *const u8, s.mzs.len() * 8) };
        w.write_all(mz_bytes).map_err(write_err)?;
        // Write intensities as raw bytes
        let int_bytes: &[u8] = unsafe {
            std::slice::from_raw_parts(s.intensities.as_ptr() as *const u8, s.intensities.len() * 4)
        };
        w.write_all(int_bytes).map_err(write_err)?;
    }

    // MS1 spectra
    for s in ms1_spectra {
        w.write_all(&s.scan_number.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.retention_time.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&(s.mzs.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        let mz_bytes: &[u8] =
            unsafe { std::slice::from_raw_parts(s.mzs.as_ptr() as *const u8, s.mzs.len() * 8) };
        w.write_all(mz_bytes).map_err(write_err)?;
        let int_bytes: &[u8] = unsafe {
            std::slice::from_raw_parts(s.intensities.as_ptr() as *const u8, s.intensities.len() * 4)
        };
        w.write_all(int_bytes).map_err(write_err)?;
    }

    w.flush().map_err(write_err)?;
    drop(w);

    osprey_core::copy_and_verify(&tmp_path, path)?;
    Ok(())
}

/// Load spectra from a binary cache file.
///
/// `source_path`, when supplied, is compared against the size + mtime
/// fingerprint stored in the header: a cache whose source mzML changed
/// since it was written is rejected (returns an error so the caller falls
/// back to re-parsing). The check is skipped when the cache recorded no
/// fingerprint (stored size == 0) or the source is unavailable for
/// comparison (e.g. a resume run whose mzML is not beside the cache) — a
/// within-run cache is trusted in that case.
pub fn load_spectra_cache(
    path: &Path,
    source_path: Option<&Path>,
) -> Result<(Vec<Spectrum>, MS1Index)> {
    let file = File::open(path).map_err(|e| {
        OspreyError::IoError(std::io::Error::other(format!(
            "Failed to open spectra cache '{}': {}",
            path.display(),
            e
        )))
    })?;
    let mut r = BufReader::new(file);

    // Validate header
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic).map_err(read_err)?;
    if &magic != MAGIC {
        return Err(OspreyError::IoError(std::io::Error::other(
            "Invalid spectra cache magic bytes",
        )));
    }

    let version = read_u32(&mut r)?;
    if version != VERSION {
        return Err(OspreyError::IoError(std::io::Error::other(format!(
            "Unsupported spectra cache version: {} (expected {})",
            version, VERSION
        ))));
    }

    // Validate source fingerprint: reject a cache whose data file changed
    // since it was written. Skipped when the cache recorded no fingerprint
    // (stored_size == 0) or the source is unavailable for comparison.
    let stored_size = read_u64(&mut r)?;
    let stored_mtime_ms = read_i64(&mut r)?;
    if stored_size != 0 {
        if let Some(src) = source_path {
            let (actual_size, actual_mtime_ms) = source_fingerprint(Some(src));
            if actual_size != 0
                && (actual_size != stored_size || actual_mtime_ms != stored_mtime_ms)
            {
                return Err(OspreyError::IoError(std::io::Error::other(format!(
                    "Spectra cache '{}' is stale: source file '{}' changed \
                     (size {} -> {}, mtime {} -> {} ms)",
                    path.display(),
                    src.display(),
                    stored_size,
                    actual_size,
                    stored_mtime_ms,
                    actual_mtime_ms
                ))));
            }
        }
    }

    let n_ms2 = read_u32(&mut r)? as usize;
    let n_ms1 = read_u32(&mut r)? as usize;

    // Read MS2 spectra
    let mut spectra = Vec::with_capacity(n_ms2);
    for _ in 0..n_ms2 {
        let scan_number = read_u32(&mut r)?;
        let retention_time = read_f64(&mut r)?;
        let precursor_mz = read_f64(&mut r)?;
        let iso_center = read_f64(&mut r)?;
        let iso_lower = read_f64(&mut r)?;
        let iso_upper = read_f64(&mut r)?;
        let n_peaks = read_u32(&mut r)? as usize;

        let mut mzs = vec![0f64; n_peaks];
        let mz_bytes: &mut [u8] =
            unsafe { std::slice::from_raw_parts_mut(mzs.as_mut_ptr() as *mut u8, n_peaks * 8) };
        r.read_exact(mz_bytes).map_err(read_err)?;

        let mut intensities = vec![0f32; n_peaks];
        let int_bytes: &mut [u8] = unsafe {
            std::slice::from_raw_parts_mut(intensities.as_mut_ptr() as *mut u8, n_peaks * 4)
        };
        r.read_exact(int_bytes).map_err(read_err)?;

        spectra.push(Spectrum {
            scan_number,
            retention_time,
            precursor_mz,
            isolation_window: IsolationWindow::new(iso_center, iso_lower, iso_upper),
            mzs,
            intensities,
        });
    }

    // Read MS1 spectra
    let mut ms1_spectra = Vec::with_capacity(n_ms1);
    for _ in 0..n_ms1 {
        let scan_number = read_u32(&mut r)?;
        let retention_time = read_f64(&mut r)?;
        let n_peaks = read_u32(&mut r)? as usize;

        let mut mzs = vec![0f64; n_peaks];
        let mz_bytes: &mut [u8] =
            unsafe { std::slice::from_raw_parts_mut(mzs.as_mut_ptr() as *mut u8, n_peaks * 8) };
        r.read_exact(mz_bytes).map_err(read_err)?;

        let mut intensities = vec![0f32; n_peaks];
        let int_bytes: &mut [u8] = unsafe {
            std::slice::from_raw_parts_mut(intensities.as_mut_ptr() as *mut u8, n_peaks * 4)
        };
        r.read_exact(int_bytes).map_err(read_err)?;

        ms1_spectra.push(MS1Spectrum {
            scan_number,
            retention_time,
            mzs,
            intensities,
        });
    }

    Ok((spectra, MS1Index::new(ms1_spectra)))
}

/// Get the spectra cache filename for a given input file, e.g.
/// `/data/sample.mzML` -> `sample.spectra.bin`. Matches the historical
/// `with_extension("spectra.bin")` result (the stem before the final
/// extension) so the produced filename is unchanged; only the directory is
/// redirected by [`spectra_cache_path_in`].
fn spectra_cache_filename(input_file: &Path) -> std::path::PathBuf {
    input_file
        .with_extension("spectra.bin")
        .file_name()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::path::PathBuf::from("spectra.bin"))
}

/// Get the spectra cache path for a given input file, written next to the
/// input (the historical default).
pub fn spectra_cache_path(input_file: &Path) -> std::path::PathBuf {
    input_file.with_extension("spectra.bin")
}

/// Get the spectra cache path for a given input file in an explicit
/// directory. The filename is unchanged (`{stem}.spectra.bin`); only the
/// directory is `cache_dir`. Pass the directory resolved by
/// [`osprey_core::resolve_cache_dir`].
pub fn spectra_cache_path_in(input_file: &Path, cache_dir: &Path) -> std::path::PathBuf {
    cache_dir.join(spectra_cache_filename(input_file))
}

fn write_err(e: std::io::Error) -> OspreyError {
    OspreyError::IoError(e)
}

fn read_err(e: std::io::Error) -> OspreyError {
    OspreyError::IoError(e)
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_f64<R: Read>(r: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(f64::from_le_bytes(buf))
}

fn read_u64<R: Read>(r: &mut R) -> Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_i64<R: Read>(r: &mut R) -> Result<i64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(i64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_spectra_cache_round_trip() {
        let spectra = vec![
            Spectrum {
                scan_number: 1,
                retention_time: 10.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::new(500.0, 1.5, 1.5),
                mzs: vec![100.0, 200.0, 300.0],
                intensities: vec![1000.0, 2000.0, 500.0],
            },
            Spectrum {
                scan_number: 2,
                retention_time: 11.0,
                precursor_mz: 600.0,
                isolation_window: IsolationWindow::new(600.0, 2.0, 3.0),
                mzs: vec![150.0, 250.0],
                intensities: vec![800.0, 1200.0],
            },
        ];

        let ms1_spectra = vec![MS1Spectrum {
            scan_number: 0,
            retention_time: 10.0,
            mzs: vec![400.0, 500.0, 600.0],
            intensities: vec![5000.0, 3000.0, 1000.0],
        }];
        let ms1_index = MS1Index::new(ms1_spectra);

        let tmp = NamedTempFile::new().unwrap();
        save_spectra_cache(tmp.path(), &spectra, &ms1_index, None).unwrap();
        let (loaded_spectra, loaded_ms1) = load_spectra_cache(tmp.path(), None).unwrap();

        assert_eq!(loaded_spectra.len(), 2);
        assert_eq!(loaded_spectra[0].scan_number, 1);
        assert!((loaded_spectra[0].retention_time - 10.5).abs() < 1e-10);
        assert_eq!(loaded_spectra[0].mzs, vec![100.0, 200.0, 300.0]);
        assert_eq!(loaded_spectra[0].intensities, vec![1000.0, 2000.0, 500.0]);
        assert!((loaded_spectra[1].isolation_window.lower_offset - 2.0).abs() < 1e-10);
        assert!((loaded_spectra[1].isolation_window.upper_offset - 3.0).abs() < 1e-10);

        assert_eq!(loaded_ms1.len(), 1);
    }

    #[test]
    fn test_empty_spectra_cache() {
        let tmp = NamedTempFile::new().unwrap();
        let ms1_index = MS1Index::new(vec![]);
        save_spectra_cache(tmp.path(), &[], &ms1_index, None).unwrap();
        let (spectra, ms1) = load_spectra_cache(tmp.path(), None).unwrap();
        assert!(spectra.is_empty());
        assert!(ms1.is_empty());
    }

    #[test]
    fn test_fingerprint_accepts_unchanged_source() {
        // Write a cache fingerprinted against a real source file, then load
        // it back with the same (unmodified) source: the fingerprint matches
        // and the cache is accepted.
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join("unchanged.mzML");
        std::fs::write(&source, b"raw mzml bytes").unwrap();
        let cache = dir.path().join("unchanged.spectra.bin");

        let ms1_index = MS1Index::new(vec![]);
        save_spectra_cache(&cache, &[], &ms1_index, Some(&source)).unwrap();
        let loaded = load_spectra_cache(&cache, Some(&source));
        assert!(loaded.is_ok(), "unchanged source should be accepted");
    }

    #[test]
    fn test_fingerprint_rejects_changed_source() {
        // Write a cache fingerprinted against a source file, then mutate the
        // source so its size (and mtime) differ: the fingerprint mismatch
        // makes the loader reject the cache.
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join("changed.mzML");
        std::fs::write(&source, b"raw mzml bytes").unwrap();
        let cache = dir.path().join("changed.spectra.bin");

        let ms1_index = MS1Index::new(vec![]);
        save_spectra_cache(&cache, &[], &ms1_index, Some(&source)).unwrap();

        // Change the source size; the fingerprint no longer matches.
        std::fs::write(&source, b"different, longer raw mzml bytes").unwrap();
        let loaded = load_spectra_cache(&cache, Some(&source));
        assert!(loaded.is_err(), "changed source should be rejected");
    }

    #[test]
    fn test_fingerprint_skipped_when_source_unavailable() {
        // A cache written with a fingerprint is still trusted when no source
        // is supplied for comparison (e.g. a resume run whose mzML is not
        // beside the cache).
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join("unavailable.mzML");
        std::fs::write(&source, b"raw mzml bytes").unwrap();
        let cache = dir.path().join("unavailable.spectra.bin");

        let ms1_index = MS1Index::new(vec![]);
        save_spectra_cache(&cache, &[], &ms1_index, Some(&source)).unwrap();

        // Delete the source, then load with no source path: trusted.
        std::fs::remove_file(&source).unwrap();
        let loaded = load_spectra_cache(&cache, None);
        assert!(loaded.is_ok(), "no source for comparison should be trusted");
    }

    #[test]
    fn test_cache_path_in_redirects_dir_keeps_filename() {
        let input = Path::new("/data/run/sample.mzML");
        let cache_dir = Path::new("/work/cache");
        assert_eq!(
            spectra_cache_path_in(input, cache_dir),
            Path::new("/work/cache/sample.spectra.bin")
        );
        // Default helper is unchanged (beside the input).
        assert_eq!(
            spectra_cache_path(input),
            Path::new("/data/run/sample.spectra.bin")
        );
    }
}
