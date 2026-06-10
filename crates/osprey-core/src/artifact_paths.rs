//! Resolution of the directory for per-file *derived* artifacts.
//!
//! Every per-file artifact path (the scores parquet and its reconciled
//! in-place rewrite, the FDR score sidecars, the `reconciliation.json`
//! boundary file, the calibration JSON, and the `.spectra.bin` cache) is
//! historically derived from the input mzML path, so all of them land next
//! to the input file. On a read-only data tree that forces the caller to
//! copy the (multi-GB) mass-spec files into a writable workdir.
//!
//! These two resolvers decouple the data location from the output location:
//!
//! - [`resolve_output_dir`] returns the directory for non-cache artifacts:
//!   the configured `output_dir` when set, else the input file's own
//!   directory (the historical default).
//! - [`resolve_cache_dir`] returns the directory for the `.spectra.bin`
//!   cache: the configured `cache_dir` when set, else beside the data file
//!   when that directory is writable (probed, not assumed), else the
//!   `output_dir`. When neither `output_dir` nor `cache_dir` is configured
//!   the input directory is returned with **no** writability probe, so a
//!   default run never touches the data directory with a temp file and the
//!   default path stays byte-identical to the historical behavior.
//!
//! Both `output_dir` and `cache_dir` default to `None`, which preserves the
//! historical behavior of writing each artifact in its input file's own
//! directory. `--work-dir` sets both; `--output-dir` / `--cache-dir` set
//! them individually. These map to the OspreySharp `ArtifactPaths`
//! process-wide holder (Track A); the byte layout of the artifacts is
//! unchanged, only the directory moves.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::OnceLock;

/// Directory for a non-cache derived artifact of `input_path`: `output_dir`
/// when set, else the input file's own directory (the historical default).
pub fn resolve_output_dir(input_path: &Path, output_dir: Option<&Path>) -> PathBuf {
    match output_dir {
        Some(dir) => dir.to_path_buf(),
        None => input_dir(input_path),
    }
}

/// Directory for the `.spectra.bin` cache of `input_path`.
///
/// Resolution order: explicit `cache_dir` -> beside the data file when that
/// directory is writable -> `output_dir`. The cache is settings-independent,
/// so a shared `cache_dir` lets many analyses (and the straight-through vs.
/// resume passes) reuse one parse.
pub fn resolve_cache_dir(
    input_path: &Path,
    output_dir: Option<&Path>,
    cache_dir: Option<&Path>,
) -> PathBuf {
    if let Some(dir) = cache_dir {
        return dir.to_path_buf();
    }
    let input = input_dir(input_path);
    // No redirection configured at all: this is exactly the historical
    // location, so skip the writability probe entirely (a default run must
    // not touch the data directory with a temp file).
    let Some(out) = output_dir else {
        return input;
    };
    // output_dir is set but no explicit cache_dir: prefer beside-data for
    // cross-analysis reuse, falling back to output_dir when the data
    // directory is read-only.
    if is_directory_writable(&input) {
        input
    } else {
        out.to_path_buf()
    }
}

fn input_dir(input_path: &Path) -> PathBuf {
    input_path
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_default()
}

// Writability is probed (an ACL check is unreliable cross-platform) and
// memoized per directory so resolution stays cheap when called for every
// input file.
fn writable_cache() -> &'static Mutex<HashMap<PathBuf, bool>> {
    static CACHE: OnceLock<Mutex<HashMap<PathBuf, bool>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

fn is_directory_writable(dir: &Path) -> bool {
    let key = if dir.as_os_str().is_empty() {
        PathBuf::from(".")
    } else {
        dir.to_path_buf()
    };
    if let Some(&cached) = writable_cache().lock().unwrap().get(&key) {
        return cached;
    }
    let result = probe_writable(&key);
    writable_cache().lock().unwrap().insert(key, result);
    result
}

fn probe_writable(dir: &Path) -> bool {
    if !dir.is_dir() {
        return false;
    }
    // create_new avoids clobbering an existing file on a (rare) name collision
    // and makes a concurrent probe in the same directory safe; retry a few times
    // on AlreadyExists before giving up.
    for _ in 0..3 {
        let probe = dir.join(format!(".{}.osprey-wtest", uuid_like()));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&probe)
        {
            Ok(_) => {
                let _ = std::fs::remove_file(&probe);
                return true;
            }
            Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(_) => return false,
        }
    }
    false
}

// A short, collision-resistant token for the probe file name. We don't need
// cryptographic randomness — just enough entropy to avoid clobbering a
// concurrent probe in the same directory.
fn uuid_like() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    format!("{:x}-{:x}", std::process::id(), nanos)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn output_dir_defaults_to_input_dir() {
        let input = Path::new("/data/run/sample1.mzML");
        // No output_dir: historical default = the input file's own dir.
        assert_eq!(resolve_output_dir(input, None), PathBuf::from("/data/run"));
    }

    #[test]
    fn output_dir_override_is_used() {
        let input = Path::new("/data/run/sample1.mzML");
        let out = Path::new("/work/out");
        assert_eq!(
            resolve_output_dir(input, Some(out)),
            PathBuf::from("/work/out")
        );
    }

    #[test]
    fn cache_dir_explicit_wins() {
        let input = Path::new("/data/run/sample1.mzML");
        let out = Path::new("/work/out");
        let cache = Path::new("/work/cache");
        assert_eq!(
            resolve_cache_dir(input, Some(out), Some(cache)),
            PathBuf::from("/work/cache")
        );
        // Explicit cache_dir wins even with no output_dir.
        assert_eq!(
            resolve_cache_dir(input, None, Some(cache)),
            PathBuf::from("/work/cache")
        );
    }

    #[test]
    fn cache_dir_default_is_input_dir_no_probe() {
        // Neither output_dir nor cache_dir configured: the cache lands in the
        // input file's own directory with NO writability probe, byte-identical
        // to the historical default. The directory here does not exist, so a
        // probe would fail — confirming no probe happened.
        let input = Path::new("/data/run/sample1.mzML");
        assert_eq!(
            resolve_cache_dir(input, None, None),
            PathBuf::from("/data/run")
        );
    }

    #[test]
    fn cache_dir_falls_back_to_output_when_input_not_writable() {
        // output_dir set, no cache_dir, and the input directory does not
        // exist (so the writability probe fails): fall back to output_dir.
        let input = Path::new("/no/such/data/dir/sample1.mzML");
        let out = Path::new("/work/out");
        assert_eq!(
            resolve_cache_dir(input, Some(out), None),
            PathBuf::from("/work/out")
        );
    }

    #[test]
    fn cache_dir_prefers_beside_data_when_writable() {
        // output_dir set, no cache_dir, and the input directory IS writable:
        // prefer beside-data for cross-analysis reuse.
        let tmp = tempfile::tempdir().unwrap();
        let input = tmp.path().join("sample1.mzML");
        let out = Path::new("/work/out");
        assert_eq!(
            resolve_cache_dir(&input, Some(out), None),
            tmp.path().to_path_buf()
        );
    }
}
