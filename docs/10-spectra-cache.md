# Spectra Cache (`.spectra.bin`)

## Purpose

Osprey loads mzML files twice during analysis:

1. **Initial search** — parse mzML, run calibration, score all candidates
2. **Post-FDR re-scoring** — reload spectra to re-score consensus and reconciliation targets

Parsing mzML XML is expensive (~3 minutes per Astral file with ~200K MS2 + ~1.2K MS1 spectra). The spectra cache eliminates redundant XML parsing by saving the parsed spectra to a binary file after the first load. The second load deserializes from the binary cache instead, which is nearly instant.

## What is cached

The cache stores **raw, uncalibrated spectra** — the direct output of mzML parsing before any mass calibration is applied. This includes:

- **MS2 spectra** (`Vec<Spectrum>`): scan number, retention time, precursor m/z, isolation window, fragment m/z array, intensity array
- **MS1 index** (`MS1Index`): MS1 spectra sorted by retention time, each with scan number, RT, m/z array, intensity array

Mass calibration (MS2 m/z offsets) is a lightweight operation that is re-applied from the calibration JSON each time spectra are loaded from cache. Storing uncalibrated spectra means the cache remains valid even if calibration parameters change between runs.

## File format

The cache file uses [bincode](https://github.com/bincode-org/bincode) (v1) serialization — a compact binary encoding of Rust's serde data model. The file is written as three sequential bincode segments:

```
┌─────────────────────────────────────────┐
│  CacheHeader (bincode-encoded)          │
│    magic: [u8; 4] = b"OSPC"            │
│    version: u32 = 1                     │
│    source_size: u64                     │
│    source_modified_secs: u64            │
├─────────────────────────────────────────┤
│  Vec<Spectrum> (bincode-encoded)        │
│    For each spectrum:                   │
│      scan_number: u32                   │
│      retention_time: f64                │
│      precursor_mz: f64                  │
│      isolation_window: IsolationWindow  │
│        center: f64                      │
│        lower_offset: f64                │
│        upper_offset: f64                │
│      mzs: Vec<f64>                      │
│      intensities: Vec<f32>              │
├─────────────────────────────────────────┤
│  MS1Index (bincode-encoded)             │
│    spectra: Vec<MS1Spectrum>            │
│      For each MS1 spectrum:             │
│        scan_number: u32                 │
│        retention_time: f64              │
│        mzs: Vec<f64>                    │
│        intensities: Vec<f32>            │
└─────────────────────────────────────────┘
```

Bincode uses little-endian encoding with variable-length integer prefixes for vectors (a `u64` length prefix followed by the packed elements). There is no compression — the data is already numeric arrays that don't compress well.

## File location and naming

The cache file is placed alongside the source mzML with the extension `.spectra.bin`:

```
data/
  sample1.mzML           # Source file
  sample1.spectra.bin    # Cache file (auto-generated)
  sample1.calibration.json
```

## Cache invalidation

The cache header records the source mzML file's **size** (bytes) and **modification time** (seconds since UNIX epoch). On load, these are compared against the current mzML file metadata. The cache is silently discarded if:

- The cache file does not exist
- The magic bytes (`OSPC`) or version number don't match
- The mzML file size has changed
- The mzML file modification time has changed
- Deserialization fails (corrupt file, incompatible version)

When the cache is discarded, Osprey falls back to full mzML parsing with no error. A stale or corrupt cache never causes a failure.

Bumping `CACHE_VERSION` in the source code invalidates all existing caches (e.g., if the `Spectrum` struct layout changes).

## Pipeline integration

```
Search phase (per file):
  1. Try load_spectra_cache(mzml_path)
     → If valid cache exists: use cached spectra (fast)
     → If not: parse mzML via load_all_spectra()
  2. save_spectra_cache(mzml_path, spectra, ms1_index)
  3. Run calibration, scoring, etc.

Post-FDR re-scoring phase (per file):
  1. FileRescoreContext::load() tries load_spectra_cache(mzml_path)
     → Cache was saved in step 2 above, so this hits
  2. Apply MS2 mass calibration from calibration JSON
  3. Re-score consensus + reconciliation targets
```

## Cleanup

Cache files can be safely deleted at any time. Osprey will regenerate them on the next run. They are intermediate files and should not be committed to version control or included in data archives.
