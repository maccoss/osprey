# 12 - Intermediate File Formats

Osprey produces several intermediate files during analysis. These files enable fast re-loading, disk-backed memory management, and calibration reuse across runs. All intermediate files are written alongside the input mzML files.

## File Overview

| File Pattern | Format | Purpose | Typical Size | Reusable? |
|---|---|---|---|---|
| `*.calibration.json` | JSON | RT and mass calibration parameters | 5–50 KB | Yes (across runs with same LC-MS setup) |
| `*.spectra.bin` | Custom binary | Decoded spectra for fast reload | 2–5 GB (larger than mzML) | Auto-created on first run |
| `*.scores.parquet` | Apache Parquet (ZSTD) | Full scored entries with features, fragments, CWT candidates | 100 MB–2 GB | Per-run, deleted after pipeline |
| `*.blib` | SQLite (BiblioSpec) | Final output for Skyline | 50–500 MB | N/A (output) |

---

## 1. Calibration JSON (`*.calibration.json`)

**Source**: `crates/osprey-chromatography/src/calibration/io.rs`

Stores the calibration parameters computed during Phase 2 (auto-calibration). If a calibration file already exists for an input file, Osprey skips the calibration phase and loads it directly.

### Structure

```json
{
  "metadata": {
    "num_confident_peptides": 4523,
    "num_sampled_precursors": 10000,
    "calibration_successful": true,
    "timestamp": "2025-01-15T10:30:00Z",
    "isolation_scheme": {
      "num_windows": 166,
      "mz_min": 401.5,
      "mz_max": 898.5,
      "typical_width": 3.0,
      "uniform_width": true,
      "windows": [[401.5, 3.0], [404.5, 3.0], "..."]
    }
  },
  "ms1_calibration": {
    "mean": -0.42,
    "median": -0.38,
    "sd": 1.85,
    "count": 4523,
    "unit": "ppm",
    "adjusted_tolerance": 5.97,
    "window_halfwidth_multiplier": 3.0,
    "histogram": { "bin_edges": [], "counts": [], "bin_width": 0.5 },
    "calibrated": true
  },
  "ms2_calibration": { "...same structure as ms1..." },
  "rt_calibration": {
    "method": "Loess",
    "residual_sd": 0.15,
    "n_points": 4523,
    "r_squared": 0.9987,
    "model_params": { "...LOESS knots and coefficients..." },
    "mad": 0.08
  },
  "second_pass_rt": null
}
```

### Key Types

| Field | Type | Description |
|---|---|---|
| `metadata` | `CalibrationMetadata` | Process info: peptide counts, timestamp, isolation scheme |
| `ms1_calibration` | `MzCalibration` | Precursor m/z error stats and adjusted tolerance |
| `ms2_calibration` | `MzCalibration` | Fragment m/z error stats and adjusted tolerance |
| `rt_calibration` | `RTCalibrationParams` | RT model (LOESS), residual SD, R² |
| `second_pass_rt` | `Option<RTCalibrationParams>` | Cross-run consensus RT recalibration (populated after reconciliation) |

### Usage

- **Reuse**: Delete the `.calibration.json` to force recalibration. Keep it to skip calibration on re-runs.
- **Visualization**: `python scripts/evaluate_calibration.py sample.calibration.json --output report.html`

---

## 2. Binary Spectra Cache (`*.spectra.bin`)

**Source**: `crates/osprey-io/src/mzml/spectra_cache.rs`

A raw binary dump of all decoded MS1 and MS2 spectra from an mzML file. Created on first parse and reloaded during reconciliation re-scoring to avoid re-parsing mzML (which takes minutes for large files).

### Format (little-endian)

```
Header (20 bytes):
  [magic:   8 bytes  "OSPRSPC\0"]
  [version: u32      currently 1]
  [n_ms2:   u32      number of MS2 spectra]
  [n_ms1:   u32      number of MS1 spectra]

For each MS2 spectrum:
  [scan_number:  u32]
  [retention_time: f64]
  [precursor_mz:   f64]
  [iso_center:     f64]
  [iso_lower:      f64]
  [iso_upper:      f64]
  [n_peaks:        u32]
  [mzs:            f64 × n_peaks]
  [intensities:    f32 × n_peaks]

For each MS1 spectrum:
  [scan_number:    u32]
  [retention_time: f64]
  [n_peaks:        u32]
  [mzs:            f64 × n_peaks]
  [intensities:    f32 × n_peaks]
```

### Per-peak storage: 12 bytes

Each peak is stored as `f64` m/z (8 bytes) + `f32` intensity (4 bytes) = **12 bytes/peak**.

### Why it's larger than mzML but faster

| Property | mzML | .spectra.bin |
|---|---|---|
| Peak encoding | Base64 + zlib compressed | Raw f64/f32 arrays |
| Bytes per peak | ~4–6 (compressed) | 12 (uncompressed) |
| Parse overhead | XML parser + Base64 decode + zlib inflate | Sequential `read_exact()` into pre-allocated `Vec` |
| Typical load time | 3–5 min (large Astral file) | 5–15 sec |
| Typical file size | 1–2 GB | 2–5 GB |

The trade-off is deliberate: disk space is cheap, but re-parsing mzML during reconciliation would add minutes per file in a multi-file experiment.

### When created/used

- **Created**: First time an mzML file is processed (end of Phase 3 scoring)
- **Used**: Reconciliation re-scoring (Phase 4) reloads spectra from cache instead of re-parsing mzML
- **Safe to delete**: Yes — will be recreated on next run

---

## 3. Scores Parquet Cache (`*.scores.parquet`)

**Source**: `crates/osprey/src/pipeline.rs` (`write_scores_parquet()`, line ~870)

Per-file cache of all scored entries (targets + decoys) with full feature vectors, fragment data, and CWT peak candidates. Written after Phase 3 scoring, read selectively during FDR and output phases.

### Schema (~40 columns)

**Identity columns (8)**:

| Column | Type | Description |
|---|---|---|
| `entry_id` | UInt32 | Library entry ID (bit 31 = decoy flag) |
| `is_decoy` | Boolean | True for decoy entries |
| `sequence` | Utf8 | Bare peptide sequence |
| `modified_sequence` | Utf8 | Peptide with modifications (e.g., `C[+57.0]`) |
| `charge` | UInt8 | Precursor charge state |
| `precursor_mz` | Float64 | Precursor m/z |
| `protein_ids` | Utf8 (nullable) | Semicolon-separated protein accessions |
| `file_name` | Utf8 | Source mzML file name |

**Boundary columns (5)**:

| Column | Type | Description |
|---|---|---|
| `scan_number` | UInt32 | Apex scan number |
| `apex_rt` | Float64 | Apex retention time (minutes) |
| `start_rt` | Float64 | Peak start RT |
| `end_rt` | Float64 | Peak end RT |
| `bounds_area` | Float64 | Integrated peak area |
| `bounds_snr` | Float64 | Signal-to-noise ratio |

**Binary columns (5)** — variable-length arrays serialized as little-endian bytes:

| Column | Type | Encoding |
|---|---|---|
| `cwt_candidates` | Binary | `[u32 count][48 bytes × count]` — each candidate: apex_rt, start_rt, end_rt, area, snr, coelution_score (all f64) |
| `fragment_mzs` | Binary | `[f64 × n_fragments]` — library fragment m/z values |
| `fragment_intensities` | Binary | `[f32 × n_fragments]` — library fragment intensities |
| `reference_xic_rts` | Binary | `[f64 × n_points]` — reference XIC retention times |
| `reference_xic_intensities` | Binary | `[f64 × n_points]` — reference XIC intensities |

**Feature columns (21)** — all Float64, one per PIN feature:

| Group | Features |
|---|---|
| Pairwise coelution (3) | `fragment_coelution_sum`, `fragment_coelution_max`, `n_coeluting_fragments` |
| Peak shape (3) | `peak_apex`, `peak_area`, `peak_sharpness` |
| Spectral at apex (3) | `xcorr`, `consecutive_ions`, `explained_intensity` |
| Mass accuracy (2) | `mass_accuracy_deviation_mean`, `abs_mass_accuracy_deviation_mean` |
| RT deviation (2) | `rt_deviation`, `abs_rt_deviation` |
| MS1 (2) | `ms1_precursor_coelution`, `ms1_isotope_cosine` |
| Median polish (2) | `median_polish_cosine`, `median_polish_residual_ratio` |
| SG-weighted multi-scan (4) | `sg_weighted_xcorr`, `sg_weighted_cosine`, `median_polish_min_fragment_r2`, `median_polish_residual_correlation` |

### Compression

ZSTD compression via Apache Parquet `WriterProperties`. Typical compression ratio: 3–5×.

### Selective Loading

Osprey never loads the full Parquet unless needed. Four specialized loaders read only the columns required:

| Loader | Columns Read | Purpose |
|---|---|---|
| `load_fdr_stubs_from_parquet()` | 9 identity + boundary columns | Convert to lightweight `FdrEntry` stubs for FDR |
| `load_pin_features_from_parquet()` | 21 feature columns + `entry_id` | SVM re-scoring with trained model |
| `load_cwt_candidates_from_parquet()` | `entry_id`, `cwt_candidates` | Reconciliation planning (peak candidates) |
| `load_scores_parquet()` | All ~40 columns | Full entry rehydration for blib output |

### Inspection

```bash
# View schema (requires pyarrow)
python -c "import pyarrow.parquet as pq; print(pq.read_schema('sample.scores.parquet'))"

# Read into pandas
python -c "import pandas as pd; df = pd.read_parquet('sample.scores.parquet'); print(df.shape); print(df.columns.tolist())"
```

---

## 4. BiblioSpec Output (`*.blib`)

The final output file in SQLite BiblioSpec format. See [08 - BiblioSpec Output Schema](08-blib-output-schema.md) for complete documentation of the schema, Osprey extension tables, and Skyline integration.

---

## Cleanup

All intermediate files can be safely deleted. They will be recreated on the next run.

```bash
# Remove all intermediate files for a sample
rm -f sample.spectra.bin sample.scores.parquet

# Remove calibration (forces recalibration on next run)
rm -f sample.calibration.json

# Remove all intermediates in a directory
rm -f *.spectra.bin *.scores.parquet
```

The `.calibration.json` files are small and worth keeping if you plan to re-run with the same LC-MS setup, since they skip the calibration phase entirely.
