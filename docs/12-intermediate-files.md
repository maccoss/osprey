# 12 - Intermediate File Formats

Osprey produces several intermediate files during analysis. These files enable fast re-loading, disk-backed memory management, calibration reuse, and intelligent cache invalidation across runs. All intermediate files are written alongside the input mzML files.

## File Overview

| File Pattern | Format | Purpose | Typical Size | Reusable? |
|---|---|---|---|---|
| `*.calibration.json` | JSON | RT and mass calibration parameters | 5–50 KB | Yes (across runs with same LC-MS setup) |
| `*.spectra.bin` | Custom binary | Decoded spectra for fast reload | 2–5 GB (larger than mzML) | Auto-created on first run |
| `*.scores.parquet` | Apache Parquet (ZSTD) | Full scored entries with features, fragments, CWT candidates | 100 MB–2 GB | Validated via parameter hash |
| `*.1st-pass.fdr_scores.bin` | Raw binary (f64 LE) | SVM discriminant scores after first-pass Percolator | ~9 MB per file | Enables skipping Percolator on rerun |
| `*.2nd-pass.fdr_scores.bin` | Raw binary (f64 LE) | SVM discriminant scores after second-pass Percolator | ~9 MB per file | Enables skipping Percolator on rerun |
| `*.blib` | SQLite (BiblioSpec) | Final output for Skyline | 50–500 MB | N/A (output) |

---

## 1. Calibration JSON (`*.calibration.json`)

**Source**: `crates/osprey-chromatography/src/calibration/io.rs`

Stores the calibration parameters computed during Phase 2 (auto-calibration). If a calibration file already exists for an input file and its `search_hash` matches the current search parameters, Osprey skips the calibration phase and loads it directly. If the hash is missing (old file) or doesn't match (parameters changed, e.g., resolution mode), the stale file is deleted and calibration is re-run.

### Structure

```json
{
  "metadata": {
    "num_confident_peptides": 4523,
    "num_sampled_precursors": 10000,
    "calibration_successful": true,
    "timestamp": "2025-01-15T10:30:00Z",
    "search_hash": "de620a21cb001009...",
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

```text
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
| `load_blib_plan_entries()` | 5 columns (entry_id, apex_rt, start_rt, end_rt, bounds_area) | Lightweight plan for streaming blib output |
| `load_scores_parquet()` | All ~40 columns | Full entry rehydration (Parquet report only) |

### Cache Invalidation via Parquet Metadata

Each `.scores.parquet` file stores key-value metadata in its footer that enables intelligent cache reuse on subsequent runs. Osprey checks this metadata before accepting a cached file.

**Metadata keys stored in every Parquet file:**

| Key | Value | Purpose |
|---|---|---|
| `osprey.version` | Osprey binary version (e.g., `0.1.0`) | Version upgrade invalidates cache (feature set may change) |
| `osprey.search_hash` | SHA-256 hex string | Hash of all search parameters that affect scoring |
| `osprey.library_hash` | SHA-256 hex string | Hash of library file path + size + modification time |
| `osprey.reconciled` | `"true"` or `"false"` | Whether this file has been through reconciliation |
| `osprey.reconciliation_hash` | SHA-256 hex string (when reconciled) | Hash of reconciliation parameters + file set |

**Search parameter hash includes**: resolution mode, fragment/precursor tolerance (value + unit), prefilter setting, decoy method, all RT calibration parameters (bandwidth, tolerance factor, min/max tolerance, sample size, retry factor), and `reconciliation.top_n_peaks`.

**Library identity hash uses**: file path + file size + filesystem modification time (fast metadata check, no content hashing).

**Reconciliation parameter hash includes**: the search hash (if search changed, reconciliation is also invalid), `reconciliation.enabled`, `reconciliation.consensus_fdr`, `run_fdr`, and sorted input file stems (changing the file set changes consensus RTs).

**Cache validity states:**

| State | Condition | Action |
|---|---|---|
| `ValidReconciled` | All hashes match, reconciled=true | Skip scoring, FDR, and reconciliation |
| `ValidFirstPass` | Search + library hashes match, not reconciled | Skip scoring; run FDR and reconciliation |
| `Stale` | Any hash mismatch or missing metadata | Delete cache and re-score from scratch |

Pre-existing Parquet files without metadata (from older Osprey versions) are treated as `Stale` and automatically re-scored.

### Inspection

```bash
# View schema (requires pyarrow)
python -c "import pyarrow.parquet as pq; print(pq.read_schema('sample.scores.parquet'))"

# Read into pandas
python -c "import pandas as pd; df = pd.read_parquet('sample.scores.parquet'); print(df.shape); print(df.columns.tolist())"

# View cache metadata
python -c "import pyarrow.parquet as pq; print(pq.read_metadata('sample.scores.parquet').metadata)"
```

---

## 4. FDR Score Sidecars (`*.1st-pass.fdr_scores.bin`, `*.2nd-pass.fdr_scores.bin`)

**Source**: `crates/osprey/src/pipeline.rs` (`write_fdr_scores_sidecar()`, `load_fdr_scores_sidecar()`)

Lightweight binary files that persist the SVM discriminant scores from Percolator training. These files enable skipping the expensive Percolator SVM training and scoring step on subsequent runs with the same parameters.

### Why Two Passes?

Osprey runs Percolator twice in a multi-file experiment:

1. **First-pass FDR** (after initial scoring): Trains SVM on all entries, produces run-level q-values used to identify passing precursors for reconciliation.
2. **Second-pass FDR** (after reconciliation): Retrains SVM restricted to first-pass passing precursors + paired decoys, incorporating re-scored reconciliation entries for final q-values.

Each pass produces different SVM scores, so both are persisted independently.

### Format

Raw little-endian `f64` array, one value per entry, in Parquet row order:

```text
[f64 score_0][f64 score_1][f64 score_2]...[f64 score_N-1]
```

- **Entry count**: Must exactly match the number of entries in the corresponding `.scores.parquet`
- **Byte size**: `N_entries × 8 bytes` (e.g., 1.1M entries = ~9 MB)
- **Validation**: On load, the file size is checked against `entries.len() * 8`. Mismatches are silently ignored (scores will be recomputed).

### When Created/Used

- **Written**: After each Percolator FDR pass completes (first-pass after Phase 4, second-pass after Phase 5 reconciliation)
- **Loaded**: On subsequent runs, during stub loading from cached Parquet. The second-pass sidecar is preferred; the first-pass sidecar is used as fallback.
- **Safe to delete**: Yes — Percolator will retrain on next run

### Skip-Percolator Optimization

When all of the following conditions are met, Osprey skips both Percolator training passes entirely:

1. All files loaded stubs from cached Parquet (no fresh scoring)
2. All files have SVM scores loaded from sidecars
3. All files are `ValidReconciled` (or reconciliation is disabled / single file)

In this fast path, Osprey runs only `compute_fdr_from_stubs()` — target-decoy competition and PEP estimation on the existing scores — which takes seconds instead of minutes.

---

## 5. BiblioSpec Output (`*.blib`)

The final output file in SQLite BiblioSpec format. See [08 - BiblioSpec Output Schema](08-blib-output-schema.md) for complete documentation of the schema, Osprey extension tables, and Skyline integration.

---

## Memory Architecture

Osprey uses a tiered memory strategy to process experiments with hundreds of files without running out of memory. Each tier uses progressively lighter data structures, loading heavy data only when needed.

### Tier 1: Full Scored Entries (~940 bytes each)

`CoelutionScoredEntry` contains all scored data: 21 PIN features, fragment m/z and intensities, reference XIC, CWT candidates, peak boundaries, and metadata. These exist in memory only during per-file scoring (Phase 3) and are immediately written to Parquet and replaced with stubs.

### Tier 2: FDR Stubs (~128 bytes each)

`FdrEntry` retains the fields needed for FDR control and reconciliation: entry_id, parquet_index, is_decoy, charge, scan_number, apex/start/end RT, coelution_sum, score, 6 q-value fields (run/experiment x precursor/peptide/protein), pep, and an `Arc<str>`-interned modified_sequence.

After first-pass FDR, stubs for non-passing precursors are dropped (compaction). The `parquet_index` field preserves the original Parquet row reference for CWT/feature lookup after compaction.

For a 240-file experiment: `271M entries x 128 bytes = ~35 GB` before compaction, reduced to `~106M x 128 bytes = ~14 GB` after compaction (drops ~21 GB of non-passing entries).

### Tier 3: Lightweight FDR Data (~48 bytes each)

`LightFdr` extracts only q-values, score, pep, charge, and modified_sequence from the FDR stubs. Created right before blib output to free the heavier FDR stubs (~19 GB) from memory.

### Tier 4: Blib Plan Entries (~96 bytes each)

`BlibPlanEntry` contains only the fields needed for blib output: entry_id (for library lookup), charge, file_name_idx, q-values, RT boundaries, bounds_area, and interned modified_sequence. These are loaded from a 5-column Parquet projection merged with LightFdr data.

For 24M passing entries: `24M × 96 bytes ≈ 2.3 GB` (vs 22 GB if full entries were loaded).

Fragment m/z, intensities, modifications, and protein mappings are looked up from the in-memory library (via entry_id) rather than reloaded from Parquet. This eliminates the need to load any full `CoelutionScoredEntry` for blib output.

### Memory Flow (240-file experiment)

| Phase | In Memory | Peak RAM |
|---|---|---|
| Scoring (per file) | 1 file's spectra + full entries | ~3-5 GB |
| First-pass FDR (all files) | FDR stubs for all files | ~35 GB |
| Post-compaction | FDR stubs (passing only) | ~14 GB |
| Reconciliation | Compacted stubs + 1 file's spectra | ~17 GB |
| Protein FDR | Compacted stubs + peptide score map | ~14 GB |
| Blib output | LightFdr + BlibPlanEntry + library | ~5 GB |

---

## Cleanup

All intermediate files can be safely deleted. They will be recreated on the next run.

```bash
# Remove all intermediate files for a sample
rm -f sample.spectra.bin sample.scores.parquet
rm -f sample.1st-pass.fdr_scores.bin sample.2nd-pass.fdr_scores.bin

# Remove calibration (forces recalibration on next run)
rm -f sample.calibration.json

# Remove all intermediates in a directory
rm -f *.spectra.bin *.scores.parquet *.fdr_scores.bin

# Remove everything (forces full re-analysis)
rm -f *.spectra.bin *.scores.parquet *.fdr_scores.bin *.calibration.json
```

The `.calibration.json` files are small and worth keeping if you plan to re-run with the same LC-MS setup, since they skip the calibration phase entirely. The `.scores.parquet` and `.fdr_scores.bin` files are worth keeping if re-running with the same parameters, as they enable skipping scoring, FDR training, and reconciliation.
