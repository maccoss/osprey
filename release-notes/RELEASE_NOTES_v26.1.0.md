# Osprey v26.1.0 Release Notes

First public release of Osprey under the `YY.feature.patch` versioning scheme, a peptide-centric DIA analysis tool for mass spectrometry proteomics. Osprey uses fragment XIC co-elution analysis to detect peptides in DIA data, with machine learning scoring and rigorous FDR control. Results integrate with Skyline for quantification.

This release consolidates all development since the initial `v0.1.0` tag, including a series of bug fixes for large-experiment memory usage and FDR scoring correctness.

## Core Analysis Pipeline

- Fragment XIC co-elution search with 21 PIN features (pairwise correlation, peak shape, spectral matching, mass accuracy, median polish, Savitzky-Golay weighted multi-scan scoring)
- Per-file RT calibration via LOESS regression with automatic retry and expanded sampling for large libraries
- Per-file MS1/MS2 mass calibration with ppm error correction
- Enzyme-aware decoy generation (sequence reversal preserving terminal residues)
- Support for DIA-NN TSV, EncyclopeDIA elib, and BiblioSpec blib spectral libraries
- YAML configuration files with CLI argument overrides
- BiblioSpec blib output with Osprey extension tables for Skyline integration

## FDR Control

- Native Rust Percolator implementation (linear SVM, cross-validated, semi-supervised)
- Two-level FDR: run-level per file, experiment-level across all files
- Dual precursor + peptide level FDR via max(precursor_qvalue, peptide_qvalue)
- Two-pass scoring: first-pass on all entries, second-pass restricted to first-pass passing precursors for concentrated target-decoy competition
- Posterior error probability (PEP) estimation via target-decoy score modeling
- External Mokapot integration with parallel workers and progress streaming
- Configurable FDR level filtering (precursor, peptide, or both) via `--fdr-level`

## Protein-Level Analysis

- Native Rust protein parsimony: bipartite graph construction, identical-set grouping, subset elimination, greedy set cover
- Razor peptide assignment for shared peptides (greedy assignment to protein with most unique peptides)
- Picked-protein FDR with target-decoy competition using 1-PEP (probability correct) scoring
- Configurable shared peptide handling (`--shared-peptides all|razor|unique`)
- Optional protein-level FDR control via `--protein-fdr`
- CSV reports for peptides and proteins

## Inter-Replicate Peak Reconciliation

- Consensus RT computation using weighted median of FDR-passing detections mapped to library RT space, weighted by coelution sum
- Peak width estimation using weighted median (weighted by coelution sum) to prevent noisy detections from inflating forced integration boundaries
- Only detections passing run-level FDR contribute to consensus (excludes noise peaks from replicates where the peptide was not confidently detected)
- Refined per-file RT calibration using consensus peptides
- Three reconciliation actions: Keep (peak already correct), UseCwtPeak (alternate CWT candidate overlaps consensus), ForcedIntegration (no peak overlaps, impute boundaries from consensus RT and median peak width)
- Gap-fill scoring for precursors passing FDR in some replicates but absent from others
- Multi-charge consensus: aligns charge states of the same peptide to the best-scoring charge state's peak

## Performance and Scalability

### Disk-Backed Caching (scales to 1000+ files)

- **Per-file Parquet caches** (`.scores.parquet`): ZSTD-compressed scored entries with 21 PIN features, fragments, CWT candidates. Eliminates redundant re-scoring on reruns.
- **Parquet metadata hashing**: SHA-256 hashes of search parameters, library identity (path + size + mtime), and Osprey version stored in Parquet key-value metadata. Stale caches from parameter or library changes are automatically detected and deleted.
- **Reconciliation metadata**: Parquet files track whether reconciliation has been applied and with which parameters. On rerun with matching parameters, reconciliation is skipped entirely.
- **FDR score sidecars** (`.1st-pass.fdr_scores.bin`, `.2nd-pass.fdr_scores.bin`): Per-file binary files storing SVM discriminant scores. Enables skipping Percolator SVM training on reruns (~9 MB per file).
- **Skip-Percolator optimization**: When all files have valid reconciled Parquet and loaded SVM scores from sidecars, Osprey skips SVM training entirely and recomputes q-values from cached scores using the same two-pass FDR structure.
- **Binary spectra cache** (`.spectra.bin`): Serialized MS2 spectra for fast reload during reconciliation re-scoring.
- **Calibration JSON cache**: Per-file RT and mass calibration parameters.
- **Binary library cache** (`.libcache`): Serialized spectral library for fast reload.

### Memory Architecture (54 GB RAM for 240-file experiments)

- **FdrEntry stubs** (~128 bytes per entry): Lightweight in-memory representation replacing full scored entries (~940 bytes). Heavy data (features, fragments, CWT candidates) stays on disk in Parquet caches and is loaded on demand.
- **Arc\<str\> string interning**: Deduplicates modified_sequence strings across files. With 240 GPF files and ~3.5M unique sequences, interning reduces heap from ~6 GB to ~90 MB.
- **BlibPlanEntry streaming** (~96 bytes per entry): Blib output uses a 5-column Parquet projection instead of loading full entries. For 24M passing entries: ~2.3 GB vs ~22 GB.
- **Sequential reconciliation re-scoring**: Files processed one at a time (not parallel) to prevent OOM from simultaneous spectra loading (~1.5 GB per file).
- **Batched CWT loading**: Dynamic batch size based on available RAM and CPU cores.
- **Reconciliation overlay persistence**: Re-scored entries written back to Parquet immediately per file, then dropped, avoiding accumulation across files.
- **LightFdr extraction**: Heavy FDR stubs dropped before blib output, replaced with lightweight structs containing only q-values needed for output filtering.

### Safe File I/O for Network Storage

- All cache writers (Parquet, spectra.bin, calibration JSON, score sidecars, blib) write to local temp directory first, then `copy_and_verify` to final destination. Prevents corrupt/zero-byte files from process kills or network interruptions on NAS/CIFS mounts.

## Bug Fixes (since v0.1.0)

### FDR Scoring

- **Fixed missing `parquet_index` population on fresh stub creation**: FdrEntry stubs created from freshly scored entries were left with `parquet_index = 0`, causing Percolator's feature lookup to return `file_features[0]` for every entry. The SVM trained on identical features and produced 0 passing precursors at 1% FDR on large datasets. Now sets `parquet_index = i` during the Vec enumeration.
- **Fixed compaction-vs-parquet-index bugs in Percolator**: The streaming and direct Percolator paths used Vec position (`local_idx`) as the Parquet row index after FDR stub compaction. After compaction the Vec is shorter than the Parquet, causing out-of-range panics (direct path) or loading the wrong features (streaming path). All three feature-loading sites now use `entry.parquet_index` for Parquet lookups.
- **Fixed gap-fill stub `parquet_index` sentinel**: Gap-fill entries (added during reconciliation) now use `parquet_index = u32::MAX` as a sentinel to indicate their features come from the overlay, not the Parquet file.
- **Fixed two-pass FDR in skip-Percolator path**: When cached SVM scores were loaded from sidecars, the q-value recomputation was running only the first-pass TDC instead of the two-pass structure (first unrestricted, then restricted to first-pass passing precursors). Experiment-level results were dropping by ~40%. Now correctly runs two-pass q-value computation in the skip path.
- **Fixed second-pass sidecar write location**: Moved the second-pass FDR score sidecar persistence out of the reconciliation scoped block so that temporary allocations are freed first. Prevented OOM on 240-file experiments.

### Reconciliation and Consensus

- **Fixed protein FDR scoring to use second-pass SVM scores**: Protein-level FDR scoring was using peptide scores from before stub compaction (first-pass scores). Reconciliation corrects peak boundaries for both targets and decoys; first-pass scores don't reflect those corrections. Now collects scores from compacted stubs after the second-pass FDR.
- **Skip multi-charge consensus when files are already reconciled**: Multi-charge consensus re-scoring was running unconditionally, even for files marked as already reconciled. It now skips the consensus step together with inter-replicate reconciliation when all files have matching reconciliation metadata.
- **Fixed consensus RT filtering**: Peak width estimation for forced integration boundaries now uses a weighted median across FDR-passing detections (weighted by coelution sum), preventing noisy detections with wide fallback peaks from inflating the forced integration window.

### Memory and I/O

- **Fixed blib output OOM on 240-file experiments**: The blib writer was loading all ~28M passing `CoelutionScoredEntry` objects from Parquet into a single Vec (~27 GB). Replaced with `BlibPlanEntry` (~96 bytes) loaded from a 5-column Parquet projection with fragments and metadata looked up from the in-memory library. Memory dropped from ~27 GB to ~2.3 GB for large experiments.
- **Fixed FDR stub memory for large experiments**: Added FDR stub compaction after first-pass FDR that drops entries for precursors that didn't pass in any replicate (saves ~21 GB for 240-file experiments). The `parquet_index` field preserves original Parquet row references for downstream lookups.
- **Fixed calibration cache invalidation**: Added `search_hash` to `CalibrationMetadata` so that stale calibration JSON files (e.g., from a different resolution mode) are detected and re-computed on rerun.

### Build and Install

- **Fixed Makefile install on WSL**: `cargo install` fails on WSL cross-filesystem builds due to a `libz-sys` cmake/zlib-ng extraction issue. The Makefile `install` target now copies the already-built release binary directly to `~/.cargo/bin/`.
- **Added `make install-clean` target**: Forces a full `cargo clean` before rebuild for cases when incremental compilation artifacts are stale.

### Regression Tests

Added tests to catch these bug classes in the future:

- `test_compaction_preserves_parquet_index`: verifies compaction keeps both targets and paired decoys with correct parquet_index
- `test_feature_lookup_uses_parquet_index_after_compaction`: documents the correct lookup pattern
- `test_gap_fill_sentinel_parquet_index`: verifies gap-fill entries use u32::MAX sentinel
- `test_to_fdr_entry_default_parquet_index_is_zero`: documents the dangerous default
- `test_fresh_stubs_must_have_correct_parquet_index`: verifies fresh stubs populate parquet_index from Vec position
- `test_buggy_pattern_collapses_parquet_index_to_zero`: explicitly catches the buggy pattern
- `test_consensus_peak_width_excludes_non_passing_detections`: verifies peak width estimation uses only FDR-passing entries

## Known Issues

- **Protein-level FDR target/decoy separation**: With the current SVM-based scoring, target and decoy protein groups separate too cleanly, producing artificially low protein FDR values. The picked-protein competition rarely selects decoy winners. This is tracked for a follow-up release; the current CSV output includes all passing target proteins with their best peptide scores, but the protein q-values should be treated with caution until the scoring metric is refined.

## Build and CI

- Makefile with `make check`, `make test`, `make install`, `make release` targets
- CI enforces `cargo fmt --check` and `cargo clippy -D warnings`
- GitHub Actions release workflow for Linux x86_64

## Compatibility

- Rust 1.75+ required
- External dependency: Mokapot (Python, `pip install mokapot`) for Mokapot FDR method only
- Input: mzML (DIA), spectral libraries (DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib)
- Output: BiblioSpec blib (for Skyline), CSV reports (peptide, protein), optional Parquet report
