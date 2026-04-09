# Osprey v26.1.0 Release Notes

First public release of Osprey, a peptide-centric DIA analysis tool for mass spectrometry proteomics. Osprey uses fragment XIC co-elution analysis to detect peptides in DIA data, with machine learning scoring and rigorous FDR control. Results integrate with Skyline for quantification.

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

## Build and CI

- Makefile with `make check`, `make test`, `make install`, `make release` targets
- CI enforces `cargo fmt --check` and `cargo clippy -D warnings`
- GitHub Actions release workflow for Linux x86_64

## Compatibility

- Rust 1.75+ required
- External dependency: Mokapot (Python, `pip install mokapot`) for Mokapot FDR method only
- Input: mzML (DIA), spectral libraries (DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib)
- Output: BiblioSpec blib (for Skyline), CSV reports (peptide, protein), optional Parquet report
