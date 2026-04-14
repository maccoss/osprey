# CLAUDE.md - Development Context for Osprey

This file provides context for Claude Code when working on the Osprey project.

## Project Overview

Osprey is a Rust-based peptide-centric DIA (Data-Independent Acquisition) analysis tool for mass spectrometry proteomics. It uses fragment XIC co-elution analysis to detect peptides in DIA data, with machine learning scoring and rigorous FDR control. Results integrate with Skyline for quantification.

**Current Status**: Working prototype with 76,169 precursors at 1% FDR from single Astral file.

## Architecture

### Workspace Structure

```
osprey/
├── Cargo.toml                    # Workspace root
├── crates/
│   ├── osprey-core/              # Core types, config, errors, traits
│   ├── osprey-io/                # File I/O (mzML, blib, DIA-NN TSV)
│   ├── osprey-chromatography/    # Peak detection, RT/mass calibration
│   ├── osprey-scoring/           # Feature extraction, decoy generation
│   ├── osprey-fdr/               # FDR control, Mokapot integration
│   ├── osprey-ml/                # Machine learning (SVM, LDA, PEP, q-values)
│   └── osprey/                   # Main library + CLI binary
├── scripts/
│   ├── evaluate_calibration.py   # Calibration report generator
│   └── inspect_mokapot_weights.py # Feature weight analysis
└── docs/                         # Algorithm documentation
```

### Key Components

- **osprey-core**: `LibraryEntry`, `Spectrum`, `OspreyConfig`, `IsolationWindow`, `CoelutionFeatureSet`, `FdrEntry`
- **osprey-io**: `MzmlReader` (uses mzdata crate), `DiannTsvLoader`, `ElibLoader`, `BlibLoader`, `BlibWriter`
- **osprey-chromatography**: `PeakDetector`, `RTCalibrator`, `MzCalibration`, `CalibrationParams`
- **osprey-scoring**: `SpectralScorer`, `DecoyGenerator`, batch scoring
- **osprey-fdr**: `FdrController`, `MokapotRunner` (PIN file + semi-supervised FDR), `protein` module (parsimony + picked-protein FDR)
- **osprey-ml**: `LinearSvmClassifier`, `PepEstimator`, `LDA`, matrix operations
- **osprey** (main): `pipeline::run_analysis`, `reconciliation`, `logging::TwoTierLogger`

## Build Commands

```bash
# Development build
cargo build

# Release build (for production runs)
cargo build --release

# Run tests
cargo test

# Run with arguments
cargo run --release -- -i sample.mzML -l library.tsv -o results.blib

# Install globally
cargo install --path crates/osprey
```

## CI Requirements

**IMPORTANT**: Before committing any Rust code, always run these checks locally:

```bash
cargo fmt          # Format code — CI rejects formatting diffs
cargo clippy --all-targets --all-features -- -D warnings  # Lint — CI treats warnings as errors
cargo test         # Run tests
```

CI runs `cargo fmt --check` and `cargo clippy -D warnings` (including test targets). Any formatting difference or clippy warning will fail the build. Always run `cargo fmt` after writing or editing code.

## Key Files

- `docs/README.md` - Pipeline overview and algorithm documentation
- `docs/07-fdr-control.md` - Two-level FDR with Mokapot
- `docs/08-blib-output-schema.md` - BiblioSpec output format and Skyline integration
- `docs/12-intermediate-files.md` - Intermediate file formats (calibration JSON, spectra cache, Parquet scores, FDR score sidecars, cache invalidation, memory architecture)
- `crates/osprey/src/pipeline.rs` - Main analysis pipeline
- `crates/osprey/src/main.rs` - CLI entry point
- `crates/osprey-fdr/src/mokapot.rs` - Mokapot integration
- `crates/osprey-io/src/output/blib.rs` - BiblioSpec blib writer
- `crates/osprey-core/src/types.rs` - CoelutionFeatureSet (~47 fields, 21 used in PIN), FdrEntry (~80 bytes inline, 13 fields)
- `crates/osprey-chromatography/src/calibration/` - RT and mass calibration

## Configuration

Osprey supports YAML configuration files:

```bash
# Generate template
osprey --generate-config config.yaml

# Run with config
osprey --config config.yaml
```

CLI arguments override config file values.

## Dependencies

- **mzdata**: mzML parsing (v0.63)
- **rusqlite**: SQLite for blib output (bundled)
- **ndarray**: Matrix operations
- **rayon**: Parallelism
- **clap**: CLI parsing
- **serde_yaml**: YAML config
- **chrono**: Timestamps

External:
- **mokapot**: Semi-supervised FDR (Python, `pip install mokapot`)

## Current Status

### Implemented (Working)

- mzML parsing with DIA isolation windows
- DIA-NN TSV library loading
- EncyclopeDIA elib library loading
- BiblioSpec blib library loading
- Enzyme-aware decoy generation (reversal)
- RT calibration with LOESS regression
- MS1/MS2 mass calibration
- Fragment XIC co-elution search with pairwise correlation scoring (21 PIN features)
- Peak detection in XIC time series
- Tukey median polish for robust peak boundaries and fragment scoring
- Two-level FDR control (run + experiment level)
- Mokapot integration with parallel workers
- Progress streaming from Mokapot
- Feature weight inspection (`--save_models`)
- blib output with Osprey extension tables (library theoretical fragments)
- YAML configuration
- CLI with all core options
- Calibration JSON save/load
- Calibration HTML report generation
- FdrEntry memory optimization with per-file Parquet caching

### Feature Set (21 PIN Features)

The `CoelutionFeatureSet` struct has ~47 fields total, but only 21 are written to the PIN file for FDR scoring. The remaining fields are computed but not used (removed during feature weight optimization).

Pairwise coelution (3): fragment_coelution_sum, fragment_coelution_max, n_coeluting_fragments
Peak shape (3): peak_apex, peak_area, peak_sharpness
Spectral at apex (3): xcorr, consecutive_ions, explained_intensity
Mass accuracy (2): mass_accuracy_deviation_mean, abs_mass_accuracy_deviation_mean
RT deviation (2): rt_deviation, abs_rt_deviation
MS1 features (2): ms1_precursor_coelution, ms1_isotope_cosine
Median polish (2): median_polish_cosine, median_polish_residual_ratio
SG-weighted multi-scan (4): sg_weighted_xcorr, sg_weighted_cosine, median_polish_min_fragment_r2, median_polish_residual_correlation

26 features removed from PIN (still in CoelutionFeatureSet struct but not used for scoring):
fragment_corr_0..5, fragment_coelution_min, n_fragment_pairs, dot_product, dot_product_smz,
dot_product_top4..6, dot_product_smz_top4..6, signal_to_noise, peak_symmetry,
elution_weighted_cosine, peptide_length, missed_cleavages, peak_width, n_scans,
fragment_coverage, hyperscore, sequence_coverage, mass_accuracy_std, median_polish_rsquared

### TODO (Future)

- EMG peak fitting (Levenberg-Marquardt)
- Feature weight optimization (remove low-value features)
- Background correction
- Two-step search strategy
- Iterative candidate expansion

## Testing

Test data should be placed in a `example_test_data/` directory (not committed):
- `*.mzML` - DIA mass spec files
- `*.tsv` - DIA-NN format spectral libraries
- `*.elib` - EncyclopeDIA spectral libraries
- `*.blib` - BiblioSpec spectral libraries

## Critical Invariants

### Cross-Validation Grouping

**CRITICAL**: Any operation that splits data into cross-validation folds or subsamples training data MUST keep these groups together:

1. **Target-decoy pairs**: A target and its paired decoy (linked by `base_id = entry_id & 0x7FFFFFFF`) must ALWAYS be in the same fold. If pairs are split across folds, unpaired targets in the training set auto-win competition, inflating the positive training set and causing the SVM to become too permissive — this silently corrupts FDR estimates.

2. **Same peptide, different charge states**: All charge states of the same peptide (same `base_id` → same target sequence) must be in the same fold. Otherwise, information about the peptide leaks between training and test sets.

3. **Subsampling**: When subsampling entries for SVM training efficiency, subsample **paired target-decoy groups** (by `base_id`), not individual entries. Subsampling should happen BEFORE fold splitting, per the Percolator paper (The et al., 2016, PMC5059416). The learned model is still applied to ALL entries for scoring.

This applies to: `percolator.rs` (fold assignment, subsampling), `calibration_ml.rs` (LDA fold assignment), and any future code that partitions entries for cross-validation.

### Protein FDR Scoring (Two-Pass Picked-Protein, Savitski 2015)

Osprey implements **true picked-protein FDR** (Savitski et al. 2015, PMC4563723) in `crates/osprey-fdr/src/protein.rs::compute_protein_fdr`. Each protein is scored by its **single best peptide** (max SVM score), not by a sum — sum aggregation is length-biased (Savitski explicitly tested and rejected it). Target and decoy sides compete via **pairwise picking**: each group produces exactly one winner. Only target winners are exposed in `group_qvalues`; decoy winners are statistical machinery for cumulative FDR.

**Two-pass architecture** (in `pipeline.rs`):

1. **First-pass protein FDR** runs AFTER first-pass peptide FDR but BEFORE compaction. Uses the **full pre-compaction peptide pool** so targets and decoys are symmetric. Writes `run_protein_qvalue` on FdrEntry stubs. This is used as a GATE for protein-aware compaction (rescue borderline peptides from strong proteins) and reconciliation consensus selection.

2. **Compaction** retains a peptide if EITHER `run_peptide_qvalue <= reconciliation_compaction_fdr` (default 0.01) OR `run_protein_qvalue <= config.protein_fdr`. The protein rule is additive and only applies when `--protein-fdr` is set.

3. **Reconciliation** uses first-pass protein q-values as an optional rescue gate in `compute_consensus_rts()`. Peptides from strong proteins can anchor consensus even if their own peptide q-value is borderline.

4. **Second-pass protein FDR** runs AFTER second-pass peptide FDR on the compacted + reconciled + second-pass-scored pool. Writes `experiment_protein_qvalue` (AUTHORITATIVE). Feeds `FdrLevel::Protein` filtering and the protein CSV report.

**CRITICAL**: The first-pass MUST see the full pre-compaction pool. v26.1.2 had a calibration bug where protein FDR only ran after compaction, giving one-sided competition because decoys paired with non-passing targets were dropped. Keeping the first-pass picked-protein FDR before compaction preserves target/decoy symmetry.

**CRITICAL: Ranking score is raw SVM discriminant — and on current Osprey output all three candidates produce a collapsed decoy null.** During development we tried three ranking scores:

- **Peptide q-value**: 99% of decoys have `q = 1.0` because they lost peptide-level TDC, leaving too few with meaningful scores. Produced 2/6102 decoy winners on Stellar — under-calibrated.
- **Peptide PEP**: `PepEstimator` is fit on TDC winners only and clamps out-of-range scores to `bins[0] ≈ 1.0`. Losing decoys all clamped to ~1.0, same failure mode. Produced 2/6102 decoy winners on Stellar.
- **SVM discriminant** (current, kept): every entry has a well-defined score on the same scale. An earlier Stellar run showed 348/5988 decoy winners (5.8% decoy rate) and was used to justify this choice, BUT that run preceded the FDR q-value mapping fix. With the bug fixed, compaction keeps the correct pool, Percolator trains on the full set, and the resulting SVM separates targets and decoys so sharply that picked-protein has no tail overlap. Post-fix Stellar produces 2/6099 decoy winners — same failure mode as q-value and PEP. SVM is kept as the least-bad option (honest monotone score over all entries), but **picked-protein FDR on Osprey's current output provides essentially no control beyond peptide-level FDR**. Report protein q-values as "protein has a peptide passing peptide-level FDR", not as independently calibrated probabilities. See `docs/07-fdr-control.md` → "Why SVM Score (Not q-Value or PEP) — and What Actually Happens in Practice" and the "Known Limitation" subsection for the full story.

**Known root cause (deferred): decoy asymmetry.** Target library entries use AI-predicted spectra and RTs (Carafe). `DecoyGenerator::Reverse` produces reversed peptide sequences but keeps the target's RT and does a mechanical reversal of fragments — the decoy side is NOT re-predicted by the same model. Every feature the SVM learns (RT residual, spectral cosine, fragment co-elution, XCorr) is systematically different for decoys than for targets because the generative processes differ. Fixing this would require running the AI predictor on the reversed sequence so target and decoy libraries come from the same model. Until then the tail overlap picked-protein needs cannot be produced. **Do not revisit picked-protein scoring without either fixing decoy generation OR running on real data first and confirming the decoy-winner fraction is meaningful (e.g., > 1% of total winners).**

**Do NOT use protein-level PEP.** Protein-level PEP is intentionally not computed (`ProteinFdrResult` has no `group_pep` field, and `write_protein_report` does not emit a `protein_pep` column). Peptide-level and precursor-level PEP are unaffected and still written to the blib.

**Do NOT revert to single-pass or composite scoring.** v26.1.2's composite log-likelihood approach produced artificially low protein q-values because the aggregation inflated target scores relative to decoys.

See `docs/16-protein-parsimony.md` for the full algorithm, `docs/07-fdr-control.md` for the pipeline integration, and the `test_picked_protein_fdr_*` regression tests in `protein.rs`.

## Code Style

- Use `log::info!`, `log::debug!` for logging
- Error handling via `thiserror` and `OspreyError`
- Prefer iterators and functional patterns
- Document public APIs with `///` comments
- Use "precursor" (not "PSM") in user-facing messages
- Precursor = peptide + charge state

## Utility Scripts

### evaluate_calibration.py

Generate HTML report from calibration JSON:
```bash
python scripts/evaluate_calibration.py calibration.json --output report.html
```

Features:
- RT calibration curve with residuals
- MS1/MS2 mass error histograms
- Candidate density heatmap (optimized with binary search)

### inspect_mokapot_weights.py

View Mokapot feature weights to identify important features:
```bash
python scripts/inspect_mokapot_weights.py mokapot.model.pkl
```

## FDR Control Logic

Osprey enforces FDR at **two levels** (run and experiment) and **two scopes** (precursor and peptide).

### Two-Level FDR

1. **Run-level FDR**: Each file is scored independently. Target-decoy competition determines q-values per file. This controls the per-file false discovery rate.

2. **Experiment-level FDR**: Takes the **single best-scoring observation** per precursor (modified_sequence + charge) across ALL files in the experiment. Target-decoy competition on this deduplicated set determines experiment-level q-values. This controls the experiment-wide false discovery rate.

### Dual Precursor + Peptide FDR

Both Percolator and Mokapot compute q-values at precursor level (modified_sequence + charge) and peptide level (modified_sequence only). A precursor must pass FDR at **both** levels.

This is enforced via: `effective_qvalue = max(precursor_qvalue, peptide_qvalue)`

- `max(0.005, 0.015) = 0.015 > 0.01` → rejected (peptide fails)
- `max(0.015, 0.005) = 0.015 > 0.01` → rejected (precursor fails)
- `max(0.003, 0.008) = 0.008 < 0.01` → accepted (both pass)

The max is applied at both run and experiment levels:
- `run_qvalue = max(run_precursor_qvalue, run_peptide_qvalue)`
- `experiment_qvalue = max(experiment_precursor_qvalue, experiment_peptide_qvalue)`

For Mokapot: `mokapot.psms.txt` provides precursor-level q-values; `mokapot.peptides.txt` provides peptide-level q-values. Both are parsed and combined via max.

### Multi-File Observation Propagation

After experiment-level FDR determines which precursors pass, **all per-file target observations** for those precursors are included in the output (blib and report). This ensures:
- Each file gets its own RT boundaries for a passing precursor
- The best experiment_qvalue is propagated to all observations of that precursor
- Skyline can use per-file peak boundaries for quantification across replicates/GPF files

### Blib RetentionTimes: Nullable retentionTime for Skyline ID Lines

In the RetentionTimes table, `retentionTime` controls Skyline's ID line display:
- **Set to apex RT**: When the precursor passes run-level FDR in that file → Skyline shows an ID line
- **Set to NULL**: When the precursor did NOT pass run-level FDR in that file (but passed experiment-level FDR via another replicate) → Skyline uses `startTime`/`endTime` for quantification boundaries without showing an ID line

This distinction is important: a NULL `retentionTime` with populated `startTime`/`endTime` tells Skyline "integrate here for quantification, but this is not an independent identification." See `docs/08-blib-output-schema.md` for full details.

### FDR Pipeline Flow (Multi-File)

1. Score each file: write full entries to per-file Parquet cache, convert to FdrEntry stubs
2. **Best-per-precursor subsampling**: Find the best-scoring observation (by coelution_sum) per base_id across all files — one target + one decoy per precursor — to maximize peptide diversity for SVM training. Subsample from this deduplicated set if > max_train.
3. Train SVM on the subsampled best-per-precursor entries (train_only mode suppresses FDR logging from training subset)
4. Score ALL entries with trained model by streaming through per-file Parquet caches
5. Compute run-level q-values per file (precursor + peptide level)
6. Select best observation per precursor across experiment
7. Compute experiment-level q-values (precursor + peptide level)
8. Apply `max(precursor, peptide)` at both run and experiment levels
9. Determine passing precursors: `experiment_qvalue <= threshold`
10. Reload full entries from Parquet only for files with passing precursors
11. Include ALL per-file observations for passing precursors in output
12. Propagate best experiment_qvalue to all observations

### Memory Architecture

Osprey uses a tiered disk-backed architecture to scale to 1000+ files:

- **FdrEntry** (~128 bytes): Lightweight stub with fields: entry_id, parquet_index, is_decoy, charge, scan_number, apex/start/end_rt, coelution_sum, score, 6 q-values (run/experiment x precursor/peptide/protein), pep, modified_sequence (`Arc<str>`). After first-pass FDR, non-passing entries are compacted out (saves ~21 GB for 240-file experiments). The `parquet_index` preserves original Parquet row references after compaction.
- **LightFdr** (~48 bytes): Extracted from FdrEntry stubs before blib output. Contains only q-values, score, pep, charge, and modified_sequence. Dropping FdrEntry stubs frees ~19 GB for 240-file experiments.
- **BlibPlanEntry** (~96 bytes): Loaded from a 5-column Parquet projection for blib output. Contains entry_id (for library lookup), RT boundaries, q-values, and file_name_idx. Fragments, modifications, and protein mappings come from the in-memory library. Avoids loading full entries (~22 GB for 24M passing entries).
- **Parquet caches** (`.scores.parquet`): Per-file ZSTD-compressed caches storing 21 PIN features, fragments, CWT candidates. Metadata footer contains SHA-256 hashes of search parameters, library identity, and reconciliation state for automatic cache invalidation.
- **FDR score sidecars** (`.1st-pass.fdr_scores.bin`, `.2nd-pass.fdr_scores.bin`): Per-file binary files storing SVM discriminant scores (~9 MB each). Enable skipping Percolator SVM training on reruns.
- **Selective loading**: `load_blib_plan_entries()` for 5-column projection, `load_fdr_stubs_from_parquet()` for FDR stubs, `load_pin_features_from_parquet()` for SVM scoring, `load_cwt_candidates_from_parquet()` for reconciliation planning
- **Sequential re-scoring**: Reconciliation re-scoring processes files sequentially (not parallel) because each file loads ~3 GB (spectra + full Parquet entries); parallel loading causes OOM on large experiments
- **Steady-state RAM**: ~80 bytes × entries × files during FDR (vs ~940 bytes without caching), ~96 bytes × passing entries during blib output, plus shared `Arc<str>` interning pool

## Recent Changes

- Removed ridge regression mode (coelution is now the sole search mode)
- Fixed blib output to write library theoretical fragments (not observed DIA peaks)
- Fixed experiment-level FDR: entries not in experiment results no longer fall back to run q-values
- Removed pearson_correlation and spearman_correlation from coelution features (poor discrimination)
- Added blib round-trip tests (fragment, modseq, protein mapping, multi-run RT)
- Replaced library-assisted median polish with Tukey median polish for peak boundaries
- Added 3 median polish features: cosine, R², residual_ratio
- Fixed all cosine/correlation scoring to include ALL library fragments (0 for unmatched)
- Added fragment co-elution features (9 features: sum, min, n_positive, per-fragment corr)
- Added elution-weighted cosine, mass accuracy features, Percolator-style features
- Added MS1 precursor co-elution and isotope cosine features
- Fixed RT calibration NaN bug from duplicate library_rts
- Fixed LDA calibration scoring: non-negative weights, early stopping
- Added `--save_models` to Mokapot for feature weight inspection
- Fixed two-level FDR to aggregate by precursor (peptide+charge), not peptide
- Single replicate now skips experiment-level FDR (uses run-level directly)
- Added parallel workers to Mokapot (auto-detected, capped at 8)
- Added progress streaming from Mokapot to console
- Enforced dual precursor + peptide level FDR via max(precursor_qvalue, peptide_qvalue)
- Added mokapot.peptides.txt parsing for peptide-level FDR in Mokapot paths
- Added multi-file observation propagation: all per-file observations for passing precursors included in output
- Nullable retentionTime in blib RetentionTimes: NULL (no ID line) for runs not passing FDR, populated for runs passing FDR
- Added binary spectra cache (.spectra.bin) for faster second-pass mzML loading
- Merged multi-charge consensus + cross-run reconciliation into single post-FDR phase (one spectra load per file)
- Skip consensus re-scoring for peptide groups where no charge state passes FDR
- FdrEntry memory architecture: lightweight stubs (~80 bytes inline, `Arc<str>`-interned `modified_sequence`) replace full scored entries (~940 bytes) in memory after per-file Parquet caching; heavy data (features, fragments, CWT candidates) reloaded on-demand from disk, enabling 1000+ file experiments without OOM
- Selective Parquet loaders: `load_cwt_candidates_from_parquet()` for reconciliation, PIN feature loading for FDR, full entry reload for blib output
- Removed `shrink_for_fdr()` in favor of FdrEntry conversion + Parquet caching
- Best-per-precursor subsampling for streaming Percolator: finds best observation per base_id across all files before subsampling, maximizing peptide diversity for SVM training on 100+ file experiments
- Added `train_only` flag to PercolatorConfig: suppresses per-file/experiment FDR logging when training on a subset (where FDR numbers are meaningless)
- Reconciliation consensus uses experiment-level FDR only (was `min(run, experiment)`), reducing consensus set to truly confident identifications
- Sequential reconciliation re-scoring: `iter_mut()` instead of `par_iter_mut()` to prevent OOM from parallel spectra loading (~3 GB per file)
- Changed `consensus_fdr` default from 0.05 to 0.01
- Removed `file_name` from FdrEntry (redundant with outer `Vec<(String, Vec<FdrEntry>)>` key) — saves ~40 bytes heap per entry
- Interned `modified_sequence` with `Arc<str>` — deduplicates peptide strings across GPF replicates, saving ~25 bytes heap per entry on large experiments
- Refactored reconciliation re-scoring to reuse `run_search()` with `boundary_overrides` parameter — eliminates separate sequential `rescore_for_reconciliation()` code path, gaining parallel window processing and shared per-window XCorr preprocessing (~50x speedup for reconciliation)
- Per-window XCorr preprocessing optimization: `preprocess_spectrum_for_xcorr()` called once per window, stored as `Vec<Vec<f32>>`, reused across all candidates via O(n_frags) bin lookup instead of redundant per-entry O(n_peaks) preprocessing
- Multi-charge consensus moved to post-FDR: uses SVM scores (not coelution_sum) to select consensus leader among FDR-passing charge states
- Consensus and reconciliation targets merged into single `run_search()` call per file (reconciliation wins on conflict)
- Added 5 boundary_overrides tests with synthetic data helpers (make_lib_entry_with_fragments, make_gaussian_spectra, make_test_config)
- Removed `rescore_for_reconciliation()`, `rescore_entry_in_window()`, `FileRescoreContext` (replaced by `run_search()` with boundary_overrides)
- Added Parquet metadata hashing for cache invalidation: `osprey.version`, `osprey.search_hash`, `osprey.library_hash`, `osprey.reconciled`, `osprey.reconciliation_hash` — SHA-256 of search parameters, library identity, and reconciliation config stored in Parquet footer
- Added `CacheValidity` enum (ValidFirstPass, ValidReconciled, Stale) — stale caches auto-detected and re-scored on parameter/library/version change
- Added FDR score sidecar files (`.1st-pass.fdr_scores.bin`, `.2nd-pass.fdr_scores.bin`) — persists SVM discriminant scores per file, enabling skip-Percolator on reruns (~9 MB per file)
- Skip-Percolator optimization: when all files have ValidReconciled Parquet + loaded SVM scores from sidecars, Osprey skips SVM training entirely and just recomputes q-values from cached scores
- Streaming blib output via `BlibPlanEntry` (~96 bytes per entry): loads 5-column Parquet projection, looks up fragments/mods/proteins from library — avoids loading full entries (~940 bytes each, ~22 GB for 24M passing entries)
- Added `LightFdr` struct (~48 bytes): extracted from FDR stubs before blib output to free ~19 GB of FDR stubs from memory
- Removed `write_blib_output()`, `build_shared_boundaries()`, `write_scored_report()` — replaced by streaming plan-based equivalents
- Safe NAS file writes: all cache writers (Parquet, spectra.bin, calibration JSON, score sidecars) write to local temp then `copy_and_verify` to final destination
- Removed unused config settings: `max_candidates_per_spectrum`, `memory_limit_gb`, `custom_bin_width`, `initial_tolerance_fraction`, `use_percentile_tolerance`
- Removed expensive unused `elution_weighted_cosine` feature computation (set to 0.0, not in PIN)
- Split FdrEntry q-value fields: `run_qvalue`/`experiment_qvalue` replaced by 6 separate fields (run/experiment x precursor/peptide/protein) with `effective_*_qvalue(FdrLevel)` methods
- Added `FdrLevel` enum (Precursor, Peptide, Protein, Both) for configurable output filtering via `--fdr-level`. Default is `Peptide`.
- Added `SharedPeptideMode` enum (All, Razor, Unique) for protein-level shared peptide handling via `--shared-peptides`
- Added native Rust protein parsimony in `osprey-fdr/src/protein.rs`: bipartite graph, identical-set grouping, subset elimination, greedy razor assignment, picked-protein FDR with TDC q-values
- Added `--protein-fdr` CLI flag and `protein_fdr` config for optional protein-level FDR control
- Protein FDR integrated into pipeline: `build_protein_parsimony` + `compute_protein_fdr` + `propagate_protein_qvalues` called after precursor/peptide FDR when `--protein-fdr` is set. Peptide scores collected before stub compaction to include decoy competition winners. CSV protein report with gene names, PEP, and q-values.
- FDR stub compaction after first-pass: drops non-passing entries, `parquet_index` field preserves original Parquet row references (~21 GB freed for 240-file experiments)
- Streaming blib output via `BlibPlanEntry` (~96 bytes): 5-column Parquet projection + library lookup, avoids loading full entries (~22 GB)
- Calibration cache validation: `search_hash` in CalibrationMetadata detects stale calibration files when parameters change (e.g., resolution mode)
- Parquet cache validation: SHA-256 hashes of search/library/reconciliation parameters in Parquet footer metadata
- Consensus RT fix: only FDR-passing detections included, weighted median peak widths by coelution_sum
- Safe NAS file writes: all cache writers use `copy_and_verify` pattern (write to local temp, verify, move to final destination)
- Fixed blib writer row misalignment after compaction: `LightFdr` now carries `parquet_index` for correct Parquet row lookup instead of zipping by Vec position
- Fixed reconciliation write-back using Vec position instead of `parquet_index` after compaction, causing re-scored entries to overwrite wrong Parquet rows
- Fixed reconciliation `determine_reconcile_action` to use apex proximity instead of boundary containment. Uses refined per-file calibration's `local_tolerance` (derived from thousands of consensus peptides) to determine whether a peak's apex is at the expected RT. CWT candidate selection also uses apex proximity, picking the closest candidate.
- Fixed reconciliation self-fulfilling tolerance: replaced `local_tolerance` interpolation with global per-file MAD-based tolerance (`3.0 * MAD * 1.4826`, floor 0.1 min). Outlier peptides no longer inflate the tolerance at their own RT.
- Protein parsimony now always runs (not gated by `--protein-fdr`). Only picked-protein q-value computation is optional. Peptide-to-protein-group mapping is always available for downstream filtering and interpretation. See `docs/16-protein-parsimony.md`.
- `SharedPeptideMode::Razor` now uses an iterative greedy set cover that picks the globally-best protein group at each step, processes its claimed shared peptides in sorted order, and repeats. Deterministic (path-independent) regardless of HashMap iteration order.
- `FdrLevel` now has a `Protein` variant (in addition to `Precursor`, `Peptide`, `Both`) and defaults to `Peptide`. When `--fdr-level protein` is selected, `passing_precursors` and `best_exp_q` read `experiment_protein_qvalue`, filtering the blib output by protein-level FDR. Requires `--protein-fdr` to be enabled; warns and falls back to peptide-level if not.
- Two-tier logging: clean `[M:SS]` terminal output + verbose log file. `--verbose` restores full detail on terminal. `TwoTierLogger` in `logging.rs` replaces `TeeWriter` + `env_logger`.
- Labeled progress bars: `run_search` accepts `search_label` parameter displayed in progress bar (Scoring, Re-scoring, Gap-fill CWT, Gap-fill forced)

## Versioning and Release Notes

Osprey uses `YY.feature.patch` versioning (e.g., `26.1.0` = 2026, first feature release, no patches).

- **Release notes** are maintained in `release-notes/RELEASE_NOTES_v{version}.md`
- See `release-notes/README.md` for the full format, conventions, and release process
- During development, maintain a draft release notes file for the next version and append entries as features and fixes land
- The workspace version in `Cargo.toml` is updated only at release time, not during development
- When adding significant features or fixes, add a brief entry to the current draft release notes file
