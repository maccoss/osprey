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
│   └── osprey/                   # Main library + CLI binary
├── scripts/
│   ├── evaluate_calibration.py   # Calibration report generator
│   └── inspect_mokapot_weights.py # Feature weight analysis
└── docs/                         # Algorithm documentation
```

### Key Components

- **osprey-core**: `LibraryEntry`, `Spectrum`, `OspreyConfig`, `IsolationWindow`, `CoelutionFeatureSet`
- **osprey-io**: `MzmlReader` (uses mzdata crate), `DiannTsvLoader`, `ElibLoader`, `BlibLoader`, `BlibWriter`
- **osprey-chromatography**: `PeakDetector`, `RTCalibrator`, `MzCalibration`, `CalibrationParams`
- **osprey-scoring**: `SpectralScorer`, `DecoyGenerator`, batch scoring
- **osprey-fdr**: `FdrController`, `MokapotRunner` (PIN file + semi-supervised FDR)

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
- `crates/osprey/src/pipeline.rs` - Main analysis pipeline
- `crates/osprey/src/main.rs` - CLI entry point
- `crates/osprey-fdr/src/mokapot.rs` - Mokapot integration
- `crates/osprey-io/src/output/blib.rs` - BiblioSpec blib writer
- `crates/osprey-core/src/types.rs` - CoelutionFeatureSet (45 features)
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
- Fragment XIC co-elution search with pairwise correlation scoring (45 features)
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

### Feature Set (45 Features for PIN)

Pairwise coelution (11):
- fragment_coelution_sum, fragment_coelution_min, fragment_coelution_max
- n_coeluting_fragments, n_fragment_pairs
- fragment_corr_0..5 (per-fragment avg correlation with other fragments)

Peak shape (7):
- peak_apex, peak_area, peak_width, peak_symmetry
- signal_to_noise, n_scans, peak_sharpness

Spectral at apex (15):
- hyperscore, xcorr, dot_product, dot_product_smz
- dot_product_top6, dot_product_top5, dot_product_top4
- dot_product_smz_top6, dot_product_smz_top5, dot_product_smz_top4
- fragment_coverage, sequence_coverage, consecutive_ions
- explained_intensity, elution_weighted_cosine

Mass accuracy (3): mass_accuracy_deviation_mean, abs_mass_accuracy_deviation_mean, mass_accuracy_std
RT deviation (2): rt_deviation, abs_rt_deviation
MS1 features (2): ms1_precursor_coelution, ms1_isotope_cosine
Peptide properties (2): peptide_length, missed_cleavages
Tukey median polish (3): median_polish_cosine, median_polish_rsquared, median_polish_residual_ratio

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

1. Train model on all data (all files combined)
2. Score all entries with trained model
3. Compute run-level q-values per file (precursor + peptide level)
4. Select best observation per precursor across experiment
5. Compute experiment-level q-values (precursor + peptide level)
6. Apply `max(precursor, peptide)` at both run and experiment levels
7. Determine passing precursors: `experiment_qvalue <= threshold`
8. Include ALL per-file observations for passing precursors in output
9. Propagate best experiment_qvalue to all observations

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
