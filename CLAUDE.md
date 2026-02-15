# CLAUDE.md - Development Context for Osprey

This file provides context for Claude Code when working on the Osprey project.

## Project Overview

Osprey is a Rust-based peptide-centric DIA (Data-Independent Acquisition) analysis tool for mass spectrometry proteomics. It supports two search modes: ridge regression (spectrum deconvolution) and coelution (DIA-NN-style fragment XIC correlation). Results integrate with Skyline for quantification.

**Current Status**: Working prototype with 76,169 precursors at 1% FDR from single Astral file.

## Architecture

### Workspace Structure

```
osprey/
├── Cargo.toml                    # Workspace root
├── crates/
│   ├── osprey-core/              # Core types, config, errors, traits
│   ├── osprey-io/                # File I/O (mzML, blib, DIA-NN TSV)
│   ├── osprey-regression/        # Spectrum regression engine
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

- **osprey-core**: `LibraryEntry`, `Spectrum`, `OspreyConfig`, `IsolationWindow`, `FeatureSet`, `CoelutionFeatureSet`, `SearchMode`
- **osprey-io**: `MzmlReader` (uses mzdata crate), `DiannTsvLoader`, `ElibLoader`, `BlibLoader`, `BlibWriter`
- **osprey-regression**: `RidgeSolver` (Cholesky), `Binner`, `DesignMatrixBuilder`, `SparseRidgeSolver`, `SparseMatrixBuilder` (HRAM)
- **osprey-chromatography**: `PeakDetector`, `RTCalibrator`, `MzCalibration`, `CalibrationParams`
- **osprey-scoring**: `SpectralScorer`, `FeatureExtractor`, `DecoyGenerator`, batch scoring
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
- `docs/fdr-control.md` - Two-level FDR with Mokapot
- `docs/blib-output-schema.md` - BiblioSpec output format and Skyline integration
- `crates/osprey/src/pipeline.rs` - Main analysis pipeline
- `crates/osprey/src/main.rs` - CLI entry point
- `crates/osprey-fdr/src/mokapot.rs` - Mokapot integration
- `crates/osprey-io/src/output/blib.rs` - BiblioSpec blib writer
- `crates/osprey-core/src/types.rs` - FeatureSet (37 regression), CoelutionFeatureSet (45 coelution)
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
- Two search modes: regression (37 features) and coelution (45 features)
- Ridge regression with NNLS solver (f32)
- HRAM sparse matrix support (ppm-based matching)
- Coelution search: DIA-NN-style fragment XIC correlation without regression
- Peak detection in coefficient/XIC time series
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

### Search Modes

Osprey supports two search modes via `SearchMode` (`--search-mode`):

- **Regression** (default): Ridge regression deconvolutes mixed DIA spectra into individual peptide contributions. 37 features extracted per precursor.
- **Coelution**: DIA-NN-style fragment XIC correlation without regression. Extracts fragment XICs, computes pairwise correlations, and scores using spectral matching at the apex. 45 features extracted per precursor.

Both modes share: calibration phase, decoy generation, library deduplication, Mokapot FDR, and blib output. The blib always stores library theoretical fragment m/z and intensities (not observed DIA peaks).

### Regression Feature Set (37 Features for Mokapot PIN)

Ridge regression (8):
- peak_apex, peak_area, peak_width, coefficient_stability
- relative_coefficient, explained_intensity, signal_to_noise, xic_signal_to_noise

Spectral matching - mixed at apex (2): xcorr, consecutive_ions
Spectral matching - deconvoluted (2): xcorr_deconv, consecutive_ions_deconv
RT deviation (1): rt_deviation

Fragment co-elution (9):
- fragment_coelution_sum, fragment_coelution_min, n_coeluting_fragments
- fragment_corr_0..5 (per-fragment correlations, top 6 by library intensity)

Elution-weighted similarity (1): elution_weighted_cosine
Mass accuracy (3): mass_accuracy_deviation_mean, abs_mass_accuracy_deviation_mean, mass_accuracy_std

Percolator-style (6):
- abs_rt_deviation, peptide_length, missed_cleavages
- ln_num_candidates, coef_zscore, coef_zscore_mean

MS1 features (2): ms1_precursor_coelution, ms1_isotope_cosine

Tukey median polish (3):
- median_polish_cosine, median_polish_rsquared, median_polish_residual_ratio

### Coelution Feature Set (45 Features for Mokapot PIN)

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

## Recent Changes

- Added coelution search mode (DIA-NN-style fragment XIC correlation, 45 features)
- Fixed coelution blib output to write library theoretical fragments (not observed DIA peaks)
- Fixed experiment-level FDR: entries not in experiment results no longer fall back to run q-values
- Removed pearson_correlation and spearman_correlation from coelution features (poor discrimination)
- Added blib round-trip tests (fragment, modseq, protein mapping, multi-run RT)
- Replaced library-assisted median polish with Tukey median polish for peak boundaries
- Added 3 median polish features: cosine, R², residual_ratio
- Fixed all cosine/correlation scoring to include ALL library fragments (0 for unmatched)
- Reduced PIN features from 30 mixed to focused 37 features with better discrimination
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
