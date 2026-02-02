# CLAUDE.md - Development Context for Osprey

This file provides context for Claude Code when working on the Osprey project.

## Project Overview

Osprey is a Rust-based peptide-centric DIA (Data-Independent Acquisition) analysis tool for mass spectrometry proteomics. It uses ridge regression to deconvolute mixed MS/MS spectra and integrates with Skyline for quantification.

## Architecture

### Workspace Structure

```
osprey/
├── Cargo.toml                    # Workspace root
├── crates/
│   ├── osprey-core/              # Core types, config, errors, traits
│   ├── osprey-io/                # File I/O (mzML, blib, DIA-NN TSV)
│   ├── osprey-regression/        # Spectrum regression engine
│   ├── osprey-chromatography/    # Peak detection, EMG fitting
│   ├── osprey-scoring/           # Feature extraction, decoy generation
│   ├── osprey-fdr/               # FDR control, q-value calculation
│   └── osprey/                   # Main library + CLI binary
```

### Key Components

- **osprey-core**: `LibraryEntry`, `Spectrum`, `OspreyConfig`, `IsolationWindow`
- **osprey-io**: `MzmlReader` (uses mzdata crate), `DiannTsvLoader`, `ElibLoader`, `BlibLoader`, `BlibWriter`
- **osprey-regression**: `RidgeSolver` (Cholesky), `Binner`, `DesignMatrixBuilder`, `SparseRidgeSolver`, `SparseMatrixBuilder` (HRAM)
- **osprey-fdr**: `FdrController`, `MokapotRunner` (semi-supervised FDR control)

## Build Commands

```bash
# Development build
cargo build

# Run tests
cargo test

# Run with arguments
cargo run -- -i sample.mzML -l library.tsv -o results.blib

# Install globally
cargo install --path crates/osprey

# Release build
cargo build --release
```

## Key Files

- `osprey_specification_v1.2.md` - Full specification document
- `crates/osprey/src/pipeline.rs` - Main analysis pipeline
- `crates/osprey/src/main.rs` - CLI entry point
- `crates/osprey-io/src/output/blib.rs` - BiblioSpec output writer
- `crates/osprey-core/src/config.rs` - YAML configuration support

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

## Current Status

### Implemented (Phase 1)

- mzML parsing
- DIA-NN TSV library loading
- EncyclopeDIA elib library loading
- BiblioSpec blib library loading
- Ridge regression with Cholesky solver
- Unit resolution binning
- HRAM sparse matrix support (ppm-based peak matching)
- Basic peak detection
- blib output with Osprey extension tables
- YAML configuration
- CLI with all core options
- Mokapot integration (PIN file generation, result parsing)

### TODO (Phase 2)

- EMG peak fitting
- Full 30+ feature extraction
- Two-step search strategy
- Background correction

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
