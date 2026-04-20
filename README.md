# Osprey

[![CI](https://github.com/maccoss/osprey/actions/workflows/ci.yml/badge.svg)](https://github.com/maccoss/osprey/actions/workflows/ci.yml)
[![Release](https://github.com/maccoss/osprey/actions/workflows/release.yml/badge.svg)](https://github.com/maccoss/osprey/actions/workflows/release.yml)
[![Latest Release](https://img.shields.io/github/v/release/maccoss/osprey?display_name=tag&sort=semver)](https://github.com/maccoss/osprey/releases/latest)
[![License: Apache 2.0](https://img.shields.io/github/license/maccoss/osprey)](LICENSE)
[![Rust 1.75+](https://img.shields.io/badge/rust-1.75%2B-blue.svg)](https://www.rust-lang.org)
[![GitHub Downloads](https://img.shields.io/github/downloads/maccoss/osprey/total)](https://github.com/maccoss/osprey/releases)
[![Last Commit](https://img.shields.io/github/last-commit/maccoss/osprey)](https://github.com/maccoss/osprey/commits)
[![GitHub Stars](https://img.shields.io/github/stars/maccoss/osprey?style=flat)](https://github.com/maccoss/osprey/stargazers)

**Open-Source Peptide-Centric Search Tool Designed for Skyline Integration**

Osprey is an open-source tool for peptide detection and quantification in data-independent acquisition (DIA) mass spectrometry data. It uses fragment XIC co-elution analysis to detect peptides in DIA data, with machine learning scoring and rigorous FDR control.

## Features

- **Peptide-centric analysis**: Directly scores peptide candidates against observed spectra
- **Fragment co-elution**: Extracts fragment XICs, computes pairwise correlations, and scores using spectral matching at the apex
- **High-resolution support**: ppm-based fragment matching for high-resolution instruments (Astral, Orbitrap)
- **RT calibration**: LOESS-based retention time calibration with stratified sampling
- **Decoy generation**: Enzyme-aware sequence reversal with fragment m/z recalculation
- **Skyline integration**: Outputs BiblioSpec (.blib) format for seamless quantification in Skyline
- **FDR control**: Built-in Percolator-style semi-supervised SVM with two-level FDR (run + experiment) at precursor and peptide levels
- **Protein FDR**: Native protein parsimony with true picked-protein FDR (Savitski 2015), shared peptide handling (All/Razor/Unique), and two-pass architecture (first-pass protein q-values gate compaction and reconciliation consensus; second-pass q-values are authoritative)
- **Cross-run reconciliation**: Consensus RT alignment across replicates with peak boundary imputation for missing detections; tolerance derived from within-peptide RT reproducibility (typically 3-5x tighter than cross-peptide calibration MAD)
- **Peptide trace diagnostics**: `OSPREY_TRACE_PEPTIDE=<sequence>` emits a detailed per-peptide log at every pipeline stage (CWT candidates, consensus, reconciliation, gap-fill, FDR) for investigating integration quality on specific peptides. See [docs/17-peptide-trace.md](docs/17-peptide-trace.md).
- **Flexible input**: Supports DIA-NN TSV, EncyclopeDIA elib, and BiblioSpec blib libraries
- **Scalable**: Disk-backed memory architecture with automatic cache invalidation; tested on 240-file experiments
- **Configurable**: YAML configuration files with CLI argument overrides for reproducible analyses

## Installation

### Prerequisites

- Rust 1.75 or later
- OpenBLAS development libraries
- CMake

### Linux (Ubuntu/Debian)

#### 1. Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

Follow the on-screen prompts (the defaults are fine). After installation, `cargo` and `rustc` will be available in your terminal.

#### 2. Install build dependencies

```bash
sudo apt-get update
sudo apt-get install build-essential libopenblas-dev cmake pkg-config
```

#### 3. Build and install

```bash
git clone https://github.com/maccoss/osprey.git
cd osprey
cargo build --release
cargo install --path crates/osprey
```

### Windows

#### 1. Install Rust

Download and run [rustup-init.exe](https://rustup.rs/). The default installation (MSVC toolchain) is recommended.

#### 2. Install Visual Studio Build Tools

Rust on Windows requires the MSVC C/C++ build tools. If you don't have Visual Studio installed, download [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/) and install the **"Desktop development with C++"** workload.

#### 3. Install CMake

Download and install [CMake](https://cmake.org/download/). During installation, select **"Add CMake to the system PATH"**.

Alternatively, install via winget:

```powershell
winget install Kitware.CMake
```

#### 4. Install OpenBLAS

Osprey requires OpenBLAS for linear algebra operations. Download pre-built binaries from the [OpenBLAS releases](https://github.com/OpenMathLib/OpenBLAS/releases) page (choose the latest `OpenBLAS-*-x64.zip`).

Extract the archive and set the environment variable so the build system can find it:

```powershell
# Example: extracted to C:\OpenBLAS
[System.Environment]::SetEnvironmentVariable("OPENBLAS_PATH", "C:\OpenBLAS", "User")
```

Alternatively, install OpenBLAS via [vcpkg](https://vcpkg.io/):

```powershell
git clone https://github.com/microsoft/vcpkg.git
cd vcpkg
.\bootstrap-vcpkg.bat
.\vcpkg install openblas:x64-windows

# Set the environment variable to the installed path
[System.Environment]::SetEnvironmentVariable("OPENBLAS_PATH", "$PWD\installed\x64-windows", "User")
```

After setting environment variables, **restart your terminal** for changes to take effect.

#### 5. Build and install

```powershell
git clone https://github.com/maccoss/osprey.git
cd osprey
cargo build --release
cargo install --path crates/osprey
```

The built binary will be at `target\release\osprey.exe`, and `cargo install` places it in `%USERPROFILE%\.cargo\bin\` (which rustup adds to PATH automatically).

### macOS

#### 1. Install Rust

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

#### 2. Install build dependencies

```bash
brew install openblas cmake
export OPENBLAS_PATH="$(brew --prefix openblas)"
```

Add the `OPENBLAS_PATH` export to your `~/.zshrc` (or `~/.bash_profile`) so it persists across sessions.

#### 3. Build and install

```bash
git clone https://github.com/maccoss/osprey.git
cd osprey
cargo build --release
cargo install --path crates/osprey
```

### Optional: Mokapot

Osprey includes a built-in Percolator-style semi-supervised SVM — no Python required. Mokapot is an optional alternative FDR engine that can be enabled with `--fdr-engine mokapot`.

```bash
# Requires Python 3.8+ and pip
# Note: mokapot 0.10.0 requires pandas 2.x and numpy <2.0
pip install mokapot 'pandas>=2.0,<3.0' 'numpy<2.0'
```

If the `mokapot` command is not found after installation, add the scripts directory to your PATH:

- **Linux**: Add `~/.local/bin` to PATH in `~/.bashrc`
- **Windows**: `pip install --user` places scripts in `%APPDATA%\Python\PythonXX\Scripts` — add this to your system PATH, or use `pip install` (without `--user`) in an activated virtual environment
- **macOS**: Add `~/Library/Python/X.Y/bin` to PATH

## Quick Start

### Basic usage

```bash
# Analyze DIA data with a spectral library (default: 10 ppm fragment tolerance)
osprey -i sample.mzML -l library.tsv -o results.blib

# Multiple input files
osprey -i *.mzML -l library.tsv -o results.blib

# With TSV report
osprey -i sample.mzML -l library.tsv -o results.blib --report results.tsv

# Custom fragment tolerance (e.g., 20 ppm)
osprey -i sample.mzML -l library.tsv -o results.blib --fragment-tolerance 20

# Unit resolution mode (uses Th/Dalton tolerances instead of ppm)
osprey -i sample.mzML -l library.tsv -o results.blib --resolution unit

# Filter output at precursor-level FDR only (less conservative than default)
osprey -i sample.mzML -l library.tsv -o results.blib --fdr-level precursor

# Protein-level FDR always runs at 0.01 by default. Change threshold or
# use razor peptide assignment:
osprey -i *.mzML -l library.tsv -o results.blib --protein-fdr 0.05 --shared-peptides razor

# Filter blib output to only peptides from protein groups passing protein FDR:
osprey -i *.mzML -l library.tsv -o results.blib --fdr-level protein

# Trace a specific peptide's journey through the pipeline (zero overhead when unset)
OSPREY_TRACE_PEPTIDE=PEPTIDEK osprey -i *.mzML -l library.tsv -o results.blib --verbose 2>&1 | tee trace.log
grep '\[trace\]' trace.log
```

### Using configuration files

```bash
# Generate a template configuration
osprey --generate-config my_analysis.yaml

# Edit the configuration file, then run
osprey --config my_analysis.yaml

# Override specific settings
osprey --config my_analysis.yaml --run-fdr 0.005
```

### Example configuration (YAML)

```yaml
input_files:
  - sample1.mzML
  - sample2.mzML

library_source:
  DiannTsv: library.tsv
  # Or: Blib: library.blib
  # Or: Elib: library.elib

output_blib: results.blib
output_report: results.tsv

# RT Calibration
rt_calibration:
  enabled: true
  loess_bandwidth: 0.3        # Fraction of data for local fits (0.2-0.5)
  min_calibration_points: 200
  rt_tolerance_factor: 3.0    # Multiplier for residual SD
  fallback_rt_tolerance: 2.0  # Used if calibration fails

# FDR control
run_fdr: 0.01
experiment_fdr: 0.01
decoy_method: Reverse  # Options: Reverse, Shuffle, FromLibrary
decoys_in_library: false
fdr_level: Precursor    # Output filtering level: Precursor (default), Peptide, Protein, or Both

# Protein-level FDR (always runs)
protein_fdr: 0.01         # Protein-level FDR threshold (default 0.01)
shared_peptides: All       # How to handle shared peptides: All (default), Razor, or Unique
```

## Command-line Options

```text
Options:
  -c, --config <CONFIG>
          Configuration file (YAML format)

      --generate-config <FILE>
          Generate a template configuration file

  -i, --input <INPUT>...
          Input mzML file(s)

  -l, --library <LIBRARY>
          Spectral library file (.tsv for DIA-NN, .blib, or .elib)

  -o, --output <OUTPUT>
          Output results file (.blib format for Skyline)

      --resolution <MODE>
          Resolution mode: unit, hram, auto [default: auto]
          - auto: Use config defaults (ppm tolerances)
          - unit: Use Th (Dalton) tolerances for low-res data
          - hram: Use ppm tolerances for high-res data

      --fragment-tolerance <VALUE>
          Fragment m/z tolerance (e.g., 10 for 10 ppm, or 0.3 for 0.3 Th)

      --fragment-unit <UNIT>
          Fragment tolerance unit: ppm, mz (Thompson)

      --precursor-tolerance <VALUE>
          Precursor m/z tolerance (e.g., 10 for 10 ppm, or 1.0 for 1.0 Th)

      --precursor-unit <UNIT>
          Precursor tolerance unit: ppm, mz (Thompson)

      --rt-tolerance <MINUTES>
          RT tolerance in minutes (fallback when calibration disabled) [default: 2]

      --no-rt-calibration
          Disable RT calibration (use fixed rt_tolerance instead)

      --run-fdr <THRESHOLD>
          Run-level FDR threshold [default: 0.01]

      --experiment-fdr <THRESHOLD>
          Experiment-level FDR threshold (for multi-file analyses) [default: 0.01]

      --fdr-method <METHOD>
          FDR method: percolator (built-in SVM, default), mokapot (external Python),
          or simple (no ML) [default: percolator]

      --fdr-level <LEVEL>
          FDR filtering level for output: precursor, peptide, protein, or both [default: precursor]
          - precursor: filter on precursor-level q-values (peptide + charge)
          - peptide: filter on peptide-level q-values (peptide sequence only)
          - protein: filter on protein-level q-values
          - both: filter on max(precursor, peptide) q-values (most conservative)
          Note: the blib output always enforces precursor-level FDR within each
          eligible peptide regardless of this setting; --fdr-level controls which
          peptide identities are admitted.

      --protein-fdr <THRESHOLD>
          Protein-level FDR threshold [default: 0.01]. Protein parsimony and
          picked-protein FDR always run; this sets the cutoff for reporting
          and for the compaction/reconciliation rescue rule. Use --fdr-level
          protein to restrict the blib output to protein-FDR-passing groups.

      --shared-peptides <MODE>
          How to handle peptides mapping to multiple protein groups: all, razor,
          or unique [default: all]
          - all: shared peptides contribute to all their protein groups
          - razor: shared peptides assigned to the group with most unique peptides
          - unique: only unique peptides used; shared peptides excluded

      --write-pin
          Write PIN files for external tools

      --no-prefilter
          Disable the coelution signal pre-filter

      --threads <N>
          Number of threads (default: all available)

      --report <FILE>
          Write TSV report to this file

  -v, --verbose
          Verbose output

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```

## Output

### BiblioSpec (.blib)

The primary output is a SQLite database in BiblioSpec format, compatible with Skyline. It includes:

- **RefSpectra**: Detected peptides with consensus spectra
- **RefSpectraPeaks**: Fragment m/z and intensities
- **Modifications**: Modification positions and masses
- **Proteins/RefSpectraProteins**: Protein accessions and mappings
- **RetentionTimes**: Per-file peak boundaries with nullable retentionTime for Skyline ID lines
- **OspreyPeakBoundaries**: Peak boundaries (StartRT, EndRT, ApexRT) per run
- **OspreyRunScores**: Run-level q-values and scores
- **OspreyExperimentScores**: Experiment-level q-values

Import directly into Skyline for quantification. See [BiblioSpec Output Schema](docs/08-blib-output-schema.md) for complete schema documentation.

### TSV Report

Optional human-readable report with peptide detections, scores, and peak boundaries.

## Algorithm Overview

### Per-File RT Calibration

Each file gets **independent calibration** because LC conditions may vary between files:

1. Sample library peptides (default 100K targets, 2D stratified on a (RT × precursor m/z) grid for uniform gradient coverage)
2. Assume linear relationship: library RT range ≈ mzML RT range
3. Wide initial tolerance (20-30% of gradient range)
4. Score peptides and detect peaks via co-elution search
5. Record (library_RT, measured_apex_RT) pairs
6. Fit LOESS calibration curve using classical Cleveland (1979) robust iterations — residuals are refreshed at each iteration so the bisquare weights progressively tighten (set `OSPREY_LOESS_CLASSICAL_ROBUST=0` for the legacy single-refresh behavior)
7. Calculate residual statistics — tight RT tolerance (3× MAD × 1.4826) used for the main search

Calibration results are saved to JSON per file for reuse on subsequent runs. See [Calibration](docs/02-calibration.md) for the full algorithm including two-pass refinement.

### Candidate Selection

Candidate selection uses:

- **Isolation window**: From mzML file (defines which precursors are fragmented)
- **RT tolerance**: Calibrated (Phase 2) or fallback (if calibration disabled/fails)

Note: No additional precursor m/z tolerance is applied - the isolation window from the mzML is sufficient.

### Decoy Generation

Following the [pyXcorrDIA](https://github.com/maccoss/pyXcorrDIA) approach:

1. **Enzyme-aware terminal preservation**:
   - Trypsin/Lys-C: Preserves C-terminal residue only (the cleavage site), reverses the rest
   - Lys-N/Asp-N: Preserves N-terminal residue only, reverses the rest
   - Example: `PEPTIDEK` (trypsin) → `EDITPEPK`

2. **Fragment m/z recalculation**:
   - b-ions become y-ions at complementary positions
   - y-ions become b-ions at complementary positions
   - Full mass recalculation with proper ion formulas

3. **Precursor mass preserved**: Same amino acid composition means same precursor m/z

### Coelution Scoring

For each candidate precursor:

1. Extract fragment XICs across the chromatographic dimension
2. Compute pairwise fragment co-elution correlations
3. Detect candidate peaks via CWT consensus (Mexican Hat wavelet across transitions)
4. Rank candidates by `coelution_score × rt_penalty × log(1 + apex_intensity)`:
   - `coelution_score`: mean pairwise Pearson correlation of fragment XICs within the peak
   - `rt_penalty`: Gaussian centered on calibration-predicted RT, sigma = 5 × MAD × 1.4826 (downweights peaks far from expected RT without eliminating them)
   - `log(1 + apex_intensity)`: intensity tiebreaker so the main peak beats its own shoulder when coelution scores are comparable
5. Score the selected peak using spectral matching (dot product, XCorr, hyperscore) at the apex
6. Extract 21 features per precursor for machine learning scoring via built-in Percolator SVM (or optionally Mokapot)

### Cross-Run Peak Reconciliation

For multi-file experiments, Osprey reconciles peak integration boundaries across replicates after the initial FDR pass:

1. Collect peptides passing run-level precursor FDR across all files (per-entry precursor q-value is a hard precondition; protein FDR can upgrade borderline peptide-level evidence but cannot override poor precursor evidence)
2. Compute consensus library RTs using weighted median of per-run detections, weighted by `sigmoid(SVM score)` so wrong-peak detections with negative scores are downweighted
3. Refit per-run LOESS calibration using consensus peptides
4. Derive `rt_tolerance` from the global median of per-peptide library-RT MADs (within-peptide reproducibility floor, typically 3-5x tighter than cross-peptide calibration MAD)
5. For each entry in each run: keep the existing peak, switch to an alternate CWT candidate at the consensus RT, or perform forced integration
6. Re-score reconciled entries and compute final experiment-level q-values

This ensures consistent quantification across replicates by aligning peak boundaries to the same chromatographic feature. See [Cross-Run Reconciliation](docs/10-cross-run-reconciliation.md) for algorithm details.

### Disk-Backed Memory Architecture

Osprey uses per-file Parquet caching to scale to 1000+ file experiments without running out of memory. After scoring each file, the full scored entries are written to ZSTD-compressed Parquet caches and replaced in memory with lightweight FdrEntry stubs (~80 bytes each). Heavy data (features, fragments, CWT candidates) is reloaded on-demand from disk only when needed. Peptide sequence strings are deduplicated using `Arc<str>` interning across files.

Each Parquet cache stores SHA-256 hashes of search parameters, library identity, and reconciliation state in its file metadata. On subsequent runs, Osprey validates these hashes and automatically invalidates stale caches when parameters change. SVM discriminant scores from Percolator are persisted in lightweight sidecar files, enabling Osprey to skip SVM training entirely on reruns with matching parameters.

See [Pipeline Overview](docs/README.md) for memory architecture details and [Intermediate File Formats](docs/12-intermediate-files.md) for cache file documentation.

## Development

### Project Structure

```text
osprey/
├── crates/
│   ├── osprey-core/          # Core types and configuration
│   ├── osprey-io/            # File I/O (mzML, blib, libraries)
│   ├── osprey-chromatography/# Peak detection, RT calibration
│   ├── osprey-scoring/       # Feature extraction, decoy generation
│   ├── osprey-fdr/           # FDR control (Percolator SVM + optional Mokapot)
│   ├── osprey-ml/            # Machine learning (SVM, PEP estimation)
│   └── osprey/               # CLI and pipeline
```

### Key Components

- **osprey-core**: `LibraryEntry`, `Spectrum`, `OspreyConfig`, `IsolationWindow`, `FdrEntry`
- **osprey-io**: `MzmlReader`, `DiannTsvLoader`, `ElibLoader`, `BlibLoader`, `BlibWriter`
- **osprey-chromatography**: `PeakDetector`, `RTCalibrator`, `RTStratifiedSampler`
- **osprey-scoring**: `DecoyGenerator`, `FeatureExtractor`
- **osprey-fdr**: `FdrController`, `MokapotRunner`
- **osprey-ml**: `LinearSvmClassifier`, `PepEstimator` (peptide-level PEP only; protein PEP is intentionally not computed — see `docs/07-fdr-control.md`)

### Build and test

A Makefile provides convenience targets that run formatting and linting before each step:

```bash
make check    # Format and lint (cargo fmt + clippy)
make test     # Format, lint, and run tests
make build    # Format, lint, and build release
make install  # Format, lint, build release, and install
```

### Building documentation

```bash
cargo doc --open
```

## Current Status

### Implemented

- mzML parsing via mzdata crate
- DIA-NN TSV library loading
- EncyclopeDIA elib library loading
- BiblioSpec blib library loading
- Fragment XIC co-elution analysis with ppm-based fragment matching
- CWT consensus peak detection (Mexican Hat wavelet, median across transitions)
- RT-penalized, intensity-weighted peak selection: score = coelution × Gaussian RT penalty × log(1 + apex_intensity); interferers at the wrong RT lose via the RT penalty, narrow low-intensity shoulders lose via the log-intensity factor
- Tukey median polish for robust peak boundaries and fragment scoring
- Signal pre-filter for ~30% speedup on HRAM data (disable with `--no-prefilter`)
- RT calibration with LOESS regression (classical Cleveland 1979 robust iterations) and (RT × precursor m/z) stratified sampling
- MS1/MS2 mass calibration
- Enzyme-aware decoy generation with fragment recalculation
- 21-feature extraction per precursor for machine learning scoring
- Built-in Percolator-style semi-supervised SVM for FDR control (no Python required)
- Two-level FDR control (run + experiment level) with dual precursor + peptide level FDR
- Two-stage blib output gate: `--fdr-level` determines eligible peptide identities; precursor-level FDR is always enforced within each eligible peptide (fallback to best charge state if no charge passes)
- Native protein parsimony with true picked-protein FDR (Savitski 2015), two-pass architecture, and protein-aware compaction + reconciliation consensus rescue
- Optional Mokapot integration with parallel workers and progress streaming
- Cross-run peak reconciliation: aligns integration boundaries across replicates using consensus RTs; sigma-clipped MAD + original-calibration cap prevent tolerance inflation from wrong-peak detections
- Multi-file observation propagation: all per-file observations for passing precursors included in output
- Disk-backed memory architecture with per-file Parquet caching for 1000+ file experiments
- Intelligent cache invalidation via SHA-256 parameter hashing in Parquet metadata
- FDR score sidecar caching — skips Percolator SVM training on reruns with same parameters
- Streaming blib output via lightweight plan entries — avoids loading full entries into memory
- Binary spectra cache for faster second-pass mzML loading
- Cross-implementation bisection diagnostics for validation against OspreySharp (env-var gated; zero overhead when unused)
- Fragment co-elution scoring
- Mass accuracy features (ppm-level)
- MS1 isotope envelope and precursor co-elution features
- blib output with Osprey extension tables (library theoretical fragments)
- YAML configuration
- CLI with all core options
- Calibration JSON save/load and HTML report generation

### TODO

- EMG peak fitting (Levenberg-Marquardt)
- Feature weight optimization (remove low-value features)
- Background correction
- Two-step search strategy
- Iterative candidate expansion

## Citation

If you use Osprey in your research, please cite:

> [Citation pending publication]

## License

Apache 2.0

## Acknowledgments

Developed by the [MacCoss Lab](https://maccosslab.org) at the University of Washington.

## Related Projects

- [Skyline](https://skyline.ms) - Targeted mass spectrometry environment
- [DIA-NN](https://github.com/vdemichev/DiaNN) - DIA data analysis
- [EncyclopeDIA](https://bitbucket.org/searleb/encyclopedia) - Library searching for DIA
- [pyXcorrDIA](https://github.com/maccoss/pyXcorrDIA) - Python DIA analysis tool
- [mokapot](https://github.com/wfondrie/mokapot) - Semi-supervised FDR control
