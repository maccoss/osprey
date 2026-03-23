# Osprey

[![GitHub Downloads](https://img.shields.io/github/downloads/maccoss/osprey/total)](https://github.com/maccoss/osprey/releases)

**Open-Source Peptide-Centric Search Tool Designed for Skyline Integration**

Osprey is an open-source tool for peptide detection and quantification in data-independent acquisition (DIA) mass spectrometry data. It uses fragment XIC co-elution analysis to detect peptides in DIA data, with machine learning scoring and rigorous FDR control.

## Features

- **Peptide-centric analysis**: Directly scores peptide candidates against observed spectra
- **Fragment co-elution**: Extracts fragment XICs, computes pairwise correlations, and scores using spectral matching at the apex
- **High-resolution support**: ppm-based fragment matching for high-resolution instruments (Astral, Orbitrap)
- **RT calibration**: LOESS-based retention time calibration with stratified sampling
- **Decoy generation**: Enzyme-aware sequence reversal with fragment m/z recalculation
- **Skyline integration**: Outputs BiblioSpec (.blib) format for seamless quantification in Skyline
- **FDR control**: Built-in Percolator-style semi-supervised SVM with target-decoy competition; optional Mokapot integration for an alternative semi-supervised approach
- **Flexible input**: Supports DIA-NN TSV, EncyclopeDIA elib, and BiblioSpec blib libraries
- **Scalable**: Disk-backed memory architecture processes 1000+ files without running out of memory
- **Configurable**: YAML configuration files for reproducible analyses

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

# Candidate selection (precursor filtering uses isolation window from mzML)
max_candidates_per_spectrum: 5250  # Use 0 for unlimited

# RT Calibration
rt_calibration:
  enabled: true
  loess_bandwidth: 0.3     # Fraction of data for local fits (0.2-0.5)
  n_rt_bins: 10            # Bins for stratified sampling
  peptides_per_bin: 100    # Peptides to sample per bin
  min_calibration_points: 50
  rt_tolerance_factor: 3.0 # Multiplier for residual SD
  fallback_rt_tolerance: 2.0  # Used if calibration fails

# FDR control
run_fdr: 0.01
experiment_fdr: 0.01
decoy_method: Reverse  # Options: Reverse, Shuffle, FromLibrary
decoys_in_library: false
```

## Command-line Options

```
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

      --max-candidates <N>
          Maximum candidates per spectrum (use 0 for unlimited) [default: 5250]

      --run-fdr <THRESHOLD>
          Run-level FDR threshold [default: 0.01]

      --experiment-fdr <THRESHOLD>
          Experiment-level FDR threshold (for multi-file analyses) [default: 0.01]

      --fdr-method <METHOD>
          FDR method: percolator (built-in SVM, default), mokapot (external Python),
          or simple (no ML) [default: percolator]

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
- **OspreyPeakBoundaries**: Peak boundaries (StartRT, EndRT, ApexRT) per run
- **OspreyRunScores**: Run-level q-values and scores
- **OspreyExperimentScores**: Experiment-level q-values

Import directly into Skyline for quantification. See [BiblioSpec Output Schema](docs/08-blib-output-schema.md) for complete schema documentation.

### TSV Report

Optional human-readable report with peptide detections, scores, and peak boundaries.

## Algorithm Overview

### Multi-File RT Calibration

When processing multiple files, Osprey uses a **multi-file calibration strategy**:

#### File 1: Calibration Discovery

1. Use **all library peptides** (no sampling)
2. Assume linear relationship: library RT range ≈ mzML RT range
3. Wide initial tolerance (20-30% of gradient range)
4. Score peptides and detect peaks
5. Record (library_RT, measured_apex_RT) pairs
6. Fit LOESS calibration curve
7. Calculate residual standard deviation

**Why all peptides?** For the first file, we don't know the RT shift, so we use all peptides with wide tolerance to maximize calibration points. The assumption is that library RTs span roughly the same range as the measured gradient.

#### Files 2-N: Calibrated Search

1. **Reuse calibration** from File 1 (same experiment → similar LC conditions)
2. Use tight RT tolerance (3× residual SD from File 1)
3. Apply LOESS to convert library RTs to predicted measured RTs
4. Run full search with calibrated candidate selection
5. Extract fragment XICs and detect peaks

**Why reuse calibration?** Files within the same experiment have similar LC conditions (same column, gradient, mobile phases). The calibration from File 1 is directly applicable to subsequent files.

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
3. Detect peaks in the co-elution signal
4. Score using spectral matching (dot product, XCorr, hyperscore) at the apex
5. Extract 21 features per precursor for machine learning scoring via built-in Percolator SVM (or optionally Mokapot)

## Development

### Project Structure

```
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
- **osprey-ml**: `LinearSvmClassifier`, `PepEstimator`

### Running tests

```bash
cargo test
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
- Peak detection with Tukey median polish peak boundaries
- RT calibration with LOESS regression and stratified sampling
- MS1/MS2 mass calibration
- Enzyme-aware decoy generation with fragment recalculation
- 21-feature extraction per precursor for machine learning scoring
- Built-in Percolator-style semi-supervised SVM for FDR control (no Python required)
- Two-level FDR control (run + experiment level) with optional Mokapot integration
- Tukey median polish for robust elution profiles and fragment scoring
- Fragment co-elution scoring, elution-weighted cosine
- Mass accuracy features (ppm-level)
- MS1 isotope envelope and precursor co-elution features
- blib output with Osprey extension tables
- YAML configuration
- CLI with all core options
- Calibration JSON save/load and HTML report generation

### TODO

- Two-step search strategy
- Background correction
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
