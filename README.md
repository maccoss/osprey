# Osprey

**Open-Source Peptide Recognition and Elution Yield**

Peptide-centric DIA analysis with Skyline integration.

Osprey is an open-source tool for peptide detection and quantification in data-independent acquisition (DIA) mass spectrometry data. It uses ridge regression to deconvolute mixed MS/MS spectra, aggregates evidence across the chromatographic dimension, and provides rigorous FDR control for peptide detections.

## Features

- **Peptide-centric analysis**: Directly scores peptide candidates against observed spectra
- **Ridge regression**: Handles overlapping isolation windows and co-eluting peptides
- **HRAM support**: Sparse matrix operations with ppm-based peak matching for high-resolution data
- **RT calibration**: LOESS-based retention time calibration with stratified sampling
- **Decoy generation**: Enzyme-aware sequence reversal with fragment m/z recalculation
- **Skyline integration**: Outputs BiblioSpec (.blib) format for seamless quantification in Skyline
- **FDR control**: Target-decoy competition with mokapot integration for semi-supervised learning
- **Flexible input**: Supports DIA-NN TSV, EncyclopeDIA elib, and BiblioSpec blib libraries
- **Configurable**: YAML configuration files for reproducible analyses

## Installation

### Prerequisites

- Rust 1.75 or later
- Python 3.8 or later (for mokapot)
- Mokapot (for FDR control)
- OpenBLAS development libraries
- OpenSSL development libraries
- CMake

On Ubuntu/Debian:
```bash
sudo apt-get update
sudo apt-get install libopenblas-dev libssl-dev cmake pkg-config python3-pip
```

Install mokapot:
```bash
# Note: mokapot 0.10.0 requires pandas 2.x and numpy <2.0
pip install mokapot 'pandas>=2.0,<3.0' 'numpy<2.0'

# Or install for your user only (recommended)
pip install --user mokapot 'pandas>=2.0,<3.0' 'numpy<2.0'
```

Verify mokapot installation:
```bash
mokapot --version
```

**Note:** If the `mokapot` command is not found after installation, add `~/.local/bin` to your PATH:
```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### Building from source

```bash
git clone https://github.com/maccoss/osprey.git
cd osprey
cargo build --release
```

### Installing

```bash
cargo install --path crates/osprey
```

## Quick Start

### Basic usage

```bash
# Analyze DIA data with a spectral library
osprey -i sample.mzML -l library.tsv -o results.blib

# Multiple input files
osprey -i *.mzML -l library.tsv -o results.blib

# With TSV report
osprey -i sample.mzML -l library.tsv -o results.blib --report results.tsv

# High-resolution mode (HRAM) - uses ppm-based fragment matching
osprey -i sample.mzML -l library.tsv -o results.blib --resolution hram --hram-tolerance 20
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

# Resolution mode (controls tolerance units and fragment matching)
resolution_mode: Auto  # Auto uses config defaults (ppm)
# Options: Auto, UnitResolution (Th units), or HRAM (ppm units)
# Ridge regression always uses ~1 Th bins; HRAM applies ppm matching in scoring

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

# Regression
regularization_lambda: CrossValidated
# Options: CrossValidated, Adaptive, or Fixed with value

# Two-step search (recommended)
two_step_search:
  enabled: true
  step1_fdr: 0.01

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
          Note: Ridge regression always uses ~1 Th bins for memory efficiency.
          HRAM precision is applied via ppm-based fragment matching in scoring.

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

      --lambda <VALUE>
          Fixed regularization parameter lambda

      --max-candidates <N>
          Maximum candidates per spectrum (use 0 for unlimited) [default: 5250]

      --run-fdr <THRESHOLD>
          Run-level FDR threshold [default: 0.01]

      --threads <N>
          Number of threads (default: all available)

      --report <FILE>
          Write TSV report to this file

      --export-coefficients
          Export coefficient time series to parquet file(s)

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

Import directly into Skyline for quantification. See [BiblioSpec Output Schema](docs/blib-output-schema.md) for complete schema documentation.

### TSV Report

Optional human-readable report with peptide detections, scores, and peak boundaries.

## Algorithm Overview

### Multi-File RT Calibration

When processing multiple files, Osprey uses a **multi-file calibration strategy**:

#### File 1: Calibration Discovery

1. Use **all library peptides** (no sampling)
2. Assume linear relationship: library RT range ≈ mzML RT range
3. Wide initial tolerance (20-30% of gradient range)
4. Run regression and detect peaks
5. Record (library_RT, measured_apex_RT) pairs
6. Fit LOESS calibration curve
7. Calculate residual standard deviation

**Why all peptides?** For the first file, we don't know the RT shift, so we use all peptides with wide tolerance to maximize calibration points. The assumption is that library RTs span roughly the same range as the measured gradient.

#### Files 2-N: Calibrated Search

1. **Reuse calibration** from File 1 (same experiment → similar LC conditions)
2. Use tight RT tolerance (3× residual SD from File 1)
3. Apply LOESS to convert library RTs to predicted measured RTs
4. Run full regression with calibrated candidate selection
5. Extract coefficient time series and detect peaks

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

### Ridge Regression

For each MS2 spectrum:
1. Select library peptides within isolation window and RT tolerance
2. Build design matrix from binned library spectra
3. Solve regularized least squares: minimize ‖Ax - b‖² + λ‖x‖²
4. Non-negative coefficients indicate peptide contributions

### Resolution Modes

**Ridge regression always uses unit resolution binning (~1 Th bins, ~2000 bins)** for memory efficiency. Using 0.02 Th bins (100K bins) would create impractically large matrices (~400MB per spectrum per thread).

Resolution mode controls where high-resolution precision is applied:

**Unit Resolution Mode** (`--resolution unit`):
- Fragment tolerances default to Th (Dalton) units
- Suitable for low-resolution instruments (e.g., Stellar)
- Fragment matching uses binned intensities

**HRAM Mode** (`--resolution hram`):
- Fragment tolerances default to ppm units
- Suitable for high-resolution instruments (e.g., Astral, Orbitrap)
- ppm-based fragment matching in:
  - Fragment pre-filter (has_top3_fragment_match)
  - Spectral scoring (dot product, XCorr)
  - Mass accuracy feature calculation
- Sparse matrix operations for efficient ppm-based scoring

**Auto Mode** (`--resolution auto`, default):
- Uses tolerance settings from config file (defaults to ppm)
- Does not auto-detect instrument type

## Development

### Project Structure

```
osprey/
├── crates/
│   ├── osprey-core/          # Core types and configuration
│   ├── osprey-io/            # File I/O (mzML, blib, libraries)
│   ├── osprey-regression/    # Ridge regression engine
│   ├── osprey-chromatography/# Peak detection, RT calibration
│   ├── osprey-scoring/       # Feature extraction, decoy generation
│   ├── osprey-fdr/           # FDR control, mokapot integration
│   └── osprey/               # CLI and pipeline
```

### Key Components

- **osprey-core**: `LibraryEntry`, `Spectrum`, `OspreyConfig`, `IsolationWindow`
- **osprey-io**: `MzmlReader`, `DiannTsvLoader`, `ElibLoader`, `BlibLoader`, `BlibWriter`
- **osprey-regression**: `RidgeSolver`, `Binner`, `DesignMatrixBuilder`, `SparseRidgeSolver`, `SparseMatrixBuilder`
- **osprey-chromatography**: `PeakDetector`, `RTCalibrator`, `RTStratifiedSampler`
- **osprey-scoring**: `DecoyGenerator`, `FeatureExtractor`
- **osprey-fdr**: `FdrController`, `MokapotRunner`

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
- Ridge regression with Cholesky solver (NNLS, f32)
- Unit resolution binning with HRAM sparse matrix support (ppm-based matching)
- Peak detection with Tukey median polish peak boundaries
- RT calibration with LOESS regression and stratified sampling
- MS1/MS2 mass calibration
- Enzyme-aware decoy generation with fragment recalculation
- 37-feature extraction per precursor for Mokapot
- Two-level FDR control (run + experiment level) via Mokapot
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
