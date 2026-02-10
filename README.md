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
pip install mokapot

# Or install for your user only (recommended)
pip install --user mokapot
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

# High-resolution mode (HRAM)
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

# Resolution mode
resolution_mode: Auto
# Options: Auto, UnitResolution, or HRAM with tolerance_ppm
# resolution_mode:
#   HRAM:
#     tolerance_ppm: 20.0

# Candidate selection (precursor filtering uses isolation window from mzML)
max_candidates_per_spectrum: 200

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
  -c, --config <CONFIG>           Configuration file (YAML format)
      --generate-config <FILE>    Generate a template configuration file
  -i, --input <INPUT>...          Input mzML file(s)
  -l, --library <LIBRARY>         Spectral library file (.tsv, .blib, or .elib)
  -o, --output <OUTPUT>           Output results file (.blib)
      --resolution <MODE>         Resolution mode: auto, unit, hram [default: auto]
      --hram-tolerance <PPM>      HRAM tolerance in ppm [default: 20.0]
      --rt-tolerance <MINUTES>    Fallback RT tolerance [default: 2.0]
      --no-rt-calibration         Disable RT calibration
      --lambda <VALUE>            Fixed regularization parameter
      --max-candidates <N>        Maximum candidates per spectrum [default: 200]
      --run-fdr <THRESHOLD>       Run-level FDR threshold [default: 0.01]
      --threads <N>               Number of threads (default: all)
      --report <FILE>             Write TSV report
  -v, --verbose                   Verbose output
  -h, --help                      Print help
  -V, --version                   Print version
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

### HRAM Mode

For high-resolution data:
- Sparse matrix operations (CSC format via sprs crate)
- ppm-based peak matching instead of unit bins
- Conjugate gradient solver for large problems (>200 candidates)
- Dense Cholesky for smaller problems

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

### Implemented (Phase 1)

- mzML parsing via mzdata crate
- DIA-NN TSV library loading
- EncyclopeDIA elib library loading
- BiblioSpec blib library loading
- Ridge regression with Cholesky solver
- Unit resolution binning
- HRAM sparse matrix support (ppm-based peak matching)
- Basic peak detection
- RT calibration with LOESS and stratified sampling
- Enzyme-aware decoy generation with fragment recalculation
- blib output with Osprey extension tables
- YAML configuration
- CLI with all core options
- Mokapot integration (PIN file generation, result parsing)

### TODO (Phase 2)

- EMG peak fitting
- Full 30+ feature extraction
- Two-step search strategy
- Background correction

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
