# Osprey Algorithm Documentation

This folder contains detailed documentation of the algorithms used in Osprey.

## Contents

- [RT Calibration](rt-calibration.md) - LibCosine-based RT calibration with target-decoy competition
- [Mass Calibration](mass-calibration.md) - MS1 (precursor) and MS2 (fragment) mass calibration
- [Spectral Scoring](spectral-scoring.md) - LibCosine, XCorr, and isotope cosine scoring
- [XCorr Scoring](xcorr-scoring.md) - Detailed Comet-style XCorr implementation
- [Decoy Generation](decoy-generation.md) - Enzyme-aware sequence reversal for FDR control
- [Ridge Regression](ridge-regression.md) - Spectral deconvolution for peptide detection
- [Peak Detection](peak-detection.md) - Chromatographic peak detection in coefficient time series

## Overview

Osprey uses a peptide-centric approach to DIA analysis with a two-phase calibration and search strategy:

```
┌─────────────────────────────────────────────────────────────────┐
│                        Input                                     │
│   mzML files (DIA data) + Spectral Library (predicted spectra)  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│              Phase 1: Calibration Discovery                      │
│   1. Generate decoys (enzyme-aware reversal)                     │
│   2. Sample ~2000 peptides for calibration                       │
│   3. Score ALL spectra via LibCosine (pyXcorrDIA-compatible)     │
│   4. Extract MS1 isotope envelopes → mass calibration            │
│   5. Target-decoy competition → confident matches at 1% FDR      │
│   6. Fit RT calibration (LOESS) from confident targets           │
│   7. Calculate MS1/MS2 mass error statistics                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Phase 2: Full Search                            │
│   For each spectrum:                                             │
│     1. Select candidates (isolation window + calibrated RT)      │
│     2. Build design matrix from library spectra                  │
│     3. Solve ridge regression                                    │
│     4. Extract non-zero coefficients                             │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Peak Detection                                  │
│   Coefficient time series → Peak boundaries → Features           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  FDR Control                                     │
│   Target-decoy competition → q-value estimation                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        Output                                    │
│   BiblioSpec (.blib) for Skyline + TSV report                   │
└─────────────────────────────────────────────────────────────────┘
```

## Key Concepts

### Peptide-Centric Analysis

Unlike spectrum-centric approaches that identify each spectrum independently, Osprey:

1. **Considers all candidate peptides simultaneously** using ridge regression
2. **Aggregates evidence across the chromatographic dimension** before scoring
3. **Handles co-eluting peptides** in overlapping isolation windows

### Calibration Discovery (pyXcorrDIA-Compatible)

The calibration phase follows the pyXcorrDIA methodology:

1. **LibCosine Scoring**: Score library vs observed spectra using:
   - SMZ preprocessing: `sqrt(intensity) × mz²`
   - ppm-based fragment matching (closest m/z, not highest intensity)
   - L2 normalization and cosine similarity

2. **Target-Decoy Competition**: Each target competes against its paired decoy
   - Higher LibCosine score wins
   - Ties go to decoy (conservative)
   - Only targets passing 1% FDR used for calibration

3. **Mass Calibration**: From FDR-filtered matches:
   - MS1 (precursor) error: Extract M+0 peak from MS1 spectrum
   - MS2 (fragment) error: Collect matched fragment mass errors

4. **RT Calibration**: LOESS fitting of (library_RT, measured_RT) pairs

### Scoring Methods

Osprey provides multiple scoring methods for different use cases:

| Score | Use Case | Method |
|-------|----------|--------|
| **LibCosine** | Calibration, FDR | Cosine similarity with SMZ preprocessing |
| **XCorr** | Secondary score | Comet-style cross-correlation with windowing |
| **Isotope Cosine** | MS1 quality | Observed vs theoretical isotope pattern |

### Why Ridge Regression?

DIA spectra are mixed spectra containing fragments from multiple precursors. Ridge regression:

- Deconvolutes overlapping signals
- Provides coefficients proportional to peptide abundance
- Handles collinearity (similar spectra) gracefully
- Is computationally efficient

### Why RT Calibration?

Library retention times may be:
- Predicted by deep learning models (e.g., Prosit, DeepLC)
- Measured on different LC systems
- Normalized to iRT scale

RT calibration converts library RTs to expected measured RTs, enabling tight RT filtering.

## Crate Architecture

```
osprey/
├── osprey-core/           # Core types, config, errors, traits
│   ├── config.rs          # YAML configuration structures
│   ├── isotope.rs         # Isotope distribution calculations
│   ├── types.rs           # LibraryEntry, Spectrum, IsolationWindow
│   └── traits.rs          # MS1SpectrumLookup, etc.
│
├── osprey-io/             # File I/O
│   ├── library/           # DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib
│   └── mzml.rs            # mzML parsing (via mzdata crate)
│
├── osprey-scoring/        # Spectral scoring
│   ├── lib.rs             # SpectralScorer (LibCosine, XCorr)
│   ├── batch.rs           # BLAS-accelerated batch scoring
│   └── decoy.rs           # Decoy generation
│
├── osprey-chromatography/ # RT processing
│   ├── calibration/       # LOESS fitting, mass calibration
│   └── peak/              # Peak detection
│
├── osprey-regression/     # Ridge regression
│   ├── ridge.rs           # Cholesky solver
│   ├── sparse.rs          # Sparse HRAM solver
│   └── binning.rs         # m/z binning
│
├── osprey-fdr/            # FDR control
│   └── controller.rs      # Target-decoy competition, q-values
│
└── osprey/                # Main binary
    ├── main.rs            # CLI entry point
    └── pipeline.rs        # Analysis pipeline
```

## Configuration

Osprey supports YAML configuration files:

```bash
# Generate template
osprey --generate-config config.yaml

# Run with config
osprey --config config.yaml
```

Example configuration:

```yaml
# Input/Output
input_files:
  - sample1.mzML
  - sample2.mzML
spectral_library: library.tsv
output_blib: results.blib

# Resolution mode
resolution_mode:
  HRAM:
    tolerance_ppm: 10.0

# RT Calibration
rt_calibration:
  enabled: true
  loess_bandwidth: 0.3
  min_calibration_points: 50
  rt_tolerance_factor: 3.0

# FDR Control
run_fdr: 0.01

# Decoy Generation
decoy_method: Reverse
```

## Debug Output

Osprey writes a `calibration_debug.csv` with paired target-decoy scores for analysis:

```csv
target_entry_id,charge,target_sequence,decoy_sequence,
target_libcosine,decoy_libcosine,target_xcorr,decoy_xcorr,
target_isotope_score,decoy_isotope_score,
target_precursor_error_ppm,decoy_precursor_error_ppm,
target_rt,decoy_rt,library_rt,expected_rt,delta_rt,
target_n_matched,decoy_n_matched
```

This enables analysis of what features best separate targets from decoys.

## References

- pyXcorrDIA: https://github.com/maccoss/pyXcorrDIA
- Comet: https://comet-ms.sourceforge.io/
- LOESS: Cleveland, W.S. (1979). "Robust locally weighted regression and smoothing scatterplots"
