# Osprey Algorithm Documentation

This folder contains detailed documentation of the algorithms used in Osprey.

## Contents

- [RT Calibration](rt-calibration.md) - LOESS-based retention time calibration with stratified sampling
- [Decoy Generation](decoy-generation.md) - Enzyme-aware sequence reversal for FDR control
- [Ridge Regression](ridge-regression.md) - Spectral deconvolution for peptide detection
- [Peak Detection](peak-detection.md) - Chromatographic peak detection in coefficient time series

## Overview

Osprey uses a peptide-centric approach to DIA analysis:

```
┌─────────────────────────────────────────────────────────────────┐
│                        Input                                     │
│   mzML files (DIA data) + Spectral Library (predicted spectra)  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Phase 1: RT Calibration                         │
│   Stratified sampling → Peak detection → LOESS fitting          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                  Decoy Generation                                │
│   Enzyme-aware reversal → Fragment m/z recalculation            │
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
