# Osprey Algorithm Documentation

This folder contains detailed documentation of the algorithms used in Osprey, a peptide-centric DIA analysis tool that uses ridge regression to deconvolute mixed MS/MS spectra.

## Current Status (Working Prototype)

Osprey has a **working prototype** that can:
- Parse mzML files with DIA data
- Load spectral libraries (DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib)
- Generate enzyme-aware decoys
- Run auto-calibration with target-decoy FDR control
- Perform ridge regression spectrum deconvolution
- Extract 30 features per precursor
- Run Mokapot semi-supervised FDR control
- Output BiblioSpec .blib files for Skyline

**Recent test run**: 76,169 precursors at 1% FDR from a single Astral file.

## Algorithm Documentation

| Document | Description |
|----------|-------------|
| [RT Calibration](rt-calibration.md) | LOESS-based RT calibration with target-decoy FDR |
| [Mass Calibration](mass-calibration.md) | MS1 (precursor) and MS2 (fragment) mass calibration |
| [Spectral Scoring](spectral-scoring.md) | LibCosine scoring (calibration phase) |
| [XCorr Scoring](xcorr-scoring.md) | Comet-style XCorr implementation |
| [Deconvolution Scoring](deconvolution-scoring.md) | Post-regression scoring: all scores computed after ridge regression |
| [Ridge Regression](ridge-regression.md) | NNLS ridge regression for spectrum deconvolution |
| [Peak Detection](peak-detection.md) | Chromatographic peak detection in coefficient time series |
| [Decoy Generation](decoy-generation.md) | Enzyme-aware sequence reversal for FDR control |
| [FDR Control](fdr-control.md) | Two-level FDR with Mokapot integration |

---

## Pipeline Overview

Osprey's pipeline has four major phases. The calibration phase uses fast spectral scoring to establish RT and mass calibration. The main search phase uses ridge regression to deconvolute mixed DIA spectra. Post-regression scoring computes features for semi-supervised FDR control via Mokapot.

```
INPUT
  mzML files (DIA data) + Spectral Library (predicted spectra)
  │
  ▼
PHASE 1: INITIALIZATION
  ├─ Load spectral library (DIA-NN TSV, elib, or blib)
  ├─ Generate decoys (enzyme-aware reversal)
  ├─ Pre-bin library spectra (f32, one-time cost)
  └─ Build m/z index for candidate lookup
  │
  ▼
PHASE 2: CALIBRATION DISCOVERY (first file only)
  ├─ Calculate library-to-measured RT mapping
  ├─ Set wide RT tolerance (20-50% of gradient)
  ├─ Score all peptides against all spectra:
  │     ├─ XCorr (BLAS-accelerated) → find best RT per peptide
  │     ├─ LibCosine (ppm matching) → score at best RT only
  │     └─ E-value from XCorr survival function
  ├─ Target-decoy competition → 1% FDR filtering
  ├─ Fit LOESS RT calibration from confident targets
  ├─ Calculate MS1/MS2 mass error statistics
  └─ Set tight RT tolerance for main search (3× residual SD)
  │
  ▼
PHASE 3: RIDGE REGRESSION (all files, per spectrum)
  ├─ Select candidates (isolation window + calibrated RT)
  ├─ Build design matrix from pre-binned library
  ├─ Solve NNLS ridge regression (f32, projected gradient)
  └─ Output: coefficients per candidate per spectrum
  │
  ▼
PHASE 4: POST-REGRESSION SCORING (per precursor)
  ├─ Aggregate coefficients into time series per precursor
  ├─ Peak detection → apex RT, boundaries
  ├─ Score at apex (all 30 features computed here):
  │     ├─ Spectral scores (dot_product, xcorr, hyperscore, etc.)
  │     ├─ Chromatographic features (peak shape, width, symmetry)
  │     └─ Contextual features (competitors, regression quality)
  ├─ Write PIN file (targets + decoys, 30 features)
  └─ Mokapot semi-supervised FDR → q-values
  │
  ▼
OUTPUT
  ├─ BiblioSpec (.blib) for Skyline
  ├─ Mokapot model weights (.pkl) for feature analysis
  └─ TSV report (optional)
```

---

## Feature Set (30 Features for Mokapot)

Osprey extracts 30 features per precursor for semi-supervised FDR control:

### Chromatographic Features (12)
| Feature | Description |
|---------|-------------|
| `peak_apex` | Peak apex coefficient (maximum deconvolution value) |
| `peak_area` | Integrated peak area (AUC of coefficients) |
| `emg_fit_quality` | EMG fit quality (R²) |
| `peak_width` | Peak width (FWHM in minutes) |
| `peak_symmetry` | Peak symmetry (leading/trailing ratio) |
| `rt_deviation` | RT deviation from prediction (minutes) |
| `rt_deviation_normalized` | Normalized RT deviation |
| `n_contributing_scans` | Number of scans contributing to peak |
| `coefficient_stability` | Coefficient variance near apex |
| `peak_sharpness` | Peak boundary sharpness |
| `peak_prominence` | Peak prominence (apex / baseline) |
| `modification_count` | Number of modifications |

### Spectral Features (13)
| Feature | Description |
|---------|-------------|
| `hyperscore` | X!Tandem-style: log(n_b!) + log(n_y!) + Σlog(I+1) |
| `xcorr` | Comet-style cross-correlation |
| `spectral_contrast_angle` | Normalized spectral contrast angle |
| `dot_product` | LibCosine with sqrt preprocessing |
| `dot_product_smz` | LibCosine with sqrt×mz² (SMZ) preprocessing |
| `pearson_correlation` | Pearson intensity correlation |
| `spearman_correlation` | Spearman rank correlation |
| `fragment_coverage` | Fraction of predicted fragments detected |
| `sequence_coverage` | Backbone coverage |
| `consecutive_ions` | Longest consecutive b/y ion run |
| `base_peak_rank` | Rank of base peak in predicted |
| `top3_matches` | Number of top-3 predicted fragments matched |
| `explained_intensity` | Fraction of observed intensity explained |

### Contextual Features (5)
| Feature | Description |
|---------|-------------|
| `n_competitors` | Number of competing candidates in regression |
| `relative_coefficient` | Coefficient relative to sum of all |
| `local_peptide_density` | Local peptide density (candidates per window) |
| `spectral_complexity` | Spectral complexity estimate |
| `regression_residual` | Regression residual (unexplained signal) |

---

## FDR Control

Osprey uses **two-level FDR control** via Mokapot:

### Run-Level FDR (Step 1)
- Per-file target-decoy competition
- All 30 features passed to Mokapot
- Mokapot trains linear SVM to combine features
- Q-values assigned per precursor

### Experiment-Level FDR (Step 2)
- Aggregate best precursor across replicates
- Only runs if multiple replicates provided
- Single replicate skips to run-level results
- Precursor = peptide + charge (not just peptide)

### Mokapot Integration
- PIN file format with 30 features
- `--save_models` flag saves feature weights
- Use `scripts/inspect_mokapot_weights.py` to view feature importance
- Parallel workers (auto-detected, capped at 8)
- Progress streaming to console

---

## Phase 2: Calibration Scoring (Quick Pass)

The calibration phase is a **fast scoring pass** designed to establish RT and mass calibration before the main search. It does NOT use ridge regression.

### How Calibration Scoring Works

Calibration uses a **hybrid XCorr + LibCosine** strategy for speed and accuracy:

1. **XCorr via BLAS** (fast, all spectra): Preprocess library and spectra into f32 vectors, score all pairs via matrix multiplication. XCorr is used to find the best-matching spectrum RT for each peptide. This is computationally cheap because it uses BLAS batched operations.

2. **LibCosine** (accurate, one spectrum): Once XCorr identifies the best RT, LibCosine is computed at that single spectrum using ppm-based fragment matching (no binning). This provides a high-quality score and mass error measurements.

3. **E-value** (significance): Calculated from the XCorr survival function (Comet-style) to estimate the probability of a match occurring by chance.

See: [Spectral Scoring](spectral-scoring.md), [XCorr Scoring](xcorr-scoring.md)

### What Calibration Produces

- **RT calibration curve** (LOESS): Maps library RT → expected measured RT
- **RT tolerance**: Residual SD × 3 (typically 0.5-2 min)
- **Mass calibration**: Mean/median PPM error for MS1 and MS2
- **Isolation scheme**: DIA window widths extracted from mzML

See: [RT Calibration](rt-calibration.md), [Mass Calibration](mass-calibration.md)

### Multi-File Strategy

- **First file**: Full calibration discovery (all peptides, wide tolerance)
- **Subsequent files**: Reuse calibration from first file (same LC column)
- Calibration saved to JSON for reuse on reruns

---

## Phase 3: Ridge Regression

Ridge regression deconvolutes each DIA spectrum into individual peptide contributions. Each spectrum contains fragments from multiple co-eluting precursors; regression separates these overlapping signals.

### Per-Spectrum Processing

For each MS2 spectrum:

1. **Candidate selection**: Find library entries whose precursor m/z falls within the isolation window AND whose calibrated RT is within tolerance
2. **Design matrix**: Columns are pre-binned library spectra (no preprocessing, raw intensity)
3. **NNLS ridge regression**: Minimize ‖Ax - b‖² + λ‖x‖² subject to x ≥ 0
4. **Output**: Non-negative coefficient per candidate (how much of the observed signal this peptide explains)

See: [Ridge Regression](ridge-regression.md)

---

## Phase 4: Post-Regression Scoring

After ridge regression assigns coefficients to all precursors across all spectra, Osprey aggregates these into coefficient time series per precursor and computes scoring features. **All scores are computed at the apex** -- the RT where the deconvoluted signal is strongest.

### Step 1: Coefficient Time Series

For each library entry, collect all (RT, coefficient) pairs from regression results across spectra. This forms a chromatographic profile of the deconvoluted peptide signal.

### Step 2: Peak Detection

The peak detector identifies elution peaks in the coefficient time series using a **threshold-crossing algorithm**:

1. Scan through the time series
2. When coefficient crosses above `min_height` (default: 0.05), mark peak start
3. Track the highest coefficient within the peak as the apex
4. When coefficient drops below `min_height`, mark peak end
5. Reject peaks narrower than `min_width` scans (default: 3)
6. If multiple peaks, select the one with the highest apex coefficient near the expected RT

See: [Peak Detection](peak-detection.md)

### Step 3: Score Computation at Apex

All spectral scores are computed at the apex spectrum (the observed MS2 spectrum closest to the coefficient apex RT). Each score captures a different aspect of the match:

| Score | Field | Description | Requires Library Intensities |
|-------|-------|-------------|------------------------------|
| LibCosine (sqrt) | `dot_product` | Cosine similarity, sqrt(intensity) preprocessing | Yes |
| LibCosine SMZ | `dot_product_smz` | Cosine similarity, sqrt(intensity)×mz² preprocessing | Yes |
| XCorr | `xcorr` | Comet-style cross-correlation | No |
| Hyperscore | `hyperscore` | X!Tandem: log(n_b!) + log(n_y!) + Σlog(I+1) | No |
| Delta RT | `rt_deviation` | Observed apex - predicted RT (minutes) | N/A |

Chromatographic features (peak shape, width, symmetry) and contextual features (regression quality, competitors) are also extracted.

See: [Deconvolution Scoring](deconvolution-scoring.md)

### Step 4: FDR Control via Mokapot

All 30 features for both targets and decoys are written to a PIN file. Mokapot performs semi-supervised learning (linear SVM) to combine features into an optimal discriminant score, then estimates q-values via target-decoy competition.

**Two-level FDR**:
1. **Run-level**: Per-file FDR control
2. **Experiment-level**: Best precursor per peptide+charge across all files

See: [FDR Control](fdr-control.md)

---

## Key Distinction: Calibration Scoring vs Post-Regression Scoring

| Aspect | Calibration (Phase 2) | Post-Regression (Phase 4) |
|--------|----------------------|---------------------------|
| **Purpose** | Establish RT/mass calibration | Determine peptide detections |
| **Uses regression?** | No | Yes (scores the deconvoluted result) |
| **Primary score** | LibCosine (with XCorr for RT) | All 30 features |
| **Fragment matching** | ppm-based (pyXcorrDIA-style) | ppm-based at apex spectrum |
| **FDR method** | Simple target-decoy (1% for calibration) | Mokapot semi-supervised ML |
| **Output** | RT curve + mass calibration | Peptide detections with q-values |
| **Scope** | First file only, wide tolerance | All files, tight calibrated tolerance |
| **Speed** | BLAS-accelerated batch scoring | Per-precursor feature extraction |

---

## Why Peptide-Centric Analysis?

Unlike spectrum-centric approaches (e.g., Comet, MSFragger) that identify each spectrum independently, Osprey:

1. **Considers all candidate peptides simultaneously** using ridge regression, handling overlapping signals in wide DIA windows
2. **Aggregates evidence across the chromatographic dimension** before scoring, using the coefficient time series rather than a single spectrum
3. **Scores at the best RT** for each peptide, where the deconvoluted signal is strongest

### Why Ridge Regression?

DIA spectra are mixed spectra containing fragments from multiple precursors. Ridge regression:
- Deconvolutes overlapping signals into individual peptide contributions
- Provides coefficients proportional to peptide abundance
- Handles collinearity (similar spectra) gracefully via L2 regularization
- Is computationally efficient with BLAS and f32 optimizations

### Why RT Calibration?

Library retention times may be predicted by deep learning models (Prosit, DeepLC), measured on different LC systems, or normalized to iRT scale. RT calibration converts library RTs to expected measured RTs, enabling tight RT filtering that reduces the candidate set and improves regression quality.

---

## Crate Architecture

```
osprey/
├── osprey-core/              # Core types, config, errors, traits
│   ├── config.rs             # YAML configuration structures
│   ├── isotope.rs            # Isotope distribution calculations
│   ├── types.rs              # LibraryEntry, Spectrum, FeatureSet, RegressionResult
│   └── traits.rs             # MS1SpectrumLookup, etc.
│
├── osprey-io/                # File I/O
│   ├── library/              # DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib
│   ├── mzml/                 # mzML parsing (via mzdata crate)
│   └── output/               # blib and report output
│
├── osprey-scoring/           # Spectral scoring and feature extraction
│   ├── lib.rs                # SpectralScorer, FeatureExtractor, DecoyGenerator
│   ├── batch.rs              # BLAS-accelerated batch scoring (calibration)
│   └── pipeline/             # Streaming pipeline components
│
├── osprey-chromatography/    # Chromatographic analysis
│   ├── lib.rs                # PeakDetector
│   └── calibration/          # LOESS RT calibration, mass calibration, I/O
│
├── osprey-regression/        # Ridge regression engine
│   ├── optimized.rs          # Pre-binned library, f32 solver, spectra cache
│   ├── ridge.rs              # Dense NNLS solver (f64)
│   ├── sparse.rs             # Sparse HRAM solver
│   └── binning.rs            # m/z binning
│
├── osprey-fdr/               # FDR control
│   ├── lib.rs                # Target-decoy competition, q-values
│   └── mokapot.rs            # PIN file output, Mokapot integration
│
└── osprey/                   # Main library + CLI binary
    ├── main.rs               # CLI entry point
    └── pipeline.rs           # Analysis pipeline orchestration
```

## Utility Scripts

```
scripts/
├── evaluate_calibration.py      # Generate HTML report from calibration JSON
├── inspect_mokapot_weights.py   # Display Mokapot feature weights/importance
└── ...
```

### evaluate_calibration.py

Generates an interactive HTML report visualizing calibration quality:
- RT calibration curve with residuals
- MS1/MS2 mass error histograms
- Candidate density heatmap (RT × m/z)

```bash
python scripts/evaluate_calibration.py calibration.json --output report.html
```

### inspect_mokapot_weights.py

Displays feature weights from trained Mokapot models to identify important features:

```bash
python scripts/inspect_mokapot_weights.py mokapot.model.pkl
```

---

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

## References

- pyXcorrDIA: https://github.com/maccoss/pyXcorrDIA
- Comet: https://comet-ms.sourceforge.io/
- Mokapot: https://github.com/wfondrie/mokapot
- LOESS: Cleveland, W.S. (1979). "Robust locally weighted regression and smoothing scatterplots"
