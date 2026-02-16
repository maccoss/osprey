# Osprey Algorithm Documentation

This folder contains detailed documentation of the algorithms used in Osprey, a peptide-centric DIA analysis tool. Osprey supports two search modes: ridge regression (spectrum deconvolution) and coelution (DIA-NN-style fragment XIC correlation).

## Current Status (Working Prototype)

Osprey has a **working prototype** that can:
- Parse mzML files with DIA data
- Load spectral libraries (DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib)
- Generate enzyme-aware decoys
- Run auto-calibration with target-decoy FDR control
- Search via ridge regression (37 features) or coelution (45 features)
- Run FDR control via native Percolator (default), external Mokapot, or simple TDC
- Output BiblioSpec .blib files for Skyline (library theoretical fragments)

## Algorithm Documentation

| Document | Description |
|----------|-------------|
| [Calibration](calibration.md) | RT, MS1, and MS2 calibration with LDA scoring and LOESS fitting |
| [Spectral Scoring](spectral-scoring.md) | XCorr + E-value scoring (calibration phase) |
| [XCorr Scoring](xcorr-scoring.md) | Comet-style XCorr implementation |
| [Deconvolution Scoring](deconvolution-scoring.md) | Post-regression scoring: all scores computed after ridge regression |
| [Ridge Regression](ridge-regression.md) | NNLS ridge regression for spectrum deconvolution |
| [Peak Detection](peak-detection.md) | Chromatographic peak detection in coefficient time series |
| [Decoy Generation](decoy-generation.md) | Enzyme-aware sequence reversal for FDR control |
| [FDR Control](fdr-control.md) | Two-level FDR: native Percolator, Mokapot, and simple TDC |
| [BiblioSpec Output Schema](blib-output-schema.md) | blib output format and Skyline integration |

---

## Pipeline Overview

Osprey's pipeline has four major phases. Phases 1 and 2 (initialization and calibration) are shared between search modes. Phase 3 diverges: regression mode deconvolutes mixed spectra via ridge regression, while coelution mode extracts and correlates fragment XICs directly. Phase 4 (scoring and FDR) follows the same structure in both modes but extracts different features.

```
INPUT
  mzML files (DIA data) + Spectral Library (predicted spectra)
  |
  v
PHASE 1: INITIALIZATION (shared)
  +- Load spectral library (DIA-NN TSV, elib, or blib)
  +- Deduplicate by (modified_sequence, charge)
  +- Generate decoys (enzyme-aware reversal)
  +- Build m/z index for candidate lookup
  |
  v
PHASE 2: CALIBRATION DISCOVERY (first file only, shared)
  +- Calculate library-to-measured RT mapping
  +- Set wide RT tolerance (20-50% of gradient)
  +- Score all peptides against all spectra:
  |     +- XCorr (unit resolution bins, BLAS sdot) -> find best RT per peptide
  |     +- Top-3 fragment matching (binary search) -> MS2 mass errors
  |     +- E-value from XCorr survival function
  +- Target-decoy competition -> 1% FDR filtering
  +- Fit LOESS RT calibration from confident targets
  +- Calculate MS1/MS2 mass error statistics
  +- Set tight RT tolerance for main search (3x residual SD)
  |
  v
PHASE 3A: RIDGE REGRESSION (--search-mode Regression)
  +- Select candidates (isolation window + calibrated RT)
  +- Build design matrix from pre-binned library
  +- Solve NNLS ridge regression (f32, projected gradient)
  +- Output: coefficients per candidate per spectrum
  |
  v
PHASE 4A: POST-REGRESSION SCORING (37 features)
  +- Aggregate coefficients into time series per precursor
  +- Peak detection -> apex RT, boundaries
  +- Tukey median polish on fragment XICs -> elution profile, FWHM
  +- Score at apex (37 features)
  +- FDR control (Percolator SVM or Mokapot) -> q-values

--- OR ---

PHASE 3B: COELUTION SEARCH (--search-mode Coelution)
  +- Extract fragment XICs per precursor within RT window
  +- Compute pairwise correlations between fragment XICs
  +- Build reference XIC (best-correlated fragment)
  +- Peak detection on reference XIC -> apex, boundaries
  |
  v
PHASE 4B: COELUTION SCORING (45 features)
  +- Pairwise coelution features (correlations between fragment XICs)
  +- Peak shape features (apex, area, width, symmetry, S/N)
  +- Spectral matching at apex (hyperscore, xcorr, dot products)
  +- Tukey median polish on fragment XICs
  +- FDR control (Percolator SVM or Mokapot) -> q-values
  |
  v
OUTPUT (shared)
  +- BiblioSpec (.blib) for Skyline (library theoretical fragments)
  +- SVM/Mokapot model weights for feature analysis
```

---

## Feature Sets

Osprey extracts features per precursor for semi-supervised FDR control via native Percolator (default) or external Mokapot. The feature set depends on the search mode: **37 features** for regression, **45 features** for coelution. All intensity-based spectral similarity scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks.

### Regression Mode (37 Features)

### Ridge Regression Features (8)

Derived from coefficient time series and chromatographic analysis:

| Feature | Description |
|---------|-------------|
| `peak_apex` | Peak apex coefficient (maximum deconvolution value) |
| `peak_area` | Integrated peak area (AUC of coefficients) |
| `peak_width` | Peak width (FWHM in minutes), from Tukey median polish elution profile |
| `coefficient_stability` | Coefficient variance near apex |
| `relative_coefficient` | Coefficient relative to sum of all candidates |
| `explained_intensity` | Fraction of observed intensity explained by matched fragments |
| `signal_to_noise` | Peak apex / noise estimate from coefficient series |
| `xic_signal_to_noise` | Signal-to-noise from the best-correlated fragment XIC |

### Spectral Matching Features - Mixed (2)

Computed on the observed (mixed) spectrum at apex:

| Feature | Description |
|---------|-------------|
| `xcorr` | Comet-style cross-correlation (unit intensity library selector) |
| `consecutive_ions` | Longest consecutive b/y ion run |

### Spectral Matching Features - Deconvoluted (2)

Computed on coefficient-weighted spectra (apex +/- 2 scans):

| Feature | Description |
|---------|-------------|
| `xcorr_deconv` | XCorr on deconvoluted spectrum |
| `consecutive_ions_deconv` | Consecutive ions on deconvoluted spectrum |

### RT Deviation (1)

| Feature | Description |
|---------|-------------|
| `rt_deviation` | Observed apex RT - predicted RT from LOESS calibration (minutes) |

### Fragment Co-elution Features (9)

Computed within peak integration boundaries. Each fragment's XIC is correlated
against the coefficient time series:

| Feature | Description |
|---------|-------------|
| `fragment_coelution_sum` | Sum of Pearson correlations across top 6 fragments |
| `fragment_coelution_min` | Minimum correlation among scored fragments |
| `n_coeluting_fragments` | Count of fragments with positive correlation |
| `fragment_corr_0` | Correlation of highest-intensity library fragment vs coefficient series |
| `fragment_corr_1` | Correlation of 2nd highest fragment |
| `fragment_corr_2` | Correlation of 3rd highest fragment |
| `fragment_corr_3` | Correlation of 4th highest fragment |
| `fragment_corr_4` | Correlation of 5th highest fragment |
| `fragment_corr_5` | Correlation of 6th highest fragment |

### Elution-Weighted Spectral Similarity (1)

| Feature | Description |
|---------|-------------|
| `elution_weighted_cosine` | LibCosine computed at each scan within peak boundaries, weighted by coefficient^2, averaged. Uses sqrt preprocessing and includes all library fragments (0 for unmatched). |

### Mass Accuracy Features (3)

| Feature | Description |
|---------|-------------|
| `mass_accuracy_deviation_mean` | Mean signed mass error (ppm) across matched fragments |
| `abs_mass_accuracy_deviation_mean` | Mean absolute mass error (ppm) |
| `mass_accuracy_std` | Standard deviation of mass errors (ppm) |

### Percolator-style Features (6)

| Feature | Description |
|---------|-------------|
| `abs_rt_deviation` | Absolute value of RT deviation |
| `peptide_length` | Number of amino acids in peptide sequence |
| `missed_cleavages` | Number of missed enzymatic cleavages |
| `ln_num_candidates` | ln(number of candidates in the regression window) |
| `coef_zscore` | Z-score of this precursor's coefficient at apex |
| `coef_zscore_mean` | Mean z-score of coefficients across the peak |

### MS1 Features (2)

HRAM only; 0.0 for unit resolution or when MS1 data is unavailable:

| Feature | Description |
|---------|-------------|
| `ms1_precursor_coelution` | Pearson correlation between coefficient series and MS1 monoisotopic precursor XIC |
| `ms1_isotope_cosine` | Cosine similarity between observed MS1 isotope envelope and theoretical distribution |

### Tukey Median Polish Features (3)

Fragment XIC matrix decomposition via Tukey median polish:

| Feature | Description |
|---------|-------------|
| `median_polish_cosine` | Cosine similarity between data-derived fragment intensities (row effects) and library intensities, sqrt preprocessing |
| `median_polish_rsquared` | R^2 of the additive model in sqrt space. Measures how well the shared elution profile + fragment intensities explain the XIC matrix. |
| `median_polish_residual_ratio` | Fraction of total signal unexplained by the model: sum|obs-pred| / sum(obs) in linear space. Lower = cleaner co-elution. |

### Coelution Mode (45 Features)

The coelution search mode extracts fragment XICs directly (no regression) and computes pairwise correlations between them. It uses a different feature set optimized for XIC-based scoring.

#### Pairwise Coelution (11)

| Feature | Description |
|---------|-------------|
| `fragment_coelution_sum` | Sum of all pairwise Pearson correlations between fragment XICs |
| `fragment_coelution_min` | Minimum pairwise correlation (worst pair) |
| `fragment_coelution_max` | Maximum pairwise correlation (best pair) |
| `n_coeluting_fragments` | Number of fragments with positive mean correlation |
| `n_fragment_pairs` | Number of fragment pairs in correlation matrix |
| `fragment_corr_0..5` | Per-fragment average correlation with all other fragments (top 6) |

#### Peak Shape (7)

Derived from the reference XIC (best-correlated fragment):

| Feature | Description |
|---------|-------------|
| `peak_apex` | Peak apex intensity (background-subtracted) |
| `peak_area` | Integrated peak area within boundaries |
| `peak_width` | Peak width (FWHM in minutes) |
| `peak_symmetry` | Leading/trailing area ratio around apex |
| `signal_to_noise` | Signal-to-noise ratio |
| `n_scans` | Number of scans within peak boundaries |
| `peak_sharpness` | Steepness of peak edges |

#### Spectral at Apex (15)

| Feature | Description |
|---------|-------------|
| `hyperscore` | X!Tandem-style hyperscore |
| `xcorr` | Comet-style cross-correlation |
| `dot_product` | Library cosine (sqrt intensity preprocessing) |
| `dot_product_smz` | Library cosine with sqrt(intensity) * mz^2 preprocessing |
| `dot_product_top6/5/4` | Library cosine using only top N fragments |
| `dot_product_smz_top6/5/4` | SMZ cosine using only top N fragments |
| `fragment_coverage` | Fraction of library fragments matched |
| `sequence_coverage` | Fraction of peptide backbone covered by b/y ions |
| `consecutive_ions` | Longest consecutive b or y ion series |
| `explained_intensity` | Fraction of observed intensity explained by matches |
| `elution_weighted_cosine` | LibCosine at each scan, weighted by reference XIC intensity^2, averaged |

#### Mass Accuracy (3), RT Deviation (2), MS1 (2), Peptide Properties (2)

Same as regression mode, except coelution includes both `rt_deviation` and `abs_rt_deviation` (2 features) and `peptide_length` + `missed_cleavages` (no ln_num_candidates, coef_zscore, coef_zscore_mean).

#### Tukey Median Polish (3)

Same as regression mode: `median_polish_cosine`, `median_polish_rsquared`, `median_polish_residual_ratio`.

---

## Tukey Median Polish

The Tukey median polish decomposes the fragment XIC matrix into an additive model:

```
ln(Observed[f,s]) = mu + alpha_f + beta_s + epsilon_fs

  f = fragment index (1..6, top 6 by library intensity)
  s = scan index (1..N)
  mu = overall effect (grand median)
  alpha_f = row effect = data-derived fragment relative intensity
  beta_s = column effect = shared elution profile shape
  epsilon_fs = residual = interference + noise
```

The algorithm iterates (max 20 iterations, convergence tolerance 1e-4):
1. **Row sweep**: subtract nanmedian of each row from residuals
2. **Column sweep**: subtract nanmedian of each column from residuals
3. Check convergence: max(|change|) < tolerance

### Uses

- **Peak boundaries**: Column effects give a robust elution profile by borrowing
  strength across all 6 transitions. The median suppresses interference from any
  single fragment. FWHM computed on this profile determines peak boundaries:
  sigma = FWHM / 2.355, boundaries at +/-1.96*sigma from apex.

- **Scoring features**: Row effects give interference-free fragment intensities
  that are scored against the library (median_polish_cosine). R^2 and residual
  ratio measure how well the additive model explains the data.

---

## FDR Control

Osprey uses **two-level FDR control** with three available methods (`--fdr-method`):

- **Percolator** (default): Native linear SVM with cross-validation. Both targets and decoys enter; Percolator does internal paired competition within each fold.
- **Mokapot**: External Python tool. Requires pre-competed PIN files (only competition winners). Uses ROC AUC to select the best feature for upstream competition.
- **Simple**: Direct target-decoy competition on the best single feature.

### Run-Level FDR
- Per-file target-decoy competition and q-value computation
- All features (37 for regression, 45 for coelution) used by SVM
- Q-values assigned at both precursor and peptide level
- Effective q-value: `max(precursor_qvalue, peptide_qvalue)` — must pass both

### Experiment-Level FDR
- Takes the **single best-scoring observation** per precursor (modified_sequence + charge) across all files
- Target-decoy competition on this deduplicated set at both precursor and peptide level
- Effective q-value: `max(precursor_qvalue, peptide_qvalue)`
- Single replicate skips experiment-level (uses run-level directly)

### Multi-File Observation Propagation
- After experiment-level FDR determines passing precursors, **all per-file target observations** for those precursors are included in the blib output
- Each file gets its own RT boundaries; best experiment_qvalue is propagated to all observations
- This enables Skyline to use per-file peak boundaries for quantification across replicates/GPF files

See [FDR Control](fdr-control.md) for full algorithm details, fold assignment, and the target-decoy competition strategy.

---

## Phase 2: Calibration Scoring (Quick Pass)

The calibration phase is a **fast scoring pass** designed to establish RT and mass calibration before the main search. It does NOT use ridge regression. Calibration always uses **unit resolution XCorr** (2001 bins, 1.0005 Da) regardless of data type, providing ~50x faster scoring for HRAM data.

### How Calibration Scoring Works

1. **Pre-filter** (top-3 fragment match): Binary search check -- at least 1 of the top-3 most intense library fragments must be present in the observed spectrum. Eliminates candidates with no spectral evidence before expensive XCorr preprocessing. O(3 x log n) per candidate.

2. **XCorr** (all passing spectra): Preprocess experimental spectrum with sqrt, windowing normalization, and flanking bin subtraction. Score against all spectra in RT window via BLAS sdot. Used to find the best-matching RT for each peptide.

3. **Top-3 fragment matching** (at best XCorr spectrum): Binary search the top-3 most intense library fragments in the best XCorr spectrum. Returns signed mass errors (ppm) for MS2 mass calibration.

4. **E-value** (significance): Calculated from the XCorr survival function (Comet-style) -- fit log-linear regression to survival counts. Used for target-decoy competition.

5. **MS1 isotope envelope** (at best XCorr spectrum): Extract monoisotopic peak from nearest MS1 spectrum. Returns MS1 mass error (ppm) for precursor mass calibration.

See: [Spectral Scoring](spectral-scoring.md), [XCorr Scoring](xcorr-scoring.md)

### What Calibration Produces

- **RT calibration curve** (LOESS): Maps library RT -> expected measured RT
- **RT tolerance**: Residual SD x 3 (typically 0.5-2 min)
- **Mass calibration**: Mean/median PPM error for MS1 and MS2
- **Isolation scheme**: DIA window widths extracted from mzML

See: [Calibration](calibration.md)

### Multi-File Strategy

- **First file**: Full calibration discovery (all peptides, wide tolerance)
- **Subsequent files**: Reuse calibration from first file (same LC column)
- Calibration saved to JSON for reuse on reruns

---

## Phase 3: Main Search

Osprey supports two search strategies. Both produce scored entries that feed into FDR control (native Percolator or Mokapot).

### Regression Mode (default)

Ridge regression deconvolutes each DIA spectrum into individual peptide contributions. Each spectrum contains fragments from multiple co-eluting precursors; regression separates these overlapping signals.

### Per-Spectrum Processing

For each MS2 spectrum:

1. **Candidate selection**: Find library entries whose precursor m/z falls within the isolation window AND whose calibrated RT is within tolerance
2. **Design matrix**: Columns are pre-binned library spectra (no preprocessing, raw intensity)
3. **NNLS ridge regression**: Minimize ||Ax - b||^2 + lambda*||x||^2 subject to x >= 0
4. **Output**: Non-negative coefficient per candidate (how much of the observed signal this peptide explains)

See: [Ridge Regression](ridge-regression.md)

### Coelution Mode

Instead of regression, the coelution search extracts fragment XICs directly from the DIA data and scores precursors based on how well their fragment ion chromatograms co-elute.

For each precursor within the RT window:

1. **Fragment XIC extraction**: Extract chromatograms for each library fragment ion at the expected m/z (ppm tolerance)
2. **Pairwise correlation**: Compute Pearson correlations between all pairs of fragment XICs. Fragments from the same peptide should co-elute
3. **Reference XIC**: Select the fragment with highest mean correlation as the reference chromatogram
4. **Peak detection**: Find the elution peak in the reference XIC
5. **Spectral scoring**: At the apex scan, score the observed spectrum against the library (hyperscore, xcorr, dot products)
6. **Tukey median polish**: Decompose fragment XIC matrix for robust peak boundaries and scoring features

---

## Phase 4: Post-Search Scoring

After the main search, Osprey aggregates results per precursor and computes scoring features. For regression mode, this means coefficient time series; for coelution mode, fragment XIC profiles.

### Step 1: Coefficient Time Series

For each library entry, collect all (RT, coefficient) pairs from regression results across spectra. This forms a chromatographic profile of the deconvoluted peptide signal.

### Step 2: Peak Detection

Two complementary peak detectors are used:

1. **Threshold-crossing** (`PeakDetector::detect()`): Scans coefficient series for contiguous regions above `min_height` (0.01). Selects the peak with highest apex coefficient near the expected RT.
2. **XIC peak detection** (`detect_xic_peak()`): Smoothed apex finding, DIA-NN-style valley-based boundary detection, and asymmetric FWHM capping. Used for signal-to-noise estimation on coefficient series.

See: [Peak Detection](peak-detection.md)

### Step 3: Tukey Median Polish and Peak Boundaries

The top 6 library fragments' XICs are decomposed via Tukey median polish. The column
effects give a robust shared elution profile (resistant to single-fragment interference).
FWHM on this profile determines peak integration boundaries at +/-1.96*sigma from apex.

### Step 4: Score Computation

All features are computed within or at the peak boundaries:

- **Fragment co-elution**: Per-fragment XIC correlations against the coefficient series
- **Elution-weighted cosine**: LibCosine at each scan, weighted by coefficient^2
- **Mass accuracy**: Mean/std of mass errors across matched fragments
- **XCorr**: Comet-style cross-correlation at the apex spectrum
- **Tukey median polish features**: Cosine, R^2, residual ratio from the XIC decomposition
- **MS1 features**: Precursor co-elution and isotope cosine (HRAM only)
- **Percolator-style**: Peptide length, missed cleavages, charge, RT deviation

### Step 5: FDR Control

All features (37 for regression, 45 for coelution) are combined into an optimal discriminant score via semi-supervised learning (linear SVM). The default native Percolator receives both targets and decoys and performs internal paired competition. For external Mokapot, targets and decoys are pre-competed using the best feature (selected by ROC AUC) and only winners are written to PIN files.

**Two-level FDR with dual precursor+peptide control**:
1. **Run-level**: Per-file FDR at both precursor and peptide level
2. **Experiment-level**: Best observation per precursor across all files, FDR at both levels
3. **Effective q-value**: `max(precursor_qvalue, peptide_qvalue)` at each level
4. **Observation propagation**: All per-file observations for passing precursors included in output

See: [FDR Control](fdr-control.md)

---

## Key Distinctions

### Calibration vs Main Search

| Aspect | Calibration (Phase 2) | Main Search (Phase 3-4) |
|--------|----------------------|---------------------------|
| **Purpose** | Establish RT/mass calibration | Determine peptide detections |
| **Method** | XCorr + E-value | Regression or coelution |
| **Features** | E-value only | 37 (regression) or 45 (coelution) |
| **FDR method** | Target-decoy on LDA score | Semi-supervised SVM (Percolator/Mokapot) |
| **Scope** | First file only, wide tolerance | All files, tight calibrated tolerance |

### Regression vs Coelution Search

| Aspect | Regression | Coelution |
|--------|-----------|-----------|
| **Core method** | Ridge regression deconvolution | Fragment XIC correlation |
| **Signal source** | Coefficient time series | Fragment XICs directly |
| **Features** | 37 (regression-specific) | 45 (XIC-specific) |
| **Unique features** | coefficient_stability, relative_coefficient, xcorr_deconv, coef_zscore | hyperscore, dot_product variants, peak_symmetry, peak_sharpness |
| **Handles interference** | Via regression (simultaneous fit) | Via pairwise correlation filtering |

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
+-- osprey-core/              # Core types, config, errors, traits
|   +-- config.rs             # YAML configuration structures
|   +-- isotope.rs            # Isotope distribution calculations
|   +-- types.rs              # LibraryEntry, Spectrum, FeatureSet, CoelutionFeatureSet, SearchMode
|   +-- traits.rs             # MS1SpectrumLookup, etc.
|
+-- osprey-io/                # File I/O
|   +-- library/              # DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib
|   +-- mzml/                 # mzML parsing (via mzdata crate)
|   +-- output/               # blib and report output
|
+-- osprey-scoring/           # Spectral scoring and feature extraction
|   +-- lib.rs                # SpectralScorer, FeatureExtractor, DecoyGenerator,
|   |                         # Tukey median polish, fragment co-elution
|   +-- batch.rs              # BLAS-accelerated batch scoring (calibration)
|   +-- calibration_ml.rs     # LDA calibration scoring
|   +-- pipeline/             # Streaming pipeline components
|
+-- osprey-chromatography/    # Chromatographic analysis
|   +-- lib.rs                # PeakDetector
|   +-- calibration/          # LOESS RT calibration, mass calibration, I/O
|
+-- osprey-regression/        # Ridge regression engine
|   +-- solver.rs             # Unified f32 CD-NNLS solver (unit res + HRAM)
|   +-- cd_nnls.rs            # Coordinate Descent NNLS algorithm
|   +-- ridge.rs              # Dense NNLS solver (f64, reference)
|   +-- binning.rs            # m/z binning (Comet BIN macro)
|   +-- matrix.rs             # Design matrix construction
|
+-- osprey-fdr/               # FDR control
|   +-- lib.rs                # Target-decoy competition, q-values
|   +-- mokapot.rs            # PIN file output, Mokapot integration
|
+-- osprey/                   # Main library + CLI binary
    +-- main.rs               # CLI entry point
    +-- pipeline.rs           # Analysis pipeline orchestration
```

## Utility Scripts

```
scripts/
+-- evaluate_calibration.py      # Generate HTML report from calibration JSON
+-- inspect_mokapot_weights.py   # Display Mokapot feature weights/importance
```

### evaluate_calibration.py

Generates an interactive HTML report visualizing calibration quality:
- RT calibration curve with residuals
- MS1/MS2 mass error histograms
- Candidate density heatmap (RT x m/z)

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

# Search mode: Regression (default) or Coelution
search_mode: Regression

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
