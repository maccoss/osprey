# Osprey Algorithm Documentation

This folder contains detailed documentation of the algorithms used in Osprey, a peptide-centric DIA analysis tool. Osprey uses a coelution search approach (DIA-NN-style fragment XIC correlation) to detect peptides in DIA data.

## Current Status (Working Prototype)

Osprey has a **working prototype** that can:
- Parse mzML files with DIA data
- Load spectral libraries (DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib)
- Generate enzyme-aware decoys
- Run auto-calibration with target-decoy FDR control
- Coelution search with 45 features extracted per precursor
- Run FDR control via native Percolator (default), external Mokapot, or simple TDC
- Output BiblioSpec .blib files for Skyline (library theoretical fragments)

## Algorithm Documentation

| Document | Description |
|----------|-------------|
| [Calibration](calibration.md) | RT, MS1, and MS2 calibration with LDA scoring and LOESS fitting |
| [Spectral Scoring](spectral-scoring.md) | XCorr + E-value scoring (calibration phase) |
| [XCorr Scoring](xcorr-scoring.md) | Comet-style XCorr implementation |
| [Peak Detection](peak-detection.md) | CWT consensus peak detection in fragment XIC time series |
| [Multi-Charge Consensus](composite-peak-selection.md) | Cross-charge-state peak boundary sharing |
| [Determinism](determinism.md) | Deterministic analysis: patterns, invariants, and maintenance |
| [Decoy Generation](decoy-generation.md) | Enzyme-aware sequence reversal for FDR control |
| [FDR Control](fdr-control.md) | Two-level FDR: native Percolator, Mokapot, and simple TDC |
| [BiblioSpec Output Schema](blib-output-schema.md) | blib output format and Skyline integration |

---

## Pipeline Overview

Osprey's pipeline has four major phases.

```
INPUT
  mzML files (DIA data) + Spectral Library (predicted spectra)
  |
  v
PHASE 1: INITIALIZATION
  +- Load spectral library (DIA-NN TSV, elib, or blib)
  +- Deduplicate by (modified_sequence, charge)
  +- Generate decoys (enzyme-aware reversal)
  +- Build m/z index for candidate lookup
  |
  v
PHASE 2: CALIBRATION DISCOVERY (first file only)
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
PHASE 3: COELUTION SEARCH
  +- Extract fragment XICs per precursor within RT window
  +- CWT consensus peak detection (Mexican Hat wavelet, median across transitions)
  +- Fallback: median polish profile or reference XIC peak detection
  +- Peak boundaries: zero-crossings ±2σ (valley guard)
  |
  v
PHASE 4: COELUTION SCORING (45 features)
  +- Pairwise coelution features (correlations between fragment XICs)
  +- Peak shape features (apex, area, width, symmetry, S/N)
  +- Spectral matching at apex (hyperscore, xcorr, dot products)
  +- Tukey median polish on fragment XICs
  +- Multi-charge consensus: all charge states of same peptide share peak boundaries
  +- FDR control (Percolator SVM or Mokapot) -> q-values
  |
  v
OUTPUT
  +- BiblioSpec (.blib) for Skyline (library theoretical fragments)
  +- SVM/Mokapot model weights for feature analysis
```

---

## Feature Set (45 Features)

Osprey extracts 45 features per precursor for semi-supervised FDR control via native Percolator (default) or external Mokapot. All intensity-based spectral similarity scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks.

### Pairwise Coelution (11)

| Feature | Description |
|---------|-------------|
| `fragment_coelution_sum` | Sum of all pairwise Pearson correlations between fragment XICs |
| `fragment_coelution_min` | Minimum pairwise correlation (worst pair) |
| `fragment_coelution_max` | Maximum pairwise correlation (best pair) |
| `n_coeluting_fragments` | Number of fragments with positive mean correlation |
| `n_fragment_pairs` | Number of fragment pairs in correlation matrix |
| `fragment_corr_0..5` | Per-fragment average correlation with all other fragments (top 6) |

### Peak Shape (7)

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

### Spectral at Apex (15)

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

### Mass Accuracy (3)

| Feature | Description |
|---------|-------------|
| `mass_accuracy_deviation_mean` | Mean signed mass error (ppm) across matched fragments |
| `abs_mass_accuracy_deviation_mean` | Mean absolute mass error (ppm) |
| `mass_accuracy_std` | Standard deviation of mass errors (ppm) |

### RT Deviation (2)

| Feature | Description |
|---------|-------------|
| `rt_deviation` | Observed apex RT - predicted RT from LOESS calibration (minutes) |
| `abs_rt_deviation` | Absolute value of RT deviation |

### MS1 Features (2)

HRAM only; 0.0 for unit resolution or when MS1 data is unavailable:

| Feature | Description |
|---------|-------------|
| `ms1_precursor_coelution` | Pearson correlation between reference XIC and MS1 monoisotopic precursor XIC |
| `ms1_isotope_cosine` | Cosine similarity between observed MS1 isotope envelope and theoretical distribution |

### Peptide Properties (2)

| Feature | Description |
|---------|-------------|
| `peptide_length` | Number of amino acids in peptide sequence |
| `missed_cleavages` | Number of missed enzymatic cleavages |

### Tukey Median Polish Features (3)

Fragment XIC matrix decomposition via Tukey median polish:

| Feature | Description |
|---------|-------------|
| `median_polish_cosine` | Cosine similarity between data-derived fragment intensities (row effects) and library intensities, sqrt preprocessing |
| `median_polish_rsquared` | R^2 of the additive model in sqrt space. Measures how well the shared elution profile + fragment intensities explain the XIC matrix. |
| `median_polish_residual_ratio` | Fraction of total signal unexplained by the model: sum|obs-pred| / sum(obs) in linear space. Lower = cleaner co-elution. |

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
- All 45 features used by SVM
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

The calibration phase is a **fast scoring pass** designed to establish RT and mass calibration before the main search. Calibration always uses **unit resolution XCorr** (2001 bins, 1.0005 Da) regardless of data type, providing ~50x faster scoring for HRAM data.

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

The coelution search extracts fragment XICs directly from the DIA data and scores precursors based on how well their fragment ion chromatograms co-elute.

For each precursor within the RT window:

1. **Fragment XIC extraction**: Extract chromatograms for each library fragment ion at the expected m/z (ppm tolerance)
2. **CWT consensus peak detection**: Convolve each fragment XIC with a Mexican Hat wavelet, compute pointwise median across transitions, find peaks in consensus signal. Boundaries via zero-crossings extended to ±2σ with valley guard. See [Peak Detection](peak-detection.md).
3. **Spectral scoring**: At the apex scan, score the observed spectrum against the library (hyperscore, xcorr, dot products)
4. **Tukey median polish**: Decompose fragment XIC matrix for robust scoring features and blib peak boundaries
5. **Multi-charge consensus**: After all windows are processed, all charge states of the same peptide are forced to share the same peak RT and integration boundaries. See [Multi-Charge Consensus](composite-peak-selection.md).

---

## Phase 4: Post-Search Scoring

After the main search, Osprey computes scoring features per precursor from fragment XIC profiles.

### Step 1: Peak Detection

CWT consensus peak detection uses Mexican Hat wavelet convolution of each fragment XIC, pointwise median across transitions, and ±2σ boundary extension with valley guard. Falls back to SG-smoothed peak detection if CWT finds no peaks.

See: [Peak Detection](peak-detection.md)

### Step 2: Tukey Median Polish and Peak Boundaries

The top 6 library fragments' XICs are decomposed via Tukey median polish. The column
effects give a robust shared elution profile (resistant to single-fragment interference).
FWHM on this profile determines peak integration boundaries at +/-1.96*sigma from apex.

### Step 3: Score Computation

All 45 features are computed within or at the peak boundaries:

- **Pairwise coelution**: Fragment XIC pairwise correlations
- **Peak shape**: Apex, area, width, symmetry, S/N, sharpness
- **Spectral at apex**: Hyperscore, xcorr, dot products, coverage, consecutive ions
- **Elution-weighted cosine**: LibCosine at each scan, weighted by reference XIC intensity^2
- **Mass accuracy**: Mean/std of mass errors across matched fragments
- **Tukey median polish features**: Cosine, R^2, residual ratio from the XIC decomposition
- **MS1 features**: Precursor co-elution and isotope cosine (HRAM only)
- **Peptide properties**: Peptide length, missed cleavages, RT deviation

### Step 4: FDR Control

All 45 features are combined into an optimal discriminant score via semi-supervised learning (linear SVM). The default native Percolator receives both targets and decoys and performs internal paired competition. For external Mokapot, targets and decoys are pre-competed using the best feature (selected by ROC AUC) and only winners are written to PIN files.

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
| **Peak detection** | CWT consensus (with fallback) | CWT consensus (with fallback) |
| **Method** | XCorr + E-value | Coelution search |
| **Features** | E-value only | 45 features |
| **FDR method** | Target-decoy on LDA score | Semi-supervised SVM (Percolator/Mokapot) |
| **Scope** | First file only, wide tolerance | All files, tight calibrated tolerance |

---

## Why Peptide-Centric Analysis?

Unlike spectrum-centric approaches (e.g., Comet, MSFragger) that identify each spectrum independently, Osprey:

1. **Aggregates evidence across the chromatographic dimension** before scoring, using fragment XIC profiles rather than a single spectrum
2. **Scores at the best RT** for each peptide, where the co-elution signal is strongest
3. **Uses pairwise fragment co-elution** to distinguish true signals from interference

### Why Coelution Search?

DIA spectra are mixed spectra containing fragments from multiple precursors. The coelution approach:
- Extracts fragment XICs and measures how well they co-elute (pairwise correlations)
- Fragments from the same peptide should rise and fall together
- Robust to interference because co-elution metrics are computed pairwise across fragments
- Scores spectral similarity at the detected apex for additional discrimination

### Why RT Calibration?

Library retention times may be predicted by deep learning models (Prosit, DeepLC), measured on different LC systems, or normalized to iRT scale. RT calibration converts library RTs to expected measured RTs, enabling tight RT filtering that reduces the candidate set and improves search quality.

---

## Crate Architecture

```
osprey/
+-- osprey-core/              # Core types, config, errors, traits
|   +-- config.rs             # YAML configuration structures
|   +-- isotope.rs            # Isotope distribution calculations
|   +-- types.rs              # LibraryEntry, Spectrum, CoelutionFeatureSet
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
|   +-- lib.rs                # PeakDetector, detect_xic_peak (legacy fallback)
|   +-- cwt.rs                # CWT consensus peak detection (Mexican Hat wavelet)
|   +-- calibration/          # LOESS RT calibration, mass calibration, I/O
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
