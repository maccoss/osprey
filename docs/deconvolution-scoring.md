# Deconvolution Scoring

After ridge regression deconvolutes each DIA spectrum into peptide contributions, Osprey computes scores for each detected peptide. These scores are features for Mokapot (semi-supervised FDR control).

This document covers the **post-deconvolution scoring pipeline**: how we select the best elution location for each peptide, and what scores are computed there.

## Pipeline Overview

```
Ridge Regression (per spectrum)
  │
  │  Output: coefficient per candidate peptide per spectrum
  ▼
Coefficient Time Series (per peptide, across RT)
  │
  │  For peptide P: [(RT₁, coef₁), (RT₂, coef₂), ..., (RTₙ, coefₙ)]
  ▼
Peak Detection
  │
  │  Find apex, boundaries, validate peak shape
  ▼
Apex Selection ← "Where is the best signal?"
  │
  │  Select the single best RT for this peptide
  ▼
Score Computation at Apex
  │
  │  Compute all features at the selected spectrum
  ▼
PIN Output (targets + decoys)
  │
  │  All features written to Mokapot PIN format
  ▼
FDR Control (Mokapot semi-supervised learning)
```

## Apex Selection: Defining "Best"

The apex is the RT where the deconvoluted signal for a peptide is strongest. This is where all spectral scores are computed.

### Current Strategy: Maximum Coefficient

```
apex_rt = argmax(coefficient[t])  over all t in coefficient series
```

The coefficient from NNLS ridge regression represents how much of the observed spectrum is explained by this peptide's library spectrum. The maximum coefficient is the point where the peptide contributes most to the observed signal.

### Why Maximum Coefficient?

1. **Physical meaning**: The coefficient represents the relative abundance of that peptide in the mixture at that time point. The maximum is the elution peak apex.

2. **Signal quality**: The spectrum at peak apex has the best signal-to-noise for that peptide, giving the most reliable spectral scores.

3. **Consistency with chromatographic peak**: The coefficient time series IS the chromatographic peak (after deconvolution), so its maximum is the natural apex.

### Alternative Strategies (Not Currently Implemented)

| Strategy | Pro | Con |
|----------|-----|-----|
| Max coefficient (current) | Simple, robust, physical meaning | May not be the best spectral match |
| Best XCorr across peak | Best spectral agreement | Computationally expensive (score every spectrum) |
| Best LibCosine across peak | Best intensity match | Same cost as XCorr approach |
| Coefficient-weighted average | Smooths noise | Loses the actual spectrum identity |
| Max explained variance | Best regression fit | Depends on all candidates, not just this one |

**Discussion**: The maximum coefficient approach selects based on **deconvolution quality** rather than spectral match quality. This is appropriate because:
- The deconvolution has already separated overlapping signals
- A high coefficient means the library spectrum explains a large fraction of the observed data
- Computing spectral scores at every time point would be expensive and redundant

If we wanted to refine apex selection, one approach would be to score the top-N coefficient spectra (e.g., apex and its two neighbors) and pick the one with the best spectral score. This would add minimal cost while potentially improving score quality.

---

## Spectral Scores

All spectral scores are computed at the apex spectrum (the observed MS2 spectrum closest in RT to the coefficient apex). Each score captures a different aspect of the match between the observed spectrum and the library prediction.

### 1. LibCosine (sqrt scaling)

**Requires library**: Yes (library fragment intensities)
**Implementation status**: Done
**FeatureSet field**: `dot_product`

LibCosine measures the cosine similarity between the observed and library spectra after sqrt intensity preprocessing.

#### Formula

```
For each matched fragment pair (obs_i, lib_i):
  obs_preprocessed_i = sqrt(obs_intensity_i)
  lib_preprocessed_i = sqrt(lib_intensity_i)

L2 normalize both vectors:
  obs_norm = obs_preprocessed / ||obs_preprocessed||₂
  lib_norm = lib_preprocessed / ||lib_preprocessed||₂

LibCosine = obs_norm · lib_norm  (dot product)
```

#### Score range: 0 to 1 (1 = perfect match)

#### Fragment matching
- PPM matching (10-20 ppm for HRAM, 0.5 Da for unit resolution)
- Closest m/z within tolerance (NOT highest intensity)
- Each library fragment matched to at most one observed peak

#### Why sqrt preprocessing?
- Down-weights dominant peaks
- Gives more influence to smaller diagnostic fragments
- Reduces impact of intensity measurement noise
- L2 normalization makes scores comparable regardless of total signal

#### Example
```
Library fragments:   y3=100  y4=80   y5=60   b3=40
Observed intensities: y3=9000 y4=7200 y5=4800  b3=?

sqrt preprocessing:
  Library:  [10.0,  8.94,  7.75,  6.32]
  Observed: [94.87, 84.85, 69.28, -   ]  (b3 not matched)

L2 normalize:
  Library:  [0.619, 0.554, 0.480, 0.392]  (only matched)
  → [0.619, 0.554, 0.480] / ||(0.619, 0.554, 0.480)||
  Observed: [0.672, 0.601, 0.491] / ||(...)||

LibCosine ≈ 0.997
```

### 2. LibCosine with sqrt(intensity) * mz² Scaling

**Requires library**: Yes (library fragment intensities)
**Implementation status**: Done
**FeatureSet field**: `dot_product_smz`

This variant weights fragment matches by their m/z value squared, giving more importance to higher m/z fragments which are more likely to be unique.

#### Formula

```
For each matched fragment pair (obs_i, lib_i) at m/z_i:
  obs_preprocessed_i = sqrt(obs_intensity_i) * mz_i²
  lib_preprocessed_i = sqrt(lib_intensity_i) * mz_i²

L2 normalize both vectors:
  obs_norm = obs_preprocessed / ||obs_preprocessed||₂
  lib_norm = lib_preprocessed / ||lib_preprocessed||₂

LibCosine_SMZ = obs_norm · lib_norm
```

#### Why m/z² weighting?
- Higher m/z fragments are more sequence-specific (longer ion series)
- Low m/z fragments (immonium ions, small b/y ions) are less diagnostic
- The m/z² weighting gives ~100× more weight to a fragment at m/z 1000 vs m/z 100
- Helps discriminate between peptides that share low-mass fragments
- This is the original SMZ (Sqrt-Mz-squared) preprocessing from spectral library searching

#### Example
```
Library fragments at m/z:  y3@348.2=100  y4@461.3=80  y5@576.3=60

sqrt * mz² preprocessing:
  Library:
    y3: sqrt(100) * 348.2² = 10.0 * 121,243 = 1,212,430
    y4: sqrt(80)  * 461.3² = 8.94 * 212,798 = 1,902,417
    y5: sqrt(60)  * 576.3² = 7.75 * 332,122 = 2,573,945

Note: y5 dominates despite lowest raw intensity, because
      it has the highest m/z → most diagnostic.
```

### 3. XCorr (Cross-Correlation)

**Requires library**: Partially (fragment m/z only, NOT intensities)
**Implementation status**: Done
**FeatureSet field**: `xcorr`

XCorr measures how well the observed spectrum matches the predicted fragment positions, independent of predicted intensities. It uses Comet-style preprocessing. Xcorr is unusual in that the preprocessing downweights peaks that don't match the predicted fragment ions used in the search. Most other scores ignore fragments that are not part of the predicted ion series.

#### Formula

```
Step 1: Bin experimental spectrum (Comet-style)
  bin_width = 1.0005079 Da
  BIN(mz) = (int)(mz / bin_width + 0.6)
  For each peak: obs_binned[BIN(mz)] += sqrt(intensity)

Step 2: Windowing normalization (10 windows)
  Divide m/z range into 10 equal windows
  For each window: normalize to max = 50.0
  Apply 5% global threshold

Step 3: Sliding window subtraction (offset = 75 bins)
  For each bin i:
    local_avg = mean(bins in [i-75, i+75], excluding i)
    xcorr_preprocessed[i] = windowed[i] - local_avg
  This removes baseline and enhances peaks above local background

Step 4: Bin theoretical spectrum (NO preprocessing!)
  For each library fragment: lib_binned[BIN(frag.mz)] = 1.0
  Unit intensity only - the library is just a "selector"

Step 5: Score
  XCorr = (Σ xcorr_preprocessed[i] * lib_binned[i]) × 0.005
```

#### Why XCorr doesn't need library intensities
- The theoretical spectrum uses unit intensity (1.0) for all fragments
- XCorr only asks: "is there signal at the predicted fragment positions?"
- This makes it robust to inaccurate library intensity predictions
- The preprocessing (windowing + flanking subtraction) removes noise and baseline

#### Score range: Unbounded, typically -0.5 to 2.0 for real matches

#### Key implementation details
- Experimental spectrum: sqrt + windowing + flanking subtraction
- Theoretical spectrum: NO preprocessing (just a selector)
- Scale factor 0.005 matches Comet's spectrum-centric scoring
- Preprocessing the theoretical spectrum would inflate scores ~50× (a common bug)

### 4. X!Tandem Hyperscore

**Requires library**: Partially (fragment m/z and ion types, NOT intensities)
**Implementation status**: Done
**FeatureSet field**: `hyperscore`

The X!Tandem hyperscore combines the number of matched fragment ions with their observed intensities, rewarding spectra with many matching fragments.

#### Formula

```
Step 1: Match predicted fragments to observed peaks
  For each predicted fragment f in library:
    Find closest observed peak within tolerance
    Record matched intensity I_f and ion type (b or y)

Step 2: Count matched ions by type
  n_b = count of matched b-ions
  n_y = count of matched y-ions

Step 3: Calculate score
  hyperscore = log(n_b!) + log(n_y!) + Σ log(I_f + 1)
  where:
    n_b! = factorial of number of matched b-ions
    n_y! = factorial of number of matched y-ions
    Σ log(I_f + 1) = sum over ALL matched fragments
```

#### Why hyperscore is useful
- Rewards having many matched ions of each type (through the factorial terms)
- log(n_b!) grows superlinearly: log(5!) = 4.79, log(10!) = 15.10
- A spectrum matching 10 b-ions and 8 y-ions scores much higher than one matching 18 b-ions and 0 y-ions
- The intensity term ensures matched peaks are real signal, not noise
- Does not depend on predicted intensities, only on fragment m/z positions

#### Score range: 0 to ~50 for typical tryptic peptides

#### Comparison to XCorr
- Hyperscore explicitly rewards having both b and y ions (factorial terms)
- XCorr treats all fragment positions equally
- Hyperscore uses raw observed intensities; XCorr uses preprocessed bins
- Both are library-intensity-independent

### 5. Delta RT (RT Deviation)

**Requires library**: Yes (expected RT from calibration)
**Implementation status**: Done
**FeatureSet field**: `rt_deviation` (minutes), `rt_deviation_normalized` (dimensionless)

The difference between the observed apex RT and the predicted RT from the LOESS calibration curve.

#### Formula

```
Step 1: Get predicted RT from calibration
  predicted_rt = LOESS_calibration.predict(library_entry.retention_time)

Step 2: Compute deviation
  rt_deviation = observed_apex_rt - predicted_rt    (in minutes)
  rt_deviation_normalized = rt_deviation / FWHM     (dimensionless)
```

#### Why it matters
- Peptides eluting far from their predicted RT are less likely to be correct
- The normalized version accounts for gradient steepness (wider peaks → more tolerance)
- Strong feature for Mokapot to discriminate targets from decoys

---

## Chromatographic Features (from Coefficient Time Series)

These features describe the shape and quality of the deconvoluted chromatographic peak.

| Feature | Field | Formula | Description |
|---------|-------|---------|-------------|
| Peak apex | `peak_apex` | max(coefficient[t]) | Maximum coefficient value |
| Peak area | `peak_area` | Σ coefficient[t] | Integrated AUC of coefficients |
| Peak width | `peak_width` | RT(right_half_max) - RT(left_half_max) | FWHM in minutes |
| Peak symmetry | `peak_symmetry` | leading_width / trailing_width | Ratio at half-max (1.0 = symmetric) |
| Contributing scans | `n_contributing_scans` | count(coefficient[t] > 0) | Non-zero coefficients |
| Coefficient stability | `coefficient_stability` | var(coefficients near apex) | Variance in top quartile |
| Peak sharpness | `peak_sharpness` | rise_rate at boundaries | Rate of coefficient change |
| Peak prominence | `peak_prominence` | apex / baseline_estimate | Signal above baseline |
| EMG fit quality | `emg_fit_quality` | estimated from symmetry | Proxy for exponentially modified Gaussian fit |

## Contextual Features (from Regression)

These features describe the regression environment at the time of deconvolution.

| Feature | Field | Formula | Description |
|---------|-------|---------|-------------|
| Competitors | `n_competitors` | max candidates per spectrum | How many peptides compete in the same spectrum |
| Relative coefficient | `relative_coefficient` | this_coef / Σ all_coefs | Fraction of total signal assigned to this peptide |
| Peptide density | `local_peptide_density` | avg competitors across spectra | Average crowding across the elution profile |
| Spectral complexity | `spectral_complexity` | estimated from regression | How complex the underlying spectra are |
| Regression residual | `regression_residual` | avg ‖Ax - b‖² | Average unexplained signal |

## Complete Scoring Flow in Code

```
score_run() in crates/osprey/src/pipeline.rs
│
├─ 1. Aggregate regression results by peptide ID
│     entry_data[lib_id] = [(RT₁, coef₁), (RT₂, coef₂), ...]
│
├─ 2. For each peptide (parallel):
│     ├─ Sort coefficient series by RT
│     ├─ Run peak detection (min_height = 0.05)
│     │     → Reject if no peaks found
│     │
│     ├─ Find apex_rt = RT at max coefficient
│     │
│     ├─ Collect peak region spectra:
│     │     Filter: isolation_window contains precursor_mz
│     │             AND RT within coefficient series range
│     │
│     ├─ Build RegressionContext:
│     │     avg_n_competitors, avg_relative_coefficient,
│     │     avg_residual, avg_explained_variance
│     │
│     ├─ Extract chromatographic features from coefficient series
│     │     peak_apex, peak_area, peak_width, peak_symmetry, ...
│     │
│     ├─ Aggregate spectra around apex (weighted by coefficients)
│     │     → Creates a single "deconvoluted" spectrum
│     │
│     ├─ Find apex spectrum (closest to apex_rt)
│     │
│     ├─ Compute spectral scores at apex spectrum:
│     │     ├─ LibCosine (sqrt preprocessing) ← primary score
│     │     ├─ XCorr (Comet-style)
│     │     ├─ Fragment coverage, explained intensity
│     │     ├─ Pearson/Spearman correlations
│     │     └─ Consecutive ions, sequence coverage, etc.
│     │
│     ├─ Apply contextual features from regression
│     │
│     └─ Return ScoredEntry with all features
│
└─ 3. Output: Vec<ScoredEntry> (targets + decoys)
       → Written to PIN file for Mokapot
```

## Score Summary Table

| Score | Field | Requires Library Intensities | Range | Status |
|-------|-------|------------------------------|-------|--------|
| LibCosine (sqrt) | `dot_product` | Yes | 0-1 | Done |
| LibCosine (sqrt*mz²) | `dot_product_smz` | Yes | 0-1 | Done |
| XCorr | `xcorr` | No (unit intensity) | -0.5 to 2.0 | Done |
| X!Tandem Hyperscore | `hyperscore` | No (fragment types only) | 0-50 | Done |
| Delta RT | `rt_deviation` | Yes (calibrated RT) | minutes | Done |
| Delta RT (normalized) | `rt_deviation_normalized` | Yes | dimensionless | Done |
| Fragment coverage | `fragment_coverage` | Yes | 0-1 | Done |
| Explained intensity | `explained_intensity` | Yes (fragment m/z) | 0-1 | Done |
| Pearson correlation | `pearson_correlation` | Yes | -1 to 1 | Done |
| Spearman correlation | `spearman_correlation` | Yes | -1 to 1 | Done |
| Consecutive ions | `consecutive_ions` | Yes (ion types) | 0-N | Done |

## Implementation Files

| File | Purpose |
|------|---------|
| `crates/osprey-scoring/src/lib.rs` | SpectralScorer (LibCosine, LibCosine SMZ, XCorr, Hyperscore), FeatureExtractor |
| `crates/osprey-scoring/src/batch.rs` | BLAS-accelerated batch scoring (calibration phase) |
| `crates/osprey/src/pipeline.rs` | score_run() - orchestrates scoring after regression |
| `crates/osprey-core/src/types.rs` | FeatureSet struct (all 30 features) |
| `crates/osprey-fdr/src/mokapot.rs` | PIN file output for Mokapot |

## TODO

1. **Consider apex refinement**: Score top-N coefficient spectra and pick best spectral match
2. **EMG peak fitting**: Replace symmetry-based heuristic with true EMG fit
