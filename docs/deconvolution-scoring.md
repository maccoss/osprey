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
**FeatureSet field**: `rt_deviation` (minutes)

The difference between the observed apex RT and the predicted RT from the LOESS calibration curve.

#### Formula

```
Step 1: Get predicted RT from calibration
  predicted_rt = LOESS_calibration.predict(library_entry.retention_time)

Step 2: Compute deviation
  rt_deviation = observed_apex_rt - predicted_rt    (in minutes)
```

#### Why it matters
- Peptides eluting far from their predicted RT are less likely to be correct
- Strong feature for Mokapot to discriminate targets from decoys

### 6. Fragment Coverage

**Requires library**: Yes (predicted fragment m/z positions)
**Implementation status**: Done
**FeatureSet field**: `fragment_coverage`

The fraction of predicted library fragments that were matched in the observed spectrum.

#### Formula

```
fragment_coverage = n_matched_fragments / n_predicted_fragments
```

Where:
- `n_matched_fragments` = count of library fragments with a matching observed peak (within tolerance)
- `n_predicted_fragments` = total fragments in the library spectrum

#### Score range: 0 to 1 (1 = all predicted fragments found)

### 7. Sequence Coverage

**Requires library**: Yes (ion series positions)
**Implementation status**: Done
**FeatureSet field**: `sequence_coverage`

The fraction of backbone cleavage sites covered by at least one matched ion (b or y).

#### Formula

```
For a peptide of length L:
  n_cleavage_sites = L - 1  (positions between each amino acid)

For each site i (1 to L-1):
  covered[i] = True if b_i OR y_{L-i} is matched

sequence_coverage = sum(covered) / n_cleavage_sites
```

#### Example
```
Peptide: PEPTIDE (7 AA, 6 cleavage sites)

Matched ions: b2, b3, y4, y5
  Site 1: b1? NO, y6? NO  → not covered
  Site 2: b2? YES         → covered
  Site 3: b3? YES         → covered
  Site 4: b4? NO, y3? NO  → not covered
  Site 5: b5? NO, y2? NO  → not covered
  Site 6: b6? NO, y1? NO  → not covered

Wait, let's recalculate with y4, y5:
  y4 covers site 7-4=3, y5 covers site 7-5=2

sequence_coverage = 2/6 = 0.33
```

### 8. Top-3 Matches

**Requires library**: Yes (fragment intensities to identify top-3)
**Implementation status**: Done
**FeatureSet field**: `top3_matches`

The number of the 3 most intense predicted fragments that were matched.

#### Formula

```
Step 1: Rank library fragments by predicted intensity (descending)
Step 2: Take top 3 most intense fragments
Step 3: Count how many of these are matched in observed spectrum

top3_matches ∈ {0, 1, 2, 3}
```

#### Why it matters
- The most intense predicted fragments are most likely to be visible
- Missing all top-3 fragments suggests a poor match or wrong peptide
- Simple, interpretable feature

### 9. Explained Intensity

**Requires library**: Yes (fragment m/z positions)
**Implementation status**: Done
**FeatureSet field**: `explained_intensity`

The fraction of observed intensity that is explained by matched library fragments.

#### Formula

```
matched_intensity = Σ intensity of observed peaks that match a library fragment
total_intensity = Σ intensity of all observed peaks (above threshold)

explained_intensity = matched_intensity / total_intensity
```

#### Score range: 0 to 1 (1 = all signal comes from predicted fragments)

---

## Deconvoluted Spectral Features

In addition to scoring the raw observed spectrum at the apex, Osprey computes a second set of spectral scores from a **deconvoluted spectrum**. This helps discriminate true peptides from false matches caused by co-eluting interference.

### What is Deconvolution?

The deconvoluted spectrum represents what the observed signal would look like if only this peptide were present, based on the regression coefficients.

```
For each scan near the apex (±2 scans):
  deconv_contribution[scan] = coefficient[scan] × library_spectrum

Aggregate across scans:
  deconvoluted_spectrum = Σ deconv_contribution[scan] (weighted by coefficient)
```

### Deconvoluted Features

These features use the same scoring algorithms as the "mixed" spectral features, but applied to the deconvoluted (cleaned) spectrum:

| Feature | Field | Description |
|---------|-------|-------------|
| Hyperscore (deconv) | `hyperscore_deconv` | X!Tandem hyperscore on deconvoluted spectrum |
| XCorr (deconv) | `xcorr_deconv` | Cross-correlation on deconvoluted spectrum |
| Dot product (deconv) | `dot_product_deconv` | LibCosine on deconvoluted spectrum |
| Dot product SMZ (deconv) | `dot_product_smz_deconv` | LibCosine with m/z² weighting on deconvoluted spectrum |
| Fragment coverage (deconv) | `fragment_coverage_deconv` | Fraction of fragments matched in deconvoluted spectrum |
| Sequence coverage (deconv) | `sequence_coverage_deconv` | Backbone coverage in deconvoluted spectrum |
| Consecutive ions (deconv) | `consecutive_ions_deconv` | Longest b/y run in deconvoluted spectrum |
| Top-3 matches (deconv) | `top3_matches_deconv` | Top-3 predicted fragments matched in deconvoluted spectrum |

### Why Both Mixed and Deconvoluted Scores?

| Scenario | Mixed Score | Deconvoluted Score |
|----------|-------------|-------------------|
| True peptide, no interference | High | High |
| True peptide, with interference | Medium (diluted) | High (cleaned) |
| False positive (interference attributed) | High | Low (doesn't match library) |
| Noise | Low | Low |

The combination of both score types helps Mokapot identify cases where a peptide appears to match the mixed spectrum but doesn't actually contribute to the signal.

---

## Chromatographic Features (from Coefficient Time Series)

These features describe the shape and quality of the deconvoluted chromatographic peak.

**Features sent to Mokapot (6 features):**

| Feature | Field | Description |
|---------|-------|-------------|
| Peak apex | `peak_apex` | Maximum coefficient value across the elution profile |
| Peak area | `peak_area` | Integrated area under the coefficient curve |
| Peak width | `peak_width` | Full width at half maximum (FWHM) in minutes |
| Contributing scans | `n_contributing_scans` | Number of scans with non-zero coefficient |
| Coefficient stability | `coefficient_stability` | Variance of coefficients near apex (lower = more stable) |
| Relative coefficient | `relative_coefficient` | This peptide's coefficient / sum of all coefficients in spectrum |

**Additional features in FeatureSet (not currently sent to Mokapot):**

| Feature | Field | Description |
|---------|-------|-------------|
| Peak symmetry | `peak_symmetry` | Leading width / trailing width at half-max |
| Peak sharpness | `peak_sharpness` | Rate of coefficient change at boundaries |
| Peak prominence | `peak_prominence` | Apex / baseline estimate |
| EMG fit quality | `emg_fit_quality` | Proxy for exponentially modified Gaussian fit |

---

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

## Mokapot PIN Feature Summary (24 Features)

The following 24 features are written to the PIN file for Mokapot:

### Chromatographic Features (7)

| # | Field | Description | Range |
|---|-------|-------------|-------|
| 1 | `peak_apex` | Maximum coefficient value | 0+ |
| 2 | `peak_area` | Integrated coefficient area | 0+ |
| 3 | `peak_width` | FWHM in minutes | 0+ |
| 4 | `n_contributing_scans` | Scans with non-zero coefficient | 0-N |
| 5 | `coefficient_stability` | Variance near apex | 0+ |
| 6 | `relative_coefficient` | Fraction of total signal | 0-1 |
| 7 | `explained_intensity` | Fraction of observed intensity matched | 0-1 |

### Spectral Features - Mixed/Observed Spectrum (8)

| # | Field | Requires Library Intensities | Range |
|---|-------|------------------------------|-------|
| 8 | `hyperscore` | No (fragment types only) | 0-50 |
| 9 | `xcorr` | No (unit intensity) | -0.5 to 2.0 |
| 10 | `dot_product` | Yes | 0-1 |
| 11 | `dot_product_smz` | Yes | 0-1 |
| 12 | `fragment_coverage` | Yes (fragment m/z) | 0-1 |
| 13 | `sequence_coverage` | Yes (ion types) | 0-1 |
| 14 | `consecutive_ions` | Yes (ion types) | 0-N |
| 15 | `top3_matches` | Yes | 0-3 |

### Spectral Features - Deconvoluted Spectrum (8)

| # | Field | Requires Library Intensities | Range |
|---|-------|------------------------------|-------|
| 16 | `hyperscore_deconv` | No | 0-50 |
| 17 | `xcorr_deconv` | No | -0.5 to 2.0 |
| 18 | `dot_product_deconv` | Yes | 0-1 |
| 19 | `dot_product_smz_deconv` | Yes | 0-1 |
| 20 | `fragment_coverage_deconv` | Yes | 0-1 |
| 21 | `sequence_coverage_deconv` | Yes | 0-1 |
| 22 | `consecutive_ions_deconv` | Yes | 0-N |
| 23 | `top3_matches_deconv` | Yes | 0-3 |

### Derived Features (1)

| # | Field | Description | Range |
|---|-------|-------------|-------|
| 24 | `rt_deviation` | Observed apex RT - predicted RT | minutes |

---

## Features in FeatureSet Not Sent to Mokapot

These features are computed and stored in `FeatureSet` but are NOT currently included in the PIN file:

| Field | Description | Reason Not Included |
|-------|-------------|---------------------|
| `emg_fit_quality` | EMG fit R² | Not yet implemented |
| `peak_symmetry` | Leading/trailing width ratio | Low discrimination power |
| `peak_sharpness` | Boundary slope | Redundant with peak_width |
| `peak_prominence` | Apex / baseline | Redundant with peak_apex |
| `rt_deviation_normalized` | RT deviation / FWHM | Redundant with rt_deviation |
| `spectral_contrast_angle` | Angle-based similarity | Redundant with dot_product |
| `pearson_correlation` | Intensity correlation | Redundant with dot_product |
| `spearman_correlation` | Rank correlation | Low discrimination power |
| `base_peak_rank` | Rank of base peak | Low discrimination power |
| `n_competitors` | Candidates per spectrum | Context feature, may add later |
| `local_peptide_density` | Average competitors | Context feature, may add later |
| `spectral_complexity` | Regression complexity | Context feature, may add later |
| `regression_residual` | Unexplained signal | Context feature, may add later |
| `precursor_intensity` | MS1 intensity | Optional, not always available |
| `modification_count` | Number of modifications | May add for PTM analysis |

## Implementation Files

| File | Purpose |
|------|---------|
| `crates/osprey-scoring/src/lib.rs` | SpectralScorer (LibCosine, XCorr, Hyperscore), FeatureExtractor |
| `crates/osprey-scoring/src/batch.rs` | BLAS-accelerated batch scoring (calibration phase) |
| `crates/osprey/src/pipeline.rs` | score_run() - orchestrates scoring after regression |
| `crates/osprey-core/src/types.rs` | FeatureSet struct (all available features) |
| `crates/osprey-fdr/src/mokapot.rs` | PIN file output (24 features), Mokapot runner |

## Feature Weight Analysis

Use `scripts/inspect_mokapot_weights.py` to see which features Mokapot finds most discriminative:

```bash
python scripts/inspect_mokapot_weights.py mokapot.model.pkl
```

This helps identify:

- Which features contribute most to target/decoy discrimination
- Features with near-zero weights (candidates for removal to speed up processing)
- Unexpected feature importance patterns

## TODO

1. **Consider apex refinement**: Score top-N coefficient spectra and pick best spectral match
2. **EMG peak fitting**: Replace symmetry-based heuristic with true Levenberg-Marquardt EMG fit
3. **Feature reduction**: Use Mokapot weights to identify and remove low-value features
4. **Context features**: Evaluate adding n_competitors, regression_residual to PIN file
