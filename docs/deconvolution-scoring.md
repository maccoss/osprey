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

## Spectral Scores (Sent to Mokapot)

All spectral scores are computed at the apex spectrum (the observed MS2 spectrum closest in RT to the coefficient apex). Only XCorr and consecutive ions are sent to the Mokapot PIN file for FDR control. Additional scores (LibCosine, hyperscore, etc.) are computed and stored in the FeatureSet for diagnostics and blib output but are not used for FDR.

**Important**: All intensity-based spectral similarity scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks. A missing peak that should be present is a strong discriminator between good and bad matches.

### 1. XCorr (Cross-Correlation)

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

### 2. Consecutive Ions

**PIN feature**: `consecutive_ions` (mixed), `consecutive_ions_deconv` (deconvoluted)

The length of the longest run of consecutive b or y ions matched in the spectrum.

#### Formula

```
For each ion series (b-ions, y-ions):
  Find the longest contiguous run of matched ions
  e.g., b3, b4, b5 matched → run length = 3

consecutive_ions = max(longest_b_run, longest_y_run)
```

#### Score range: 0 to n (peptide length - 1)

#### Why it matters
- True peptide matches tend to have long runs of sequential fragment ions
- Random noise matches tend to have scattered, non-consecutive ion matches
- Simple integer feature that captures ion series completeness

### 3. Explained Intensity

**PIN feature**: `explained_intensity`

The fraction of total observed intensity in the spectrum that is accounted for by matched library fragments.

#### Formula

```
matched_intensity = Σ intensity of observed peaks that match a library fragment
total_intensity = Σ intensity of all observed peaks (above threshold)

explained_intensity = matched_intensity / total_intensity
```

#### Score range: 0 to 1 (1 = all signal comes from predicted fragments)

### 4. Delta RT (RT Deviation)

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

---

## Deconvoluted Spectral Features

In addition to scoring the raw observed spectrum at the apex, Osprey computes spectral scores from a **deconvoluted spectrum**. This helps discriminate true peptides from false matches caused by co-eluting interference.

### What is Deconvolution?

The deconvoluted spectrum represents what the observed signal would look like if only this peptide were present, based on the regression coefficients.

```
For each scan near the apex (±2 scans):
  deconv_contribution[scan] = coefficient[scan] × library_spectrum

Aggregate across scans:
  deconvoluted_spectrum = Σ deconv_contribution[scan] (weighted by coefficient)
```

### Features Sent to Mokapot (Deconvoluted)

| Feature | Field | Description |
|---------|-------|-------------|
| XCorr (deconv) | `xcorr_deconv` | Cross-correlation on deconvoluted spectrum |
| Consecutive ions (deconv) | `consecutive_ions_deconv` | Longest b/y run in deconvoluted spectrum |

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

**Features sent to Mokapot (8 features):**

| Feature | Field | Description |
|---------|-------|-------------|
| Peak apex | `peak_apex` | Maximum coefficient value across the elution profile |
| Peak area | `peak_area` | Integrated area under the coefficient curve |
| Peak width | `peak_width` | FWHM from Tukey median polish elution profile (minutes) |
| Coefficient stability | `coefficient_stability` | CV of coefficients near apex |
| Relative coefficient | `relative_coefficient` | This peptide's coefficient / sum of all coefficients in spectrum |
| Explained intensity | `explained_intensity` | Fraction of observed intensity explained by matched fragments |
| Signal to noise | `signal_to_noise` | Peak apex / noise estimate from coefficient time series |
| XIC signal to noise | `xic_signal_to_noise` | S/N from best-correlated fragment XIC |

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
│     ├─ Extract fragment XICs (top 6 by library intensity)
│     │
│     ├─ Run Tukey median polish on 6×N XIC matrix
│     │     → Column effects = elution profile (robust peak shape)
│     │     → Row effects = interference-free fragment intensities
│     │     → FWHM from elution profile → peak boundaries
│     │
│     ├─ Compute fragment co-elution features within peak boundaries
│     │
│     ├─ Compute elution-weighted cosine (libcosine at each scan, weighted by coef²)
│     │
│     ├─ Compute MS1 features (precursor co-elution, isotope envelope)
│     │
│     ├─ Find apex spectrum (closest to apex_rt)
│     │
│     ├─ Compute spectral scores at apex spectrum:
│     │     ├─ XCorr (Comet-style, mixed and deconvoluted)
│     │     ├─ Consecutive ions (mixed and deconvoluted)
│     │     └─ Explained intensity
│     │
│     ├─ Compute Tukey median polish scores:
│     │     ├─ median_polish_cosine (row effects vs library, sqrt space)
│     │     ├─ median_polish_rsquared (R² in sqrt space)
│     │     └─ median_polish_residual_ratio (residual ratio in linear space)
│     │
│     ├─ Compute mass accuracy features (ppm errors)
│     │
│     ├─ Compute Percolator-style features (peptide_length, missed_cleavages, etc.)
│     │
│     └─ Return ScoredEntry with all features
│
└─ 3. Output: Vec<ScoredEntry> (targets + decoys)
       → 37 features written to PIN file for Mokapot
```

## Mokapot PIN Feature Summary (37 Features)

See [docs/README.md](README.md) for the complete 37-feature table organized by category.

The features fall into these groups:
- **Ridge regression (8)**: peak_apex, peak_area, peak_width, coefficient_stability, relative_coefficient, explained_intensity, signal_to_noise, xic_signal_to_noise
- **Spectral matching - mixed (2)**: xcorr, consecutive_ions
- **Spectral matching - deconvoluted (2)**: xcorr_deconv, consecutive_ions_deconv
- **RT deviation (1)**: rt_deviation
- **Fragment co-elution (9)**: coelution sum/min/n_positive, 6 per-fragment correlations
- **Elution-weighted cosine (1)**: cosine at each scan, weighted by coefficient^2
- **Mass accuracy (3)**: mean/abs_mean/std of fragment mass errors (ppm)
- **Percolator-style (6)**: abs_rt_deviation, peptide_length, missed_cleavages, ln_num_candidates, coef_zscore, coef_zscore_mean
- **MS1 features (2)**: ms1_precursor_coelution, ms1_isotope_cosine
- **Tukey median polish (3)**: median_polish_cosine, median_polish_rsquared, median_polish_residual_ratio

### Design Principles

- All intensity-based scores (cosine, Pearson, Spearman) include ALL library fragments within the spectrum's mass range, using 0 for unmatched peaks
- sqrt preprocessing (Poisson noise model) used for cosine scoring
- Fragment co-elution features capture whether individual fragment XICs track the coefficient time series
- Tukey median polish provides robust peak boundaries and interference-free fragment scoring

## Implementation Files

| File | Purpose |
|------|---------|
| `crates/osprey-scoring/src/lib.rs` | SpectralScorer (LibCosine, XCorr), FeatureExtractor, Tukey median polish, fragment co-elution |
| `crates/osprey-scoring/src/batch.rs` | BLAS-accelerated batch scoring (calibration phase) |
| `crates/osprey/src/pipeline.rs` | score_run() - orchestrates scoring after regression |
| `crates/osprey-core/src/types.rs` | FeatureSet struct (37 PIN features + additional computed fields) |
| `crates/osprey-fdr/src/mokapot.rs` | PIN file output (37 features), Mokapot runner |

## Feature Weight Analysis

Use `scripts/inspect_mokapot_weights.py` to see which features Mokapot finds most discriminative:

```bash
python scripts/inspect_mokapot_weights.py mokapot.model.pkl
```

This helps identify:

- Which features contribute most to target/decoy discrimination
- Features with near-zero weights (candidates for removal to speed up processing)
- Unexpected feature importance patterns
