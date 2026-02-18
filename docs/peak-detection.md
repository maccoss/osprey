# Peak Detection

Osprey uses a multi-stage peak detection system that operates on fragment ion chromatograms (XICs). Peak detection determines the elution apex, integration boundaries, and chromatographic features that feed into FDR control.

## Overview

Osprey uses **DIA-NN-style XIC peak detection** (`detect_xic_peak()`) with smoothed apex finding, valley-based boundary detection, and asymmetric FWHM capping. This is complemented by **Tukey median polish**, which decomposes fragment XIC matrices to produce robust elution profiles and refined peak boundaries.

```
Peak Detection Pipeline
────────────────────────

  Fragment XICs → pairwise correlation → reference XIC (best-correlated fragment)
  Reference XIC → detect_xic_peak()   → apex, boundaries, area, S/N
  Fragment XICs → Tukey median polish  → refined FWHM boundaries
                                      → scoring features
```

---

## XIC Peak Detection

### Algorithm: `detect_xic_peak()`

The primary peak detector used for fragment XICs and calibration scoring. Combines smoothing, adaptive apex finding, valley-based boundary detection, and asymmetric FWHM capping.

```
detect_xic_peak(xic, min_height, peak_boundary, expected_rt)

  Step 1: Smooth XIC with weighted moving average
  Step 2: Find apex (local maximum nearest to expected RT)
  Step 3: Walk boundaries outward using valley detection
  Step 4: Cap boundaries using asymmetric FWHM
  Step 5: Compute signal-to-noise from background flanking regions
```

#### Step 1: Smoothing

Applies a 5-point Savitzky-Golay quadratic filter: `[-3, 12, 17, 12, -3] / 35`

```
Interior (i >= 2 and i < n-2):
  smoothed[i] = (-3*v[i-2] + 12*v[i-1] + 17*v[i] + 12*v[i+1] - 3*v[i+2]) / 35

Endpoints (first 2 and last 2 points): left unsmoothed
Short series (< 5 points): returned unsmoothed
Clamping: negative values from the filter are clamped to zero
```

The SG quadratic filter fits a local parabola at each point, preserving peak position and shape better than triangular or moving average kernels. It is the same filter used in calibration coelution scoring (`batch.rs`).

#### Step 2: Apex Finding

```
find_apex(smoothed, min_height, expected_rt, xic):
  1. Collect all local maxima above min_height:
     - smoothed[i] >= smoothed[i-1] AND smoothed[i] >= smoothed[i+1]
     - Also check endpoints as potential maxima
  2. If expected_rt provided: return local maximum nearest to expected RT
  3. Otherwise: return global maximum
```

Using the nearest-to-expected-RT strategy ensures the correct peak is selected when multiple chromatographic peaks are present (e.g., from co-eluting interferences).

#### Step 3: Valley-Based Boundary Detection

Two independent boundary walks from the apex, each stopping at the earlier of:

1. **Intensity threshold**: `smoothed[i] < apex_intensity / peak_boundary`
   - Default `peak_boundary = 5.0` → stops at 20% of apex intensity
2. **Valley detection**: A local minimum that satisfies both:
   - `< 50%` of apex intensity
   - `< 50%` of the neighboring rising peak

```
Walk left from apex:
  For i = apex-1 down to 0:
    If smoothed[i] < threshold → stop (threshold boundary)
    If smoothed[i] < smoothed[i-1]:  # valley detected
      If smoothed[i] < 0.5 * apex AND smoothed[i] < 0.5 * smoothed[i-1]:
        → stop (valley boundary)

Walk right from apex: (symmetric, scanning right)
```

Valley detection prevents over-extension into adjacent peaks by recognizing the saddle point between two overlapping peaks.

#### Step 4: Asymmetric FWHM Capping

**Problem**: Chromatographic tailing (slow exponential decay on the trailing edge) causes valley detection to extend too far because the tail decays slowly but never reaches a valley or the 20% threshold.

**Solution**: Compute separate left and right half-widths at half-maximum, then cap boundaries at `2.0 x half_width`:

```
compute_asymmetric_half_widths(smoothed, xic, apex_idx):

  half_max = smoothed[apex_idx] / 2.0

  Scan left from apex to find where smoothed first crosses half_max:
    Linear interpolate between bracketing points for exact crossing RT
    left_half_width = apex_rt - crossing_rt

  Scan right similarly:
    right_half_width = crossing_rt - apex_rt

  Return (left_half_width, right_half_width)
```

**Applying the cap** (factor = 2.0, chosen for ~98% Gaussian area coverage):

```
If valley_start_rt < apex_rt - 2.0 * left_half_width:
  Tighten start_rt to apex_rt - 2.0 * left_half_width

If valley_end_rt > apex_rt + 2.0 * right_half_width:
  Tighten end_rt to apex_rt + 2.0 * right_half_width
```

**Key insight**: Separate left/right half-widths naturally capture chromatographic tailing where `right_hw >> left_hw`. A symmetric FWHM would either over-extend the leading edge or under-capture the trailing edge.

```
Intensity
  |       ╭──╮
  |      ╱    ╲      ← Tailing peak: fast rise, slow decay
  |     ╱      ╲╲
  |    ╱        ╲╲╲
  |───╱──────────╲╲╲╲───── valley too far right
  |  ╱            ╲╲╲╲╲
  +──|──|─────|──────|────→ RT
     ^  ^     ^      ^
     |  apex  |      valley boundary (too wide)
     |        FWHM cap (2× right_hw)
     FWHM cap (2× left_hw)
```

#### Step 5: Signal-to-Noise

```
compute_snr(raw_intensities, apex_idx, start_idx, end_idx):
  background_left  = 5 points immediately before peak start
  background_right = 5 points immediately after peak end
  bg_mean = mean(background_left + background_right)
  bg_sd   = sd(background_left + background_right)
  SNR     = (apex_intensity - bg_mean) / bg_sd
```

Uses **raw (unsmoothed) intensities** for honest noise estimation.

### Output

```rust
pub struct XICPeakBounds {
    pub apex_rt: f64,           // Peak apex time (minutes)
    pub apex_intensity: f64,    // Intensity at apex (smoothed)
    pub apex_index: usize,      // Index in XIC array
    pub start_rt: f64,          // Start boundary time
    pub end_rt: f64,            // End boundary time
    pub start_index: usize,     // Index at start boundary
    pub end_index: usize,       // Index at end boundary
    pub area: f64,              // Trapezoidal area within boundaries
    pub signal_to_noise: f64,   // Background-subtracted S/N
}
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_height` | 0.01 | Minimum smoothed intensity for apex |
| `peak_boundary` | 5.0 | Apex/boundary intensity ratio (5.0 = stop at 20%) |
| `expected_rt` | Library RT | Hint for apex selection among multiple maxima |
| FWHM cap factor | 2.0 | Boundary cap at 2x half-width from apex |

---

## Tukey Median Polish

After initial peak detection, Osprey uses Tukey median polish to decompose the fragment XIC matrix into additive components. This provides both refined peak boundaries and scoring features.

### Mathematical Model

The matrix of fragment XIC intensities (6 fragments x N scans) is decomposed in log space:

```
ln(Observed[f,s]) = μ + α_f + β_s + ε_fs

  f = fragment index (top 6 by library intensity)
  s = scan index
  μ     = overall effect (grand median)
  α_f   = row effects (data-derived fragment relative intensities)
  β_s   = column effects (shared elution profile shape)
  ε_fs  = residuals (interference + noise)
```

### Algorithm

```
Initialize: matrix[f,s] = ln(intensity_fs), NaN for zeros (missing data)
            row_effects = [0; n_fragments]
            col_effects = [0; n_scans]
            overall = 0

For iteration = 1..20:
  Row sweep:
    For each fragment f:
      m = nanmedian(residuals[f, :])
      residuals[f, :] -= m
      row_effects[f] += m

  Column sweep:
    For each scan s:
      m = nanmedian(residuals[:, s])
      residuals[:, s] -= m
      col_effects[s] += m

  Update overall:
    m = median(col_effects)
    col_effects -= m
    overall += m

  If max(|change|) < 1e-4: converged, stop

Elution profile: (RT_s, exp(overall + col_effects[s]))  for each scan s
```

**NaN handling**: Zero-intensity cells (undetected fragment at that scan) are set to NaN and skipped by `nanmedian`. This naturally handles missing data without imputation.

### How Column Effects Determine Peak Boundaries

Column effects (β_s) represent the **shared elution profile** across all fragments. Because the median operation suppresses interference on individual transitions, this profile is more robust than any single fragment XIC.

```
From Polish elution profile:
  1. Find FWHM via linear interpolation at half-maximum
  2. Convert to sigma: σ = FWHM / 2.355
  3. Boundaries: [apex - 1.96σ, apex + 1.96σ]  (~95% Gaussian interval)
```

These boundaries are used as **fragment peak bounds** in the blib output for Skyline. Skyline uses these per-file RT boundaries to define integration regions for quantification.

### Scoring Features from Median Polish

Three features are extracted from the decomposition:

| Feature | Computation | Interpretation |
|---------|-------------|----------------|
| `median_polish_cosine` | cosine(sqrt(exp(μ + α_f)), sqrt(library_intensity_f)) | Do row effects match the library? High for targets. |
| `median_polish_rsquared` | 1 - SS_residual/SS_total in sqrt space | How well does the additive model fit? Higher = cleaner co-elution. |
| `median_polish_residual_ratio` | Σ\|obs - pred\| / Σ obs in linear space | Fraction of signal unexplained by model. Lower = better. |

**R² computation details:**
- Uses sqrt preprocessing (appropriate for Poisson-like count data)
- Predicted values: `exp(μ + α_f + β_s)` converted to linear space, then sqrt
- Zero observed cells contribute `sqrt(pred)²` to residual (penalizes hallucinated fragments)

---

## Data Flow

Peak detection operates on **fragment ion chromatograms** (XICs) extracted directly from DIA spectra:

```
Fragment XIC Extraction
  → Extract top 6 fragment XICs within RT window (ppm tolerance)
  → Compute pairwise Pearson correlations between all fragment pairs

Reference XIC Selection
  → Fragment with highest mean correlation = reference XIC

Peak Detection on Reference XIC
  → detect_xic_peak(reference_xic, expected_rt)
  → Returns: apex, boundaries, area, S/N

Feature Extraction (within peak boundaries)
  → Pairwise coelution features from fragment correlations
  → Peak shape features from reference XIC
  → Spectral matching at apex scan
  → Tukey median polish on fragment XIC matrix
  → MS1 features from nearest MS1 spectrum
```

### Peak Shape Features

| Feature | Description |
|---------|-------------|
| `peak_apex` | Reference XIC intensity at apex (background-subtracted) |
| `peak_area` | Integrated area within boundaries |
| `peak_width` | end_rt - start_rt (minutes) |
| `peak_symmetry` | Leading area / trailing area around apex (capped at 10.0) |
| `signal_to_noise` | S/N from `detect_xic_peak()` |
| `n_scans` | Number of spectrum points within peak boundaries |
| `peak_sharpness` | Average absolute slope at peak edges (intensity/minute) |

### Peak Symmetry

```
peak_symmetry = left_area / right_area

Where:
  left_area  = trapezoidal area from start_index to apex_index
  right_area = trapezoidal area from apex_index to end_index
```

- Symmetric peak: ~1.0
- Fronting peak: < 1.0
- Tailing peak: > 1.0 (capped at 10.0)

### Peak Sharpness

```
left_slope  = (apex_intensity - start_intensity) / (apex_rt - start_rt)
right_slope = (apex_intensity - end_intensity) / (end_rt - apex_rt)
peak_sharpness = (|left_slope| + |right_slope|) / 2
```

Steeper slopes indicate a sharper, better-defined peak. Broad, noisy peaks have low sharpness.

---

## Fragment Co-Elution Analysis

After peak detection determines the integration boundaries, fragment co-elution features are computed **within the peak region only**:

Pairwise correlations between fragment XICs:

```
For each pair of top 6 fragment XICs (within peak window):
  Compute Pearson correlation
Per-fragment: average correlation with all other fragments
Return: sum, min, max, n_coeluting, n_pairs, per-fragment averages
```

---

## Peak Boundaries for Blib Output

The final peak boundaries written to the blib file for Skyline come from **Tukey median polish** (when successful) or fall back to the initial detection boundaries:

```
Priority:
  1. Tukey median polish FWHM → σ = FWHM/2.355, boundaries = apex ± 1.96σ
  2. detect_xic_peak() boundaries (valley detection + FWHM capping)
```

**Important**: Do NOT use integrated area for quantification. XIC intensities are relative measures. Skyline handles quantification using the extracted ion chromatograms within the peak boundaries provided by Osprey.

---

## Calibration Scoring

During RT calibration discovery, the same XIC peak detection is used to score calibration matches:

```
run_coelution_calibration_scoring():
  For each sampled library entry:
    Extract top-6 fragment XICs from candidate spectra
    Compute pairwise correlations
    detect_xic_peak() on reference XIC → apex, boundaries, S/N
    Score at apex: libcosine, hyperscore, mass errors
    → CalibrationMatch with correlation_score, S/N, etc.
```

This provides the confident matches used for RT and mass calibration fitting.

---

## Example

```
Peptide: PEPTIDEK (charge 2+)
Library RT: 25.3 min
Expected RT: 24.8 min

Fragment XICs extracted (top 6 by library intensity):
  y5+1, y6+1, y7+1, b3+1, b4+1, y4+1

Pairwise correlations:
  y5-y6: 0.95, y5-y7: 0.92, y6-y7: 0.94, ...
  Reference XIC: y6+1 (highest mean correlation: 0.93)

detect_xic_peak(reference_xic, expected_rt=24.8):
  Smoothing: [0.25, 0.5, 0.25] kernel
  Apex: RT=24.9, intensity=4500 (nearest local max to 24.8)
  Valley walk: start=24.1, end=25.8
  FWHM cap: left_hw=0.35, right_hw=0.45
    → cap_start = 24.9 - 0.70 = 24.2 (no change, valley tighter)
    → cap_end   = 24.9 + 0.90 = 25.8 (no change, valley matches)
  Final boundaries: [24.1, 25.8]
  SNR: 12.5

Features:
  peak_apex: 4500, peak_area: 28500, peak_width: 1.7 min
  peak_symmetry: 0.85, n_scans: 14, peak_sharpness: 3200
  signal_to_noise: 12.5
```

---

## Multiple Peaks

A peptide may have multiple peaks due to:
- **Isomers** eluting at different RTs
- **False positive signals** from co-eluting peptides with similar spectra
- **Peak splitting** from chromatographic artifacts

`detect_xic_peak()` selects the local maximum nearest to the expected RT. If no peak meets the minimum height threshold, the peptide is not reported.

---

## Testing

Test coverage in `crates/osprey-chromatography/src/lib.rs`:

| Test | Description |
|------|-------------|
| `test_peak_detector` | Simple Gaussian peak detection via threshold-crossing |
| `test_find_best_peak` | Peak selection within RT tolerance |
| `test_fwhm_cap_symmetric_peak` | Symmetric peak: FWHM cap should not tighten boundaries |
| `test_fwhm_cap_tailing_peak` | Tailing peak: FWHM cap tightens trailing boundary |
| `test_fwhm_cap_preserves_valley` | Adjacent peaks: valley boundary preserved (not capped) |
| `test_fwhm_cap_fallback` | Graceful fallback when FWHM cannot be computed |
| `test_asymmetric_half_widths` | Half-width computation for asymmetric peaks |

Test helpers generate synthetic data:
- `make_gaussian_xic(center, sigma, n)`: N(center, sigma) peak
- `make_tailing_xic(center, sigma, tail_factor, n)`: Gaussian rise + exponential decay

---

## Implementation

| File | Purpose |
|------|---------|
| `crates/osprey-chromatography/src/lib.rs` | `PeakDetector::detect()`, `find_best_peak()`, `detect_xic_peak()`, `smooth_weighted_avg()`, `find_apex()`, `walk_boundary_left/right()`, `compute_asymmetric_half_widths()`, `compute_snr()`, `compute_fwhm_interpolated()` |
| `crates/osprey-core/src/types.rs` | `PeakBoundaries`, `PeakQuality`, `XICPeakBounds` structs |
| `crates/osprey-scoring/src/lib.rs` | `tukey_median_polish()`, `TukeyMedianPolishResult`, `compute_fragment_coelution()`, `FeatureExtractor`, scoring features from peak shape |
| `crates/osprey/src/pipeline.rs` | Peak detection calls in coelution pipeline, Tukey median polish integration, boundary propagation to blib output |
