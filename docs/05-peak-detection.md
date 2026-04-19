# Peak Detection

Osprey uses a multi-stage peak detection system that operates on fragment ion chromatograms (XICs). Peak detection determines the elution apex, integration boundaries, and chromatographic features that feed into FDR control.

## Overview

Osprey's primary peak detection method is **CWT (Continuous Wavelet Transform) consensus peak detection** — a multi-transition approach that convolves each fragment XIC with a Mexican Hat wavelet, then takes the pointwise median across all transitions to produce a consensus signal. Peaks are detected in the consensus signal, which is naturally robust to single-fragment interference.

CWT consensus is complemented by **Tukey median polish**, which decomposes the fragment XIC matrix into additive components for robust scoring features and peak boundary refinement.

```text
Peak Detection Pipeline
────────────────────────

  Fragment XICs → Mexican Hat wavelet convolution (per fragment)
                → Pointwise median (consensus CWT signal)
                → Local maxima → peak apices
                → Zero-crossings ±2σ extension (with valley guard) → boundaries
                → Reference signal (sum of raw XICs) → area, S/N

  Fragment XICs → Tukey median polish → scoring features
                                      → peak boundaries (blib output)
```

---

## CWT Consensus Peak Detection

### Motivation

The challenge in DIA peak detection is distinguishing true co-eluting fragment signals from interference caused by co-isolated peptides. Individual fragment XICs may contain interfering peaks from other peptides. A multi-transition consensus approach — where a peak is only called if the **majority** of fragments simultaneously show peak-like shapes — naturally rejects single-fragment interference.

The Mexican Hat wavelet is a natural matched filter for Gaussian-like chromatographic peaks: it responds positively to peak shapes and negatively to flat or monotonic regions. Taking the pointwise median of CWT coefficients across fragments means the consensus is only high where most fragments agree.

### Algorithm: `detect_cwt_consensus_peaks()`

```text
detect_cwt_consensus_peaks(xics, min_consensus_height)

  Step 1: Estimate scale parameter (σ) from fragment XICs
  Step 2: Generate Mexican Hat wavelet kernel
  Step 3: Convolve each fragment XIC with the kernel
  Step 4: Compute pointwise median consensus across all transitions
  Step 5: Find local maxima in consensus signal
  Step 6: Define boundaries via zero-crossings with ±2σ extension
  Step 7: Compute area and S/N from reference signal (sum of raw XICs)
```

#### Step 1: Scale Estimation

The CWT scale parameter (σ) controls which peak widths the wavelet is tuned to detect. Osprey estimates σ from the data:

```text
estimate_cwt_scale(xics):
  For each fragment XIC with detectable signal:
    Find apex and measure FWHM via linear interpolation at half-maximum
  sigma = median(FWHM values) / 2.355    (Gaussian: FWHM = 2.355σ)
  Clamp to [2.0, 20.0] scan units
  Fallback: 4.0 scans if no FWHM can be estimated
```

#### Step 2: Mexican Hat Wavelet Kernel

The Mexican Hat (Ricker) wavelet is the negative normalized second derivative of a Gaussian:

```text
ψ(t) = (2 / sqrt(3σ) π^(1/4)) × (1 - (t/σ)²) × exp(-t²/(2σ²))
```

Properties:

- **Positive center, negative tails**: Responds positively to peaks, negatively to flanking regions
- **Zero-mean**: After discretization, a DC offset correction ensures zero response to constant signals
- **Symmetric**: Equal sensitivity to leading and trailing edges
- **Kernel radius**: `ceil(5σ)` points on each side, capturing >99.99% of wavelet energy

#### Step 3: Wavelet Convolution

Each fragment XIC is convolved independently with the Mexican Hat kernel using direct "same"-size convolution (zero-padded at edges). The output length equals the input length.

Direct convolution is O(N×K) per fragment, but with typical XIC lengths (30-200 scans) and kernel sizes (11-51 points), this is trivially fast.

#### Step 4: Pointwise Median Consensus

At each scan position, the median CWT coefficient across all fragments is computed:

```text
For each scan s:
  consensus[s] = median(cwt_coeffs[f][s] for f in 0..n_fragments)
```

The **median** (not mean) is the key to interference rejection:

- If 5/6 fragments have a peak shape at scan s, the median is positive
- If only 1/6 fragments has interference at scan s, the median suppresses it
- A peak must be supported by the **majority** of transitions to survive

#### Step 5: Apex Detection

Local maxima in the consensus signal above `min_consensus_height` are identified:

- Interior: `consensus[i] > consensus[i-1]` AND `consensus[i] > consensus[i+1]`
- Endpoints checked as potential maxima
- Sorted by consensus coefficient descending (strongest peaks first)

#### Step 6: Boundary Extension (±2σ with Valley Guard)

Peak boundaries are determined in two stages:

**Stage 1 — Zero-crossings (±σ):**
Walk outward from the apex to where the consensus signal crosses zero. For a Gaussian peak, zero-crossings occur at approximately ±σ from the apex, capturing ~68% of the peak area.

**Stage 2 — Extend to ±2σ (~95% coverage):**
The zero-crossing distance gives an asymmetric estimate of σ (left_σ = apex - left_zc, right_σ = right_zc - apex). Boundaries are extended to ±2σ for ~95% Gaussian area coverage:

```text
target_start = apex - 2 × left_σ
target_end   = apex + 2 × right_σ
```

**Valley guard:** During extension beyond the zero-crossing, the raw reference signal (sum of fragment intensities) is monitored. If the signal rises more than 5% of the apex value above a running minimum, we've entered a neighboring peak — stop at the valley minimum:

```text
Extend left from zero-crossing toward target_start:
  Track running_min of raw reference signal
  If ref_signal[i] - running_min > 0.05 × apex_intensity:
    Stop at valley minimum (running_min position)
  Otherwise: continue extending

Extend right: symmetric procedure
```

The valley guard prevents the ±2σ extension from bleeding into adjacent peaks while still capturing the full width of isolated peaks.

```text
Consensus CWT Signal
  |
  |           ╭──╮
  |          ╱    ╲
  |    ─────╱      ╲─────────    ← zero line
  |        ╱ˉˉˉˉˉˉˉˉ╲
  |       ╱            ╲
  +──|──|──|────|────|──|──|──→ scan index
     ^  ^  ^    ^    ^  ^  ^
     |  |  zc   apex zc |  |
     |  |               |  |
     |  target          target
     |  (2σ)            (2σ)
     valley guard
     stopped here
```

#### Step 7: Area and S/N from Reference Signal

The consensus CWT coefficients are used only for peak detection and boundary definition. Quantitative measures come from the **reference signal** (sum of raw fragment XIC intensities):

- **Apex**: The maximum of the reference signal within the boundary region
- **Area**: Trapezoidal integration of the reference signal within boundaries
- **S/N**: Computed from 5 flanking points on each side of the peak in the reference signal

### Output

```rust
pub struct XICPeakBounds {
    pub apex_rt: f64,           // Peak apex time (minutes)
    pub apex_intensity: f64,    // Intensity at apex (from reference signal)
    pub apex_index: usize,      // Index in XIC array
    pub start_rt: f64,          // Start boundary time
    pub end_rt: f64,            // End boundary time
    pub start_index: usize,     // Index at start boundary
    pub end_index: usize,       // Index at end boundary
    pub area: f64,              // Trapezoidal area within boundaries
    pub signal_to_noise: f64,   // Background-subtracted S/N
}
```

### Fallback Chain

If CWT consensus finds no peaks (e.g., too few fragments, very low signal), Osprey falls back to simpler methods:

```text
Priority:
  1. CWT consensus peak detection (multi-transition)
  2. Tukey median polish elution profile peak detection (detect_all_xic_peaks)
  3. Reference XIC peak detection (detect_all_xic_peaks on best-correlated fragment)
```

The fallback uses `detect_all_xic_peaks()`, which applies Savitzky-Golay smoothing, valley-based boundary detection, and asymmetric FWHM capping on a single signal.

### Where CWT Is Used

CWT consensus peak detection is used in both search phases:

| Phase | File | Purpose |
|-------|------|---------|
| Calibration | `osprey-scoring/src/batch.rs` | Peak detection for calibration matches |
| Main search | `osprey/src/pipeline.rs` | Peak detection for all precursors |

---

## Legacy: SG-Smoothed Peak Detection

### Algorithm: `detect_xic_peak()`

Used as a fallback when CWT consensus finds no peaks, and as the underlying detector in `detect_all_xic_peaks()`. Combines smoothing, adaptive apex finding, valley-based boundary detection, and asymmetric FWHM capping.

```text
detect_xic_peak(xic, min_height, peak_boundary, expected_rt)

  Step 1: Smooth XIC with weighted moving average
  Step 2: Find apex (local maximum nearest to expected RT)
  Step 3: Walk boundaries outward using valley detection
  Step 4: Cap boundaries using asymmetric FWHM
  Step 5: Compute signal-to-noise from background flanking regions
```

#### Step 1: Smoothing

Applies a 5-point Savitzky-Golay quadratic filter: `[-3, 12, 17, 12, -3] / 35`

```text
Interior (i >= 2 and i < n-2):
  smoothed[i] = (-3*v[i-2] + 12*v[i-1] + 17*v[i] + 12*v[i+1] - 3*v[i+2]) / 35

Endpoints (first 2 and last 2 points): left unsmoothed
Short series (< 5 points): returned unsmoothed
Clamping: negative values from the filter are clamped to zero
```

#### Step 2: Apex Finding

```text
find_apex(smoothed, min_height, expected_rt, xic):
  1. Collect all local maxima above min_height:
     - smoothed[i] >= smoothed[i-1] AND smoothed[i] >= smoothed[i+1]
     - Also check endpoints as potential maxima
  2. If expected_rt provided: return local maximum nearest to expected RT
  3. Otherwise: return global maximum
```

#### Step 3: Valley-Based Boundary Detection

Two independent boundary walks from the apex, each stopping at the earlier of:

1. **Intensity threshold**: `smoothed[i] < apex_intensity / peak_boundary`
   - Default `peak_boundary = 5.0` → stops at 20% of apex intensity
2. **Valley detection**: A local minimum that satisfies both:
   - `< 50%` of apex intensity
   - `< 50%` of the neighboring rising peak

#### Step 4: Asymmetric FWHM Capping

Compute separate left and right half-widths at half-maximum, then cap boundaries at `2.0 × half_width`:

```text
If valley_start_rt < apex_rt - 2.0 * left_half_width:
  Tighten start_rt to apex_rt - 2.0 * left_half_width

If valley_end_rt > apex_rt + 2.0 * right_half_width:
  Tighten end_rt to apex_rt + 2.0 * right_half_width
```

Separate left/right half-widths naturally capture chromatographic tailing where `right_hw >> left_hw`.

#### Step 5: Signal-to-Noise

```text
compute_snr(raw_intensities, apex_idx, start_idx, end_idx):
  background_left  = 5 points immediately before peak start
  background_right = 5 points immediately after peak end
  bg_mean = mean(background_left + background_right)
  bg_sd   = sd(background_left + background_right)
  SNR     = (apex_intensity - bg_mean) / bg_sd
```

Uses **raw (unsmoothed) intensities** for honest noise estimation.

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_height` | 0.01 | Minimum smoothed intensity for apex |
| `peak_boundary` | 5.0 | Apex/boundary intensity ratio (5.0 = stop at 20%) |
| `expected_rt` | Library RT | Hint for apex selection among multiple maxima |
| FWHM cap factor | 2.0 | Boundary cap at 2x half-width from apex |

---

## Tukey Median Polish

After peak detection, Osprey uses Tukey median polish to decompose the fragment XIC matrix into additive components. This provides both scoring features and peak boundaries for the blib output.

### Mathematical Model

The matrix of fragment XIC intensities (6 fragments x N scans) is decomposed in log space:

```text
ln(Observed[f,s]) = μ + α_f + β_s + ε_fs

  f = fragment index (top 6 by library intensity)
  s = scan index
  μ     = overall effect (grand median)
  α_f   = row effects (data-derived fragment relative intensities)
  β_s   = column effects (shared elution profile shape)
  ε_fs  = residuals (interference + noise)
```

### Algorithm

```text
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

```text
From Polish elution profile:
  1. Find FWHM via linear interpolation at half-maximum
  2. Convert to sigma: σ = FWHM / 2.355
  3. Boundaries: [apex - 1.96σ, apex + 1.96σ]  (~95% Gaussian interval)
```

These boundaries are used as **fragment peak bounds** in the blib output for Skyline.

### Scoring Features from Median Polish

Three features are extracted from the decomposition:

| Feature | Computation | Interpretation |
|---------|-------------|----------------|
| `median_polish_cosine` | cosine(sqrt(exp(μ + α_f)), sqrt(library_intensity_f)) | Do row effects match the library? High for targets. |
| `median_polish_rsquared` | 1 - SS_residual/SS_total in sqrt space | How well does the additive model fit? Higher = cleaner co-elution. |
| `median_polish_residual_ratio` | Σ\|obs - pred\| / Σ obs in linear space | Fraction of signal unexplained by model. Lower = better. |

---

## Data Flow

Peak detection operates on **fragment ion chromatograms** (XICs) extracted directly from DIA spectra:

```text
Fragment XIC Extraction
  → Extract top 6 fragment XICs within RT window (ppm tolerance)

CWT Consensus Peak Detection
  → Mexican Hat wavelet convolution (per fragment)
  → Pointwise median consensus
  → Local maxima → apex candidates
  → ±2σ boundary extension (with valley guard) → boundaries
  → Reference signal → area, S/N

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
| `signal_to_noise` | S/N from peak detection |
| `n_scans` | Number of spectrum points within peak boundaries |
| `peak_sharpness` | Average absolute slope at peak edges (intensity/minute) |

---

## Peak Boundaries for Blib Output

The final peak boundaries written to the blib file for Skyline come from **Tukey median polish** (when successful) or fall back to the CWT detection boundaries:

```text
Priority:
  1. Tukey median polish FWHM → σ = FWHM/2.355, boundaries = apex ± 1.96σ
  2. CWT consensus boundaries (zero-crossing ±2σ with valley guard)
  3. detect_xic_peak() boundaries (valley detection + FWHM capping)
```

**Important**: Do NOT use integrated area for quantification. XIC intensities are relative measures. Skyline handles quantification using the extracted ion chromatograms within the peak boundaries provided by Osprey.

---

## Peak Selection: Choosing Among Multiple Candidates

Peak detection and peak selection are **two separate stages**. CWT consensus detects candidate peaks; a pairwise correlation score selects the best one.

### Stage 1: Candidate Detection (CWT Consensus Median)

CWT consensus returns all detected peaks sorted by **consensus CWT coefficient** descending. The consensus coefficient reflects how strongly the majority of fragment transitions simultaneously exhibit peak-like shapes at that location. This determines **which regions are candidate peaks** — it filters the chromatographic landscape down to regions where most fragments agree a peak exists.

A peptide may have multiple candidate peaks due to:

- **Isomers** eluting at different RTs
- **False positive signals** from co-eluting peptides with similar spectra
- **Peak splitting** from chromatographic artifacts

### Stage 2: Peak Selection (Pairwise Fragment Correlation)

When there are **multiple candidate peaks**, each is scored by computing the **mean pairwise Pearson correlation** of fragment intensities within that peak's boundaries. The peak with the highest mean pairwise correlation wins:

```text
For each candidate peak:
  Extract fragment intensities within [start_index, end_index]
  For each pair of fragments (i, j):
    Compute Pearson correlation r(fragment_i, fragment_j)
  coelution_score = mean of all pairwise correlations

Select peak with highest coelution_score
```

If only one candidate exists, it is used directly without scoring.

### Why Correlation, Not CWT Amplitude

The CWT consensus coefficient tells you **how peak-shaped** the consensus signal is at a location, but it is influenced by absolute signal intensity — a large interference peak could produce a higher CWT coefficient than a smaller true peak. The pairwise correlation score measures **how well fragments co-elute within the peak boundaries**, which is the defining characteristic of a true peptide signal:

- **True peak**: All fragments co-elute with high pairwise correlation
- **Interference peak**: Fragments from different peptides have uncorrelated intensities
- **CWT coefficient**: Could be high for either, if intensity is high enough

The CWT consensus median already filters out obvious interference (regions where only a minority of fragments have signal), so the candidates that survive to Stage 2 all have reasonable multi-transition support. Among these plausible candidates, pairwise correlation is the strongest discriminator.

### Minimum Correlation Threshold

After selecting the best peak, a minimum co-elution threshold is applied. In the calibration path (`batch.rs`), this threshold is `coelution_sum >= 0.5`. If no candidate meets this threshold, the precursor is rejected entirely — this prevents accepting peaks where fragments are present but not co-eluting.

### Summary

```text
Fragment XICs
  → CWT consensus (Mexican Hat wavelet + pointwise median)
  → Candidate peaks (ranked by CWT coefficient)
  → Score each candidate by mean pairwise fragment correlation
  → Select candidate with highest correlation
  → Apply minimum correlation threshold
  → Accepted peak (or reject if threshold not met)
```

---

## Testing

Test coverage:

### CWT tests (`crates/osprey-chromatography/src/cwt.rs`)

| Test | Description |
|------|-------------|
| `test_kernel_zero_mean` | Mexican Hat kernel sums to zero |
| `test_kernel_symmetric` | Kernel is symmetric around center |
| `test_kernel_positive_center_negative_tails` | Correct wavelet shape |
| `test_kernel_size` | Kernel length = 2 × radius + 1 |
| `test_convolve_same_length` | Output length equals input length |
| `test_convolve_delta_function` | Delta convolved with kernel reproduces kernel |
| `test_convolve_gaussian_response` | CWT of Gaussian is positive at center |
| `test_estimate_scale_known_peak` | Scale estimation from known FWHM |
| `test_estimate_scale_fallback` | Falls back to 4.0 for all-zero XICs |
| `test_consensus_single_gaussian_peak` | Finds peak with correct apex and boundaries |
| `test_consensus_interference_rejection` | Median suppresses single-fragment interference |
| `test_consensus_two_separated_peaks` | Resolves two peaks at different RTs |
| `test_consensus_degenerate_single_xic` | Returns empty for < 2 XICs |
| `test_consensus_degenerate_short_xic` | Returns empty for < 5 scans |
| `test_consensus_all_zero` | Returns empty for no signal |
| `test_consensus_peak_bounds_valid` | Boundaries are ordered and within array |
| `test_consensus_noise_robustness` | Finds correct peak despite noise |

### Legacy peak detection tests (`crates/osprey-chromatography/src/lib.rs`)

| Test | Description |
|------|-------------|
| `test_peak_detector` | Simple Gaussian peak detection |
| `test_find_best_peak` | Peak selection within RT tolerance |
| `test_fwhm_cap_symmetric_peak` | FWHM cap does not tighten symmetric peaks |
| `test_fwhm_cap_tailing_peak` | FWHM cap tightens trailing boundary |
| `test_fwhm_cap_preserves_valley` | Valley boundary preserved |
| `test_fwhm_cap_fallback` | Graceful fallback when FWHM cannot be computed |
| `test_asymmetric_half_widths` | Half-width computation for asymmetric peaks |

---

## Implementation

| File | Purpose |
|------|---------|
| `crates/osprey-chromatography/src/cwt.rs` | `detect_cwt_consensus_peaks()`, `mexican_hat_kernel()`, `estimate_cwt_scale()`, `convolve_same()` |
| `crates/osprey-chromatography/src/lib.rs` | `PeakDetector::detect()`, `detect_xic_peak()`, `detect_all_xic_peaks()`, `smooth_weighted_avg()`, `find_apex()`, `walk_boundary_left/right()`, `compute_asymmetric_half_widths()`, `compute_snr()` |
| `crates/osprey-core/src/types.rs` | `PeakBoundaries`, `PeakQuality`, `XICPeakBounds` structs |
| `crates/osprey-scoring/src/lib.rs` | `tukey_median_polish()`, `TukeyMedianPolishResult`, `compute_fragment_coelution()`, scoring features from peak shape |
| `crates/osprey/src/pipeline.rs` | Peak detection calls in coelution pipeline, Tukey median polish integration, boundary propagation to blib output |
