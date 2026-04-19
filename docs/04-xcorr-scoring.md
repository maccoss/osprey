# XCorr Scoring

XCorr (cross-correlation) is a spectral similarity score originally developed for Comet and SEQUEST. This document describes the exact implementation used in Osprey, which matches Comet's calculation.

## Overview

```text
XCorr Calculation Flow:

  Experimental Spectrum              Theoretical Spectrum
         │                                   │
         ▼                                   ▼
    ┌─────────┐                        ┌─────────┐
    │ Binning │                        │ Binning │
    │ + sqrt  │                        │ (unit   │
    └────┬────┘                        │  1.0)   │
         │                             └────┬────┘
         ▼                                  │
    ┌─────────┐                             │
    │Windowing│                             │
    │ (10 win,│                             │
    │ max=50) │                             │
    └────┬────┘                             │
         │                                  │
         ▼                                  │
    ┌─────────┐                             │
    │Flanking │                             │
    │Subtract │                             │
    │(off=75) │                             │
    └────┬────┘                             │
         │                                  │
         ▼                                  ▼
    ┌────────────────────────────────────────┐
    │  XCorr = Σ exp_preprocessed[frag_bins] │
    │        × 0.005                         │
    └────────────────────────────────────────┘
```

**Key insight**: The theoretical spectrum is NOT preprocessed with windowing. Comet simply looks up values from the preprocessed experimental spectrum at fragment bin positions and sums them.

## Binning

Both experimental and theoretical spectra are binned using Comet's BIN macro:

```text
BIN(mass) = (int)(mass / bin_width + (1 - offset))
          = (int)(mass / bin_width + 0.6)    // with offset = 0.4
```

### Parameters

| Parameter | Unit Resolution | HRAM |
|-----------|----------------|------|
| bin_width | 1.0005079 Da | 0.02 Da |
| offset | 0.4 | 0.0 |
| m/z range | 0-2000 | 0-2000 |
| num_bins | ~2000 | ~100,000 |

### Experimental Spectrum Binning

```rust
// For each peak in the spectrum
let bin = ((mz / bin_width) + 0.6) as usize;
if bin < num_bins {
    binned[bin] += intensity.sqrt();  // Apply sqrt transformation
}
```

### Theoretical Spectrum Binning

```rust
// For each fragment in the library entry
let bin = ((fragment_mz / bin_width) + 0.6) as usize;
if bin < num_bins {
    binned[bin] = 1.0;  // Unit intensity, NO sqrt, NO library intensity
}
```

## Experimental Spectrum Preprocessing

The experimental spectrum undergoes three preprocessing steps:

### Step 1: Square Root Transformation

Applied during binning (see above). This dampens intensity differences.

### Step 2: Windowing Normalization

Comet's `MakeCorrData` function divides the spectrum into 10 windows and normalizes each to max=50:

```rust
fn apply_windowing_normalization(spectrum: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; spectrum.len()];
    let num_windows = 10;
    let window_size = (spectrum.len() / num_windows) + 1;

    // Find global max for threshold (5% of base peak)
    let global_max = spectrum.iter().cloned().fold(0.0, f64::max);
    let threshold = global_max * 0.05;

    for window_idx in 0..num_windows {
        let start = window_idx * window_size;
        let end = ((window_idx + 1) * window_size).min(spectrum.len());

        // Find max in this window
        let window_max = spectrum[start..end].iter().cloned().fold(0.0, f64::max);

        // Normalize this window to 50.0
        if window_max > 0.0 {
            let norm_factor = 50.0 / window_max;
            for i in start..end {
                if spectrum[i] > threshold {
                    result[i] = spectrum[i] * norm_factor;
                }
            }
        }
    }

    result
}
```

**Why windowing?**

- Prevents high-intensity peaks from dominating the score
- Gives equal weight to different mass regions
- The 5% threshold removes noise peaks

### Step 3: Flanking Bin Subtraction (Fast XCorr)

Comet's fast XCorr preprocessing subtracts the local average to create an autocorrelation-like effect.

**Naive O(n × 151) implementation** (for reference):

```rust
fn apply_sliding_window_naive(spectrum: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; spectrum.len()];
    let offset: i64 = 75;  // Comet's iXcorrProcessingOffset
    let window_range = 2 * offset + 1;  // 151 bins total
    let norm_factor = 1.0 / (window_range - 1) as f64;  // Exclude center

    for i in 0..spectrum.len() {
        let mut sum = 0.0;
        for j in (i as i64 - offset)..=(i as i64 + offset) {
            if j >= 0 && j < spectrum.len() as i64 && j != i as i64 {
                sum += spectrum[j as usize];
            }
        }
        result[i] = spectrum[i] - sum * norm_factor;
    }
    result
}
```

**Optimized O(n) prefix-sum implementation** (used in Osprey):

The naive approach is O(n × 151) — fine for 2K unit-resolution bins but ~14M ops for 95K HRAM bins per spectrum. Osprey uses a prefix-sum approach:

```rust
fn apply_sliding_window(spectrum: &[f64]) -> Vec<f64> {
    let n = spectrum.len();
    let offset: usize = 75;
    let norm_factor = 1.0 / (2 * offset) as f64;  // Exclude center from 151 bins

    // Build prefix sum: prefix[i] = spectrum[0] + spectrum[1] + ... + spectrum[i-1]
    let mut prefix = vec![0.0; n + 1];
    for i in 0..n {
        prefix[i + 1] = prefix[i] + spectrum[i];
    }

    // For each bin, compute window sum using prefix sums
    let mut result = vec![0.0; n];
    for i in 0..n {
        let left = if i >= offset { i - offset } else { 0 };
        let right = (i + offset + 1).min(n);

        // Window sum = prefix[right] - prefix[left], then subtract center
        let window_sum = prefix[right] - prefix[left] - spectrum[i];
        result[i] = spectrum[i] - window_sum * norm_factor;
    }
    result
}
```

This reduces sliding window from O(n × 151) to O(n), critical for HRAM bins where n ≈ 100,000.

**Why flanking subtraction?**

- Enhances peaks that stand out from their local background
- Reduces contribution from broad, noisy regions
- Creates positive values where peaks are higher than surroundings
- Creates negative values where peaks are lower than surroundings

## Theoretical Spectrum (NO Preprocessing)

**Critical**: The theoretical spectrum does NOT undergo windowing or flanking subtraction. It is simply a vector of unit intensities (1.0) at fragment bin positions.

```rust
fn preprocess_library_for_xcorr(entry: &LibraryEntry) -> Vec<f64> {
    let mut binned = vec![0.0; num_bins];

    for fragment in &entry.fragments {
        let bin = ((fragment.mz / bin_width) + 0.6) as usize;
        if bin < num_bins {
            binned[bin] = 1.0;  // Unit intensity only
        }
    }

    // Return directly - NO windowing, NO flanking subtraction
    binned
}
```

**Why no preprocessing for theoretical?**

- Comet's XCorr is designed as a lookup: sum the preprocessed experimental values at fragment positions
- Applying windowing would multiply values by ~50, inflating scores
- The theoretical spectrum serves as a "selector" for which bins to sum

## XCorr Calculation

The final XCorr score is a simple dot product with scaling:

```rust
fn xcorr_from_preprocessed(
    experimental_preprocessed: &[f64],
    theoretical: &[f64]
) -> f64 {
    let raw: f64 = experimental_preprocessed
        .iter()
        .zip(theoretical.iter())
        .map(|(exp, theo)| exp * theo)
        .sum();

    raw * 0.005  // Comet's scaling factor
}
```

Since the theoretical spectrum has unit intensities (1.0) at fragment positions and zeros elsewhere, this is equivalent to:

```rust
// Comet's actual implementation (simplified)
let mut xcorr = 0.0;
for fragment in &entry.fragments {
    let bin = ((fragment.mz / bin_width) + 0.6) as usize;
    if bin < num_bins {
        xcorr += experimental_preprocessed[bin];
    }
}
xcorr *= 0.005;
```

### Scaling Factor (0.005)

The 0.005 scaling factor comes from Comet's normalization:

- Windowing normalizes to max=50 per window
- 0.005 = 1/(50 × 4) approximately normalizes scores to a reasonable range
- Typical good XCorr values: 0.5-5.0
- Typical noise XCorr values: < 0.3

## BLAS Vectorization

For efficiency, Osprey uses BLAS-accelerated dot products for XCorr scoring.

### Calibration XCorr (per-peptide scoring)

During calibration, each peptide is scored against all spectra in its RT window. The experimental spectrum is preprocessed once, and XCorr is computed via `ndarray::ArrayView1::dot()` which dispatches to BLAS `sdot` for contiguous f32 slices:

```rust
// Preprocess experimental spectrum ONCE
let exp_preprocessed = scorer.preprocess_spectrum_for_xcorr(&spectrum);

// Preprocess library entry (just binning, no windowing)
let lib_preprocessed = scorer.preprocess_library_for_xcorr(&entry);

// XCorr via BLAS sdot (dispatched by ndarray)
let xcorr = exp_view.dot(&lib_view) * 0.005;
```

### Feature Extraction XCorr (O(n_fragments) optimization)

During feature extraction, XCorr is computed using an O(n_fragments) shortcut. Since the theoretical spectrum has unit intensities at fragment bin positions, the dot product reduces to summing the preprocessed experimental values at those positions:

```rust
// O(n_fragments) instead of O(n_bins)
let mut xcorr = 0.0;
for fragment in &entry.fragments {
    let bin = bin_config.mz_to_bin(fragment.mz);
    if bin < preprocessed.len() {
        xcorr += preprocessed[bin];
    }
}
xcorr *= 0.005;
```

This avoids creating a dense theoretical vector and is faster when n_fragments << n_bins (typical: 20-50 fragments vs 2K-100K bins).

## Per-Window Preprocessing Optimization

During the coelution search (Phase 3), `run_search()` processes each DIA isolation window in parallel. Within a window, the same set of spectra is shared by hundreds of precursor candidates. Without optimization, each candidate would independently preprocess every spectrum in its RT range — redundantly repeating the same sqrt + windowing + flanking subtraction work.

### The Problem

A typical DIA window contains ~200 spectra and hundreds of candidate precursors. Naive per-candidate preprocessing would preprocess the same spectrum hundreds of times:

```text
Window [400, 403] → 200 spectra, 300 candidates
Without optimization:  300 × 200 = 60,000 preprocessing calls
With optimization:                  200 preprocessing calls
```

### The Solution

All spectra in a window are preprocessed **once** at the start of the per-window loop, before iterating over candidates:

```rust
let preprocessed_xcorr: Vec<Vec<f32>> = window_spectra
    .iter()
    .map(|s| scorer.preprocess_spectrum_for_xcorr(s))
    .collect();
```

This `Vec<Vec<f32>>` is indexed by position in the window's spectrum list. Each candidate's `FeatureComputeContext` receives a reference to this shared array along with a mapping from the candidate's local spectrum indices to the window-level indices:

```rust
struct FeatureComputeContext<'a> {
    // ... other fields ...
    preprocessed_xcorr: Option<&'a [Vec<f32>]>,  // pre-preprocessed for all window spectra
    cand_window_local: Option<&'a [usize]>,       // candidate spectra → window index mapping
}
```

Inside `compute_features_at_peak()`, XCorr computation becomes an O(n_fragments) lookup per spectrum rather than O(n_bins) preprocessing + O(n_fragments) scoring:

```rust
// Look up the pre-preprocessed vector for this spectrum
let preprocessed = &ctx.preprocessed_xcorr[window_idx];

// O(n_fragments) bin lookups — no preprocessing needed
let mut xcorr = 0.0;
for fragment in &entry.fragments {
    if let Some(bin) = bin_config.mz_to_bin(fragment.mz) {
        xcorr += preprocessed[bin] as f64;
    }
}
xcorr *= 0.005;
```

### Memory Impact

Each preprocessed vector has `n_bins` entries of 4 bytes (f32):

- **Unit resolution**: 2,001 bins × 4 bytes = ~8 KB per spectrum
- **HRAM**: 100,000 bins × 4 bytes = ~400 KB per spectrum

For a window with 200 spectra: ~80 MB (HRAM) or ~1.6 MB (unit resolution). This is allocated once per window and freed when the window loop moves on.

## Calibration: Always Unit Resolution Bins

Calibration XCorr always uses **unit resolution bins** (2001 bins, 1.0005 Da) regardless of whether the data is unit resolution or HRAM. This provides ~50x faster scoring for HRAM data while maintaining sufficient discriminating power for calibration.

HRAM-resolution bins (0.02 Da, 100K bins) are used during the main search phase (Phase 3) for fragment matching, not for calibration XCorr.

## Comparison: XCorr vs LibCosine

| Aspect | XCorr | LibCosine |
|--------|-------|-----------|
| **Preprocessing** | Binning + windowing + flanking | sqrt(intensity) only |
| **Normalization** | Window max=50 | L2 normalization |
| **Matching** | Bin position lookup | ppm tolerance matching |
| **Score range** | 0-10 typical | 0-1 (cosine similarity) |
| **Used in** | Calibration (all spectra) + feature extraction | Feature extraction only |
| **Library intensity** | Ignored (unit=1.0) | sqrt-transformed |

## Example

```text
Peptide: PEPTIDEK (charge 2+)
Fragments: y3(391.2), y4(504.3), y5(605.3), y6(702.4), y7(831.4)

Experimental spectrum at RT=25.3 min after preprocessing:
  Bin 391: +12.3  (peak above local average)
  Bin 504: +8.7   (peak above local average)
  Bin 605: +15.2  (strong peak)
  Bin 702: -2.1   (below local average - noise)
  Bin 831: +22.4  (strong peak)

Theoretical spectrum (unit intensities):
  Bin 391: 1.0
  Bin 504: 1.0
  Bin 605: 1.0
  Bin 702: 1.0
  Bin 831: 1.0

XCorr calculation:
  raw = 12.3 + 8.7 + 15.2 + (-2.1) + 22.4 = 56.5
  xcorr = 56.5 × 0.005 = 0.283

Note: The negative value at bin 702 reduces the score,
penalizing the missing/weak y6 ion.
```

## Common Mistakes

### 1. Applying windowing to theoretical spectrum (WRONG)

```rust
// WRONG - inflates scores by ~50x
let theoretical_windowed = apply_windowing_normalization(&theoretical);
```

### 2. Using library intensities instead of unit intensities (WRONG)

```rust
// WRONG - should use 1.0, not library intensity
binned[bin] = fragment.relative_intensity;
```

### 3. Starting m/z range at 200 instead of 0 (WRONG)

```rust
// WRONG - Comet starts at 0 m/z
let bin = ((mz - 200.0) / bin_width) as usize;

// CORRECT - use Comet's BIN macro
let bin = ((mz / bin_width) + 0.6) as usize;
```

### 4. Forgetting the 0.005 scaling factor (WRONG)

```rust
// WRONG - scores will be ~200x too large
return raw_dot_product;

// CORRECT
return raw_dot_product * 0.005;
```

## Implementation

Key files:

- `crates/osprey-scoring/src/lib.rs` - SpectralScorer with XCorr methods
- `crates/osprey/src/pipeline.rs` - Integration in calibration workflow

### API

```rust
use osprey_scoring::SpectralScorer;
use osprey_core::BinConfig;

// Create scorer — unit resolution (default, used for calibration)
let scorer = SpectralScorer::new();  // 2001 bins, 1.0005 Da

// Create scorer — HRAM (for feature extraction)
let scorer = SpectralScorer::with_bin_config(BinConfig::hram());  // 100K bins, 0.02 Da

// Single spectrum scoring
let score = scorer.xcorr(&spectrum, &library_entry);

// Batch preprocessing for BLAS sdot
let exp_preprocessed = scorer.preprocess_spectrum_for_xcorr(&spectrum);
let theo_preprocessed = scorer.preprocess_library_for_xcorr(&entry);
let xcorr = SpectralScorer::xcorr_from_preprocessed(&exp_preprocessed, &theo_preprocessed);
```

## References

- Comet source code: `CometPreprocess.cpp` (MakeCorrData, FastXcorrPreprocessing)
- Comet source code: `CometSearch.cpp` (XcorrScore function)
- Eng, J.K., et al. (2013). "Comet: An open-source MS/MS sequence database search tool"
- pyXcorrDIA: <https://github.com/maccoss/pyXcorrDIA>
