# XCorr Scoring

XCorr (cross-correlation) is a spectral similarity score originally developed for Comet and SEQUEST. This document describes the exact implementation used in Osprey, which matches Comet's calculation.

## Overview

```
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

```
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

Comet's fast XCorr preprocessing subtracts the local average to create an autocorrelation-like effect:

```rust
fn apply_sliding_window(spectrum: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; spectrum.len()];
    let offset: i64 = 75;  // Comet's iXcorrProcessingOffset
    let window_range = 2 * offset + 1;  // 151 bins total
    let norm_factor = 1.0 / (window_range - 1) as f64;  // Exclude center

    for i in 0..spectrum.len() {
        let mut sum = 0.0;

        // Sum values in window, excluding center
        for j in (i as i64 - offset)..=(i as i64 + offset) {
            if j >= 0 && j < spectrum.len() as i64 && j != i as i64 {
                sum += spectrum[j as usize];
            }
        }

        // Subtract local average from center value
        result[i] = spectrum[i] - sum * norm_factor;
    }

    result
}
```

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

For efficiency, Osprey preprocesses all spectra and library entries once per isolation window, then uses matrix multiplication:

```rust
// Preprocess all experimental spectra ONCE
let experimental_matrix: Vec<Vec<f64>> = spectra
    .iter()
    .map(|s| preprocess_spectrum_for_xcorr(s))
    .collect();

// Preprocess all library entries ONCE (just binning, no windowing)
let theoretical_matrix: Vec<Vec<f64>> = entries
    .iter()
    .map(|e| preprocess_library_for_xcorr(e))
    .collect();

// Matrix multiply for all-vs-all scoring
// xcorr_matrix[i][j] = dot(theoretical[i], experimental[j]) * 0.005
```

This allows scoring thousands of peptide-spectrum pairs with a single BLAS call.

## Comparison: XCorr vs LibCosine

| Aspect | XCorr | LibCosine |
|--------|-------|-----------|
| **Preprocessing** | Binning + windowing + flanking | sqrt(intensity) only |
| **Normalization** | Window max=50 | L2 normalization |
| **Matching** | Bin position lookup | ppm tolerance matching |
| **Score range** | 0-10 typical | 0-1 (cosine similarity) |
| **Best for** | RT selection (all spectra) | Fragment scoring (best RT) |
| **Library intensity** | Ignored (unit=1.0) | sqrt-transformed |

## Example

```
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

// Create scorer (unit resolution)
let scorer = SpectralScorer::unit_resolution();

// Single spectrum scoring
let score = scorer.xcorr(&spectrum, &library_entry);

// Batch preprocessing for BLAS
let exp_preprocessed = scorer.preprocess_spectrum_for_xcorr(&spectrum);
let theo_preprocessed = scorer.preprocess_library_for_xcorr(&entry);
let xcorr = SpectralScorer::xcorr_from_preprocessed(&exp_preprocessed, &theo_preprocessed);
```

## References

- Comet source code: `CometPreprocess.cpp` (MakeCorrData, FastXcorrPreprocessing)
- Comet source code: `CometSearch.cpp` (XcorrScore function)
- Eng, J.K., et al. (2013). "Comet: An open-source MS/MS sequence database search tool"
- pyXcorrDIA: https://github.com/maccoss/pyXcorrDIA
