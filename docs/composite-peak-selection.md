# Composite Peak Selection for Better Apex RT

## Problem

When multiple chromatographic peaks are detected for the same precursor, selecting the correct apex RT is critical. If the wrong RT is picked, ALL downstream features (spectral scores, mass accuracy, explained intensity, coelution) are scored at the wrong location.

### Current Flow

```
Library + Spectra
      |
[1] CANDIDATE SELECTION (mz_window + RT_tolerance)
      |
[2] FRAGMENT XIC EXTRACTION
      |
[3] PEAK DETECTION on reference XIC
    |-- Pairwise fragment correlations
    |-- Select reference XIC (best-correlated fragment)
    |-- detect_xic_peak() → apex, boundaries
      |
[4] FEATURE EXTRACTION (all 45 features scored at locked apex_rt)
      |
[5] FDR (Percolator/Mokapot)
```

### How DIA-NN Does It

DIA-NN uses a much more careful approach:
1. Detects multiple peak candidates via fragment XIC inter-correlation
2. For each peak, computes ~40 features (spectral, chromatographic, MS1, RT)
3. Uses a learned linear classifier (Percolator-trained weights) to pick the best peak
4. Fragment inter-correlation is the key signal -- peaks where fragments co-elute strongly are preferred
5. Apex must be at 99% of local max (rejects shoulders/noise spikes)

## Proposed Approach: Lightweight Composite Peak Scoring

For each detected peak candidate, compute a cheap composite score from orthogonal signals, then select the peak with the highest composite score. This happens BEFORE the expensive full feature extraction.

### The Signals

| Signal | Weight | Description |
|--------|--------|-------------|
| Fragment coelution | 0.35 | Per-peak pairwise correlation between fragment XICs. Strongest discriminator: interference peaks have poor fragment coelution. |
| RT agreement | 0.25 | `1.0 - |peak_apex_rt - expected_rt| / rt_tolerance`, clamped to [0,1]. Peaks near calibrated expected RT preferred. Essentially free. |
| Peak shape | 0.20 | Combine: (a) n_scans in peak (log-scaled), (b) signal-to-noise, (c) relative area. Real peaks vs noise spikes. |
| Peak intensity | 0.20 | `log(apex_intensity)` normalized across peaks. The simplest signal, but now one component rather than sole determinant. |

Each component is normalized to [0, 1] across the peaks for this precursor, then combined with fixed weights.

### Cost Analysis

- Fragment XICs are extracted once for the full series, then sliced per peak -- same binary search cost as existing scoring
- Per-peak Pearson correlation on 5-15 points: negligible
- RT score and peak shape: arithmetic only
- Most peptides have 1 peak (no extra work); ~10-20% have 2-3 peaks
- **Expected overhead: ~5% on scoring phase, ~1-3% on total pipeline**

## Implementation Plan

### 1. `crates/osprey-scoring/src/lib.rs` -- Add `score_peak_candidates()`

New public function:

```rust
pub fn score_peak_candidates(
    peaks: &[XICPeakBounds],
    fragment_xics: &[Vec<(f64, f64)>],
    expected_rt: f64,
    rt_tolerance: f64,
) -> Vec<(usize, f64)>  // (peak_index, composite_score)
```

Implementation:
- For each peak, slice fragment XICs to `[start_rt, end_rt]`
- Compute lightweight pairwise fragment coelution on the sub-series
- Compute RT score: `(1.0 - |apex_rt - expected_rt| / rt_tolerance).clamp(0.0, 1.0)`
- Compute shape score from: `n_scans.ln()`, signal-to-noise, `area / max_area`
- Compute intensity score: `apex_intensity.ln()` normalized across peaks
- Normalize each component to [0, 1] across peaks, combine with fixed weights

### 2. `crates/osprey/src/pipeline.rs` -- Use composite scoring in coelution search

In the per-peptide closure, when multiple peaks are detected from `detect_xic_peak()`:

- Score each peak candidate using `score_peak_candidates()`
- Select the peak with highest composite score
- Pass the selected peak's boundaries to feature extraction

### 3. Logging -- Peak selection diagnostics

In the scoring phase, after scoring all peptides, log:
- How many precursors had multiple peak candidates
- How often the non-highest-intensity peak was selected
- Distribution of peak counts

This directly measures whether the composite scoring is changing decisions.

## Existing Code to Reuse

| Function | Location | Purpose |
|----------|----------|---------|
| `pearson_correlation_raw()` | `osprey-scoring/src/lib.rs` | Raw Pearson correlation between two f64 slices |
| `extract_fragment_xics()` | `osprey-scoring/src/lib.rs` | Extract XICs for top N fragments from spectra |
| `detect_xic_peak()` | `osprey-chromatography/src/lib.rs` | Returns `XICPeakBounds` with apex, boundaries, area, S/N |

## Verification

1. `cargo fmt && cargo clippy --all-targets --all-features -- -D warnings && cargo test` -- all pass
2. `cargo build --release`
3. Run on test data and verify:
   - Log shows peak selection statistics
   - Compare FDR results: expect same or better precursor count at 1% FDR
4. Manual spot-check: find cases where non-max-intensity peak was selected, verify the selected peak is more plausible (closer to expected RT, better fragment coelution)

## Future Extension: Learned Weights

The fixed weights (0.35, 0.25, 0.20, 0.20) can later be replaced with learned weights using the same LDA pattern from `calibration_ml.rs`. After Percolator/Mokapot identifies true positives at stringent FDR, peak selection weights could be optimized to maximize agreement with confirmed results. Not needed for initial implementation.
