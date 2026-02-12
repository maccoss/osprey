# Composite Peak Selection for Better Apex RT

## Problem

Osprey currently picks the apex RT as simply "the spectrum with the highest ridge regression coefficient." This is fragile because:

- A high coefficient doesn't mean the peptide is really there — interference from co-eluting peptides sharing fragment ions can inflate coefficients at the wrong RT
- If the wrong RT is picked, ALL downstream features (spectral scores, mass accuracy, explained intensity, coelution, z-scores) are scored at the wrong location
- Multiple peaks ARE detected by `PeakDetector::detect()` but only the max coefficient is used; other peaks are discarded

### Current Flow

```
Library + Spectra
      |
[1] CANDIDATE SELECTION (mz_window + RT_tolerance + top-2-of-6 fragment filter)
      |
[2] RIDGE REGRESSION (solve for coefficients in each spectrum)
      |
[3] APEX SELECTION  <-- THE PROBLEM
    |-- Aggregate (RT, coefficient) pairs
    |-- Detect peaks (PeakDetector)
    |-- SELECT RT with max coefficient  <-- Too simple!
      |
[4] FEATURE EXTRACTION (all 28+ features scored at locked apex_rt)
      |
[5] FDR (Mokapot)
```

### How DIA-NN Does It

DIA-NN uses a much more careful approach:
1. Detects multiple peak candidates via fragment XIC inter-correlation
2. For each peak, computes ~40 features (spectral, chromatographic, MS1, RT)
3. Uses a learned linear classifier (Percolator-trained weights) to pick the best peak
4. Fragment inter-correlation is the key signal — peaks where fragments co-elute strongly are preferred
5. Apex must be at 99% of local max (rejects shoulders/noise spikes)

## Proposed Approach: Lightweight Composite Peak Scoring

For each detected peak, compute a cheap composite score from 4 orthogonal signals, then select the peak with the highest composite score. This happens BEFORE the expensive full feature extraction.

### The 4 Signals

| Signal | Weight | Description |
|--------|--------|-------------|
| Fragment coelution | 0.35 | Per-peak Pearson correlation between top-6 fragment XICs and coefficient sub-series. Strongest discriminator: interference peaks have high coefficients but poor fragment coelution. |
| RT agreement | 0.25 | `1.0 - |peak_apex_rt - expected_rt| / rt_tolerance`, clamped to [0,1]. Peaks near calibrated expected RT preferred. Essentially free. |
| Peak shape | 0.20 | Combine: (a) n_scans in peak (log-scaled), (b) coefficient SNR (apex/median), (c) relative area. Real peaks vs noise spikes. |
| Coefficient amplitude | 0.20 | `log(apex_coefficient)` normalized across peaks. The current signal, but now one component rather than sole determinant. |

Each component is normalized to [0, 1] across the peaks for this precursor, then combined with fixed weights.

### Peak Boundaries and Background

Pass both the full coefficient series AND the selected peak's boundaries to feature extraction:
- **Peak-focused features** (FWHM, peak_area, peak_width, peak_symmetry) use only data within the peak boundaries — no confusion from neighboring peaks
- **Background/noise features** (SNR, peak_prominence, coefficient_stability) use data outside the peak boundaries for noise estimation
- The `extract_with_deconvolution` function gets the selected peak's `start_rt`/`end_rt` as an additional parameter

### Cost Analysis

- Fragment XICs are extracted once for the full series, then sliced per peak — same binary search cost as existing `compute_fragment_coelution` call
- Per-peak Pearson correlation on 5-15 points: negligible
- RT score and peak shape: arithmetic only
- Most peptides have 1 peak (no extra work); ~10-20% have 2-3 peaks
- **Expected overhead: ~5% on scoring phase, ~1-3% on total pipeline**

## Implementation Plan

### 1. `crates/osprey-scoring/src/lib.rs` — Add `score_peak_candidates()`

New public function:

```rust
pub fn score_peak_candidates(
    peaks: &[PeakBoundaries],
    coefficient_series: &[(f64, f64)],
    spectra: &[&Spectrum],
    library_entry: &LibraryEntry,
    expected_rt: f64,
    rt_tolerance: f64,
    tolerance_da: f64,
    tolerance_ppm: f64,
) -> Vec<(usize, f64)>  // (peak_index, composite_score)
```

Implementation:
- For each peak, slice `coefficient_series` to `[start_rt, end_rt]`
- Filter `spectra` to same RT range
- Compute lightweight fragment coelution on the sub-series (reuse XIC extraction + Pearson correlation from `compute_fragment_coelution`, restricted to peak window)
- Compute RT score: `(1.0 - |apex_rt - expected_rt| / rt_tolerance).clamp(0.0, 1.0)`
- Compute shape score from: `n_scans.ln()`, `apex_coef / median_coef`, `area / max_area`
- Compute coefficient score: `apex_coefficient.ln()` normalized across peaks
- Normalize each component to [0, 1] across peaks, combine with fixed weights

**Optimization**: Extract fragment XICs once for the full spectra set, then slice per peak by RT range. Avoids redundant m/z binary searches.

### 2. `crates/osprey/src/pipeline.rs` — Use composite scoring in `score_run()`

In the per-peptide closure (~line 1965-1976), replace simple `max_by(coefficient)`:

**Before:**
```rust
let peaks = peak_detector.detect(&rt_coef_pairs);
if peaks.is_empty() { return None; }
let apex_rt = sorted_data.iter().max_by(|a, b| a.1.total_cmp(&b.1))...;
```

**After:**
```rust
let peaks = peak_detector.detect(&rt_coef_pairs);
if peaks.is_empty() { return None; }

// Move peak_region_spectra computation here (before apex selection)
let peak_region_spectra = ...;  // existing code, moved earlier

// Score each peak candidate
let peak_scores = osprey_scoring::score_peak_candidates(
    &peaks, &rt_coef_pairs, &peak_region_spectra,
    entry, entry.retention_time, rt_tolerance,
    spectral_scorer.tolerance_da(), spectral_scorer.tolerance_ppm(),
);

// Select best peak
let best_peak_idx = peak_scores.iter()
    .max_by(|a, b| a.1.total_cmp(&b.1))
    .map(|(idx, _)| *idx)
    .unwrap_or(0);
let best_peak = &peaks[best_peak_idx];
let apex_rt = best_peak.apex_rt;
```

Pass both `rt_coef_pairs` (full series) and `best_peak` boundaries to `extract_with_deconvolution` so it can focus peak features on the selected peak while using surrounding data for background estimation.

Update `apex_scan_number` to use `best_peak.apex_rt` rather than global coefficient max.

### 3. `crates/osprey-core/src/types.rs` — Add `n_peaks_detected` to FeatureSet

```rust
/// Number of chromatographic peaks detected for this precursor
pub n_peaks_detected: u32,
```

Precursors with many detected peaks are more likely to have interference — useful as a mokapot feature.

### 4. `crates/osprey-fdr/src/mokapot.rs` — Add `n_peaks_detected` to PIN

Add to header and format_features (29 features total).

### 5. Logging — Peak selection diagnostics

In `score_run()`, after scoring all peptides, log:
- How many precursors had multiple peaks
- How often the non-highest-coefficient peak was selected
- Distribution of peak counts

This directly measures whether the composite scoring is changing decisions.

## Existing Code to Reuse

| Function | Location | Purpose |
|----------|----------|---------|
| `compute_fragment_coelution()` | `osprey-scoring/src/lib.rs:266` | Fragment XIC extraction + Pearson correlation (reuse XIC extraction logic) |
| `pearson_correlation_raw()` | `osprey-scoring/src/lib.rs` | Raw Pearson correlation between two f64 slices |
| `extract_fragment_xics()` | `osprey-scoring/src/lib.rs:444` | Extract XICs for top N fragments from spectra |
| `PeakDetector::detect()` | `osprey-chromatography/src/lib.rs:88` | Returns `Vec<PeakBoundaries>` with start_rt, end_rt, apex_rt, apex_coefficient, integrated_area |
| `PeakDetector::find_best_peak()` | `osprey-chromatography/src/lib.rs:157` | Current best-peak selection (to be replaced) |
| `compute_fwhm_interpolated()` | `osprey-scoring/src/lib.rs:659` | FWHM computation from coefficient series |

## Verification

1. `cargo fmt && cargo clippy --all-targets --all-features -- -D warnings && cargo test` — all pass
2. `cargo build --release`
3. Run on test data and verify:
   - PIN file has 29 features (was 28) + 5 charge columns
   - `n_peaks_detected` values are positive integers (typically 1-5)
   - Log shows peak selection statistics
   - Compare FDR results: expect same or better precursor count at 1% FDR
4. Manual spot-check: find cases where non-max-coefficient peak was selected, verify the selected peak is more plausible (closer to expected RT, better fragment coelution)

## Future Extension: Learned Weights

The fixed weights (0.35, 0.25, 0.20, 0.20) can later be replaced with learned weights using the same LDA pattern from `calibration_ml.rs`. After mokapot identifies true positives at stringent FDR, peak selection weights could be optimized to maximize agreement with mokapot-confirmed results. Not needed for initial implementation.