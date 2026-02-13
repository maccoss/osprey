# Peak Detection

After ridge regression, each peptide has a **coefficient time series** -- a sequence of (RT, coefficient) pairs across all MS2 spectra where that peptide had a non-zero regression coefficient. Peak detection identifies where the peptide elutes and determines the apex RT where all spectral scores will be computed.

## Overview

```
Coefficient
    |
1.0 |           .---.
    |          /     \
0.5 |         /       \
    |        /         \
0.05|......./.............\.......  ← min_height threshold
0.0 |------/---------------\------
    +-----|-----|-----|-----|------→ RT
         start apex  end
```

The peak detector uses a **threshold-crossing algorithm**: it scans through the coefficient time series and identifies contiguous regions where the coefficient stays above a minimum height. The highest point within each region is the apex.

## Algorithm

### Threshold-Crossing Peak Detection

The actual algorithm implemented in `PeakDetector::detect()`:

```
Input:  series[] = [(RT₁, coef₁), (RT₂, coef₂), ..., (RTₙ, coefₙ)]
Params: min_height = 0.05, min_width = 3 scans

Initialize:
  in_peak = false
  peaks = []

For each point (RT_i, coef_i) in series:

  If coef_i >= min_height:
    If NOT in_peak:
      # Entering a peak region
      in_peak = true
      peak_start = i
      apex_idx = i
      apex_value = coef_i
    Else if coef_i > apex_value:
      # Update apex within current peak
      apex_idx = i
      apex_value = coef_i

  Else if in_peak:
    # Leaving a peak region (coefficient dropped below threshold)
    peak_end = i - 1
    in_peak = false

    If (peak_end - peak_start + 1) >= min_width:
      area = sum(coef[peak_start..peak_end])
      peaks.append(Peak{
        start_rt = series[peak_start].RT,
        end_rt   = series[peak_end].RT,
        apex_rt  = series[apex_idx].RT,
        apex_coefficient = apex_value,
        integrated_area  = area
      })

# Handle peak at end of series
If in_peak AND (len - peak_start + 1) >= min_width:
  peaks.append(...)

Return peaks
```

### Key Properties

- **Threshold-based**: A point is "in a peak" if its coefficient >= `min_height`. This means the boundary is defined by where the signal crosses the threshold, not by local minima or inflection points.
- **Single pass**: O(n) scan through the time series.
- **No smoothing**: Operates directly on raw regression coefficients.
- **Multiple peaks**: Can detect multiple separate peaks for the same peptide (e.g., from isomers or false positives).

### Parameters

| Parameter | Default | Production | Description |
|-----------|---------|------------|-------------|
| `min_height` | 0.01 | 0.05 | Minimum coefficient value for a point to be part of a peak |
| `min_width` | 3 | 3 | Minimum number of scans (data points) in a valid peak |

The production value of `min_height = 0.05` is set in `score_run()` in `pipeline.rs`. This means the regression must assign at least 5% of the observed spectrum's signal to a peptide for it to be considered part of that peptide's chromatographic peak.

## Selecting the Best Peak

When multiple peaks are detected for the same peptide, the algorithm selects the best one:

```
find_best_peak(peaks, expected_rt, rt_tolerance):
  1. Filter to peaks where |apex_rt - expected_rt| <= rt_tolerance
  2. Among those, select the peak with the highest apex_coefficient
  3. Return None if no peaks are within tolerance
```

The expected RT comes from the LOESS calibration curve (`RTCalibration.predict(library_rt)`). This ensures that among multiple coefficient peaks, the one closest to the expected elution time with the strongest signal is chosen.

## Coefficient Time Series Construction

Before peak detection, the coefficient time series is built by aggregating regression results:

```
For each library entry (peptide):
  series = []
  For each regression result (one per spectrum):
    If this peptide has a non-zero coefficient:
      series.append((result.retention_time, coefficient))
  Sort series by RT
```

Points where the peptide was NOT a candidate (not in isolation window, not in RT tolerance) or had a zero coefficient are simply absent from the series. The series only contains positive evidence.

## Peak Boundaries

```rust
pub struct PeakBoundaries {
    pub start_rt: f64,           // First RT above threshold
    pub end_rt: f64,             // Last RT above threshold
    pub apex_rt: f64,            // RT of highest coefficient
    pub apex_coefficient: f64,   // Maximum coefficient value
    pub integrated_area: f64,    // Sum of all coefficients in peak
    pub peak_quality: PeakQuality,
}
```

### Integrated Area

The integrated area is the sum of all regression coefficients within the peak boundaries. This is used for:
- Peak ranking within a peptide (when multiple peaks exist)
- Quality control

**Important**: Do NOT use integrated area for quantification. Regression coefficients are relative within each spectrum and are not directly comparable across spectra. Skyline handles quantification using the extracted ion chromatogram within the peak boundaries provided by Osprey.

## Chromatographic Features from Peak Shape

After peak detection, chromatographic features are extracted from the coefficient time series for Mokapot scoring. These are computed by `FeatureExtractor` in `osprey-scoring`:

**Features sent to Mokapot PIN:**

| Feature | Field | Description |
|---------|-------|-------------|
| Peak apex | `peak_apex` | Maximum coefficient value |
| Peak area | `peak_area` | Integrated coefficient sum |
| Peak width | `peak_width` | FWHM from Tukey median polish elution profile (minutes) |
| Coefficient stability | `coefficient_stability` | CV of coefficients near apex |
| Signal to noise | `signal_to_noise` | Peak apex / noise estimate |
| XIC signal to noise | `xic_signal_to_noise` | S/N from best-correlated fragment XIC |

## Example

```
Peptide: PEPTIDEK (charge 2+)
Library RT: 25.3 min
Calibrated expected RT: 24.8 min

Coefficient time series (from ridge regression):
  RT 22.0: 0.00  (not in series - zero coefficient)
  RT 22.5: 0.00
  RT 23.0: 0.02  (below min_height, not in peak)
  RT 23.5: 0.08  ← peak start (first point >= 0.05)
  RT 24.0: 0.25
  RT 24.5: 0.58
  RT 25.0: 0.89  ← apex (highest coefficient)
  RT 25.5: 0.72
  RT 26.0: 0.35
  RT 26.5: 0.12
  RT 27.0: 0.03  ← peak ends here (drops below 0.05)
  RT 27.5: 0.00

Peak detected:
  start_rt: 23.5 min
  apex_rt:  25.0 min
  end_rt:   26.5 min
  apex_coefficient: 0.89
  integrated_area:  3.00
  width (scans):    7 (passes min_width=3)

All spectral scores computed at RT=25.0 (apex spectrum).
Delta RT = 25.0 - 24.8 = 0.2 min (observed - predicted).
```

## Multiple Peaks

A peptide may have multiple peaks due to:
- **Isomers** eluting at different RTs
- **False positive signals** from co-eluting peptides with similar spectra
- **Peak splitting** from chromatographic artifacts

Osprey selects the peak with the highest apex coefficient that is within RT tolerance of the expected RT. If no peak is within tolerance, the peptide is not reported.

## Peak Boundaries via Tukey Median Polish

After initial peak detection identifies the apex, Osprey refines peak boundaries using Tukey median polish on fragment XICs:

1. Extract XICs for top 6 library fragments
2. Run Tukey median polish on the 6×N XIC matrix (log space)
3. Column effects give a robust shared elution profile
4. Compute FWHM on the elution profile
5. Convert to boundaries: `σ = FWHM / 2.355`, `[apex - 1.96σ, apex + 1.96σ]`

The Tukey median polish approach is more robust than using any single fragment XIC because the median across 6 fragments suppresses interference on individual transitions.

See [docs/README.md](README.md) for full details on the Tukey median polish algorithm.

## Implementation

| File | Purpose |
|------|---------|
| `crates/osprey-chromatography/src/lib.rs` | `PeakDetector::detect()`, `find_best_peak()` |
| `crates/osprey-core/src/types.rs` | `PeakBoundaries`, `PeakQuality` structs |
| `crates/osprey-scoring/src/lib.rs` | `FeatureExtractor` (chromatographic features from peak shape) |
| `crates/osprey/src/pipeline.rs` | `score_run()` (calls peak detection, sets `min_height=0.05`) |
