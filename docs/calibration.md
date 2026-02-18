# Calibration

Osprey performs joint calibration of retention time, MS1 (precursor) m/z, and MS2 (fragment) m/z from a single set of high-confidence peptide matches. This enables tighter tolerance windows for the main search and more accurate quantification.

## Overview

All three calibrations share the same calibration discovery phase: a co-elution search scored with LDA-based machine learning, followed by paired target-decoy competition at 1% FDR. The confident matches provide (library_RT, measured_RT) pairs for RT calibration and mass error measurements for m/z calibration.

```
Calibration Discovery Workflow
──────────────────────────────

1. Sample library entries (up to 100K targets + paired decoys)
2. Co-elution scoring: fragment XIC correlation, spectral similarity
3. LDA scoring: learn optimal feature weights via 3-fold cross-validation
4. Paired target-decoy competition → q-values
5. Filter to targets at ≤1% FDR with S/N ≥ 5.0
6. From confident matches, simultaneously extract:
   ├── (library_RT, measured_RT) pairs  → LOESS RT calibration
   ├── MS1 errors (M+0 isotope peak)   → precursor m/z calibration
   └── MS2 errors (matched fragments)  → fragment m/z calibration
7. Save calibration JSON to disk for caching
```

## Calibration Discovery

### Library Sampling

For large libraries (millions of entries), scoring all entries is unnecessary. Osprey samples a subset with a retry loop:

- **Attempt 1**: Sample `calibration_sample_size` targets (default: 100,000) plus their paired decoys
- **Attempt 2**: Expand by `retry_factor` (default: 2x)
- **Attempt 3**: Use ALL library entries (guaranteed fallback)

Matches accumulate across attempts — FDR runs on the combined set. If the library has fewer targets than the sample size, all entries are used and retries are skipped.

### Co-Elution Scoring

Each sampled peptide is scored using fragment XIC co-elution, not simple XCorr:

1. Map library RT to expected measured RT (linear mapping if RT scales differ, identity if similar)
2. Set initial RT tolerance (20% of gradient range if scales similar, 50% if linear mapping needed)
3. For each library entry, extract fragment XICs from spectra within tolerance
4. Detect candidate peaks using **CWT consensus peak detection** (Mexican Hat wavelet convolution of each fragment XIC, pointwise median across transitions — see [Peak Detection](peak-detection.md))
5. For each candidate peak, score pairwise fragment correlation within its boundaries
6. Select the peak with the highest co-elution score
7. At apex, measure spectral similarity (libcosine), fragment matches, and signal-to-noise

If CWT consensus finds no peaks, falls back to SG-smoothed peak detection on the reference XIC (`detect_all_xic_peaks`).

### LDA Scoring

A Linear Discriminant Analysis model learns optimal feature weights from target-decoy pairs using Percolator-style iterative training with 3-fold cross-validation.

**4 features** (normalized to 0–1 range):

| Feature | Normalization | Description |
|---------|---------------|-------------|
| `correlation` | `score / 6.0` | XIC co-elution sum across fragments |
| `libcosine` | Already 0–1 | Spectral similarity at apex |
| `top6_matched` | `count / 6.0` | Number of top-6 library fragments matched at apex |
| `signal_to_noise` | `ln(S/N + 1) / 5.0` | Peak signal-to-noise ratio (log-transformed) |

Hyperscore and isotope cosine are intentionally excluded from the LDA — hyperscore dominates target-decoy discrimination but doesn't correlate with RT quality, and isotope cosine shows negative correlation.

**Training procedure:**

1. Evaluate each single feature independently, select the best as baseline
2. For up to 3 iterations:
   - 3-fold stratified CV (grouped by peptide sequence to keep charge states together)
   - In each fold: select targets at q ≤ 1% from the train set, train LDA on those targets + all decoys
   - Average weights across folds (consensus), clip negative weights to zero, renormalize
   - Score all data with consensus weights
3. Track best iteration — revert if LDA degrades, stop early after 2 consecutive non-improvements

### Target-Decoy Competition

Each target competes against its paired decoy (linked by `base_id = entry_id & 0x7FFFFFFF`):
- Higher discriminant score wins
- Ties go to the decoy (conservative for FDR estimation)
- Only winners enter the FDR walk
- Q-values computed using the conservative `(decoys + 1) / targets` formula

See [FDR Control](fdr-control.md) for full details on the target-decoy strategy.

### Quality Filtering

From the FDR-filtered winners, an additional quality filter removes matches with signal-to-noise < 5.0. This prevents peptides with good spectral scores but poor peak quality from corrupting RT calibration.

## RT Calibration

### LOESS Fitting

RT calibration uses LOESS (Locally Estimated Scatterplot Smoothing) to fit (library_RT, measured_RT) pairs from confident matches.

**Parameters:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bandwidth` | 0.3 | Fraction of data used for each local fit (30%) |
| `degree` | 1 | Local polynomial degree (linear) |
| `min_points` | 20 | Minimum calibration points required |
| `robustness_iter` | 2 | Bisquare robustness iterations |

**Algorithm:**

1. Sort calibration points by library RT
2. For each point, find k nearest neighbors (k = bandwidth × n)
3. Weight neighbors by tricube kernel: `w(u) = (1 - u³)³` for normalized distance u
4. Fit weighted linear regression (2×2 system)
5. Compute residuals
6. Robustness iterations: downweight outliers using bisquare function
   - Scale: `s = 6 × MAD` (median absolute deviation)
   - Weight: `(1 - (r/s)²)²` if |r/s| < 1, else 0
7. Final refit with robustness weights

### Prediction and Interpolation

For prediction at a query library RT:
- **Within range**: Binary search to find bracketing calibration points, linear interpolation between fitted values
- **Outside range**: Linear extrapolation using the nearest non-duplicate pair

**Duplicate library RT handling**: When two calibration points share the same library_rt (|x1 - x0| < 1e-12), the fitted values are averaged instead of interpolated. This prevents division-by-zero that would produce NaN values.

### Local RT Tolerance

Rather than a single global tolerance, Osprey computes a **local** RT tolerance that varies across the gradient:

1. At each calibration point, store the absolute residual |measured - fitted|
2. For a query RT, interpolate absolute residuals from neighboring calibration points
3. Smooth over ±2 neighbors to reduce noise from individual outliers
4. Tolerance = `max(smoothed_residual × factor, min_tolerance)`

Defaults: `factor = 3.0`, `min_tolerance = 0.1 min` (6 seconds)

This means well-calibrated regions of the gradient get tighter RT windows while regions with larger calibration residuals get wider windows.

### RT Calibration Statistics

After fitting, Osprey reports:

```
RT calibration: 500 peptides at 1% FDR
RT calibration fit: R²=0.9987, residual_SD=0.32 min
Calibrated RT tolerance: 0.96 min (3× residual SD)
```

- **n_points**: Number of peptides passing FDR and quality filters
- **R²**: Coefficient of determination (goodness of fit)
- **residual_sd**: Standard deviation of (measured - predicted) residuals
- **Global tolerance**: `max(residual_sd × rt_tolerance_factor, min_rt_tolerance)` used as fallback

## Mass Calibration

### Unit-Aware Calibration

Mass calibration is unit-aware, supporting both instrument types:
- **HRAM instruments** (e.g., Orbitrap, Astral): Errors measured in ppm
- **Unit resolution instruments** (e.g., Stellar): Errors measured in Th (Thomson)

Both MS1 and MS2 errors use the same unit within a run.

### MS1 (Precursor) Calibration

MS1 errors are collected from the monoisotopic (M+0) peak in the nearest MS1 spectrum:

```
For each confident match:
  1. Find nearest MS1 spectrum to the peak apex RT
  2. Calculate expected isotope m/z values:
     - Isotope spacing = 1.002868 Da / charge
     - Expected: M-1, M+0, M+1, M+2, M+3
  3. Extract M+0 peak within precursor tolerance
  4. Calculate error: ((observed - theoretical) / theoretical) × 10⁶  (ppm)
```

If no MS1 spectra are available, the isolation window center is used as a fallback for the observed precursor m/z.

### MS2 (Fragment) Calibration

MS2 errors are collected from all matched fragments at the peak apex spectrum:

```
For each confident match:
  1. At the best apex spectrum, match library fragments to observed peaks
  2. For each match, calculate error in the configured unit (ppm or Th)
  3. Add all individual fragment errors to the MS2 error collection
```

### Error Statistics

From the collected errors, Osprey calculates:
- **Mean**: Systematic mass offset (used for correction)
- **Median**: Robust estimate of systematic error
- **SD**: Random error component (sample standard deviation, n-1 denominator)
- **Adjusted tolerance**: |mean| + 3 × SD

Histograms are generated with bin widths of 1.0 ppm (HRAM) or 0.01 Th (unit resolution).

```
MS1 calibration: mean=-2.50 ppm, median=-2.40 ppm, SD=0.80 ppm (n=500)
MS2 calibration: mean=1.20 ppm, median=1.10 ppm, SD=1.00 ppm (n=3000)
```

## Applying Calibration in the Main Search

### MS2 Spectrum Correction

After calibration, systematic m/z offset is removed from all MS2 spectra before the main search:

- **PPM correction**: `corrected_mz = observed_mz - observed_mz × mean / 10⁶`
- **Th correction**: `corrected_mz = observed_mz - mean`

### Calibrated Fragment Tolerance

After correcting the systematic offset, the fragment matching tolerance is tightened:

```
calibrated_tolerance = max(3 × SD, minimum_floor)
```

| Instrument Type | Minimum Floor |
|-----------------|---------------|
| HRAM (ppm)      | 1.0 ppm       |
| Unit resolution (Th) | 0.05 Th  |

If calibration is not available, the base fragment tolerance from configuration is used.

### RT Candidate Selection

During the main search, candidates are selected using the calibrated RT model:

1. For each library entry: `expected_rt = LOESS.predict(library_rt)`
2. Compute local tolerance: `local_tolerance(library_rt, factor=3.0, min=0.1 min)`
3. Filter candidates by: `|observed_rt - expected_rt| ≤ local_tolerance`

## Calibration Caching

Calibration results are saved to a JSON file alongside the input mzML for reuse:

```
sample.mzML → sample.osprey_calibration.json
```

On subsequent runs, Osprey checks for a cached calibration file and reuses it if valid, skipping the calibration discovery phase.

### JSON Format

```json
{
  "metadata": {
    "num_confident_peptides": 500,
    "num_sampled_precursors": 50000,
    "calibration_successful": true,
    "timestamp": "2025-01-15T10:30:00Z",
    "isolation_scheme": {
      "num_windows": 60,
      "mz_min": 380.0,
      "mz_max": 980.0,
      "typical_width": 10.0,
      "uniform_width": true,
      "windows": [[385.0, 10.0], ...]
    }
  },
  "ms1_calibration": {
    "mean": -2.5,
    "median": -2.4,
    "sd": 0.8,
    "count": 500,
    "unit": "ppm",
    "adjusted_tolerance": 4.9,
    "window_halfwidth_multiplier": 3.0,
    "calibrated": true
  },
  "ms2_calibration": {
    "mean": 1.2,
    "median": 1.1,
    "sd": 1.0,
    "count": 3000,
    "unit": "ppm",
    "adjusted_tolerance": 4.2,
    "window_halfwidth_multiplier": 3.0,
    "calibrated": true
  },
  "rt_calibration": {
    "method": "LOESS",
    "residual_sd": 0.32,
    "n_points": 500,
    "r_squared": 0.9987,
    "model_params": {
      "library_rts": [1.0, 2.0, ...],
      "fitted_rts": [1.1, 2.2, ...],
      "abs_residuals": [0.05, 0.08, ...]
    }
  }
}
```

The `model_params` field stores the fitted LOESS curve as (library_rt, fitted_rt, abs_residual) triples, enabling reconstruction of the calibration via interpolation without re-fitting.

## Configuration

```yaml
rt_calibration:
  enabled: true                      # Enable/disable RT calibration
  loess_bandwidth: 0.3               # Fraction of data for local LOESS fits
  min_calibration_points: 200        # Target minimum calibration points
  rt_tolerance_factor: 3.0           # Multiplier for residual SD
  fallback_rt_tolerance: 2.0         # Used if calibration fails (minutes)
  min_rt_tolerance: 0.1              # Minimum RT tolerance floor (minutes, 6 sec)
  calibration_sample_size: 100000    # Target library entries per attempt
  calibration_retry_factor: 2.0      # Sample size multiplier for retry

# Mass tolerance settings
precursor_tolerance_ppm: 10.0        # For MS1 isotope extraction
resolution_mode:
  HRAM:
    tolerance_ppm: 10.0              # Base fragment tolerance (before calibration)
```

## Fallback Behavior

If calibration fails (insufficient confident matches):

```
RT calibration failed: Insufficient calibration points: 23 < 50 absolute minimum
Using fallback tolerance: 2.0 min
```

The absolute minimum for LOESS fitting is 50 points (hardcoded), regardless of the `min_calibration_points` config setting. If fewer than 50 points pass FDR and quality filtering after all retry attempts, calibration fails and the fallback tolerance is used with library RTs directly.

## Implementation

### Key Files

| File | Description |
|------|-------------|
| `crates/osprey-chromatography/src/calibration/mod.rs` | Core types: `CalibrationParams`, `MzCalibration`, `RTCalibrationParams` |
| `crates/osprey-chromatography/src/calibration/rt.rs` | LOESS fitting, `RTCalibration`, local tolerance |
| `crates/osprey-chromatography/src/calibration/mass.rs` | m/z calibration: `MzQCData`, error statistics, spectrum correction |
| `crates/osprey-chromatography/src/calibration/io.rs` | JSON serialization and disk caching |
| `crates/osprey-scoring/src/calibration_ml.rs` | LDA training with cross-validation |
| `crates/osprey-scoring/src/batch.rs` | Co-elution scoring, `CalibrationMatch` extraction |
| `crates/osprey/src/pipeline.rs` | Calibration discovery workflow and pipeline integration |
| `crates/osprey-core/src/config.rs` | `RTCalibrationConfig` with defaults |

## Troubleshooting

| Symptom | Possible Cause | Solution |
|---------|----------------|----------|
| Low confident count | Poor library match or mismatched RT scales | Check library source and RT units |
| High RT residual SD | Non-linear RT shift or gradient issues | Check LC gradient; LOESS should handle nonlinearity |
| Low R² | Outliers or multiple RT modes | Investigate gradient problems |
| Calibration fails | Too few confident matches | Widen initial tolerance, check library quality |
| Low MS1 count | Missing MS1 spectra in mzML | Verify mzML contains MS1 scans |
| High mass SD | Interference or noise peaks | May indicate poor spectral quality |
| Large systematic mass error | Instrument miscalibration | Correction is applied automatically |
| NaN in RT predictions | Duplicate library RTs | Fixed: averaged fitted values at duplicates |
