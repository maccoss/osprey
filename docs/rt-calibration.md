# RT Calibration

Osprey uses LOESS (Locally Estimated Scatterplot Smoothing) regression to calibrate library retention times against measured retention times. This is essential because library RTs may come from different sources (predicted, measured on different systems, or iRT-normalized).

## Overview

RT calibration uses a **multi-file strategy** where the first file establishes the calibration and subsequent files reuse it:

```
File 1 (First file in experiment):
  └─→ Use ALL library peptides
  └─→ Assume linear mapping: library RT range ≈ mzML RT range
  └─→ Wide initial tolerance (20-30% of gradient range)
  └─→ Peak detection to find measured RTs
  └─→ LOESS curve fitting → calibration + residual SD

Files 2-N (Subsequent files):
  └─→ Start with calibration from File 1
  └─→ Use residual SD from File 1 as tolerance
  └─→ Refine calibration if needed (optional)
```

**Rationale**: Files within the same experiment have similar LC conditions, so the calibration from the first file provides a good starting point for subsequent files.

## First File Calibration

### Initial Assumptions

For the first file, we assume:

1. **Library RT range spans mzML RT range**: The gradient range in the library should roughly match the actual run
2. **Linear relationship**: Before LOESS fitting, assume linear mapping with wide tolerance
3. **Tolerance = 20-30% of gradient range**: Wide enough to catch all peptides regardless of RT shift

```
Library RT range: 10 - 60 min (50 min total)
Initial tolerance: 50 × 0.25 = 12.5 min

Any peptide within ±12.5 min of its library RT will be considered
```

### Processing All Peptides

For the first file, we use **all library peptides** (no stratified sampling):

```python
def calibrate_first_file(library, spectra):
    # Calculate initial wide tolerance (25% of RT range)
    library_rts = [e.retention_time for e in library]
    rt_range = max(library_rts) - min(library_rts)
    initial_tolerance = rt_range * 0.25

    calibration_points = []

    for entry in library:
        # Find spectra where:
        # 1. Precursor falls within isolation window
        # 2. Spectrum RT is within wide tolerance of library RT
        matching_spectra = [
            s for s in spectra
            if s.isolation_window.contains(entry.precursor_mz)
            and abs(s.retention_time - entry.retention_time) <= initial_tolerance
        ]

        # Run regression, detect peak
        coefficients = run_regression(entry, matching_spectra)
        peak = detect_peak(coefficients)

        if peak:
            calibration_points.append((entry.retention_time, peak.apex_rt))

    # Fit LOESS
    calibration = fit_loess(calibration_points)

    return calibration
```

### Why Not Stratified Sampling for First File?

| Approach | Use Case |
|----------|----------|
| All peptides (first file) | Need maximum calibration points, unknown RT shift |
| Stratified sampling | Faster, but assumes reasonable RT alignment already |

For the first file, we don't know how shifted the RTs are, so we use all peptides with a wide tolerance to ensure we capture enough calibration points.

## LOESS Fitting

### What is LOESS?

LOESS (Local Polynomial Regression) fits the data using weighted local polynomials:

1. For each query point x₀, find nearby training points
2. Weight points by distance using tricube kernel
3. Fit weighted polynomial (linear or quadratic)
4. Return predicted value

### Tricube Kernel

```
w(d) = (1 - |d|³)³  for |d| < 1
w(d) = 0            for |d| ≥ 1
```

### Bandwidth Parameter

The bandwidth (0-1) controls how much data to use for each local fit:

| Bandwidth | Effect |
|-----------|--------|
| 0.2 | Very local, captures fine structure, may overfit |
| 0.3 | **Default**: Good balance |
| 0.5 | Smoother, more robust |
| 0.8 | Nearly global, may underfit |

### Algorithm

```python
def loess_predict(library_rt, training_library_rts, training_measured_rts, bandwidth=0.3):
    n = len(training_library_rts)
    k = int(bandwidth * n)  # Number of neighbors

    # Find k nearest neighbors
    distances = [abs(rt - library_rt) for rt in training_library_rts]
    sorted_indices = argsort(distances)[:k]

    # Calculate weights (tricube)
    max_dist = distances[sorted_indices[-1]]
    weights = []
    for idx in sorted_indices:
        d = distances[idx] / max_dist
        w = (1 - d**3)**3
        weights.append(w)

    # Fit weighted linear regression
    X = [training_library_rts[i] for i in sorted_indices]
    Y = [training_measured_rts[i] for i in sorted_indices]

    # Weighted least squares: predict at library_rt
    predicted = weighted_linear_fit(X, Y, weights, library_rt)

    return predicted
```

## Calibration Statistics

After fitting, Osprey reports:

```
RT calibration (file 1): 2847 points, R²=0.9987, residual_SD=0.32 min
Calibrated RT tolerance: 0.96 min (3× residual SD)
```

- **n_points**: Number of peptides with detected peaks
- **R²**: Goodness of fit (1.0 = perfect)
- **residual_std**: Standard deviation of (measured - predicted)

## Subsequent Files

### Reusing First File Calibration

For files 2-N in the same experiment:

```python
def calibrate_subsequent_file(library, spectra, first_file_calibration):
    # Use calibration from first file
    calibration = first_file_calibration

    # Use residual SD from first file as tolerance
    tolerance = calibration.residual_std * 3.0  # 3σ

    # Optional: refine calibration with this file's data
    # (usually not necessary within same experiment)

    return calibration, tolerance
```

### Why Reuse Calibration?

Files within the same experiment typically have:
- Same LC column
- Same gradient program
- Same mobile phases
- Similar RT behavior

The calibration from the first file should be directly applicable. Only if there's significant drift between runs would recalibration be needed.

### Adaptive Tolerance

The calibrated RT tolerance uses the first file's residual SD:

```
tolerance = residual_std × rt_tolerance_factor
tolerance = max(tolerance, 0.5)  # Minimum 0.5 min
```

With typical values:
- residual_std = 0.3 min
- factor = 3.0
- tolerance = 0.9 min

This captures ~99.7% of true peptides (3σ for normal distribution).

## Multi-Experiment Handling

If processing files from **different experiments** (different LC conditions):

1. Group files by experiment
2. Calibrate first file of each experiment separately
3. Apply per-experiment calibration

```yaml
# Future configuration option
experiments:
  - name: "experiment_1"
    files: [sample1.mzML, sample2.mzML, sample3.mzML]
  - name: "experiment_2"
    files: [sample4.mzML, sample5.mzML]
```

## Extrapolation

For RTs outside the training range:

```
                LOESS fit
               ___________
              /           \
    Linear   /             \   Linear
    extrap. /               \ extrap.
         __/                 \__
        |                       |
     min_RT                  max_RT
```

Osprey uses linear extrapolation from the nearest training points.

## Fallback Behavior

If calibration fails (< min_calibration_points detected):

```
RT calibration failed: Insufficient calibration points: 23 < 50 required
Using fallback tolerance: 2.0 min
```

The fallback uses library RTs directly with a wider tolerance.

## Configuration Reference

```yaml
rt_calibration:
  enabled: true              # Enable/disable RT calibration
  loess_bandwidth: 0.3       # Fraction of data for local fits
  min_calibration_points: 50 # Minimum required for fitting
  rt_tolerance_factor: 3.0   # Multiplier for residual SD
  fallback_rt_tolerance: 2.0 # Used if calibration fails (minutes)
  initial_tolerance_fraction: 0.25  # 25% of RT range for first file
```

## Example Walkthrough

```
Experiment: 3 replicate injections
Library: 5000 peptides, RT range 10-60 min

=== File 1 (sample1.mzML) ===
Initial tolerance: (60-10) × 0.25 = 12.5 min
Processing all 5000 peptides...
Detected peaks for 3847 peptides
Fitting LOESS...
Result: R²=0.9991, residual_SD=0.28 min
Calibrated tolerance: 0.28 × 3 = 0.84 min

=== File 2 (sample2.mzML) ===
Using calibration from File 1
RT tolerance: 0.84 min
Processing with calibrated RTs...

=== File 3 (sample3.mzML) ===
Using calibration from File 1
RT tolerance: 0.84 min
Processing with calibrated RTs...
```

## Implementation

Key files:
- `crates/osprey-chromatography/src/calibration.rs` - LOESS fitting
- `crates/osprey/src/pipeline.rs` - Multi-file processing
- `crates/osprey-core/src/config.rs` - RTCalibrationConfig
