# RT Calibration

Osprey uses LOESS (Locally Estimated Scatterplot Smoothing) regression to calibrate library retention times against measured retention times. The calibration phase uses XCorr scoring with E-value-based target-decoy competition to identify high-confidence matches for fitting.

## Overview

RT calibration follows the pyXcorrDIA methodology:

```
┌─────────────────────────────────────────────────────────────────┐
│              RT Calibration Workflow                             │
│                                                                  │
│  1. Generate decoys (enzyme-aware reversal)                     │
│  2. Sample ~2000 peptides for calibration                       │
│  3. Score spectra via XCorr (unit resolution bins, E-value)     │
│  4. Target-decoy competition → winners                          │
│  5. Filter to targets at ≤1% FDR                                │
│  6. Fit LOESS on (library_RT, measured_RT) pairs                │
│  7. Calculate residual SD → RT tolerance for main search        │
└─────────────────────────────────────────────────────────────────┘
```

**Key insight**: Only high-confidence matches (targets passing 1% FDR) are used for LOESS fitting. This ensures the calibration is not corrupted by false matches.

## Calibration Discovery Phase

### Peptide Sampling

For large libraries (millions of entries), scoring all entries would be too slow. Osprey samples a representative subset:

```rust
// Default: 2000 peptides (targets + their paired decoys)
// Samples across the precursor m/z range for coverage
let calibration_library = sample_calibration_peptides(&library, 2000);
```

### XCorr Scoring

Each sampled peptide is scored against all spectra within RT tolerance using XCorr:

```
For each library entry:
  1. Filter spectra by RT tolerance + top-3 fragment match
  2. Calculate XCorr against all passing spectra (unit resolution bins: 2001 bins)
  3. Select best spectrum by XCorr
  4. Calculate E-value from XCorr survival function
  5. Collect top-3 fragment mass errors at best spectrum
```

Calibration always uses unit resolution bins (1.0005 Da, 2001 bins) regardless of data type for ~50x faster scoring on HRAM data.

### Target-Decoy Competition

Each target competes against its paired decoy using E-values:

```
Target-decoy pairing: decoy_id = target_id | 0x80000000

For each target-decoy pair:
  if target_evalue < decoy_evalue:  # Lower E-value is better
    winner = target
  else:  # Including ties - conservative
    winner = decoy

Sort winners by E-value (ascending)
FDR = cumulative_decoys / cumulative_targets
Filter to targets where FDR ≤ 0.01
```

### LOESS Fitting

From the FDR-filtered confident matches, extract (library_RT, measured_RT) pairs:

```python
calibration_points = []
for match in confident_targets:
    calibration_points.append((match.library_rt, match.measured_rt))

rt_calibration = fit_loess(calibration_points)
```

## LOESS Algorithm

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

### Implementation

```rust
pub struct RTCalibrator {
    bandwidth: f64,      // Default: 0.3
    degree: usize,       // Default: 1 (linear)
    min_points: usize,   // Default: 50
    robustness_iter: usize,  // Default: 2
}

impl RTCalibrator {
    pub fn fit(&self, library_rts: &[f64], measured_rts: &[f64]) -> Result<RTCalibration>;
    pub fn predict(&self, library_rt: f64) -> f64;
}
```

## Calibration Statistics

After fitting, Osprey reports:

```
RT calibration: 150 peptides at 1% FDR (from 1847 target wins, 153 decoy wins)
RT calibration fit: R²=0.9987, residual_SD=0.32 min
Calibrated RT tolerance: 0.96 min (3× residual SD)
```

- **n_points**: Number of peptides passing 1% FDR
- **R²**: Goodness of fit (1.0 = perfect)
- **residual_sd**: Standard deviation of (measured - predicted)
- **rt_tolerance**: 3 × residual_sd (captures ~99.7% of true peptides)

## Multi-File Processing

### First File

The first file establishes calibration:

```
File 1 (sample1.mzML):
  1. Run full calibration discovery
  2. Score ~2000 sampled peptides
  3. Target-decoy competition → 1% FDR
  4. Fit LOESS → calibration parameters
  5. Save calibration for subsequent files
```

### Subsequent Files

Subsequent files reuse the calibration:

```
Files 2-N:
  1. Load calibration from File 1
  2. Use calibrated RT tolerance: 3 × residual_sd
  3. Run main search with tight RT window
  (Optional: refine calibration with this file's data)
```

**Rationale**: Files within the same experiment have similar LC conditions, so the calibration from the first file is directly applicable.

## Configuration

```yaml
rt_calibration:
  enabled: true              # Enable/disable RT calibration
  loess_bandwidth: 0.3       # Fraction of data for local fits
  min_calibration_points: 50 # Minimum required for fitting
  rt_tolerance_factor: 3.0   # Multiplier for residual SD
  fallback_rt_tolerance: 2.0 # Used if calibration fails (minutes)
```

## Extrapolation

For RTs outside the training range, LOESS uses linear extrapolation:

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

## Fallback Behavior

If calibration fails (insufficient points at 1% FDR):

```
RT calibration failed: Insufficient calibration points: 23 < 50 required
Using fallback tolerance: 2.0 min
```

The fallback uses library RTs directly with a wider tolerance.

## Implementation Files

Key files:
- `crates/osprey-chromatography/src/calibration/rt.rs` - LOESS fitting
- `crates/osprey-chromatography/src/calibration/mod.rs` - CalibrationParams
- `crates/osprey-scoring/src/batch.rs` - XCorr calibration scoring, E-value calculation
- `crates/osprey-fdr/src/controller.rs` - Target-decoy competition
- `crates/osprey/src/pipeline.rs` - Calibration workflow

### Pipeline Integration

```rust
// Run calibration scoring (XCorr + top-3 fragment matching, unit resolution bins)
let matches = run_xcorr_calibration_scoring(
    &calibration_library,
    &spectra,
    Some(&ms1_index),
    fragment_tolerance,
    precursor_tolerance_ppm,
    initial_rt_tolerance,
);

// Target-decoy competition
let fdr_controller = FdrController::new(0.01);  // 1% FDR
let competition_result = fdr_controller.compete_and_filter(matches);

// Extract calibration points from passing targets
let mut library_rts = Vec::new();
let mut measured_rts = Vec::new();
for m in &competition_result.passing_targets {
    library_rts.push(m.library_rt);
    measured_rts.push(m.measured_rt);
}

// Fit LOESS
let calibrator = RTCalibrator::with_config(config);
let rt_calibration = calibrator.fit(&library_rts, &measured_rts)?;
```

## Example Walkthrough

```
Experiment: DIA analysis with predicted library

=== Calibration Discovery ===
Library: 100,000 peptides
Sampled for calibration: 2,000 targets + 2,000 decoys

XCorr scoring (unit resolution bins):
  Initial RT tolerance: 5.0 min (wide)
  Fragment tolerance: 10 ppm
  Matched peptides: 3,847

Target-decoy competition:
  Target wins: 1,847
  Decoy wins: 153
  Total winners: 2,000

FDR filtering (1%):
  Confident targets: 150
  (at score 0.72, FDR = 0.98%)

=== LOESS Fitting ===
Input: 150 (library_RT, measured_RT) pairs
Bandwidth: 0.3
Result:
  R² = 0.9987
  Residual SD = 0.32 min
  RT tolerance = 0.96 min (3× SD)

=== Calibrated Search ===
For main search, use:
  - RT tolerance: 0.96 min (instead of 5.0 min)
  - RT prediction: LOESS(library_RT) → expected_RT
```

## Troubleshooting

| Symptom | Possible Cause | Solution |
|---------|----------------|----------|
| Low confident count | Poor library match | Check library source |
| High residual SD | Non-linear RT shift | Increase LOESS bandwidth |
| Low R² | Outliers or multiple modes | Check for gradient issues |
| Calibration fails | Too few confident matches | Widen initial RT tolerance |

## Best Practices

1. **Sample diverse m/z range**: Calibration peptides should span the full precursor range
2. **Use 1% FDR**: Stricter FDR gives cleaner calibration points
3. **Check residual SD**: Should be < 1 min for good chromatography
4. **Monitor R²**: Should be > 0.99 for predictable RT behavior
5. **Save calibration JSON**: Enables reproducibility and troubleshooting
