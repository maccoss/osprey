# Mass Calibration

Osprey performs mass calibration for both MS1 (precursor) and MS2 (fragment) measurements during the calibration discovery phase. This enables tighter tolerance windows for the main search and more accurate quantification.

## Overview

Mass calibration follows the pyXcorrDIA methodology:

```
┌─────────────────────────────────────────────────────────────────┐
│              Mass Calibration Workflow                           │
│                                                                  │
│  1. Score library vs spectra using XCorr (unit resolution bins) │
│  2. Target-decoy competition (E-value) → filter to 1% FDR      │
│  3. Collect mass errors ONLY from confident matches:            │
│     - MS1: Extract M+0 peak from MS1 spectrum                   │
│     - MS2: Top-3 fragment mass errors (binary search)           │
│  4. Calculate statistics: mean, median, SD                      │
│  5. Set adjusted tolerance: |mean| + 3×SD                       │
└─────────────────────────────────────────────────────────────────┘
```

**Key insight**: Mass calibration uses ONLY high-confidence matches (targets passing 1% FDR), not all matches. This ensures the calibration statistics are not contaminated by false positives.

## MS1 (Precursor) Calibration

### Isotope Envelope Extraction

MS1 calibration is based on extracting the monoisotopic (M+0) peak from MS1 spectra:

```
For each confident match:
  1. Find nearest MS1 spectrum to the best MS2 retention time
  2. Calculate expected isotope m/z values:
     - Isotope spacing = 1.002868 Da / charge
     - Expected m/z: M-1, M+0, M+1, M+2, M+3
  3. Search for peaks within ppm tolerance
  4. Record M+0 observed m/z and intensity
  5. Calculate PPM error: ((observed - theoretical) / theoretical) × 10⁶
```

### MS1 Error Statistics

From FDR-filtered confident matches:

```
MS1 calibration: mean=-2.50 ppm, median=-2.40 ppm, SD=0.80 ppm (n=150)
Adjusted MS1 tolerance: 4.90 ppm (|mean| + 3×SD)
```

**Output fields:**
- `mean`: Systematic mass error (can be used for correction)
- `median`: Robust estimate of systematic error
- `sd`: Random error component
- `count`: Number of observations
- `adjusted_tolerance`: Recommended search tolerance

### MS1 Histogram

Osprey optionally generates an ASCII histogram of MS1 errors:

```
--- MS1 (precursor) Mass Error Histogram (ppm) ---
   -4.5 |    12 | ████
   -3.5 |    28 | █████████
   -2.5 |    67 | ██████████████████████
   -1.5 |    31 | ██████████
   -0.5 |     8 | ██
    0.5 |     4 | █
        | mean=-2.50, median=-2.40, SD=0.80 (n=150)
```

## MS2 (Fragment) Calibration

### Fragment Mass Error Collection

MS2 calibration uses mass errors from top-3 fragment matches at the best XCorr spectrum:

```
For each confident match:
  Top-3 fragments (by library intensity) are matched via binary search:
    1. Find closest observed peak within tolerance
    2. Calculate error: ((observed_mz - library_mz) / library_mz) × 10⁶ (ppm)
    3. Add to MS2 error collection

Typical yield: 1-3 mass errors per peptide (up to 3 from top-3 fragments)
With ~10K calibration peptides → ~20-30K MS2 data points
```

**Note**: Only the top-3 most intense library fragments are matched (using binary search, O(3 × log n)). This is much faster than full fragment matching while providing sufficient data for calibration.

### MS2 Error Statistics

```
MS2 calibration: mean=1.20 ppm, median=1.10 ppm, SD=1.00 ppm (n=500)
Adjusted MS2 tolerance: 4.20 ppm (|mean| + 3×SD)
```

## CalibrationParams Structure

Mass calibration results are stored in the `CalibrationParams` structure:

```rust
pub struct CalibrationParams {
    pub metadata: CalibrationMetadata,
    pub ms1_calibration: MzCalibration,
    pub ms2_calibration: MzCalibration,
    pub rt_calibration: RTCalibrationParams,
}

pub struct MzCalibration {
    pub mean: f64,              // Mean PPM error
    pub median: f64,            // Median PPM error
    pub sd: f64,                // Standard deviation
    pub count: usize,           // Number of observations
    pub unit: String,           // Always "ppm"
    pub adjusted_tolerance: Option<f64>,
    pub histogram: Option<MzHistogram>,
    pub calibrated: bool,
}
```

## MzQCData Collection

During calibration, errors are collected in the `MzQCData` structure:

```rust
pub struct MzQCData {
    ms1_errors: Vec<f64>,
    ms2_errors: Vec<f64>,
}

impl MzQCData {
    pub fn add_ms1_error(&mut self, error: f64);
    pub fn add_ms2_error(&mut self, error: f64);
    pub fn n_ms1(&self) -> usize;
    pub fn n_ms2(&self) -> usize;
}
```

## Calibration Calculation

The `calculate_mz_calibration` function computes statistics from collected errors:

```rust
pub fn calculate_mz_calibration(qc_data: &MzQCData) -> (MzCalibration, MzCalibration) {
    let ms1_calibration = if qc_data.n_ms1() > 0 {
        let mean = mean(&qc_data.ms1_errors);
        let median = median(&qc_data.ms1_errors);
        let sd = std_dev(&qc_data.ms1_errors);

        MzCalibration {
            mean,
            median,
            sd,
            count: qc_data.n_ms1(),
            unit: "ppm".to_string(),
            adjusted_tolerance: Some(mean.abs() + 3.0 * sd),
            histogram: Some(build_histogram(&qc_data.ms1_errors)),
            calibrated: true,
        }
    } else {
        MzCalibration::uncalibrated()
    };

    // Similar for MS2...
    (ms1_calibration, ms2_calibration)
}
```

## Using Calibration in Search

After calibration, the adjusted tolerances can be used for the main search:

```rust
// Get effective tolerance (uses calibrated if available)
let ms1_tolerance = calibration.ms1_calibration.effective_tolerance(config.precursor_tolerance_ppm);
let ms2_tolerance = calibration.ms2_calibration.effective_tolerance(config.fragment_tolerance_ppm);

// Optionally apply mass correction
let corrected_mz = observed_mz - (calibration.ms1_calibration.mean * observed_mz / 1e6);
```

## JSON Output

Calibration results are saved to a JSON file for reproducibility:

```json
{
  "metadata": {
    "num_confident_peptides": 150,
    "num_sampled_precursors": 2000,
    "calibration_successful": true,
    "timestamp": "2024-01-15T10:30:00Z"
  },
  "ms1_calibration": {
    "mean": -2.5,
    "median": -2.4,
    "sd": 0.8,
    "count": 150,
    "unit": "ppm",
    "adjusted_tolerance": 4.9,
    "calibrated": true
  },
  "ms2_calibration": {
    "mean": 1.2,
    "median": 1.1,
    "sd": 1.0,
    "count": 500,
    "unit": "ppm",
    "adjusted_tolerance": 4.2,
    "calibrated": true
  },
  "rt_calibration": {
    "method": "LOESS",
    "residual_sd": 0.8,
    "n_points": 150,
    "r_squared": 0.98
  }
}
```

## Configuration

```yaml
# Precursor tolerance for MS1 isotope extraction
precursor_tolerance_ppm: 10.0

# Fragment tolerance (from resolution mode)
resolution_mode:
  HRAM:
    tolerance_ppm: 10.0

# Calibration is computed from matches at 1% FDR
# (No separate configuration needed)
```

## Implementation

Key files:
- `crates/osprey-chromatography/src/calibration/mod.rs` - CalibrationParams, MzCalibration
- `crates/osprey-chromatography/src/calibration/io.rs` - JSON serialization
- `crates/osprey/src/pipeline.rs` - Calibration workflow integration
- `crates/osprey-core/src/types.rs` - IsotopeEnvelope extraction

### Calibration Workflow in Pipeline

```rust
// Collect errors ONLY from FDR-passing targets
let mut mz_qc_data = MzQCData::new();

for m in &competition_result.passing_targets {
    // MS1 error from isotope envelope
    if let Some(ms1_error) = m.ms1_ppm_error {
        mz_qc_data.add_ms1_error(ms1_error);
    }

    // MS2 errors from matched fragments
    for &ms2_error in &m.ms2_mass_errors {
        mz_qc_data.add_ms2_error(ms2_error);
    }
}

// Calculate calibration statistics
let (ms1_calibration, ms2_calibration) = calculate_mz_calibration(&mz_qc_data);

log::info!(
    "MS1 calibration: mean={:.2} ppm, SD={:.2} ppm (from {} observations)",
    ms1_calibration.mean,
    ms1_calibration.sd,
    mz_qc_data.n_ms1()
);
```

## Best Practices

1. **FDR filtering is critical**: Always use only high-confidence matches for calibration
2. **Check observation count**: Low counts may indicate poor MS1 extraction
3. **Monitor systematic error**: Large mean values may indicate instrument miscalibration
4. **Use adjusted tolerance**: `|mean| + 3×SD` captures ~99.7% of true matches
5. **Save calibration JSON**: Enables reproducibility and troubleshooting

## Troubleshooting

| Symptom | Possible Cause | Solution |
|---------|----------------|----------|
| Low MS1 count | Missing MS1 spectra in mzML | Check mzML contains MS1 |
| High MS1 SD | Interference or noise | Use stricter precursor tolerance |
| Large systematic error | Instrument miscalibration | Consider applying correction |
| No calibration | Too few confident matches | Check library or FDR threshold |
