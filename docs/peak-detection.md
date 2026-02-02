# Peak Detection

After ridge regression, each peptide has a **coefficient time series** - a sequence of coefficients across all MS2 spectra. Peak detection identifies where the peptide elutes and extracts chromatographic features.

## Overview

```
Regression coefficients across RT:

Coefficient
    │
1.0 │           ╭───╮
    │          ╱     ╲
0.5 │         ╱       ╲
    │        ╱         ╲
0.0 │───────╱───────────╲───────
    └─────────────────────────────► RT
           start apex  end
```

## Coefficient Time Series

For each library entry, we collect (RT, coefficient) pairs:

```python
# After processing all spectra
entry_series = {}

for result in regression_results:
    for lib_id, coefficient in zip(result.library_ids, result.coefficients):
        if lib_id not in entry_series:
            entry_series[lib_id] = []
        entry_series[lib_id].append((result.retention_time, coefficient))

# Sort by RT
for lib_id in entry_series:
    entry_series[lib_id].sort(key=lambda x: x[0])
```

## Peak Finding Algorithm

### Simple Peak Detection

```python
def detect_peaks(time_series, min_height=0.05, min_width=3):
    peaks = []
    n = len(time_series)

    for i in range(1, n - 1):
        rt, coef = time_series[i]

        # Check if local maximum
        if coef <= time_series[i-1][1] or coef <= time_series[i+1][1]:
            continue

        # Check minimum height
        if coef < min_height:
            continue

        # Find peak boundaries
        start = find_peak_start(time_series, i)
        end = find_peak_end(time_series, i)

        # Check minimum width
        if end - start + 1 < min_width:
            continue

        peak = Peak(
            apex_idx=i,
            apex_rt=rt,
            apex_coefficient=coef,
            start_idx=start,
            end_idx=end,
            start_rt=time_series[start][0],
            end_rt=time_series[end][0]
        )
        peaks.append(peak)

    return peaks
```

### Finding Peak Boundaries

```python
def find_peak_start(series, apex_idx):
    """Walk left until coefficient drops to baseline or rises again."""
    baseline = 0.05 * series[apex_idx][1]  # 5% of apex

    for i in range(apex_idx - 1, -1, -1):
        if series[i][1] < baseline:
            return i + 1
        if series[i][1] > series[i + 1][1]:  # Rising again = different peak
            return i + 1

    return 0

def find_peak_end(series, apex_idx):
    """Walk right until coefficient drops to baseline or rises again."""
    baseline = 0.05 * series[apex_idx][1]

    for i in range(apex_idx + 1, len(series)):
        if series[i][1] < baseline:
            return i - 1
        if series[i][1] > series[i - 1][1]:
            return i - 1

    return len(series) - 1
```

## Selecting the Best Peak

When multiple peaks are detected, select the one closest to the expected RT:

```python
def find_best_peak(peaks, expected_rt, rt_tolerance):
    valid_peaks = [
        p for p in peaks
        if abs(p.apex_rt - expected_rt) <= rt_tolerance
    ]

    if not valid_peaks:
        return None

    # Return peak closest to expected RT
    return min(valid_peaks, key=lambda p: abs(p.apex_rt - expected_rt))
```

## Peak Quality Metrics

### Peak Boundaries

```python
@dataclass
class PeakBoundaries:
    start_rt: float      # Peak start (minutes)
    end_rt: float        # Peak end (minutes)
    apex_rt: float       # Peak apex (minutes)
    apex_coefficient: float  # Coefficient at apex
    integrated_area: float   # Sum of coefficients (not for quant)
```

### Peak Quality Flags

```python
@dataclass
class PeakQuality:
    emg_fit_r2: float      # R² of EMG fit (0-1)
    is_split: bool         # Peak appears split
    is_truncated: bool     # Truncated at gradient edge
    has_shoulder: bool     # Shoulder detected
    width_percentile: float  # Width relative to other peaks
```

## EMG Peak Fitting (TODO)

Exponentially Modified Gaussian (EMG) models chromatographic peak shape:

```
EMG(t) = (A * σ * √(2π)) / (2τ) * exp((σ²)/(2τ²) - (t-μ)/τ) * erfc((σ/τ - (t-μ)/σ) / √2)

Where:
  A = amplitude
  μ = Gaussian mean (RT)
  σ = Gaussian standard deviation
  τ = exponential decay constant (tailing)
```

EMG fitting provides:
- More accurate apex RT estimation
- Peak asymmetry quantification
- Robust boundary determination

## Integrated Area

The integrated area (sum of coefficients) is used for:
- Peak ranking within a spectrum
- Quality control

**Important**: Do NOT use integrated area for quantification! The coefficients are relative within each spectrum and not comparable across spectra.

## Peak Width

```python
def calculate_fwhm(series, apex_idx):
    """Full Width at Half Maximum."""
    half_max = series[apex_idx][1] / 2

    # Find left half-max point
    left = apex_idx
    for i in range(apex_idx - 1, -1, -1):
        if series[i][1] < half_max:
            # Interpolate
            left = i + (half_max - series[i][1]) / (series[i+1][1] - series[i][1])
            break

    # Find right half-max point
    right = apex_idx
    for i in range(apex_idx + 1, len(series)):
        if series[i][1] < half_max:
            right = i - (half_max - series[i][1]) / (series[i-1][1] - series[i][1])
            break

    fwhm = series[int(right)][0] - series[int(left)][0]
    return fwhm
```

## Chromatographic Features

Features extracted from peaks for scoring:

| Feature | Description |
|---------|-------------|
| `peak_apex` | Maximum coefficient value |
| `peak_area` | Integrated coefficient area |
| `peak_width` | FWHM in minutes |
| `peak_symmetry` | Leading/trailing ratio |
| `n_contributing_scans` | Number of non-zero coefficients |
| `coefficient_stability` | Variance near apex |
| `peak_sharpness` | Rate of rise/fall |
| `peak_prominence` | Apex / baseline ratio |

## Configuration

```yaml
# Peak detection is controlled by code defaults
# Future: expose these in config

min_peak_height: 0.05    # Minimum coefficient for peak
min_peak_width: 3        # Minimum scans for valid peak
```

## Implementation

Key files:
- `crates/osprey-chromatography/src/lib.rs` - PeakDetector
- `crates/osprey-core/src/types.rs` - PeakBoundaries, PeakQuality
- `crates/osprey/src/pipeline.rs` - Peak detection in calibration discovery

## Example

```
Peptide: PEPTIDEK
Library RT: 25.3 min

Coefficient time series:
  RT 22.0: 0.00
  RT 22.5: 0.00
  RT 23.0: 0.02
  RT 23.5: 0.08  ← start_idx
  RT 24.0: 0.25
  RT 24.5: 0.58
  RT 25.0: 0.89  ← apex_idx
  RT 25.5: 0.72
  RT 26.0: 0.35
  RT 26.5: 0.12
  RT 27.0: 0.03  ← end_idx
  RT 27.5: 0.00

Peak detected:
  start_rt: 23.5
  apex_rt: 25.0
  end_rt: 27.0
  apex_coefficient: 0.89
  integrated_area: 3.04
  fwhm: 1.8 min
```

## Multiple Peaks

A peptide may have multiple peaks due to:
- Isomers eluting at different RTs
- False positive signals
- Peak splitting

Osprey selects the peak closest to the expected (calibrated) RT.
