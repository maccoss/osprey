# Spectral Scoring

Osprey uses spectral similarity scoring to assess whether an observed spectrum matches its predicted library spectrum. This is the primary score used for target-decoy FDR computation during calibration.

## Overview

The scoring system follows the pyXcorrDIA methodology for calibration:

```
Workflow:
  1. For each library entry, find best matching spectrum (within RT tolerance)
  2. Score using LibCosine (primary) and XCorr (secondary)
  3. Extract MS1 isotope envelope → isotope cosine score
  4. Use LibCosine for target-decoy competition
  5. FDR filtering → calibration points
```

## Scoring Methods

### LibCosine (Primary Score)

LibCosine is a cosine similarity score with SMZ preprocessing (sqrt intensity × m/z²):

```
Preprocessing (per fragment):
  value = sqrt(intensity) × m/z²

Scoring:
  1. Match library fragments to observed peaks (within ppm tolerance)
  2. Apply SMZ preprocessing to both matched vectors
  3. L2 normalize both vectors
  4. Compute cosine similarity (dot product of normalized vectors)

Score range: 0-1 (1 = perfect match)
```

**Why SMZ preprocessing?**
- Square root dampens intensity differences (accounts for measurement variability)
- m/z² weighting emphasizes higher mass fragments (more informative)
- L2 normalization enables meaningful comparison regardless of total intensity

### XCorr (Secondary Score)

XCorr uses Comet-style preprocessing for cross-correlation scoring. See [XCorr Scoring](xcorr-scoring.md) for detailed implementation.

```
Experimental Spectrum Preprocessing:
  1. Bin into discrete bins (0-2000 m/z, bin_width = 1.0005079 Da)
  2. Apply sqrt transformation
  3. Windowing normalization (10 windows, normalize each to max=50, threshold 5%)
  4. Flanking bin subtraction (offset=75, removes local average)

Theoretical Spectrum (NO preprocessing - this is critical!):
  1. Bin fragment m/z values with unit intensity (1.0)
  2. NO windowing, NO sqrt, NO flanking subtraction
  3. The theoretical spectrum is just a "selector" for which bins to sum

Scoring:
  - XCorr = Σ experimental_preprocessed[fragment_bins] × 0.005
  - Equivalent to dot product since theoretical has unit intensities
```

**Key insight**: Comet does NOT preprocess the theoretical spectrum. It simply looks up values from the preprocessed experimental spectrum at each fragment bin position and sums them. Applying windowing to the theoretical spectrum would inflate scores by ~50x.

**Batch XCorr Preprocessing (pyXcorrDIA approach):**

For efficiency, Osprey preprocesses all spectra and library entries once per isolation window, then looks up XCorr at the best LibCosine match:

```rust
// Preprocess all spectra for XCorr ONCE per window
let spec_xcorr_preprocessed: Vec<Vec<f64>> = window_spectra
    .iter()
    .map(|spec| xcorr_scorer.preprocess_spectrum_for_xcorr(spec))
    .collect();

// Preprocess all library entries for XCorr ONCE per window
let entry_xcorr_preprocessed: HashMap<u32, Vec<f64>> = window_entries
    .iter()
    .map(|entry| (entry.id, xcorr_scorer.preprocess_library_for_xcorr(entry)))
    .collect();

// After finding best LibCosine match, lookup XCorr:
let xcorr_score = SpectralScorer::xcorr_from_preprocessed(
    &spec_xcorr_preprocessed[best_local_idx],
    &entry_xcorr_preprocessed[entry.id],
);
```

### Isotope Cosine Score (MS1 Quality)

The isotope cosine score compares the observed MS1 isotope envelope to the theoretical distribution calculated from the peptide's exact elemental composition:

```
1. Calculate elemental composition from peptide sequence:
   - Sum amino acid compositions (C, H, N, O, S)
   - Add terminal H2O

2. Calculate theoretical isotope distribution:
   - Use natural isotope abundances (C13=1.084%, N15=0.364%, etc.)
   - Polynomial expansion method for [M+0, M+1, M+2, M+3, M+4]

3. Extract observed envelope from MS1 spectrum:
   - Find peaks at M-1, M+0, M+1, M+2, M+3 (isotope spacing = 1.002868 / charge)
   - Match within ppm tolerance

4. Compute cosine similarity:
   - Align M+0..M+3 from both distributions
   - cosine = dot(observed, theoretical) / (||observed|| × ||theoretical||)
```

**Why exact elemental composition?**
- More accurate than averagine approximation
- Sulfur-containing peptides (Cys, Met) have distinct M+2 patterns
- Enables confident MS1-level validation

## Fragment Matching

Library fragments are matched to observed peaks using **closest m/z** within tolerance (pyXcorrDIA approach):

```rust
// Find closest matching experimental peak within tolerance
// Select peak with smallest m/z difference, NOT highest intensity
let mut best_mz_diff = f64::INFINITY;

for (&exp_mz, &exp_intensity) in spectrum.mzs.iter().zip(spectrum.intensities.iter()) {
    if tolerance.within_tolerance(lib_mz, exp_mz) {
        let mz_diff = (exp_mz - lib_mz).abs();
        if mz_diff < best_mz_diff {
            best_mz_diff = mz_diff;
            best_intensity = exp_intensity;
            best_mz = Some(exp_mz);
        }
    }
}
```

**Tolerance modes:**
- **ppm tolerance**: 10-20 ppm (HRAM instruments - Orbitrap, TOF)
- **Da tolerance**: 0.5 Da (unit resolution instruments)

**Metrics computed:**
- `n_matched`: Number of library fragments with matches
- `n_library`: Total library fragments
- `mass_errors`: Vector of (observed - theoretical) errors for calibration

## Target-Decoy Competition

During calibration, each target competes against its paired decoy:

```
Target-decoy pairing: decoy_id = target_id | 0x80000000

Competition rules (pyXcorrDIA):
  1. Compare LibCosine scores for target and decoy
  2. Higher score wins
  3. Ties go to decoy (conservative)
  4. Winner enters FDR calculation

FDR calculation:
  1. Sort winners by score (descending)
  2. FDR = cumulative_decoys / cumulative_targets
  3. Filter to targets at ≤1% FDR
```

## Configuration

```yaml
# Fragment matching tolerance (based on resolution mode)
resolution_mode: unit  # Uses 0.5 Da
# OR
resolution_mode:
  HRAM:
    tolerance_ppm: 10.0  # Uses 10 ppm

# Precursor tolerance for MS1 isotope extraction
precursor_tolerance_ppm: 10.0
```

## Feature Set

Spectral scoring populates these fields in CalibrationMatch:

| Field | Source | Description |
|-------|--------|-------------|
| `score` | LibCosine | Cosine similarity score (0-1) |
| `xcorr_score` | XCorr | Comet-style cross-correlation |
| `isotope_cosine_score` | MS1 | Isotope envelope match (0-1) |
| `n_matched_fragments` | Matching | Number of library fragments matched |
| `ms2_mass_errors` | Matching | Fragment mass errors (ppm) |
| `ms1_ppm_error` | MS1 | Precursor mass error (ppm) |

## Implementation

Key files:
- `crates/osprey-scoring/src/lib.rs` - SpectralScorer (XCorr preprocessing)
- `crates/osprey-scoring/src/batch.rs` - LibCosineScorer, batch scoring
- `crates/osprey-core/src/isotope.rs` - Isotope distribution calculations
- `crates/osprey/src/pipeline.rs` - Integration in calibration workflow

### LibCosineScorer API

```rust
use osprey_scoring::batch::LibCosineScorer;

// Create scorer with ppm tolerance
let scorer = LibCosineScorer::hram(10.0);  // 10 ppm

// Score with mass error collection
let (score, mass_errors) = scorer.score_with_errors(&library_entry, &spectrum);
println!("LibCosine: {:.3}", score);
println!("Matched fragments: {}", mass_errors.len());

// Fragment matching details
let match_result = scorer.match_fragments(&library_entry, &spectrum);
println!("Matched {}/{} fragments", match_result.n_matched, library_entry.fragments.len());
```

### SpectralScorer API (XCorr)

```rust
use osprey_scoring::SpectralScorer;

let xcorr_scorer = SpectralScorer::new();

// Preprocess spectrum for XCorr
let spec_preprocessed = xcorr_scorer.preprocess_spectrum_for_xcorr(&spectrum);

// Preprocess library entry for XCorr
let lib_preprocessed = xcorr_scorer.preprocess_library_for_xcorr(&entry);

// Compute XCorr from preprocessed vectors
let xcorr = SpectralScorer::xcorr_from_preprocessed(&spec_preprocessed, &lib_preprocessed);
```

### Isotope Scoring API

```rust
use osprey_core::{IsotopeEnvelope, peptide_isotope_cosine};

// Extract isotope envelope from MS1 spectrum
let envelope = IsotopeEnvelope::extract(
    ms1_spectrum,
    precursor_mz,
    charge,
    tolerance_ppm,
);

// Calculate isotope cosine score
if envelope.has_m0() {
    let isotope_score = peptide_isotope_cosine(&sequence, &envelope.intensities);
    println!("Isotope cosine: {:?}", isotope_score);
}
```

## Example

```
Peptide: PEPTIDEK (charge 2+)
Library fragments: y3(348.2), y4(461.3), y5(576.3), b3(357.2), b4(458.2)

Observed spectrum at apex RT=25.3 min:
  m/z: [348.201, 461.299, 576.305, 400.1, 550.2, ...]
  intensity: [1000, 800, 600, 200, 150, ...]

Fragment matching (10 ppm tolerance):
  y3: 348.2 → 348.201 ✓ (error: +2.9 ppm, closest m/z)
  y4: 461.3 → 461.299 ✓ (error: -2.2 ppm)
  y5: 576.3 → 576.305 ✓ (error: +8.7 ppm)
  b3: 357.2 → no match ✗
  b4: 458.2 → no match ✗

Metrics:
  n_matched: 3
  n_library: 5
  avg_ms2_error: +3.1 ppm

LibCosine scoring:
  Library SMZ: [sqrt(100)×348.2², sqrt(80)×461.3², sqrt(60)×576.3²]
  Observed SMZ: [sqrt(1000)×348.2², sqrt(800)×461.3², sqrt(600)×576.3²]
  (normalized and dot product computed)
  LibCosine: 0.95

XCorr scoring:
  (preprocessed via windowing + flanking subtraction)
  XCorr: 0.42

Isotope cosine (from MS1):
  Theoretical: [0.55, 0.30, 0.12, 0.03]  (from C41H68N10O15)
  Observed:    [0.52, 0.31, 0.13, 0.04]
  Isotope cosine: 0.998
```

## Notes

1. **Score interpretation**: LibCosine > 0.7 is typically a reasonable match
2. **Closest m/z matching**: Ensures accurate mass error calculation for calibration
3. **Batch preprocessing**: XCorr is precomputed for ALL spectrum-library pairs, not just the best match
4. **Isotope scoring**: Only available when MS1 spectra are present in the mzML file
5. **Multiple scores**: LibCosine for FDR, but all scores written to debug CSV for analysis
