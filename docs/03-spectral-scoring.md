# Spectral Scoring

Osprey uses spectral similarity scoring to assess whether an observed spectrum matches its predicted library spectrum. This is the primary score used for target-decoy FDR computation during calibration.

## Overview

The scoring system uses XCorr with E-value for calibration:

```
Calibration Workflow:
  1. For each library entry, calculate XCorr against ALL spectra in RT window
  2. Select best spectrum by XCorr (RT selection)
  3. Collect MS2 mass errors from top-3 fragment matches at best spectrum
  4. Calculate E-value from XCorr survival function (Comet-style)
  5. Extract MS1 isotope envelope → isotope cosine score
  6. Use E-value for target-decoy competition
  7. FDR filtering → calibration points
```

## Scoring Methods

### XCorr (Primary Calibration Score)

XCorr uses Comet-style preprocessing for cross-correlation scoring. It is calculated for ALL spectra in the RT tolerance window and used to select the best retention time. See [XCorr Scoring](04-xcorr-scoring.md) for detailed implementation.

Calibration XCorr always uses **unit resolution bins** (2001 bins, 1.0005 m/z) regardless of whether the data is unit resolution or HRAM. This provides ~50x faster scoring for HRAM data while maintaining sufficient discriminating power for calibration.

```
Experimental Spectrum Preprocessing:
  1. Bin into discrete bins (Comet BIN macro: bin_width = 1.0005079 Da, offset = 0.4)
  2. Apply sqrt transformation
  3. Windowing normalization (10 windows, normalize each to max=50, threshold 5%)
  4. Flanking bin subtraction (prefix-sum O(n), offset=75, removes local average)

Theoretical Spectrum (NO preprocessing - this is critical!):
  1. Bin fragment m/z values with unit intensity (1.0)
  2. NO windowing, NO sqrt, NO flanking subtraction
  3. The theoretical spectrum is just a "selector" for which bins to sum

Scoring:
  - XCorr = Σ experimental_preprocessed[fragment_bins] × 0.005
  - Computed via BLAS sdot (ndarray ArrayView1::dot)
```

**Key insight**: Comet does NOT preprocess the theoretical spectrum. It simply looks up values from the preprocessed experimental spectrum at each fragment bin position and sums them.

### E-Value (Significance Score)

Osprey computes Comet-style E-values from the XCorr survival function for each peptide:

```
E-value calculation:
  1. Collect ALL XCorr scores in RT window (~700 spectra typically)
  2. Sort scores descending, build survival function
  3. Fit linear regression: log10(survival_count) = intercept - slope × xcorr
  4. E-value = 10^(intercept - slope × best_xcorr)
```

Lower E-value = better statistical significance. E-value is used for target-decoy competition.

### Top-3 Fragment Matching (MS2 Calibration)

After XCorr selects the best spectrum, Osprey collects MS2 mass errors from the top-3 most intense library fragments using binary search:

```
For the best XCorr spectrum:
  1. Find top 3 library fragments by intensity
  2. Binary search each in sorted observed m/z array
  3. Find closest peak within tolerance window
  4. Compute signed mass error (ppm or Th)

Returns: (has_match: bool, mass_errors: Vec<f64>)
Cost: O(3 × log n_peaks) per peptide
```

These mass errors feed into MS2 mass calibration for the main search phase.

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

## Pre-Filtering

Before XCorr scoring, candidates are filtered using the `has_top3_fragment_match` function — a fast binary search check that requires at least 1 of the top-3 most intense library fragments to be present in the observed spectrum. This eliminates candidates with no spectral evidence before expensive XCorr preprocessing.

**Tolerance modes:**
- **ppm tolerance**: 10-20 ppm (HRAM instruments - Orbitrap, TOF)
- **Da tolerance**: 0.5 Da (unit resolution instruments)

## Target-Decoy Competition

During calibration, each target competes against its paired decoy:

```
Target-decoy pairing: decoy_id = target_id | 0x80000000

Competition rules:
  1. Compare E-values for target and decoy
  2. Lower E-value wins
  3. Ties go to decoy (conservative)
  4. Winner enters FDR calculation

FDR calculation:
  1. Sort winners by E-value (ascending - lower is better)
  2. FDR = cumulative_decoys / cumulative_targets
  3. Filter to targets at ≤1% FDR
```

## Scoring Flow Summary

```
Per peptide (target or decoy):
┌─────────────────────────────────────────────────────────────┐
│ 1. Filter by RT tolerance + top-3 fragment match            │
│                                                             │
│ 2. Calculate XCorr for ALL passing spectra                  │
│    → Collect all scores for E-value calculation             │
│    → Track best XCorr and its spectrum index                │
│                                                             │
│ 3. At best XCorr spectrum:                                  │
│    → Collect MS2 mass errors from top-3 fragments           │
│    → Extract MS1 isotope envelope                           │
│                                                             │
│ 4. Calculate E-value from XCorr survival function           │
│                                                             │
│ 5. Output: E-value, XCorr, MS2 errors, isotope_cosine      │
└─────────────────────────────────────────────────────────────┘
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
| `evalue` | XCorr survival | Comet-style E-value (lower = better) |
| `score` | XCorr | Same as xcorr_score (primary score) |
| `xcorr_score` | XCorr | Comet-style cross-correlation |
| `isotope_cosine_score` | MS1 | Isotope envelope match (0-1) |
| `n_matched_fragments` | Top-3 matching | Number of top-3 fragments matched (0-3) |
| `ms2_mass_errors` | Top-3 matching | Fragment mass errors (ppm or Th) |
| `ms1_error` | MS1 | Precursor mass error (ppm or Th) |

## Implementation

Key files:
- `crates/osprey-scoring/src/lib.rs` - SpectralScorer (XCorr preprocessing), `top3_fragment_match_with_errors`
- `crates/osprey-scoring/src/batch.rs` - E-value calculation, calibration batch scoring
- `crates/osprey-core/src/isotope.rs` - Isotope distribution calculations
- `crates/osprey/src/pipeline.rs` - Integration in calibration workflow

### SpectralScorer API (XCorr)

```rust
use osprey_scoring::SpectralScorer;

// Create scorer (always unit resolution for calibration)
let xcorr_scorer = SpectralScorer::new();

// Preprocess spectrum for XCorr
let spec_preprocessed = xcorr_scorer.preprocess_spectrum_for_xcorr(&spectrum);

// Preprocess library entry for XCorr
let lib_preprocessed = xcorr_scorer.preprocess_library_for_xcorr(&entry);

// Compute XCorr from preprocessed vectors (BLAS sdot)
let xcorr = SpectralScorer::xcorr_from_preprocessed(&spec_preprocessed, &lib_preprocessed);
```

### Top-3 Fragment Matching API

```rust
use osprey_scoring::top3_fragment_match_with_errors;

// Get both filter result and mass errors in one call
let (has_match, ms2_errors) = top3_fragment_match_with_errors(
    &entry.fragments,
    &spectrum.mzs,
    tolerance,
    unit,
);

// has_match: true if at least 1 of top-3 matched
// ms2_errors: signed mass errors for each matched fragment
```

### E-Value Calculation API

```rust
use osprey_scoring::batch::calculate_evalue_from_xcorr_distribution;

// Collect ALL XCorr scores in RT window
let all_xcorr_scores: Vec<f64> = candidate_spectra
    .iter()
    .map(|spec| SpectralScorer::xcorr_from_preprocessed(&spec_preprocessed, &lib_preprocessed))
    .collect();

// Calculate E-value from survival function
let evalue = calculate_evalue_from_xcorr_distribution(&all_xcorr_scores, best_xcorr);
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
Top-3 by intensity: y3, y4, y5

Step 1: Pre-filter (top-3 check)
  y3: 348.2 → binary search → 348.201 found ✓ → passes filter

Step 2: XCorr against ALL spectra in RT window (20% tolerance)
  700 spectra preprocessed (unit resolution: 2001 bins)
  Best XCorr: 0.42 at RT=25.3 min (scan 1234)
  All XCorr scores collected for E-value calculation

Step 3: Top-3 fragment errors at best XCorr spectrum (RT=25.3 min)
  y3: 348.2 → 348.201 ✓ (error: +2.9 ppm)
  y4: 461.3 → 461.299 ✓ (error: -2.2 ppm)
  y5: 576.3 → 576.305 ✓ (error: +8.7 ppm)
  → 3 MS2 mass errors for calibration

Step 4: E-value from XCorr survival function
  700 XCorr scores → fit survival function → E-value: 1.2e-8

Step 5: Isotope cosine (from MS1)
  Theoretical: [0.55, 0.30, 0.12, 0.03]  (from C41H68N10O15)
  Observed:    [0.52, 0.31, 0.13, 0.04]
  Isotope cosine: 0.998

Final scores:
  E-value: 1.2e-8 (used for target-decoy competition)
  XCorr: 0.42 (primary score)
  MS2 errors: [+2.9, -2.2, +8.7] ppm (for MS2 calibration)
  Isotope cosine: 0.998
```

## Main Search Scoring (Phase 3-4)

The main search uses a fundamentally different scoring approach than calibration. Instead of XCorr with E-value, the main search uses **fragment XIC co-elution analysis** with 45 features.

### Peak Detection

The main search uses **CWT consensus peak detection** to find peaks in fragment XICs. See [Peak Detection](05-peak-detection.md) for details.

### Spectral Scores at Apex (15 features)

At the detected peak apex, library fragments are matched against the observed spectrum. All intensity-based scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks.

| Score | Description |
|-------|-------------|
| `hyperscore` | X!Tandem-style hyperscore |
| `xcorr` | Comet-style cross-correlation (same as calibration, but at coelution-detected apex) |
| `dot_product` | Library cosine with sqrt intensity preprocessing |
| `dot_product_smz` | Library cosine with sqrt(intensity) * mz^2 preprocessing |
| `dot_product_top6/5/4` | Library cosine using only top N fragments |
| `dot_product_smz_top6/5/4` | SMZ cosine using only top N fragments |
| `fragment_coverage` | Fraction of library fragments matched |
| `sequence_coverage` | Fraction of peptide backbone covered by b/y ions |
| `consecutive_ions` | Longest consecutive b or y ion series |
| `explained_intensity` | Fraction of observed intensity explained by matches |
| `elution_weighted_cosine` | LibCosine at each scan within peak, weighted by reference XIC intensity^2, averaged |

### Co-Elution Features (11 features)

Pairwise Pearson correlations between fragment XICs within the peak boundaries:

| Feature | Description |
|---------|-------------|
| `fragment_coelution_sum` | Sum of all pairwise correlations |
| `fragment_coelution_min` | Minimum pairwise correlation (worst pair) |
| `fragment_coelution_max` | Maximum pairwise correlation (best pair) |
| `n_coeluting_fragments` | Number of fragments with positive mean correlation |
| `n_fragment_pairs` | Number of fragment pairs |
| `fragment_corr_0..5` | Per-fragment average correlation with all other fragments |

### Multi-Charge Consensus

After the main search, all charge states of the same peptide are forced to share the same peak RT and integration boundaries. See [Multi-Charge Consensus](06-multi-charge-consensus.md) for details.

### Feature Computation

All 45 features are computed via `compute_features_at_peak()`, a reusable function that takes a `FeatureComputeContext` (read-only references to library entry, XICs, spectra, calibration, etc.) and an `XICPeakBounds`. This function is called both during the initial per-precursor search and during multi-charge consensus re-scoring.

---

## Notes

1. **XCorr for everything**: XCorr is both the RT selection score and the primary calibration score
2. **Unit resolution bins for calibration**: Always 2001 bins regardless of data type (~50x faster for HRAM)
3. **E-value for competition**: E-value (from XCorr survival function) is used for target-decoy competition
4. **Top-3 for MS2 calibration**: MS2 mass errors come from top-3 fragment matching (binary search), not full fragment matching
5. **Isotope scoring**: Only available when MS1 spectra are present in the mzML file
6. **LibCosine in main search**: LibCosine (ppm fragment matching) is used during feature extraction (Phase 4), not during calibration
7. **Multiple scores**: All scores written to debug CSV for analysis
8. **CWT peak detection**: Both calibration and main search use CWT consensus peak detection for finding peaks in fragment XICs
