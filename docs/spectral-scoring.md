# Spectral Scoring

Osprey uses spectral similarity scoring to assess whether an observed spectrum matches its predicted library spectrum. This is the primary score used for target-decoy FDR computation.

## Overview

After ridge regression produces coefficient time series for each peptide, spectral scoring compares the observed spectrum at the peak apex to the library spectrum. This determines if the peptide's fragments are actually present in the data.

```
Workflow:
  1. Ridge regression → coefficient time series per peptide
  2. Peak detection → identify apex scan
  3. Spectral scoring → compare apex spectrum to library
  4. FDR computation → use spectral score for target-decoy competition
```

## Scoring Methods

### LibCosine (Primary Score)

LibCosine is a cosine similarity score with SMZ preprocessing (sqrt intensity × m/z²):

```
Preprocessing (per fragment):
  value = sqrt(intensity) × m/z²

Scoring:
  1. Match library fragments to observed peaks (within tolerance)
  2. Apply SMZ preprocessing to both
  3. L2 normalize both vectors
  4. Compute cosine similarity (dot product of normalized vectors)

Score range: 0-1 (1 = perfect match)
```

**Why SMZ preprocessing?**
- Square root dampens intensity differences (accounts for measurement variability)
- m/z² weighting emphasizes higher mass fragments (more informative)
- L2 normalization enables meaningful comparison regardless of total intensity

### XCorr (Secondary Score)

XCorr uses Comet-style preprocessing for cross-correlation scoring:

```
Preprocessing:
  1. Bin spectrum into discrete bins
  2. Apply sqrt transformation
  3. Windowing normalization (10 windows, normalize each to max=50)
  4. Sliding window subtraction (offset=75, removes local average)

Scoring:
  - Dot product between preprocessed observed and library spectra
  - Scaled by 0.005 (spectrum-centric scaling)
```

**When to use XCorr?**
- XCorr is more robust to noise and interfering ions
- May perform better in complex DIA spectra
- Currently available but not used as primary score

## Fragment Matching

Library fragments are matched to observed peaks using tolerance:

```
Da tolerance: 0.5 (unit resolution)
ppm tolerance: 20 (HRAM)

For each library fragment:
  Find closest observed peak within tolerance
  Record match if found
```

**Metrics computed:**
- `n_matched`: Number of library fragments with matches
- `n_library`: Total library fragments
- `fragment_coverage`: n_matched / n_library
- `explained_intensity`: Sum of matched intensities / total observed intensity

## Integration with FDR

Spectral scores are used for target-decoy FDR computation:

```python
def compute_fdr(targets, decoys):
    # Primary score is LibCosine
    for peptide in targets + decoys:
        apex_spectrum = get_apex_spectrum(peptide)
        score = lib_cosine(apex_spectrum, peptide.library)

    # Sort by score (higher = better)
    all_peptides = sort_by_score(targets + decoys)

    # Compute q-values using target-decoy competition
    for peptide in all_peptides:
        q_value = count_decoys_above / count_targets_above

    # Filter to FDR threshold
    return [p for p in targets if p.q_value <= fdr_threshold]
```

## Configuration

```yaml
# Spectral scoring is automatic when apex spectrum is available

# Fragment matching tolerance (based on resolution mode)
resolution_mode: unit  # Uses 0.5 Da
# OR
resolution_mode:
  HRAM:
    tolerance_ppm: 20.0  # Uses 20 ppm
```

## Feature Set

Spectral scoring populates these FeatureSet fields:

| Field | Source | Description |
|-------|--------|-------------|
| `dot_product` | LibCosine | Cosine similarity score (0-1) |
| `spectral_contrast_angle` | LibCosine | acos(LibCosine) in degrees |
| `hyperscore` | XCorr | Cross-correlation score |
| `fragment_coverage` | Matching | Fraction of library fragments matched |
| `explained_intensity` | Matching | Fraction of observed intensity explained |

## Implementation

Key files:
- `crates/osprey-scoring/src/lib.rs` - SpectralScorer, LibCosine, XCorr
- `crates/osprey/src/pipeline.rs` - Integration in compute_fdr_and_filter()

### SpectralScorer API

```rust
use osprey_scoring::SpectralScorer;

let scorer = SpectralScorer::new()
    .with_tolerance_da(0.5)      // Unit resolution
    .with_tolerance_ppm(20.0);   // HRAM

// Compute LibCosine score
let score = scorer.lib_cosine(&observed_spectrum, &library_entry);
println!("LibCosine: {}", score.lib_cosine);
println!("Matched {}/{} fragments", score.n_matched, score.n_library);

// Compute XCorr score
let score = scorer.xcorr(&observed_spectrum, &library_entry);
println!("XCorr: {}", score.xcorr);
```

## Example

```
Peptide: PEPTIDEK
Library fragments: y3(348.2), y4(461.3), y5(576.3), b3(357.2), b4(458.2)

Observed spectrum at apex RT=25.3 min:
  m/z: [348.2, 461.3, 576.4, 400.1, 550.2, ...]
  intensity: [1000, 800, 600, 200, 150, ...]

Fragment matching (0.5 Da tolerance):
  y3: 348.2 → 348.2 ✓ (error: 0.0 Da)
  y4: 461.3 → 461.3 ✓ (error: 0.0 Da)
  y5: 576.3 → 576.4 ✓ (error: 0.1 Da)
  b3: 357.2 → no match ✗
  b4: 458.2 → no match ✗

Metrics:
  n_matched: 3
  n_library: 5
  fragment_coverage: 0.60

LibCosine scoring:
  Library SMZ: [sqrt(100)×348.2², sqrt(80)×461.3², sqrt(60)×576.3²]
  Observed SMZ: [sqrt(1000)×348.2², sqrt(800)×461.3², sqrt(600)×576.4²]
  (normalized and dot product computed)
  LibCosine: 0.95
```

## Notes

1. **Score interpretation**: LibCosine > 0.9 is typically a good match
2. **Missing fragments**: Low fragment coverage may indicate:
   - Peptide not present
   - Co-eluting interference
   - Ion suppression
3. **Decoy scoring**: Decoys should have lower spectral scores on average
4. **Future work**: Combine multiple scores with machine learning (mokapot)
