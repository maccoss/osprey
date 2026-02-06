# FDR Control

Osprey uses two-level FDR (False Discovery Rate) control to ensure high-quality peptide identifications. The primary method is semi-supervised learning via Mokapot, which combines multiple scoring features into an optimal discriminant score.

## Overview

```
FDR Control Workflow:
  1. Extract 30 features per precursor (targets + decoys)
  2. Write PIN file (Percolator/Mokapot input format)
  3. Run Mokapot with semi-supervised learning
  4. Two-level FDR:
     - Run-level: Per-file FDR control
     - Experiment-level: Aggregate best precursor across replicates
  5. Output: q-values per precursor
```

## Key Terminology

| Term | Definition |
|------|------------|
| **Precursor** | Peptide + charge state (e.g., PEPTIDEK+2) |
| **Target** | Real peptide from the spectral library |
| **Decoy** | Reversed/shuffled peptide for FDR estimation |
| **q-value** | Minimum FDR at which a result would pass |
| **Run-level FDR** | FDR controlled within a single file |
| **Experiment-level FDR** | FDR controlled across all replicates |

## Two-Level FDR Strategy

### Step 1: Run-Level FDR

Each mzML file is processed independently:

1. **Ridge regression** produces coefficients for all candidate precursors
2. **Peak detection** finds the apex for each precursor's coefficient time series
3. **Feature extraction** computes 30 features at the apex spectrum
4. **Mokapot** combines features via linear SVM and computes run-level q-values

```
Per-file processing:
  ┌─────────────────────────────────────────────┐
  │ file1.mzML                                  │
  │   ├─ Ridge regression → coefficients        │
  │   ├─ Peak detection → apex RTs              │
  │   ├─ Feature extraction → 30 features       │
  │   └─ Mokapot → run_qvalue per precursor     │
  └─────────────────────────────────────────────┘
```

### Step 2: Experiment-Level FDR

If multiple replicates are provided, Osprey aggregates results:

1. **Aggregate by precursor** (peptide + charge, not just peptide)
2. **Keep best score** per precursor across all replicates
3. **Re-run Mokapot** on the aggregated set
4. **Assign experiment-level q-values**

```
Experiment-level aggregation:
  ┌─────────────────────────────────────────────┐
  │ Best precursor per peptide+charge:          │
  │   PEPTIDEK+2: file1 score=0.95 (best)       │
  │   PEPTIDEK+3: file2 score=0.87 (best)       │
  │   ANOTHERK+2: file1 score=0.91 (best)       │
  │   ...                                       │
  │                                             │
  │ Re-run Mokapot → experiment_qvalue          │
  └─────────────────────────────────────────────┘
```

### Single Replicate Handling

When only one file is provided, Step 2 is **skipped**:
- Run-level q-values are copied to experiment-level q-values
- No redundant Mokapot run
- Faster processing

## Mokapot Integration

### What is Mokapot?

Mokapot is a Python tool for semi-supervised learning in proteomics FDR control. It:
- Uses linear SVM (like Percolator) to combine features
- Learns to distinguish targets from decoys
- Provides calibrated q-values via target-decoy competition

### PIN File Format

Osprey writes a PIN (Percolator INput) file with 30 features:

```
SpecId  Label  ScanNr  ChargeState  peak_apex  peak_area  ...  Peptide  Proteins
psm_001    1    1234           2       0.85      12.3  ...  -.PEPTIDEK.-  PROT1
psm_002   -1    1234           2       0.12       1.1  ...  -.KEDITPEP.-  DECOY_PROT1
```

- `Label`: 1 for target, -1 for decoy
- `SpecId`: Unique identifier (file_precursorId)
- 30 feature columns (see Feature Set below)
- `Peptide`: Flanking.SEQUENCE.Flanking format
- `Proteins`: Protein accessions

### Mokapot Command

Osprey runs Mokapot with these settings:

```bash
mokapot input.pin \
  --dest_dir output/ \
  --train_fdr 0.01 \
  --test_fdr 0.01 \
  --max_iter 10 \
  --seed 42 \
  --num_workers 8 \
  --verbosity 1 \
  --save_models
```

| Flag | Description |
|------|-------------|
| `--train_fdr` | FDR threshold for training set |
| `--test_fdr` | FDR threshold for test set |
| `--max_iter` | Maximum SVM iterations |
| `--num_workers` | Parallel cross-validation workers (auto-detected, max 8) |
| `--save_models` | Save trained model for feature weight inspection |

### Feature Weights

Mokapot saves model weights to `mokapot.model.pkl`. Use the inspection script to see which features are most important:

```bash
python scripts/inspect_mokapot_weights.py mokapot.model.pkl
```

Output:
```
FEATURE WEIGHTS (ranked by absolute importance)
======================================================================
Rank   Feature                                  Weight     |Weight|
----------------------------------------------------------------------
1      dot_product                              0.8532      0.8532  ++++++++
2      xcorr                                    0.7891      0.7891  +++++++
3      peak_apex                                0.6234      0.6234  ++++++
4      rt_deviation_normalized                 -0.4521      0.4521  ----
5      hyperscore                               0.3987      0.3987  +++
...
```

Features with very small weights may be candidates for removal to speed up processing.

## Feature Set (30 Features)

### Chromatographic Features (12)

| Feature | Description | Range |
|---------|-------------|-------|
| `peak_apex` | Maximum coefficient value | 0-∞ |
| `peak_area` | Integrated area under curve | 0-∞ |
| `emg_fit_quality` | EMG fit R² | 0-1 |
| `peak_width` | FWHM in minutes | 0-∞ |
| `peak_symmetry` | Leading/trailing ratio | 0-∞ |
| `rt_deviation` | Observed - predicted RT (min) | -∞ to ∞ |
| `rt_deviation_normalized` | RT deviation / RT tolerance | -∞ to ∞ |
| `n_contributing_scans` | Scans with non-zero coefficient | 1-∞ |
| `coefficient_stability` | CV of coefficients near apex | 0-∞ |
| `peak_sharpness` | Boundary steepness | 0-∞ |
| `peak_prominence` | Apex / baseline ratio | 0-∞ |
| `modification_count` | Number of modifications | 0-∞ |

### Spectral Features (13)

| Feature | Description | Range |
|---------|-------------|-------|
| `hyperscore` | X!Tandem-style score | 0-∞ |
| `xcorr` | Comet-style cross-correlation | 0-∞ |
| `spectral_contrast_angle` | Normalized spectral angle | 0-1 |
| `dot_product` | LibCosine (sqrt preprocessing) | 0-1 |
| `dot_product_smz` | LibCosine (sqrt×mz² preprocessing) | 0-1 |
| `pearson_correlation` | Pearson intensity correlation | -1 to 1 |
| `spearman_correlation` | Spearman rank correlation | -1 to 1 |
| `fragment_coverage` | Matched/predicted fragments | 0-1 |
| `sequence_coverage` | Backbone coverage | 0-1 |
| `consecutive_ions` | Longest b/y ion run | 0-n |
| `base_peak_rank` | Base peak rank in predicted | 1-n |
| `top3_matches` | Top-3 fragments matched | 0-3 |
| `explained_intensity` | Explained/total intensity | 0-1 |

### Contextual Features (5)

| Feature | Description | Range |
|---------|-------------|-------|
| `n_competitors` | Candidates in regression | 1-∞ |
| `relative_coefficient` | Coefficient / sum(all) | 0-1 |
| `local_peptide_density` | Candidates per window | 1-∞ |
| `spectral_complexity` | Spectral complexity estimate | 0-∞ |
| `regression_residual` | Unexplained signal fraction | 0-1 |

## Target-Decoy Competition

### Decoy Generation

Osprey generates decoys using enzyme-aware sequence reversal:

```
Target:  K.PEPTIDEK.R  (tryptic)
Decoy:   K.KEDITPEP.R  (reversed, preserving cleavage sites)
```

Each target has a paired decoy with the same:
- Precursor m/z
- Charge state
- Fragment m/z values (but from reversed sequence)
- Library RT (same elution window)

### Competition Rules

For each target-decoy pair:
1. Both are scored independently
2. Higher score wins the competition
3. **Ties go to decoy** (conservative)
4. Only winners enter FDR calculation

### FDR Calculation

After competition, winners are sorted by score (descending):

```
FDR at position i = (decoy wins so far) / (target wins so far)

Position  Score  Winner   Decoy_Wins  Target_Wins  FDR
   1      0.98   target        0           1       0.00
   2      0.95   target        0           2       0.00
   3      0.92   decoy         1           2       0.50
   4      0.90   target        1           3       0.33
   5      0.88   target        1           4       0.25
   ...
```

Find the maximum number of targets where FDR ≤ threshold (e.g., 1%):
- Walk down the ranked list
- Track cumulative targets and decoys
- Find the deepest position where FDR ≤ threshold
- Return all targets up to that position

## Output

### Mokapot Output Files

| File | Description |
|------|-------------|
| `mokapot.psms.txt` | PSM-level results with q-values |
| `mokapot.peptides.txt` | Peptide-level results |
| `mokapot.model.pkl` | Trained model (with `--save_models`) |

### Result Fields

Each passing precursor has:

| Field | Description |
|-------|-------------|
| `run_qvalue` | q-value within the single file |
| `experiment_qvalue` | q-value across all replicates |
| `mokapot_score` | Combined discriminant score |
| `pep` | Posterior error probability |

## Configuration

```yaml
# FDR thresholds
run_fdr: 0.01           # 1% run-level FDR
experiment_fdr: 0.01    # 1% experiment-level FDR

# Mokapot settings
mokapot:
  train_fdr: 0.01
  test_fdr: 0.01
  max_iter: 10
  num_workers: 8        # Auto-detected if not specified
```

## Implementation

Key files:
- `crates/osprey-fdr/src/lib.rs` - FdrController, target-decoy competition
- `crates/osprey-fdr/src/mokapot.rs` - MokapotRunner, PIN file generation
- `crates/osprey/src/pipeline.rs` - Two-level FDR orchestration
- `scripts/inspect_mokapot_weights.py` - Feature weight analysis

### API Example

```rust
use osprey_fdr::{FdrController, MokapotRunner, PsmFeatures};

// Simple FDR control (no Mokapot)
let fdr = FdrController::new(0.01);  // 1% FDR
let passing = fdr.filter_with_competition(&targets, &decoys, |entry| entry.score);

// Mokapot integration
let runner = MokapotRunner::new()
    .with_train_fdr(0.01)
    .with_num_workers(8);

// Write PIN file
runner.write_pin(&psm_features, "output.pin")?;

// Run Mokapot
let results = runner.run("output.pin", "output_dir/")?;

// Results contain q-values
for result in results {
    println!("{}: q-value = {:.4}", result.psm_id, result.q_value);
}
```

## Troubleshooting

### Mokapot Not Found

```
Error: Failed to run mokapot: No such file or directory
```

Install Mokapot:
```bash
pip install mokapot
```

### Too Few Targets Passing

If FDR filtering is too aggressive:
1. Check calibration quality (RT tolerance may be too tight)
2. Verify decoy generation is working (should be ~50% decoy wins at random)
3. Inspect feature weights - some features may be hurting performance
4. Try increasing `train_fdr` for Mokapot training

### Slow Mokapot Processing

Mokapot cross-validation can be slow with many PSMs:
1. Ensure `--num_workers` is set (Osprey auto-detects CPUs)
2. Consider reducing feature set to important features only
3. Use `--save_models` to inspect which features matter

## References

- Mokapot: https://github.com/wfondrie/mokapot
- Percolator: https://github.com/percolator/percolator
- Target-decoy approach: Elias & Gygi (2007) Nature Methods
