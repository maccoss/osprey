# FDR Control

Osprey uses two-level FDR (False Discovery Rate) control to ensure high-quality peptide identifications. The primary method is semi-supervised learning via Mokapot, which combines multiple scoring features into an optimal discriminant score.

## Overview

```
FDR Control Workflow:
  1. Extract 37 features per precursor (targets + decoys)
  2. Write PIN files (one per mzML file)
  3. Run Mokapot CLI with two-step analysis:
     - Step 1: Train joint model, report per-file q-values
     - Step 2: Reuse model, report experiment-level q-values
  4. Output: run-level and experiment-level q-values per precursor
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

Osprey uses Mokapot's CLI with `--save_models` and `--load_models` for efficient two-step analysis:

### Step 1: Run-Level FDR (Joint Model Training)

All mzML files are processed together to train a single joint model:

```bash
mokapot file1.pin file2.pin file3.pin \
  --dest_dir mokapot/run_level/ \
  --save_models \
  --train_fdr 0.01 \
  --test_fdr 0.01
```

**Without `--aggregate`**, Mokapot:
- Trains ONE model on PSMs from all files
- Reports separate results for each file
- Saves the trained model for reuse

```
mokapot/run_level/
├── file1.mokapot.psms.txt    # Run-level q-values for file1
├── file2.mokapot.psms.txt    # Run-level q-values for file2
├── file3.mokapot.psms.txt    # Run-level q-values for file3
└── mokapot.model             # Trained model for Step 2
```

### Step 2: Experiment-Level FDR (Model Reuse)

The saved model is reused with `--aggregate` for experiment-level results:

```bash
mokapot file1.pin file2.pin file3.pin \
  --dest_dir mokapot/experiment_level/ \
  --load_models mokapot/run_level/mokapot.model \
  --aggregate \
  --test_fdr 0.01
```

**With `--aggregate`**, Mokapot:
- Loads the model from Step 1 (no retraining)
- Combines all files for experiment-level FDR
- Reports a single combined result

```
mokapot/experiment_level/
└── mokapot.psms.txt    # Experiment-level q-values
```

### Single File Handling

When only one mzML file is provided:
- Only Step 1 runs (no aggregation needed)
- Run-level q-values are used as experiment-level q-values
- Faster processing

## Mokapot Integration

### What is Mokapot?

Mokapot is a Python tool for semi-supervised learning in proteomics FDR control. It:
- Uses linear SVM (like Percolator) to combine features
- Learns to distinguish targets from decoys
- Provides calibrated q-values via target-decoy competition

### PIN File Format

Osprey writes one PIN file per mzML file with 37 features:

```
SpecId  Label  ScanNr  ChargeState  peak_apex  peak_area  ...  Peptide  Proteins
psm_001    1    1234           2       0.85      12.3  ...  -.PEPTIDEK.-  PROT1
psm_002   -1    1234           2       0.12       1.1  ...  -.KEDITPEP.-  DECOY_PROT1
```

- `Label`: 1 for target, -1 for decoy
- `SpecId`: Unique identifier for the PSM
- 37 feature columns (see Feature Set below)
- `Peptide`: Flanking.SEQUENCE.Flanking format
- `Proteins`: Protein accessions

### Mokapot CLI Options

Osprey uses these CLI options:

| Flag | Description |
|------|-------------|
| `--dest_dir` | Output directory for results |
| `--train_fdr` | FDR threshold for training set (default: 0.01) |
| `--test_fdr` | FDR threshold for test set (default: 0.01) |
| `--max_iter` | Maximum SVM iterations (default: 10) |
| `--seed` | Random seed for reproducibility (default: 42) |
| `--max_workers` | Parallel cross-validation workers (auto-detected, max 8) |
| `--subset_max_train` | Max PSMs for training (memory-aware, auto-calculated) |
| `--save_models` | Save trained model for reuse |
| `--load_models` | Load a previously trained model |
| `--aggregate` | Combine all files for experiment-level FDR |

### Memory-Aware Training

For large datasets, Osprey automatically limits the number of PSMs used for training via `--subset_max_train`. This prevents mokapot from running out of memory:

- Osprey detects available system memory at runtime
- Estimates ~50 KB per PSM for training memory usage
- Reserves 2 GB for system overhead
- Uses 70% of remaining memory as a safety margin
- If the calculated limit exceeds 1 million PSMs, no limit is applied

The subset is used only for model training. All PSMs are still scored and assigned q-values using the trained model.

### Feature Weights

Mokapot saves model weights to `mokapot.model`. Use the inspection script to see which features are most important:

```bash
python scripts/inspect_mokapot_weights.py mokapot/run_level/mokapot.model
```

Output:
```
FEATURE WEIGHTS (ranked by absolute importance)
======================================================================
Rank   Feature                                  Weight     |Weight|
----------------------------------------------------------------------
1      elution_weighted_cosine                  0.8532      0.8532  ++++++++
2      xcorr                                    0.7891      0.7891  +++++++
3      peak_apex                                0.6234      0.6234  ++++++
4      rt_deviation                            -0.4521      0.4521  ----
5      median_polish_cosine                     0.3987      0.3987  +++
...
```

## Feature Set (37 Features)

All intensity-based spectral similarity scores include ALL library fragments
within the spectrum's mass range, using 0 intensity for unmatched peaks.

### Ridge Regression Features (8)

| Feature | Description | Range |
|---------|-------------|-------|
| `peak_apex` | Maximum coefficient value | 0+ |
| `peak_area` | Integrated area under curve | 0+ |
| `peak_width` | FWHM from Tukey median polish elution profile | 0+ |
| `coefficient_stability` | CV of coefficients near apex | 0+ |
| `relative_coefficient` | Coefficient / sum(all) | 0-1 |
| `explained_intensity` | Explained/total intensity | 0-1 |
| `signal_to_noise` | Peak apex / noise estimate | 0+ |
| `xic_signal_to_noise` | S/N from best-correlated fragment XIC | 0+ |

### Spectral Matching - Mixed (2)

| Feature | Description | Range |
|---------|-------------|-------|
| `xcorr` | Comet-style cross-correlation | 0+ |
| `consecutive_ions` | Longest b/y ion run | 0-n |

### Spectral Matching - Deconvoluted (2)

| Feature | Description | Range |
|---------|-------------|-------|
| `xcorr_deconv` | XCorr on deconvoluted spectrum | 0+ |
| `consecutive_ions_deconv` | Consecutive ions on deconvoluted spectrum | 0-n |

### RT Deviation (1)

| Feature | Description | Range |
|---------|-------------|-------|
| `rt_deviation` | Observed - predicted RT (min) | -inf to inf |

### Fragment Co-elution (9)

| Feature | Description | Range |
|---------|-------------|-------|
| `fragment_coelution_sum` | Sum of per-fragment correlations | 0+ |
| `fragment_coelution_min` | Minimum fragment correlation | -1 to 1 |
| `n_coeluting_fragments` | Fragments with positive correlation | 0-6 |
| `fragment_corr_0..5` | Per-fragment correlation (ranked by library intensity) | -1 to 1 |

### Elution-Weighted Similarity (1)

| Feature | Description | Range |
|---------|-------------|-------|
| `elution_weighted_cosine` | LibCosine at each scan, weighted by coef^2, averaged | 0-1 |

### Mass Accuracy (3)

| Feature | Description | Range |
|---------|-------------|-------|
| `mass_accuracy_deviation_mean` | Mean signed mass error (ppm) | -inf to inf |
| `abs_mass_accuracy_deviation_mean` | Mean |mass error| (ppm) | 0+ |
| `mass_accuracy_std` | Std dev of mass errors (ppm) | 0+ |

### Percolator-style (6)

| Feature | Description | Range |
|---------|-------------|-------|
| `abs_rt_deviation` | |RT deviation| | 0+ |
| `peptide_length` | Number of amino acids | 5-50 |
| `missed_cleavages` | Missed enzymatic cleavages | 0+ |
| `ln_num_candidates` | ln(candidates in regression) | 0+ |
| `coef_zscore` | Z-score of coefficient at apex | -inf to inf |
| `coef_zscore_mean` | Mean z-score across peak | -inf to inf |

### MS1 Features (2)

| Feature | Description | Range |
|---------|-------------|-------|
| `ms1_precursor_coelution` | Coefficient vs MS1 XIC correlation | -1 to 1 |
| `ms1_isotope_cosine` | Observed vs theoretical isotope envelope | 0-1 |

### Tukey Median Polish (3)

| Feature | Description | Range |
|---------|-------------|-------|
| `median_polish_cosine` | Row effects vs library, sqrt preprocessing | 0-1 |
| `median_polish_rsquared` | Additive model R^2 in sqrt space | 0-1 |
| `median_polish_residual_ratio` | sum|obs-pred| / sum(obs), linear | 0+ |

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

### Directory Structure

```
output_dir/
└── mokapot/
    ├── file1.pin              # PIN file for file1
    ├── file2.pin              # PIN file for file2
    ├── run_level/             # Step 1 outputs
    │   ├── file1.mokapot.psms.txt
    │   ├── file2.mokapot.psms.txt
    │   └── mokapot.model
    └── experiment_level/      # Step 2 outputs
        └── mokapot.psms.txt
```

### Result Fields

Each passing precursor has:

| Field | Description |
|-------|-------------|
| `run_qvalue` | q-value within the single file |
| `experiment_qvalue` | q-value across all replicates |
| `mokapot_score` | Combined discriminant score |
| `pep` | Posterior error probability |

### Terminal Output

Osprey reports precursors and peptides at each level:

```
=== Per-file results (1% FDR) ===
  file1.mzML: 1234 precursors, 987 peptides
  file2.mzML: 1456 precursors, 1123 peptides

=== Experiment-level results (1% FDR) ===
  Experiment: 1567 precursors, 1234 peptides
```

## Configuration

```yaml
# FDR thresholds
run_fdr: 0.01           # 1% run-level FDR
experiment_fdr: 0.01    # 1% experiment-level FDR
```

CLI options are passed directly to Mokapot. Osprey auto-detects the number of workers (capped at 8).

## Implementation

Key files:
- `crates/osprey-fdr/src/lib.rs` - FdrController, target-decoy competition
- `crates/osprey-fdr/src/mokapot.rs` - MokapotRunner, PIN file generation, CLI execution
- `crates/osprey/src/pipeline.rs` - Two-level FDR orchestration
- `scripts/inspect_mokapot_weights.py` - Feature weight analysis

### API Example

```rust
use osprey_fdr::{MokapotRunner, PsmFeatures};
use std::collections::HashMap;

// Build mokapot runner
let runner = MokapotRunner::new()
    .with_test_fdr(0.01)
    .with_num_workers(8);

// Write PIN files (one per mzML file)
let psms_by_file: HashMap<String, Vec<PsmFeatures>> = collect_psm_features();
let pin_files = runner.write_pin_files(&psms_by_file, "mokapot/")?;

// Run two-step analysis
let (per_file_results, experiment_results) =
    runner.run_two_step_analysis(&pin_files, "mokapot/")?;

// Per-file results (run-level q-values)
for (file_name, results) in per_file_results {
    println!("{}: {} PSMs passing FDR", file_name, results.len());
    for result in results {
        println!("  {}: q={:.4}", result.psm_id, result.q_value);
    }
}

// Experiment-level results
println!("Experiment: {} PSMs passing FDR", experiment_results.len());
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
1. Ensure `--max_workers` is set (Osprey auto-detects CPUs, max 8)
2. Consider reducing feature set to important features only
3. Use feature weight inspection to identify low-importance features

### Model Loading Errors

If Step 2 fails to load the model:
1. Check that Step 1 completed successfully
2. Verify `mokapot.model` exists in `run_level/` directory
3. Ensure mokapot version matches between steps

## References

- Mokapot: https://github.com/wfondrie/mokapot
- Mokapot CLI documentation: https://mokapot.readthedocs.io/en/latest/cli.html
- Percolator: https://github.com/percolator/percolator
- Target-decoy approach: Elias & Gygi (2007) Nature Methods
