# FDR Control

Osprey uses two-level FDR (False Discovery Rate) control to ensure high-quality peptide identifications. Three FDR methods are available: a native Percolator implementation (default), external Mokapot, and simple target-decoy competition.

## Overview

```
FDR Control Workflow:
  1. Extract features per precursor (37 for regression, 45 for coelution)
  2. Dispatch to selected FDR method:
     - Percolator (default): native linear SVM with cross-validation
     - Mokapot: write PIN files, run external mokapot CLI
     - Simple: direct target-decoy FDR on a single feature
  3. Two-level q-values: per-run (within file) and experiment-level (across files)
  4. Posterior error probabilities (PEP) via KDE + isotonic regression
  5. Blib output: only precursors passing experiment-level FDR
```

## Key Terminology

| Term | Definition |
|------|------------|
| **Precursor** | Peptide + charge state (e.g., PEPTIDEK+2) |
| **Target** | Real peptide from the spectral library |
| **Decoy** | Reversed peptide for FDR estimation (enzyme-aware reversal) |
| **base_id** | `entry_id & 0x7FFFFFFF` — links each target to its paired decoy |
| **q-value** | Minimum FDR at which a result would pass |
| **PEP** | Posterior error probability — probability that a single result is incorrect |
| **Run-level FDR** | FDR controlled within a single file |
| **Experiment-level FDR** | FDR controlled across all replicates |

## Target-Decoy Strategy

### Decoy Generation

Osprey generates decoys using enzyme-aware sequence reversal:

```
Target:  K.PEPTIDEK.R  (tryptic)
Decoy:   K.KEDITPEP.R  (reversed, preserving cleavage sites)
```

Each target has a paired decoy with the same precursor m/z, charge state, and library RT. Decoys are identified by `entry_id | 0x80000000` (high bit set), so `entry_id & 0x7FFFFFFF` yields the shared base_id linking a target-decoy pair.

### Target-Decoy Competition (TDC)

The target-decoy competition approach assumes that decoys model the distribution of incorrect matches. Targets are a mixture of correct and incorrect matches. By competing each target against its paired decoy, incorrect targets are removed at roughly the same rate as decoys, and the remaining decoys estimate the residual false discovery rate:

```
FDR = (decoys + 1) / targets       (conservative formula)
```

The +1 in the numerator provides a conservative estimate (Levitsky et al., 2017).

### When TDC Happens

The timing of TDC differs between FDR methods:

| FDR Method | When TDC occurs | What enters the FDR engine |
|------------|-----------------|---------------------------|
| **Percolator** (native) | After SVM scoring, inside Percolator | Both targets and decoys — Percolator does internal paired competition |
| **Mokapot** (external) | Before PIN writing, using best feature by ROC AUC | Only competition winners — mokapot assumes upstream TDC |
| **Simple** | After scoring | Competition on final scores |

This distinction is critical: mokapot's FDR estimates are only valid when it receives pre-competed data (winners only). The native Percolator handles paired competition internally within each cross-validation fold.

### Feature Selection for Pre-Competition (Mokapot)

When pre-competing for mokapot, the best-separating feature is selected by **ROC AUC** (area under the receiver operating characteristic curve), computed via the Mann-Whitney U statistic:

```
AUC = P(target_score > decoy_score)
```

- AUC > 0.5: higher values are more target-like
- AUC < 0.5: lower values are more target-like (direction flipped)
- The feature with highest |AUC - 0.5| is selected

For each target-decoy pair, the entry with the better value on the selected feature wins. Ties go to the decoy (conservative).

## Native Percolator (Default)

### Algorithm Overview

The native Percolator implements the semi-supervised learning approach of Käll et al. (2007), using a linear SVM trained iteratively via cross-validation. Both targets and their paired decoys enter the algorithm — no upstream competition is performed.

```
Native Percolator Workflow:
  1. Standardize all features to zero mean, unit variance
  2. Assign folds: group by target peptide (base_id), keeping pairs together
  3. Find best initial feature (test both ascending and descending directions)
  4. For each CV fold (trained in parallel):
     a. Grid search for SVM cost parameter C (once, on initial labels)
     b. Iterative SVM training (up to 10 iterations):
        - Select positive training set: targets passing FDR threshold
        - Train linear SVM on selected targets + all decoys
        - Score training set, check for improvement
        - Stop after 2 consecutive non-improvements
     c. Score held-out test set with best model
  5. Calibrate scores between folds (Granholm et al. 2012)
  6. Compute PEP via KDE + isotonic regression
  7. Compute per-run and experiment-level q-values
```

### Fold Assignment

Cross-validation fold assignment must keep target-decoy pairs in the same fold. If pairs are split, unpaired targets in the training set auto-win competition, inflating the positive training set and causing the SVM to become too permissive.

Fold assignment groups all entries by their **target peptide** via `base_id`:

1. Build a mapping from `base_id` → target peptide sequence
2. Each entry (target or decoy) looks up its target peptide via `base_id`
3. All entries sharing a target peptide (all charge states + their decoys) → same fold
4. Round-robin assignment across folds, sorted alphabetically for determinism

```
Example with 3 folds:
  PEPTIDEK z=2 (target, base_id=1)  → group "PEPTIDEK" → fold 0
  KEDITPEP z=2 (decoy,  base_id=1)  → group "PEPTIDEK" → fold 0
  PEPTIDEK z=3 (target, base_id=2)  → group "PEPTIDEK" → fold 0
  KEDITPEP z=3 (decoy,  base_id=2)  → group "PEPTIDEK" → fold 0
  ANOTHERONE z=2 (target, base_id=3) → group "ANOTHERONE" → fold 1
  ENOREHTONA z=2 (decoy, base_id=3) → group "ANOTHERONE" → fold 1
```

This matches the approach of both mokapot (groups by spectrum) and C++ Percolator (groups by scan number). In precursor-centric DIA, the equivalent of a spectrum is the target-decoy pair identified by `base_id`.

### Initial Feature Selection

Before SVM training, the single best-discriminating feature is found by testing **both ascending and descending** directions for each feature:

- For each feature: score all entries using the raw standardized value
- Also test negated values (descending — lower is better)
- Count targets passing the FDR threshold after paired competition
- The feature + direction with the most passing targets is selected

This ensures features like `rt_deviation` (lower is better) are not missed.

### Iterative SVM Training

For each cross-validation fold:

1. **Select positive training set**: targets passing FDR threshold on current scores (initial feature or previous SVM scores). If fewer than 50 targets pass, progressively relax to 5%, 10%, 25%, 50% FDR.

2. **Build SVM training set**: selected targets (positive class) + all decoys (negative class).

3. **Grid search for C** (first iteration only): 6 values × 3 inner folds = 18 SVMs. The best C is reused for all subsequent iterations.

4. **Train linear SVM**: L2-regularized linear SVM with dual coordinate descent.

5. **Score and evaluate**: score the entire training set, count targets passing FDR. Track the best model across iterations.

6. **Convergence**: stop after 2 consecutive iterations without improvement in passing target count.

### Score Calibration Between Folds

Each fold's SVM produces scores on a different scale. The Granholm et al. (2012) method normalizes them with a linear transform per fold:

- Score at the FDR threshold → 0
- Median decoy score → -1

This ensures scores from different folds are comparable for the final q-value computation.

### Posterior Error Probability (PEP)

After score calibration, PEP is estimated on competition winners:

1. Fit target and decoy score distributions using kernel density estimation (KDE)
2. PEP = P(decoy | score) = p(score | decoy) × π₀ / p(score)
3. Apply isotonic regression to enforce monotonicity (PEP decreases as score increases)

### Q-Value Computation

Q-values are computed at four levels, all using the conservative `(decoys + 1) / targets` formula:

| Level | Scope |
|-------|-------|
| Per-run precursor | Within each file, per precursor |
| Per-run peptide | Within each file, best precursor per peptide |
| Experiment precursor | Across all files, per precursor |
| Experiment peptide | Across all files, best precursor per peptide |

For single-file runs, experiment-level q-values equal run-level q-values (no separate aggregation step).

### Dual Precursor + Peptide FDR

A precursor must pass FDR at **both** the precursor level (modified_sequence + charge) and the peptide level (modified_sequence only). This is enforced by taking the max of the two q-values:

```
effective_qvalue = max(precursor_qvalue, peptide_qvalue)
```

This is applied at both run and experiment levels:

```
run_qvalue        = max(run_precursor_qvalue,        run_peptide_qvalue)
experiment_qvalue = max(experiment_precursor_qvalue,  experiment_peptide_qvalue)
```

**Examples:**

| Precursor q | Peptide q | Effective q | Result at 1% FDR |
|-------------|-----------|-------------|-------------------|
| 0.005 | 0.015 | 0.015 | Rejected (peptide fails) |
| 0.015 | 0.005 | 0.015 | Rejected (precursor fails) |
| 0.003 | 0.008 | 0.008 | Accepted (both pass) |

**Why dual control matters**: Without peptide-level FDR, multiple charge states of the same peptide can each pass precursor-level FDR independently, even if the underlying peptide identification is weak. Peptide-level FDR ensures the peptide itself is confidently identified before any of its charge states are reported.

**Implementation by FDR method:**

- **Native Percolator**: `run_percolator()` computes all four q-value levels internally. When mapping results back to entries, the max is applied directly.
- **Mokapot**: `mokapot.psms.txt` provides precursor-level q-values; `mokapot.peptides.txt` provides peptide-level q-values. Both files are parsed and combined via max. The peptide file uses `-.PEPTIDEK.-` format with flanking characters that are stripped for matching against `modified_sequence`.

### Configuration

```yaml
fdr_method: percolator          # Default
run_fdr: 0.01                   # 1% run-level FDR threshold
experiment_fdr: 0.01            # 1% experiment-level FDR threshold
```

The native Percolator uses these defaults:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_iterations` | 10 | Maximum SVM training iterations per fold |
| `n_folds` | 3 | Cross-validation folds |
| `seed` | 42 | Random seed for reproducibility |
| `c_values` | [0.001, 0.01, 0.1, 1.0, 10.0, 100.0] | Grid search C values |
| `train_fdr` | 0.01 | FDR threshold for positive training set selection |

## Mokapot (External)

### Overview

Mokapot is an external Python tool that implements the same Percolator algorithm. Osprey writes PIN files and invokes the mokapot CLI. Mokapot assumes that target-decoy competition has already been performed — PIN files contain only competition winners.

### Two-Level FDR Strategy

Osprey uses mokapot's CLI with `--save_models` and `--load_models` for two-step analysis:

#### Step 1: Run-Level FDR (Joint Model Training)

All mzML files are processed together to train a single joint model:

```bash
mokapot file1.pin file2.pin file3.pin \
  --dest_dir mokapot/run_level/ \
  --save_models \
  --train_fdr 0.01 \
  --test_fdr 0.01
```

Without `--aggregate`, mokapot:
- Trains ONE model on PSMs from all files
- Reports separate results for each file
- Saves the trained model for reuse

#### Step 2: Experiment-Level FDR (Model Reuse)

The saved model is reused with `--aggregate` for experiment-level results:

```bash
mokapot file1.pin file2.pin file3.pin \
  --dest_dir mokapot/experiment_level/ \
  --load_models mokapot/run_level/mokapot.model \
  --aggregate \
  --test_fdr 0.01
```

With `--aggregate`, mokapot:
- Loads the model from Step 1 (no retraining)
- Combines all files for experiment-level FDR
- Reports a single combined result

#### Single File Handling

When only one mzML file is provided:
- Only Step 1 runs (no aggregation needed)
- Run-level q-values are used as experiment-level q-values

### PIN File Format

Osprey writes one PIN file per mzML file. The number of feature columns depends on the search mode: 37 for regression, 45 for coelution. PIN files contain **only competition winners** (pre-competed using the best feature by ROC AUC).

```
SpecId  Label  ScanNr  ChargeState  peak_apex  peak_area  ...  Peptide  Proteins
psm_001    1    1234           2       0.85      12.3  ...  -.PEPTIDEK.-  PROT1
psm_002   -1    5678           3       0.12       1.1  ...  -.KEDITPEP.-  DECOY_PROT1
```

- `Label`: 1 for target, -1 for decoy
- `SpecId`: Unique identifier for the PSM
- `Peptide`: Flanking.SEQUENCE.Flanking format
- `Proteins`: Protein accessions

### Mokapot CLI Options

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

For large datasets, Osprey automatically limits the number of PSMs used for training via `--subset_max_train`:

- Detects available system memory at runtime
- Estimates ~50 KB per PSM for training memory usage
- Reserves 2 GB for system overhead, uses 70% of remaining as safety margin
- If the calculated limit exceeds 1 million PSMs, no limit is applied

The subset is used only for model training. All PSMs are still scored using the trained model.

## Simple FDR

The simplest method applies target-decoy competition directly on the best single feature (selected by ROC AUC) without SVM reranking. Useful as a baseline or when the feature set is small.

## Feature Sets

All intensity-based spectral similarity scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks.

### Regression Mode (37 Features)

#### Ridge Regression Features (8)

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

### Coelution Mode (45 Features)

The coelution search mode does not use ridge regression. Instead, it extracts fragment XICs and scores based on pairwise correlations and spectral matching.

| Category | Count | Features |
|----------|-------|----------|
| Pairwise coelution | 11 | fragment_coelution_sum/min/max, n_coeluting_fragments, n_fragment_pairs, fragment_corr_0..5 |
| Peak shape | 7 | peak_apex, peak_area, peak_width, peak_symmetry, signal_to_noise, n_scans, peak_sharpness |
| Spectral at apex | 15 | hyperscore, xcorr, dot_product, dot_product_smz, dot_product_top6/5/4, dot_product_smz_top6/5/4, fragment_coverage, sequence_coverage, consecutive_ions, explained_intensity, elution_weighted_cosine |
| Mass accuracy | 3 | mass_accuracy_deviation_mean, abs_mass_accuracy_deviation_mean, mass_accuracy_std |
| RT deviation | 2 | rt_deviation, abs_rt_deviation |
| MS1 | 2 | ms1_precursor_coelution, ms1_isotope_cosine |
| Peptide properties | 2 | peptide_length, missed_cleavages |
| Tukey median polish | 3 | median_polish_cosine, median_polish_rsquared, median_polish_residual_ratio |

**Key differences from regression mode**: No coefficient_stability, relative_coefficient, xic_signal_to_noise, xcorr_deconv, consecutive_ions_deconv, ln_num_candidates, coef_zscore, coef_zscore_mean. Added hyperscore, dot_product variants (top6/5/4), fragment_coelution_max, n_fragment_pairs, peak_symmetry, peak_sharpness, n_scans, fragment_coverage, sequence_coverage.

## Feature Weight Inspection

Both FDR methods produce SVM/model weights. Use the inspection script:

```bash
python scripts/inspect_mokapot_weights.py mokapot/run_level/mokapot.model
```

For the native Percolator, weights per fold are logged to the console during the run.

## Output

### Result Fields

Each passing precursor has:

| Field | Description |
|-------|-------------|
| `run_qvalue` | max(run precursor q, run peptide q) within the single file |
| `experiment_qvalue` | max(experiment precursor q, experiment peptide q) across all replicates |
| `score` | Combined SVM discriminant score |
| `pep` | Posterior error probability |

### Multi-File Observation Propagation

After experiment-level FDR determines which precursors pass, **all per-file target observations** for those passing precursors are included in the blib output. This is critical for multi-file experiments (replicates, GPF runs):

```
FDR Pipeline Flow (Multi-File):

  1. Train model on all data (all files combined)
  2. Score all entries with trained model
  3. Compute run-level q-values per file (precursor + peptide level)
  4. Select best observation per precursor across experiment
  5. Compute experiment-level q-values (precursor + peptide level)
  6. Apply max(precursor, peptide) at both run and experiment levels
  7. Determine passing precursors: experiment_qvalue <= threshold
  8. Include ALL per-file observations for passing precursors in output
  9. Propagate best experiment_qvalue to all observations
```

**Why this matters:**
- Each file gets its own RT boundaries for a passing precursor
- Skyline uses per-file peak boundaries for quantification across replicates
- A precursor detected in 3 out of 5 files produces 3 blib entries with the same experiment_qvalue
- The best experiment_qvalue is propagated to all observations of that precursor

**Example:**

```
Precursor PEPTIDEK+2 detected in 3 files:
  file1.mzML: score=2.5, experiment_qvalue=0.003  ← best observation
  file2.mzML: score=1.8, experiment_qvalue=0.012
  file3.mzML: score=2.1, experiment_qvalue=0.007

Experiment-level competition uses best (score=2.5):
  → experiment_qvalue = 0.003

All 3 observations written to blib with experiment_qvalue = 0.003
  Each with its own apex_rt and peak boundaries from that file
```

### Terminal Output

```
=== Per-file results (1% FDR) ===
  file1.mzML: 1234 precursors, 987 peptides
  file2.mzML: 1456 precursors, 1123 peptides

=== Experiment-level results (1% FDR) ===
  Experiment: 1567 precursors, 1234 peptides
```

The per-file counts reflect dual precursor+peptide FDR (max of both q-values). The experiment-level count is always <= the union of per-file counts because only the best observation per precursor enters experiment-level competition.

## Configuration

```yaml
fdr_method: percolator          # percolator (default), mokapot, or simple
run_fdr: 0.01                   # 1% run-level FDR
experiment_fdr: 0.01            # 1% experiment-level FDR
```

## Implementation

Key files:

| File | Description |
|------|-------------|
| `crates/osprey-fdr/src/percolator.rs` | Native Percolator: SVM training, fold assignment, q-values, PEP |
| `crates/osprey-fdr/src/mokapot.rs` | MokapotRunner, PIN file generation, CLI execution |
| `crates/osprey-fdr/src/lib.rs` | FdrController, target-decoy competition |
| `crates/osprey-ml/src/svm.rs` | Linear SVM (dual CD), grid search, feature standardization |
| `crates/osprey-ml/src/pep.rs` | PEP estimation via KDE + isotonic regression |
| `crates/osprey/src/pipeline.rs` | FDR orchestration, pre-competition for mokapot, ROC AUC feature selection |
| `scripts/inspect_mokapot_weights.py` | Feature weight analysis |

## Troubleshooting

### Mokapot Not Found

```
Error: Failed to run mokapot: No such file or directory
```

Install Mokapot: `pip install mokapot`

### Too Few Targets Passing

If FDR filtering is too aggressive:
1. Check calibration quality (RT tolerance may be too tight)
2. Verify decoy generation is working (should be ~50% decoy wins at random)
3. Inspect feature weights — some features may be hurting performance
4. Try increasing `train_fdr` for training

### Inflated FDR Results

If FDR results are significantly higher than expected:
1. Verify that target-decoy pairs are being kept in the same fold (check Percolator logs)
2. For mokapot: ensure pre-competition is running before PIN writing
3. Check that the best separating feature has high ROC AUC (logged at competition time)

## References

- Percolator: Käll L, Canterbury J, Weston J, Noble WS, MacCoss MJ. Nat Methods. 2007;4(11):923-925.
- Score calibration: Granholm V, Noble WS, Käll L. BMC Bioinformatics. 2012;13 Suppl 16:S3.
- Mokapot: Fondrie WE, Noble WS. J Proteome Res. 2021;20(4):1966-1971.
- Target-decoy approach: Elias JE, Gygi SP. Nat Methods. 2007;4(3):207-214.
- Conservative FDR: Levitsky LI, et al. J Proteome Res. 2017;16(2):689-694.
