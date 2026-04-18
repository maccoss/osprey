# FDR Control

Osprey uses two-level FDR (False Discovery Rate) control to ensure high-quality peptide identifications. Three FDR methods are available: a native Percolator implementation (default), external Mokapot, and simple target-decoy competition.

## Overview

```
FDR Control Workflow:
  1. Extract features per precursor (21 PIN features written to per-file Parquet cache)
  2. Convert to FdrEntry stubs (~128 bytes each, Arc<str>-interned peptide sequences)
  3. Dispatch to selected FDR method:
     - Percolator (default): native linear SVM with cross-validation
       (loads PIN features on-demand from Parquet)
     - Mokapot: write PIN files, run external mokapot CLI
     - Simple: direct target-decoy FDR on a single feature
  4. Two-level q-values: per-run (within file) and experiment-level (across files)
  5. Posterior error probabilities (PEP) via KDE + isotonic regression
  6. Compact stubs: drop non-passing entries to free ~21 GB (240-file experiments)
  7. Cross-run reconciliation + second-pass FDR
  8. Protein FDR (optional): picked-protein competition using peptide SVM discriminant (Savitski 2015)
  9. Blib output: load lightweight plan entries from Parquet projection
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

**CRITICAL INVARIANT**: Target-decoy pairs and same-peptide charge states must ALWAYS stay together in any data partitioning — fold assignment, subsampling, or any other split. If pairs are separated:
- Unpaired targets in a training fold auto-win competition, inflating the positive training set
- The SVM becomes too permissive, silently corrupting FDR estimates
- This invariant applies to ALL cross-validation code: `percolator.rs`, `calibration_ml.rs`, and any future partitioning

### Training Set Subsampling

For large datasets (millions of entries), SVM training can be very slow (O(n²) with the number of training examples). Following The et al. (2016, PMC5059416), entries can be subsampled before fold splitting:

1. Group all entries by target peptide sequence (via `base_id`)
2. Randomly sample N peptide groups (keeping all entries in each group: target + decoy + all charge states)
3. Split the subsampled set into cross-validation folds
4. Train the SVM on the subsampled folds
5. Score ALL original entries with the trained model

The default training cap is 300,000 total entries. Subsampling operates on peptide groups to preserve target-decoy pairs and charge state groupings.

#### Best-Per-Precursor Subsampling (Multi-File Streaming)

For multi-file experiments (100+ files), naive peptide-group subsampling can be catastrophically bad: with 240 files, each peptide has ~480 entries (target + decoy × 240 files), so 300K entries yields only ~625 unique peptides — far too few for SVM training.

The streaming Percolator path uses **best-per-precursor subsampling**:

1. Find the best-scoring observation (by `coelution_sum`) per `base_id` across ALL files — one target and one decoy per precursor
2. This produces ~500K deduplicated entries with maximum peptide diversity
3. If this deduplicated set exceeds `max_train`, subsample peptide groups from it
4. Train SVM on the subsampled best-per-precursor entries (`train_only` mode suppresses FDR logging since q-values from the training subset are meaningless)
5. Score ALL entries across all files with the trained model by streaming through per-file Parquet caches

This ensures ~100K+ unique peptides in the training set regardless of file count, vs ~625 with naive subsampling.

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

Osprey writes one PIN file per mzML file with 21 feature columns. PIN files contain **only competition winners** (pre-competed using the best feature by ROC AUC).

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

## Feature Set (21 PIN Features)

All intensity-based spectral similarity scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks.

Osprey computes ~47 features per precursor in the `CoelutionFeatureSet` struct, but only **21 are written to the PIN file** for scoring. The remaining features were removed during feature weight optimization (they provided little discriminative value or were redundant). Osprey extracts fragment XICs and scores based on pairwise correlations and spectral matching.

| Category | Count | PIN Features |
|----------|-------|--------------|
| Pairwise coelution | 3 | fragment_coelution_sum, fragment_coelution_max, n_coeluting_fragments |
| Peak shape | 3 | peak_apex, peak_area, peak_sharpness |
| Spectral at apex | 3 | xcorr, consecutive_ions, explained_intensity |
| Mass accuracy | 2 | mass_accuracy_deviation_mean, abs_mass_accuracy_deviation_mean |
| RT deviation | 2 | rt_deviation, abs_rt_deviation |
| MS1 | 2 | ms1_precursor_coelution, ms1_isotope_cosine |
| Tukey median polish | 2 | median_polish_cosine, median_polish_residual_ratio |
| SG-weighted multi-scan | 4 | sg_weighted_xcorr, sg_weighted_cosine, median_polish_min_fragment_r2, median_polish_residual_correlation |

See [Pipeline Overview](README.md#feature-set-21-pin-features) for the full feature set documentation including removed features.

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

  1. Score each file: write full entries to Parquet, convert to FdrEntry stubs
  2. Load PIN features from Parquet, train model on all files
  3. Score all entries with trained model, update FdrEntry stubs
  4. Compute run-level q-values per file (precursor + peptide level)
  5. Select best observation per precursor across experiment
  6. Compute experiment-level q-values (precursor + peptide level)
  7. Apply max(precursor, peptide) at both run and experiment levels
  8. Determine passing precursors: experiment_qvalue <= threshold
  9. Reload full entries from Parquet only for files with passing precursors
  10. Include ALL per-file observations for passing precursors in output
  11. Propagate best experiment_qvalue to all observations
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

## FDR Filtering Level

Osprey computes q-values at precursor level (modified_sequence + charge), peptide level (modified_sequence), and protein level (when `--protein-fdr` is enabled), at both run and experiment scope. The `fdr_level` setting controls which level is used for output filtering.

### Configuration

CLI: `--fdr-level <precursor|peptide|protein|both>`

YAML:

```yaml
fdr_level: Peptide  # default
```

### Modes

| Mode | Filter applied | Description |
|------|---------------|-------------|
| `Precursor` | precursor q-value only | Filters on modified_sequence + charge. Least conservative. |
| `Peptide` (default) | peptide q-value only | Filters on modified_sequence. Typically the most biologically meaningful — a peptide with multiple charges gets only one competition opportunity at peptide level, making it stricter than precursor. |
| `Protein` | protein q-value only | Filters on protein group. Requires `--protein-fdr` to be enabled. |
| `Both` | max(precursor, peptide) | A precursor must pass at both levels. Most conservative dual FDR. |

The default is `Peptide` because peptide-level FDR produces biologically meaningful identifications — you want to know which peptides are present in the sample, not which (peptide, charge) combinations. `Both` is the most conservative and matches the behavior of Percolator/Mokapot when enforcing dual precursor + peptide FDR. Use `Precursor` if you want charge-state-specific filtering without the peptide-level constraint. Use `Protein` when you need all output to be gated by protein-level FDR (e.g., for downstream protein quantification workflows).

### How q-values are stored

Each `FdrEntry` stores six q-value fields:

- `run_precursor_qvalue`, `run_peptide_qvalue`, `run_protein_qvalue`
- `experiment_precursor_qvalue`, `experiment_peptide_qvalue`, `experiment_protein_qvalue`

The `effective_run_qvalue(level)` and `effective_experiment_qvalue(level)` methods return the appropriate value based on the configured `fdr_level`.

## Protein Parsimony and Protein-Level FDR

**Protein parsimony always runs** — bipartite graph construction, identical-set merging, subset elimination, and shared-peptide assignment are inexpensive deterministic steps that produce the authoritative peptide-to-protein-group mapping. Protein q-values (picked-protein FDR) are computed only when `--protein-fdr` is set.

See [16-protein-parsimony.md](16-protein-parsimony.md) for the full parsimony algorithm, including the iterative greedy set cover used for `Razor` mode.

### Configuration

CLI:

```bash
# Parsimony runs automatically; shared peptide mode defaults to All
osprey -i *.mzML -l library.tsv -o results.blib

# Enable protein FDR at 1% (parsimony still runs; picked-protein q-values computed)
osprey -i *.mzML -l library.tsv -o results.blib --protein-fdr 0.01

# With razor peptide assignment
osprey -i *.mzML -l library.tsv -o results.blib --protein-fdr 0.01 --shared-peptides razor

# Filter the blib output by protein-level FDR (requires --protein-fdr)
osprey -i *.mzML -l library.tsv -o results.blib --protein-fdr 0.01 --fdr-level protein
```

YAML:

```yaml
protein_fdr: 0.01       # Optional: enable picked-protein FDR at 1%
shared_peptides: All     # Options: All (default), Razor, Unique
fdr_level: Peptide       # Default; also accepts Precursor, Protein, Both
```

When `protein_fdr` is not set, parsimony still runs (producing the peptide-to-protein-group mapping) but no protein q-values are computed and the blib output is gated by the selected `fdr_level` (defaults to peptide-level).

### Two-Pass Picked-Protein FDR (Savitski 2015)

Osprey implements **true picked-protein FDR** (Savitski et al. 2015, PMC4563723). The algorithm runs twice in a two-pass architecture: first-pass is used for protein-aware gating during compaction and reconciliation, second-pass is authoritative for output.

**First-pass protein FDR** (runs before compaction, when `--protein-fdr` is set):

1. Build the parsimony graph from peptides passing first-pass peptide FDR.
2. Call `collect_best_peptide_scores()` on the **full pre-compaction** `per_file_entries`. This gives symmetric target + decoy peptide pools for picked-protein competition. No compaction survivorship bias.
3. Run `compute_protein_fdr()` with `qvalue_gate = config.run_fdr` (Savitski's 1% convention).
4. Propagate protein q-values into `FdrEntry.run_protein_qvalue` via `propagate_protein_qvalues()`.

**Protein-aware compaction** (uses first-pass protein q-values):

A peptide survives compaction if EITHER its peptide-level q-value is ≤ `reconciliation_compaction_fdr` (default 0.05, loosened from the old 0.01 run_fdr) OR its first-pass protein group q-value is ≤ `config.protein_fdr`. Rule (b) rescues borderline peptides whose protein has strong evidence. The loosened peptide gate alone already broadens the compaction pool; the protein rule is additive.

**Reconciliation consensus** (uses first-pass protein q-values):

`compute_consensus_rts()` includes a peptide in the consensus RT computation if it passes peptide-level FDR directly OR if its first-pass protein group passes protein-level FDR. Peptides from strong proteins can anchor the consensus even if their individual peptide q-values are borderline.

**Second-pass protein FDR** (runs after second-pass peptide FDR — AUTHORITATIVE):

1. Build the parsimony graph from peptides passing second-pass peptide FDR (at experiment level).
2. Call `collect_best_peptide_scores()` on the compacted + reconciled + second-pass-scored `per_file_entries`. These are the reconciliation-corrected scores.
3. Run `compute_protein_fdr()` with `qvalue_gate = config.run_fdr`.
4. Propagate into `FdrEntry.experiment_protein_qvalue` (leaves `run_protein_qvalue` alone from first-pass).
5. Write the protein CSV report.

### Picked-Protein Algorithm (Savitski 2015)

For each protein group:

1. **Score by best peptide**: `target_score = max(SVM score over target peptides passing gate)`, `decoy_score = max(SVM score over DECOY_-prefixed peptides passing gate)`. Uses the **single best peptide** per side, not a sum — Savitski explicitly rejected sum aggregation as length-biased.

2. **Pairwise picking**: each group produces exactly one winner. `target_score >= decoy_score` → target wins; otherwise → decoy wins. Groups with only a target side win as target; only a decoy side win as decoy; no peptides → skip.

3. **Cumulative FDR on winners**: sort by score descending (tiebreak: group_id ascending), `q = cum_decoys / max(1, cum_targets)` at each position, backward sweep for monotonicity.

4. **Target winners only** are exposed in `group_qvalues`. Decoy winners are statistical machinery for the FDR computation and are not written to output.

5. **Peptide propagation**: each peptide's q-value is the best (lowest) across the protein groups it belongs to. Important for shared peptides in `SharedPeptideMode::All`.

### Why Picked-Protein?

Classical protein-level TDC suffers from decoy over-representation: as genuine targets dominate the pool, random decoy matches accumulate disproportionately in the low-scoring region. **Pairwise picking eliminates this through structural symmetry** — each target-decoy pair produces exactly one winner, so the pool of winners is balanced by construction, and classical cumulative FDR on the winner list is well-calibrated.

Osprey's v26.1.2 used a DIA-NN-inspired composite log-likelihood approach which was found to produce artificially low protein q-values (most target groups passing at 1% FDR with very few decoy wins). The switch to true picked-protein corrects this calibration.

### Why SVM Score (Not q-Value or PEP) — and What Actually Happens in Practice

The ranking score used inside picked-protein is the **maximum peptide SVM discriminant** across the group, after gating targets on peptide q-value. We arrived at SVM by elimination, with an important empirical wrinkle we document honestly at the end.

- **Peptide q-value**: ranking by `min(run_peptide_qvalue)` looks like a literal reading of Savitski's text. It collapses the decoy null distribution — about 99% of decoy peptides have `q = 1.0` because they lost peptide-level TDC, so only the ~1% of decoys that won peptide TDC contribute meaningful scores. On Stellar 3-file HeLa this produced 2 decoy winners out of 6102 (~0.03%); cumulative FDR pinned at ~0.0003 and every target trivially passed 1% FDR.

- **Peptide PEP**: PEP is bounded `[0, 1]` but continuous, so in principle every decoy should contribute a meaningful score. In Osprey it collapses the null distribution too: [`PepEstimator`](../crates/osprey-ml/src/pep.rs) is fit on TDC winners only, its bins cover `[min_winner_score, max_winner_score]`, and `posterior_error()` clamps out-of-range queries to the nearest bin. Losing decoys have SVM scores below the winner range, so they all clamp to `bins[0] ≈ 1.0` regardless of their true density. Applying the estimator to all entries in the percolator code does not help because the clamping is intrinsic to the lookup. Fitting a second PepEstimator over the full entry range was considered but adds a separate PEP model just for protein FDR, and since PEP is a monotone function of SVM score the ranking would be identical to raw SVM anyway.

- **SVM discriminant** (current choice): every entry — target or decoy, TDC winner or loser — has a well-defined SVM score on the same scale. On paper this preserves the decoy null end-to-end.

- **Empirical outcome on Stellar**: prior to the RT-penalized peak selection fix, SVM-based picked-protein collapsed to 2 decoy winners out of 6099, producing pinned-at-zero cumulative FDR. The root cause was wrong-peak selection: interferer peaks at the wrong RT were selected because CWT candidates were ranked purely by coelution_score. Wrong-peak features inflated target SVM scores relative to decoys, eliminating the tail overlap picked-protein depends on. After adding the [Gaussian RT penalty to CWT peak selection](10-cross-run-reconciliation.md#rt-penalized-peak-selection), Stellar produces **57 decoy winners out of 6441 (0.89% decoy rate)** and 6384 target groups at 1% FDR — well calibrated. The RT penalty prevents wrong peaks from being selected, which restores natural target-decoy score overlap at the protein level.

The target-side **gate** still uses peptide q-value (`best_qvalue <= qvalue_gate`) so the reportable set matches Savitski's convention; only the ranking within that gate uses SVM.

### Sensitivity to Peak Selection Quality

The protein-level decoy-winner rate is a useful diagnostic of overall peak selection quality. When the RT penalty is effective and peaks are selected correctly, the target-decoy score distributions overlap naturally and picked-protein produces a meaningful null (~1% decoy winners at 1% FDR). If the decoy-winner rate drops to near zero on a dataset, it may indicate that wrong peaks are being selected (inflating target scores) or that the target-decoy asymmetry in the spectral library is too large.

**Decoy asymmetry remains a concern for future work.** Target library entries use AI-predicted spectra and retention times (Carafe), while decoy entries are produced by `DecoyGenerator::Reverse` with mechanical fragment reversal and the target's RT carried over unchanged. This asymmetry means the SVM can learn features that distinguish targets from decoys beyond what random matching would produce. Re-predicting decoy spectra and RTs from the same model used for targets would make the two pools more comparable and may further improve protein FDR calibration.

### Decoy Protein Pairing

Osprey's `DecoyGenerator` prefixes each decoy protein accession with `DECOY_` (e.g., target `P12345` becomes decoy `DECOY_P12345`). The picked-protein algorithm pairs target and decoy at the **peptide** level: for each peptide in the parsimony graph, it looks up the `DECOY_`-prefixed version in the peptide score map to compute the decoy-side score for each protein group the peptide belongs to.

### Protein Report Output

When `--protein-fdr` is set, Osprey writes a CSV protein report (`{output}.proteins.csv`) with columns:

| Column | Description |
|--------|-------------|
| `protein_group` | Numeric group ID |
| `protein_accessions` | Semicolon-separated protein accessions |
| `gene_names` | Semicolon-separated gene names (from library) |
| `protein_qvalue` | Protein group q-value from picked-protein FDR (second-pass, authoritative) |
| `best_peptide_score` | Maximum peptide SVM discriminant score among the group's peptides (the picked-protein ranking score) |
| `n_unique_peptides` | Number of unique (proteotypic) peptides |
| `n_shared_peptides` | Number of shared peptides |
| `unique_peptides` | Semicolon-separated unique peptide sequences |
| `shared_peptides` | Semicolon-separated shared peptide sequences |

Protein q-values are also propagated to each FdrEntry stub's `experiment_protein_qvalue` field for downstream filtering (used by `--fdr-level protein`).

> **Note**: Protein-level posterior error probability (PEP) is intentionally not computed. Peptide-level and precursor-level PEP (used internally by Percolator) are unaffected and remain available.

## References

- Percolator: Kall L, Canterbury J, Weston J, Noble WS, MacCoss MJ. Nat Methods. 2007;4(11):923-925.
- Score calibration: Granholm V, Noble WS, Kall L. BMC Bioinformatics. 2012;13 Suppl 16:S3.
- Mokapot: Fondrie WE, Noble WS. J Proteome Res. 2021;20(4):1966-1971.
- Target-decoy approach: Elias JE, Gygi SP. Nat Methods. 2007;4(3):207-214.
- Conservative FDR: Levitsky LI, et al. J Proteome Res. 2017;16(2):689-694.
- Picked protein FDR: Savitski MM, et al. Mol Cell Proteomics. 2015;14(9):2394-2404.
- Picked protein group FDR: The M, Kall L. J Proteome Res. 2016;15(4):1456-1461.
