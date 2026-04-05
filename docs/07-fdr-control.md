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
  8. Protein FDR (optional): picked-protein competition using peptide PEP
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

Osprey computes q-values at both precursor level (modified_sequence + charge) and peptide level (modified_sequence only), at both run and experiment scope. The `fdr_level` setting controls which level is used for output filtering.

### Configuration

CLI: `--fdr-level <precursor|peptide|both>`

YAML:

```yaml
fdr_level: Both  # default
```

### Modes

| Mode | Filter applied | Description |
|------|---------------|-------------|
| `Precursor` | precursor q-value only | Filters on modified_sequence + charge. Less conservative. |
| `Peptide` | peptide q-value only | Filters on modified_sequence (ignoring charge). |
| `Both` (default) | max(precursor, peptide) | A precursor must pass at both levels. Most conservative. |

The `Both` mode is the most conservative and matches the behavior of Percolator and Mokapot, which enforce dual precursor + peptide level FDR. Use `Precursor` if you want charge-state-specific filtering without the peptide-level constraint.

### How q-values are stored

Each `FdrEntry` stores six q-value fields:

- `run_precursor_qvalue`, `run_peptide_qvalue`, `run_protein_qvalue`
- `experiment_precursor_qvalue`, `experiment_peptide_qvalue`, `experiment_protein_qvalue`

The `effective_run_qvalue(level)` and `effective_experiment_qvalue(level)` methods return the appropriate value based on the configured `fdr_level`.

## Protein-Level FDR

Osprey supports optional protein-level FDR control using native Rust protein parsimony and the picked-protein approach.

### Configuration

CLI:

```bash
# Enable protein FDR at 1%
osprey -i *.mzML -l library.tsv -o results.blib --protein-fdr 0.01

# With razor peptide assignment
osprey -i *.mzML -l library.tsv -o results.blib --protein-fdr 0.01 --shared-peptides razor
```

YAML:

```yaml
protein_fdr: 0.01       # Enable protein-level FDR at 1%
shared_peptides: All     # Options: All (default), Razor, Unique
```

When `protein_fdr` is not set (the default), protein-level FDR is disabled and all peptides passing precursor/peptide FDR are included in the output.

### Protein Parsimony Algorithm

Protein inference runs once on the spectral library at startup. Protein accessions come from `LibraryEntry.protein_ids` (loaded from DIA-NN TSV, elib, or blib libraries). No FASTA file is required.

```
Protein Parsimony Steps:
  1. Build bipartite graph: peptide <-> protein from target library entries
  2. Group proteins with identical peptide sets into protein groups
  3. Subset elimination: remove groups whose peptides are a strict subset of another
  4. Classify peptides as unique (1 group) or shared (2+ groups)
  5. Apply shared peptide mode (All, Razor, or Unique)
```

### Shared Peptide Modes

| Mode | Behavior | Use case |
|------|----------|----------|
| `All` (default) | Shared peptides contribute to all their protein groups. Each peptide inherits the best (lowest) protein q-value among its groups. | Maximum sensitivity. |
| `Razor` | Shared peptides assigned exclusively to the protein group with the most unique peptides (tiebreak: lowest group ID). After assignment, treated as unique. | Balanced approach. Matches MaxQuant's razor peptide logic. |
| `Unique` | Shared peptides excluded from protein scoring and output entirely. Only unique peptides are used. | Most conservative. |

### Picked-Protein FDR

After SVM scoring, protein-level FDR is computed using the picked-protein approach (Savitski et al., 2015; The and Kall, 2016):

```
Picked-Protein FDR Steps:
  1. Collect best peptide PEP per modified_sequence BEFORE stub compaction
     (must include decoy competition winners that are dropped during compaction)
  2. Score each target protein group: best (lowest) target peptide PEP
  3. Score each decoy "shadow" group: best (lowest) decoy peptide PEP
     (decoy peptides are DECOY_ prefixed versions of the target group's peptides)
  4. Picked competition: for each group, target vs decoy shadow -- lower PEP wins
  5. Compute q-values via standard TDC on picked proteins
  6. Estimate protein PEP via kernel density estimation on picked winners
  7. Propagate protein q-values to peptides (best q-value among groups)
```

**Why PEP instead of SVM score**: The SVM discriminant score is trained to maximize separation between targets and decoys. At the protein level, the best target peptide's SVM score almost always exceeds the best decoy peptide's score, producing nearly zero decoy wins and degenerate FDR. PEP (posterior error probability) is calibrated as a probability with meaningful overlap between target and decoy distributions, enabling proper picked-protein competition.

**Pre-compaction collection**: After first-pass FDR, non-passing stubs are compacted to reduce memory. This drops decoy entries that won their peptide-level competition (their base_ids aren't in `first_pass_base_ids`). Peptide scores must be collected before compaction to include these decoy winners.

Protein FDR is computed at experiment level (across all files) and the protein q-values are propagated to the `experiment_protein_qvalue` field on each FdrEntry stub.

### Decoy Protein Pairing

Osprey's `DecoyGenerator` prefixes each decoy protein accession with `DECOY_` (e.g., target `P12345` becomes decoy `DECOY_P12345`). The picked-protein approach pairs target and decoy proteins by stripping this prefix, so each target protein group competes against its decoy counterpart.

### Protein Report Output

When `--protein-fdr` is set, Osprey writes a CSV protein report (`{output}.proteins.csv`) with columns:

| Column | Description |
|--------|-------------|
| `protein_group` | Numeric group ID |
| `protein_accessions` | Semicolon-separated protein accessions |
| `gene_names` | Semicolon-separated gene names (from library) |
| `protein_qvalue` | Protein group q-value from picked-protein TDC |
| `protein_pep` | Protein group posterior error probability |
| `best_score` | Best peptide score (negative PEP) for the group |
| `n_unique_peptides` | Number of unique (proteotypic) peptides |
| `n_shared_peptides` | Number of shared peptides |
| `unique_peptides` | Semicolon-separated unique peptide sequences |
| `shared_peptides` | Semicolon-separated shared peptide sequences |

Protein q-values are also propagated to each FdrEntry stub's `experiment_protein_qvalue` field for downstream filtering.

## References

- Percolator: Kall L, Canterbury J, Weston J, Noble WS, MacCoss MJ. Nat Methods. 2007;4(11):923-925.
- Score calibration: Granholm V, Noble WS, Kall L. BMC Bioinformatics. 2012;13 Suppl 16:S3.
- Mokapot: Fondrie WE, Noble WS. J Proteome Res. 2021;20(4):1966-1971.
- Target-decoy approach: Elias JE, Gygi SP. Nat Methods. 2007;4(3):207-214.
- Conservative FDR: Levitsky LI, et al. J Proteome Res. 2017;16(2):689-694.
- Picked protein FDR: Savitski MM, et al. Mol Cell Proteomics. 2015;14(9):2394-2404.
- Picked protein group FDR: The M, Kall L. J Proteome Res. 2016;15(4):1456-1461.
