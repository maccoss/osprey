# Multi-Charge-State Peak Consensus

## Problem

A peptide elutes at one retention time regardless of its ionization charge state. In DIA, different charge states of the same peptide (e.g., PEPTIDEK+2 and PEPTIDEK+3) have different precursor m/z values and therefore fall into **different DIA isolation windows**. Because isolation windows are processed by independent parallel tasks, each charge state finds its peak independently — they may end up selecting different chromatographic peaks.

This is physically inconsistent: if a peptide is detected at charge 2+ at RT=25.3 and at charge 3+ at RT=27.1, one of them is wrong. Forcing all charge states to share the same peak RT and integration boundaries improves both accuracy and downstream quantification consistency.

## Solution: Post-FDR Consensus

Cross-charge consensus **cannot** happen within the per-window parallel loop because different charge states are in different windows. Osprey adds a **post-FDR** step that identifies charge states at different peaks and re-scores them at the consensus position using the same `run_search()` function with boundary overrides.

### Execution Order

```
1. Parallel per-window search (run_search) → all_entries
2. Overlapping-window dedup → deduped
3. First-pass FDR → run-level q-values
   ↓
4. select_post_fdr_consensus() → identifies entries needing re-scoring
5. Re-scoring via run_search() with boundary_overrides
   (merged with reconciliation targets if multi-file)
   ↓
6. Second-pass FDR → final q-values
   ↓
7. Cross-run reconciliation (multi-file only)
   → consensus RT computation, calibration refit, second-pass FDR
   See: Cross-Run Reconciliation (10-cross-run-reconciliation.md)
```

---

## Algorithm

### Step 1: Group by Modified Sequence

Entries are grouped by `modified_sequence`. This naturally separates:
- **Different charge states** of the same peptide (same modified_sequence, different charge)
- **Targets from decoys** (decoys have `DECOY_` prefix in modified_sequence)
- **Different peptides** (different modified_sequence)

### Step 2: Select Consensus Peak

For each group with multiple charge states:

1. The charge state with the **highest SVM score** among those passing the FDR threshold defines the consensus peak
   - SVM score is preferred over `coelution_sum` because it incorporates all 21 features and reflects the trained model's assessment
   - Only FDR-passing charge states are considered as consensus leaders
   - This charge state's peak boundaries (apex_rt, start_rt, end_rt) become the consensus

2. For each other charge state in the group:
   - If its `apex_rt` is within half the consensus peak width (minimum 0.1 min) of the consensus apex → **already agrees**, keep as-is
   - Otherwise → **mark for re-scoring** at the consensus boundaries

### Step 3: Re-Score via `run_search()` with Boundary Overrides

Entries marked for re-scoring are collected as `boundary_overrides` — a `HashMap<u32, (f64, f64, f64)>` mapping `entry_id` to `(apex_rt, start_rt, end_rt)`. These are passed to `run_search()` along with a subset library containing only the entries to re-score.

Inside `run_search()`, entries with boundary overrides:
1. **Skip the pre-filter** (we're told to score here regardless of initial signal)
2. **Skip CWT peak detection** (boundaries are already determined from the consensus leader)
3. **Map consensus RT boundaries to XIC scan indices** via binary search
4. **Compute all features at the consensus boundaries** via `compute_features_at_peak()`
5. If no signal is found at the consensus RT (too few spectra, no fragment evidence): **drop the entry**

This reuses the same parallel window processing and per-window XCorr preprocessing as the first-pass search, so consensus re-scoring has the same performance characteristics. See [Boundary Overrides & Forced Integration](13-boundary-overrides.md) for details on the override mechanism.

### Step 4: Merge Results

The final result combines:
- Entries that were kept as-is (single charge state, or already agreed with consensus)
- Re-scored entries with features computed at consensus boundaries
- Dropped entries are excluded (no evidence at consensus RT)

Results are re-sorted by `(entry_id, scan_number)` for deterministic output.

---

## Why "Best SVM Score Wins"

The charge state with the highest SVM score among FDR-passing entries is most likely to have found the true elution peak because:
- **SVM score** integrates all 21 discriminative features (coelution, spectral, mass accuracy, RT deviation, etc.) into a single score optimized by the semi-supervised training
- A charge state at an interference peak will have a lower SVM score than one at the true peak
- Using SVM score (available after the first-pass FDR) gives better consensus selection than raw coelution_sum

Alternative strategies considered:
- **Highest coelution_sum**: Good discriminator but doesn't use the full feature set; used in the pre-FDR design
- **Highest apex intensity**: Biased by ionization efficiency and window interference
- **Closest to expected RT**: Penalizes real peaks displaced by calibration error
- **Majority vote**: Complex, and most peptides have only 2-3 charge states

---

## Performance

- **Consensus selection**: O(n) HashMap grouping — negligible
- **Re-scoring**: Only for entries that selected a different peak than consensus. Expected ~2-6% of entries. Re-scoring uses `run_search()` with boundary overrides, so it benefits from parallel window processing and shared XCorr preprocessing.
- **Combined with reconciliation**: In multi-file experiments, consensus and reconciliation targets are merged into a single `run_search()` call per file, avoiding redundant spectra loading.

---

## Example

```
Peptide: PEPTIDEK
  Charge 2+: precursor m/z = 458.25 → DIA window [455, 465]
  Charge 3+: precursor m/z = 305.83 → DIA window [300, 310]

After first-pass search and FDR:
  Charge 2+: apex_rt = 25.3, svm_score = 2.1, boundaries = [24.8, 25.9], qvalue = 0.002
  Charge 3+: apex_rt = 27.1, svm_score = 0.8, boundaries = [26.6, 27.5], qvalue = 0.03

Post-FDR consensus selection:
  Best SVM score among FDR-passing: charge 2+ (2.1 > 0.8)
  Consensus: apex = 25.3, boundaries = [24.8, 25.9]
  Peak width = 1.1, tolerance = max(0.55, 0.1) = 0.55 min

  Charge 2+: apex_diff = 0.0 ≤ 0.55 → keep as-is
  Charge 3+: apex_diff = 1.8 > 0.55 → needs re-scoring

Re-scoring charge 3+ via run_search() with boundary_overrides:
  boundary_overrides = {entry_id_3plus: (25.3, 24.8, 25.9)}
  run_search() processes window [300, 310] with override:
    Skip pre-filter, skip CWT
    Map [24.8, 25.9] to XIC scan indices
    Compute all 21 features at these boundaries

Result:
  Both charge states now share RT = 25.3 and boundaries [24.8, 25.9]
  Each is scored independently at those shared boundaries
```

---

## Scope

- **Post-FDR**: Consensus selection happens after the first-pass FDR so that SVM scores are available and only FDR-passing charge states can be consensus leaders
- **Per-file**: Consensus is computed within each file independently (different files may have different chromatographic conditions)
- **Targets and decoys separately**: Grouped by `modified_sequence` which naturally separates targets from decoys (decoys have `DECOY_` prefix)
- **Merged with reconciliation**: In multi-file experiments, consensus targets and reconciliation targets are merged into a single `run_search()` call per file. If both consensus and reconciliation want to re-score the same entry, reconciliation wins (it has refined calibration + consensus RT).

### What stays unchanged

- CWT peak detection in the first-pass per-precursor loop
- All 21 PIN features and Percolator/Mokapot scoring
- FDR control, blib output, everything downstream
- Calibration path (`batch.rs`) — does not need multi-charge consensus

---

## Determinism

See [Determinism](09-determinism.md) for comprehensive determinism documentation.

The multi-charge consensus step follows the project's determinism patterns:
- After HashMap-based grouping, results are sorted by `(entry_id, scan_number)` before return
- Re-scoring goes through `run_search()`, which has its own determinism guarantees (sorted output, deterministic window processing)

---

## Implementation

| File | Function | Purpose |
|------|----------|---------|
| `crates/osprey/src/pipeline.rs` | `select_post_fdr_consensus()` | Group by modified_sequence, find best FDR-passing charge state, mark re-score targets |
| `crates/osprey/src/pipeline.rs` | `run_search()` | Re-scores entries at consensus boundaries via `boundary_overrides` |
| `crates/osprey/src/pipeline.rs` | `compute_features_at_peak()` | Reusable feature computation (shared with initial search) |
| `crates/osprey/src/pipeline.rs` | `FeatureComputeContext` | Context struct bundling read-only references for feature computation |

### Tests

| Test | Description |
|------|-------------|
| `test_consensus_single_charge_no_change` | Single charge state: kept, no re-scoring |
| `test_consensus_two_charges_same_peak` | Two charges, same apex RT: both kept |
| `test_consensus_two_charges_different_peaks` | Two charges, different apex RTs: best kept, other marked |
| `test_consensus_three_charges_two_agree` | Three charges, two agree + one disagrees: two kept, one marked |
| `test_consensus_decoys_separate_from_targets` | Decoys grouped separately from targets |
| `test_consensus_multiple_peptides_independent` | Different peptides grouped independently |
