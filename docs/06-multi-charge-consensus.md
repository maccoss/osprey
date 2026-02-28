# Multi-Charge-State Peak Consensus

## Problem

A peptide elutes at one retention time regardless of its ionization charge state. In DIA, different charge states of the same peptide (e.g., PEPTIDEK+2 and PEPTIDEK+3) have different precursor m/z values and therefore fall into **different DIA isolation windows**. Because isolation windows are processed by independent parallel tasks, each charge state finds its peak independently — they may end up selecting different chromatographic peaks.

This is physically inconsistent: if a peptide is detected at charge 2+ at RT=25.3 and at charge 3+ at RT=27.1, one of them is wrong. Forcing all charge states to share the same peak RT and integration boundaries improves both accuracy and downstream quantification consistency.

## Solution: Post-Search Consensus

Cross-charge consensus **cannot** happen within the per-window parallel loop because different charge states are in different windows. Instead, Osprey adds a post-processing step after the main search, before FDR control.

### Execution Order in `run_search()`

```
1. Parallel per-window search (existing) → all_entries
2. Overlapping-window dedup → deduped
3. Multi-charge consensus selection → identifies entries needing re-scoring
4. Re-scoring at consensus boundaries → replaces entries with corrected features
5. Return final results
   ↓
6. First-pass FDR control (caller) → run-level q-values
   ↓
7. Cross-run reconciliation (caller, multi-file only)
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

1. The charge state with the **highest `coelution_sum`** defines the consensus peak
   - `coelution_sum` = sum of all pairwise fragment correlations within the peak
   - Higher coelution_sum means stronger, more reliable fragment co-elution evidence
   - This charge state's peak boundaries (apex_rt, start_rt, end_rt) become the consensus

2. For each other charge state in the group:
   - If its `apex_rt` is within half the consensus peak width (minimum 0.1 min) of the consensus apex → **already agrees**, keep as-is
   - Otherwise → **mark for re-scoring** at the consensus boundaries

### Step 3: Re-Score at Consensus Boundaries

For entries that selected a different peak than the consensus:

1. Look up the library entry and find the correct DIA isolation window for this precursor's m/z
2. Gather spectra from that window near the consensus apex RT
3. Extract fragment XICs from those spectra
4. Map the consensus RT boundaries to XIC scan indices (find closest RT in XIC time grid)
5. Compute all 45 features at the consensus boundaries via `compute_features_at_peak()`
6. If no signal is found at the consensus RT (too few spectra, no fragment evidence): **drop the entry** — this charge state has no evidence at the peptide's true elution time

### Step 4: Merge Results

The final result combines:
- Entries that were kept as-is (single charge state, or already agreed with consensus)
- Re-scored entries with features computed at consensus boundaries
- Dropped entries are excluded (no evidence at consensus RT)

Results are re-sorted by `(entry_id, scan_number)` for deterministic output.

---

## Why "Best Coelution Sum Wins"

The charge state with the highest pairwise fragment correlation is most likely to have found the true elution peak because:
- **Fragment co-elution** is the strongest discriminator between true peaks and interference
- A charge state that happens to pick an interference peak will have low co-elution (fragments don't co-elute at a false peak)
- The best-co-eluting charge state has the most reliable boundaries for defining where the peptide actually elutes

Alternative strategies considered:
- **Highest apex intensity**: Biased by ionization efficiency and window interference
- **Closest to expected RT**: Penalizes real peaks displaced by calibration error
- **Majority vote**: Complex, and most peptides have only 2-3 charge states

---

## Performance

- **Consensus selection**: O(n) HashMap grouping — negligible
- **Re-scoring**: Only for entries that selected a different peak than consensus. Expected ~2-6% of entries (most multi-charge peptides already agree on their peak). Each re-scoring has the same cost as initial scoring.
- **Net overhead**: ~2-6% additional scoring time in `run_search`

---

## Example

```
Peptide: PEPTIDEK
  Charge 2+: precursor m/z = 458.25 → DIA window [455, 465]
  Charge 3+: precursor m/z = 305.83 → DIA window [300, 310]

After independent parallel search:
  Charge 2+: apex_rt = 25.3, coelution_sum = 8.5, boundaries = [24.8, 25.9]
  Charge 3+: apex_rt = 27.1, coelution_sum = 3.2, boundaries = [26.6, 27.5]

Consensus selection:
  Best coelution_sum: charge 2+ (8.5 > 3.2)
  Consensus: apex = 25.3, boundaries = [24.8, 25.9]
  Peak width = 1.1, tolerance = max(0.55, 0.1) = 0.55 min

  Charge 2+: apex_diff = 0.0 ≤ 0.55 → keep as-is
  Charge 3+: apex_diff = 1.8 > 0.55 → needs re-scoring

Re-scoring charge 3+ at consensus boundaries [24.8, 25.9]:
  Find spectra from window [300, 310] near RT = 25.3
  Extract fragment XICs
  Map [24.8, 25.9] to XIC scan indices
  Compute all 45 features at these boundaries

Result:
  Both charge states now share RT = 25.3 and boundaries [24.8, 25.9]
  Each is scored independently at those shared boundaries
```

---

## Scope

- **Main search only**: Multi-charge consensus is applied in `pipeline.rs::run_search()`, not during calibration
- **Per-file**: Consensus is computed within each file independently (different files may have different chromatographic conditions)
- **Targets and decoys separately**: Grouped by `modified_sequence` which naturally separates targets from decoys (decoys have `DECOY_` prefix)

### What stays unchanged

- CWT peak detection within the per-precursor loop
- All 45 PIN features and Mokapot/Percolator scoring
- FDR control, blib output, everything downstream
- Calibration path (`batch.rs`) — does not need multi-charge consensus
- Existing blib `build_shared_boundaries()` function (uses same grouping pattern for output)

---

## Determinism

See [Determinism](09-determinism.md) for comprehensive determinism documentation.

The multi-charge consensus step follows the project's determinism patterns:
- After HashMap-based grouping, results are sorted by `(entry_id, scan_number)` before return
- Re-scoring is done sequentially (not `par_iter`) to avoid nondeterministic floating-point accumulation order

---

## Implementation

| File | Function | Purpose |
|------|----------|---------|
| `crates/osprey/src/pipeline.rs` | `select_consensus_peaks()` | Group by modified_sequence, identify consensus, mark re-score targets |
| `crates/osprey/src/pipeline.rs` | `rescore_at_consensus()` | Re-extract XICs and compute features at consensus boundaries |
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
