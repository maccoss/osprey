# Cross-Run Peak Reconciliation

Cross-run reconciliation aligns peak integration boundaries across replicate files so that the same peptide is quantified at a consistent chromatographic position in every run. It runs after the initial per-run FDR pass, before the final experiment-level FDR, and only applies to multi-file experiments.

## Problem

Without reconciliation, each file's per-window search finds peaks independently using CWT peak detection and pairwise coelution scoring. For a given peptide, the search may confidently find the correct peak in most runs but select a completely different peak — a co-eluting interferer, a noise spike, or an isomer — in one or more runs where the signal is weaker or the interference is stronger.

```
Peptide PEPTIDEK:
  file1.mzML: coelution_sum=8.5  → apex=25.3 (correct peak)
  file2.mzML: coelution_sum=7.9  → apex=25.4 (correct peak)
  file3.mzML: coelution_sum=3.1  → apex=31.7 (wrong peak — low-signal run picked an interferer)
```

In this case, file3 has detected a completely different peptide or artifact. Cross-run reconciliation uses the high-confidence runs (file1, file2) to establish where PEPTIDEK actually elutes, then goes back to file3 and either finds an alternate CWT candidate at that RT or imputes the integration boundaries at the expected position.

---

## Algorithm Overview

```
After per-run FDR (first pass):
  1. collect_consensus_peptides() — targets at consensus_fdr experiment-level (from FdrEntry stubs)
  2. compute_consensus_rts()      — weighted median library RT across runs (from FdrEntry stubs)
  3. refit_calibration()          — tighter per-run LOESS from consensus points
  4. plan_reconciliation()        — Keep | UseCwtPeak | ForcedIntegration per entry
       (CWT candidates loaded selectively from per-file Parquet caches)
  5. Re-score via run_search()    — boundary_overrides skip CWT; parallel window processing
  6. Second-pass FDR              — final experiment-level q-values
```

---

## Step 1: Consensus Peptide Selection

Target peptides qualify for consensus RT computation if they pass `consensus_fdr` (default 1%) at **run-level peptide FDR** in at least one replicate. Their paired decoys (DECOY_ prefix) are included so target-decoy competition remains fair in the second FDR pass.

```
targets: {PEPTIDEK, ANOTHERR, ...}     ← pass peptide-level FDR
decoys:  {DECOY_PEPTIDEK, DECOY_ANOTHERR, ...}  ← paired decoys included
```

Both sides get independent consensus RTs — decoy consensus RTs are computed from decoy detections (not from target detections) to avoid any information leakage.

### Protein FDR Rescue Gate

When `--protein-fdr` is set, `compute_consensus_rts()` accepts an additional `protein_fdr_threshold` parameter. A target peptide qualifies for consensus RT computation if EITHER:

1. Its peptide-level q-value is ≤ `consensus_fdr`, OR
2. Its first-pass protein group's q-value is ≤ `protein_fdr_threshold`

Rule (2) lets peptides from strong proteins contribute to consensus RT computation even when their individual peptide q-values are borderline. This is particularly valuable for peptides that are weakly scored on their own but belong to a protein with many other strong peptides — the protein evidence "rescues" them into the consensus anchor set.

First-pass protein FDR runs before compaction and reconciliation (see [07-fdr-control.md](07-fdr-control.md#two-pass-picked-protein-fdr-savitski-2015)), so `run_protein_qvalue` is already populated by the time `compute_consensus_rts()` is called. When `--protein-fdr` is not set, rule (2) is disabled (`protein_fdr_threshold = 0.0`).

---

## Step 2: Consensus Library RT

For each consensus peptide, its detections across all files are collected:

```
(file_name, apex_rt, coelution_sum, peak_width)
```

Each measured `apex_rt` is mapped back to **library RT space** using the run's inverse RT calibration:

```
library_rt = RTCalibration::inverse_predict(apex_rt)
```

The **consensus library RT** is the **weighted median** of these library RT values, with weights equal to `coelution_sum` (pairwise fragment correlation sum). Higher coelution_sum means more reliable detection → stronger weight in the consensus.

```
PeptideConsensusRT {
    modified_sequence: String,
    is_decoy: bool,
    consensus_library_rt: f64,   // weighted median in library RT space
    median_peak_width: f64,       // median peak width across detections (minutes)
    n_runs_detected: usize,
}
```

Working in library RT space is important: it decouples the consensus from run-specific RT drift. The consensus library RT can then be mapped into any run's measured RT space using that run's calibration.

---

## Step 3: Calibration Refit

For each run, the LOESS RT calibration is refit using `(consensus_library_rt → measured_apex_rt)` pairs from **target** peptides detected in that run at experiment-level `consensus_fdr`.

This produces a tighter calibration than the initial discovery pass because:
- The calibration points are high-confidence FDR-controlled detections, not noisy initial matches
- Outlier removal is disabled (`outlier_retention=1.0`) — the LOESS robustness iterations still downweight genuine outliers
- The consensus library RTs are more consistent across the gradient than first-pass matches

The refined calibration is stored per-run and used in Step 4. If too few consensus points exist for a given run (< 20), that run falls back to the original calibration.

---

## Step 4: Reconciliation Planning

For each scored entry (target or decoy), the expected measured RT is predicted from the refined calibration:

```
expected_rt = refined_calibration.predict(consensus_library_rt)
```

Three actions are possible based on **apex proximity** to the expected RT:

| Action | Condition | What happens |
|--------|-----------|--------------|
| **Keep** | `|apex_rt - expected_rt| <= rt_tolerance` | No change |
| **UseCwtPeak** | A CWT candidate has `|cand.apex_rt - expected_rt| <= rt_tolerance` | Switch to the closest-apex candidate |
| **ForcedIntegration** | Neither current peak nor any CWT candidate has apex within tolerance | Integrate at `expected_rt ± median_peak_width/2` |

### RT Tolerance: Global MAD, Not Local Interpolation

The `rt_tolerance` is computed **once per file** as:

```
rt_tolerance = max(0.1, 3.0 × MAD × 1.4826)
```

where `MAD` is the median absolute residual from the refined calibration's training points (thousands of consensus peptides across the gradient). This is a global measure of how tightly the calibration fits, and it is robust to individual outliers because it uses the median, not per-point residuals.

**Why not per-point local tolerance?** An earlier version used `local_tolerance(query_rt)` which interpolated the absolute residual from the nearest training points at the query RT. This created a self-fulfilling prophecy: a peptide with a wrong apex RT contributed a huge residual to the calibration training set, which inflated the local tolerance at that RT, which then allowed the wrong-RT detection to pass the apex proximity check. Using a global MAD eliminates this feedback loop — one bad peptide's residual barely moves the median.

A minimum floor of 0.1 min prevents being too strict when the calibration is exceptionally tight (e.g., very short gradients with highly reproducible LC).

### Apex Proximity vs Boundary Containment

Earlier versions used **boundary containment** (`start_rt <= expected_rt <= end_rt`) to decide Keep vs re-score. That check was too permissive: a peak detected at the wrong RT (e.g., apex 1.2 min off) but with a wide trailing edge that happened to span the expected RT was classified as Keep. The apex proximity check rejects such peaks regardless of their boundary width.

For CWT candidate selection, apex proximity also picks the candidate whose apex is **closest** to the expected RT, not just any candidate whose boundaries overlap. If a wrong peak's boundaries span the expected RT but a better candidate's apex is closer, the better candidate wins.

The CWT candidates are stored during the initial search (configurable, default 5 per precursor) and persisted in the per-file Parquet score cache. During reconciliation planning, only the CWT candidate column is loaded selectively via `load_cwt_candidates_from_parquet()`, avoiding the cost of reloading full entries with all features and fragment data. This lets reconciliation switch to an alternate peak that was already detected but not selected as the best, without requiring a full re-extraction.

`ForcedIntegration` is a last resort — it integrates at the expected position even if no peak (current or CWT candidate) has its apex within tolerance, using the median peak width from all detections as the integration window.

---

## Step 5: Re-Scoring via `run_search()` with Boundary Overrides

Entries with `UseCwtPeak` or `ForcedIntegration` actions are re-scored by calling the same `run_search()` function used in the first-pass search, with `boundary_overrides` that direct it to score at specific RT boundaries instead of running CWT peak detection.

### How It Works

1. **Merge targets**: Multi-charge consensus targets and reconciliation targets are merged into a single set per file. If a reconciliation action conflicts with a consensus action for the same entry, reconciliation wins (it has more information: refined calibration + consensus RT).

2. **Build boundary overrides**: Each `ReconcileAction` is converted to a `(apex_rt, start_rt, end_rt)` tuple and stored in a `HashMap<u32, (f64, f64, f64)>` keyed by `entry_id`:
   - `UseCwtPeak` → uses the stored CWT candidate's `(apex_rt, start_rt, end_rt)`
   - `ForcedIntegration` → `(expected_rt, expected_rt - half_width, expected_rt + half_width)`

3. **Subset library**: Only library entries corresponding to entries needing re-scoring are included (subset by `entry_id`).

4. **Call `run_search()`**: The subset library, loaded spectra, calibration, and `Some(&boundary_overrides)` are passed to `run_search()`. Inside `run_search()`, entries with overrides:
   - **Skip the pre-filter** (we're told to score here regardless of signal)
   - **Skip CWT peak detection** (boundaries are already determined)
   - **Map target RTs to XIC scan indices** via binary search (`partition_point`)
   - **Compute all 21 PIN features** at the override boundaries via `compute_features_at_peak()`

5. **Update FdrEntry stubs**: Re-scored entries replace the originals. Entries that cannot be scored (no spectral data at the expected RT) are dropped.

### Why `run_search()` Instead of Sequential Re-Scoring

Reconciliation previously used a dedicated sequential re-scoring loop that processed entries one at a time. This was ~50x slower than the first-pass search because it missed the parallelism of `run_search()`:

- **First-pass search**: Isolation windows processed in parallel via `rayon`, spectra grouped by window, XCorr preprocessing done once per window and reused across all candidates
- **Old reconciliation**: Each entry loaded spectra, extracted XICs, and computed features independently — no window-level parallelism, redundant XCorr preprocessing

By reusing `run_search()` with `boundary_overrides`, reconciliation gets the same performance characteristics as the first pass:
- **Parallel window processing** via rayon
- **Per-window XCorr preprocessing** done once and shared across all entries in the window (see [XCorr Scoring](04-xcorr-scoring.md#per-window-preprocessing-optimization))
- **Shared spectra grouping and indexing** across entries

### Memory: Sequential Per-File Processing

Files are still processed **sequentially** (not in parallel) because each file requires loading spectra (~1–2 GB) and building the window index. Parallel file loading would OOM on large experiments. The sequencing is at the file level; within each file, `run_search()` parallelizes across isolation windows.

---

## Step 6: Second-Pass FDR

After reconciliation, FDR is recomputed on the updated entry set:
- Features have been recomputed at consensus-aligned boundaries
- The same Percolator/Mokapot/Simple FDR pipeline applies
- This produces the **final experiment-level q-values** written to the blib output

---

## Configuration

```yaml
reconciliation:
  enabled: true         # Set false to skip reconciliation entirely
  consensus_fdr: 0.01   # FDR threshold for selecting consensus peptides
```

CLI flag: `--no-reconciliation` disables it entirely.

The number of CWT candidates stored per precursor is configured separately:

```yaml
# (set in search config — affects memory for all search results)
n_cwt_candidates: 5
```

---

## When Reconciliation Is Skipped

- **Single-file experiments**: No cross-run consensus possible; experiment-level FDR equals run-level FDR directly.
- **Disabled via config**: `reconciliation.enabled: false` or `--no-reconciliation`.
- **Too few consensus points per run**: If fewer than 20 consensus points exist for a run, that run uses its original calibration but still participates in reconciliation planning.

---

## Example

```
Peptide PEPTIDEK detected in 3 files:
  file1.mzML: apex_rt=25.3, coelution_sum=8.5, peak_width=0.9, library_rt_mapped=32.1
  file2.mzML: apex_rt=25.4, coelution_sum=7.9, peak_width=1.0, library_rt_mapped=32.0
  file3.mzML: apex_rt=31.7, coelution_sum=3.1, peak_width=1.2, library_rt_mapped=38.5
              ^^^ wrong peak: low coelution score, maps to wrong library RT

Weighted median library RT (weight = coelution_sum):
  weights = [8.5, 7.9, 3.1], values = [32.1, 32.0, 38.5]
  consensus_library_rt ≈ 32.1  (file3's outlier point doesn't move the weighted median)

Per-run refined calibration → expected measured RT (tolerance ~0.4 min from global MAD):
  file1: predict(32.1) = 25.3 → |25.3 - 25.3| = 0.00 ≤ 0.4 → Keep
  file2: predict(32.1) = 25.4 → |25.4 - 25.4| = 0.00 ≤ 0.4 → Keep
  file3: predict(32.1) = 25.2 → |31.7 - 25.2| = 6.5 > 0.4 → NOT Keep
           → check stored CWT candidates by apex proximity:
             candidate 0 apex=25.2 → |25.2 - 25.2| = 0.00 ≤ 0.4 → closest
             candidate 1 apex=31.7 → |31.7 - 25.2| = 6.5 > 0.4 → rejected
           → UseCwtPeak (candidate 0, apex=25.2)
             (CWT had found the correct peak but scored it lower than the interferer)

Result:
  file1, file2: unchanged — kept their correct peaks
  file3: switched to the stored CWT candidate at RT≈25.2 (the true elution time)
         features recomputed at new boundaries
```

If file3 had no CWT candidate with its apex within tolerance (e.g., the peptide was truly not present or below detection), the outcome would be `ForcedIntegration` — boundaries placed at `expected_rt ± median_peak_width/2`. This provides a quantification window for Skyline even when the peptide was not detected in that run.

---

## Determinism

- `compute_consensus_rts()` sorts output by `(is_decoy, modified_sequence)` before returning
- `plan_reconciliation()` iterates entries in deterministic order (per-file, per-entry index)
- Re-scoring is done sequentially (`iter_mut()`, not `par_iter_mut()`) to avoid both nondeterministic floating-point accumulation and OOM from parallel spectra loading

---

## Implementation

| File | Function | Purpose |
|------|----------|---------|
| `crates/osprey/src/reconciliation.rs` | `compute_consensus_rts()` | Collect detections from FdrEntry stubs, compute weighted median library RTs |
| `crates/osprey/src/reconciliation.rs` | `refit_calibration_with_consensus()` | Refit per-run LOESS from consensus points |
| `crates/osprey/src/reconciliation.rs` | `plan_reconciliation()` | Determine Keep/UseCwtPeak/ForcedIntegration per entry |
| `crates/osprey/src/reconciliation.rs` | `determine_reconcile_action()` | Single-entry action determination |
| `crates/osprey/src/pipeline.rs` | `run_search()` | Reused for re-scoring with `boundary_overrides` parameter |
| `crates/osprey/src/pipeline.rs` | `select_post_fdr_consensus()` | Multi-charge consensus target selection (merged with reconciliation targets) |
| `crates/osprey/src/pipeline.rs` | `load_cwt_candidates_from_parquet()` | Selective CWT-only Parquet column reader |
| `crates/osprey/src/pipeline.rs` | reconciliation orchestration | Build overrides, load spectra, call `run_search()`, second-pass FDR |

### Tests

| Test | Description |
|------|-------------|
| `test_weighted_median_single` | Single detection → returns that value |
| `test_weighted_median_equal_weights` | Equal weights → regular median |
| `test_weighted_median_skewed_weights` | Heavy weight pulls consensus toward that run |
| `test_simple_median_odd/even` | Median computation for odd/even counts |

---

## References

- Cross-run alignment strategy follows the general approach of DIA-NN and Spectronaut
- Weighted median for robust RT aggregation: more reliable than arithmetic mean in the presence of outlier detections
