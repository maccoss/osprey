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
  5. Re-score reconciled entries  — reload full entries from Parquet, recompute features
  6. Second-pass FDR              — final experiment-level q-values
```

---

## Step 1: Consensus Peptide Selection

Target peptides passing `consensus_fdr` (default 1%) at **experiment-level FDR** are collected. Their paired decoys (DECOY_ prefix) are included so target-decoy competition remains fair in the second FDR pass.

```
targets: {PEPTIDEK, ANOTHERR, ...}     ← pass experiment-level FDR
decoys:  {DECOY_PEPTIDEK, DECOY_ANOTHERR, ...}  ← paired decoys included
```

Both sides get independent consensus RTs — decoy consensus RTs are computed from decoy detections (not from target detections) to avoid any information leakage.

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

Three actions are possible:

| Action | Condition | What happens |
|--------|-----------|--------------|
| **Keep** | Current peak `[start_rt, end_rt]` contains `expected_rt` | No change |
| **UseCwtPeak** | An alternate stored CWT candidate contains `expected_rt` | Switch to that candidate's boundaries |
| **ForcedIntegration** | No CWT candidate contains `expected_rt` | Integrate at `expected_rt ± median_peak_width/2` |

The CWT candidates are stored during the initial search (configurable, default 5 per precursor) and persisted in the per-file Parquet score cache. During reconciliation planning, only the CWT candidate column is loaded selectively via `load_cwt_candidates_from_parquet()`, avoiding the cost of reloading full entries with all features and fragment data. This lets reconciliation switch to an alternate peak that was already detected but not selected as the best, without requiring a full re-extraction.

`ForcedIntegration` is a last resort — it integrates at the expected position even if no CWT peak was found there, using the median peak width from all detections as the integration window.

---

## Step 5: Re-Scoring

Entries with `UseCwtPeak` or `ForcedIntegration` actions are re-scored at their new boundaries:
- Full `CoelutionScoredEntry` data is reloaded from the per-file Parquet cache
- Files are processed **sequentially** (not in parallel) to limit peak memory — each file loads spectra (~1-2 GB) and full Parquet entries (~1.4 GB), so parallel loading would OOM on large experiments
- All features are recomputed at the new peak boundaries via `compute_features_at_peak()`
- Updated entries are converted back to FdrEntry stubs for the second-pass FDR
- Entries that cannot be re-scored (no spectral data at the expected RT) are dropped

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

Per-run refined calibration → expected measured RT:
  file1: predict(32.1) = 25.3 → current peak [25.0, 25.9] contains 25.3 → Keep
  file2: predict(32.1) = 25.4 → current peak [25.0, 25.9] contains 25.4 → Keep
  file3: predict(32.1) = 25.2 → current peak [31.2, 32.1] does NOT contain 25.2
           → check stored CWT candidates...
             candidate at [25.0, 25.7] contains 25.2 → UseCwtPeak
             (CWT had found the correct peak but scored it lower than the interferer)

Result:
  file1, file2: unchanged — kept their correct peaks
  file3: switched to the stored CWT candidate at RT≈25.2 (the true elution time)
         features recomputed at new boundaries
```

If file3 had no CWT candidate near the expected RT (e.g., the peptide was truly not present or below detection), the outcome would be `ForcedIntegration` — boundaries placed at `expected_rt ± median_peak_width/2`. This provides a quantification window for Skyline even when the peptide was not detected in that run.

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
| `crates/osprey/src/pipeline.rs` | `load_cwt_candidates_from_parquet()` | Selective CWT-only Parquet column reader |
| `crates/osprey/src/pipeline.rs` | reconciliation orchestration | Load data from Parquet, re-score, second-pass FDR |

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
