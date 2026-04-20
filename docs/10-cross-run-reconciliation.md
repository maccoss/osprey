# Cross-Run Peak Reconciliation

Cross-run reconciliation aligns peak integration boundaries across replicate files so that the same peptide is quantified at a consistent chromatographic position in every run. It runs after the initial per-run FDR pass, before the final experiment-level FDR, and only applies to multi-file experiments.

## Problem

Without reconciliation, each file's per-window search finds peaks independently using CWT peak detection and pairwise coelution scoring. For a given peptide, the search may confidently find the correct peak in most runs but select a completely different peak — a co-eluting interferer, a noise spike, or an isomer — in one or more runs where the signal is weaker or the interference is stronger.

### RT-Penalized Peak Selection

To mitigate wrong-peak selection, CWT candidate peaks are ranked by `coelution_score * rt_penalty * intensity_weight` rather than by `coelution_score` alone.

**RT penalty** — Gaussian centered on the calibration-predicted RT:

```text
rt_penalty = exp(-residual^2 / (2 * sigma^2))
sigma      = max(5 * MAD * 1.4826, 0.1)   // from the per-file calibration
residual   = |peak_apex - expected_rt|
```

A peak at the expected position gets penalty 1.0; a peak 1 sigma away gets ~0.61; a peak 3 sigma away gets ~0.011; a peak 5 sigma away is effectively zero. At 5-sigma width the penalty is gentle on peaks with slight RT deviations from the calibration (e.g., 0.3 min off with Stellar's MAD ~0.145 gives sigma ~1.07 and penalty ~0.96) but still strongly rejects genuine wrong-peak selections (1.0 min off gives penalty ~0.65, 2.0 min off gives ~0.17).

An earlier 3-sigma version was too aggressive — it could downweight the correct peak enough that a nearby shoulder won. The 5-sigma widening in v26.3.1 preserves the interferer-rejection benefit while leaving room for legitimate peptide-specific RT deviations from the LOESS prediction.

**Intensity weight** — log-scaled peak height:

```text
intensity_weight = log(1 + apex_intensity)
```

Peak intensity is a multiplicative factor in the score, always contributing alongside coelution and the RT penalty. Log scaling compresses the dynamic range (a 100× intensity difference is only ~1.5× in score) so intensity doesn't dominate the coelution ranking, but it's still meaningful enough to disambiguate candidates that would otherwise score similarly — most importantly, a narrow low-intensity shoulder of a main peak.

Without these modifiers, a common failure mode on multi-replicate experiments is: an interferer with slightly better fragment co-elution than the correct peak gets selected in multiple replicates, corrupting the downstream consensus RT and causing reconciliation to "correct" the good replicates to the wrong position.

```text
Peptide PEPTIDEK (with RT penalty + intensity tiebreaker):
  file1.mzML: coelution=8.5, rt_penalty=1.0, intensity=high → wins at apex=25.3 (correct)
  file2.mzML: coelution=7.9, rt_penalty=1.0, intensity=high → wins at apex=25.4 (correct)
  file3.mzML: correct peak   coelution=3.1, rt_penalty=1.0,  intensity=medium → score dominates
              interferer     coelution=4.2, rt_penalty=0.17, intensity=medium → score=~0.7
              → correct peak wins despite lower coelution
```

The stored CWT candidates for reconciliation retain their **raw** `coelution_score` (not the RT-penalized version) since reconciliation has its own consensus-based RT tolerance logic via `plan_reconciliation`.

Cross-run reconciliation uses the high-confidence runs to establish where PEPTIDEK actually elutes, then goes back to files where the selected peak is wrong and either finds an alternate CWT candidate at that RT or imputes the integration boundaries at the expected position.

---

## Algorithm Overview

```text
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

A target **detection** (not just the peptide as a whole) qualifies for consensus RT computation when ALL of the following hold:

1. Its **per-entry precursor q-value** is ≤ `consensus_fdr` (default 0.01). This is a **hard gate** — protein FDR cannot rescue poor precursor evidence.
2. EITHER its peptide-level q-value is ≤ `consensus_fdr`, OR its first-pass protein group q-value is ≤ `protein_fdr_threshold`.

Paired decoys (DECOY_ prefix, any matching modified_sequence) are included regardless of their own q-values so the second FDR pass sees fair competition:

```text
targets: {PEPTIDEK, ANOTHERR, ...}     ← at least one detection with prec_q ≤ fdr
decoys:  {DECOY_PEPTIDEK, DECOY_ANOTHERR, ...}  ← all paired decoys
```

Both sides get independent consensus RTs — decoy consensus RTs are computed from decoy detections (not from target detections) to avoid any information leakage.

### Why Precursor-Level Evidence Is a Hard Gate

Consensus RT is driven by each surviving entry's own `apex_rt`. If the only thing qualifying a detection is its protein group being strong (but the entry itself has a weak/wrong-peak apex), that wrong apex gets pulled into the weighted median and can shift the consensus toward an interferer.

**Regression case (Stellar, DAQVVGMTTTGAAK):** this peptide's z=3 first-pass CWT consistently picked an interferer at ~8.46 in three replicates (scores −4.1 / −3.6 / −0.26, precursor q-values 0.099 / 0.091 / 0.016). The correct z=2 peak was at ~8.67 in two replicates. Under the looser earlier gate, the three z=3 wrong-peak detections were rescued into consensus by their protein's picked-protein FDR and — with `coelution_sum` as the median weight — just barely outweighed the two correct z=2 detections (total weight 18.62 vs 17.64). The consensus collapsed to the interferer position. Reconciliation then force-integrated the two correct z=2 observations to the interferer RT, and experiment-level FDR rejected them as "wrong peak."

The tightened gate drops all three z=3 wrong-peak detections (precursor q > 0.01) and the two z=2 correct detections alone anchor the consensus at 8.68. Reconciliation now keeps the two correct observations and the third sample is gap-filled (via gap-fill CWT) to the correct 8.66 peak. See the `test_consensus_rejects_low_precursor_q_despite_protein_rescue` regression test in [reconciliation.rs](../crates/osprey/src/reconciliation.rs).

### Protein FDR Rescue Gate (borderline peptide evidence)

The protein-FDR rescue path lets peptides from strong proteins enter consensus when their *peptide*-level q-value is borderline but their own detection is decent. Specifically: a detection with `run_precursor_qvalue ≤ consensus_fdr` AND `run_peptide_qvalue > consensus_fdr` can still qualify if its first-pass protein group passes (`run_protein_qvalue ≤ protein_fdr_threshold`). This retains the original intent — boost borderline peptides from strong proteins — without allowing protein FDR to override individual detection quality. Protein FDR runs before reconciliation (see [07-fdr-control.md](07-fdr-control.md)), so `run_protein_qvalue` is already populated. When `--protein-fdr` is disabled (threshold 0.0), the rescue branch is off.

---

## Step 2: Consensus Library RT

For each consensus peptide, its qualifying detections across all files are collected:

```text
(file_name, apex_rt, score, peak_width, coelution_sum)
```

Each measured `apex_rt` is mapped back to **library RT space** using the run's inverse RT calibration:

```text
library_rt = RTCalibration::inverse_predict(apex_rt)
```

The **consensus library RT** is the **weighted median** of these library RT values, with weights derived from the SVM discriminant score:

```text
weight = max(1e-6, 1 / (1 + exp(-score)))     // sigmoid(score), floor prevents zero
```

The SVM score is our strongest per-detection quality signal — it's what FDR itself ranks on. Mapping it through `sigmoid` produces a smoothly monotonic weight in `(0, 1)`: a detection with positive score (target-like) dominates the weighted median; a detection with negative score (noise/interferer-like) contributes near-zero weight and cannot poison the median.

Per-peptide peak widths are aggregated the same way (weighted median with sigmoid-score weights), so the `median_peak_width` used for gap-fill / forced integration reflects high-quality detections' peak shapes and not noisy wide fallbacks.

An additional `coelution_sum > 0` sanity filter rejects anti-correlated "noise integration" detections (e.g., forced integrations at an empty RT window where fragments happen to anti-correlate). They never reach the weighted median.

```text
PeptideConsensusRT {
    modified_sequence: String,
    is_decoy: bool,
    consensus_library_rt: f64,    // weighted median in library RT space
    median_peak_width: f64,        // weighted median peak width across detections (minutes)
    n_runs_detected: usize,
    apex_library_rt_mad: Option<f64>,  // per-peptide MAD in library RT space (≥ 3 detections)
}
```

Working in library RT space is important: it decouples the consensus from run-specific RT drift. The consensus library RT can then be mapped into any run's measured RT space using that run's calibration.

### Why SVM Score, Not `coelution_sum`

`coelution_sum` (sum of pairwise fragment Pearson correlations at the selected peak) is a single feature. It captures whether *the chosen apex has co-eluting fragments*, but says nothing about whether *this apex is the correct peak for this peptide*. A strong interferer with 5-6 coincidentally co-eluting fragments can have `coelution_sum > 5` at a completely wrong RT. The SVM is trained on all 21 PIN features (coelution, mass accuracy, xcorr, peak shape, RT deviation, median polish, etc.) explicitly to separate correct detections from wrong-peak targets — so it is by construction a better weight.

**Regression case (Stellar, DAQVVGMTTTGAAK):** under coelution_sum weighting, three z=3 wrong-peak detections had weights 5.48 / 7.11 / 6.03 (total 18.62), barely outweighing the two correct z=2 detections with weights 10.25 and 7.39 (total 17.64). Under sigmoid-score weighting (given the z=3 scores of −4.1 / −3.6 / −0.26 and z=2 scores of +0.32 / +1.61), the wrong side totals ~0.46 and the correct side totals ~1.41 — a 3x margin. See the `test_consensus_weighting_downweights_negative_score_detections` regression test.

### Multi-Charge Reconciliation Is Inherently Cross-Run

Consensus grouping keys are `(modified_sequence, is_decoy)` — **charge is not part of the key**. All detections of a peptide across every file and every charge state contribute to one shared `PeptideConsensusRT`. This has two important consequences:

1. **A peptide's charge states reinforce each other across files.** If PEPTIDEK is detected as 2⁺ in file1 and 3⁺ in file2, both detections contribute (weighted by their respective `coelution_sum`) to a single consensus library RT. During planning, each `(peptide, charge, file)` triple then derives its own expected measured RT from that shared consensus via the file's refined calibration.

2. **Gas-phase fractionation (GPF) is handled automatically.** In GPF each file covers a disjoint precursor m/z range, so a peptide's different charge states land in different files — its 2⁺ may only be visible in file1 and its 3⁺ only in file3. Because consensus is charge-agnostic, file1's 2⁺ detection still anchors the consensus RT that file3's 3⁺ detection is checked against, and vice versa. There is no separate "intra-file multi-charge consensus" step in GPF mode — the cross-run consensus IS the multi-charge alignment.

The per-file `select_post_fdr_consensus()` function that handles multi-charge alignment within a single replicate (standard DIA case) is effectively a no-op in GPF because each file holds at most one charge state per peptide. That is correct behavior, not a bug: the inter-replicate consensus path covers the GPF case by construction.

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

```text
expected_rt = refined_calibration.predict(consensus_library_rt)
```

Three actions are possible based on **apex proximity** to the expected RT:

| Action | Condition | What happens |
|--------|-----------|--------------|
| **Keep** | `\|apex_rt - expected_rt\| <= rt_tolerance` | No change |
| **UseCwtPeak** | A CWT candidate has `\|cand.apex_rt - expected_rt\| <= rt_tolerance` | Switch to the closest-apex candidate |
| **ForcedIntegration** | Neither current peak nor any CWT candidate has apex within tolerance | Integrate at `expected_rt ± median_peak_width/2` |

### RT Tolerance: Global Within-Peptide MAD

The `rt_tolerance` used by `plan_reconciliation` is now driven by the **within-peptide** RT reproducibility of the experiment, not the cross-peptide calibration residuals:

```text
global_within_peptide_mad = median over all target peptides (with n_runs_detected ≥ 3) of:
    median_i( | library_rt_i − consensus_library_rt | )

rt_tolerance = max(0.1, 3.0 × global_within_peptide_mad × 1.4826)
rt_tolerance = rt_tolerance.min(file_calibration_ceiling)  // safety ceiling
```

Each peptide's own `apex_library_rt_mad` is the median absolute deviation of its library-space apex RTs across replicates. For a well-aligned peptide detected 3+ times, this MAD reflects only the LC/instrument reproducibility floor (drift between replicates), not any cross-peptide calibration residual. The **global** tolerance is the median across all peptides' MADs — a robust experiment-wide estimate that individual noisy peptides can barely move. Applied uniformly to every peptide.

**Why global, not per-peptide.** Individual per-peptide MADs from 3-5 replicates are very noisy estimators; the median across thousands of peptides is much more stable. After LOESS alignment, within-peptide scatter is approximately peptide-independent — it reflects the instrument/LC floor, not anything unique to the peptide — so the global median is an appropriate summary. Per-peptide MADs are still computed and stored on `PeptideConsensusRT` (and included in trace output for diagnostic use) but not used as the primary tolerance.

**Why not cross-peptide calibration MAD.** The per-file LOESS fit's MAD conflates two very different things: (1) genuine LC drift between replicates (what we care about for reconciliation tolerance), and (2) peptide-to-peptide deviation from the LOESS model at each RT (which has nothing to do with whether a single peptide's apex drifted between runs). On a well-aligned experiment the calibration MAD is typically 3-5x larger than the within-peptide MAD. Using it as the tolerance makes reconciliation blind to wrong-peak selections that are within 0.2-0.3 min of the correct peak, because those residuals fall within the calibration's cross-peptide noise envelope.

**Concrete example.** On a Stellar 3-replicate HeLa experiment the cross-peptide calibration MAD was ~0.07 min, giving a tolerance of ~0.31 min. The global within-peptide MAD was **0.0207 min**, giving a tolerance of ~0.1 min (after the floor). The previous tolerance of 0.31 min Kept the DAQVVGMTTTGAAK wrong-peak detections 0.21 min off-consensus; the new tolerance of 0.1 min correctly flags them for re-scoring.

**The per-file calibration MAD is retained as a safety ceiling.** If a pathological global MAD is somehow tighter than scan resolution (e.g., degenerate dataset with all apexes rounded to the same scan), the calibration-MAD-based ceiling guards against a sub-threshold tolerance that would force every peak to get re-scored.

A minimum floor of 0.1 min prevents being too strict when the global MAD is exceptionally tight. (This floor is intentionally coarse for now — it exceeds the raw global MAD × 3-sigma on well-aligned experiments, so in practice the floor is often the binding constraint. Tightening it is a future refinement if warranted by specific datasets.)

**No per-point local tolerance.** An earlier version used `local_tolerance(query_rt)` which interpolated the absolute residual from the nearest training points at the query RT. That created a self-fulfilling prophecy: a peptide with a wrong apex RT contributed a huge residual to the calibration training set, which inflated the local tolerance at that RT, which then allowed the wrong-RT detection to pass the apex proximity check. Using a single global number eliminates this feedback loop.

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

## Gap-Fill

Reconciliation planning only acts on entries that were actually scored in each file during the first pass. If a precursor passed FDR in one replicate but was never scored in another (either because its signal fell below the pre-filter or because the scoring pass simply missed it), reconciliation alone would leave that file with no integration for that precursor, and the blib output would be missing a per-file boundary. Gap-fill closes this hole.

### Identification

`identify_gap_fill_targets()` walks the set of precursors `(modified_sequence, charge)` passing run-level FDR in at least one replicate, and for each file determines which of those precursors have no corresponding FdrEntry. For each missing precursor in each file it emits a `GapFillTarget` carrying:

- `target_entry_id` / `decoy_entry_id` — library IDs for the target and its paired decoy
- `expected_rt` — consensus library RT mapped through this file's refined calibration
- `half_width` — half the consensus `median_peak_width`
- `modified_sequence`, `charge` — for downstream routing

Both target and decoy are included together so the second-pass FDR pass sees symmetric competition.

### Isolation Window m/z Filter (GPF-Aware)

Naively forcing an integration at every missing precursor is wrong for gas-phase fractionation. In GPF each file covers a disjoint precursor m/z range, so for many "missing" precursors the peptide is not merely below the detection threshold — the file **physically cannot observe** that m/z because no isolation window selects it for MS2. Forcing integration there would land on noise and inflate the null distribution with meaningless features.

Gap-fill therefore filters candidates by per-file isolation window coverage. Each file's isolation scheme (captured from its first MS2 cycle during the calibration step) is stored as a list of `[lower_mz, upper_mz]` intervals. A candidate survives the filter only if its library precursor m/z falls inside at least one of the target file's intervals:

```text
precursor_mz = library[target_entry_id].precursor_mz
in_range = any(lo <= precursor_mz < hi for (lo, hi) in file_windows)
if not in_range:
    skip this candidate (filtered_by_mz += 1)
```

Skipped candidates are counted and logged at INFO level:

```text
Gap-fill: 18423 candidates skipped because precursor m/z is outside
the file's isolation windows (GPF or disjoint m/z ranges)
```

When a file has no stored isolation scheme (e.g., older calibration caches), filtering is disabled for that file as a graceful fallback — the pipeline reverts to the old "fill every missing precursor" behavior rather than silently dropping everything.

The same filter also helps non-GPF acquisitions where different replicates happen to have partially disjoint isolation schemes: a precursor whose m/z lands in a gap between the target file's windows is correctly skipped even in standard DIA runs.

### Scoring

Surviving gap-fill targets are re-scored via `run_search()` with `boundary_overrides` set from `(expected_rt - half_width, expected_rt, expected_rt + half_width)`, using the same parallel-window path as reconciliation re-scoring. New entries are written back into `per_file_entries` and feed into the second-pass FDR computation.

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

```text
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
| `crates/osprey/src/reconciliation.rs` | `compute_consensus_rts()` | Collect detections from FdrEntry stubs, compute weighted median library RTs (charge-agnostic grouping) |
| `crates/osprey/src/reconciliation.rs` | `refit_calibration_with_consensus()` | Refit per-run LOESS from consensus points |
| `crates/osprey/src/reconciliation.rs` | `plan_reconciliation()` | Determine Keep/UseCwtPeak/ForcedIntegration per entry |
| `crates/osprey/src/reconciliation.rs` | `determine_reconcile_action()` | Single-entry action determination |
| `crates/osprey/src/reconciliation.rs` | `identify_gap_fill_targets()` | Find precursors missing from each file; filter by per-file isolation window m/z coverage |
| `crates/osprey/src/pipeline.rs` | `run_search()` | Reused for re-scoring with `boundary_overrides` parameter |
| `crates/osprey/src/pipeline.rs` | `select_post_fdr_consensus()` | Intra-file multi-charge consensus (no-op in GPF; cross-run consensus handles GPF instead) |
| `crates/osprey/src/pipeline.rs` | `load_cwt_candidates_from_parquet()` | Selective CWT-only Parquet column reader |
| `crates/osprey/src/pipeline.rs` | reconciliation orchestration | Build overrides, load spectra, call `run_search()`, second-pass FDR |

### Tests

| Test | Description |
|------|-------------|
| `test_weighted_median_single` | Single detection → returns that value |
| `test_weighted_median_equal_weights` | Equal weights → regular median |
| `test_weighted_median_skewed_weights` | Heavy weight pulls consensus toward that run |
| `test_simple_median_odd/even` | Median computation for odd/even counts |
| `test_gap_fill_identifies_missing_precursors` | Precursor in file1 but not file2 → file2 gets a gap-fill target |
| `test_gap_fill_filters_precursors_outside_isolation_windows` | GPF disjoint ranges: precursor visible in file1 is NOT gap-filled in file2 when its m/z is outside file2's isolation windows |
| `test_gap_fill_allows_precursors_inside_isolation_windows` | Overlapping replicate-style ranges: same precursor IS gap-filled in file2 when m/z falls inside the shared window |
| `test_consensus_rejects_low_precursor_q_despite_protein_rescue` | Three z=3 wrong-peak detections (precursor q > 0.01) rescued by protein FDR no longer enter consensus; two correct z=2 detections alone anchor the consensus on the correct peak |
| `test_consensus_weighting_downweights_negative_score_detections` | Two detections with identical precursor q but scores −4 vs +1.5 — weighted median collapses to the high-score detection (sigmoid-score weight, not coelution_sum weight) |
| `test_plan_reconciliation_global_mad_is_median_across_peptides` | Five peptides with MADs 0.01-0.50 yield median=0.05 and tolerance=0.222 min; an outlier peptide (MAD 0.50) does not inflate the global tolerance |
| `test_plan_reconciliation_uses_global_within_peptide_mad_tighter_than_calibration` | Peptide MAD (0.02) much tighter than file calibration MAD (0.1) — tolerance driven by global MAD, not calibration ceiling |
| `test_plan_reconciliation_tolerance_matches_expected_mad_formula` | Global MAD of 0.1 produces tolerance 3 × 0.1 × 1.4826 ≈ 0.445 min; entry 0.4 min off passes, entry 0.5 min off is re-scored |

---

## References

- Cross-run alignment strategy follows the general approach of DIA-NN and Spectronaut
- Weighted median for robust RT aggregation: more reliable than arithmetic mean in the presence of outlier detections
