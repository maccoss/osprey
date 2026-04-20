# Peptide Trace Diagnostics

Osprey includes an env-var gated trace facility for diagnosing per-peptide behavior at every major pipeline stage. It's useful for investigating single-peptide questions that would otherwise be buried in GB-scale Parquet caches: why a specific peptide got integrated at the wrong RT, whether consensus was poisoned, which CWT candidates were available, how FDR ranked a given charge state, and so on.

## Enabling

Set `OSPREY_TRACE_PEPTIDE` to a modified_sequence (or a comma-separated list) before running Osprey:

```bash
OSPREY_TRACE_PEPTIDE=DAQVVGMTTTGAAK osprey -i *.mzML -l library.tsv -o results.blib --verbose 2>&1 | tee trace.log
grep '\[trace\]' trace.log
```

Multiple peptides can be traced in one run:

```bash
OSPREY_TRACE_PEPTIDE=AMLNNVVRPR,DAQVVGMTTTGAAK,CPLNEEVIVQAR osprey ...
```

The match is on the **bare modified_sequence** (no `DECOY_` prefix, no flanking residues). Paired decoys (`DECOY_<target>`) are traced automatically when the target is listed.

Trace is zero overhead when the env var is not set (`OnceLock`-guarded probe).

## Invalidating Caches for a Clean Run

The trace only fires on code paths that actually execute. Parquet `.scores.parquet` caches, FDR score sidecars (`.1st-pass.fdr_scores.bin`, `.2nd-pass.fdr_scores.bin`), and the skip-Percolator / skip-reconciliation fast paths will all shortcut past the traced logic. To get a full trace, delete the per-file caches first:

```bash
rm -v <prefix>.scores.parquet
rm -v <prefix>.{1st,2nd}-pass.fdr_scores.bin   # optional but recommended
```

Calibration JSON (`<prefix>.calibration.json`) and binary spectra cache (`<prefix>.spectra.bin`) can be kept — they do not affect the trace path.

## Trace Format

All trace lines are emitted at `log::info!` and share a `[trace]` prefix. They're interleaved with normal log output; grep `'\[trace\]'` to isolate them. Trace lines appear at five pipeline stages:

### 1. First-pass CWT peak scoring

For every entry of a traced peptide in every file, after CWT candidates are sorted by penalized score:

```text
[trace] <modified_sequence> z=<charge> <TARGET|DECOY> <file> expected_rt=<val> rt_sigma=<val> — CWT candidates (sorted by penalized score):
[trace]   rank 0: apex=<val> (rt_res=<val>) coelution=<val> rt_penalty=<val> intensity=<val> int_weight=<val> penalized=<val>
[trace]   rank 1: ...
...
```

`rank 0` is the winner — it becomes the entry's `apex_rt` in first-pass. `rank 1`-`rank N` are the other stored top-N candidates that reconciliation may later pick from.

Fields:

- `expected_rt` — first-pass (discovery) calibration prediction for the library RT
- `rt_sigma` — width of the Gaussian RT penalty (5 × calibration MAD × 1.4826)
- `apex` — candidate apex RT in measured space
- `rt_res` — signed residual `apex - expected_rt`
- `coelution` — mean pairwise fragment correlation at the candidate peak (unpenalized)
- `rt_penalty` — `exp(-rt_res² / (2·rt_sigma²))`
- `intensity` — raw apex intensity
- `int_weight` — `log(1 + intensity)` (intensity tiebreaker)
- `penalized` — `coelution × rt_penalty × int_weight` (the ranking score)

Use this to answer "did CWT actually find the correct peak?" — if so, it'll appear as rank ≥ 1 with visible `coelution`/`rt_residual` context; if the winner (rank 0) was a wrong peak, the trace shows exactly how much higher its penalized score was and which factor dominated (coelution, intensity, or RT proximity).

### 2. `compute_consensus_rts`

Once per qualifying peptide (target), plus paired decoy, after weighted-median consensus is computed:

```text
[trace] consensus <modified_sequence> <TARGET|DECOY> — detections (<n> runs contributing):
[trace]   file=<file> apex_rt=<val> lib_rt=<val> score=<val> weight=<val> coelution_sum=<val> peak_width=<val>
[trace]   ...
[trace] consensus <modified_sequence> <TARGET|DECOY> → consensus_library_rt=<val>, median_peak_width=<val>
```

`weight` is `sigmoid(score)` clamped to `[1e-6, 1)` — this is the weight used in the weighted median. `score` is the first-pass SVM discriminant. Use this to see whether consensus was built from the right detections and whether a single strong detection or a group of weaker ones anchored it.

### 3. `plan_reconciliation`

Once per entry of a traced peptide/charge/file that survives to reconciliation planning:

```text
[trace] plan <modified_sequence> z=<charge> <TARGET|DECOY> file=<file> current_apex=<val> expected_rt=<val> residual=<val> tolerance=<val> (global_MAD=<val>, this_peptide_MAD=<val> from n=<val>, ceiling=<val>) n_cwt=<val> → <action>
[trace]   cwt[0]: apex=<val> (res=<val>) coelution=<val> area=<val>
[trace]   cwt[1]: ...
...
```

`action` is one of:

- `Keep` — apex within tolerance of consensus-predicted RT; entry unchanged
- `UseCwtPeak(idx=<i>, apex=<val>)` — a stored CWT candidate's apex is closer to expected_rt; reconciliation re-scores at that candidate
- `ForcedIntegration(center=<val>, half_width=<val>)` — no suitable CWT candidate exists; integrate at `expected_rt ± half_width`

Fields:

- `current_apex` — entry's existing apex_rt going into reconciliation
- `expected_rt` — refined calibration's prediction from `consensus_library_rt`
- `residual` — signed `current_apex - expected_rt`
- `tolerance` — the per-peptide RT tolerance (derived from `global_MAD`, floored at 0.1 min, capped at `ceiling`)
- `global_MAD` — median of within-peptide MADs across the experiment (library RT space)
- `this_peptide_MAD` — this peptide's own library-RT MAD (diagnostic; not used for tolerance derivation)
- `ceiling` — per-file calibration-MAD-based safety ceiling
- `n_cwt` — number of stored CWT candidates available for UseCwtPeak
- `cwt[i]` — each stored candidate's apex, residual, coelution, area

Use this to see exactly why an entry was Kept vs re-scored and, for non-Keep entries, whether a stored CWT candidate matched the consensus RT.

### 4. Gap-fill identification

Emitted once per (peptide, charge, file) where the precursor passed FDR in some other replicate but was missing from this file:

```text
[trace] gap-fill: <modified_sequence> z=<charge> file=<file> expected_rt=<val> half_width=<val> (missing from this file, passed in another)
```

If the peptide's m/z is filtered out (GPF isolation window mismatch) or the precursor wasn't missing, no gap-fill line is emitted for that file. Subsequent behavior (gap-fill CWT Pass 1 finding a natural peak vs Pass 2 forced integration) is reflected in the second-pass FDR trace.

### 5. First-pass and second-pass FDR q-values

After each FDR pass completes, every traced entry in every file gets one line:

```text
[trace] fdr(<stage>) <modified_sequence> z=<charge> <TARGET|DECOY> file=<file> apex=<val> score=<val> pep=<val> run_q(prec/pep/eff)=<v>/<v>/<v> exp_q(prec/pep/eff)=<v>/<v>/<v>
```

where `<stage>` is `first-pass` or `second-pass`. `run_q` triple is (precursor q, peptide q, effective = max of both). `exp_q` is the same at experiment scope.

Use this to see the peptide's post-FDR disposition at both passes, and to identify entries whose q-values changed between passes (typically because reconciliation moved them to a different peak with better / worse features).

## Typical Diagnostic Workflows

### "Why was this peptide integrated at RT X in file Y when the real peak is at Z?"

1. Enable trace for the peptide, clear Parquet caches, rerun.
2. Check first-pass CWT trace for file Y — if rank 0 = wrong peak but rank 1-4 include the right one, it's an RT penalty / coelution problem at first pass.
3. Check consensus trace — did only wrong detections qualify? Did the weighted median land on the wrong RT?
4. Check plan trace for file Y — was the wrong apex within tolerance (Kept) or was reconciliation force-integrating at a bad consensus?

### "Why did this peptide fail FDR?"

1. Check fdr(first-pass) and fdr(second-pass) lines — if score went from positive to negative between passes, reconciliation moved the entry to a location with worse features (a symptom of poisoned consensus).
2. Look at consensus trace: was it built from only 1-2 detections, or did the weighted median land on a bad peak?

### "Why is this peptide inconsistent across replicates in the blib?"

1. Check plan traces across all files — mix of Keep / UseCwtPeak / ForcedIntegration?
2. Check gap-fill trace — files that didn't have first-pass detections get forced-integrated at `expected_rt`; if that's wrong, gap-fill CWT Pass 1 may have picked a natural peak at a different RT instead.
3. Cross-reference with the second-pass fdr trace `apex` values — they show where each observation actually ended up.

## Implementation

| File | Purpose |
|------|---------|
| `crates/osprey/src/trace.rs` | `trace_set()` (OnceLock of env-var-split targets) and `is_traced(modified_sequence)` predicate; also the `log_fdr_qvalues()` helper called after each FDR pass |
| `crates/osprey/src/pipeline.rs` (inside `run_search` scoring loop, after `scored_candidates.sort_by`) | First-pass CWT candidate trace |
| `crates/osprey/src/reconciliation.rs` (`compute_consensus_rts`) | Consensus detection + weight trace |
| `crates/osprey/src/reconciliation.rs` (`plan_reconciliation`) | Per-entry plan trace with tolerance derivation and stored CWT candidates |
| `crates/osprey/src/reconciliation.rs` (`identify_gap_fill_targets`) | Gap-fill target trace |
| `crates/osprey/src/pipeline.rs` (after each FDR pass) | `trace::log_fdr_qvalues(per_file_entries, "first-pass" / "second-pass")` |

Trace output lives alongside normal logging (terminal + log file) so it's preserved across sessions. Zero allocations or function calls for non-traced peptides.
