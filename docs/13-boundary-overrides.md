# Boundary Overrides & Forced Integration

This document describes the `boundary_overrides` mechanism in `run_search()` that enables multi-charge consensus and cross-run reconciliation to re-score entries at specific RT boundaries without running CWT peak detection.

## Overview

`run_search()` accepts an optional `boundary_overrides: Option<&HashMap<u32, (f64, f64, f64)>>` parameter. When present, entries whose `entry_id` appears in the map skip peak detection and instead score at the given `(apex_rt, start_rt, end_rt)` boundaries.

```rust
fn run_search(
    library: &[LibraryEntry],
    spectra: &[Spectrum],
    ms1_index: &MS1Index,
    calibration: Option<&CalibrationParams>,
    rt_calibration: Option<&RTCalibration>,
    config: &OspreyConfig,
    file_name: &str,
    boundary_overrides: Option<&HashMap<u32, (f64, f64, f64)>>,
) -> Result<Vec<CoelutionScoredEntry>>
```

This design means the first-pass search and all re-scoring paths (multi-charge consensus, cross-run reconciliation) use **exactly the same function** — the only difference is whether boundaries come from CWT peak detection or from the caller.

## Three Override Types

### 1. Multi-Charge Consensus (Forced Boundaries)

When a peptide has multiple charge states at different peaks, the best SVM-scoring FDR-passing charge state defines the consensus. Other charge states are re-scored at the consensus boundaries.

```
Source:    select_post_fdr_consensus()
Boundaries: (consensus_apex, consensus_start, consensus_end)
```

### 2. UseCwtPeak (Stored CWT Candidate)

During reconciliation, if the current peak doesn't contain the expected RT but a stored CWT candidate does, the system switches to that candidate's boundaries.

```
Source:    plan_reconciliation() → ReconcileAction::UseCwtPeak
Boundaries: (candidate_apex_rt, candidate_start_rt, candidate_end_rt)
```

CWT candidates are stored during the initial search (up to `n_cwt_candidates` per precursor, default 5) and persisted in the per-file Parquet cache. They are loaded selectively via `load_cwt_candidates_from_parquet()`.

### 3. ForcedIntegration (Imputed Boundaries)

When no CWT candidate contains the expected RT, boundaries are imputed from the consensus RT and median peak width:

```
Source:    plan_reconciliation() → ReconcileAction::ForcedIntegration
Boundaries: (expected_rt, expected_rt - half_width, expected_rt + half_width)
```

Where:
- `expected_rt = refined_calibration.predict(consensus_library_rt)` — the predicted measured RT in this run
- `half_width = median_peak_width / 2` — half the median peak width across all confident detections of this peptide

This is a last resort — it integrates at the expected position even if no chromatographic peak was detected there.

## How Overrides Work Inside `run_search()`

### Step 1: Override Detection

For each candidate entry in a window, the override is checked:

```rust
let has_override = boundary_overrides
    .and_then(|m| m.get(&entry.id))
    .copied();
```

### Step 2: RT Range Selection

Entries with overrides use the override boundaries ± a margin for spectrum gathering (instead of the expected RT ± RT tolerance):

```rust
if let Some((_, target_start, target_end)) = has_override {
    let peak_width = (target_end - target_start).max(0.1);
    let margin = peak_width.max(0.2);
    let rt_lo = target_start - margin;
    let rt_hi = target_end + margin;
    // Gather spectra in [rt_lo, rt_hi]
}
```

The margin ensures enough spectra are gathered around the boundaries for XIC extraction.

### Step 3: Pre-Filter Skip

Override entries skip the signal pre-filter entirely — the caller has determined these boundaries should be scored regardless of initial signal:

```rust
if has_override.is_none() && config.prefilter_enabled {
    // Pre-filter logic (only for non-override entries)
}
```

### Step 4: RT → XIC Index Mapping

Instead of CWT peak detection, override entries map the target RTs to XIC scan indices using binary search:

```rust
let start_index = ref_xic
    .partition_point(|(rt, _)| *rt < target_start)
    .saturating_sub(1)
    .min(ref_xic.len() - 1);

let end_index = ref_xic
    .partition_point(|(rt, _)| *rt < target_end)
    .min(ref_xic.len() - 1);

let apex_index = ref_xic
    .partition_point(|(rt, _)| *rt < target_apex)
    .min(ref_xic.len() - 1);

// Refine apex to nearest neighbor
let apex_index = if apex_index > 0
    && (ref_xic[apex_index - 1].0 - target_apex).abs()
        < (ref_xic[apex_index].0 - target_apex).abs()
{
    apex_index - 1
} else {
    apex_index
};
let apex_index = apex_index.clamp(start_index, end_index);
```

This gives O(log n) mapping from RT to XIC indices, compared to CWT which is O(n × scales).

### Step 5: Feature Computation

With peak bounds determined, `compute_features_at_peak()` computes all 21 PIN features — identical to the first-pass path. The `FeatureComputeContext` includes the per-window preprocessed XCorr vectors, so no redundant preprocessing occurs.

### Early Return on Insufficient Data

If the mapped boundaries are too narrow (≤ 1 scan between start and end) or the XIC has < 3 points, the entry is dropped:

```rust
if ref_xic.len() < 3 { return None; }
if end_index <= start_index + 1 { return None; }
```

## Merging Consensus and Reconciliation Targets

In multi-file experiments, both multi-charge consensus and cross-run reconciliation may identify entries for re-scoring. These are merged into a single `boundary_overrides` map per file:

1. Multi-charge consensus targets are collected first
2. Reconciliation targets are added — **reconciliation wins on conflict** (it has more information: refined calibration + consensus RT from all runs)
3. A single `run_search()` call per file processes all override entries together

This avoids loading spectra twice per file.

## Why This Design?

### Same Code Path = Same Features

By using the same `run_search()` and `compute_features_at_peak()` for both first-pass and re-scoring, the features are computed on exactly the same scale. This is critical because the SVM/Percolator model trained on first-pass features is also used to score re-scored entries. If re-scoring used a different code path, subtle differences in feature computation could bias the scores.

### Parallel Window Processing

Re-scoring benefits from `run_search()`'s parallel window processing via rayon. Each isolation window is processed independently, and entries in the same window share:
- Grouped spectra
- Pre-extracted fragment XICs
- Per-window XCorr preprocessing (see [XCorr Scoring](04-xcorr-scoring.md#per-window-preprocessing-optimization))

### Performance vs Previous Design

The previous design used a dedicated sequential re-scoring loop (`rescore_for_reconciliation()`) that processed entries one at a time. This was ~50x slower because it missed window-level parallelism and redundant XCorr preprocessing. The current design reusing `run_search()` with overrides matches first-pass performance.

---

## Implementation

| File | Function | Purpose |
|------|----------|---------|
| `crates/osprey/src/pipeline.rs` | `run_search()` | Main search function with optional boundary_overrides |
| `crates/osprey/src/pipeline.rs` | `select_post_fdr_consensus()` | Identify multi-charge consensus targets |
| `crates/osprey/src/pipeline.rs` | `compute_features_at_peak()` | Feature computation (shared by all paths) |
| `crates/osprey/src/pipeline.rs` | `FeatureComputeContext` | Read-only context for feature computation |
| `crates/osprey/src/reconciliation.rs` | `plan_reconciliation()` | Generate ReconcileAction per entry |
| `crates/osprey/src/reconciliation.rs` | `determine_reconcile_action()` | Keep / UseCwtPeak / ForcedIntegration decision |

### Tests

| Test | Description |
|------|-------------|
| `test_run_search_override_returns_entry_at_specified_boundaries` | Override produces result at given boundaries |
| `test_run_search_override_skips_prefilter` | Override bypasses pre-filter for low-signal entries |
| `test_run_search_mixed_override_and_cwt` | Override and CWT entries coexist in same window |
| `test_run_search_no_override_normal_detection` | `None` overrides use normal CWT path |
| `test_run_search_override_narrow_boundaries` | Very tight override boundaries handled correctly |

---

## References

- See [Cross-Run Reconciliation](10-cross-run-reconciliation.md) for the full reconciliation algorithm
- See [Multi-Charge Consensus](06-multi-charge-consensus.md) for the consensus selection logic
- See [RT Alignment & Consensus Library RT](14-rt-alignment.md) for how expected RTs are computed
- See [XCorr Scoring](04-xcorr-scoring.md#per-window-preprocessing-optimization) for per-window preprocessing
