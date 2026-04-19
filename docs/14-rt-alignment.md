# RT Alignment & Consensus Library RT

This document describes how Osprey aligns retention times across runs and computes a consensus library RT for each peptide, which is the foundation for cross-run reconciliation and peak imputation.

## Problem

Retention times vary across runs due to LC system drift, gradient differences, column conditioning, and sample matrix effects. A peptide might elute at RT=25.3 in file1 and RT=25.8 in file2. To reconcile peak boundaries across runs, Osprey needs a **run-independent representation** of where each peptide elutes — the consensus library RT.

## RT Spaces

Osprey works with three RT spaces:

| Space | Description | Example |
|-------|-------------|---------|
| **Library RT** | RT from the spectral library (predicted or measured on a reference system) | 32.1 min |
| **Measured RT** | RT observed in a specific run | 25.3 min (file1), 25.8 min (file2) |
| **Consensus library RT** | Weighted median library RT across all high-confidence detections | 32.05 min |

The key insight: **library RT space is run-independent**. By mapping measured RTs back to library RT space, aggregating, and then mapping back to each run's measured RT space, we can predict where a peptide should elute in any run.

## LOWESS RT Calibration

### Initial Calibration (Phase 2)

During calibration discovery, Osprey fits a LOESS (locally estimated scatterplot smoothing) curve mapping library RT → measured RT:

```text
calibration.predict(library_rt) → expected_measured_rt
```

This is fit from high-confidence calibration targets (LDA scoring with target-decoy competition, 1% FDR) with:

- Bandwidth: 0.3 (default)
- Minimum points: 50
- Outlier retention: based on residual SD

### Inverse Calibration

The LOESS calibrator also supports **inverse prediction** — mapping measured RT back to library RT space:

```text
calibration.inverse_predict(measured_rt) → estimated_library_rt
```

This uses LOESS-interpolated inverse mapping (the calibration stores both the forward and inverse mappings). Inverse prediction is critical for consensus RT computation: it converts each run's measured apex RT into a comparable library RT value.

### Refined Calibration (Phase 5)

During reconciliation, the calibration is **refit** using only FDR-controlled consensus peptides:

```rust
refit_calibration_with_consensus(
    consensus_rts: &[PeptideConsensusRT],
    fdr_entries: &[(String, Vec<FdrEntry>)],
    original_calibration: &RTCalibration,
    file_name: &str,
) → RTCalibration
```

Parameters:

- Points: `(consensus_library_rt → measured_apex_rt)` for target peptides passing FDR in this run
- LOESS bandwidth: 0.3
- Outlier retention: 1.0 (no outlier removal — all points are FDR-controlled)
- Minimum points: 20 (falls back to original calibration if fewer)

This produces a tighter calibration because:

1. All calibration points are high-confidence FDR-controlled detections
2. The consensus library RTs are more consistent than first-pass library RTs
3. LOESS robustness iterations still downweight genuine outliers

## Consensus Library RT

### Algorithm

For each peptide passing experiment-level FDR at `consensus_fdr` (default 1%):

1. **Collect detections** across all files:

   ```text
   (file_name, apex_rt, coelution_sum, peak_width)
   ```

2. **Map to library RT space** using each run's inverse calibration:

   ```text
   library_rt = cal.inverse_predict(apex_rt)
   ```

3. **Compute weighted median** with `coelution_sum` as weights:

   ```text
   consensus_library_rt = weighted_median([(library_rt, coelution_sum), ...])
   ```

4. **Compute median peak width** across all detections:

   ```text
   median_peak_width = simple_median([peak_width_1, peak_width_2, ...])
   ```

### Why Weighted Median?

The **weighted median** is robust to outlier detections — if one run picked the wrong peak and mapped to a wildly different library RT, the weighted median ignores it as long as most runs agree. The `coelution_sum` weighting ensures that high-quality detections (strong fragment co-elution) have more influence than dubious ones.

Example:

```text
File1: library_rt=32.1, coelution_sum=8.5  ← strong detection
File2: library_rt=32.0, coelution_sum=7.9  ← strong detection
File3: library_rt=38.5, coelution_sum=3.1  ← wrong peak (outlier)

Weighted median ≈ 32.1  (file3's outlier doesn't move the median)
Arithmetic mean = 34.2  (outlier pulls the mean away)
```

### Output

```rust
PeptideConsensusRT {
    modified_sequence: String,
    is_decoy: bool,
    consensus_library_rt: f64,    // weighted median in library RT space
    median_peak_width: f64,        // median peak width across detections (minutes)
    n_runs_detected: usize,        // number of runs where this peptide was detected
}
```

### Target-Decoy Separation

Targets and decoys get independent consensus RTs. Decoy consensus RTs are computed from decoy detections only (not from target detections) to avoid information leakage. This ensures fair target-decoy competition in the second-pass FDR.

## Predicting Expected Measured RT

Once a consensus library RT exists, Osprey predicts where a peptide should elute in any run:

```text
expected_measured_rt = refined_calibration.predict(consensus_library_rt)
```

This is the RT used for reconciliation planning — the system asks "is the current peak's apex within `rt_tolerance` of this expected RT?" If not, it looks for an alternate CWT candidate whose apex is close to the expected RT, or imputes boundaries. See [Boundary Overrides & Forced Integration](13-boundary-overrides.md).

## RT Tolerance from Global MAD

The reconciliation apex-proximity check uses a **global per-file tolerance** derived from the refined calibration's residuals:

```text
rt_tolerance = max(0.1, 3.0 × MAD × 1.4826)
```

where `MAD` is the median absolute residual from the refined calibration's training points (thousands of consensus peptides). This is a single value per file, applied to all peptides in that file during reconciliation planning.

**Why global, not per-point?** The `RTCalibration::local_tolerance()` method interpolates absolute residuals at a query RT, giving a position-specific tolerance. This is appropriate for the initial search (where tolerance varies along the gradient). But for reconciliation it creates a self-fulfilling prophecy: a peptide with a wrong apex RT contributes a large residual to the training set, which inflates the local tolerance at that RT, which then allows the wrong-RT detection to pass the apex proximity check.

Using a global MAD eliminates this feedback loop. One bad peptide's residual barely moves the median over thousands of points, so the tolerance stays tight even if a few training points are outliers. The factor of 3.0 gives approximately 3-sigma coverage (99.7%) under a normal assumption, and the 0.1 min minimum floor prevents being overly strict when the calibration is exceptionally tight.

On a 3-replicate Astral HeLa dataset, the refined calibration MAD was ~0.09 min, giving `rt_tolerance = 3.0 × 0.09 × 1.4826 ≈ 0.40 min`. This flags ~0.3% of peptides (those with apex deviations beyond 0.4 min from the consensus) for re-scoring, which is consistent with the empirical distribution where the 99.7th percentile of apex deviations is near this value.

## Peak Width for Imputation

The `median_peak_width` serves as the integration window for forced integration when no CWT candidate exists at the expected RT:

```text
forced_start = expected_rt - median_peak_width / 2
forced_end   = expected_rt + median_peak_width / 2
```

Using the median peak width from confident detections ensures the imputed window matches the typical chromatographic peak width for this peptide.

---

## Implementation

| File | Function | Purpose |
|------|----------|---------|
| `crates/osprey-chromatography/src/calibration/rt.rs` | `RTCalibration::predict()` | Library RT → measured RT |
| `crates/osprey-chromatography/src/calibration/rt.rs` | `RTCalibration::inverse_predict()` | Measured RT → library RT |
| `crates/osprey/src/reconciliation.rs` | `compute_consensus_rts()` | Collect detections, compute weighted median library RTs |
| `crates/osprey/src/reconciliation.rs` | `refit_calibration_with_consensus()` | Refit per-run LOESS from consensus points |
| `crates/osprey/src/reconciliation.rs` | `weighted_median()` | Weighted median computation |

---

## References

- LOESS: Cleveland, W.S. (1979). "Robust locally weighted regression and smoothing scatterplots"
- RT alignment approach follows concepts from DIA-NN (Demichev et al. 2020) and Spectronaut (Bruderer et al. 2015)
- Weighted median for robust aggregation in the presence of outlier detections
