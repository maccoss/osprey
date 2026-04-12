# Osprey v26.1.3 Release Notes

Patch release fixing a self-fulfilling tolerance bug in reconciliation. Peaks detected at the wrong RT were inflating the local RT tolerance at that gradient position, which then allowed those same wrong-RT detections to pass the apex proximity check. Replaced per-point local tolerance with a global MAD-based tolerance robust to outliers.

## Bug Fixes

- Fixed reconciliation self-fulfilling tolerance. The apex proximity check in `plan_reconciliation` previously used `RTCalibration::local_tolerance(query_rt, 3.0, 0.1)`, which interpolates the absolute residual from the nearest training points at the query RT. This created a feedback loop: a peptide with a wrong apex RT contributed a large residual to the refined calibration's training set, which inflated the local tolerance at that RT, which then allowed the wrong-RT detection to pass the apex proximity check.

  Replaced with a **global per-file RT tolerance** derived from the refined calibration's MAD (median absolute deviation): `rt_tolerance = max(0.1, 3.0 × MAD × 1.4826)`. Since MAD is a median over thousands of training points, a single outlier's residual barely moves it. The tolerance is now computed once per file and applied to all peptides in that file during reconciliation planning.

  On a 3-replicate Astral HeLa dataset, this catches outliers like IPEAPAGPPSDFGLFLSDDDPK in file 49 (apex 19.78 min when consensus expects 18.60 min), which v26.1.2 incorrectly classified as "Keep" because the bad peptide's own residual had inflated the local tolerance to ~3.6 min. With the global MAD (~0.09 min), tolerance is ~0.4 min, and the 1.2 min deviation is correctly flagged for re-scoring.

  Added 4 regression tests:
  - `test_plan_reconciliation_tolerance_from_global_mad`
  - `test_plan_reconciliation_tolerance_minimum_floor`
  - `test_plan_reconciliation_tolerance_matches_expected_mad_formula`
  - `test_plan_reconciliation_outliers_in_training_do_not_inflate_tolerance`

## Documentation

- Updated `docs/10-cross-run-reconciliation.md` with the new apex-proximity-with-global-MAD logic and explanation of why per-point local tolerance was insufficient.
- Updated `docs/14-rt-alignment.md` with the RT Tolerance from Global MAD section, including empirical validation on the Astral dataset (0.3% of peptides flagged for re-scoring at 3.0 × MAD × 1.4826 = ~0.4 min tolerance).

## Upgrade Notes

Users running v26.1.2 should delete all `.scores.parquet` cache files before re-running with v26.1.3. The old caches may contain outlier peaks that the v26.1.2 reconciliation missed. Calibration JSON files and spectra caches (`.spectra.bin`) do not need to be deleted.
