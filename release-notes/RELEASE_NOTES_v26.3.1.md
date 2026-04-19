# Osprey v26.3.1 Release Notes

Post-v26.3.0 refinements to peak selection and RT calibration. Headline changes: classical Cleveland 1979 robust LOESS is now default, RT penalty sigma widened from 3x to 5x for better handling of peptides with slight calibration residuals, intensity tiebreaker added to CWT peak selection, plus cross-implementation diagnostics and regression tests from Brendan MacLean.

## New Features

- **Widened RT penalty sigma from 3x to 5x calibration MAD.** The 3-sigma Gaussian RT penalty introduced in v26.3.0 was too aggressive for peptides whose true elution time deviates slightly from the LOESS prediction. A peptide with its correct peak 0.3 min from the predicted RT could lose to a shoulder at the predicted position because the 10% RT penalty (at 3-sigma) was enough to flip the selection when coelution scores were comparable. At 5-sigma, the penalty at 0.3 min is ~4% (essentially negligible), while the penalty at 1.0 min is still ~37% (substantial for genuine wrong-peak selections). This preserves the interferer-rejection benefit while avoiding over-penalization of peaks with slight RT deviations from the calibration.

- **Added intensity tiebreaker to CWT peak selection.** When two CWT candidates from the same chromatographic peak (main peak vs shoulder) have nearly identical coelution scores, the more intense peak now wins. The selection score is multiplied by `log(1 + apex_intensity)`, which keeps intensity as a secondary factor that breaks ties without dominating the coelution-based ranking. This fixes cases where a narrow shoulder near the calibration-predicted RT beat the main peak because mean pairwise correlation is similar for both (they share the same eluting species) and the RT penalty gave the shoulder a slight edge.

- **Switched default LOESS to classical Cleveland 1979 robust iterations** (thanks to Brendan MacLean, [PR #10](https://github.com/maccoss/osprey/pull/10)). The legacy RT calibration reused residuals from the initial LOESS fit across every robustness iteration, giving only a single refinement regardless of `robustness_iter`. The classical Cleveland (1979) algorithm refreshes residuals from the current fit at the top of each iteration, producing tighter calibration curves on datasets with outliers. This is now the default; set `OSPREY_LOESS_CLASSICAL_ROBUST=0` to revert to the legacy single-refresh behavior for comparison testing.

- **Added cross-implementation bisection diagnostics** (thanks to Brendan MacLean, [PR #9](https://github.com/maccoss/osprey/pull/9)). New diagnostic modules in `osprey-core`, `osprey-scoring`, and `osprey` dump intermediate state at key pipeline stages (LOESS input points, XCorr per scan, median polish decomposition, CWT candidate scoring) to text files for byte-exact comparison against a parallel C# implementation (OspreySharp). Controlled by environment variables so there is zero overhead when not actively bisecting. Includes early-exit hooks (`should_exit_after_calibration`, `should_exit_after_scoring`) to stop the pipeline after a target stage without running the full analysis.

- **Added cross-implementation regression tests** (thanks to Brendan MacLean, [PR #11](https://github.com/maccoss/osprey/pull/11)). Five new tests in `osprey-scoring/src/lib.rs` cover XCorr and median polish scoring on fixed inputs with known reference outputs from the C# implementation. These guard against silent regressions when porting scoring changes between the two implementations.

- **Added reconciliation tolerance regression tests.** Three tests in `reconciliation.rs` guard against the reconciliation failures discovered on the Seer/Floyd and Stellar datasets: (1) `test_plan_reconciliation_tolerance_capped_by_original_calibration` verifies a contaminated refined calibration is capped at the original calibration's tolerance, (2) `test_plan_reconciliation_wrong_peak_1min_off_is_caught` verifies that a peak 1.0 min from the consensus-predicted RT is flagged for re-scoring, (3) `test_plan_reconciliation_correct_peak_near_expected_is_kept` verifies that a peak 0.05 min from expected is classified as Keep (not over-aggressively re-scored).

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
