# Osprey v26.6.1 Release Notes

Bug-fix patch: cross-run reconciliation now pairs library-supplied decoys to their targets via `base_id` instead of stripping a `DECOY_` prefix from `modified_sequence`. Without this fix, Carafe-predicted / FDRBench-manifest-paired decoys were silently excluded from reconciliation (`0 paired decoy peptides included`), leaving second-pass FDR trained on asymmetric pools.

## New Features

<!-- none yet -->

## Bug Fixes

- **Reconciliation now pairs library-supplied decoys to their targets by `base_id`, not by stripping a `DECOY_` prefix.** Both `compute_consensus_rts` and `plan_reconciliation` previously identified paired decoys by stripping a `DECOY_` prefix from the decoy's `modified_sequence` and matching the stripped form against the target set. That works for Osprey-generated decoys (modseq `DECOY_PEPTIDEK` → `PEPTIDEK`) but silently misses library-supplied decoys (Carafe / FDRBench-manifest-paired) whose modseq is the raw scrambled sequence with no prefix. Symptom in the log: `Inter-replicate reconciliation: 0 paired decoy peptides included` after a non-zero target count. With the bug, target boundaries got reconciliation-improved while decoy boundaries stayed at first-pass; the second-pass SVM then trained on asymmetric pools (targets enjoying the lift, decoys not), which could bias second-pass FDR optimistic. Empirically, FDRBench peptide-level FDR on Carafe-predicted decoys remained well-calibrated, suggesting the score lift from re-scoring is modest in aggregate, but the asymmetry was a correctness defect. After this fix, paired decoys are recognised via `base_id = entry_id & 0x7FFFFFFF` (the linkage Osprey establishes at library load via the FDRBench manifest or composition fallback), so reconciliation includes them regardless of decoy origin. Two new regression tests cover the case (`test_consensus_rts_pairs_library_decoy_by_base_id`, `test_plan_reconciliation_includes_library_decoy_via_base_id`).

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
