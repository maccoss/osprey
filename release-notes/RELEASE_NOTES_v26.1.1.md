# Osprey v26.1.1 Release Notes

Patch release fixing critical `parquet_index` bugs in the blib output writer and reconciliation write-back that caused incorrect retention times in blib files and silently disabled multi-charge consensus and inter-replicate reconciliation. All users of v26.1.0 should upgrade and re-run affected analyses from scratch (delete `.scores.parquet` caches before re-running).

## Bug Fixes

- Fixed blib output writing wrong RT/boundaries per peptide after first-pass FDR compaction. The blib plan entry construction zipped compacted LightFdr data (shorter Vec after non-passing entries removed) against the full Parquet projection by position. After compaction, Vec position != Parquet row index, so each peptide received the apex_rt/start_rt/end_rt of a different peptide. This produced blib files where peptides appeared at wildly wrong RTs across replicates (e.g., 9-minute span for a peptide that should cluster within 0.6 min). Fixed by adding `parquet_index` to LightFdr and using it to index into the Parquet projection instead of zipping by position. Added 2 regression tests.

- Fixed reconciliation/consensus write-back overwriting wrong Parquet row after compaction. When multi-charge consensus or inter-replicate reconciliation re-scored an entry, the updated `CoelutionScoredEntry` was written back to `full_entries[idx]` where `idx` was the compacted FdrEntry Vec position, not the actual Parquet row. This corrupted an unrelated entry's data and left the target entry unchanged in the Parquet cache, causing multi-charge consensus and inter-replicate reconciliation to have no effect on Parquet-cached scores. The second-pass SVM trained on corrupted features, producing inflated precursor counts. Fixed by using `fdr_entries[idx].parquet_index` for the Parquet row lookup. Added 2 regression tests.

- Fixed gap-fill `parquet_index` update during reconciliation Parquet write-back. Gap-fill entries (added when a precursor passes FDR in some replicates but is missing from others) had their features appended to the end of the Parquet file but the corresponding FdrEntry stubs were left with `parquet_index = u32::MAX` (the "no Parquet row" sentinel). When the second-pass FDR ran, the Phase 2 feature loader filtered out these entries, producing fewer feature vectors than entries and causing a matrix shape mismatch panic. Now the write-back assigns sequential `parquet_index` values to gap-fill stubs matching the rows they occupy in the updated Parquet. Added 4 regression tests.

## Regression Tests

Added 8 regression tests to catch the `parquet_index` vs Vec-position bug class:

- `test_blib_plan_uses_parquet_index_not_vec_position`: verifies blib plan entry construction uses parquet_index for Parquet row lookup after compaction
- `test_blib_plan_with_gap_fill_entries`: verifies gap-fill entries in the blib plan use their updated parquet_index values
- `test_rescore_writeback_uses_parquet_index_not_vec_position`: verifies reconciliation write-back uses parquet_index to locate the correct Parquet row
- `test_rescore_writeback_multi_charge_consensus_scenario`: end-to-end test of the multi-charge consensus write-back with realistic entry layout
- `test_compaction_preserves_parquet_index`: verifies compaction keeps both targets and paired decoys with correct parquet_index
- `test_feature_lookup_uses_parquet_index_after_compaction`: documents the correct lookup pattern
- `test_gap_fill_sentinel_parquet_index`: verifies gap-fill entries use u32::MAX sentinel
- `test_reconcile_writeback_updates_gap_fill_parquet_index`: verifies gap-fill stubs get real parquet_index after write-back

## Upgrade Notes

Users running v26.1.0 must delete all `.scores.parquet` cache files and `.fdr_scores.bin` sidecar files before re-running with v26.1.1. The old caches contain corrupted reconciliation data from the write-back bug. Calibration JSON files and spectra caches (`.spectra.bin`) do not need to be deleted.
