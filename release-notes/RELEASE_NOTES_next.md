# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

<!-- none yet -->

## Bug Fixes

- Fixed gap-fill `parquet_index` update during reconciliation Parquet write-back. Gap-fill entries (added when a precursor passes FDR in some replicates but is missing from others) had their features appended to the end of the Parquet file but the corresponding FdrEntry stubs were left with `parquet_index = u32::MAX` (the "no Parquet row" sentinel). When the second-pass FDR ran, the Phase 2 feature loader filtered out these entries, producing fewer feature vectors than entries and causing a matrix shape mismatch panic. Now the write-back assigns sequential `parquet_index` values to gap-fill stubs matching the rows they occupy in the updated Parquet. Added 4 regression tests to catch this bug class.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->

