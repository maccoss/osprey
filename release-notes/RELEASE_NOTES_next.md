# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

<!-- none yet -->

## Bug Fixes

- Fixed missing ID lines in Skyline for precursors passing experiment-level FDR. After two-pass FDR with reconciliation, second-pass run q-values can shift slightly above the run FDR threshold (e.g., 0.0107 vs 0.01) even though the precursor passes experiment-level FDR. Previously, these precursors had NULL `retentionTime` in all runs, producing no ID line in Skyline despite being confidently identified at the experiment level. Now the best run (lowest run q-value) gets an ID line as a fallback, guaranteeing at least one ID line per precursor in the blib. Affected ~17% of precursors (4,940/28,636) in a 3-replicate Stellar dataset.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->

