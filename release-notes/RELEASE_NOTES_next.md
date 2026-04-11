# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

<!-- none yet -->

## Bug Fixes

- Fixed missing ID lines in Skyline for precursors passing experiment-level FDR. After two-pass FDR with reconciliation, second-pass run q-values can shift slightly above the run FDR threshold (e.g., 0.0107 vs 0.01) even though the precursor passes experiment-level FDR. Previously, these precursors had NULL `retentionTime` in all runs, producing no ID line in Skyline despite being confidently identified at the experiment level. Now the best run (lowest run q-value) gets an ID line as a fallback, guaranteeing at least one ID line per precursor in the blib. Affected ~17% of precursors (4,940/28,636) in a 3-replicate Stellar dataset.

- Fixed reconciliation incorrectly keeping peaks at wrong RTs. The `determine_reconcile_action` function used boundary containment (`start_rt <= expected <= end_rt`) to decide whether a peak was at the correct RT. Peaks detected at the wrong RT (e.g., apex 1+ min off consensus) but with wide trailing edges that happened to span the expected RT were classified as "Keep" and never re-scored. Replaced with apex proximity logic using the refined per-file calibration's `local_tolerance`, which reflects actual run-to-run RT variation at each gradient position derived from thousands of consensus peptides. CWT candidate selection also changed from boundary overlap to apex proximity, picking the candidate with the closest apex. Added 2 regression tests.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->

