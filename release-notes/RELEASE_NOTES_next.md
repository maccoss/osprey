# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

- **`--fdrbench <path>` flag**: writes an [FDRBench](https://github.com/Noble-Lab/FDRBench)-compatible input TSV directly from the pre-compaction first-pass FDR pool. The output includes **every scored target** (not just FDR-passing ones), with the raw SVM discriminant as `score` and the q-value matching `--fdr-level` (Both emits precursor-level; Protein writes the picked-protein TSV from the first-pass parsimony result). Decoys are excluded per FDRBench convention; entrapment sequences (carrying `_p_target` on protein accessions) pass through unchanged. Add `--fdrbench-per-run` to emit one row per (precursor, run) using run-level q-values. Invoke FDRBench with `-score 'score:1'`. This supersedes the lossy blib-based path in `scripts/build_fdrbench_input.py` for Osprey results (the blib only contains FDR-passing entries and has no SVM-score column).

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
