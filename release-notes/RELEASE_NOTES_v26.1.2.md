# Osprey v26.1.2 Release Notes

Patch release fixing a critical reconciliation bug and several related issues. The reconciliation fix significantly improves RT consistency across replicates on large datasets, and increases second-pass precursor counts by preserving correct first-pass detections that were being damaged by unnecessary re-scoring. All users should upgrade and re-run affected analyses from scratch (delete `.scores.parquet` caches before re-running).

## Bug Fixes

- Fixed reconciliation incorrectly keeping peaks at wrong RTs. The `determine_reconcile_action` function used boundary containment (`start_rt <= expected <= end_rt`) to decide whether a peak was at the correct RT. Peaks detected at the wrong RT (e.g., apex 1+ min off consensus) but with wide trailing edges that happened to span the expected RT were classified as "Keep" and never re-scored. Replaced with apex proximity logic using the refined per-file calibration's `local_tolerance`, which reflects actual run-to-run RT variation at each gradient position derived from thousands of consensus peptides. CWT candidate selection also changed from boundary overlap to apex proximity, picking the candidate with the closest apex. Added 2 regression tests.

- Fixed missing ID lines in Skyline for precursors passing experiment-level FDR. After two-pass FDR with reconciliation, second-pass run q-values can shift slightly above the run FDR threshold (e.g., 0.0107 vs 0.01) even though the precursor passes experiment-level FDR. Previously, these precursors had NULL `retentionTime` in all runs, producing no ID line in Skyline despite being confidently identified at the experiment level. Now the best run (lowest run q-value) gets an ID line as a fallback, guaranteeing at least one ID line per precursor in the blib. Affected ~17% of precursors (4,940/28,636) in a 3-replicate Stellar dataset.

- Fixed log file not being fully written. The `BufWriter` wrapping the log file was not flushed on exit, causing the last ~8KB of log output to be lost. Now flushes after every log line and at program exit.

## New Features

### Two-tier logging

Terminal output is now clean and scannable, while the log file captures full detail:

- Terminal shows info-level only, formatted as `[M:SS] message` with no module paths, matching DIA-NN's style. Warn/Error get a level prefix.
- Log file always captures debug-level with full timestamps and module paths.
- `--verbose` flag makes the terminal show debug-level too (restores old behavior).

Implemented via a new `TwoTierLogger` in `crates/osprey/src/logging.rs` that replaces the previous `TeeWriter` + `env_logger` setup. Downgraded ~60 verbose `log::info!` calls to `log::debug!` across pipeline, percolator, calibration_ml, and calibration modules.

### Labeled progress bars

Progress bars now include descriptive labels so users can identify which step is running:

- `Scoring` for the initial coelution search
- `Re-scoring` for multi-charge consensus and inter-replicate reconciliation
- `Gap-fill CWT` for the CWT pass of gap-fill
- `Gap-fill forced` for the forced-integration pass
- `SVM training` for Percolator cross-validation training

Added a new progress bar for Percolator SVM training showing fold-iteration progress across the 3 parallel folds (up to 10 iterations each, stops early when no improvement).

### Clearer calibration and FDR log messages

- Calibration summary now shows mass shift alongside tolerances: `Calibration: 2758 peptides, RT 0.54 min, MS1 shift -0.92 ppm (tol 5.20), MS2 shift -1.03 ppm (tol 5.96)`
- LDA scoring line shows the S/N filter impact: `LDA scoring: 4865 at 1% FDR, 2758/4865 with S/N >= 5`
- First-pass FDR message clarifies the distinction from experiment-level: `First-pass: 147694 unique precursors pass run-level 1% FDR across all files (compacting for reconciliation)`
- Blib output uses "observations" instead of "passing plan entries" for clarity

## Upgrade Notes

Users running v26.1.1 should delete all `.scores.parquet` cache files before re-running with v26.1.2. The old caches were produced with the boundary-containment reconciliation logic, which left some peaks at incorrect RTs. Calibration JSON files and spectra caches (`.spectra.bin`) do not need to be deleted.

## Results Impact

On a 3-replicate Stellar HeLa dataset at 1% FDR:

| Metric | v26.1.1 | v26.1.2 | Change |
|--------|---------|---------|--------|
| Experiment precursors | 28,636 | 34,941 | +22% |
| Experiment peptides | 25,412 | 31,770 | +25% |
| Protein groups | 5,117 | 5,967 | +17% |

The increase comes from preserving correct first-pass detections. The old boundary-containment reconciliation was too aggressive, re-scoring entries whose apex was at the correct RT but whose boundaries didn't quite contain the expected RT. Those entries were re-scored at forced boundaries, producing worse features than the original, and they failed second-pass FDR. The apex-proximity logic correctly keeps these entries.
