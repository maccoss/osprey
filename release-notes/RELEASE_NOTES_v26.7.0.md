# Osprey v26.7.0 Release Notes

Hardens FDR scoring against high-intensity interferences, decouples input
directories from derived output, and lands a large batch of cross-implementation
bit-parity and HPC-distribution correctness fixes validated against OspreySharp.

## Upgrade Note

This release changes how three scoring features are computed, so existing
per-file `.scores.parquet` caches from 26.6.x are treated as stale and are
transparently re-scored on the first 26.7.0 run. This is automatic: the parquet
footer records the Osprey version, and the cache validator re-scores any file
whose cached version does not match the running binary. Cross-implementation
`--join-only` runs abort (rather than silently mixing scales) if pointed at
parquets scored by a different major/minor version, so upgrade Rust and
OspreySharp together.

## New Features

- **`--work-dir` / `--output-dir` / `--cache-dir`: input data can stay
  read-only (#47).** Per-file derived artifacts (scores parquet, calibration
  JSON, FDR score sidecars, `reconciliation.json`, the `.spectra.bin` cache, and
  the library `.libcache`) were previously written next to the input mzML and
  library, forcing a copy of read-only data just to run an analysis. They now
  resolve through `resolve_output_dir` / `resolve_cache_dir`, so input data can
  stay in place while only derived output is written elsewhere. `--work-dir`
  sets both `--output-dir` and `--cache-dir`; explicit flags override. The
  `.spectra.bin` header gains a source size + mtime fingerprint (format
  version 2 to 3), byte-compatible with OspreySharp so the two tools can share
  one cache; the default path stays byte-identical.

## Scoring and FDR

- **Intensity magnitude features are log10-conditioned to prevent a Percolator
  hijack.** `peak_apex`, `peak_area`, and `peak_sharpness` reached the
  first-pass Percolator SVM as raw, linear, heavy-tailed values. Because a
  single experiment-wide standardizer and one linear model score every peak
  across all runs, a lone high-intensity DIA interference standardized to a
  z-score of 100-300 that dominated the linear discriminant, letting intensity
  outliers (a random mix of target / decoy / entrapment) hijack the top of the
  score ranking. On SEA-AD Astral-DIA entrapment runs this floored the
  achievable q-value, produced a top-of-ranking FDP spike, and at 1:1
  entrapment collapsed a whole run to zero IDs. All three features are now
  conditioned with `log10(max(0, x) + 1)` at their single computation point,
  matching Skyline mProphet's `MQuestIntensityCalc`. The transform is monotonic
  (within-feature ranking is preserved) and the `max(0, x)` floor keeps the
  feature finite even when `peak_sharpness` goes negative. `peak.area` itself
  stays raw for quantification; only the PIN feature is conditioned. Mirrors the
  paired ProteoWizard/pwiz change to preserve C#/Rust bit parity.

## Bug Fixes

- **Calibration pass-2 refit no longer discards a good pass-1 calibration
  (#52).** Pass 2 reused pass 1's minimum-points floor, so a file whose pass-1
  count fell in the `[ABSOLUTE_MIN_CALIBRATION_POINTS, min_calibration_points)`
  band could clear the absolute guard on a narrowed-window refit only to trip
  the calibrator's own `min_points` check, discarding a perfectly good pass-1
  calibration and running the file uncalibrated. The refit now builds its own
  adaptive floor, and graduated RT fit tiers step from full-bandwidth LOESS to a
  stiffer LOESS to a global linear fit as the point set thins.
- **MS1 envelope: non-ppm precursor tolerance now defaults to 10 ppm (#41).**
  A precursor tolerance specified in non-ppm units no longer degrades the MS1
  isotope-envelope match; it falls back to a 10 ppm default.
- **LDA side-effect re-sort no longer corrupts the MS2 calibration accumulator
  (#44).** An LDA re-sort mutated the ordering the MS2 calibration accumulator
  depended on. Added the `dump_ms2_cal_errors` diagnostic for verifying MS2
  calibration error distributions.
- **Reconciliation gap-fill excludes decoys (#45).** `identify_gap_fill_targets`
  seeded gap-fill with both target and decoy entry IDs; decoys already have a
  first-pass parquet row, so appending a gap-fill decoy duplicated it and skewed
  competition. Gap-fill now targets only real targets. Closes an Astral 3-file
  `group_qvalue` divergence.
- **Gap-fill isolation m/z filter is applied on the join and cached paths
  (#50).** The isolation-window m/z filter used during gap-fill scoring was not
  applied consistently when reading cached or joined data.
- **HPC chain Stage 7 correctness at `--join-at-pass=2` (#40).** Fixed a set of
  HPC-distribution-path bugs surfaced by a strict bit-parity test comparing the
  4-step HPC chain against a straight-through run: second-pass Percolator now
  runs at `--join-at-pass=2` when sidecars are missing, the reconciliation hash
  is `file_stems`-aware, and Stage 6 hard-error gates were tightened.
- **Cross-impl Stage 5/6/7 parity (#43).** Aligned Percolator sort order, dedup,
  `parquet_index` remap, and experiment q-value propagation with OspreySharp.
- **`reconciliation.json` carries the join-wide first-pass base_id set
  (v3, #48).** The boundary file now records `first_pass_base_ids` (format
  version 2 to 3), the join-wide set of base_ids that survived first-pass
  compaction, so a per-file HPC rescore worker can compact to exactly the set
  the in-memory pipeline used. Format stays byte-identical across
  implementations.
- **Calibration determinism (#39, #46).** Replaced the Welford running mean with
  a `sum / n` computation in `calculate_single_calibration` and applied
  ULP-scale tweaks so Rust and OspreySharp produce bit-identical calibrations.
- **Calibration pass 2 diagnostics (#42).** Refreshed the LOESS dump and added
  `num_confident_peptides` metadata.
- **Library-decoy mode: the "no decoys detected" hard error is deferred until
  after the manifest pass (#35).** The pipeline previously bailed out
  immediately if the prefix scan and `Decoy` column flagged zero entries,
  defeating the manifest-authoritative recovery path for predictor-stripped
  libraries (Carafe runs where the generator dropped the `decoy_`/`rev_` prefix
  on protein accessions). The check now runs after manifest application and
  composition pairing have had a chance to flip `is_decoy`, and the error
  message mentions the manifest as a third way to supply decoys. Behavior is
  unchanged when a library has any prefix-matched or column-flagged decoys.

## Diagnostics

- Env-var-gated dumps for cross-impl bisection (#38) and buffered diagnostic
  dump writers that flush before an early exit (#36), plus tests for the
  manifest decoy-recovery path.

## Breaking Changes

- **Relicensed from Apache-2.0 to LGPL-3.0-only.** The `LICENSE` file now
  contains the GNU LGPL v3 text and `LICENSE.GPL` ships the GPL-3.0 text it
  incorporates by reference. `NOTICE.md` is unchanged; its Apache-2.0 references
  describe upstream dependencies, which keep their original licensing. The
  README now points users to OspreySharp, the ProteoWizard port maintained by
  the Skyline development team, which gives byte-identical results with direct
  vendor RAW access; the Osprey Rust tool is maintained as the reference
  implementation for cross-impl bit parity.
