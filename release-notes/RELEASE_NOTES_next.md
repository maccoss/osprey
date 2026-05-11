# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

<!-- none yet -->

## Bug Fixes

- **`[unsorted-spectrum]` log demoted from `info` to `debug` (per-scan); aggregated end-of-load summary added.** The PR #33 centroid-sort fix that landed in v26.5.1 emits a log line for each MS2 spectrum whose centroids needed reordering — typically ~10–30 lines on a HeLa Astral 3 m/z DIA run (~0.07% of MS2 spectra), enough to clutter the default-verbosity log without adding actionable signal at that volume. The per-scan log is now `debug` level (visible only with `--verbose`), and both bulk mzML loaders (`load_all_spectra` and `load_ms1_spectra`) emit a single info-level summary at end-of-load: `"Sorted N spectra with non-monotonic centroids (run with --verbose to list scan_numbers)"`. The malformed-mzML length-mismatch case stays at `warn` since it's a data-integrity signal. `ensure_sorted` now returns `(mzs, intensities, did_sort)` so the bulk loaders can aggregate; the existing four unit tests grew explicit `did_sort` assertions.

- **CI: opted v4 GitHub Actions into Node 24 via `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24=true`.** Silences five Node 20 deprecation warnings on `actions/checkout@v4`, `actions/cache@v4`, `actions/upload-artifact@v4`, and `actions/download-artifact@v4` (visible on the CI run for v26.5.1). GitHub forces Node 24 for all JS actions on 2026-06-02; this opts in early without bumping each pinned v4 action's version tag. Applied to both `.github/workflows/ci.yml` and `.github/workflows/release.yml`. Pure CI hygiene — no behavioral change to the build itself. Comment in the workflows explains the transitional intent so a future maintainer knows when this env var can come out.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
