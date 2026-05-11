# Osprey v26.5.3 Release Notes

CI-only patch release that supersedes v26.5.2's `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24` workaround with proper version pins to GitHub Actions that natively ship Node 24. No behavioral changes to the analysis pipeline.

## New Features

<!-- none yet -->

## Bug Fixes

- **CI: bumped v4 GitHub Actions to Node 24-native majors.** v26.5.2 silenced the "deprecated, may not work as expected" warning by setting `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24=true`, but the release CI run still emitted a milder "Node 20 is deprecated, actions are being forced to run on Node 24" reminder because the underlying `action.yml` files still declared `runs.using: node20`. The proper fix is version pins to majors that ship Node 24 natively: `actions/checkout` v4 → v6, `actions/cache` v4 → v5, `actions/upload-artifact` v4 → v7, `actions/download-artifact` v4 → v8. All four current majors use `runs.using: node24` so the reminder disappears entirely. The `FORCE_JAVASCRIPT_ACTIONS_TO_NODE24` env var added in v26.5.2 is now redundant and removed from both `.github/workflows/ci.yml` and `.github/workflows/release.yml`.

- **CI: trimmed redundant `brew install cmake openssl` on macOS runners.** `cmake` and `openssl@3` are pre-installed on the `macos-latest` runner image, so `brew install` was emitting "already installed" notices on every release build. Only `openblas` actually needs installing now.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
