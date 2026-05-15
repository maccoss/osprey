# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

<!-- none yet -->

## Bug Fixes

- **Library-decoy mode: the "no decoys detected" hard error is now
  deferred until after the manifest pass.** Previously the pipeline
  bailed out immediately after `apply_library_decoy_marking` if the
  prefix scan and `Decoy` column together flagged zero entries.
  That defeated the manifest-authoritative recovery path from
  d23d496, where the manifest's `peptide_type=decoy` rows are the
  source of truth for predictor-stripped libraries (Carafe runs
  where the library generator dropped the `decoy_`/`rev_` prefix on
  protein accessions). The check now runs after manifest application
  and composition pairing have had a chance to flip `is_decoy` on
  entries the prefix scan missed; the error message has been
  updated to mention the manifest path as a third option for
  supplying decoys. Behavior is unchanged when a library has any
  prefix-matched or column-flagged decoys.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
