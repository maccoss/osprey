# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

- **`--fdrbench <path>` flag** (native, Rust): writes an [FDRBench](https://github.com/Noble-Lab/FDRBench)-compatible input TSV directly from the pre-compaction first-pass FDR pool. The output includes **every scored target** (not just FDR-passing ones), with the raw SVM discriminant as `score` and the q-value matching `--fdr-level` (Both emits precursor-level; Protein writes the picked-protein TSV from the first-pass parsimony result). Decoys are excluded per FDRBench convention; entrapment sequences (carrying `_p_target` on protein accessions) pass through unchanged. Add `--fdrbench-per-run` to emit one row per (precursor, run) using run-level q-values. Invoke FDRBench with `-score 'score:1'`. Verified on Astral entrapment data: 2.5M target rows emitted (20x the blib's FDR-passing count), 1.23M entrapment rows preserved, zero decoys, spearman(score, q_value) = -0.97 within the FDR-passing subset.

- **`scripts/build_fdrbench_input.py`**: companion utility for the cases where you don't have access to the source data Osprey scored from. Reads either a DIA-NN `report.parquet` or an Osprey `.blib` (+ companion `.proteins.csv`) and emits the same FDRBench TSV format at peptide/precursor/protein level. `score` is the raw upstream discriminant (DIA-NN `Evidence`, Osprey blib `DiscriminantScore` with `1 - q_value` fallback when the blib doesn't carry it, `best_peptide_score` for protein groups from the Osprey proteins CSV). The native `--fdrbench` flag is preferred for Osprey results (the blib only contains FDR-passing entries); the script is the right choice for DIA-NN parquets and for after-the-fact analysis of existing blibs.

- **Library-supplied decoys** (`--decoys-in-library`): Osprey now correctly trusts decoys already present in the spectral library instead of generating them from targets. Library entries whose protein accession matches one of the configured `decoy_prefixes` (default: `DECOY_`, `rev_`, `decoy_`, case-insensitive) get `is_decoy = true` and the high bit set on their `id`. Set via the new `--decoys-in-library` CLI flag or `decoys_in_library: true` in YAML; customize the prefix list via `decoy_prefixes` in YAML. Osprey now bails with a clear error if the flag is set but no entries match.

- **Target-decoy pairing for library decoys**: marking alone is not enough — the SVM, LDA calibration, and CV fold splitting all rely on each decoy sharing a `base_id` with its target. Osprey now runs one of two pairing strategies after marking. Preferred path is the new `--decoy-pairing-manifest <PATH>` flag, which reads a FDRBench-style 5-column TSV (`sequence`, `decoy`, `proteins`, `peptide_type`, `peptide_pair_index`) and links target↔decoy and `p_target`↔`p_decoy` within each `peptide_pair_index` group, honoring charge state. Without a manifest, Osprey falls back to composition-based pairing: for each decoy, strip the configured prefix from each protein accession to recover the target accession, then find a target peptide with the same `(stripped_accession, charge, sorted_amino_acid_composition)`. Within a composition group on a single protein, entries are sorted by `(sequence, id)` and zipped 1:1 deterministically. Osprey logs a pairing summary and **bails with a hard error if fewer than `decoy_pair_min_fraction` (default 80%) of decoys pair successfully** — FDR estimates would be unreliable otherwise. 14 new unit tests cover composition pairing edge cases (multi-peptide composition groups, shared peptides, charge-state separation, determinism) and manifest reader correctness.

## Bug Fixes

- **`decoys_in_library: true` and `decoy_method: FromLibrary` were silently broken.** The flag skipped Osprey's `DecoyGenerator` (as intended), but no code path detected decoys already present in the library — all four library loaders (DIA-NN TSV, blib, elib, libcache) hardcoded `is_decoy: false`. The `DecoyMethod::FromLibrary` enum variant additionally fell through to `Reverse` generation at the call site. Result: target-decoy competition silently treated every entry as a target and reported FDR estimates were meaningless. Fixed by the library-decoy support above; the two settings are now equivalent and both produce correctly-flagged decoys.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
