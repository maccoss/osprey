# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

- **Protein-level FDR now runs by default at 0.01.** Protein parsimony and picked-protein FDR (Savitski 2015) are always on — no more `--protein-fdr` flag required. The default threshold of 0.01 controls reporting and the compaction/reconciliation rescue rule. Override with `--protein-fdr <threshold>` (e.g., `--protein-fdr 0.05` for a looser cutoff). The blib output is gated by precursor-level FDR by default; to restrict output to peptides from protein-FDR-passing groups, use `--fdr-level protein`. Config type changed from `Option<f64>` to `f64`; existing YAML configs with `protein_fdr: null` or omitting the field will pick up the new default.

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

- **`protein_fdr` config field is now a required `f64` (default 0.01), not `Option<f64>`.** Existing YAML configs with the field omitted or set to `null` will load with the new 0.01 default. CLI behavior is backward compatible: `--protein-fdr <value>` still overrides the default, and omitting the flag now means "use the default 0.01" instead of "disable protein FDR".
