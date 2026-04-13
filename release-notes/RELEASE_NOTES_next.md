# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

- **Rewrote protein FDR as true picked-protein (Savitski 2015).** v26.1.2's composite log-likelihood scoring was length-biased and produced artificially low protein q-values. The new implementation follows Savitski et al. (PMC4563723) exactly: each protein is scored by its **single best peptide** (not a sum), target and decoy sides compete via **pairwise picking** (target wins its pair if `target_score >= decoy_score`, else decoy wins), and cumulative FDR is computed on the winner list (`q = cum_decoys / cum_targets`, backward sweep for monotonicity). Target winners only are written into output q-values; decoy winners are statistical machinery. Added 10 regression tests including `test_picked_protein_fdr_best_not_sum` (verifies max-not-sum aggregation) and `test_picked_protein_fdr_cumulative` (verifies cumulative FDR math).

- **Two-pass protein FDR architecture.** Picked-protein FDR now runs twice:
  - **First pass** (before compaction, when `--protein-fdr` is set): uses the full pre-compaction peptide pool so target + decoy sides are symmetric (no compaction survivorship bias). Writes `run_protein_qvalue` on FdrEntry stubs. Used as a gate for compaction and reconciliation.
  - **Second pass** (after second-pass peptide FDR, authoritative): uses the reconciled + second-pass-scored peptides. Writes `experiment_protein_qvalue`. Feeds `FdrLevel::Protein` filtering and the protein CSV report.

- **Protein-aware compaction.** The compaction filter now rescues borderline peptides from strong proteins. A peptide survives compaction if EITHER its peptide-level q-value ≤ `reconciliation_compaction_fdr` (default 0.05, loosened from the old 0.01 `run_fdr`) OR its first-pass protein group q-value ≤ `config.protein_fdr`. The loosened peptide gate alone broadens the pool fed into reconciliation and second-pass FDR; the protein rule is additive and only applies when `--protein-fdr` is set.

- **Reconciliation consensus protein rescue.** `compute_consensus_rts()` accepts an optional `protein_fdr_threshold`. A target peptide qualifies as a consensus anchor if it passes peptide-level FDR directly OR its first-pass protein group passes protein-level FDR. Peptides from strong proteins can anchor consensus RT computation even when their own peptide q-value is borderline.

- **New config field**: `reconciliation_compaction_fdr` (default 0.05). Peptide q-value threshold for compaction.

- **Changed peptide q-value gate in protein FDR**: from `2.0 × run_fdr` (DIA-NN convention) to `1.0 × run_fdr` (Savitski convention). Matches the peptide-level FDR threshold exactly.

- Added `Protein` variant to `FdrLevel` and changed the default from `Both` to `Peptide`. The `--fdr-level` flag now accepts `precursor`, `peptide` (default), `protein`, or `both`.

- Protein parsimony runs on every analysis regardless of whether `--protein-fdr` is set.

- `SharedPeptideMode::Razor` uses a proper iterative greedy set cover (path-independent, byte-identical across runs).

- Clarified FDR log output to distinguish precursor, peptide, and protein level counts.

## Bug Fixes

- Fixed `--protein-fdr` not actually filtering the blib output. Previously, setting `--protein-fdr 0.01` computed protein q-values and wrote a protein report CSV, but the blib output filter always used `experiment_precursor_qvalue`, ignoring the protein q-values. Now, setting `--fdr-level protein` (along with `--protein-fdr`) correctly gates the blib output by protein-level FDR.

- Fixed compaction survivorship bias corrupting protein FDR. v26.1.2 ran protein FDR on compacted stubs, which retained only decoys paired with passing targets. This gave picked-protein a one-sided comparison and produced nonsense q-values (~all target groups passing at 1% FDR). The new first-pass protein FDR runs on the full pre-compaction pool so decoys paired with non-passing targets are still available as competitors. This removes the v26.1.2 "known limitation" warning about cleanly-separated protein score distributions.

## Breaking Changes

- **Protein PEP removed from output.** The `protein_pep` column has been removed from the protein CSV report, and `ProteinFdrResult.group_pep` has been dropped. Savitski's picked-protein algorithm does not compute protein-level PEP. Peptide-level and precursor-level PEP are unaffected.

- **`ProteinFdrResult` struct simplified.** Removed `group_pep` field. Downstream code that reads `fdr_result.group_pep` will not compile. Use `fdr_result.group_qvalues` for the authoritative q-value; PEP is no longer provided at the protein level.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
