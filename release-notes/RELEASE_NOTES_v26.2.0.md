# Osprey v26.2.0 Release Notes

Second feature release of 2026. Headline changes: true picked-protein FDR (Savitski 2015) with a two-pass architecture, protein-aware compaction and reconciliation, gas-phase-fractionation awareness in reconciliation and gap-fill, and a documented "picked-protein null collapses on Percolator-SVM output" known limitation with a path forward via decoy re-prediction.

## New Features

- **Rewrote protein FDR as true picked-protein (Savitski 2015).** v26.1.2's composite log-likelihood scoring was length-biased and produced artificially low protein q-values. The new implementation follows Savitski et al. (PMC4563723) exactly: each protein is scored by its **single best peptide SVM discriminant** (not a sum), target and decoy sides compete via **pairwise picking** (target wins its pair if `target_score >= decoy_score`, else decoy wins), and cumulative FDR is computed on the winner list (`q = cum_decoys / cum_targets`, backward sweep for monotonicity). Target winners only are written into output q-values; decoy winners are statistical machinery. The target-side **gate** still uses peptide q-value (`best_qvalue <= run_fdr`) per Savitski's convention, so only "reportable" target proteins enter the ranking; the ranking score itself is raw SVM discriminant. Added 10 regression tests including `test_picked_protein_fdr_best_not_sum` (verifies max-not-sum aggregation), `test_picked_protein_fdr_cumulative` (verifies cumulative FDR math), and `test_picked_protein_fdr_decoys_not_gated` (verifies asymmetric gating preserves the decoy null).

- **Protein FDR ranking score investigation — all three candidates produce a collapsed decoy null on this pipeline.** During development we tried three ranking scores and documented the outcome of each:
  - **Peptide q-value**: ~99% of decoy peptides have `q = 1.0` because they lost peptide-level TDC, so only the ~1% of decoys that won peptide TDC contribute a meaningful score. On Stellar 3-file HeLa this produced **2 decoy winners out of 6102** → cumulative FDR pinned at ~0.0003 and every target trivially passed 1% FDR. Structurally under-calibrated.
  - **Peptide PEP**: Osprey's `PepEstimator` is fit on TDC winners only, its bins cover `[min_winner_score, max_winner_score]`, and `posterior_error()` clamps out-of-range scores to `bins[0] ≈ 1.0`. Losing decoys all clamped to ~1.0. Same failure mode, **2 decoy winners out of 6102**.
  - **SVM discriminant (current)**: chosen because the earlier pre-fix run had given ~5.8% decoy winners (5640 target + 348 decoy), which looked well calibrated. On the post-fix rebuild of the same dataset it produces **2 decoy winners out of 6099** — same failure mode as the other two. The 348-decoy result was an **artifact of the q-value mapping bug** (see Bug Fixes): with `run_peptide_qvalue` pinned at 1.0, compaction rule (a) rejected everything, the compaction pool was tiny, second-pass Percolator trained on a smaller set, and the resulting SVM was less discriminating. Fixing the bug restored the correct compaction pool, which trained a sharper SVM, which eliminated the protein-level target/decoy overlap picked-protein depends on.
  - **Accepted state (option A)**: we kept SVM as the ranking score and documented that, on Osprey's Percolator output, picked-protein protein FDR provides essentially no additional control beyond peptide-level FDR. The 6097 groups passing at 1% on Stellar are exactly the target protein groups with at least one peptide passing peptide-level FDR. That is a defensible statement ("trust the peptide classifier and the protein FDR collapses to ~0") but it means the reported `protein_qvalue` column should be read as "protein has a qualifying peptide", not as an independently calibrated probability. The full decision history lives in the `compute_protein_fdr` doc comment in `crates/osprey-fdr/src/protein.rs`, in `docs/07-fdr-control.md` ("Why SVM Score (Not q-Value or PEP)?"), and in `CLAUDE.md` (protein FDR invariants).
  - **Known open issue — decoy asymmetry.** Target library entries use AI-predicted spectra and AI-predicted retention times (Carafe). Osprey's decoy generator produces reversed peptide sequences but **keeps the target's predicted RT and does not re-predict the decoy spectrum** — the decoy uses a reversed-fragment spectrum, not an AI-predicted spectrum for the reversed sequence. This asymmetry (targets: learned by Carafe; decoys: mechanical reversal) is a plausible reason the SVM can separate targets from decoys so sharply that picked-protein has no tail overlap. Re-predicting decoy spectra and RTs from Carafe for the reversed sequences would make the target and decoy pools comparable and may restore a meaningful protein-level null. Deferred to a future investigation; noted here so we can come back to it.

- **Two-pass protein FDR architecture.** Picked-protein FDR now runs twice:
  - **First pass** (before compaction, when `--protein-fdr` is set): uses the full pre-compaction peptide pool so target + decoy sides are symmetric (no compaction survivorship bias). Writes `run_protein_qvalue` on FdrEntry stubs. Used as a gate for compaction and reconciliation.
  - **Second pass** (after second-pass peptide FDR, authoritative): uses the reconciled + second-pass-scored peptides. Writes `experiment_protein_qvalue`. Feeds `FdrLevel::Protein` filtering and the protein CSV report.

- **Protein-aware compaction.** The compaction filter rescues borderline peptides from strong proteins. A peptide survives compaction if EITHER its peptide-level q-value ≤ `reconciliation_compaction_fdr` (default 0.01, matching `run_fdr`) OR its first-pass protein group q-value ≤ `config.protein_fdr`. The protein rule is additive and only applies when `--protein-fdr` is set. (An earlier draft loosened the peptide gate to 0.05 to broaden the reconciliation pool, but this inflated second-pass precursor counts because Percolator re-trained on an enriched set and produced sharper-than-warranted discrimination; reverted to 0.01.)

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

- **Fixed FDR result mapping dropping peptide-level q-values on FdrEntry stubs.** The Percolator and Mokapot result-mapping blocks in `pipeline.rs` were collapsing `max(precursor_q, peptide_q)` into `run_precursor_qvalue` and never writing `run_peptide_qvalue` or `experiment_peptide_qvalue`, leaving those fields at their init value of 1.0. This was a silent landmine after the 6-qvalue split in v26.x: anything downstream that read `run_peptide_qvalue` saw every entry as failing. In particular, the first-pass protein FDR built its detected-peptide set from `run_peptide_qvalue <= run_fdr`, which matched zero peptides, producing an empty parsimony graph, an empty compaction pool, an empty reconciliation pool, and an empty blib. Per-file FDR log lines ("N precursors, N peptides") still looked fine because they read from the local `PercolatorResult`, not the FdrEntry stub. Fix: write all four q-values (`run/experiment x precursor/peptide`) through separately in both Percolator and Mokapot mapping paths. The dual-FDR gate is applied downstream via `effective_*_qvalue(FdrLevel)`.

- **Filtered gap-fill by per-file isolation window m/z coverage (GPF fix).** Gap-fill previously generated a forced-integration target in every file for every precursor that passed run-level FDR in any replicate. That is wrong for gas-phase fractionation: each GPF file covers a disjoint precursor m/z range, so forcing an integration at a precursor outside the file's isolation windows lands on noise where that charge state physically cannot be observed. The same issue applies, more subtly, to any multi-replicate acquisition with partially disjoint isolation schemes. Gap-fill now takes a per-file map of isolation window intervals and a library `entry_id → precursor_mz` lookup; a candidate is skipped unless at least one of the target file's windows contains its precursor m/z. When windows are missing for a file the filter is disabled as a graceful fallback. Skipped candidates are counted and logged at INFO level. This also fixed a related GPF issue where nearly all inter-replicate reconciliation was failing silently because the gap-fill pool was fabricating bogus second-pass decoys from out-of-range precursors. Added two regression tests (`test_gap_fill_filters_precursors_outside_isolation_windows`, `test_gap_fill_allows_precursors_inside_isolation_windows`).

- **Multi-charge reconciliation is inherently cross-run for GPF, not per-file.** Clarified (in code and docs, no behavior change) that `select_post_fdr_consensus()` being a no-op in GPF is correct — per-file grouping has no multi-charge pairs to reconcile in GPF because each file holds at most one charge state per peptide. The cross-run consensus path (`compute_consensus_rts()` → `plan_reconciliation()`) is charge-agnostic: its keys are `(modified_sequence, is_decoy)` without charge, so a peptide's 2⁺ in one file and 3⁺ in another both contribute to a single shared consensus library RT, and each `(peptide, charge, file)` triple derives its expected measured RT from that consensus via the file's refined calibration. See `docs/10-cross-run-reconciliation.md` for the full write-up.

## Breaking Changes

- **Protein PEP removed from output.** The `protein_pep` column has been removed from the protein CSV report, and `ProteinFdrResult.group_pep` has been dropped. Savitski's picked-protein algorithm does not compute protein-level PEP. Peptide-level and precursor-level PEP are unaffected.

- **`ProteinFdrResult` struct simplified.** Removed `group_pep` field. Downstream code that reads `fdr_result.group_pep` will not compile. Use `fdr_result.group_qvalues` for the authoritative q-value; PEP is no longer provided at the protein level.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
