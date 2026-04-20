# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

- **Protein-level FDR now runs by default at 0.01.** Protein parsimony and picked-protein FDR (Savitski 2015) are always on — no more `--protein-fdr` flag required. The default threshold of 0.01 controls reporting and the compaction/reconciliation rescue rule. Override with `--protein-fdr <threshold>` (e.g., `--protein-fdr 0.05` for a looser cutoff). The blib output is gated by precursor-level FDR by default; to restrict output to peptides from protein-FDR-passing groups, use `--fdr-level protein`. Config type changed from `Option<f64>` to `f64`; existing YAML configs with `protein_fdr: null` or omitting the field will pick up the new default.

- **Reconciliation tolerance now derived from within-peptide RT reproducibility, not cross-peptide calibration residuals.** The RT tolerance used by `plan_reconciliation` (to decide Keep / UseCwtPeak / ForcedIntegration) is now driven by the **global median of per-peptide apex-RT MADs in library space**, rather than the per-file LOESS calibration MAD. Once cross-run RT alignment has been applied, within-peptide scatter reflects the LC/instrument reproducibility floor and is typically 3-5x tighter than the cross-peptide calibration MAD. For well-aligned Stellar data this shrinks tolerance from ~0.3-0.6 min to ~0.03-0.1 min, so reconciliation now catches wrong-peak selections (e.g. 0.2 min off the correct peak) that previously squeaked through as "Keep". The per-file calibration MAD is retained as a safety ceiling in case a pathological global MAD arises. A log line at the start of reconciliation reports the global MAD value and the number of peptides (with ≥ 3 detections) that contributed.

- **Consensus RT weighting now uses sigmoid(SVM score) instead of coelution_sum.** The weighted median that computes each peptide's `consensus_library_rt` now weights each detection by `1 / (1 + exp(-score))`, where `score` is the SVM discriminant from first-pass FDR. A wrong-peak detection with negative SVM score gets near-zero weight (e.g., score=-4 → weight ~0.018) and cannot poison the median, while a strong detection dominates (score=+1.5 → weight ~0.818). Previously `coelution_sum` was used as the weight, which does not discriminate wrong-peak-but-high-coelution detections from true peaks. The `coelution_sum > 0` sanity filter is retained to reject anti-correlated noise integrations.

- **Protein-FDR rescue gate tightened: per-entry precursor q-value is now a hard precondition.** `compute_consensus_rts` previously allowed protein-FDR rescue to override any per-entry q-value, so a peptide with a wrong-peak detection (poor precursor q-value) could still be pulled into the consensus if its protein passed. Now the rescue can only upgrade borderline peptide-level evidence — the entry's own precursor q-value must still be ≤ `consensus_fdr`. This prevents a strong protein's weak or wrong-peak charge-state detections from poisoning the charge-agnostic consensus RT. Observed effect: on Stellar, target SVM scores shift higher because reconciliation no longer drags correct-peak observations to incorrect consensus positions, improving target/decoy separation in second-pass FDR.

- **Peptide trace diagnostic via `OSPREY_TRACE_PEPTIDE` environment variable.** Set `OSPREY_TRACE_PEPTIDE=<modified_sequence>` (comma-separated for multiple) to emit detailed per-peptide trace lines at `log::info!` level at five pipeline stages: first-pass CWT candidate scoring, `compute_consensus_rts`, `plan_reconciliation` (including per-entry tolerance derivation and stored CWT candidates), gap-fill target identification, and first/second-pass FDR q-values. Matches both target and `DECOY_`-prefixed paired decoy. Zero overhead when the env var is unset (`OnceLock` guarded). See [docs/17-peptide-trace.md](docs/17-peptide-trace.md) for the full trace format and typical diagnostic workflows.

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

- **`protein_fdr` config field is now a required `f64` (default 0.01), not `Option<f64>`.** Existing YAML configs with the field omitted or set to `null` will load with the new 0.01 default. CLI behavior is backward compatible: `--protein-fdr <value>` still overrides the default, and omitting the flag now means "use the default 0.01" instead of "disable protein FDR".
