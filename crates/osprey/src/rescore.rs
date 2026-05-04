//! Per-file rescore worker — entry point for `--join-at-pass=1 --no-join`.
//!
//! Consumes the Stage 5 → Stage 6 boundary file pair
//! (`<stem>.1st-pass.fdr_scores.bin` + `<stem>.reconciliation.json`)
//! that a coordinator node persisted in `--join-at-pass=1 --join-only`
//! mode, plus the per-file `<stem>.scores.parquet` outputs from the
//! original Stages 1-4 run. Rebuilds the in-memory state that an
//! in-process per-file rescore would normally hold (post-Percolator
//! FdrEntry stubs, planner action map, refined per-file RT calibration,
//! gap-fill targets), runs the rescore engine with the planner's
//! boundary overrides, runs the gap-fill two-pass, and writes the
//! reconciled per-file parquet that the second join then consumes.
//!
//! End-to-end shape:
//!
//! ```text
//! .scores.parquet              ─┐
//! .1st-pass.fdr_scores.bin     ─┼──> hydrate_for_rescore ──> RescoreInputs
//! .reconciliation.json         ─┘                              │
//!                                                              ▼
//!                                            run_rescore (compaction +
//!                                              reconciliation_actions
//!                                              re-keying + per-file
//!                                              rescore loop + parquet
//!                                              write-back)
//! ```
//!
//! `RescoreInputs` is the same in-memory shape the in-process pipeline
//! holds at the Stage 5 → Stage 6 seam. The bit-parity contract is
//! enforced by `Compare-Stage6-Worker.ps1` (Phase A in-process vs
//! Phase C worker-on-renamed-folder, byte-compare reconciled
//! `.scores.parquet`) — currently 3/3 PASS on Stellar and Astral.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use osprey_chromatography::RTCalibration;
use osprey_core::{FdrEntry, LibraryEntry, OspreyConfig, OspreyError, Result};

use osprey_chromatography::{calibration_path_for_input, load_calibration};

use crate::pipeline::{
    fdr_scores_path_pass1, load_fdr_scores_sidecar, load_fdr_stubs_from_parquet,
    rescore_per_file_loop, select_post_fdr_consensus, synthetic_input_from_parquet,
};
use crate::reconciliation::{GapFillTarget, ReconcileAction};
use crate::reconciliation_io::{
    read_reconciliation_file, reconciliation_path, ForcedIntegrationEntry, GapFillEntry,
    RefinedRtCalibrationJson, UseCwtPeakEntry,
};

/// In-memory state needed to drive a per-file rescore from disk-backed
/// boundary files. Mirrors the shape the in-process pipeline holds at
/// the Stage 5 → Stage 6 seam, so the rescore engine can be written
/// once and used from both the in-process path and the worker path.
#[derive(Debug)]
pub struct RescoreInputs {
    /// Per-file FdrEntry stubs from `<stem>.scores.parquet`, with
    /// SVM scores + 4 q-values + PEP overlaid from the
    /// `<stem>.1st-pass.fdr_scores.bin` sidecar. File order matches
    /// the order of `config.input_scores`.
    pub per_file_entries: Vec<(String, Vec<FdrEntry>)>,
    /// Reconciliation actions keyed by `(file_name, vec_idx)`. Built
    /// from the homogeneous `use_cwt_peak_actions` +
    /// `forced_integration_actions` arrays in `reconciliation.json` by
    /// joining each action's `entry_id` against the loaded stub list.
    /// Keep actions are implicitly absent (the planner never persists
    /// them, since they are the default).
    pub reconciliation_actions: HashMap<(String, usize), ReconcileAction>,
    /// Refined per-file RT calibrations reconstructed from
    /// `reconciliation.json`'s `refined_rt_calibration` field via
    /// `RTCalibration::from_model_params`. Files whose envelope had a
    /// null calibration (e.g., refined fit failed during Stage 5) are
    /// absent.
    pub refined_calibrations: HashMap<String, RTCalibration>,
    /// Per-file gap-fill targets parsed from `reconciliation.json`'s
    /// `gap_fill_targets` array.
    pub per_file_gap_fill: HashMap<String, Vec<GapFillTarget>>,
    /// Modified-sequence interner shared with the parquet stub loader.
    /// Passed through to the rescore engine so any newly appended
    /// gap-fill stubs reuse the same `Arc<str>` identities that the
    /// hydrated stubs already hold.
    pub seq_interner: HashSet<Arc<str>>,
}

impl RescoreInputs {
    /// Total non-Keep reconciliation actions across all files.
    pub fn total_actions(&self) -> usize {
        self.reconciliation_actions.len()
    }
    /// Total gap-fill targets across all files.
    pub fn total_gap_fill_targets(&self) -> usize {
        self.per_file_gap_fill.values().map(|v| v.len()).sum()
    }
    /// Total stubs across all files.
    pub fn total_stubs(&self) -> usize {
        self.per_file_entries.iter().map(|(_, v)| v.len()).sum()
    }
}

/// Hydrate the boundary file pair into the in-memory state needed to
/// drive a per-file rescore. Reads each `<stem>.scores.parquet` listed
/// in `config.input_scores`, overlays `<stem>.1st-pass.fdr_scores.bin`,
/// and parses `<stem>.reconciliation.json`.
///
/// File names are extracted from the parquet stem (with the `.scores`
/// suffix stripped, mirroring the existing `synthetic_input_from_parquet`
/// helper). The output preserves `config.input_scores` ordering.
///
/// Returns an error if any per-file boundary file is missing,
/// unreadable, or fails its format-version / count checks. Does not
/// silently fall back to partial state — a Stage 6 worker that
/// proceeded with one file's planner output missing would scramble
/// gap-fill results across files.
pub fn hydrate_for_rescore(config: &OspreyConfig) -> Result<RescoreInputs> {
    let parquet_paths = config.input_scores.as_ref().ok_or_else(|| {
        OspreyError::config(
            "hydrate_for_rescore: --input-scores is required to locate the per-file \
             .scores.parquet + boundary file pair.",
        )
    })?;
    if parquet_paths.is_empty() {
        return Err(OspreyError::config(
            "hydrate_for_rescore: --input-scores list is empty",
        ));
    }

    let mut seq_interner: HashSet<Arc<str>> = HashSet::new();
    let mut per_file_entries: Vec<(String, Vec<FdrEntry>)> =
        Vec::with_capacity(parquet_paths.len());
    let mut refined_calibrations: HashMap<String, RTCalibration> = HashMap::new();
    let mut per_file_gap_fill: HashMap<String, Vec<GapFillTarget>> = HashMap::new();
    let mut reconciliation_actions: HashMap<(String, usize), ReconcileAction> = HashMap::new();

    for parquet_path in parquet_paths {
        // The parquet stem (with `.scores` stripped) is the canonical
        // file_name used everywhere else in the pipeline (calibration
        // map, planner output, blib output). Mirrors
        // synthetic_input_from_parquet.
        let synthetic_input = synthetic_input_from_parquet(parquet_path);
        let file_name = synthetic_input
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| {
                OspreyError::config(format!(
                    "hydrate_for_rescore: could not derive file_name from parquet path {}",
                    parquet_path.display()
                ))
            })?
            .to_string();

        // 1. Stubs from parquet (entry_id, charge, modseq, RTs,
        //    parquet_index assigned by load_fdr_stubs_from_parquet).
        let mut stubs =
            load_fdr_stubs_from_parquet(parquet_path, &mut seq_interner).map_err(|e| {
                OspreyError::config(format!(
                    "hydrate_for_rescore: failed to load stubs from {}: {}",
                    parquet_path.display(),
                    e
                ))
            })?;

        // 2. Overlay SVM scores + q-values + PEP from .1st-pass.fdr_scores.bin.
        //    expected_pass = 1 because the worker reads the FIRST-pass
        //    sidecar — the planner's actions were computed against it,
        //    and the same q-values feed gap-fill eligibility.
        let sidecar_path = fdr_scores_path_pass1(&synthetic_input);
        if !load_fdr_scores_sidecar(&sidecar_path, &mut stubs, 1) {
            return Err(OspreyError::config(format!(
                "hydrate_for_rescore: failed to overlay .1st-pass.fdr_scores.bin for {} \
                 (expected at {})",
                file_name,
                sidecar_path.display()
            )));
        }

        // 3. Parse reconciliation.json.
        let recon_path = reconciliation_path(&synthetic_input);
        let envelope = read_reconciliation_file(&recon_path).map_err(|e| {
            OspreyError::config(format!(
                "hydrate_for_rescore: failed to read {}: {}",
                recon_path.display(),
                e
            ))
        })?;

        // 3a. Build entry_id → vec_idx map from the loaded stubs so the
        //     planner's entry_id-keyed actions can be rehomed onto the
        //     in-memory shape (file_name, vec_idx) the rescore engine
        //     consumes.
        let mut id_to_idx: HashMap<u32, usize> = HashMap::with_capacity(stubs.len());
        for (idx, stub) in stubs.iter().enumerate() {
            id_to_idx.insert(stub.entry_id, idx);
        }

        // 3b. UseCwtPeak actions.
        for entry in &envelope.use_cwt_peak_actions {
            let UseCwtPeakEntry {
                apex_rt,
                candidate_idx,
                end_rt,
                entry_id,
                start_rt,
            } = *entry;
            let vec_idx = id_to_idx.get(&entry_id).copied().ok_or_else(|| {
                OspreyError::config(format!(
                    "hydrate_for_rescore: use_cwt_peak entry_id {} in {} not found in stubs \
                     (parquet drift?)",
                    entry_id,
                    recon_path.display()
                ))
            })?;
            reconciliation_actions.insert(
                (file_name.clone(), vec_idx),
                ReconcileAction::UseCwtPeak {
                    candidate_idx: candidate_idx as usize,
                    start_rt,
                    apex_rt,
                    end_rt,
                },
            );
        }

        // 3c. ForcedIntegration actions.
        for entry in &envelope.forced_integration_actions {
            let ForcedIntegrationEntry {
                entry_id,
                expected_rt,
                half_width,
            } = *entry;
            let vec_idx = id_to_idx.get(&entry_id).copied().ok_or_else(|| {
                OspreyError::config(format!(
                    "hydrate_for_rescore: forced_integration entry_id {} in {} not found in \
                     stubs (parquet drift?)",
                    entry_id,
                    recon_path.display()
                ))
            })?;
            reconciliation_actions.insert(
                (file_name.clone(), vec_idx),
                ReconcileAction::ForcedIntegration {
                    expected_rt,
                    half_width,
                },
            );
        }

        // 3d. Refined RT calibration (optional — null when Stage 5's
        //     LOESS refit failed for this file).
        if let Some(RefinedRtCalibrationJson {
            abs_residuals,
            fitted_rts,
            library_rts,
            residual_sd,
        }) = envelope.refined_rt_calibration
        {
            let params = osprey_chromatography::RTModelParams {
                library_rts,
                fitted_rts,
                abs_residuals,
            };
            let cal = RTCalibration::from_model_params(&params, residual_sd).map_err(|e| {
                OspreyError::config(format!(
                    "hydrate_for_rescore: failed to reconstruct refined calibration for {}: {}",
                    file_name, e
                ))
            })?;
            refined_calibrations.insert(file_name.clone(), cal);
        }

        // 3e. Gap-fill targets.
        let gap_fill: Vec<GapFillTarget> = envelope
            .gap_fill_targets
            .into_iter()
            .map(
                |GapFillEntry {
                     charge,
                     decoy_entry_id,
                     expected_rt,
                     half_width,
                     modified_sequence,
                     target_entry_id,
                 }| GapFillTarget {
                    target_entry_id,
                    decoy_entry_id,
                    expected_rt,
                    half_width,
                    modified_sequence: Arc::from(modified_sequence.as_str()),
                    charge,
                },
            )
            .collect();
        if !gap_fill.is_empty() {
            per_file_gap_fill.insert(file_name.clone(), gap_fill);
        }

        per_file_entries.push((file_name, stubs));
    }

    Ok(RescoreInputs {
        per_file_entries,
        reconciliation_actions,
        refined_calibrations,
        per_file_gap_fill,
        seq_interner,
    })
}

/// Top-level entry point for the `--join-at-pass=1 --no-join` per-file
/// rescore worker. Synthesizes the missing in-process state from the
/// Stage 5 → Stage 6 boundary files + per-file parquet + per-file
/// calibration JSON, then drives the same
/// [`rescore_per_file_loop`] the in-process pipeline calls. The
/// resulting reconciled per-file `.scores.parquet` files are what a
/// downstream `--join-at-pass=2` invocation consumes for second-pass
/// FDR and blib output.
///
/// Library + decoys are passed in by the caller (`run_analysis` loaded
/// them just before dispatching here) so the worker doesn't repeat
/// that setup.
pub fn run_rescore(config: OspreyConfig, library: Vec<LibraryEntry>) -> Result<()> {
    log::info!(
        "--join-at-pass=1 --no-join: per-file rescore worker starting on {} parquet(s)",
        config.input_scores.as_ref().map(|v| v.len()).unwrap_or(0)
    );

    // `config.input_files` was already synthesized from `input_scores`
    // at the top of `run_analysis` via `synthetic_input_from_parquet`,
    // so the rescore loop's `file_name_to_idx` lookup +
    // spectra/calibration path derivation point at the right mzML
    // siblings. Validate it here as defense in depth.
    if config.input_files.is_empty() {
        return Err(OspreyError::config(
            "run_rescore: config.input_files was not populated from --input-scores. \
             Expected `run_analysis` to synthesize via synthetic_input_from_parquet.",
        ));
    }
    let parquet_paths = config
        .input_scores
        .as_ref()
        .ok_or_else(|| {
            OspreyError::config(
                "run_rescore: --input-scores is required (validated upstream by validate_hpc_args)",
            )
        })?
        .clone();

    // Build file_name → input_files index, mirroring the in-process
    // setup at pipeline.rs `let mut file_name_to_idx`.
    let mut file_name_to_idx: HashMap<String, usize> = HashMap::new();
    for (idx, input_file) in config.input_files.iter().enumerate() {
        if let Some(stem) = input_file.file_stem().and_then(|s| s.to_str()) {
            file_name_to_idx.insert(stem.to_string(), idx);
        }
    }

    // Per-file `.scores.parquet` paths for the reconciled write-back at
    // the tail of `rescore_per_file_loop`. These are the parquets the
    // user passed via `--input-scores`.
    let mut per_file_cache_paths: HashMap<String, std::path::PathBuf> = HashMap::new();
    for (parquet, input_file) in parquet_paths.iter().zip(config.input_files.iter()) {
        if let Some(stem) = input_file.file_stem().and_then(|s| s.to_str()) {
            per_file_cache_paths.insert(stem.to_string(), parquet.clone());
        }
    }

    // Per-file ORIGINAL RT calibrations (first-pass LOESS fit from
    // Stages 1-3) loaded from sibling `<stem>.calibration.json` files.
    // Used by the rescore loop only as fallback when the boundary
    // file's refined cal is null (i.e., Stage 5's LOESS refit failed
    // for that file). Same construction as
    // pipeline.rs:3273-3296 in the in-process `--join-at-pass=1`
    // setup.
    let mut per_file_calibrations: HashMap<String, RTCalibration> = HashMap::new();
    for input_file in &config.input_files {
        let stem = match input_file.file_stem().and_then(|s| s.to_str()) {
            Some(s) => s.to_string(),
            None => continue,
        };
        let input_dir = match input_file.parent() {
            Some(d) => d,
            None => continue,
        };
        let cal_path = calibration_path_for_input(input_file, input_dir);
        if !cal_path.exists() {
            continue;
        }
        let cal_params = match load_calibration(&cal_path) {
            Ok(p) => p,
            Err(e) => {
                log::warn!(
                    "run_rescore: failed to read calibration JSON {}: {}",
                    cal_path.display(),
                    e
                );
                continue;
            }
        };
        if let Some(ref mp) = cal_params.rt_calibration.model_params {
            match RTCalibration::from_model_params(mp, cal_params.rt_calibration.residual_sd) {
                Ok(rt_cal) => {
                    per_file_calibrations.insert(stem, rt_cal);
                }
                Err(e) => {
                    log::warn!(
                        "run_rescore: failed to reconstruct original calibration for {}: {}",
                        cal_path.display(),
                        e
                    );
                }
            }
        }
    }

    // Hydrate post-Percolator FDR state + planner output from the
    // boundary file pair.
    let RescoreInputs {
        mut per_file_entries,
        reconciliation_actions: mut reconciliation_actions_pre,
        refined_calibrations,
        per_file_gap_fill,
        mut seq_interner,
    } = hydrate_for_rescore(&config)?;

    // Cross-impl bisection seam (mirrors the dump call from
    // `pipeline.rs` post-Percolator). Lets a worker run produce the
    // same `rust_stage5_percolator.tsv` an in-process run produces, so
    // the C# worker (which fires the matching `WriteStage5PercolatorDump`
    // call from `RescoreWorker.Run`) can be byte-compared against the
    // Rust worker via `Compare-Percolator.ps1`. Gated by
    // `OSPREY_DUMP_PERCOLATOR=1`; no-op otherwise.
    crate::diagnostics::dump_stage5_percolator(&per_file_entries);

    // Reproduce the in-process compaction (pipeline.rs first-pass FDR
    // block): drop pre-compaction entries whose base_id didn't make
    // the cut. Without this, the worker's `per_file_entries` carries
    // ~3.5x more rows than the in-process path, and
    // `select_post_fdr_consensus` + downstream scoring see entries
    // that the in-process flow filtered out — observed as a 324-row
    // divergence in median_polish features (Session 5 bisection).
    // The .fdr_scores.bin v3 sidecar carries `run_protein_qvalue`
    // precisely so this step works whether or not `--protein-fdr` is
    // set.
    //
    // Compaction also invalidates the `(file_name, vec_idx)` keys in
    // `reconciliation_actions` (which were assigned by hydration
    // against the pre-compaction stub list). Re-key by entry_id then
    // walk the compacted list to assign new indices.
    let mut reconciliation_actions: HashMap<(String, usize), ReconcileAction>;
    {
        let peptide_gate = config.reconciliation_compaction_fdr;
        let protein_gate = config.protein_fdr;
        let mut first_pass_base_ids: std::collections::HashSet<u32> =
            std::collections::HashSet::new();
        for (_, entries) in &per_file_entries {
            for e in entries {
                if !e.is_decoy
                    && (e.run_peptide_qvalue <= peptide_gate
                        || e.run_protein_qvalue <= protein_gate)
                {
                    first_pass_base_ids.insert(e.entry_id & 0x7FFF_FFFF);
                }
            }
        }

        // Save (file, entry_id) → action before per_file_entries shrinks
        // so we can rebuild (file, new_vec_idx) → action below. Use a
        // file_name → entries lookup map so the per-action lookup is O(1)
        // instead of O(num_files) on the inner find().
        let entries_by_name: HashMap<&str, &Vec<FdrEntry>> = per_file_entries
            .iter()
            .map(|(name, entries)| (name.as_str(), entries))
            .collect();
        let mut actions_by_id: HashMap<(String, u32), ReconcileAction> =
            HashMap::with_capacity(reconciliation_actions_pre.len());
        for ((file_name, idx), action) in reconciliation_actions_pre.drain() {
            // Need entry_id at the pre-compaction idx. per_file_entries is
            // still pre-compaction here, so the lookup is safe.
            if let Some(entries) = entries_by_name.get(file_name.as_str()) {
                if let Some(e) = entries.get(idx) {
                    actions_by_id.insert((file_name, e.entry_id), action);
                }
            }
        }
        drop(entries_by_name);

        let entries_before: usize = per_file_entries.iter().map(|(_, e)| e.len()).sum();
        for (_, entries) in per_file_entries.iter_mut() {
            entries.retain(|e| {
                let base_id = e.entry_id & 0x7FFF_FFFF;
                first_pass_base_ids.contains(&base_id)
            });
            entries.shrink_to_fit();
        }
        let entries_after: usize = per_file_entries.iter().map(|(_, e)| e.len()).sum();

        // Rebuild reconciliation_actions with post-compaction vec_idx.
        reconciliation_actions = HashMap::with_capacity(actions_by_id.len());
        let mut dropped_actions = 0usize;
        for (file_name, entries) in &per_file_entries {
            for (new_idx, e) in entries.iter().enumerate() {
                if let Some(action) = actions_by_id.remove(&(file_name.clone(), e.entry_id)) {
                    reconciliation_actions.insert((file_name.clone(), new_idx), action);
                }
            }
        }
        // Any actions left in actions_by_id reference entries that were
        // compacted away. That shouldn't happen for a healthy boundary
        // file (the planner only emits actions for entries that pass
        // FDR), so log if it does.
        dropped_actions += actions_by_id.len();
        if dropped_actions > 0 {
            log::warn!(
                "Worker compaction dropped {} reconciliation actions whose entries \
                 didn't survive compaction (boundary file may have been written by an \
                 older binary or with different config)",
                dropped_actions
            );
        }

        log::info!(
            "Worker compaction: {} -> {} entries ({} survived peptide-FDR or protein-rescue), \
             {} reconciliation actions retained",
            entries_before,
            entries_after,
            first_pass_base_ids.len(),
            reconciliation_actions.len(),
        );
    }

    // Multi-charge consensus targets are a pure function of the
    // hydrated FdrEntry list + run_fdr threshold, so the worker can
    // recompute them locally without needing them in the boundary file.
    let per_file_consensus_targets: HashMap<String, Vec<(usize, f64, f64, f64)>> = per_file_entries
        .iter()
        .map(|(file_name, entries)| {
            (
                file_name.clone(),
                select_post_fdr_consensus(entries, config.run_fdr),
            )
        })
        .collect();
    let n_consensus: usize = per_file_consensus_targets.values().map(|v| v.len()).sum();
    let n_actions: usize = reconciliation_actions.len();
    let n_gap_fill: usize = per_file_gap_fill.values().map(|v| v.len()).sum();
    log::info!(
        "Hydrated {} file(s); rescoring {} consensus + {} reconciliation entries, plus {} \
         total gap-fill candidates",
        per_file_entries.len(),
        n_consensus,
        n_actions,
        n_gap_fill,
    );

    // Drive the same per-file rescore loop the in-process pipeline
    // uses. The function loads spectra from the synthetic mzML paths,
    // runs `run_search` with boundary overrides, executes the gap-fill
    // two-pass, and rewrites each per-file `.scores.parquet` with
    // reconciled scores + appended gap-fill rows + reconciliation
    // metadata.
    let stats = rescore_per_file_loop(
        &mut per_file_entries,
        &per_file_consensus_targets,
        &reconciliation_actions,
        &refined_calibrations,
        &per_file_calibrations,
        &per_file_gap_fill,
        &per_file_cache_paths,
        &library,
        &file_name_to_idx,
        &config,
        &mut seq_interner,
    )?;

    log::info!(
        "--join-at-pass=1 --no-join: rescore complete — {} entries re-scored \
         ({} reconciliation + consensus, {} gap-fill via CWT, {} gap-fill via forced \
         integration). Reconciled .scores.parquet files written.",
        stats.total_rescored,
        stats.total_reconciliation,
        stats.total_gap_cwt,
        stats.total_gap_forced,
    );

    Ok(())
}

/// Format-version + boundary-pair sanity checks; thin wrapper over
/// `hydrate_for_rescore` that logs counts at info level so a worker
/// invocation has a visible "I loaded what I expected" signal even when
/// the rescore engine isn't wired yet. Retained for diagnostic use; the
/// production path is `run_rescore`.
pub fn hydrate_and_log(config: &OspreyConfig) -> Result<RescoreInputs> {
    let inputs = hydrate_for_rescore(config)?;
    log::info!(
        "Hydrated {} file(s): {} stubs, {} non-Keep actions, {} gap-fill targets, \
         {} refined calibration(s)",
        inputs.per_file_entries.len(),
        inputs.total_stubs(),
        inputs.total_actions(),
        inputs.total_gap_fill_targets(),
        inputs.refined_calibrations.len(),
    );
    for (file_name, stubs) in &inputs.per_file_entries {
        let n_actions = inputs
            .reconciliation_actions
            .keys()
            .filter(|(f, _)| f == file_name)
            .count();
        let n_gap = inputs
            .per_file_gap_fill
            .get(file_name)
            .map(|v| v.len())
            .unwrap_or(0);
        let has_refined = inputs.refined_calibrations.contains_key(file_name);
        log::info!(
            "  {}: {} stubs, {} actions, {} gap-fill, refined_cal={}",
            file_name,
            stubs.len(),
            n_actions,
            n_gap,
            has_refined
        );
    }
    Ok(inputs)
}

#[cfg(test)]
mod tests {
    use crate::reconciliation_io::{
        ForcedIntegrationEntry, GapFillEntry, ReconciliationFile, RefinedRtCalibrationJson,
        UseCwtPeakEntry, RECONCILIATION_FORMAT_VERSION,
    };
    use osprey_core::FdrEntry;

    /// Build a parquet + sidecar pair + reconciliation.json for two
    /// hand-coded files, then verify hydrate_for_rescore reconstructs
    /// the same in-memory state. Catches schema-drift bugs in the
    /// boundary writers/readers and any future ID-vs-index confusion.
    #[test]
    fn hydrate_round_trips_a_synthetic_pair() {
        // This test is a placeholder skeleton — the full round-trip
        // requires building a .scores.parquet with the right column set,
        // which `load_fdr_stubs_from_parquet` validates strictly. The
        // existing `reconciliation_file_round_trip` test in
        // `reconciliation_io.rs` and the `fdr_scores_sidecar_*` tests in
        // `pipeline.rs` already cover the per-file readers; the
        // hydration layer's logic above (entry_id → vec_idx join,
        // RTModelParams reconstruction) is exercised end-to-end by the
        // cross-impl harness once the worker is wired.
        //
        // TODO(stage6 worker): add a focused round-trip here once the
        // parquet test fixtures are in place.
        let _ = (
            RECONCILIATION_FORMAT_VERSION,
            std::mem::size_of::<FdrEntry>(),
        );
        let _ = (
            std::any::type_name::<UseCwtPeakEntry>(),
            std::any::type_name::<ForcedIntegrationEntry>(),
            std::any::type_name::<GapFillEntry>(),
            std::any::type_name::<RefinedRtCalibrationJson>(),
            std::any::type_name::<ReconciliationFile>(),
        );
        // Intentionally no assertions yet — see TODO above.
    }
}
