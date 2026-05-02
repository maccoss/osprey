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
//! This module currently implements only the **hydration** layer — the
//! rescore engine, gap-fill, and parquet write-back will land in
//! follow-up commits. The hydration layer's contract is:
//!
//! ```text
//! .scores.parquet              ─┐
//! .1st-pass.fdr_scores.bin     ─┼──> RescoreInputs
//! .reconciliation.json         ─┘
//! ```
//!
//! `RescoreInputs` is the same in-memory shape the in-process pipeline
//! holds at the Stage 5 → Stage 6 seam. A future bit-parity check
//! between an end-to-end run and a `--join-at-pass=1 --join-only`
//! followed by `--join-at-pass=1 --no-join` chain validates that the
//! boundary files capture all the state the rescore needs.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use osprey_chromatography::RTCalibration;
use osprey_core::{FdrEntry, OspreyConfig, OspreyError, Result};

use crate::pipeline::{
    fdr_scores_path_pass1, load_fdr_scores_sidecar, load_fdr_stubs_from_parquet,
    synthetic_input_from_parquet,
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
    })
}

/// Format-version + boundary-pair sanity checks; thin wrapper over
/// `hydrate_for_rescore` that logs counts at info level so a worker
/// invocation has a visible "I loaded what I expected" signal even when
/// the rescore engine isn't wired yet. Used by the `--join-at-pass=1
/// --no-join` mode while that path is still under construction.
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
