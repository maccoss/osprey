//! Runtime peptide trace for debugging reconciliation and peak selection.
//!
//! Enable by setting the environment variable `OSPREY_TRACE_PEPTIDE` to the
//! modified_sequence of interest (e.g. `DAQVVGMTTTGAAK`). Multiple peptides
//! can be traced by separating with commas. All matches are logged at
//! `log::info!` level with a `[trace]` prefix so they appear in both the
//! terminal and log file.

use std::sync::OnceLock;

static TRACE_SET: OnceLock<Option<Vec<String>>> = OnceLock::new();

fn trace_set() -> Option<&'static [String]> {
    TRACE_SET
        .get_or_init(|| {
            std::env::var("OSPREY_TRACE_PEPTIDE").ok().map(|v| {
                v.split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect::<Vec<_>>()
            })
        })
        .as_deref()
}

/// Returns true when tracing is enabled for `modified_sequence`.
///
/// Matches on suffix after the optional `DECOY_` prefix so that a trace
/// target set to the target sequence also catches its paired decoy.
pub fn is_traced(modified_sequence: &str) -> bool {
    let Some(targets) = trace_set() else {
        return false;
    };
    let stripped = modified_sequence
        .strip_prefix("DECOY_")
        .unwrap_or(modified_sequence);
    targets
        .iter()
        .any(|t| t == modified_sequence || t == stripped)
}

/// Log final q-values for all traced peptides after an FDR pass.
///
/// `stage` should be a short label like "first-pass" or "second-pass" to
/// distinguish the caller in the trace output.
pub fn log_fdr_qvalues(per_file_entries: &[(String, Vec<osprey_core::FdrEntry>)], stage: &str) {
    if trace_set().is_none() {
        return;
    }
    use osprey_core::FdrLevel;
    for (file_name, entries) in per_file_entries {
        for entry in entries {
            if !is_traced(&entry.modified_sequence) {
                continue;
            }
            let kind = if entry.is_decoy { "DECOY" } else { "TARGET" };
            log::info!(
                "[trace] fdr({}) {} z={} {} file={} apex={:.4} score={:.4} pep={:.4} run_q(prec/pep/eff)={:.4}/{:.4}/{:.4} exp_q(prec/pep/eff)={:.4}/{:.4}/{:.4}",
                stage,
                entry.modified_sequence,
                entry.charge,
                kind,
                file_name,
                entry.apex_rt,
                entry.score,
                entry.pep,
                entry.run_precursor_qvalue,
                entry.run_peptide_qvalue,
                entry.effective_run_qvalue(FdrLevel::Both),
                entry.experiment_precursor_qvalue,
                entry.experiment_peptide_qvalue,
                entry.effective_experiment_qvalue(FdrLevel::Both),
            );
        }
    }
}
