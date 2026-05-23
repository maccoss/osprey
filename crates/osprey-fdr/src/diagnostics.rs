//! Cross-impl bisection dumps for protein-FDR aggregation.
//!
//! All functions in this module are env-var-gated via `is_dump_enabled`
//! and no-op in production runs. They write small TSVs to the current
//! working directory so they can be diffed against the OspreySharp port's
//! matching dumps in `OspreyDiagnostics.cs`. See
//! `ai/.tmp/diff_winners.py` and `ai/.tmp/diff_best_peptide_scores.py`
//! for the numerical comparison drivers.

use osprey_core::diagnostics::{format_f64_roundtrip, is_dump_enabled};
use std::collections::HashMap;
use std::io::Write;
use std::sync::Arc;

use crate::protein::PeptideScore;

/// Cheap predicate the caller can use to gate any pre-dump allocation
/// (e.g. extracting `&[f64]` / `&[bool]` slices from a private struct).
/// Mirrors the env var checked by [`dump_stage7_winners`] so the caller
/// and the dump function stay in lockstep.
pub fn stage7_winners_enabled() -> bool {
    is_dump_enabled("OSPREY_DUMP_STAGE7_WINNERS")
}

/// Dump the full cumulative-FDR winners list (target + decoy together).
///
/// The existing `dump_stage7_protein_fdr` (in `crates/osprey/src/diagnostics.rs`)
/// emits only target winners' scores; decoy-winner scores fall to NaN
/// because `ProteinFdrResult.group_scores` is target-only. Cross-impl
/// divergences driven by decoy-winner scores or sort-position interleaving
/// in the cumulative-FDR sweep are invisible without exposing them.
///
/// Gated by `OSPREY_DUMP_STAGE7_WINNERS=1`. Caller passes the per-winner
/// arrays in sort order so the dump's `rank` column matches the order
/// driving the sweep. Use [`stage7_winners_enabled`] to short-circuit
/// the slice construction in the disabled-dump (production) path.
pub fn dump_stage7_winners(
    winner_scores: &[f64],
    winner_is_decoys: &[bool],
    raw_qvalues: &[f64],
    monotonic_qvalues: &[f64],
) {
    if !is_dump_enabled("OSPREY_DUMP_STAGE7_WINNERS") {
        return;
    }
    let path = "rust_stage7_winners.tsv";
    let Ok(file) = std::fs::File::create(path) else {
        log::warn!("Could not create {}", path);
        return;
    };
    let mut f = std::io::BufWriter::with_capacity(8 << 20, file);
    let _ = writeln!(f, "rank\tscore\tis_decoy\traw_qvalue\tmonotonic_qvalue");
    for i in 0..winner_scores.len() {
        let _ = writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            i,
            format_f64_roundtrip(winner_scores[i]),
            winner_is_decoys[i],
            format_f64_roundtrip(raw_qvalues[i]),
            format_f64_roundtrip(monotonic_qvalues[i]),
        );
    }
    let _ = f.flush();
    log::info!("Wrote {} ({} winners)", path, winner_scores.len());
}

/// Dump the per-modseq aggregated best-score map from
/// `collect_best_peptide_scores`.
///
/// Surfaces the protein-FDR input set so divergences in upstream
/// aggregation (e.g. different per-peptide max scores cross-impl from
/// post-compaction entry-list asymmetry) can be diffed directly,
/// instead of being inferred from later-stage dumps.
///
/// Gated by `OSPREY_DUMP_BEST_PEPTIDE_SCORES=1`. Rows are sorted by
/// modified_sequence so the cross-impl diff is stable.
pub fn dump_best_peptide_scores(best: &HashMap<Arc<str>, PeptideScore>) {
    if !is_dump_enabled("OSPREY_DUMP_BEST_PEPTIDE_SCORES") {
        return;
    }
    let path = "rust_best_peptide_scores.tsv";
    let Ok(file) = std::fs::File::create(path) else {
        log::warn!("Could not create {}", path);
        return;
    };
    let mut f = std::io::BufWriter::with_capacity(8 << 20, file);
    let _ = writeln!(f, "modified_sequence\tscore\tis_decoy\tbest_qvalue");
    let mut rows: Vec<(&str, &PeptideScore)> = best.iter().map(|(k, v)| (k.as_ref(), v)).collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    for (seq, ps) in rows {
        let _ = writeln!(
            f,
            "{}\t{}\t{}\t{}",
            seq,
            format_f64_roundtrip(ps.score),
            ps.is_decoy,
            format_f64_roundtrip(ps.best_qvalue),
        );
    }
    let _ = f.flush();
    log::info!("Wrote {} ({} peptides)", path, best.len());
}
