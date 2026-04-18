//! Cross-implementation bisection diagnostic dumps for the pipeline stage.
//!
//! Each function below is gated by an `OSPREY_DUMP_*` or `OSPREY_DIAG_*`
//! env var and writes a `rust_*.txt` file matching the format of the
//! corresponding C# OspreySharp dump (`cs_*.txt`). Dumps that have an
//! `_ONLY` companion env var exit the process after writing for fast
//! cycle-time during bisection.
//!
//! See `osprey-core::diagnostics` for shared primitives and
//! [`ai/scripts/OspreySharp/DIAGNOSTICS.md`] in the workspace for the full
//! env-var reference.

use osprey_core::diagnostics::{exit_if_only, is_dump_enabled};
use osprey_core::{LibraryEntry, Spectrum, XICPeakBounds};
use osprey_scoring::{SpectralScorer, TukeyMedianPolishResult};
use std::io::Write;

/// Dump the (lib_rt, measured_rt) pairs fed to LOESS to
/// `rust_loess_input.txt`, sorted by lib_rt ascending at `{:.17}`
/// full-bit precision so cross-impl diff catches any input drift before
/// LOESS fitting noise enters.
///
/// Gated by `OSPREY_DUMP_LOESS_INPUT=1`. When `OSPREY_LOESS_INPUT_ONLY=1`
/// is also set, exits the process after writing.
pub fn dump_loess_input(library_rts_detected: &[f64], measured_rts_detected: &[f64]) {
    if !is_dump_enabled("OSPREY_DUMP_LOESS_INPUT") {
        return;
    }

    if let Ok(mut f) = std::fs::File::create("rust_loess_input.txt") {
        writeln!(f, "idx\tlib_rt\tmeasured_rt").ok();
        let mut pairs: Vec<(f64, f64)> = library_rts_detected
            .iter()
            .zip(measured_rts_detected.iter())
            .map(|(&x, &y)| (x, y))
            .collect();
        pairs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        for (i, (lib_rt, meas_rt)) in pairs.iter().enumerate() {
            writeln!(f, "{}\t{:.17}\t{:.17}", i, lib_rt, meas_rt).ok();
        }
        log::info!(
            "Wrote LOESS input dump: rust_loess_input.txt ({} pairs)",
            pairs.len()
        );
    }

    exit_if_only("OSPREY_LOESS_INPUT_ONLY", "LOESS input dump");
}

/// Append a per-scan XCorr diagnostic block to `rust_xcorr_diag.txt`.
///
/// Gated by `OSPREY_DIAG_XCORR_SCAN=<scan_number>` matching
/// `apex_spectrum.scan_number`. Writes the preprocessed XCorr vector's
/// first 20 nonzero bins, sum, and per-fragment bin lookups so a
/// divergence in the slow-vs-fast XCorr code path can be localized to a
/// specific (scan, fragment, bin) tuple.
pub fn dump_xcorr_scan(
    apex_spectrum: &Spectrum,
    entry: &LibraryEntry,
    pre_vec: &[f32],
    lib_preprocessed: &[f32],
    xcorr_scaled: f64,
    scorer: &SpectralScorer,
) {
    let Ok(diag_scan) = std::env::var("OSPREY_DIAG_XCORR_SCAN") else {
        return;
    };
    if format!("{}", apex_spectrum.scan_number) != diag_scan {
        return;
    }
    let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open("rust_xcorr_diag.txt")
    else {
        return;
    };
    writeln!(
        f,
        "# XCORR DIAG scan={} entry={}",
        apex_spectrum.scan_number, entry.modified_sequence
    )
    .ok();
    writeln!(f, "# nbins={} xcorr_scaled={}", pre_vec.len(), xcorr_scaled).ok();
    let psum: f64 = pre_vec.iter().map(|&v| v as f64).sum();
    let pnz = pre_vec.iter().filter(|&&v| v != 0.0).count();
    writeln!(f, "# preprocessed_sum={} nonzero={}", psum, pnz).ok();
    // First 20 nonzero preprocessed bins
    let mut dumped = 0;
    for (i, &v) in pre_vec.iter().enumerate() {
        if v != 0.0 {
            writeln!(f, "pre\t{}\t{}", i, v).ok();
            dumped += 1;
        }
        if dumped >= 20 {
            break;
        }
    }
    // Fragment bin lookups + library preprocessed values
    writeln!(f, "# fragment_bins").ok();
    for (fi, frag) in entry.fragments.iter().enumerate() {
        let bin = scorer.bin_config().mz_to_bin(frag.mz);
        let pre_val = bin
            .map(|b| pre_vec.get(b).copied().unwrap_or(0.0))
            .unwrap_or(0.0);
        let lib_val = bin
            .map(|b| lib_preprocessed.get(b).copied().unwrap_or(0.0))
            .unwrap_or(0.0);
        writeln!(
            f,
            "frag\t{}\tmz={}\tbin={:?}\tpre_val={}\tlib_val={}",
            fi, frag.mz, bin, pre_val, lib_val
        )
        .ok();
    }
}

/// Write a median polish diagnostic block to `rust_mp_diag.txt`.
///
/// Gated by `OSPREY_DIAG_MP_SCAN=<scan_number>` matching
/// `apex_spectrum.scan_number`. Currently scoped to the
/// `DECOY_ALQFAQWWK` peptide (target hardcoded at the call site for a
/// specific historical bisection); the env var only triggers the dump
/// when both the scan and modified-sequence match.
///
/// Writes the row/col effects, grand mean, convergence info, and the
/// input matrix that the polish was run on.
#[allow(clippy::too_many_arguments)]
pub fn dump_mp_diag(
    apex_spectrum: &Spectrum,
    entry: &LibraryEntry,
    mp: &TukeyMedianPolishResult,
    peak: &XICPeakBounds,
    cos: f64,
    res: f64,
    min_r2: f64,
    resid_corr: f64,
    peak_xics: &[(usize, Vec<(f64, f64)>)],
) {
    let Ok(diag_scan) = std::env::var("OSPREY_DIAG_MP_SCAN") else {
        return;
    };
    let scan_str = format!("{}", apex_spectrum.scan_number);
    if scan_str != diag_scan || !entry.modified_sequence.contains("DECOY_ALQFAQWWK") {
        return;
    }
    let Ok(mut f) = std::fs::File::create("rust_mp_diag.txt") else {
        return;
    };
    writeln!(
        f,
        "# Median polish diagnostic for {} scan={}",
        entry.modified_sequence, apex_spectrum.scan_number
    )
    .ok();
    writeln!(
        f,
        "# peak range: start={} apex={} end={} len={}",
        peak.start_index,
        peak.apex_index,
        peak.end_index,
        peak.end_index - peak.start_index + 1
    )
    .ok();
    writeln!(
        f,
        "# mp_cosine={:.10} mp_rr={:.10} mp_r2={:.10} mp_rc={:.10}",
        cos, res, min_r2, resid_corr
    )
    .ok();
    writeln!(f, "# ELUTION PROFILE (col_effects)").ok();
    for (i, v) in mp.col_effects.iter().enumerate() {
        writeln!(f, "elution\t{}\t{:.10}", i, v).ok();
    }
    writeln!(f, "# FRAGMENT EFFECTS (row_effects)").ok();
    for (i, v) in mp.row_effects.iter().enumerate() {
        writeln!(f, "frag_effect\t{}\t{:.10}", i, v).ok();
    }
    writeln!(f, "# grand_mean={:.10}", mp.overall).ok();
    writeln!(
        f,
        "# n_iterations={} converged={}",
        mp.n_iterations, mp.converged
    )
    .ok();
    writeln!(f, "# INPUT MATRIX (frag_idx, scan_idx, value)").ok();
    for (xi, (_, xic_data)) in peak_xics.iter().enumerate() {
        for (s, (_, v)) in xic_data.iter().enumerate() {
            writeln!(f, "input\t{}\t{}\t{:.10}", xi, s, v).ok();
        }
    }
    log::info!("[BISECT] Wrote median polish diagnostic: rust_mp_diag.txt");
    let _ = exit_if_only; // not used here -- mp_diag has no _ONLY variant
}

// --------------------------------------------------------------------------
// Per-entry main-search XIC dump
// --------------------------------------------------------------------------

/// State for the per-entry main-search XIC diagnostic dump.
///
/// Reads `OSPREY_DIAG_SEARCH_ENTRY_IDS=<id1,id2,...>` once at construction.
/// Unlike [`super::diagnostics::CalXicEntryDump`], this dump does NOT exit
/// after writing — it accumulates rust_search_xic_entry_<id>.txt files for
/// every listed entry id encountered during the main search and lets the
/// pipeline run to completion so downstream analysis is also available.
pub struct SearchXicDump {
    target_ids: Option<std::collections::HashSet<u32>>,
}

impl Default for SearchXicDump {
    fn default() -> Self {
        Self::new()
    }
}

impl SearchXicDump {
    pub fn new() -> Self {
        let target_ids: Option<std::collections::HashSet<u32>> =
            std::env::var("OSPREY_DIAG_SEARCH_ENTRY_IDS").ok().map(|s| {
                let ids: std::collections::HashSet<u32> = s
                    .split(',')
                    .filter_map(|p| p.trim().parse::<u32>().ok())
                    .collect();
                log::info!(
                    "[BISECT] OSPREY_DIAG_SEARCH_ENTRY_IDS: will dump {} entries",
                    ids.len()
                );
                ids
            });
        Self { target_ids }
    }

    fn is_active_for(&self, entry_id: u32) -> bool {
        self.target_ids
            .as_ref()
            .is_some_and(|ids| ids.contains(&entry_id))
    }

    /// Write the header section of the search XIC dump (CANDIDATES +
    /// EXTRACTED XICS) to `rust_search_xic_entry_<id>.txt`. No-op if the
    /// entry is not in the target list.
    pub fn dump_header(
        &self,
        entry: &LibraryEntry,
        expected_rt: f64,
        rt_tolerance: f64,
        cand_spectra: &[&Spectrum],
        xics: &[(usize, Vec<(f64, f64)>)],
    ) {
        if !self.is_active_for(entry.id) {
            return;
        }
        let dump_path = format!("rust_search_xic_entry_{}.txt", entry.id);
        let Ok(mut f) = std::fs::File::create(&dump_path) else {
            return;
        };
        writeln!(f, "# search XIC dump for entry_id={}", entry.id).ok();
        writeln!(
            f,
            "# {} ({}, charge={}, lib_rt={:.10}, mz={:.10})",
            entry.modified_sequence,
            entry.sequence,
            entry.charge,
            entry.retention_time,
            entry.precursor_mz
        )
        .ok();
        writeln!(f, "# is_decoy={}", if entry.is_decoy { 1 } else { 0 }).ok();
        writeln!(f, "# expected_rt={:.10}", expected_rt).ok();
        writeln!(f, "# rt_tolerance={:.10}", rt_tolerance).ok();
        writeln!(
            f,
            "# scan_range=[0..{}] n_scans={}",
            cand_spectra.len().saturating_sub(1),
            cand_spectra.len()
        )
        .ok();
        writeln!(f, "# CANDIDATES (scan_idx, scan_number, rt)").ok();
        writeln!(f, "candidate\tscan_idx\tscan_number\trt").ok();
        for (i, spec) in cand_spectra.iter().enumerate() {
            writeln!(
                f,
                "candidate\t{}\t{}\t{:.10}",
                i, spec.scan_number, spec.retention_time
            )
            .ok();
        }
        writeln!(f, "# EXTRACTED XICS (lib_idx, scan_idx, rt, intensity)").ok();
        writeln!(f, "xic\tlib_idx\tscan_idx\trt\tintensity").ok();
        for (frag_idx, xic_data) in xics {
            for (i, (rt, intensity)) in xic_data.iter().enumerate() {
                writeln!(f, "xic\t{}\t{}\t{:.10}\t{:.10}", frag_idx, i, rt, intensity).ok();
            }
        }
        log::info!(
            "[BISECT] Search XIC dump for entry {}: {} xics, {} scans -> {}",
            entry.id,
            xics.len(),
            cand_spectra.len(),
            dump_path
        );
    }

    /// Append the CWT PEAKS + BEST PEAK section to the search XIC dump
    /// for the given entry. No-op if the entry is not in the target list.
    /// `scored_candidates` is `(bounds, raw_coelution_score, rt_penalized_score)`
    /// ordered by `rt_penalized_score` descending (the sort key used by
    /// `pipeline.rs`). The dump writes the raw score so the column matches
    /// what fork/C# bisection tooling produces; candidate *ordering* in the
    /// dump reflects the RT-penalized sort that upstream `main` uses.
    pub fn dump_peaks(
        &self,
        entry: &LibraryEntry,
        scored_candidates: &[(&XICPeakBounds, f64, f64)],
    ) {
        if !self.is_active_for(entry.id) {
            return;
        }
        if scored_candidates.is_empty() {
            return;
        }
        let dump_path = format!("rust_search_xic_entry_{}.txt", entry.id);
        let Ok(mut f) = std::fs::OpenOptions::new().append(true).open(&dump_path) else {
            return;
        };
        writeln!(f, "# CWT PEAKS: {} candidates", scored_candidates.len()).ok();
        writeln!(f, "peak\tidx\tstart\tapex\tend\tcorr_score").ok();
        for (pi, (bp, raw_score, _penalized_score)) in scored_candidates.iter().enumerate() {
            writeln!(
                f,
                "peak\t{}\t{}\t{}\t{}\t{:.10}",
                pi, bp.start_index, bp.apex_index, bp.end_index, raw_score
            )
            .ok();
        }
        let best_bp = scored_candidates[0].0;
        writeln!(
            f,
            "# BEST PEAK: idx=0 start={} apex={} end={}",
            best_bp.start_index, best_bp.apex_index, best_bp.end_index
        )
        .ok();
    }
}
