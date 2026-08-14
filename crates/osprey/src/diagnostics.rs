//! Cross-implementation bisection diagnostic dumps for the pipeline stage.
//!
//! Each function below is gated by an `OSPREY_DUMP_*` or `OSPREY_DIAG_*`
//! env var and writes a `rust_*.txt` file matching the format of the
//! corresponding C# OspreySharp dump (`cs_*.txt`). Dumps that have an
//! `_ONLY` companion env var exit the process after writing for fast
//! cycle-time during bisection.
//!
//! See `osprey-core::diagnostics` for shared primitives and the env-var
//! gating convention.

use osprey_core::diagnostics::{exit_if_only, format_f64_roundtrip, is_dump_enabled};
use osprey_core::{LibraryEntry, Spectrum, XICPeakBounds};
use osprey_fdr::protein::{PeptideScore, ProteinFdrResult, ProteinParsimonyResult};
use osprey_scoring::{SpectralScorer, TukeyMedianPolishResult};
use std::collections::HashMap;
use std::io::Write;
use std::sync::Arc;

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

    let n_lib = library_rts_detected.len();
    let n_meas = measured_rts_detected.len();
    if n_lib != n_meas {
        log::warn!(
            "LOESS input length mismatch: library_rts={} measured_rts={} (zip will truncate to {})",
            n_lib,
            n_meas,
            n_lib.min(n_meas)
        );
    }

    if let Ok(mut f) = std::fs::File::create("rust_loess_input.txt") {
        writeln!(f, "# n_library_rts={} n_measured_rts={}", n_lib, n_meas).ok();
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
    writeln!(
        f,
        "# nbins={} xcorr_scaled={:.10}",
        pre_vec.len(),
        xcorr_scaled
    )
    .ok();
    let psum: f64 = pre_vec.iter().map(|&v| v as f64).sum();
    let pnz = pre_vec.iter().filter(|&&v| v != 0.0).count();
    writeln!(f, "# preprocessed_sum={:.10} nonzero={}", psum, pnz).ok();
    // First 20 nonzero preprocessed bins
    let mut dumped = 0;
    for (i, &v) in pre_vec.iter().enumerate() {
        if v != 0.0 {
            writeln!(f, "pre\t{}\t{:.10}", i, v).ok();
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
            "frag\t{}\tmz={:.10}\tbin={:?}\tpre_val={:.10}\tlib_val={:.10}",
            fi, frag.mz, bin, pre_val, lib_val
        )
        .ok();
    }
}

/// Write a median polish diagnostic block to `rust_mp_diag.txt`.
///
/// Gated by `OSPREY_DIAG_MP_SCAN=<scan_number>` matching
/// `apex_spectrum.scan_number`. Writes the row/col effects, grand mean,
/// convergence info, and the input matrix that the polish was run on.
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
    if format!("{}", apex_spectrum.scan_number) != diag_scan {
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
}

// --------------------------------------------------------------------------
// Bisection: per-entry inputs to tukey_median_polish
// --------------------------------------------------------------------------
//
// Captures the (frag_idx, scan_idx, rt, intensity) matrix that goes into
// every `tukey_median_polish` call from the rescore path, keyed by
// (entry_id, apex_scan_number). Gated by
// `OSPREY_DUMP_MP_INPUTS=<path>`; thread-safe append + lazy file open
// behind a `Mutex<File>` so it can fire from parallel rayon scoring.
//
// Workflow:
//   OSPREY_DUMP_MP_INPUTS=/tmp/inproc_mp_inputs.tsv  osprey ... # in-process
//   OSPREY_DUMP_MP_INPUTS=/tmp/worker_mp_inputs.tsv  osprey ... # worker
//   sort -k1,1n -k2,2n -k3,3n /tmp/inproc_mp_inputs.tsv > /tmp/inproc.sorted
//   sort -k1,1n -k2,2n -k3,3n /tmp/worker_mp_inputs.tsv > /tmp/worker.sorted
//   diff /tmp/inproc.sorted /tmp/worker.sorted
//
// PASS = byte-identical sorted dumps. FAIL = the very first divergence
// between in-process and worker at the median-polish input layer.

use std::sync::{Mutex, OnceLock};

static MP_INPUTS_DUMP: OnceLock<Option<Mutex<std::fs::File>>> = OnceLock::new();

fn mp_inputs_writer() -> Option<&'static Mutex<std::fs::File>> {
    MP_INPUTS_DUMP
        .get_or_init(|| {
            let path = std::env::var("OSPREY_DUMP_MP_INPUTS").ok()?;
            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .truncate(true)
                .write(true)
                .open(&path)
                .ok()?;
            // Header comment; tools that sort the rows by the first three
            // numeric columns (entry_id, apex_scan, frag_pos) stay
            // hash-stable across runs regardless of rayon scheduling.
            writeln!(
                file,
                "# entry_id\tapex_scan\tfrag_pos\tfrag_idx\tscan_idx\trt\tintensity"
            )
            .ok();
            log::info!(
                "[BISECT] OSPREY_DUMP_MP_INPUTS active: writing tukey_median_polish inputs to {}",
                path
            );
            Some(Mutex::new(file))
        })
        .as_ref()
}

/// Dump the inputs to a single `tukey_median_polish` call for cross-
/// implementation bisection. `frag_pos` is the position of each
/// fragment within `peak_xics` (so order-sensitive bugs surface);
/// `frag_idx` is the library fragment index (so the same fragment can
/// be matched between dumps regardless of `peak_xics` ordering).
pub fn dump_mp_inputs(entry_id: u32, apex_scan: u32, peak_xics: &[(usize, Vec<(f64, f64)>)]) {
    let Some(writer) = mp_inputs_writer() else {
        return;
    };
    let Ok(mut file) = writer.lock() else {
        return;
    };
    // Build the dump for this call into a single buffer so the lock is
    // held for one write, not 6 × N_scans writes per entry.
    let mut buf = String::new();
    for (frag_pos, (frag_idx, xic)) in peak_xics.iter().enumerate() {
        for (scan_idx, (rt, intensity)) in xic.iter().enumerate() {
            buf.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                entry_id,
                apex_scan,
                frag_pos,
                frag_idx,
                scan_idx,
                format_f64_roundtrip(*rt),
                format_f64_roundtrip(*intensity),
            ));
        }
    }
    let _ = file.write_all(buf.as_bytes());
}

// --------------------------------------------------------------------------
// Bisection: per-call inputs + outputs to RTCalibration::predict
// --------------------------------------------------------------------------
//
// Captures (entry_id, library_rt_input, expected_rt_output) for every
// `cal.predict(entry.retention_time)` call from the rescore search
// path, plus a one-time dump per (file_name) of the calibration's
// `library_rts` and `fitted_values` arrays. Together these tell us:
//
// - If the per-call output differs between in-process and worker for
//   the SAME entry_id (same library_rt input, since both sides load
//   the same library), then the calibration object differs.
// - If the per-file array dump differs between sides, the JSON
//   round-trip through `format_f64_roundtrip` is not preserving f64
//   bits exactly for some values — fix the serialization.
// - If arrays match but per-call output differs, some hidden state in
//   `RTCalibration` is affecting `predict()` beyond what
//   `from_model_params` reconstructs.
//
// Gated by `OSPREY_DUMP_PREDICT_RT=<path>`. Two-section TSV output
// sharing the same five-column header. The third column ("array_or_apex")
// is the array kind for cal_arrays rows and a literal `-` placeholder
// for predict_calls rows (reserved for a future apex_scan extension if
// the per-call site ever has it in scope):
//
//     # SECTION cal_arrays
//     #   "cal_arrays" | file_name | array_kind ("library_rts" | "fitted_values") | idx | value
//     # SECTION predict_calls
//     #   "predict_calls" | entry_id | "-" | library_rt | expected_rt
//
// Both sections are append-only (the cal_arrays writer does not dedup
// per file_name; if the rescore loop reuses a file's calibration twice,
// expect duplicate cal_arrays rows). The bisection workflow handles
// duplicates with `sort -u`:
//
//     OSPREY_DUMP_PREDICT_RT=/tmp/inproc_rt.tsv  osprey ...   # in-process
//     OSPREY_DUMP_PREDICT_RT=/tmp/worker_rt.tsv  osprey ...   # worker
//     sort -u /tmp/inproc_rt.tsv > /tmp/i.sorted
//     sort -u /tmp/worker_rt.tsv > /tmp/w.sorted
//     diff /tmp/i.sorted /tmp/w.sorted

static PREDICT_RT_DUMP: OnceLock<Option<Mutex<std::fs::File>>> = OnceLock::new();

fn predict_rt_writer() -> Option<&'static Mutex<std::fs::File>> {
    PREDICT_RT_DUMP
        .get_or_init(|| {
            let path = std::env::var("OSPREY_DUMP_PREDICT_RT").ok()?;
            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .truncate(true)
                .write(true)
                .open(&path)
                .ok()?;
            writeln!(
                file,
                "# section\tfile_name_or_entry_id\tarray_or_apex\tidx_or_lib_rt\tvalue_or_expected_rt"
            )
            .ok();
            log::info!(
                "[BISECT] OSPREY_DUMP_PREDICT_RT active: writing predict() inputs/outputs to {}",
                path
            );
            Some(Mutex::new(file))
        })
        .as_ref()
}

/// Dump the cal's `library_rts` and `fitted_values` arrays once per
/// file, full f64 precision. Append-only — the writer does NOT dedup
/// per file_name, so repeated calls for the same file produce
/// repeated rows. The bisection workflow uses `sort -u` post-hoc to
/// collapse duplicates. Call from the top of the rescore loop where
/// the cal is first selected.
pub fn dump_predict_rt_arrays(file_name: &str, library_rts: &[f64], fitted_values: &[f64]) {
    let Some(writer) = predict_rt_writer() else {
        return;
    };
    let Ok(mut file) = writer.lock() else {
        return;
    };
    let mut buf = String::new();
    for (i, v) in library_rts.iter().enumerate() {
        buf.push_str(&format!(
            "cal_arrays\t{}\tlibrary_rts\t{}\t{}\n",
            file_name,
            i,
            format_f64_roundtrip(*v)
        ));
    }
    for (i, v) in fitted_values.iter().enumerate() {
        buf.push_str(&format!(
            "cal_arrays\t{}\tfitted_values\t{}\t{}\n",
            file_name,
            i,
            format_f64_roundtrip(*v)
        ));
    }
    let _ = file.write_all(buf.as_bytes());
}

/// Dump a single `predict()` input + output pair. Called from the
/// rescore search path right after `expected_rt = cal.predict(...)`.
/// The same entry_id may appear many times (once per per-window
/// candidate scoring iteration); deduplicate post-hoc with
/// `sort -u`.
pub fn dump_predict_rt_call(entry_id: u32, library_rt: f64, expected_rt: f64) {
    let Some(writer) = predict_rt_writer() else {
        return;
    };
    let Ok(mut file) = writer.lock() else {
        return;
    };
    let line = format!(
        "predict_calls\t{}\t-\t{}\t{}\n",
        entry_id,
        format_f64_roundtrip(library_rt),
        format_f64_roundtrip(expected_rt),
    );
    let _ = file.write_all(line.as_bytes());
}

// --------------------------------------------------------------------------
// Per-(file, entry) CWT path summary for cross-impl bisection
// --------------------------------------------------------------------------
//
// Captures one row per scoring call describing the CWT detection /
// fallback / apex-acceptance pipeline outcome. Ten tab-separated
// columns: file_name, entry_id, n_cwt_peaks, n_final_peaks,
// n_scored, scored, sigma, consensus_l1, consensus_max_abs,
// consensus_argmax. The four consensus-signal stats (sigma + L1 +
// max-abs + argmax of the median CWT consensus across fragment XICs)
// give a single-number signature that distinguishes "convolve /
// median diverged" from "downstream peak finder diverged" without
// dumping the full per-scan consensus signal. Use to localize where
// the rescore set diverges between Rust and the C# port.
//
//   OSPREY_DUMP_CWT_PATH=/tmp/rust_cwt.tsv  osprey ...   # Rust
//   OSPREY_DUMP_CWT_PATH=1                  OspreySharp  # C# (writes cs_stage6_cwt_path.tsv)

static CWT_PATH_DUMP: OnceLock<Option<Mutex<std::fs::File>>> = OnceLock::new();

fn cwt_path_writer() -> Option<&'static Mutex<std::fs::File>> {
    CWT_PATH_DUMP
        .get_or_init(|| {
            let path = std::env::var("OSPREY_DUMP_CWT_PATH").ok()?;
            let mut file = std::fs::OpenOptions::new()
                .create(true)
                .truncate(true)
                .write(true)
                .open(&path)
                .ok()?;
            writeln!(
                file,
                "file_name\tentry_id\tn_cwt_peaks\tn_final_peaks\tn_scored\tscored\tsigma\tconsensus_l1\tconsensus_max_abs\tconsensus_argmax"
            )
            .ok();
            log::info!(
                "[BISECT] OSPREY_DUMP_CWT_PATH active: writing CWT path summary to {}",
                path
            );
            Some(Mutex::new(file))
        })
        .as_ref()
}

/// Append one per-(file, entry) CWT path summary row. The four
/// consensus-signal scalars (sigma, l1, max_abs, argmax) are computed
/// inside this function from `xics` so the production hot path
/// doesn't pay for the recompute when `OSPREY_DUMP_CWT_PATH` is
/// unset (the OnceLock writer returns None and the function bails
/// before touching xics).
pub fn dump_cwt_path(
    file_name: &str,
    entry_id: u32,
    n_cwt_peaks: usize,
    n_final_peaks: usize,
    n_scored: usize,
    scored: bool,
    xics: &[(usize, Vec<(f64, f64)>)],
) {
    let Some(writer) = cwt_path_writer() else {
        return;
    };
    let (sigma, l1, max_abs, argmax) = match osprey_chromatography::cwt::get_consensus_signal(xics)
    {
        Some((cons, sig)) => {
            let mut l1 = 0.0f64;
            let mut max_abs = 0.0f64;
            let mut argmax: i32 = -1;
            for (i, v) in cons.iter().enumerate() {
                let a = v.abs();
                l1 += a;
                if a > max_abs {
                    max_abs = a;
                    argmax = i as i32;
                }
            }
            (sig, l1, max_abs, argmax)
        }
        None => (0.0, 0.0, 0.0, -1),
    };
    let Ok(mut file) = writer.lock() else {
        return;
    };
    let line = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        file_name,
        entry_id,
        n_cwt_peaks,
        n_final_peaks,
        n_scored,
        if scored { 1 } else { 0 },
        format_f64_roundtrip(sigma),
        format_f64_roundtrip(l1),
        format_f64_roundtrip(max_abs),
        argmax
    );
    let _ = file.write_all(line.as_bytes());
}

// --------------------------------------------------------------------------
// Per-entry main-search XIC dump
// --------------------------------------------------------------------------

/// State for the per-entry main-search XIC diagnostic dump.
///
/// Reads `OSPREY_DIAG_SEARCH_ENTRY_IDS=<id1,id2,...>` once at construction.
/// Unlike `osprey_scoring::diagnostics::CalXicEntryDump`, this dump does
/// NOT exit after writing -- it accumulates rust_search_xic_entry_<id>.txt
/// files for every listed entry id encountered during the main search and
/// lets the pipeline run to completion so downstream analysis is also
/// available.
/// Module-level entry point for the per-fragment apex-match dump used in
/// the Astral Mode A bisection. Lazy-inits a `SearchXicDump` so the
/// `OSPREY_DIAG_SEARCH_ENTRY_IDS` env-var parse happens at most once per
/// process, then forwards to `SearchXicDump::dump_fragment_match`.
pub fn dump_fragment_match(
    entry: &LibraryEntry,
    apex_spectrum: &Spectrum,
    tol_da: f64,
    tol_ppm: f64,
) {
    use std::sync::OnceLock;
    static DUMPER: OnceLock<SearchXicDump> = OnceLock::new();
    DUMPER.get_or_init(SearchXicDump::new).dump_fragment_match(
        entry,
        apex_spectrum,
        tol_da,
        tol_ppm,
    );
}

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

        // CWT CONSENSUS: per-scan median consensus value across the
        // fragment CWT coefficients. Cross-impl diff at this section
        // pinpoints the first scan where the consensus signal diverges
        // -- the seam upstream of peak detection. Use
        // format_f64_roundtrip so f64 bits compare exactly between
        // Rust and C#.
        if let Some((consensus, sigma)) = osprey_chromatography::cwt::get_consensus_signal(xics) {
            writeln!(f, "# CWT CONSENSUS (sigma, scan_idx, value)").ok();
            writeln!(f, "# sigma={}", format_f64_roundtrip(sigma)).ok();
            writeln!(f, "consensus\tscan_idx\tvalue").ok();
            for (i, v) in consensus.iter().enumerate() {
                writeln!(f, "consensus\t{}\t{}", i, format_f64_roundtrip(*v)).ok();
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

    /// Cross-impl bisection helper: dump per-fragment match info at the
    /// apex spectrum to `rust_fragmatch_entry_<id>_scan_<scan>.txt`. Mirrors
    /// the C# `cs_fragmatch_entry_<id>_scan_<scan>.txt` format so the two
    /// can be diffed directly. Each row is the library fragment annotation
    /// (ion_type/ordinal/charge), library m/z, tolerance window, and the
    /// matched flag (1 if at least one observed peak falls in the window).
    /// For each fragment, also writes one `peak` line per observed peak in
    /// the tolerance window so we can inspect what spectrum content the
    /// scorer actually saw — and whether it differs cross-impl after
    /// per-file MS2 calibration is applied.
    pub fn dump_fragment_match(
        &self,
        entry: &LibraryEntry,
        apex_spectrum: &Spectrum,
        tol_da: f64,
        tol_ppm: f64,
    ) {
        if !self.is_active_for(entry.id) {
            return;
        }
        let dump_path = format!(
            "rust_fragmatch_entry_{}_scan_{}.txt",
            entry.id, apex_spectrum.scan_number
        );
        let Ok(mut f) = std::fs::File::create(&dump_path) else {
            return;
        };
        let unit = if tol_ppm > 0.0 { "Ppm" } else { "Mz" };
        let tol_val = if tol_ppm > 0.0 { tol_ppm } else { tol_da };
        writeln!(
            f,
            "# entry_id={} scan={} modseq={} seq={} is_decoy={} n_fragments={} n_peaks={} tol={} {}",
            entry.id,
            apex_spectrum.scan_number,
            entry.modified_sequence,
            entry.sequence,
            if entry.is_decoy { 1 } else { 0 },
            entry.fragments.len(),
            apex_spectrum.mzs.len(),
            tol_val,
            unit
        )
        .ok();
        writeln!(
            f,
            "ion_type\tordinal\tcharge\tlib_mz\ttol_da\tlower\tupper\tmatched"
        )
        .ok();
        for frag in &entry.fragments {
            let frag_tol = tol_da.max(frag.mz * tol_ppm / 1e6);
            let lower = frag.mz - frag_tol;
            let upper = frag.mz + frag_tol;
            // Mirror match_fragments inclusion logic: any peak in [lower, upper] counts.
            let matched = apex_spectrum
                .mzs
                .iter()
                .any(|&mz| mz >= lower && mz <= upper);
            writeln!(
                f,
                "{:?}\t{}\t{}\t{:.10}\t{:.10}\t{:.10}\t{:.10}\t{}",
                frag.annotation.ion_type,
                frag.annotation.ordinal,
                frag.annotation.charge,
                frag.mz,
                frag_tol,
                lower,
                upper,
                if matched { 1 } else { 0 }
            )
            .ok();
            for (i, &mz) in apex_spectrum.mzs.iter().enumerate() {
                if mz < lower {
                    continue;
                }
                if mz > upper {
                    break;
                }
                let inten = apex_spectrum.intensities.get(i).copied().unwrap_or(0.0_f32);
                writeln!(
                    f,
                    "  peak\t{:?}\t{}\t{:.17}\t{:.6}",
                    frag.annotation.ion_type, frag.annotation.ordinal, mz, inten
                )
                .ok();
            }
        }
    }
}

/// Dump per-precursor Stage 5 (Percolator FDR) state to
/// `rust_stage5_percolator.tsv` so cross-impl parity can be checked at
/// end-of-first-pass-FDR, before compaction.
///
/// Gated by `OSPREY_DUMP_PERCOLATOR=1`. When `OSPREY_PERCOLATOR_ONLY=1`
/// is also set, exits the process after writing. Columns:
/// `file_name, entry_id, charge, modified_sequence, is_decoy, score, pep,
/// run_precursor_q, run_peptide_q, run_protein_q, experiment_precursor_q,
/// experiment_peptide_q`. Rows sorted by `(file_name, entry_id)` for
/// stable human inspection; `Compare-Percolator.ps1` hash-joins on the
/// composite key and is sort-order-agnostic.
///
/// `run_protein_qvalue` is the default `1.0` when this dump fires from
/// `pipeline.rs` (the dump runs BEFORE first-pass protein FDR populates
/// real values), and the real persisted value when this dump fires from
/// `rescore::run_rescore` (the worker hydrates the v3 sidecar which
/// carries post-protein-FDR values). The C# `WriteStage5PercolatorDump`
/// has the same dual-call shape and the same column ordering.
pub fn dump_stage5_percolator(per_file_entries: &[(String, Vec<osprey_core::FdrEntry>)]) {
    if !is_dump_enabled("OSPREY_DUMP_PERCOLATOR") {
        return;
    }

    let path = "rust_stage5_percolator.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Without this the ~74 MB TSV took ~14 minutes to write on WSL
    // drvfs/9P (one syscall per row); buffering completes in seconds.
    // Diagnostic-only output; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "file_name\tentry_id\tcharge\tmodified_sequence\tis_decoy\tscore\tpep\trun_precursor_q\trun_peptide_q\trun_protein_q\texperiment_precursor_q\texperiment_peptide_q"
    )
    .ok();

    let mut rows: Vec<(&str, &osprey_core::FdrEntry)> = per_file_entries
        .iter()
        .flat_map(|(file_name, entries)| entries.iter().map(move |e| (file_name.as_str(), e)))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0).then(a.1.entry_id.cmp(&b.1.entry_id)));

    let mut n_written = 0usize;
    for (file_name, e) in &rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            file_name,
            e.entry_id,
            e.charge,
            e.modified_sequence,
            if e.is_decoy { "true" } else { "false" },
            format_f64_roundtrip(e.score),
            format_f64_roundtrip(e.pep),
            format_f64_roundtrip(e.run_precursor_qvalue),
            format_f64_roundtrip(e.run_peptide_qvalue),
            format_f64_roundtrip(e.run_protein_qvalue),
            format_f64_roundtrip(e.experiment_precursor_qvalue),
            format_f64_roundtrip(e.experiment_peptide_qvalue),
        )
        .ok();
        n_written += 1;
    }
    log::info!(
        "Wrote Stage 5 Percolator dump: {} ({} rows)",
        path,
        n_written
    );

    // Flush BufWriter before exit_if_only -- process::exit(0) skips
    // destructors so the 8 MB buffer would otherwise truncate the file.
    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_PERCOLATOR_ONLY", "Stage 5 Percolator dump");
}

/// Dump per-precursor state to `rust_stage6_rescored.tsv` AFTER the
/// per-file rescore loop completes. Same column shape as
/// `dump_stage5_percolator` so `Compare-Percolator.ps1` can be reused
/// for the diff. Catches divergence in the boundary-overrides
/// rescore + gap-fill output BEFORE drilling into the inner-loop
/// `OSPREY_DUMP_MP_INPUTS` / `OSPREY_DUMP_PREDICT_RT` bisection
/// ladder.
///
/// Gated by `OSPREY_DUMP_RESCORED=1`. When `OSPREY_RESCORED_ONLY=1`
/// is also set, exits the process after writing. Fired from BOTH
/// the in-process pipeline (after `rescore_per_file_loop` in
/// `pipeline.rs::run_analysis`) and the worker
/// (`rescore::run_rescore`); the OspreySharp side has the matching
/// `WriteStage6RescoredDump` call from `AnalysisPipeline.Run` and
/// `AnalysisPipeline.RunWorker`.
pub fn dump_stage6_rescored(per_file_entries: &[(String, Vec<osprey_core::FdrEntry>)]) {
    if !is_dump_enabled("OSPREY_DUMP_RESCORED") {
        return;
    }

    let path = "rust_stage6_rescored.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "file_name\tentry_id\tcharge\tmodified_sequence\tis_decoy\tscore\tpep\trun_precursor_q\trun_peptide_q\trun_protein_q\texperiment_precursor_q\texperiment_peptide_q"
    )
    .ok();

    let mut rows: Vec<(&str, &osprey_core::FdrEntry)> = per_file_entries
        .iter()
        .flat_map(|(file_name, entries)| entries.iter().map(move |e| (file_name.as_str(), e)))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0).then(a.1.entry_id.cmp(&b.1.entry_id)));

    let mut n_written = 0usize;
    for (file_name, e) in &rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            file_name,
            e.entry_id,
            e.charge,
            e.modified_sequence,
            if e.is_decoy { "true" } else { "false" },
            format_f64_roundtrip(e.score),
            format_f64_roundtrip(e.pep),
            format_f64_roundtrip(e.run_precursor_qvalue),
            format_f64_roundtrip(e.run_peptide_qvalue),
            format_f64_roundtrip(e.run_protein_qvalue),
            format_f64_roundtrip(e.experiment_precursor_qvalue),
            format_f64_roundtrip(e.experiment_peptide_qvalue),
        )
        .ok();
        n_written += 1;
    }
    log::info!("Wrote Stage 6 rescored dump: {} ({} rows)", path, n_written);

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_RESCORED_ONLY", "Stage 6 rescored dump");
}

/// Dump the per-peptide consensus RT planning state to
/// `rust_stage6_consensus.tsv` for cross-impl parity at the start of
/// Stage 6, before any per-file rescoring. Mirrors
/// OspreySharp.FDR.Reconciliation.ConsensusRts.Compute.
///
/// Gated by `OSPREY_DUMP_CONSENSUS=1`. When `OSPREY_CONSENSUS_ONLY=1`
/// is also set, exits the process after writing. Columns:
/// `is_decoy, modified_sequence, consensus_library_rt,
/// median_peak_width, n_runs_detected, apex_library_rt_mad`.
/// Rows are emitted in the order produced by `compute_consensus_rts`,
/// which sorts by `(is_decoy, modified_sequence)` for deterministic
/// output. `apex_library_rt_mad` is empty when fewer than 3 detections
/// contributed.
pub fn dump_stage6_consensus(consensus: &[crate::reconciliation::PeptideConsensusRT]) {
    if !is_dump_enabled("OSPREY_DUMP_CONSENSUS") {
        return;
    }

    let path = "rust_stage6_consensus.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "is_decoy\tmodified_sequence\tconsensus_library_rt\tmedian_peak_width\tn_runs_detected\tapex_library_rt_mad"
    )
    .ok();

    for c in consensus {
        let mad_str = c
            .apex_library_rt_mad
            .map(format_f64_roundtrip)
            .unwrap_or_default();
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            if c.is_decoy { "true" } else { "false" },
            c.modified_sequence,
            format_f64_roundtrip(c.consensus_library_rt),
            format_f64_roundtrip(c.median_peak_width),
            c.n_runs_detected,
            mad_str,
        )
        .ok();
    }
    log::info!(
        "Wrote Stage 6 consensus dump: {} ({} rows)",
        path,
        consensus.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_CONSENSUS_ONLY", "Stage 6 consensus dump");
}

/// Dump the per-file multi-charge consensus rescore targets to
/// `rust_stage6_multicharge.tsv` for cross-impl parity. Mirrors
/// `OspreySharp.FDR.Reconciliation.MultiChargeConsensus.SelectRescoreTargets`.
///
/// Gated by `OSPREY_DUMP_MULTICHARGE=1`. When
/// `OSPREY_MULTICHARGE_ONLY=1` is also set, exits the process after
/// writing. Columns: `file_name, entry_id, consensus_apex,
/// consensus_start, consensus_end`. Rows sorted by
/// `(file_name, entry_id)` for stable diff. The dump uses the stable
/// library entry id (as written to the parquet `entry_id` column)
/// instead of the per-file Vec position, so cross-impl comparison is
/// invariant to whether the implementation has compacted the
/// per-file FdrEntry list before computing multi-charge consensus.
pub fn dump_stage6_multicharge(
    per_file_entries: &[(String, Vec<osprey_core::FdrEntry>)],
    per_file_consensus_targets: &HashMap<String, Vec<(usize, f64, f64, f64)>>,
) {
    if !is_dump_enabled("OSPREY_DUMP_MULTICHARGE") {
        return;
    }

    let path = "rust_stage6_multicharge.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "file_name\tentry_id\tconsensus_apex\tconsensus_start\tconsensus_end"
    )
    .ok();

    // Resolve idx -> entry_id via the matching per_file_entries list so the
    // dump key is the stable library entry id rather than the (compaction-
    // dependent) Vec position.
    let entries_by_file: HashMap<&str, &Vec<osprey_core::FdrEntry>> = per_file_entries
        .iter()
        .map(|(name, entries)| (name.as_str(), entries))
        .collect();

    let mut rows: Vec<(&str, u32, f64, f64, f64)> = per_file_consensus_targets
        .iter()
        .flat_map(|(file_name, targets)| {
            let entries = entries_by_file.get(file_name.as_str()).copied();
            targets.iter().filter_map(move |&(idx, apex, start, end)| {
                entries
                    .and_then(|e| e.get(idx))
                    .map(|entry| (file_name.as_str(), entry.entry_id, apex, start, end))
            })
        })
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0).then(a.1.cmp(&b.1)));

    for (file_name, entry_id, apex, start, end) in &rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            file_name,
            entry_id,
            format_f64_roundtrip(*apex),
            format_f64_roundtrip(*start),
            format_f64_roundtrip(*end),
        )
        .ok();
    }
    log::info!(
        "Wrote Stage 6 multi-charge dump: {} ({} rows)",
        path,
        rows.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_MULTICHARGE_ONLY", "Stage 6 multi-charge dump");
}

/// Append the loaded calibration arrays (library_rts + fitted_values)
/// for one file to `rust_stage6_calibration.tsv` for cross-impl
/// JSON-decode parity check. Each call writes one row per
/// `(library_rts[i], fitted_values[i])` pair, prefixed with the file
/// name. The header is written on the first call (file does not yet
/// exist) and subsequent calls append.
///
/// Diffing this dump against `cs_stage6_calibration.tsv` localizes
/// whether cross-impl divergence in `inverse_predict` enters at the
/// JSON `f64` parser (decimal-to-binary) or inside the LOESS
/// interpolation arithmetic.
///
/// Gated by `OSPREY_DUMP_CALIBRATION=1`. The companion
/// `OSPREY_CALIBRATION_ONLY=1` env-var is checked after the
/// `--join-only` calibration-load loop completes — when all files
/// have appended their rows — and exits the process cleanly.
pub fn dump_stage6_calibration(file_name: &str, library_rts: &[f64], fitted_values: &[f64]) {
    if !is_dump_enabled("OSPREY_DUMP_CALIBRATION") {
        return;
    }

    let path = "rust_stage6_calibration.tsv";
    let header_needed = !std::path::Path::new(path).exists();

    let Ok(mut f) = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
    else {
        log::warn!("Could not open {}", path);
        return;
    };

    if header_needed {
        writeln!(f, "file_name\tidx\tlibrary_rt\tfitted_value").ok();
    }

    let n = library_rts.len().min(fitted_values.len());
    for i in 0..n {
        writeln!(
            f,
            "{}\t{}\t{}\t{}",
            file_name,
            i,
            format_f64_roundtrip(library_rts[i]),
            format_f64_roundtrip(fitted_values[i]),
        )
        .ok();
    }
    log::info!(
        "Appended {} calibration rows for {} to {}",
        n,
        file_name,
        path
    );
}

/// Dump per-detection (apex_rt, library_rt) pairs going into
/// `compute_consensus_rts` to `rust_stage6_inv_predict.tsv` for
/// cross-impl ULP bisection. apex_rt is the FdrEntry value loaded
/// from Parquet; library_rt is `inverse_predict(apex_rt)` against the
/// per-file refined RT calibration. weight is the sigmoid-of-SVM-score
/// Stage 6 contribution weight.
///
/// Diffing this dump against `cs_stage6_inv_predict.tsv` localizes
/// whether cross-impl divergence in `consensus_library_rt` enters at
/// Parquet f64 decode (apex_rt diverges) or LOESS interpolation
/// (only library_rt diverges).
///
/// Each record is `(file_name, modified_sequence, is_decoy, apex_rt,
/// library_rt, weight)`.
///
/// Gated by `OSPREY_DUMP_INV_PREDICT=1`. When `OSPREY_INV_PREDICT_ONLY=1`
/// is also set, exits the process after writing.
pub type InvPredictRecord = (String, String, bool, f64, f64, f64);

pub fn dump_stage6_inv_predict(records: &[InvPredictRecord]) {
    if !is_dump_enabled("OSPREY_DUMP_INV_PREDICT") {
        return;
    }

    let path = "rust_stage6_inv_predict.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "file_name\tis_decoy\tmodified_sequence\tapex_rt\tlibrary_rt\tweight"
    )
    .ok();

    let mut sorted: Vec<&(String, String, bool, f64, f64, f64)> = records.iter().collect();
    // Tuple is (file_name, modseq, is_decoy, apex_rt, library_rt, weight).
    // Sort by (is_decoy, modseq, file_name) and tiebreak on (apex_rt, library_rt)
    // because a peptide with multiple charge-state detections in one file
    // produces multiple rows that tie on the (is_decoy, modseq, file_name) key.
    sorted.sort_by(|a, b| {
        a.2.cmp(&b.2)
            .then(a.1.cmp(&b.1))
            .then(a.0.cmp(&b.0))
            .then(a.3.total_cmp(&b.3))
            .then(a.4.total_cmp(&b.4))
    });

    for (file_name, modseq, is_decoy, apex_rt, library_rt, weight) in &sorted {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            file_name,
            if *is_decoy { "true" } else { "false" },
            modseq,
            format_f64_roundtrip(*apex_rt),
            format_f64_roundtrip(*library_rt),
            format_f64_roundtrip(*weight),
        )
        .ok();
    }
    log::info!(
        "Wrote Stage 6 inverse-predict dump: {} ({} rows)",
        path,
        sorted.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_INV_PREDICT_ONLY", "Stage 6 inverse-predict dump");
}

/// Dump the per-peptide first-pass protein FDR state to
/// `rust_stage6_protein_fdr.tsv`. Mirrors
/// `OspreyDiagnostics.WriteStage6ProteinFdrDump`.
///
/// One row per peptide that appears in `best_scores` (the union of
/// target + decoy modified_sequences seen across all per-file
/// FdrEntry stubs at the moment first-pass protein FDR runs).
/// Columns: `is_decoy, modified_sequence, best_qvalue, score,
/// protein_qvalue`. `best_qvalue` is the input gate (peptide-level
/// run q-value, min across files). `score` is the input ranking
/// (max SVM discriminant across files). `protein_qvalue` is the
/// propagated output -- the value `propagate_protein_qvalues` will
/// write to `FdrEntry.run_protein_qvalue` (1.0 if the peptide is
/// not in `protein_fdr.peptide_qvalues`, matching the
/// `propagate_protein_qvalues` default). Rows sorted by
/// `(is_decoy, modified_sequence)` for stable diff.
///
/// Gated by `OSPREY_DUMP_PROTEIN_FDR=1`. When `OSPREY_PROTEIN_FDR_ONLY=1`
/// is also set, exits the process after writing.
pub fn dump_stage6_protein_fdr(
    best_scores: &HashMap<Arc<str>, PeptideScore>,
    peptide_qvalues: &HashMap<String, f64>,
) {
    if !is_dump_enabled("OSPREY_DUMP_PROTEIN_FDR") {
        return;
    }

    let path = "rust_stage6_protein_fdr.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(
        f,
        "is_decoy\tmodified_sequence\tbest_qvalue\tscore\tprotein_qvalue"
    )
    .ok();

    let mut rows: Vec<(&Arc<str>, &PeptideScore)> = best_scores.iter().collect();
    rows.sort_by(|a, b| {
        a.1.is_decoy
            .cmp(&b.1.is_decoy)
            .then_with(|| a.0.as_ref().cmp(b.0.as_ref()))
    });

    for (modseq, ps) in &rows {
        let q = peptide_qvalues.get(modseq.as_ref()).copied().unwrap_or(1.0);
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            if ps.is_decoy { "true" } else { "false" },
            modseq,
            format_f64_roundtrip(ps.best_qvalue),
            format_f64_roundtrip(ps.score),
            format_f64_roundtrip(q),
        )
        .ok();
    }
    log::info!(
        "Wrote Stage 6 first-pass protein FDR dump: {} ({} rows)",
        path,
        rows.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_PROTEIN_FDR_ONLY", "Stage 6 protein FDR dump");
}

/// Dump the sorted set of detected target peptides handed to
/// `build_protein_parsimony` for Stage 7. Mirrors the OspreySharp
/// `WriteStage7DetectedPeptidesDump` so the protein-FDR input set can be
/// diffed cross-impl before debugging downstream divergence.
///
/// Gated by `OSPREY_DUMP_DETECTED_PEPTIDES=1`. Streams via `BufWriter` so
/// peak memory stays bounded on large datasets (the sorted Vec of
/// references is the only working set beyond the input HashSet).
pub fn dump_stage7_detected_peptides(detected_peptides: &std::collections::HashSet<String>) {
    if !is_dump_enabled("OSPREY_DUMP_DETECTED_PEPTIDES") {
        return;
    }
    let mut sorted: Vec<&String> = detected_peptides.iter().collect();
    sorted.sort();
    let path = "rust_stage7_detected_peptides.txt";
    let Ok(file) = std::fs::File::create(path) else {
        log::warn!("Could not create {}", path);
        return;
    };
    let mut f = std::io::BufWriter::with_capacity(8 << 20, file);
    for s in &sorted {
        if writeln!(f, "{}", s).is_err() {
            log::warn!("Failed writing {}", path);
            return;
        }
    }
    let _ = f.flush();
    log::info!("[DIAG] Wrote {} ({} entries)", path, sorted.len());
}

/// Dump per-protein-group state at the end of second-pass picked-protein
/// FDR (Stage 7 authoritative protein FDR) to `rust_stage7_protein_fdr.tsv`.
/// Mirrors what the OspreySharp port of `compute_protein_fdr` emits.
///
/// One row per protein group present in the parsimony result. The
/// `is_target_winner` column captures the pairwise pick outcome:
/// `true` if the target side scored at or above the decoy side for this
/// group (or the group has only a target side); `false` otherwise. Groups
/// with no winner -- both sides absent because no peptide passed the
/// gate -- are emitted with `group_qvalue = 1.0`, `best_peptide_score = NaN`,
/// `is_target_winner = false`.
///
/// Sort order: `(is_target_winner DESC, group_qvalue ASC, accessions ASC)`.
/// Targets-first keeps the calibration-relevant rows (target winners) at
/// the top of the file; secondary keys make the diff stable. The numeric
/// `group_id` is intentionally NOT a sort key (or a column) -- see below.
///
/// Columns: `accessions`, `n_unique`, `n_shared`, `best_peptide_score`,
/// `group_qvalue`, `is_target_winner`. The numeric `group_id` is omitted
/// because `build_protein_parsimony` assigns it in HashMap iteration order,
/// which would inject per-run noise into a cross-impl bisection. Joining on
/// `accessions` keeps the dump stable.
///
/// Gated by `OSPREY_DUMP_STAGE7_PROTEIN_FDR=1`. When
/// `OSPREY_STAGE7_PROTEIN_FDR_ONLY=1` is also set, exits the process after
/// writing for fast bisection cycle.
pub fn dump_stage7_protein_fdr(parsimony: &ProteinParsimonyResult, fdr_result: &ProteinFdrResult) {
    if !is_dump_enabled("OSPREY_DUMP_STAGE7_PROTEIN_FDR") {
        return;
    }

    let path = "rust_stage7_protein_fdr.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    // `accessions` is the canonical protein-group identifier in the dump;
    // `build_protein_parsimony` assigns numeric `group_id` in HashMap
    // iteration order, which would inject per-run noise into a cross-impl
    // bisection. Joining on accessions instead keeps the dump stable.
    writeln!(
        f,
        "accessions\tn_unique\tn_shared\tbest_peptide_score\tgroup_qvalue\tis_target_winner"
    )
    .ok();

    struct Row {
        accessions: String,
        n_unique: usize,
        n_shared: usize,
        best_peptide_score: f64,
        group_qvalue: f64,
        is_target_winner: bool,
    }

    let mut rows: Vec<Row> = parsimony
        .groups
        .iter()
        .map(|g| {
            // group_qvalues / group_scores contain only TARGET winners
            // (per ProteinFdrResult doc). Presence in either map ⇔ target
            // won the pair for this group.
            let is_target_winner = fdr_result.group_qvalues.contains_key(&g.id);
            let group_qvalue = fdr_result.group_qvalues.get(&g.id).copied().unwrap_or(1.0);
            let best_peptide_score = fdr_result
                .group_scores
                .get(&g.id)
                .copied()
                .unwrap_or(f64::NAN);

            // Stable accession ordering: parsimony.accessions is built by
            // iterating a HashMap, so sort here for deterministic dump.
            let mut accs = g.accessions.clone();
            accs.sort();
            let accessions = accs.join(";");

            Row {
                accessions,
                n_unique: g.unique_peptides.len(),
                n_shared: g.shared_peptides.len(),
                best_peptide_score,
                group_qvalue,
                is_target_winner,
            }
        })
        .collect();

    rows.sort_by(|a, b| {
        // Target winners first (DESC on bool ≡ true before false). Final
        // tiebreak is the sorted-accessions string, NOT group_id, because
        // `build_protein_parsimony` assigns group IDs in HashMap iteration
        // order (random per run); accession lists are stable across runs.
        b.is_target_winner
            .cmp(&a.is_target_winner)
            .then(a.group_qvalue.total_cmp(&b.group_qvalue))
            .then_with(|| a.accessions.cmp(&b.accessions))
    });

    for r in &rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            r.accessions,
            r.n_unique,
            r.n_shared,
            format_f64_roundtrip(r.best_peptide_score),
            format_f64_roundtrip(r.group_qvalue),
            if r.is_target_winner { "true" } else { "false" },
        )
        .ok();
    }

    log::info!(
        "Wrote Stage 7 second-pass protein FDR dump: {} ({} rows)",
        path,
        rows.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_STAGE7_PROTEIN_FDR_ONLY", "Stage 7 protein FDR dump");
}

/// Dump the per-point LOESS fit state of every refit RTCalibration to
/// `rust_stage6_loess_fit.tsv`. Mirrors
/// `OspreyDiagnostics.WriteStage6LoessFitDump`.
///
/// One row per `(file_name, idx)` into the refit's `library_rts` +
/// `fitted_values` + `abs_residuals` arrays. Rows sorted by
/// `(file_name, idx)` for stable diff. The refit dump captures scalar
/// stats (R²/SD/MAD); this dump captures the LOESS curve itself so a
/// stats-vs-smoother bisection is possible: if `(library_rt,
/// fitted_value, abs_residual)` match cross-impl, the divergence is
/// in the stats computation. If `fitted_value` diverges, the LOESS
/// smoother arithmetic itself differs.
///
/// Gated by `OSPREY_DUMP_LOESS_FIT=1`. When `OSPREY_LOESS_FIT_ONLY=1`
/// is also set, exits the process after writing.
pub fn dump_stage6_loess_fit(
    refined_calibrations: &HashMap<String, osprey_chromatography::RTCalibration>,
) {
    if !is_dump_enabled("OSPREY_DUMP_LOESS_FIT") {
        return;
    }

    let path = "rust_stage6_loess_fit.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(f, "file_name\tidx\tlibrary_rt\tfitted_value\tabs_residual").ok();

    let mut file_names: Vec<&String> = refined_calibrations.keys().collect();
    file_names.sort();

    let mut total_rows = 0usize;
    for file_name in &file_names {
        let cal = &refined_calibrations[file_name.as_str()];
        let params = cal.export_model_params();
        let n = params
            .library_rts
            .len()
            .min(params.fitted_rts.len())
            .min(params.abs_residuals.len());
        for i in 0..n {
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}",
                file_name,
                i,
                format_f64_roundtrip(params.library_rts[i]),
                format_f64_roundtrip(params.fitted_rts[i]),
                format_f64_roundtrip(params.abs_residuals[i]),
            )
            .ok();
        }
        total_rows += n;
    }
    log::info!(
        "Wrote Stage 6 LOESS fit dump: {} ({} rows across {} files)",
        path,
        total_rows,
        file_names.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_LOESS_FIT_ONLY", "Stage 6 LOESS fit dump");
}

/// Dump the per-file refined-calibration statistics produced by
/// `refit_calibration_with_consensus` to `rust_stage6_refit.tsv`.
/// Mirrors `OspreySharp.FDR.Reconciliation.CalibrationRefit.Refit`.
///
/// Gated by `OSPREY_DUMP_REFIT=1`. When `OSPREY_REFIT_ONLY=1` is also
/// set, exits the process after writing. Columns: `file_name, n_points,
/// r_squared, residual_sd, mad`. Rows are emitted only for files where
/// the refit produced a calibration (insufficient-points failures are
/// absent), sorted by `file_name` for stable diff.
pub fn dump_stage6_refit(
    refined_calibrations: &HashMap<String, osprey_chromatography::RTCalibration>,
) {
    if !is_dump_enabled("OSPREY_DUMP_REFIT") {
        return;
    }

    let path = "rust_stage6_refit.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    writeln!(f, "file_name\tn_points\tr_squared\tresidual_sd\tmad").ok();

    let mut rows: Vec<(&str, _)> = refined_calibrations
        .iter()
        .map(|(name, cal)| (name.as_str(), cal.stats()))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));

    for (file_name, stats) in &rows {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            file_name,
            stats.n_points,
            format_f64_roundtrip(stats.r_squared),
            format_f64_roundtrip(stats.residual_std),
            format_f64_roundtrip(stats.mad),
        )
        .ok();
    }
    log::info!("Wrote Stage 6 refit dump: {} ({} rows)", path, rows.len());

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_REFIT_ONLY", "Stage 6 refit dump");
}

/// Dump the per-(file, entry) `ReconcileAction` map produced by
/// `plan_reconciliation` to `rust_stage6_reconciliation.tsv`. Pairs with
/// the C# `WriteStage6ReconciliationDump` for cross-impl byte-parity at
/// the planning step, before per-file rescoring runs.
///
/// Gated by `OSPREY_DUMP_RECONCILIATION=1`. When `OSPREY_RECONCILIATION_ONLY=1`
/// is also set, exits the process after writing. Columns:
/// `file_name, entry_id, action, apex_or_expected_rt, start_rt, end_rt,
/// half_width, candidate_index`. Only non-Keep actions are emitted (matches
/// the C# planner's contract -- Keep entries are absent from the map).
/// Rows sorted by `(file_name, entry_id)` for stable diff.
pub fn dump_stage6_reconciliation(
    actions: &HashMap<(String, usize), crate::reconciliation::ReconcileAction>,
    per_file_entries: &[(String, Vec<osprey_core::FdrEntry>)],
) {
    if !is_dump_enabled("OSPREY_DUMP_RECONCILIATION") {
        return;
    }

    let path = "rust_stage6_reconciliation.tsv";
    // 8 MB BufWriter so each writeln! is a memcpy, not a syscall.
    // Critical on WSL drvfs/9P where unbuffered per-row writes turn
    // a 1-minute dump into 10+ minutes. Diagnostic-only; bytes unchanged.
    let Ok(mut f) =
        std::fs::File::create(path).map(|file| std::io::BufWriter::with_capacity(8 << 20, file))
    else {
        log::warn!("Could not create {}", path);
        return;
    };

    // Map (file_name -> &[FdrEntry]) so we can resolve list-index to entry_id.
    let entries_by_file: HashMap<&str, &[osprey_core::FdrEntry]> = per_file_entries
        .iter()
        .map(|(name, e)| (name.as_str(), e.as_slice()))
        .collect();

    use crate::reconciliation::ReconcileAction;
    let mut rows: Vec<(&str, u32, &ReconcileAction)> = Vec::with_capacity(actions.len());
    for ((file, idx), action) in actions {
        if matches!(action, ReconcileAction::Keep) {
            continue;
        }
        let Some(entries) = entries_by_file.get(file.as_str()) else {
            continue;
        };
        if *idx >= entries.len() {
            continue;
        }
        rows.push((file.as_str(), entries[*idx].entry_id, action));
    }
    rows.sort_by(|a, b| a.0.cmp(b.0).then(a.1.cmp(&b.1)));

    writeln!(
        f,
        "file_name\tentry_id\taction\tapex_or_expected_rt\tstart_rt\tend_rt\thalf_width\tcandidate_index"
    )
    .ok();
    for (file_name, entry_id, action) in &rows {
        match action {
            ReconcileAction::UseCwtPeak {
                candidate_idx,
                start_rt,
                apex_rt,
                end_rt,
            } => {
                writeln!(
                    f,
                    "{}\t{}\tuse_cwt_peak\t{}\t{}\t{}\t\t{}",
                    file_name,
                    entry_id,
                    format_f64_roundtrip(*apex_rt),
                    format_f64_roundtrip(*start_rt),
                    format_f64_roundtrip(*end_rt),
                    candidate_idx,
                )
                .ok();
            }
            ReconcileAction::ForcedIntegration {
                expected_rt,
                half_width,
            } => {
                writeln!(
                    f,
                    "{}\t{}\tforced_integration\t{}\t\t\t{}\t",
                    file_name,
                    entry_id,
                    format_f64_roundtrip(*expected_rt),
                    format_f64_roundtrip(*half_width),
                )
                .ok();
            }
            ReconcileAction::Keep => unreachable!("filtered above"),
        }
    }
    log::info!(
        "Wrote Stage 6 reconciliation dump: {} ({} rows)",
        path,
        rows.len()
    );

    let _ = f.flush();
    drop(f);
    exit_if_only("OSPREY_RECONCILIATION_ONLY", "Stage 6 reconciliation dump");
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::FdrEntry;
    use std::sync::Arc;

    fn make_entry(id: u32, score: f64, is_decoy: bool) -> FdrEntry {
        FdrEntry {
            entry_id: id,
            parquet_index: id,
            is_decoy,
            charge: 2,
            scan_number: 0,
            apex_rt: 0.0,
            start_rt: 0.0,
            end_rt: 0.0,
            coelution_sum: 0.0,
            score,
            run_precursor_qvalue: 0.1,
            run_peptide_qvalue: 0.1,
            run_protein_qvalue: 1.0,
            experiment_precursor_qvalue: 0.1,
            experiment_peptide_qvalue: 0.1,
            experiment_protein_qvalue: 1.0,
            pep: 0.5,
            experiment_aggregate_score: 0.0,
            modified_sequence: Arc::from("PEPTIDE"),
        }
    }

    #[test]
    fn dump_is_noop_without_env_var() {
        // Hermetic: chdir into a fresh temp dir so a stale
        // rust_stage5_percolator.tsv from a previous run (or another
        // test) can't produce a spurious pass. The Rust-standard
        // tempfile crate is a dev-dep of this workspace.
        let tmp = tempfile::tempdir().expect("tempdir");
        let prev_cwd = std::env::current_dir().expect("cwd");
        std::env::set_current_dir(tmp.path()).expect("chdir");

        std::env::remove_var("OSPREY_DUMP_PERCOLATOR");
        let entries = vec![("f1".to_string(), vec![make_entry(0, 1.0, false)])];
        dump_stage5_percolator(&entries);
        let written = std::path::Path::new("rust_stage5_percolator.tsv").exists();

        // Restore cwd before asserting so a failure doesn't strand the
        // process in the temp dir for other tests.
        std::env::set_current_dir(prev_cwd).expect("restore cwd");
        assert!(
            !written,
            "dump wrote a file even though the gate was not set"
        );
    }
}
