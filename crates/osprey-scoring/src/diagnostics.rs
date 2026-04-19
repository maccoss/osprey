//! Cross-implementation bisection diagnostic dumps for the scoring stage.
//!
//! Each function below is gated by an `OSPREY_DUMP_*` env var and writes a
//! `rust_*.txt` file matching the format of the corresponding C# OspreySharp
//! dump (`cs_*.txt`). Dumps that have an `_ONLY` companion env var exit
//! the process after writing for fast cycle-time during bisection.
//!
//! See `osprey-core::diagnostics` for shared primitives and the env-var
//! gating convention.

use crate::batch::CalibrationMatch;
use osprey_core::diagnostics::{exit_if_only, is_dump_enabled};
use osprey_core::LibraryEntry;
use std::collections::HashMap;
use std::io::Write;

/// Dump per-entry calibration match info to `rust_cal_match.txt`.
///
/// Writes a TSV row per library entry (target or decoy) with match status,
/// apex scan number when matched, and the four LDA features (correlation,
/// libcosine, top6, xcorr) plus signal-to-noise. Sorted by entry_id for a
/// stable diff against C#'s `cs_cal_match.txt`. The `scan` column is a
/// Rust-only extension: C# leaves it blank, Rust populates it from
/// `CalibrationMatch::scan_number`.
///
/// Gated by `OSPREY_DUMP_CAL_MATCH=1`. When `OSPREY_CAL_MATCH_ONLY=1` is
/// also set, the process exits after writing for fast cycle-time during
/// calibration bisection.
pub fn dump_cal_match(library: &[LibraryEntry], results: &[CalibrationMatch]) {
    if !is_dump_enabled("OSPREY_DUMP_CAL_MATCH") {
        return;
    }

    let dump_path = "rust_cal_match.txt";

    // Build match lookup keyed by entry_id
    let by_id: HashMap<u32, &CalibrationMatch> = results.iter().map(|m| (m.entry_id, m)).collect();

    if let Ok(mut f) = std::fs::File::create(dump_path) {
        // 11 columns. snr is the signal-to-noise the S/N filter gates on
        // (matches entering LOESS must have snr >= MIN_SNR_FOR_RT_CAL=5.0).
        writeln!(
            f,
            "entry_id\tis_decoy\tcharge\thas_match\tscan\tapex_rt\tcorrelation\tlibcosine\ttop6\txcorr\tsnr"
        )
        .ok();

        let mut entries: Vec<&LibraryEntry> = library.iter().collect();
        entries.sort_by_key(|e| e.id);

        let mut n_matched = 0usize;
        let mut n_unmatched = 0usize;

        for entry in entries {
            if let Some(m) = by_id.get(&entry.id) {
                // Use :.10 everywhere so we don't hit banker's vs round-
                // half-up rounding differences between Rust and C#.
                writeln!(
                    f,
                    "{}\t{}\t{}\t1\t{}\t{:.10}\t{:.10}\t{:.10}\t{}\t{:.10}\t{:.10}",
                    entry.id,
                    if entry.is_decoy { 1 } else { 0 },
                    entry.charge,
                    m.scan_number,
                    m.measured_rt,
                    m.correlation_score,
                    m.libcosine_apex,
                    m.top6_matched_apex,
                    m.xcorr_score,
                    m.signal_to_noise
                )
                .ok();
                n_matched += 1;
            } else {
                writeln!(
                    f,
                    "{}\t{}\t{}\t0\t\t\t\t\t\t\t",
                    entry.id,
                    if entry.is_decoy { 1 } else { 0 },
                    entry.charge
                )
                .ok();
                n_unmatched += 1;
            }
        }
        log::info!(
            "Wrote calibration match dump: {} ({} matched, {} unmatched)",
            dump_path,
            n_matched,
            n_unmatched
        );
    }

    exit_if_only("OSPREY_CAL_MATCH_ONLY", "match dump");
}

/// Dump per-entry LDA discriminant + q-value to `rust_lda_scores.txt`,
/// sorted by entry_id for a stable diff against `cs_lda_scores.txt`.
///
/// Uses `{:.10}` formatting to avoid banker's-vs-half-up text rounding
/// mismatches with C#.
///
/// Gated by `OSPREY_DUMP_LDA_SCORES=1`. When `OSPREY_LDA_SCORES_ONLY=1`
/// is also set, exits the process after writing.
pub fn dump_lda_scores(matches: &[CalibrationMatch]) {
    if !is_dump_enabled("OSPREY_DUMP_LDA_SCORES") {
        return;
    }

    if let Ok(mut f) = std::fs::File::create("rust_lda_scores.txt") {
        writeln!(f, "entry_id\tis_decoy\tdiscriminant\tq_value").ok();
        let mut indices: Vec<usize> = (0..matches.len()).collect();
        indices.sort_by_key(|&i| matches[i].entry_id);
        for i in indices {
            let m = &matches[i];
            writeln!(
                f,
                "{}\t{}\t{:.10}\t{:.10}",
                m.entry_id,
                if m.is_decoy { 1 } else { 0 },
                m.discriminant_score,
                m.q_value
            )
            .ok();
        }
        log::info!(
            "Wrote LDA scores dump: rust_lda_scores.txt ({} entries)",
            matches.len()
        );
    }

    exit_if_only("OSPREY_LDA_SCORES_ONLY", "LDA dump");
}

/// Dump the calibration sampling grid setup to `rust_cal_scalars.txt` and
/// `rust_cal_grid.txt`. Writes the scalar parameters (target/decoy counts,
/// RT and m/z bounds and bin widths, occupied-cell count, per-cell quota,
/// seed) and the per-cell content.
///
/// Called *before* the deterministic stride sampling runs; pair with
/// [`dump_cal_sample_entries`] which writes the sampled entries afterward.
/// Gated by `OSPREY_DUMP_CAL_SAMPLE=1`.
#[allow(clippy::too_many_arguments)]
pub fn dump_cal_sample_setup(
    targets: &[&LibraryEntry],
    decoys: &[&LibraryEntry],
    bins_per_axis: usize,
    rt_min: f64,
    rt_max: f64,
    mz_min: f64,
    mz_max: f64,
    rt_range: f64,
    mz_range: f64,
    rt_bin_width: f64,
    mz_bin_width: f64,
    n_occupied: usize,
    per_cell: usize,
    seed: u64,
    grid: &[Vec<Vec<usize>>],
) {
    if !is_dump_enabled("OSPREY_DUMP_CAL_SAMPLE") {
        return;
    }

    if let Ok(mut f) = std::fs::File::create("rust_cal_scalars.txt") {
        writeln!(f, "n_targets\t{}", targets.len()).ok();
        writeln!(f, "n_decoys\t{}", decoys.len()).ok();
        writeln!(f, "bins_per_axis\t{}", bins_per_axis).ok();
        writeln!(f, "rt_min\t{:.17}", rt_min).ok();
        writeln!(f, "rt_max\t{:.17}", rt_max).ok();
        writeln!(f, "mz_min\t{:.17}", mz_min).ok();
        writeln!(f, "mz_max\t{:.17}", mz_max).ok();
        writeln!(f, "rt_range\t{:.17}", rt_range).ok();
        writeln!(f, "mz_range\t{:.17}", mz_range).ok();
        writeln!(f, "rt_bin_width\t{:.17}", rt_bin_width).ok();
        writeln!(f, "mz_bin_width\t{:.17}", mz_bin_width).ok();
        writeln!(f, "n_occupied\t{}", n_occupied).ok();
        writeln!(f, "per_cell\t{}", per_cell).ok();
        writeln!(f, "seed\t{}", seed).ok();
    }
    // Dump full grid: rt_bin, mz_bin, count, target_ids (sorted ascending
    // within cell, to make comparison order-independent).
    if let Ok(mut f) = std::fs::File::create("rust_cal_grid.txt") {
        writeln!(f, "rt_bin\tmz_bin\tcount\ttarget_ids").ok();
        for (r, row) in grid.iter().enumerate().take(bins_per_axis) {
            for (c, cell) in row.iter().enumerate().take(bins_per_axis) {
                if cell.is_empty() {
                    continue;
                }
                let mut ids: Vec<u32> = cell.iter().map(|&i| targets[i].id).collect();
                ids.sort();
                let ids_str: Vec<String> = ids.iter().map(|i| i.to_string()).collect();
                writeln!(f, "{}\t{}\t{}\t{}", r, c, cell.len(), ids_str.join(",")).ok();
            }
        }
    }
}

/// Dump the post-sampling target entries to `rust_cal_sample.txt`,
/// sorted by stringified row for stable diff against `cs_cal_sample.txt`.
///
/// Pair with [`dump_cal_sample_setup`] which dumps the grid before sampling.
/// Gated by `OSPREY_DUMP_CAL_SAMPLE=1`. When `OSPREY_CAL_SAMPLE_ONLY=1` is
/// also set, exits the process after writing.
pub fn dump_cal_sample_entries(sampled: &[LibraryEntry]) {
    if !is_dump_enabled("OSPREY_DUMP_CAL_SAMPLE") {
        return;
    }

    let dump_path = "rust_cal_sample.txt";
    if let Ok(mut f) = std::fs::File::create(dump_path) {
        let mut tuples: Vec<String> = sampled
            .iter()
            .filter(|e| !e.is_decoy)
            .map(|e| {
                format!(
                    "{}\t{}\t{}\t{:.4}\t{:.4}",
                    e.id, e.modified_sequence, e.charge, e.precursor_mz, e.retention_time
                )
            })
            .collect();
        tuples.sort();
        writeln!(f, "id\tmodseq\tcharge\tmz\trt").ok();
        for t in &tuples {
            writeln!(f, "{}", t).ok();
        }
        log::info!(
            "Wrote calibration sample to {} ({} targets)",
            dump_path,
            tuples.len()
        );
    }

    exit_if_only("OSPREY_CAL_SAMPLE_ONLY", "sample dump");
}

// --------------------------------------------------------------------------
// Calibration window/prefilter dumps share a Mutex<Vec<String>> collector
// pattern: rows are produced inside a parallel scoring loop and written
// (sorted, with header) once scoring completes.
// --------------------------------------------------------------------------

/// Allocate a row collector for the per-(entry, window) calibration window
/// dump, gated by `OSPREY_DUMP_CAL_WINDOWS=1`. Returns `None` if the env
/// var isn't set so the parallel loop can skip row formatting entirely.
pub fn dump_cal_windows_init(library_len: usize) -> Option<std::sync::Mutex<Vec<String>>> {
    if !is_dump_enabled("OSPREY_DUMP_CAL_WINDOWS") {
        return None;
    }
    Some(std::sync::Mutex::new(Vec::with_capacity(library_len * 2)))
}

/// Format a single row for the calibration window dump. Called from within
/// the parallel scoring loop where the local variables already exist.
pub fn dump_cal_windows_row(
    entry: &LibraryEntry,
    iso_lower: f64,
    iso_upper: f64,
    expected_rt: f64,
    rt_window_start: f64,
    rt_window_end: f64,
) -> String {
    format!(
        "{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
        entry.id,
        if entry.is_decoy { 1 } else { 0 },
        entry.charge,
        entry.precursor_mz,
        entry.retention_time,
        iso_lower,
        iso_upper,
        expected_rt,
        rt_window_start,
        rt_window_end,
    )
}

/// Write the collected rows to `rust_cal_windows.txt` (sorted) once
/// scoring completes. Exits the process if `OSPREY_CAL_WINDOWS_ONLY=1`.
pub fn dump_cal_windows_finalize(window_dump: Option<std::sync::Mutex<Vec<String>>>) {
    let Some(mtx) = window_dump else {
        return;
    };
    let mut rows = mtx.into_inner().unwrap_or_default();
    rows.sort();
    if let Ok(mut f) = std::fs::File::create("rust_cal_windows.txt") {
        writeln!(
            f,
            "entry_id\tis_decoy\tcharge\tprecursor_mz\tlibrary_rt\tiso_lower\tiso_upper\texpected_rt\trt_window_start\trt_window_end"
        )
        .ok();
        for r in &rows {
            writeln!(f, "{}", r).ok();
        }
        log::info!(
            "Wrote calibration windows dump: rust_cal_windows.txt ({} rows)",
            rows.len()
        );
    }
    exit_if_only("OSPREY_CAL_WINDOWS_ONLY", "window dump");
}

/// Allocate a row collector for the per-(entry, window) calibration
/// prefilter dump, gated by `OSPREY_DUMP_CAL_PREFILTER=1`.
///
/// Rust-only diagnostic; OspreySharp does not implement this dump.
pub fn dump_cal_prefilter_init(library_len: usize) -> Option<std::sync::Mutex<Vec<String>>> {
    if !is_dump_enabled("OSPREY_DUMP_CAL_PREFILTER") {
        return None;
    }
    Some(std::sync::Mutex::new(Vec::with_capacity(library_len * 2)))
}

/// Format a single row for the calibration prefilter dump.
pub fn dump_cal_prefilter_row(entry: &LibraryEntry, candidate_scans: &[u32]) -> String {
    let scans_str = candidate_scans
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{}\t{}\t{}\t{:.6}\t{}\t{}",
        entry.id,
        if entry.is_decoy { 1 } else { 0 },
        entry.charge,
        entry.precursor_mz,
        candidate_scans.len(),
        scans_str
    )
}

/// Write the collected rows to `rust_cal_prefilter.txt` (sorted). Exits
/// the process if `OSPREY_CAL_PREFILTER_ONLY=1`.
pub fn dump_cal_prefilter_finalize(prefilter_dump: Option<std::sync::Mutex<Vec<String>>>) {
    let Some(mtx) = prefilter_dump else {
        return;
    };
    let mut rows = mtx.into_inner().unwrap_or_default();
    rows.sort();
    if let Ok(mut f) = std::fs::File::create("rust_cal_prefilter.txt") {
        writeln!(
            f,
            "entry_id\tis_decoy\tcharge\tprecursor_mz\tn_candidates\tscan_numbers"
        )
        .ok();
        for r in &rows {
            writeln!(f, "{}", r).ok();
        }
        log::info!(
            "Wrote calibration prefilter dump: rust_cal_prefilter.txt ({} rows)",
            rows.len()
        );
    }
    exit_if_only("OSPREY_CAL_PREFILTER_ONLY", "prefilter dump");
}

// --------------------------------------------------------------------------
// Per-entry calibration XIC dump
// --------------------------------------------------------------------------

/// State for the per-entry calibration XIC diagnostic dump.
///
/// Reads env vars once at construction and checks them on every entry
/// inside the parallel scoring loop without per-entry env-var lookups.
/// Used by `OSPREY_DIAG_XIC_ENTRY_ID` (which entry id to dump in detail)
/// + `OSPREY_DIAG_XIC_PASS` (1 or 2; default 1).
///
/// On first match the dump writes `rust_xic_entry_<id>.txt` (header +
/// LOESS stats + PASS calculations + CANDIDATES + TOP-6 fragments via
/// [`Self::dump_header`], then EXTRACTED XICS via
/// [`Self::dump_xics_and_exit`]) and exits the process so the cycle
/// time is bounded by the time to first match (not full pipeline).
pub struct CalXicEntryDump {
    target_id: Option<u32>,
    current_pass: u32,
    active_for_pass: bool,
}

impl CalXicEntryDump {
    /// Create a new dump state for the given calibration pass (1 or 2).
    /// Reads `OSPREY_DIAG_XIC_ENTRY_ID` and `OSPREY_DIAG_XIC_PASS` once.
    pub fn new(current_pass: u32) -> Self {
        let target_id: Option<u32> = std::env::var("OSPREY_DIAG_XIC_ENTRY_ID")
            .ok()
            .and_then(|s| s.parse::<u32>().ok());
        let target_pass: u32 = std::env::var("OSPREY_DIAG_XIC_PASS")
            .ok()
            .and_then(|s| s.parse::<u32>().ok())
            .unwrap_or(1);
        Self {
            target_id,
            current_pass,
            active_for_pass: current_pass == target_pass,
        }
    }

    /// Cheap per-entry check: returns true only if this entry's id matches
    /// `OSPREY_DIAG_XIC_ENTRY_ID` AND the configured pass matches.
    pub fn is_active_for(&self, entry_id: u32) -> bool {
        self.active_for_pass && self.target_id == Some(entry_id)
    }

    /// Write the header section (LOESS stats, PASS calculations,
    /// CANDIDATES, TOP-6 fragments) for this entry. No-op if the entry
    /// is not the diagnostic target.
    pub fn dump_header(
        &self,
        entry: &LibraryEntry,
        expected_rt: f64,
        rt_tolerance: f64,
        candidate_spectra: &[&osprey_core::Spectrum],
        rt_calibration_for_diag: Option<&osprey_chromatography::RTCalibration>,
    ) {
        if !self.is_active_for(entry.id) {
            return;
        }
        let dump_path = format!("rust_xic_entry_{}.txt", entry.id);
        let Ok(mut f) = std::fs::File::create(&dump_path) else {
            return;
        };
        writeln!(
            f,
            "# per-entry chromatogram dump for entry_id={} (pass {})",
            entry.id, self.current_pass
        )
        .ok();
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

        // Pass 2 calculations block: dump the inputs that feed into XIC
        // extraction so if the XICs don't match, we already have the
        // intermediate values to localize the divergence (LOESS model
        // stats, predicted RT, refined tolerance, selection window).
        if let Some(rt_cal) = rt_calibration_for_diag {
            let loess_stats = rt_cal.stats();
            writeln!(f, "# LOESS MODEL (pass 2 RT calibration)").ok();
            writeln!(f, "# loess.n_points={}", loess_stats.n_points).ok();
            writeln!(f, "# loess.r_squared={:.10}", loess_stats.r_squared).ok();
            writeln!(f, "# loess.residual_sd={:.10}", loess_stats.residual_std).ok();
            writeln!(f, "# loess.mean_residual={:.10}", loess_stats.mean_residual).ok();
            writeln!(f, "# loess.max_residual={:.10}", loess_stats.max_residual).ok();
            writeln!(
                f,
                "# loess.p20_abs_residual={:.10}",
                loess_stats.p20_abs_residual
            )
            .ok();
            writeln!(
                f,
                "# loess.p80_abs_residual={:.10}",
                loess_stats.p80_abs_residual
            )
            .ok();
            writeln!(f, "# loess.mad={:.10}", loess_stats.mad).ok();
        }

        writeln!(f, "# PASS CALCULATIONS").ok();
        writeln!(f, "# pass.library_rt={:.10}", entry.retention_time).ok();
        writeln!(f, "# pass.expected_rt={:.10}", expected_rt).ok();
        writeln!(f, "# pass.tolerance={:.10}", rt_tolerance).ok();
        writeln!(f, "# pass.rt_window_lo={:.10}", expected_rt - rt_tolerance).ok();
        writeln!(f, "# pass.rt_window_hi={:.10}", expected_rt + rt_tolerance).ok();
        // Rust's pipeline uses an identity linear mapping (slope=1,
        // intercept=0) when ranges are similar. We emit these so the
        // diff column-for-column matches the C# header layout.
        writeln!(f, "# pass.rt_slope={:.10}", 1.0_f64).ok();
        writeln!(f, "# pass.rt_intercept={:.10}", 0.0_f64).ok();

        writeln!(
            f,
            "# n_post_prefilter_candidates={}",
            candidate_spectra.len()
        )
        .ok();
        writeln!(f, "# CANDIDATES (post-prefilter, sorted by RT)").ok();
        writeln!(
            f,
            "candidate\tscan_idx\tscan_number\trt\tiso_lower\tiso_upper"
        )
        .ok();
        for (i, s) in candidate_spectra.iter().enumerate() {
            writeln!(
                f,
                "candidate\t{}\t{}\t{:.10}\t{:.10}\t{:.10}",
                i,
                s.scan_number,
                s.retention_time,
                s.isolation_window.lower_bound(),
                s.isolation_window.upper_bound()
            )
            .ok();
        }
        // Also list top-6 fragments with their m/z so the diff is interpretable.
        let top_indices = crate::get_top_n_fragment_indices(&entry.fragments, 6);
        writeln!(f, "# TOP-6 FRAGMENTS (selected by intensity desc)").ok();
        writeln!(f, "topfrag\ttop_idx\tlib_idx\tlib_mz\tlib_intensity").ok();
        for (rank, &fi) in top_indices.iter().enumerate() {
            let f_obj = &entry.fragments[fi];
            writeln!(
                f,
                "topfrag\t{}\t{}\t{:.10}\t{:.10}",
                rank, fi, f_obj.mz, f_obj.relative_intensity
            )
            .ok();
        }
    }

    /// Append the EXTRACTED XICS section, then exit the process. Called
    /// after [`Self::dump_header`] once the per-fragment XICs have been
    /// extracted. No-op if the entry is not the diagnostic target.
    pub fn dump_xics_and_exit(&self, entry: &LibraryEntry, xics: &[(usize, Vec<(f64, f64)>)]) {
        if !self.is_active_for(entry.id) {
            return;
        }
        let dump_path = format!("rust_xic_entry_{}.txt", entry.id);
        if let Ok(mut f) = std::fs::OpenOptions::new().append(true).open(&dump_path) {
            writeln!(f, "# EXTRACTED XICS (lib_idx, scan_idx, rt, intensity)").ok();
            writeln!(f, "xic\tlib_idx\tscan_idx\trt\tintensity").ok();
            // Use 10 decimal places — wide enough that half-way rounding
            // on values like 1756.6640625 (exact f32 mantissa) doesn't
            // produce text diffs from banker's-vs-half-up rounding
            // mode differences between Rust and C#.
            for (lib_idx, xic) in xics {
                for (i, (rt, intensity)) in xic.iter().enumerate() {
                    writeln!(f, "xic\t{}\t{}\t{:.10}\t{:.10}", lib_idx, i, rt, intensity).ok();
                }
            }
        }
        log::info!(
            "[BISECT] OSPREY_DIAG_XIC_ENTRY_ID matched on pass {} - wrote {} and exiting",
            self.current_pass,
            dump_path
        );
        std::process::exit(0);
    }
}
