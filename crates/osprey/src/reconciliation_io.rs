//! Per-file `<stem>.reconciliation.json` Stage 5 → Stage 6 boundary file.
//!
//! Produced at the end of Stage 5 (first-pass FDR + multi-charge
//! consensus + cross-run consensus RTs + per-file calibration refit +
//! reconciliation planning) for each input file. Carries everything a
//! Stage 6 worker needs to do per-file rescore + gap-fill + reconciled
//! parquet write-back without re-running any of the joined work: the
//! non-Keep `ReconcileAction`s for the file's entries (split into two
//! homogeneous arrays — `use_cwt_peak_actions` and
//! `forced_integration_actions` — so the JSON has no discriminator
//! field gymnastics), the `GapFillTarget`s for the file (precursors
//! that passed FDR in a sibling replicate but missed in this one), and
//! the refined RT calibration (LOESS refit on consensus peptides).
//!
//! Field order is alphabetical at every nesting level — both at the
//! top-level envelope and within each entry — so the JSON output matches
//! the C# emitter byte-for-byte. Numeric fields are routed through
//! `osprey_core::diagnostics::format_f64_roundtrip` on both sides
//! (rather than each language's default shortest-roundtrip path), giving
//! a single canonical fixed-point decimal form (no scientific notation)
//! that is byte-identical across runtimes for every f64 input. The
//! Newtonsoft `R`/Grisu and Rust `ryu` formatters disagree on the
//! decimal-vs-scientific threshold for small values, so neither default
//! is suitable for cross-impl byte parity.
//!
//! Mirrors `OspreySharp.IO.ReconciliationFile` on the C# side. Cross-impl
//! byte parity is verified by a sibling test in both languages.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use anyhow::{anyhow, Result};
use osprey_chromatography::RTCalibration;
use osprey_core::diagnostics::format_f64_roundtrip;
use osprey_core::FdrEntry;
use serde::{Deserialize, Serialize};
use serde_json::ser::{Formatter, PrettyFormatter};

use crate::reconciliation::{GapFillTarget, ReconcileAction};

/// JSON formatter that delegates layout to `serde_json`'s
/// `PrettyFormatter` but routes every f64 through
/// `format_f64_roundtrip`. Together with the matching
/// `RoundtripJsonConverter` on the C# side this guarantees that every
/// f64 in the reconciliation envelope is emitted as the same shortest-
/// roundtrip fixed-point decimal on both runtimes — sidestepping the
/// ryu-vs-Grisu threshold disagreement on values like `4.58e-5` that
/// would otherwise silently break cross-impl byte parity.
struct RoundtripPrettyFormatter<'a> {
    inner: PrettyFormatter<'a>,
}

impl<'a> RoundtripPrettyFormatter<'a> {
    fn new() -> Self {
        Self {
            inner: PrettyFormatter::new(),
        }
    }
}

impl<'a> Formatter for RoundtripPrettyFormatter<'a> {
    fn write_f64<W>(&mut self, writer: &mut W, value: f64) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        // format_f64_roundtrip returns "NaN" / "inf" / "-inf" for non-
        // finite values which are not valid JSON numbers; reconciliation
        // arrays should never carry those, but guard anyway by failing
        // the write rather than emitting a malformed token.
        if !value.is_finite() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("non-finite f64 in reconciliation JSON: {}", value),
            ));
        }
        writer.write_all(format_f64_roundtrip(value).as_bytes())
    }

    // Delegate the indent-aware methods to the wrapped PrettyFormatter
    // so 2-space indentation and key/value spacing are preserved.
    fn begin_array<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.begin_array(writer)
    }
    fn end_array<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.end_array(writer)
    }
    fn begin_array_value<W>(&mut self, writer: &mut W, first: bool) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.begin_array_value(writer, first)
    }
    fn end_array_value<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.end_array_value(writer)
    }
    fn begin_object<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.begin_object(writer)
    }
    fn end_object<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.end_object(writer)
    }
    fn begin_object_key<W>(&mut self, writer: &mut W, first: bool) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.begin_object_key(writer, first)
    }
    fn begin_object_value<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        self.inner.begin_object_value(writer)
    }
    fn end_object_value<W>(&mut self, writer: &mut W) -> std::io::Result<()>
    where
        W: ?Sized + std::io::Write,
    {
        // PrettyFormatter tracks `has_value` here to decide whether
        // end_object should emit a newline before the closing brace; if
        // we don't delegate, every object collapses to `<value>}` with
        // no trailing newline.
        self.inner.end_object_value(writer)
    }
}

/// Top-level JSON envelope. Field declaration order is alphabetical so
/// `serde_json::to_writer_pretty` emits keys in alphabetical order and
/// matches the C# emitter byte-for-byte.
///
/// `file_stems` (added 2026-05-19, format_version 2) carries the full
/// set of file stems that participated in the first-join phase. The
/// per-file rescore worker (`--join-at-pass=1 --no-join`) consumes only
/// its own parquet, so its `config.input_files` reflects a single
/// stem. To produce a reconciled parquet that the downstream
/// `--join-at-pass=2` merge node can validate (its
/// `reconciliation_parameter_hash` is computed over all parquets it
/// receives), the worker uses this list to compute the same multi-stem
/// hash the original join used. v1 files (no `file_stems` field) are
/// rejected at deserialization — re-run the first-join phase to
/// regenerate a v2 envelope. The pipeline is not in active production
/// use, so backward compat is not maintained; v1 acceptance previously
/// silently fell back to a single-stem hash that the merge node would
/// reject anyway.
///
/// `first_pass_base_ids` (added format_version 3) carries the join-wide
/// set of base_ids that survived first-pass compaction. The first-join
/// phase computes it with every file in memory (the distinct base_ids
/// remaining across all files after compaction); a per-file rescore
/// worker, which only holds its own file, uses it to compact to exactly
/// the set the in-memory straight-through pipeline used, instead of a
/// per-file subset that drops cross-file entries. v2 files (no
/// `first_pass_base_ids` field) are rejected at deserialization for the
/// same no-backward-compat reason as v1. Mirrors `first_pass_base_ids` on
/// the C# `ReconciliationFile` for cross-impl byte parity.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReconciliationFile {
    /// Sorted list of file stems that participated in the first-join
    /// phase. Required (no `#[serde(default)]`): a missing field
    /// triggers a serde deserialization error so legacy v1 envelopes
    /// are rejected at parse time with a clear error rather than
    /// silently falling back to single-stem hashes.
    pub file_stems: Vec<String>,
    /// Join-wide set of first-pass passing base_ids, sorted ascending for
    /// deterministic byte-parity output. Required (no `#[serde(default)]`):
    /// a missing field triggers a serde deserialization error so legacy
    /// v1/v2 envelopes are rejected at parse time rather than a per-file
    /// worker silently recomputing a divergent per-file subset.
    pub first_pass_base_ids: Vec<u32>,
    pub forced_integration_actions: Vec<ForcedIntegrationEntry>,
    pub format_version: u32,
    pub gap_fill_targets: Vec<GapFillEntry>,
    pub library_hash: String,
    pub refined_rt_calibration: Option<RefinedRtCalibrationJson>,
    pub search_hash: String,
    pub use_cwt_peak_actions: Vec<UseCwtPeakEntry>,
}

/// Wire form of a `ReconcileAction::UseCwtPeak`. Field order alphabetical.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UseCwtPeakEntry {
    pub apex_rt: f64,
    pub candidate_idx: u32,
    pub end_rt: f64,
    pub entry_id: u32,
    pub start_rt: f64,
}

/// Wire form of a `ReconcileAction::ForcedIntegration`. Field order alphabetical.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ForcedIntegrationEntry {
    pub entry_id: u32,
    pub expected_rt: f64,
    pub half_width: f64,
}

/// Wire form of a `GapFillTarget`. Field order alphabetical.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GapFillEntry {
    pub charge: u8,
    pub decoy_entry_id: u32,
    pub expected_rt: f64,
    pub half_width: f64,
    pub modified_sequence: String,
    pub target_entry_id: u32,
}

/// Wire form of the refined per-file RT calibration. Carries the LOESS
/// model parameters; Stage 6 workers reconstruct an `RTCalibration` via
/// `RTCalibration::from_model_params`.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RefinedRtCalibrationJson {
    pub abs_residuals: Vec<f64>,
    pub fitted_rts: Vec<f64>,
    pub library_rts: Vec<f64>,
    pub residual_sd: f64,
}

/// Current format version. Bump on incompatible schema changes.
///
/// v1: initial format.
/// v2: added `file_stems` so per-file rescore workers can compute the
///     reconciliation parameter hash that the downstream
///     `--join-at-pass=2` merge node expects (which is computed over
///     all files participating in the join, not the worker's single
///     parquet). Required (no `#[serde(default)]`): v1 envelopes are
///     rejected at deserialization rather than falling back to a
///     single-stem hash the merge node would reject anyway.
/// v3: added `first_pass_base_ids`, the join-wide set of base_ids that
///     survived first-pass compaction. A per-file rescore worker compacts
///     to exactly this set instead of recomputing a per-file subset that
///     drops cross-file entries, so the HPC chain matches the in-memory
///     pipeline. Required in v3 (no `#[serde(default)]`): v2 envelopes are
///     rejected at parse time.
pub const RECONCILIATION_FORMAT_VERSION: u32 = 3;

impl ReconciliationFile {
    /// Build the wire envelope for a single file. Filters
    /// `reconciliation_actions` (which is keyed by `(file_name, vec_idx)`
    /// across all files) down to entries belonging to `file_name`,
    /// resolves each `vec_idx` to its `entry_id` via `file_entries`, and
    /// splits non-Keep actions into the two homogeneous arrays.
    ///
    /// O(num_actions) per call; callers writing many files in a row
    /// should prefer [`Self::from_planner_output_pre_grouped`] together
    /// with a single up-front grouping to keep the total cost
    /// O(num_actions) rather than O(num_files * num_actions).
    #[allow(clippy::too_many_arguments)]
    pub fn from_planner_output(
        file_name: &str,
        file_entries: &[FdrEntry],
        reconciliation_actions: &HashMap<(String, usize), ReconcileAction>,
        gap_fill_targets: &[GapFillTarget],
        refined_calibration: Option<&RTCalibration>,
        search_hash: &str,
        library_hash: &str,
        file_stems: &[String],
        first_pass_base_ids: &HashSet<u32>,
    ) -> Self {
        let file_actions: Vec<(usize, &ReconcileAction)> = reconciliation_actions
            .iter()
            .filter_map(|((fname, vec_idx), action)| {
                if fname == file_name {
                    Some((*vec_idx, action))
                } else {
                    None
                }
            })
            .collect();
        Self::from_planner_output_pre_grouped(
            file_entries,
            &file_actions,
            gap_fill_targets,
            refined_calibration,
            search_hash,
            library_hash,
            file_stems,
            first_pass_base_ids,
        )
    }

    /// Build the wire envelope when the caller has already grouped
    /// reconciliation actions by file. `file_actions` is the slice of
    /// `(vec_idx, action)` entries that belong to this file (Keep actions
    /// included or not — they are filtered here either way). Use this
    /// when emitting envelopes for many files so the per-file cost stays
    /// O(actions_for_this_file) rather than O(total_actions).
    #[allow(clippy::too_many_arguments)]
    pub fn from_planner_output_pre_grouped(
        file_entries: &[FdrEntry],
        file_actions: &[(usize, &ReconcileAction)],
        gap_fill_targets: &[GapFillTarget],
        refined_calibration: Option<&RTCalibration>,
        search_hash: &str,
        library_hash: &str,
        file_stems: &[String],
        first_pass_base_ids: &HashSet<u32>,
    ) -> Self {
        // Sort and deduplicate the file_stems for deterministic output
        // — matches the sort+dedup used by reconciliation_parameter_hash
        // on the read side, so the worker hash matches the merge node hash.
        let mut stems: Vec<String> = file_stems.to_vec();
        stems.sort();
        stems.dedup();
        // Sort the join-wide passing base_ids ascending for deterministic,
        // byte-parity output (matches the Array.Sort on the C# side).
        let mut base_ids: Vec<u32> = first_pass_base_ids.iter().copied().collect();
        base_ids.sort_unstable();
        let mut use_cwt: Vec<UseCwtPeakEntry> = Vec::new();
        let mut forced: Vec<ForcedIntegrationEntry> = Vec::new();
        for (vec_idx, action) in file_actions {
            let entry_id = match file_entries.get(*vec_idx) {
                Some(e) => e.entry_id,
                None => continue,
            };
            match action {
                ReconcileAction::Keep => {}
                ReconcileAction::UseCwtPeak {
                    candidate_idx,
                    start_rt,
                    apex_rt,
                    end_rt,
                } => {
                    use_cwt.push(UseCwtPeakEntry {
                        apex_rt: *apex_rt,
                        candidate_idx: *candidate_idx as u32,
                        end_rt: *end_rt,
                        entry_id,
                        start_rt: *start_rt,
                    });
                }
                ReconcileAction::ForcedIntegration {
                    expected_rt,
                    half_width,
                } => {
                    forced.push(ForcedIntegrationEntry {
                        entry_id,
                        expected_rt: *expected_rt,
                        half_width: *half_width,
                    });
                }
            }
        }
        // Sort each array by entry_id for deterministic output.
        use_cwt.sort_by_key(|e| e.entry_id);
        forced.sort_by_key(|e| e.entry_id);

        let gap: Vec<GapFillEntry> = {
            let mut v: Vec<GapFillEntry> = gap_fill_targets
                .iter()
                .map(|g| GapFillEntry {
                    charge: g.charge,
                    decoy_entry_id: g.decoy_entry_id,
                    expected_rt: g.expected_rt,
                    half_width: g.half_width,
                    modified_sequence: g.modified_sequence.to_string(),
                    target_entry_id: g.target_entry_id,
                })
                .collect();
            v.sort_by_key(|e| e.target_entry_id);
            v
        };

        let cal = refined_calibration.map(|c| RefinedRtCalibrationJson {
            abs_residuals: c.abs_residuals().to_vec(),
            fitted_rts: c.fitted_values().to_vec(),
            library_rts: c.library_rts().to_vec(),
            residual_sd: c.residual_std(),
        });

        Self {
            file_stems: stems,
            first_pass_base_ids: base_ids,
            forced_integration_actions: forced,
            format_version: RECONCILIATION_FORMAT_VERSION,
            gap_fill_targets: gap,
            library_hash: library_hash.to_string(),
            refined_rt_calibration: cal,
            search_hash: search_hash.to_string(),
            use_cwt_peak_actions: use_cwt,
        }
    }
}

/// Compute the per-file reconciliation JSON path in the resolved output
/// directory: `<output_dir>/<stem>.reconciliation.json`. Mirrors the
/// existing pattern used for `<stem>.calibration.json`,
/// `<stem>.scores.parquet`, etc. The filename is unchanged; only the
/// directory is redirected. Pass the directory from
/// [`osprey_core::OspreyConfig::resolve_output_dir`].
pub fn reconciliation_path(input_path: &Path, output_dir: &Path) -> std::path::PathBuf {
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    output_dir.join(format!("{}.reconciliation.json", stem))
}

/// Write the boundary file as pretty JSON with 2-space indent + LF line
/// endings. Stages through a sibling `.tmp` file in the same directory
/// and renames into place; this avoids leaving a partially-written
/// destination on serializer failure, but the rename is not strictly
/// atomic when the destination already exists (the existing file is
/// removed first). A crash between the remove and the rename leaves
/// the `.tmp` next to the missing destination, which the next run
/// either overwrites or — if the writer fails identically — leaves
/// recoverable by hand.
pub fn write_reconciliation_file(path: &Path, file: &ReconciliationFile) -> Result<()> {
    let parent = path.parent().unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent).map_err(|e| {
        anyhow!(
            "Failed to ensure reconciliation output dir {}: {}",
            parent.display(),
            e
        )
    })?;
    let tmp_path = parent.join(format!(
        ".{}.tmp",
        path.file_name()
            .unwrap_or(std::ffi::OsStr::new("reconciliation.json"))
            .to_string_lossy()
    ));
    {
        let f = std::fs::File::create(&tmp_path)
            .map_err(|e| anyhow!("Failed to create {}: {}", tmp_path.display(), e))?;
        let mut writer = std::io::BufWriter::new(f);
        let formatter = RoundtripPrettyFormatter::new();
        let mut serializer = serde_json::Serializer::with_formatter(&mut writer, formatter);
        file.serialize(&mut serializer)
            .map_err(|e| anyhow!("Failed to serialize reconciliation JSON: {}", e))?;
        // The serializer does not write a trailing newline; emit one
        // explicitly so the file ends with `}\n` (POSIX text file
        // convention) and the C# side mirrors it.
        use std::io::Write;
        writer.write_all(b"\n")?;
    }
    if path.exists() {
        std::fs::remove_file(path).ok();
    }
    std::fs::rename(&tmp_path, path).map_err(|e| {
        anyhow!(
            "Failed to rename {} to {}: {}",
            tmp_path.display(),
            path.display(),
            e
        )
    })?;
    Ok(())
}

/// Read and parse the boundary file. Validates `format_version` against
/// `RECONCILIATION_FORMAT_VERSION` and returns an error on mismatch.
pub fn read_reconciliation_file(path: &Path) -> Result<ReconciliationFile> {
    let f = std::fs::File::open(path)
        .map_err(|e| anyhow!("Failed to open {}: {}", path.display(), e))?;
    let reader = std::io::BufReader::new(f);
    let parsed: ReconciliationFile = serde_json::from_reader(reader).map_err(|e| {
        anyhow!(
            "Failed to parse reconciliation JSON {}: {}",
            path.display(),
            e
        )
    })?;
    if parsed.format_version != RECONCILIATION_FORMAT_VERSION {
        return Err(anyhow!(
            "Reconciliation file {} has unsupported format_version {} (expected {})",
            path.display(),
            parsed.format_version,
            RECONCILIATION_FORMAT_VERSION
        ));
    }
    Ok(parsed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_file() -> ReconciliationFile {
        // NOTE: keep these hardcoded inputs byte-identical to the C#
        // sibling `MakeSampleReconciliationFile` so the cross-impl
        // byte-parity hook below (OSPREY_CROSS_IMPL_RECONCILIATION_OUT)
        // compares like against like.
        ReconciliationFile {
            file_stems: vec!["round_trip".to_string()],
            first_pass_base_ids: vec![3, 100, 101, 200, 201],
            forced_integration_actions: vec![
                ForcedIntegrationEntry {
                    entry_id: 200,
                    expected_rt: 41.125,
                    half_width: 0.075,
                },
                ForcedIntegrationEntry {
                    entry_id: 201,
                    expected_rt: 18.5,
                    half_width: 0.05,
                },
            ],
            format_version: RECONCILIATION_FORMAT_VERSION,
            gap_fill_targets: vec![GapFillEntry {
                charge: 2,
                decoy_entry_id: 0x80000003,
                expected_rt: 33.5,
                half_width: 0.08,
                modified_sequence: "PEPTIDE".to_string(),
                target_entry_id: 3,
            }],
            library_hash: "lib-hash-abc".to_string(),
            refined_rt_calibration: Some(RefinedRtCalibrationJson {
                abs_residuals: vec![0.01, 0.02, 0.015],
                fitted_rts: vec![10.5, 20.5, 30.5],
                library_rts: vec![10.0, 20.0, 30.0],
                residual_sd: 0.123,
            }),
            search_hash: "search-hash-xyz".to_string(),
            use_cwt_peak_actions: vec![
                UseCwtPeakEntry {
                    apex_rt: 23.45,
                    candidate_idx: 1,
                    end_rt: 23.80,
                    entry_id: 100,
                    start_rt: 23.10,
                },
                UseCwtPeakEntry {
                    apex_rt: 8.07,
                    candidate_idx: 0,
                    end_rt: 8.20,
                    entry_id: 101,
                    start_rt: 7.95,
                },
            ],
        }
    }

    /// Regression guard for `serde_json`'s default f64 parser, which is
    /// not always correctly rounded — for the shortest-roundtrip string
    /// `"3.1575921556296254"` it returns bits `0x...878` while
    /// `f64::from_str` returns the correct `0x...877`. Our boundary file
    /// round-trips f64 calibration arrays through this parser, and a
    /// 1-ULP drift propagates to per-entry `rt_deviation` divergence
    /// at Stage 6 (Session 5 bisection). The fix is the
    /// `float_roundtrip` feature in the workspace `Cargo.toml`. This
    /// test asserts every f64 it touches round-trips bit-exactly so a
    /// future feature-flag regression or serde_json upgrade is caught.
    #[test]
    fn serde_json_f64_roundtrip_is_bit_exact() {
        // Values that diverged in the cross-impl bisection of file_20's
        // refined RT calibration's `fitted_values` array.
        let candidates = [
            3.1575921556296254_f64,
            3.1585537473115846_f64,
            3.1741410573166426_f64,
            // Stress: subnormal-near, integer-valued, and a typical
            // q-value so the test fires for the broader Stage-N JSON
            // path beyond just calibration arrays.
            f64::from_bits(0x0010_0000_0000_0000), // smallest normal
            1.0_f64,
            0.0123_f64,
        ];
        for &v in &candidates {
            let s = osprey_core::diagnostics::format_f64_roundtrip(v);
            let parsed: f64 = serde_json::from_str(&s).unwrap_or_else(|e| {
                panic!("serde_json failed to parse {:?}: {}", s, e);
            });
            assert_eq!(
                v.to_bits(),
                parsed.to_bits(),
                "f64 round-trip lost bits: input={:?} ({:#x}) -> str={:?} -> parsed={:?} ({:#x}). \
                 Likely cause: serde_json `float_roundtrip` feature is not enabled (see workspace Cargo.toml).",
                v,
                v.to_bits(),
                s,
                parsed,
                parsed.to_bits()
            );
        }
    }

    #[test]
    fn reconciliation_file_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("round_trip.reconciliation.json");
        let file = sample_file();

        write_reconciliation_file(&path, &file).unwrap();

        // Cross-impl byte-parity hook: when the harness runs this test
        // with `OSPREY_CROSS_IMPL_RECONCILIATION_OUT=<path>` set, copy
        // our output to that path so a sibling test on the OspreySharp
        // side (using identical hardcoded inputs) can be byte-compared.
        if let Ok(out) = std::env::var("OSPREY_CROSS_IMPL_RECONCILIATION_OUT") {
            std::fs::copy(&path, &out).unwrap();
        }

        let parsed = read_reconciliation_file(&path).unwrap();

        assert_eq!(parsed.format_version, file.format_version);
        assert_eq!(parsed.search_hash, file.search_hash);
        assert_eq!(parsed.library_hash, file.library_hash);
        assert_eq!(
            parsed.use_cwt_peak_actions.len(),
            file.use_cwt_peak_actions.len()
        );
        assert_eq!(
            parsed.forced_integration_actions.len(),
            file.forced_integration_actions.len()
        );
        assert_eq!(parsed.gap_fill_targets.len(), file.gap_fill_targets.len());
        assert_eq!(parsed.first_pass_base_ids, file.first_pass_base_ids);
        // Spot-check bit-exact f64 round-trip on a non-trivial value.
        let cwt_orig = &file.use_cwt_peak_actions[0];
        let cwt_got = &parsed.use_cwt_peak_actions[0];
        assert_eq!(cwt_orig.apex_rt.to_bits(), cwt_got.apex_rt.to_bits());
        assert_eq!(cwt_orig.start_rt.to_bits(), cwt_got.start_rt.to_bits());
        assert_eq!(cwt_orig.end_rt.to_bits(), cwt_got.end_rt.to_bits());
    }

    #[test]
    fn reconciliation_file_format_version_mismatch_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad_version.reconciliation.json");
        // Hand-craft a file with an out-of-range format_version. Field
        // order matches the alphabetical struct declaration so the
        // reader doesn't trip over deny_unknown_fields. All required
        // fields (file_stems, first_pass_base_ids) are present so serde
        // parse succeeds and the runtime format_version check fires.
        let bad = r#"{
  "file_stems": ["a", "b"],
  "first_pass_base_ids": [1, 2],
  "forced_integration_actions": [],
  "format_version": 99,
  "gap_fill_targets": [],
  "library_hash": "x",
  "refined_rt_calibration": null,
  "search_hash": "y",
  "use_cwt_peak_actions": []
}
"#;
        std::fs::write(&path, bad).unwrap();
        let err = read_reconciliation_file(&path).unwrap_err().to_string();
        assert!(err.contains("unsupported format_version"), "got: {}", err);
    }

    /// v1 envelopes (no `file_stems` field) are rejected at parse time:
    /// the missing-field serde error fires before the runtime
    /// `format_version` check. This locks in the cross-impl alignment
    /// with OspreySharp which also hard-rejects v1.
    #[test]
    fn reconciliation_file_v1_envelope_rejected_at_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("v1.reconciliation.json");
        let v1 = r#"{
  "forced_integration_actions": [],
  "format_version": 1,
  "gap_fill_targets": [],
  "library_hash": "x",
  "refined_rt_calibration": null,
  "search_hash": "y",
  "use_cwt_peak_actions": []
}
"#;
        std::fs::write(&path, v1).unwrap();
        let err = read_reconciliation_file(&path).unwrap_err().to_string();
        assert!(
            err.contains("missing field `file_stems`"),
            "expected v1 rejection at serde parse, got: {}",
            err
        );
    }

    /// v2 envelopes (have `file_stems` but no `first_pass_base_ids`) are
    /// rejected at parse time: the missing-field serde error fires before
    /// the runtime `format_version` check. Mirrors the OspreySharp side,
    /// which hard-fails when `first_pass_base_ids` is absent so a per-file
    /// worker never silently recomputes a divergent per-file base_id subset.
    #[test]
    fn reconciliation_file_v2_envelope_rejected_at_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("v2.reconciliation.json");
        let v2 = r#"{
  "file_stems": ["a", "b"],
  "forced_integration_actions": [],
  "format_version": 2,
  "gap_fill_targets": [],
  "library_hash": "x",
  "refined_rt_calibration": null,
  "search_hash": "y",
  "use_cwt_peak_actions": []
}
"#;
        std::fs::write(&path, v2).unwrap();
        let err = read_reconciliation_file(&path).unwrap_err().to_string();
        assert!(
            err.contains("missing field `first_pass_base_ids`"),
            "expected v2 rejection at serde parse, got: {}",
            err
        );
    }
}
