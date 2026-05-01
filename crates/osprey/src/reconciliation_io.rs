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

use std::collections::HashMap;
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
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReconciliationFile {
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
pub const RECONCILIATION_FORMAT_VERSION: u32 = 1;

impl ReconciliationFile {
    /// Build the wire envelope for a single file. Filters
    /// `reconciliation_actions` (which is keyed by `(file_name, vec_idx)`
    /// across all files) down to entries belonging to `file_name`,
    /// resolves each `vec_idx` to its `entry_id` via `file_entries`, and
    /// splits non-Keep actions into the two homogeneous arrays.
    pub fn from_planner_output(
        file_name: &str,
        file_entries: &[FdrEntry],
        reconciliation_actions: &HashMap<(String, usize), ReconcileAction>,
        gap_fill_targets: &[GapFillTarget],
        refined_calibration: Option<&RTCalibration>,
        search_hash: &str,
        library_hash: &str,
    ) -> Self {
        let mut use_cwt: Vec<UseCwtPeakEntry> = Vec::new();
        let mut forced: Vec<ForcedIntegrationEntry> = Vec::new();
        for ((fname, vec_idx), action) in reconciliation_actions {
            if fname != file_name {
                continue;
            }
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

/// Compute the per-file reconciliation JSON path: sibling to the input
/// mzML at `<dir>/<stem>.reconciliation.json`. Mirrors the existing
/// pattern used for `<stem>.calibration.json`, `<stem>.spectra.bin`,
/// etc.
pub fn reconciliation_path(input_path: &Path) -> std::path::PathBuf {
    let stem = input_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown");
    let parent = input_path.parent().unwrap_or(Path::new("."));
    parent.join(format!("{}.reconciliation.json", stem))
}

/// Write the boundary file as pretty JSON with 2-space indent + LF line
/// endings. Atomic write via temp file in the same directory.
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
    std::fs::rename(&tmp_path, path)
        .map_err(|e| anyhow!("Failed to atomically rename to {}: {}", path.display(), e))?;
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
        ReconciliationFile {
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
        // reader doesn't trip over deny_unknown_fields.
        let bad = r#"{
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
}
