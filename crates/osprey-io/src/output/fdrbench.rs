//! FDRBench input TSV writer for peptide / precursor level.
//!
//! Emits a tab-separated input file in the format expected by
//! [FDRBench](https://github.com/Noble-Lab/FDRBench) for evaluating FDR control
//! quality via entrapment counting. The output contains every target Osprey
//! scored (regardless of q-value) so FDRBench can compute true-FDR across the
//! full ranking. Decoys are excluded per FDRBench convention; entrapment
//! sequences (marked by `_p_target` in protein accessions) are targets in
//! Osprey's world and pass through unchanged.
//!
//! Output columns (per-precursor, default):
//!
//! ```text
//! peptide<TAB>mod_peptide<TAB>charge<TAB>q_value<TAB>score<TAB>protein
//! ```
//!
//! With `per_run = true`, an additional `run` column is appended.
//!
//! `score` is the raw SVM discriminant from [`FdrEntry::score`], i.e. the
//! upstream score that drives the q-value. Invoke FDRBench with
//! `-score 'score:1'` (higher is better).
//!
//! The `protein` column is capped at [`MAX_PROTEIN_FIELD_BYTES`] (under the
//! 4096-char per-column limit in FDRBench's Univocity CSV parser). Lists
//! that exceed the cap, e.g. peptides shared across large paralog families
//! such as zinc fingers, are truncated with a trailing `;...+N_more` marker
//! indicating how many IDs were dropped. The canonical blib and CSV outputs
//! still carry the full protein-ID list.

use osprey_core::config::FdrLevel;
use osprey_core::types::{FdrEntry, LibraryEntry};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Mask for extracting the target-side base_id from `FdrEntry.entry_id`.
///
/// The high bit is set for decoys; clearing it yields the library entry id
/// shared by a target and its paired decoy.
const BASE_ID_MASK: u32 = 0x7FFF_FFFF;

/// Maximum width of the `protein` column in the FDRBench TSV, in bytes.
///
/// FDRBench's bundled Univocity CSV parser caps each column at 4096
/// characters and aborts the run with an `ArrayIndexOutOfBoundsException`
/// when that limit is exceeded. Shared peptides (e.g. ZNF / olfactory
/// receptor families) can legitimately map to hundreds of protein IDs
/// joined with `;`, which trips the cap. We keep a margin under 4096 so
/// the surrounding tabs and a truncation marker fit comfortably.
const MAX_PROTEIN_FIELD_BYTES: usize = 4000;

/// Format a protein-ID list for the FDRBench TSV, truncating when the joined
/// string would exceed [`MAX_PROTEIN_FIELD_BYTES`].
///
/// On truncation, returns `<kept_ids_joined_with_';'>;...+N_more` where `N`
/// is the number of IDs dropped. Returns `(field, true)` when truncation
/// occurred. If even the first ID is longer than the budget (pathological),
/// it is hard-truncated to the budget and the marker is appended.
fn format_protein_field(ids: &[String]) -> (String, bool) {
    let joined = ids.join(";");
    if joined.len() <= MAX_PROTEIN_FIELD_BYTES {
        return (joined, false);
    }
    // Reserve worst-case marker width (uses the initial dropped count, which
    // shrinks monotonically as we keep more IDs).
    let marker_reserve = format!(";...+{}_more", ids.len()).len();
    let budget = MAX_PROTEIN_FIELD_BYTES.saturating_sub(marker_reserve);

    let mut kept = 0usize;
    let mut len = 0usize;
    for (i, id) in ids.iter().enumerate() {
        let sep_len = if i == 0 { 0 } else { 1 };
        if len + sep_len + id.len() > budget {
            break;
        }
        len += sep_len + id.len();
        kept = i + 1;
    }

    if kept == 0 {
        // First ID alone exceeds the budget. Hard-truncate it (byte-safe via
        // chars().take()) and append the marker.
        let truncated: String = ids[0].chars().take(budget).collect();
        return (format!("{}...+{}_more", truncated, ids.len()), true);
    }
    let dropped = ids.len() - kept;
    let mut s = ids[..kept].join(";");
    s.push_str(&format!(";...+{}_more", dropped));
    (s, true)
}

/// Write the FDRBench peptide / precursor-level input TSV.
///
/// Every non-decoy entry in `per_file_entries` contributes at least one row.
/// With `per_run = false`, rows are deduplicated to one per precursor (level
/// `Precursor` / `Both`) or one per peptide (level `Peptide`), keeping the
/// minimum experiment-level q-value (ties broken by highest discriminant).
/// With `per_run = true`, every (precursor, file) observation is emitted with
/// the per-run q-value and a `run` column populated from the file key.
///
/// The `protein` column is the `;`-joined `LibraryEntry.protein_ids` for the
/// library entry referenced by `FdrEntry.entry_id`. Entries whose library
/// lookup fails (which should not happen in practice) get an empty protein
/// field and are still emitted.
///
/// # Arguments
///
/// * `path` - Destination TSV path. Parent directories must already exist.
/// * `per_file_entries` - Pre-compaction FDR stubs grouped by source file.
/// * `library` - The full library slice; used to resolve `entry_id` to
///   protein accessions and to confirm the entry's modified sequence.
/// * `fdr_level` - Drives which q-value field is emitted and the dedup key.
///   `Both` collapses to precursor-level for output (it preserves charge).
///   `Protein` falls back to precursor here since this writer is for the
///   peptide-level FDRBench TSV.
/// * `per_run` - If true, emit one row per (precursor, file); otherwise dedup
///   across runs to the best observation per precursor.
///
/// # Errors
///
/// Returns an `io::Error` if the file cannot be created or written.
pub fn write_fdrbench_peptide_input(
    path: &Path,
    per_file_entries: &[(String, Vec<FdrEntry>)],
    library: &[LibraryEntry],
    fdr_level: FdrLevel,
    per_run: bool,
) -> std::io::Result<()> {
    let effective_level = match fdr_level {
        FdrLevel::Protein | FdrLevel::Both => FdrLevel::Precursor,
        other => other,
    };

    // Build entry_id -> &LibraryEntry lookup. Library entry IDs are not
    // guaranteed to match Vec positions across loaders, so use a HashMap.
    let library_by_id: HashMap<u32, &LibraryEntry> = library.iter().map(|e| (e.id, e)).collect();

    let mut file = std::fs::File::create(path)?;
    if per_run {
        writeln!(
            file,
            "peptide\tmod_peptide\tcharge\tq_value\tscore\tprotein\trun"
        )?;
    } else {
        writeln!(
            file,
            "peptide\tmod_peptide\tcharge\tq_value\tscore\tprotein"
        )?;
    }

    let mut n_rows: usize = 0;
    let mut n_missing_lib: usize = 0;
    let mut n_truncated_protein: usize = 0;

    if per_run {
        for (run_name, entries) in per_file_entries {
            for entry in entries {
                if entry.is_decoy {
                    continue;
                }
                let base_id = entry.entry_id & BASE_ID_MASK;
                let (peptide, protein) = match library_by_id.get(&base_id) {
                    Some(lib) => {
                        let (field, truncated) = format_protein_field(&lib.protein_ids);
                        if truncated {
                            n_truncated_protein += 1;
                        }
                        (lib.sequence.as_str(), field)
                    }
                    None => {
                        n_missing_lib += 1;
                        ("", String::new())
                    }
                };
                let q = entry.effective_run_qvalue(effective_level);
                writeln!(
                    file,
                    "{peptide}\t{mod_seq}\t{charge}\t{q:.10e}\t{score:.10e}\t{protein}\t{run_name}",
                    peptide = peptide,
                    mod_seq = entry.modified_sequence,
                    charge = entry.charge,
                    q = q,
                    score = entry.score,
                    protein = protein,
                    run_name = run_name,
                )?;
                n_rows += 1;
            }
        }
    } else {
        // Dedup across files: keep best (min q-value, ties by max score) per dedup key.
        // Dedup key uses the modified sequence directly (no charge for peptide-level).
        // We store the entry reference and its q-value so we can write at the end.
        let mut best: HashMap<DedupKey, BestRow> = HashMap::new();
        for (_run_name, entries) in per_file_entries {
            for entry in entries {
                if entry.is_decoy {
                    continue;
                }
                let q = entry.effective_experiment_qvalue(effective_level);
                let key = DedupKey::new(entry, effective_level);
                let candidate = BestRow {
                    q_value: q,
                    score: entry.score,
                    entry,
                };
                best.entry(key)
                    .and_modify(|cur| {
                        if candidate.q_value < cur.q_value
                            || (candidate.q_value == cur.q_value && candidate.score > cur.score)
                        {
                            *cur = candidate.clone();
                        }
                    })
                    .or_insert(candidate);
            }
        }

        for row in best.into_values() {
            let base_id = row.entry.entry_id & BASE_ID_MASK;
            let (peptide, protein) = match library_by_id.get(&base_id) {
                Some(lib) => {
                    let (field, truncated) = format_protein_field(&lib.protein_ids);
                    if truncated {
                        n_truncated_protein += 1;
                    }
                    (lib.sequence.as_str(), field)
                }
                None => {
                    n_missing_lib += 1;
                    ("", String::new())
                }
            };
            writeln!(
                file,
                "{peptide}\t{mod_seq}\t{charge}\t{q:.10e}\t{score:.10e}\t{protein}",
                peptide = peptide,
                mod_seq = row.entry.modified_sequence,
                charge = row.entry.charge,
                q = row.q_value,
                score = row.score,
                protein = protein,
            )?;
            n_rows += 1;
        }
    }

    log::info!(
        "Wrote FDRBench {level:?} input ({mode}) to {path}: {n_rows} rows",
        level = effective_level,
        mode = if per_run { "per-run" } else { "per-precursor" },
        path = path.display(),
        n_rows = n_rows,
    );
    if n_missing_lib > 0 {
        log::warn!(
            "{} FDRBench rows had no library entry for entry_id; protein column left blank",
            n_missing_lib
        );
    }
    if n_truncated_protein > 0 {
        // FDRBench's CSV parser caps each column at 4096 chars; we truncate
        // longer protein-ID lists (e.g. ZNF families) to keep the run alive.
        // The `;...+N_more` marker preserves the dropped-ID count.
        log::warn!(
            "{} FDRBench rows had protein-ID lists exceeding {} bytes; truncated with ';...+N_more' marker",
            n_truncated_protein,
            MAX_PROTEIN_FIELD_BYTES,
        );
    }
    Ok(())
}

#[derive(Clone)]
struct BestRow<'a> {
    q_value: f64,
    score: f64,
    entry: &'a FdrEntry,
}

#[derive(PartialEq, Eq, Hash)]
enum DedupKey {
    /// Precursor level: (modified_sequence, charge)
    Precursor(String, u8),
    /// Peptide level: modified_sequence only
    Peptide(String),
}

impl DedupKey {
    fn new(entry: &FdrEntry, level: FdrLevel) -> Self {
        match level {
            FdrLevel::Peptide => DedupKey::Peptide(entry.modified_sequence.to_string()),
            // Precursor, Both, Protein all dedup by (modseq, charge) here.
            _ => DedupKey::Precursor(entry.modified_sequence.to_string(), entry.charge),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use tempfile::tempdir;

    fn make_entry(
        entry_id: u32,
        is_decoy: bool,
        modseq: &str,
        charge: u8,
        run_q_precursor: f64,
        exp_q_precursor: f64,
        score: f64,
    ) -> FdrEntry {
        FdrEntry {
            entry_id,
            parquet_index: 0,
            is_decoy,
            charge,
            scan_number: 0,
            apex_rt: 0.0,
            start_rt: 0.0,
            end_rt: 0.0,
            coelution_sum: 0.0,
            score,
            run_precursor_qvalue: run_q_precursor,
            run_peptide_qvalue: run_q_precursor,
            experiment_precursor_qvalue: exp_q_precursor,
            experiment_peptide_qvalue: exp_q_precursor,
            experiment_protein_qvalue: 1.0,
            pep: 0.0,
            experiment_aggregate_score: 0.0,
            modified_sequence: Arc::from(modseq),
        }
    }

    fn make_lib(
        id: u32,
        sequence: &str,
        modseq: &str,
        charge: u8,
        proteins: Vec<&str>,
    ) -> LibraryEntry {
        LibraryEntry {
            id,
            sequence: sequence.to_string(),
            modified_sequence: modseq.to_string(),
            modifications: Vec::new(),
            charge,
            precursor_mz: 0.0,
            retention_time: 0.0,
            rt_calibrated: false,
            fragments: Vec::new(),
            protein_ids: proteins.iter().map(|s| s.to_string()).collect(),
            gene_names: Vec::new(),
            is_decoy: false,
        }
    }

    fn read_lines(path: &Path) -> Vec<String> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .map(String::from)
            .collect()
    }

    #[test]
    fn precursor_dedup_keeps_best_qvalue() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let lib = vec![make_lib(1, "PEPTIDE", "PEPTIDE", 2, vec!["P1"])];
        // Same precursor scored in two runs; second is better.
        let per_file = vec![
            (
                "run_a".to_string(),
                vec![make_entry(1, false, "PEPTIDE", 2, 0.05, 0.05, 1.2)],
            ),
            (
                "run_b".to_string(),
                vec![make_entry(1, false, "PEPTIDE", 2, 0.01, 0.01, 3.4)],
            ),
        ];
        write_fdrbench_peptide_input(&path, &per_file, &lib, FdrLevel::Precursor, false).unwrap();
        let lines = read_lines(&path);
        assert_eq!(lines.len(), 2, "header + 1 deduped row");
        assert!(lines[1].contains("PEPTIDE"));
        assert!(lines[1].contains("\tP1\t") || lines[1].ends_with("\tP1"));
        assert!(
            lines[1].contains("3.4"),
            "expected best score (3.4), got: {}",
            lines[1]
        );
    }

    #[test]
    fn peptide_level_dedups_across_charges() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let lib = vec![
            make_lib(1, "PEPTIDE", "PEPTIDE", 2, vec!["P1"]),
            make_lib(2, "PEPTIDE", "PEPTIDE", 3, vec!["P1"]),
        ];
        let per_file = vec![(
            "run_a".to_string(),
            vec![
                make_entry(1, false, "PEPTIDE", 2, 0.05, 0.05, 1.0),
                make_entry(2, false, "PEPTIDE", 3, 0.02, 0.02, 2.0),
            ],
        )];
        write_fdrbench_peptide_input(&path, &per_file, &lib, FdrLevel::Peptide, false).unwrap();
        let lines = read_lines(&path);
        assert_eq!(
            lines.len(),
            2,
            "header + 1 peptide row (dedup across charges)"
        );
    }

    #[test]
    fn decoys_are_excluded() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let lib = vec![make_lib(1, "PEPTIDE", "PEPTIDE", 2, vec!["P1"])];
        let per_file = vec![(
            "run_a".to_string(),
            vec![
                make_entry(1, false, "PEPTIDE", 2, 0.01, 0.01, 3.0),
                // High bit set = decoy id; is_decoy=true.
                make_entry(1 | 0x8000_0000, true, "DECOY_PEPTIDE", 2, 0.5, 0.5, -1.0),
            ],
        )];
        write_fdrbench_peptide_input(&path, &per_file, &lib, FdrLevel::Precursor, false).unwrap();
        let lines = read_lines(&path);
        assert_eq!(lines.len(), 2);
        assert!(!lines[1].contains("DECOY_"));
    }

    #[test]
    fn per_run_emits_row_per_observation() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let lib = vec![make_lib(1, "PEPTIDE", "PEPTIDE", 2, vec!["P1"])];
        let per_file = vec![
            (
                "run_a".to_string(),
                vec![make_entry(1, false, "PEPTIDE", 2, 0.05, 0.05, 1.2)],
            ),
            (
                "run_b".to_string(),
                vec![make_entry(1, false, "PEPTIDE", 2, 0.01, 0.01, 3.4)],
            ),
        ];
        write_fdrbench_peptide_input(&path, &per_file, &lib, FdrLevel::Precursor, true).unwrap();
        let lines = read_lines(&path);
        assert_eq!(lines.len(), 3, "header + one row per observation");
        assert!(lines[0].ends_with("run"));
        assert!(lines.iter().any(|l| l.ends_with("\trun_a")));
        assert!(lines.iter().any(|l| l.ends_with("\trun_b")));
    }

    #[test]
    fn format_protein_field_passes_short_lists_through() {
        let ids = vec!["P1".to_string(), "P2".to_string(), "P3".to_string()];
        let (field, truncated) = format_protein_field(&ids);
        assert_eq!(field, "P1;P2;P3");
        assert!(!truncated);
    }

    #[test]
    fn format_protein_field_truncates_oversize_list() {
        // 200 IDs of ~30 chars each = ~6000 bytes joined, well over 4000.
        let ids: Vec<String> = (0..200)
            .map(|i| format!("sp|A{:08}|PROT_{:04}", i, i))
            .collect();
        let (field, truncated) = format_protein_field(&ids);
        assert!(truncated, "expected truncation for oversize list");
        assert!(
            field.len() <= MAX_PROTEIN_FIELD_BYTES,
            "field len {} exceeds budget {}",
            field.len(),
            MAX_PROTEIN_FIELD_BYTES
        );
        assert!(
            field.contains(";...+") && field.ends_with("_more"),
            "expected ';...+N_more' marker; got tail: {}",
            &field[field.len().saturating_sub(40)..]
        );
        // First ID must survive intact.
        assert!(field.starts_with("sp|A00000000|PROT_0000"));
    }

    #[test]
    fn format_protein_field_handles_oversize_single_id() {
        // One pathological ID longer than the budget on its own.
        let ids = vec!["X".repeat(MAX_PROTEIN_FIELD_BYTES + 500)];
        let (field, truncated) = format_protein_field(&ids);
        assert!(truncated);
        assert!(field.len() <= MAX_PROTEIN_FIELD_BYTES);
        assert!(field.ends_with("...+1_more"));
    }

    #[test]
    fn protein_column_joins_with_semicolons() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("out.tsv");
        let lib = vec![make_lib(
            42,
            "PEPTIDE",
            "PEPTIDE",
            2,
            vec!["P1", "P2_p_target"],
        )];
        let per_file = vec![(
            "run_a".to_string(),
            vec![make_entry(42, false, "PEPTIDE", 2, 0.01, 0.01, 3.0)],
        )];
        write_fdrbench_peptide_input(&path, &per_file, &lib, FdrLevel::Precursor, false).unwrap();
        let lines = read_lines(&path);
        assert!(
            lines[1].contains("P1;P2_p_target"),
            "expected semicolon-joined accessions including entrapment marker; got: {}",
            lines[1]
        );
    }
}
