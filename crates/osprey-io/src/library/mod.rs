//! Spectral library loading
//!
//! This module provides parsers for various spectral library formats:
//! - DIA-NN TSV format
//! - BiblioSpec blib format
//! - EncyclopeDIA elib format

mod blib;
mod cache;
mod diann;
mod elib;

pub use blib::BlibLoader;
pub use cache::{load_library_cache, save_library_cache};
pub use diann::DiannTsvLoader;
pub use elib::ElibLoader;

use osprey_core::{LibraryEntry, LibrarySource, OspreyError, Result};
use std::collections::HashMap;

/// Shortest peptide a library may contain.
///
/// Deliberately a constant rather than a setting: a tryptic peptide below this length is not
/// specific enough to identify, and the field has converged on excluding them (MaxQuant uses
/// 6; the Carafe libraries the regression runs against start at 7).
pub const MIN_PEPTIDE_LENGTH: usize = 6;

/// Reject a library peptide shorter than [`MIN_PEPTIDE_LENGTH`].
///
/// A hard error rather than a skip: a library carrying peptides this short is almost certainly
/// built wrong, and silently dropping them would change the searched set without saying so.
/// Downstream code is allowed to rely on the bound. `DecoyGenerator::is_candidate_acceptable`
/// is the case that matters today: its fragment-overlap ratio carries a structural 1/(n-1)
/// floor, which is 0.2 against a 0.4 threshold at the enforced minimum but 0.5 at length 3,
/// where every candidate decoy would be rejected and the peptide would vanish from the search.
///
/// Called by every format loader. An invariant enforced on only one of the paths into the
/// library is not an invariant, and the decoy generator cannot tell which loader produced its
/// input.
pub fn validate_peptide_length(sequence: &str) -> Result<()> {
    let len = sequence.chars().count();
    if len >= MIN_PEPTIDE_LENGTH {
        return Ok(());
    }
    // Name the offending peptide: on a multi-million-row library an error that does not is
    // unactionable.
    Err(OspreyError::InvalidLibraryEntry(format!(
        "Library peptide '{sequence}' has {len} residues; the minimum supported length is \
         {MIN_PEPTIDE_LENGTH}. Peptides this short are not specific enough to identify and are \
         excluded by convention from spectral libraries. Rebuild the library with a longer \
         minimum."
    )))
}

/// Load a library from any supported source.
///
/// Checks for a binary cache file (`.libcache`) alongside the source library.
/// If the cache exists and is newer than the source, loads from cache directly
/// (skipping TSV/blib/elib parsing and deduplication). Otherwise, loads from
/// source, deduplicates, and saves the cache for future runs.
pub fn load_library(
    source: &LibrarySource,
    cache_dir: &std::path::Path,
) -> Result<Vec<LibraryEntry>> {
    let source_path = source.path();
    // Library cache filename ("<library-name>.libcache"), placed in cache_dir
    // (--cache-dir / --work-dir) rather than beside a possibly read-only library.
    // Only the directory moves; the filename is unchanged.
    let cache_name = source_path
        .file_name()
        .map(|n| format!("{}.libcache", n.to_string_lossy()))
        .unwrap_or_else(|| "library.libcache".to_string());
    let cache_path = cache_dir.join(cache_name);

    // Try loading from cache if it exists and is newer than the source file
    if cache_path.exists() {
        let cache_ok = match (
            std::fs::metadata(source_path),
            std::fs::metadata(&cache_path),
        ) {
            (Ok(src_meta), Ok(cache_meta)) => match (src_meta.modified(), cache_meta.modified()) {
                (Ok(src_time), Ok(cache_time)) => cache_time >= src_time,
                _ => false,
            },
            _ => false,
        };

        if cache_ok {
            match cache::load_library_cache(&cache_path) {
                Ok(entries) => {
                    log::info!(
                        "Loaded {} library entries from cache '{}'",
                        entries.len(),
                        cache_path.display()
                    );
                    return Ok(entries);
                }
                Err(e) => {
                    log::debug!("Library cache not usable, falling back to source: {}", e);
                }
            }
        }
    }

    // Load from source format
    let entries = match source {
        LibrarySource::DiannTsv(path) => {
            let loader = DiannTsvLoader::new();
            loader.load(path)
        }
        LibrarySource::Blib(path) => {
            let loader = BlibLoader::new();
            loader.load(path)
        }
        LibrarySource::Elib(path) => {
            let loader = ElibLoader::new();
            loader.load(path)
        }
        LibrarySource::SkylineDocument(path) => {
            // TODO: Implement Skyline document extraction
            Err(OspreyError::LibraryLoadError(format!(
                "Skyline document format not yet implemented: {}",
                path.display()
            )))
        }
    }?;

    let deduped = deduplicate_library(entries);

    // Save cache for future runs
    match cache::save_library_cache(&deduped, &cache_path) {
        Ok(()) => {
            log::info!(
                "Saved library cache ({} entries) to '{}'",
                deduped.len(),
                cache_path.display()
            );
        }
        Err(e) => {
            log::warn!("Failed to save library cache: {}", e);
        }
    }

    Ok(deduped)
}

/// Deduplicate library entries by (modified_sequence, charge).
///
/// When multiple entries share the same peptide+charge, keeps the one with the
/// most fragments (ties broken by highest total fragment intensity). Averages
/// retention_time across all duplicates and merges protein_ids and gene_names.
/// Re-assigns sequential IDs.
fn deduplicate_library(entries: Vec<LibraryEntry>) -> Vec<LibraryEntry> {
    let original_count = entries.len();

    // Group entries by (modified_sequence, charge)
    let mut groups: HashMap<(String, u8), Vec<LibraryEntry>> = HashMap::new();
    for entry in entries {
        groups
            .entry((entry.modified_sequence.clone(), entry.charge))
            .or_default()
            .push(entry);
    }

    let mut deduped: Vec<LibraryEntry> = Vec::with_capacity(groups.len());

    for (_key, mut group) in groups {
        if group.len() == 1 {
            deduped.push(group.pop().unwrap());
        } else {
            // Average retention time across all duplicates
            let avg_rt = group.iter().map(|e| e.retention_time).sum::<f64>() / group.len() as f64;

            // Merge protein_ids and gene_names from all entries
            let mut all_proteins: Vec<String> = group
                .iter()
                .flat_map(|e| e.protein_ids.iter().cloned())
                .collect();
            all_proteins.sort();
            all_proteins.dedup();

            let mut all_genes: Vec<String> = group
                .iter()
                .flat_map(|e| e.gene_names.iter().cloned())
                .collect();
            all_genes.sort();
            all_genes.dedup();

            // Pick the best entry: most fragments, then highest total intensity
            group.sort_by(|a, b| {
                let frag_cmp = b.fragments.len().cmp(&a.fragments.len());
                if frag_cmp != std::cmp::Ordering::Equal {
                    return frag_cmp;
                }
                let sum_a: f64 = a
                    .fragments
                    .iter()
                    .map(|f| f.relative_intensity as f64)
                    .sum();
                let sum_b: f64 = b
                    .fragments
                    .iter()
                    .map(|f| f.relative_intensity as f64)
                    .sum();
                sum_b.total_cmp(&sum_a)
            });

            let mut best = group.into_iter().next().unwrap();
            best.retention_time = avg_rt;
            best.protein_ids = all_proteins;
            best.gene_names = all_genes;
            deduped.push(best);
        }
    }

    // Sort deterministically before assigning IDs: HashMap iteration order is random,
    // so without sorting, entries get different IDs across runs — causing non-deterministic
    // calibration sampling, tiebreakers, and FDR results.
    deduped.sort_by(|a, b| {
        a.modified_sequence
            .cmp(&b.modified_sequence)
            .then(a.charge.cmp(&b.charge))
    });

    // Re-assign sequential IDs
    for (i, entry) in deduped.iter_mut().enumerate() {
        entry.id = i as u32;
    }

    let removed = original_count - deduped.len();
    if removed > 0 {
        log::info!(
            "Library deduplication: {} → {} entries ({} duplicates removed)",
            original_count,
            deduped.len(),
            removed
        );
    }

    deduped
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{FragmentAnnotation, IonType, LibraryFragment};

    fn make_fragment(mz: f64, intensity: f32) -> LibraryFragment {
        LibraryFragment {
            mz,
            relative_intensity: intensity,
            annotation: FragmentAnnotation {
                ion_type: IonType::Y,
                ordinal: 1,
                charge: 1,
                neutral_loss: None,
            },
        }
    }

    #[test]
    fn test_deduplicate_no_duplicates() {
        let entries = vec![
            LibraryEntry::new(0, "AAA".into(), "AAA".into(), 2, 300.0, 10.0),
            LibraryEntry::new(1, "BBB".into(), "BBB".into(), 2, 400.0, 20.0),
        ];
        let result = deduplicate_library(entries);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_deduplicate_removes_duplicates() {
        let mut e1 = LibraryEntry::new(0, "AAA".into(), "AAA".into(), 2, 300.0, 10.0);
        e1.fragments = vec![make_fragment(200.0, 100.0), make_fragment(300.0, 80.0)];
        e1.protein_ids = vec!["P1".into()];

        let mut e2 = LibraryEntry::new(1, "AAA".into(), "AAA".into(), 2, 300.0, 12.0);
        e2.fragments = vec![
            make_fragment(200.0, 90.0),
            make_fragment(300.0, 70.0),
            make_fragment(400.0, 50.0),
        ];
        e2.protein_ids = vec!["P2".into()];

        let mut e3 = LibraryEntry::new(2, "BBB".into(), "BBB".into(), 3, 500.0, 15.0);
        e3.fragments = vec![make_fragment(250.0, 100.0)];

        let result = deduplicate_library(vec![e1, e2, e3]);
        assert_eq!(result.len(), 2);

        // Find the AAA entry - should have 3 fragments (e2 was best)
        let aaa = result.iter().find(|e| e.sequence == "AAA").unwrap();
        assert_eq!(aaa.fragments.len(), 3);
        // Average RT of 10.0 and 12.0
        assert!((aaa.retention_time - 11.0).abs() < 0.01);
        // Merged proteins
        assert_eq!(aaa.protein_ids.len(), 2);
    }

    #[test]
    fn test_deduplicate_different_charges_are_separate() {
        let e1 = LibraryEntry::new(0, "AAA".into(), "AAA".into(), 2, 300.0, 10.0);
        let e2 = LibraryEntry::new(1, "AAA".into(), "AAA".into(), 3, 200.0, 10.0);
        let result = deduplicate_library(vec![e1, e2]);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_deduplicate_sequential_ids() {
        let e1 = LibraryEntry::new(100, "AAA".into(), "AAA".into(), 2, 300.0, 10.0);
        let e2 = LibraryEntry::new(200, "AAA".into(), "AAA".into(), 2, 300.0, 12.0);
        let e3 = LibraryEntry::new(300, "BBB".into(), "BBB".into(), 2, 400.0, 15.0);
        let result = deduplicate_library(vec![e1, e2, e3]);
        assert_eq!(result.len(), 2);
        // IDs should be 0 and 1 (sequential)
        let mut ids: Vec<u32> = result.iter().map(|e| e.id).collect();
        ids.sort();
        assert_eq!(ids, vec![0, 1]);
    }

    #[test]
    fn test_min_peptide_length_is_six() {
        // Pinned deliberately: this is a judgement encoded as a constant, not a tuning knob.
        // Changing it changes which libraries osprey will accept, so a change here should be
        // a conscious decision that trips this test first.
        assert_eq!(MIN_PEPTIDE_LENGTH, 6);
    }

    #[test]
    fn test_validate_peptide_length_enforces_the_bound() {
        let err = validate_peptide_length("PEPTK").unwrap_err();
        let msg = err.to_string();
        // The message has to identify WHICH peptide, or the error is unactionable on a
        // multi-million-row library.
        assert!(
            msg.contains("PEPTK"),
            "message did not name the peptide: {msg}"
        );
        assert!(msg.contains('6'), "message did not state the bound: {msg}");

        // Boundary: exactly MIN_PEPTIDE_LENGTH must pass. A test that only checked the
        // rejection case would pass just as well with an off-by-one that rejected 6.
        assert!(validate_peptide_length("PEPTIK").is_ok());
    }
}
