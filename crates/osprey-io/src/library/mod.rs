//! Spectral library loading
//!
//! This module provides parsers for various spectral library formats:
//! - DIA-NN TSV format
//! - BiblioSpec blib format
//! - EncyclopeDIA elib format

mod blib;
mod diann;
mod elib;

pub use blib::BlibLoader;
pub use diann::DiannTsvLoader;
pub use elib::ElibLoader;

use osprey_core::{LibraryEntry, LibrarySource, OspreyError, Result};
use std::collections::HashMap;

/// Load a library from any supported source
pub fn load_library(source: &LibrarySource) -> Result<Vec<LibraryEntry>> {
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

    Ok(deduplicate_library(entries))
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
}
