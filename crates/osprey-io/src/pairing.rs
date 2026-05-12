//! Reader for FDRBench-style decoy pairing manifests.
//!
//! [FDRBench](https://github.com/Noble-Lab/FDRBench) emits a 5-column TSV
//! ("entrapment peptides" file) when it generates an entrapment library.
//! Each row describes one peptide and its membership in a `peptide_pair_index`
//! group that contains exactly four entries: a regular target, an entrapment
//! target (`p_target`), the regular target's randomized decoy, and the
//! entrapment target's randomized decoy.
//!
//! ```text
//! sequence                      decoy   proteins                                    peptide_type   peptide_pair_index
//! PSLDQLAAHPWMLGADGGVPESCDLR    No      sp|Q86V86|PIM3_HUMAN                         target         0
//! PALLAVGGADSLLEDGHQPCSWDMPR    No      sp|Q86V86_p_target|PIM3_HUMAN_p_target       p_target       0
//! CLGDLGLWAQPSHEDSVAPMADLGPR    Yes     rev_sp|Q86V86|PIM3_HUMAN                     decoy          0
//! LWLGSGAGDQPEAVLHPDLCMAPDSR    Yes     rev_sp|Q86V86_p_target|PIM3_HUMAN_p_target   p_decoy        0
//! ```
//!
//! From this we recover two target-decoy pairs per `peptide_pair_index`:
//! `(target, decoy)` and `(p_target, p_decoy)`. The pairing is used to rewrite
//! library decoy ids so each decoy shares a `base_id` with its target — what
//! the SVM/LDA target-decoy competition relies on.

use osprey_core::{LibraryEntry, PairingStats, DECOY_ID_BIT};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// One of the four peptide kinds in a FDRBench manifest row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PeptideKind {
    /// Original target peptide from a real protein.
    Target,
    /// Decoy of a real-protein target (randomized sequence).
    Decoy,
    /// Entrapment target peptide (a synthetic target).
    PTarget,
    /// Decoy of an entrapment target.
    PDecoy,
}

impl PeptideKind {
    fn parse(s: &str) -> Option<Self> {
        match s {
            "target" => Some(Self::Target),
            "decoy" => Some(Self::Decoy),
            "p_target" => Some(Self::PTarget),
            "p_decoy" => Some(Self::PDecoy),
            _ => None,
        }
    }

    /// Which "side" of a pair this kind sits on (true = target-like, false = decoy-like).
    fn is_target_side(self) -> bool {
        matches!(self, Self::Target | Self::PTarget)
    }

    /// Pairing partition: target↔decoy and p_target↔p_decoy are independent
    /// pairings within the same `peptide_pair_index`.
    fn partition(self) -> u8 {
        match self {
            Self::Target | Self::Decoy => 0,
            Self::PTarget | Self::PDecoy => 1,
        }
    }
}

/// In-memory representation of a FDRBench pairing manifest.
pub struct DecoyPairingManifest {
    /// Lookup by peptide sequence: gives `(kind, pair_index)`.
    ///
    /// Each sequence appears in exactly one manifest row, so this is a 1:1 map.
    seq_to_info: HashMap<String, (PeptideKind, u32)>,
}

impl DecoyPairingManifest {
    /// Parse a FDRBench `*_pep.txt` pairing manifest from disk.
    ///
    /// Expected header (tab-separated):
    /// `sequence  decoy  proteins  peptide_type  peptide_pair_index`
    ///
    /// Returns an error if the header is missing or column names don't match.
    /// Rows with an unknown `peptide_type` or a non-numeric `peptide_pair_index`
    /// are skipped with a `log::warn!`; this keeps a partially-malformed
    /// manifest usable rather than failing the whole pipeline on a typo.
    pub fn from_tsv(path: &Path) -> std::io::Result<Self> {
        let f = File::open(path)?;
        let mut reader = BufReader::new(f);

        let mut header = String::new();
        reader.read_line(&mut header)?;
        let header = header.trim_end_matches(['\r', '\n']);
        let cols: Vec<&str> = header.split('\t').collect();
        let i_seq = cols.iter().position(|c| *c == "sequence");
        let i_type = cols.iter().position(|c| *c == "peptide_type");
        let i_pair = cols.iter().position(|c| *c == "peptide_pair_index");
        let (i_seq, i_type, i_pair) = match (i_seq, i_type, i_pair) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!(
                        "FDRBench manifest header missing required columns \
                         (need sequence, peptide_type, peptide_pair_index). Got: {}",
                        header
                    ),
                ));
            }
        };

        let mut map: HashMap<String, (PeptideKind, u32)> = HashMap::new();
        let mut n_skipped: usize = 0;
        for (line_no, line) in reader.lines().enumerate() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() <= i_pair.max(i_type).max(i_seq) {
                n_skipped += 1;
                continue;
            }
            let kind = match PeptideKind::parse(fields[i_type]) {
                Some(k) => k,
                None => {
                    if n_skipped < 5 {
                        log::warn!(
                            "FDRBench manifest line {}: unknown peptide_type '{}', skipping",
                            line_no + 2,
                            fields[i_type]
                        );
                    }
                    n_skipped += 1;
                    continue;
                }
            };
            let pair_index = match fields[i_pair].parse::<u32>() {
                Ok(v) => v,
                Err(_) => {
                    n_skipped += 1;
                    continue;
                }
            };
            map.insert(fields[i_seq].to_string(), (kind, pair_index));
        }
        if n_skipped > 0 {
            log::info!(
                "FDRBench manifest: {} rows skipped (malformed or unknown peptide_type)",
                n_skipped
            );
        }
        log::info!(
            "Loaded FDRBench pairing manifest from {}: {} sequences",
            path.display(),
            map.len()
        );
        Ok(Self { seq_to_info: map })
    }

    /// Number of sequences indexed.
    pub fn len(&self) -> usize {
        self.seq_to_info.len()
    }

    /// True if no sequences are indexed.
    pub fn is_empty(&self) -> bool {
        self.seq_to_info.is_empty()
    }

    /// Apply the manifest to a library: for each pair_index/partition/charge
    /// group, identify target-side and decoy-side entries and rewrite each
    /// decoy's `id` so its `base_id` matches a target's `id`. Returns
    /// `PairingStats` reflecting how many decoys were successfully paired.
    ///
    /// Library entries whose sequence is not in the manifest are left
    /// untouched (and counted as unpaired if they have `is_decoy = true`).
    /// Within a single (pair_index, partition, charge) bucket, multiple
    /// target/decoy entries (e.g. modification variants) are sorted by
    /// `(sequence, id)` and zipped 1:1 deterministically.
    pub fn apply_to_library(&self, library: &mut [LibraryEntry]) -> PairingStats {
        // Group library entries by (pair_index, partition, charge, is_target_side).
        // Each bucket holds Vec<library_index>.
        let mut buckets: HashMap<(u32, u8, u8, bool), Vec<usize>> = HashMap::new();
        let mut n_targets: usize = 0;
        let mut n_decoys: usize = 0;
        for (idx, entry) in library.iter().enumerate() {
            if entry.is_decoy {
                n_decoys += 1;
            } else {
                n_targets += 1;
            }
            if let Some(&(kind, pair_index)) = self.seq_to_info.get(&entry.sequence) {
                buckets
                    .entry((
                        pair_index,
                        kind.partition(),
                        entry.charge,
                        kind.is_target_side(),
                    ))
                    .or_default()
                    .push(idx);
            }
        }

        // Walk every target-side bucket; find its decoy-side counterpart and
        // zip them 1:1 in deterministic order.
        let mut pairings: Vec<(usize, usize)> = Vec::new();
        let mut claimed_targets = std::collections::HashSet::new();
        // Collect target-side keys first so we can iterate deterministically.
        let mut target_keys: Vec<&(u32, u8, u8, bool)> = buckets.keys().filter(|k| k.3).collect();
        target_keys.sort();
        // Clone the keys so we can mutably look up the decoy side without
        // holding a borrow of `buckets`.
        let target_keys: Vec<(u32, u8, u8, bool)> = target_keys.iter().map(|k| **k).collect();
        for tkey in target_keys {
            let dkey = (tkey.0, tkey.1, tkey.2, false);
            let target_idxs = buckets.get(&tkey).cloned().unwrap_or_default();
            let decoy_idxs = buckets.get(&dkey).cloned().unwrap_or_default();
            // Deterministic sort.
            let mut t_sorted = target_idxs;
            let mut d_sorted = decoy_idxs;
            t_sorted.sort_by(|&a, &b| {
                library[a]
                    .sequence
                    .cmp(&library[b].sequence)
                    .then(library[a].id.cmp(&library[b].id))
            });
            d_sorted.sort_by(|&a, &b| {
                library[a]
                    .sequence
                    .cmp(&library[b].sequence)
                    .then(library[a].id.cmp(&library[b].id))
            });
            for (t_idx, d_idx) in t_sorted.iter().zip(d_sorted.iter()) {
                pairings.push((*d_idx, *t_idx));
                claimed_targets.insert(*t_idx);
            }
        }

        // Apply rewrites.
        let n_paired = pairings.len();
        for (decoy_idx, target_idx) in pairings {
            let target_id = library[target_idx].id;
            library[decoy_idx].id = target_id | DECOY_ID_BIT;
        }

        PairingStats {
            n_targets,
            n_decoys,
            n_paired,
            n_unpaired_decoys: n_decoys - n_paired,
            n_unpaired_targets: n_targets.saturating_sub(claimed_targets.len()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_entry(id: u32, sequence: &str, charge: u8, is_decoy: bool) -> LibraryEntry {
        let mut e = LibraryEntry::new(id, sequence.into(), sequence.into(), charge, 500.0, 10.0);
        e.is_decoy = is_decoy;
        if is_decoy {
            e.id |= DECOY_ID_BIT;
        }
        e
    }

    fn write_manifest(rows: &[&str]) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(
            f,
            "sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index"
        )
        .unwrap();
        for row in rows {
            writeln!(f, "{}", row).unwrap();
        }
        f
    }

    #[test]
    fn parses_minimal_manifest() {
        let f = write_manifest(&[
            "PEPTIDEA\tNo\tprotA\ttarget\t0",
            "PEPTIDEB\tNo\tprotA_p_target\tp_target\t0",
            "AEPTPIDE\tYes\trev_protA\tdecoy\t0",
            "BPETPIDE\tYes\trev_protA_p_target\tp_decoy\t0",
        ]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        assert_eq!(m.len(), 4);
    }

    #[test]
    fn header_with_extra_columns_is_accepted() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(
            f,
            "sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index\textra"
        )
        .unwrap();
        writeln!(f, "PEPTIDE\tNo\tprotA\ttarget\t0\tignored").unwrap();
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        assert_eq!(m.len(), 1);
    }

    #[test]
    fn missing_required_column_errors() {
        let mut f = NamedTempFile::new().unwrap();
        writeln!(f, "sequence\tdecoy\tproteins").unwrap();
        let err = DecoyPairingManifest::from_tsv(f.path()).err().unwrap();
        assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    #[test]
    fn pairs_target_with_decoy_and_p_target_with_p_decoy() {
        let f = write_manifest(&[
            "PEPTIDEA\tNo\tprotA\ttarget\t0",
            "PEPTIDEB\tNo\tprotA_p_target\tp_target\t0",
            "AEPTPIDE\tYes\trev_protA\tdecoy\t0",
            "BPETPIDE\tYes\trev_protA_p_target\tp_decoy\t0",
        ]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        let mut lib = vec![
            make_entry(1, "PEPTIDEA", 2, false),
            make_entry(2, "PEPTIDEB", 2, false),
            make_entry(3, "AEPTPIDE", 2, true), // decoy of PEPTIDEA
            make_entry(4, "BPETPIDE", 2, true), // p_decoy of PEPTIDEB
        ];
        let stats = m.apply_to_library(&mut lib);
        assert_eq!(stats.n_paired, 2);
        assert_eq!(stats.n_unpaired_decoys, 0);
        // The decoy of PEPTIDEA should base_id-match target id 1.
        let aeptpide = lib.iter().find(|e| e.sequence == "AEPTPIDE").unwrap();
        let bpetpide = lib.iter().find(|e| e.sequence == "BPETPIDE").unwrap();
        assert_eq!(aeptpide.id & 0x7FFF_FFFF, 1);
        assert_eq!(bpetpide.id & 0x7FFF_FFFF, 2);
    }

    #[test]
    fn target_decoy_paired_per_charge_state() {
        let f = write_manifest(&[
            "PEPTIDEA\tNo\tprotA\ttarget\t0",
            "AEPTPIDE\tYes\trev_protA\tdecoy\t0",
        ]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        let mut lib = vec![
            make_entry(1, "PEPTIDEA", 2, false),
            make_entry(2, "PEPTIDEA", 3, false),
            make_entry(3, "AEPTPIDE", 2, true),
            make_entry(4, "AEPTPIDE", 3, true),
        ];
        let stats = m.apply_to_library(&mut lib);
        assert_eq!(stats.n_paired, 2);
        // Charge 2 decoy pairs to charge 2 target; charge 3 to charge 3.
        let d_c2 = lib
            .iter()
            .find(|e| e.sequence == "AEPTPIDE" && e.charge == 2)
            .unwrap();
        let d_c3 = lib
            .iter()
            .find(|e| e.sequence == "AEPTPIDE" && e.charge == 3)
            .unwrap();
        let t_c2 = lib
            .iter()
            .find(|e| e.sequence == "PEPTIDEA" && e.charge == 2)
            .unwrap();
        let t_c3 = lib
            .iter()
            .find(|e| e.sequence == "PEPTIDEA" && e.charge == 3)
            .unwrap();
        assert_eq!(d_c2.id & 0x7FFF_FFFF, t_c2.id);
        assert_eq!(d_c3.id & 0x7FFF_FFFF, t_c3.id);
    }

    #[test]
    fn unmatched_sequences_are_unpaired() {
        let f = write_manifest(&["PEPTIDEA\tNo\tprotA\ttarget\t0"]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        // Library has a decoy not in the manifest.
        let mut lib = vec![
            make_entry(1, "PEPTIDEA", 2, false),
            make_entry(2, "UNKNOWNK", 2, true), // not in manifest
        ];
        let stats = m.apply_to_library(&mut lib);
        assert_eq!(stats.n_paired, 0);
        assert_eq!(stats.n_unpaired_decoys, 1);
    }
}
