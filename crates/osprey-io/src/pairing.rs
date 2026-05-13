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

use osprey_core::{LibraryEntry, PairingState, DECOY_ID_BIT};
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

/// Stats returned by `DecoyPairingManifest::apply_to_library`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ManifestApplyStats {
    /// Number of decoys paired with a target on this call (decoy.id rewritten
    /// so its `base_id` matches the target's `id`).
    pub n_paired: usize,
    /// Number of library entries whose `is_decoy` was flipped from `false`
    /// to `true` because the manifest classified them as `decoy` or
    /// `p_decoy`. The protein-prefix scan in
    /// `apply_library_decoy_marking` misses these when the library
    /// predictor strips the decoy prefix from protein accessions; the
    /// manifest's sequence-based classification is authoritative.
    pub n_newly_marked_decoy: usize,
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
    /// decoy's `id` so its `base_id` matches a target's `id`. Returns the
    /// number of decoys paired on this call.
    ///
    /// Library entries whose sequence is not in the manifest are left
    /// untouched. Within a single (pair_index, partition, charge) bucket,
    /// multiple target/decoy entries (e.g. modification variants) are sorted
    /// by `(sequence, id)` and zipped 1:1 deterministically.
    ///
    /// `state` carries already-claimed targets and already-paired decoys so
    /// this function can be chained with a composition-based fallback pass
    /// without claiming the same target twice. New pairings made here are
    /// added to `state`.
    ///
    /// **Authoritative classification by manifest.** Any library entry whose
    /// sequence appears in the manifest as `decoy` or `p_decoy` gets
    /// `is_decoy = true` and `DECOY_ID_BIT` set on its `id`, even if it
    /// doesn't pair. This catches the common case where a library predictor
    /// (e.g. Carafe) stripped the protein-accession decoy prefix during
    /// processing, leaving Osprey's prefix scan unable to recognise the
    /// decoys. The manifest's `peptide_type` column is taken as the source
    /// of truth.
    pub fn apply_to_library(
        &self,
        library: &mut [LibraryEntry],
        state: &mut PairingState,
    ) -> ManifestApplyStats {
        // Group library entries by (pair_index, partition, charge, is_target_side).
        // Each bucket holds Vec<library_index>. Skip entries already paired
        // (decoys) or claimed (targets) by an earlier pass.
        //
        // Also collect the list of decoy-side library indices identified by
        // the manifest, so we can stamp `is_decoy = true` + `DECOY_ID_BIT`
        // on them afterwards (catches predictors that stripped the
        // protein-accession prefix from decoy entries).
        let mut buckets: HashMap<(u32, u8, u8, bool), Vec<usize>> = HashMap::new();
        let mut decoy_side_indices: Vec<usize> = Vec::new();
        for (idx, entry) in library.iter().enumerate() {
            if entry.is_decoy && state.paired_decoys.contains(&idx) {
                continue;
            }
            if !entry.is_decoy && state.claimed_targets.contains(&idx) {
                continue;
            }
            if let Some(&(kind, pair_index)) = self.seq_to_info.get(&entry.sequence) {
                let is_target_side = kind.is_target_side();
                buckets
                    .entry((pair_index, kind.partition(), entry.charge, is_target_side))
                    .or_default()
                    .push(idx);
                if !is_target_side {
                    decoy_side_indices.push(idx);
                }
            }
        }

        // Mark manifest-identified decoys BEFORE pairing. Idempotent: entries
        // already flagged keep their flag and high bit. Pairing will overwrite
        // the low 31 bits to the target's base_id for paired decoys; unpaired
        // decoys keep their original id with just the high bit set.
        let mut n_newly_marked_decoy: usize = 0;
        for &idx in &decoy_side_indices {
            if !library[idx].is_decoy {
                library[idx].is_decoy = true;
                n_newly_marked_decoy += 1;
            }
            if library[idx].id & DECOY_ID_BIT == 0 {
                library[idx].id |= DECOY_ID_BIT;
            }
        }

        // Walk every target-side bucket; find its decoy-side counterpart and
        // zip them 1:1 in deterministic order.
        let mut pairings: Vec<(usize, usize)> = Vec::new();
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
                state.claimed_targets.insert(*t_idx);
                state.paired_decoys.insert(*d_idx);
            }
        }

        // Apply pairing rewrites. For each paired decoy, set its id so that
        // base_id matches the target.
        let n_paired = pairings.len();
        for (decoy_idx, target_idx) in pairings {
            let target_id = library[target_idx].id;
            library[decoy_idx].id = target_id | DECOY_ID_BIT;
        }
        ManifestApplyStats {
            n_paired,
            n_newly_marked_decoy,
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
        let mut state = PairingState::new();
        let stats = m.apply_to_library(&mut lib, &mut state);
        assert_eq!(stats.n_paired, 2);
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
        let mut state = PairingState::new();
        let stats = m.apply_to_library(&mut lib, &mut state);
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
        let mut state = PairingState::new();
        let stats = m.apply_to_library(&mut lib, &mut state);
        assert_eq!(stats.n_paired, 0);
        // The unpaired decoy is not in state.paired_decoys.
        assert!(!state.paired_decoys.contains(&1));
    }

    /// Real-world failure mode: a library predictor (Carafe) strips the
    /// `decoy_` / `rev_` prefix from protein accessions, so the library
    /// loads with every decoy looking like a target (`is_decoy = false`).
    /// The manifest's sequence-based classification is authoritative and
    /// must flip those entries to `is_decoy = true`.
    #[test]
    fn manifest_flips_unflagged_decoys_to_is_decoy() {
        let f = write_manifest(&[
            "PEPTIDEA\tNo\tprotA\ttarget\t0",
            "AEPTPIDE\tYes\trev_protA\tdecoy\t0",
        ]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        // Library has the decoy peptide but with is_decoy=false (predictor
        // stripped the rev_ prefix). Both entries look like targets to Osprey.
        let mut lib = vec![
            make_entry(1, "PEPTIDEA", 2, false),
            make_entry(2, "AEPTPIDE", 2, false), // <-- looked like target
        ];
        assert!(!lib[1].is_decoy);
        assert_eq!(lib[1].id & DECOY_ID_BIT, 0);

        let mut state = PairingState::new();
        let stats = m.apply_to_library(&mut lib, &mut state);

        assert_eq!(stats.n_paired, 1);
        assert_eq!(stats.n_newly_marked_decoy, 1);
        // The manifest flipped the AEPTPIDE entry.
        let decoy = lib.iter().find(|e| e.sequence == "AEPTPIDE").unwrap();
        assert!(decoy.is_decoy, "manifest should flip is_decoy=true");
        assert_ne!(decoy.id & DECOY_ID_BIT, 0, "decoy id should have high bit");
        // Paired: base_id matches the target.
        assert_eq!(decoy.id & 0x7FFF_FFFF, 1);
    }

    /// An unpaired manifest-flagged decoy (no target-side counterpart in the
    /// library) still gets `is_decoy = true` and the high bit on its `id`,
    /// even though pairing skips it.
    #[test]
    fn manifest_marks_unpaired_decoy_without_pairing_it() {
        let f = write_manifest(&[
            // Target-side counterpart missing from the library entirely.
            "PEPTIDEA\tNo\tprotA\ttarget\t0",
            "AEPTPIDE\tYes\trev_protA\tdecoy\t0",
        ]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        let mut lib = vec![
            // Only the decoy peptide is in the library; its target sibling
            // PEPTIDEA isn't.
            make_entry(7, "AEPTPIDE", 2, false),
        ];
        let mut state = PairingState::new();
        let stats = m.apply_to_library(&mut lib, &mut state);
        assert_eq!(stats.n_paired, 0); // no target -> no pair
        assert_eq!(stats.n_newly_marked_decoy, 1);
        assert!(lib[0].is_decoy);
        assert_ne!(lib[0].id & DECOY_ID_BIT, 0);
        // Unpaired: low 31 bits keep the original id.
        assert_eq!(lib[0].id & 0x7FFF_FFFF, 7);
    }

    #[test]
    fn manifest_then_composition_pairs_unmatched_decoys() {
        // Hybrid mode: manifest covers PEPTIDEA's pair; composition covers
        // PEPTIDEB↔EPBPTIDE (a permutation pair not in the manifest).
        use osprey_core::pair_library_decoys_by_composition;
        let f = write_manifest(&[
            "PEPTIDEA\tNo\tprotA\ttarget\t0",
            "AEPTPIDE\tYes\trev_protA\tdecoy\t0",
        ]);
        let m = DecoyPairingManifest::from_tsv(f.path()).unwrap();
        let mut lib = vec![
            // Pair 0: in the manifest.
            {
                let mut e = make_entry(1, "PEPTIDEA", 2, false);
                e.protein_ids = vec!["protA".into()];
                e
            },
            {
                let mut e = make_entry(2, "AEPTPIDE", 2, true);
                e.protein_ids = vec!["rev_protA".into()];
                e
            },
            // Pair 1: NOT in the manifest but composition-pairable
            // (PEPTIDEB and EPBPTIDE share AA composition; same protein pair).
            {
                let mut e = make_entry(3, "PEPTIDEB", 2, false);
                e.protein_ids = vec!["protB".into()];
                e
            },
            {
                let mut e = make_entry(4, "EPBPTIDE", 2, true);
                e.protein_ids = vec!["rev_protB".into()];
                e
            },
        ];
        let mut state = PairingState::new();
        let n_manifest = m.apply_to_library(&mut lib, &mut state);
        assert_eq!(n_manifest.n_paired, 1);
        // Run composition fallback. Should pair the remaining (PEPTIDEB, EPBPTIDE).
        let n_composition =
            pair_library_decoys_by_composition(&mut lib, &["rev_".to_string()], &mut state);
        assert_eq!(n_composition, 1);
        // Verify both pairs have matching base_ids.
        let dec_a = lib.iter().find(|e| e.sequence == "AEPTPIDE").unwrap();
        let dec_b = lib.iter().find(|e| e.sequence == "EPBPTIDE").unwrap();
        assert_eq!(dec_a.id & 0x7FFF_FFFF, 1);
        assert_eq!(dec_b.id & 0x7FFF_FFFF, 3);
    }
}
