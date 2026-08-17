//! Normalization of Carafe's per-peptide pseudo-protein accessions.
//!
//! Mirrors the C# Osprey `CarafeProteinIdNormalizer` (ProteoWizard/pwiz) so the two
//! implementations report the same protein identities on a Carafe-built library.

use osprey_core::LibraryEntry;
use std::collections::{HashMap, HashSet};

/// The token Carafe appends to a source accession, followed by the peptide ordinal.
///
/// Used both as the cheap pre-filter and as the literal matched by [`strip_pep_suffixes`].
/// The scan is O(entries x accessions) over a multi-million-entry library, and a `contains`
/// on a 4-byte needle is far cheaper than the character walk that follows it.
const PEP_TOKEN: &str = "_pep";

/// Strip the per-peptide `_pepNNNNN` pseudo-protein suffix Carafe writes into a spectral
/// library's ProteinID column, so protein-level statistics see real proteins. Returns the
/// number of entries rewritten (0 when the library is already clean, which is the no-op every
/// non-Carafe library takes).
///
/// Carafe emits one synthetic accession PER PEPTIDE - `sp|O95139_pep00019|NDUB6_HUMAN` - so
/// every peptide becomes its own protein. Measured on the 3-file Stellar entrapment library:
/// 359,656 distinct target accessions that collapse to 23,874 real ones.
///
/// The damage is not confined to the protein report. `protein-compact` second-pass q-value
/// estimation builds its stratum from proteins with >= 2 DETECTED peptides, and a protein that
/// owns exactly one peptide can never qualify - so the stratum collapses, nothing recompetes,
/// and every reported q-value is the carried-over first-pass value with no sign of failure in
/// the output.
///
/// **Call before `deduplicate_library`.** Dedup groups on (modified sequence, charge)
/// and merges the group's accessions through a sort + dedup; stripping afterwards would leave
/// two entries that differed only by `_pepNNNNN` as duplicate identical accessions in that
/// merge. Stripping first lets the existing merge collapse them, which is also what
/// pre-stripping the source TSV produces - and the two must agree.
///
/// Entrapment markers are preserved: only the `_pep` + digits token is removed, so
/// `sp|Q9ULK4_p_target_pep00052|MED23_HUMAN_p_target` becomes
/// `sp|Q9ULK4_p_target|MED23_HUMAN_p_target` - exactly the form the pairing manifest uses - and
/// a `decoy_` prefix is likewise untouched.
pub fn normalize_carafe_protein_ids(entries: &mut [LibraryEntry]) -> usize {
    // One strip per DISTINCT accession rather than one per (entry, accession) pair. The map is
    // NOT small on the libraries this exists for - a per-peptide accession is distinct by
    // construction, so it holds roughly one entry per peptide (759,805 on the 3-file Stellar
    // entrapment library). It still pays for itself: an accession shared across charge states
    // and duplicate rows is stripped once, and the alternative is re-scanning every string.
    let mut cleaned: HashMap<String, String> = HashMap::new();
    let mut real_accessions: HashSet<String> = HashSet::new();
    let mut example: Option<(String, String)> = None;
    for entry in entries.iter() {
        for id in &entry.protein_ids {
            if !id.contains(PEP_TOKEN) || cleaned.contains_key(id.as_str()) {
                continue;
            }
            // None when the accession carries a bare "_pep" with no digits after it, which is
            // a real accession and must be left alone.
            let Some(clean) = strip_pep_suffixes(id) else {
                continue;
            };
            real_accessions.insert(clean.clone());
            if example.is_none() {
                example = Some((id.clone(), clean.clone()));
            }
            cleaned.insert(id.clone(), clean);
        }
    }
    if cleaned.is_empty() {
        return 0;
    }

    // Rewrite into a sorted, distinct list. Distinct matters: once the suffix is gone, several
    // of an entry's accessions can collapse onto one. Only TOUCHED entries are rewritten, so a
    // clean entry keeps the accession order its loader produced.
    let mut n_entries = 0usize;
    for entry in entries.iter_mut() {
        if !entry
            .protein_ids
            .iter()
            .any(|id| cleaned.contains_key(id.as_str()))
        {
            continue;
        }
        let mut mapped: Vec<String> = entry
            .protein_ids
            .iter()
            .map(|id| cleaned.get(id.as_str()).unwrap_or(id).clone())
            .collect();
        mapped.sort();
        mapped.dedup();
        entry.protein_ids = mapped;
        n_entries += 1;
    }

    let (example_id, example_clean) = example.unwrap_or_default();
    log::warn!(
        "Spectral library protein accessions carry Carafe's per-peptide '_pepNNNNN' suffix: \
         {} distinct accessions on {} entries collapse to {} real proteins (e.g. '{}' -> '{}'). \
         Stripping it - left in place every peptide is its own protein, which breaks protein \
         parsimony and picked-protein FDR, and leaves protein-compact second-pass q-values with \
         an empty stratum because no protein can reach 2 detected peptides.",
        cleaned.len(),
        n_entries,
        real_accessions.len(),
        example_id,
        example_clean
    );
    n_entries
}

/// Remove every `_pep` + digits token from one accession, or `None` when it carries none.
///
/// Hand-rolled rather than a regex: the pattern is a fixed ASCII token followed by ASCII
/// digits, and osprey-io has no regex dependency to add for it. Both the token and the digits
/// are ASCII, so every slice boundary below is a char boundary.
fn strip_pep_suffixes(id: &str) -> Option<String> {
    let bytes = id.as_bytes();
    let token = PEP_TOKEN.as_bytes();
    let mut out = String::with_capacity(id.len());
    let mut copied = 0usize;
    let mut i = 0usize;
    let mut stripped = false;
    while i + token.len() <= bytes.len() {
        if &bytes[i..i + token.len()] == token {
            let mut j = i + token.len();
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            // Digits are required, so a real accession containing a bare "_pep" is left alone.
            if j > i + token.len() {
                out.push_str(&id[copied..i]);
                copied = j;
                i = j;
                stripped = true;
                continue;
            }
        }
        i += 1;
    }
    if !stripped {
        return None;
    }
    out.push_str(&id[copied..]);
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_entry(id: u32, protein_ids: &[&str]) -> LibraryEntry {
        make_precursor(id, "PEPTIDEK", 2, protein_ids)
    }

    /// Distinct (sequence, charge) so a test can model several peptides of ONE protein -
    /// entries that share a precursor would be collapsed by `deduplicate_library` downstream.
    fn make_precursor(id: u32, sequence: &str, charge: u8, protein_ids: &[&str]) -> LibraryEntry {
        let mut entry = LibraryEntry::new(
            id,
            sequence.to_string(),
            sequence.to_string(),
            charge,
            500.0,
            10.0,
        );
        entry.protein_ids = protein_ids.iter().map(|s| s.to_string()).collect();
        entry
    }

    #[test]
    fn test_strips_pep_suffix_and_preserves_markers() {
        assert_eq!(
            strip_pep_suffixes("sp|O95139_pep00019|NDUB6_HUMAN").as_deref(),
            Some("sp|O95139|NDUB6_HUMAN")
        );
        // Entrapment and decoy markers survive: only the _pep + digits token goes.
        assert_eq!(
            strip_pep_suffixes("sp|Q9ULK4_p_target_pep00052|MED23_HUMAN_p_target").as_deref(),
            Some("sp|Q9ULK4_p_target|MED23_HUMAN_p_target")
        );
        assert_eq!(
            strip_pep_suffixes("decoy_sp|O95139_pep00019|NDUB6_HUMAN").as_deref(),
            Some("decoy_sp|O95139|NDUB6_HUMAN")
        );
        // Digits are required - a bare "_pep" is a real accession.
        assert_eq!(strip_pep_suffixes("sp|Q9NP61_pep|ARFG3_HUMAN"), None);
        assert_eq!(strip_pep_suffixes("sp|Q9NP61|ARFG3_HUMAN"), None);
    }

    #[test]
    fn test_normalize_collapses_pseudo_proteins() {
        let mut entries = vec![
            make_entry(0, &["sp|O95139_pep00019|NDUB6_HUMAN"]),
            make_entry(1, &["sp|O95139_pep00020|NDUB6_HUMAN"]),
        ];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 2);
        // Two pseudo-proteins become the SAME real one - the collapse the stratum needs.
        assert_eq!(entries[0].protein_ids, vec!["sp|O95139|NDUB6_HUMAN"]);
        assert_eq!(entries[1].protein_ids, entries[0].protein_ids);
    }

    #[test]
    fn test_normalize_dedups_within_one_entry() {
        // Once the suffix is gone, several of an entry's own accessions collapse onto one.
        let mut entries = vec![make_entry(
            0,
            &[
                "sp|O95139_pep00019|NDUB6_HUMAN",
                "sp|O95139_pep00020|NDUB6_HUMAN",
                "sp|Q9NP61|ARFG3_HUMAN",
            ],
        )];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 1);
        assert_eq!(
            entries[0].protein_ids,
            vec!["sp|O95139|NDUB6_HUMAN", "sp|Q9NP61|ARFG3_HUMAN"]
        );
    }

    #[test]
    fn test_normalize_leaves_clean_libraries_untouched() {
        // The no-op every non-Carafe library takes: the map is empty, so the function
        // returns before the rewrite loop. That early return is why this test canNOT also
        // prove order preservation - see
        // test_normalize_preserves_clean_entry_order_in_a_dirty_library for that arm.
        let mut entries = vec![make_entry(
            0,
            &["sp|Q9NP61|ARFG3_HUMAN", "sp|Q04726|TLE3_HUMAN"],
        )];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 0);
        assert_eq!(
            entries[0].protein_ids,
            vec!["sp|Q9NP61|ARFG3_HUMAN", "sp|Q04726|TLE3_HUMAN"]
        );
    }

    #[test]
    fn test_normalize_counts_only_touched_entries() {
        let mut entries = vec![
            make_entry(0, &["sp|O95139_pep00019|NDUB6_HUMAN"]),
            make_entry(1, &["sp|Q9NP61|ARFG3_HUMAN"]),
        ];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 1);
        assert_eq!(entries[1].protein_ids, vec!["sp|Q9NP61|ARFG3_HUMAN"]);
    }

    #[test]
    fn test_normalize_preserves_clean_entry_order_in_a_dirty_library() {
        // The DISCRIMINATING arm for the touched-only guard. A wholly clean library returns
        // early on the empty map, so it cannot detect an unconditional sort - deleting the
        // guard and sorting every entry leaves such a test green. Here the library IS dirty,
        // so the rewrite loop runs, and the clean entry's deliberately un-sorted accessions
        // must come through in their original order.
        let mut entries = vec![
            make_entry(0, &["sp|O95139_pep00019|NDUB6_HUMAN"]),
            make_entry(1, &["sp|Q9NP61|ARFG3_HUMAN", "sp|Q04726|TLE3_HUMAN"]),
        ];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 1);
        assert_eq!(entries[0].protein_ids, vec!["sp|O95139|NDUB6_HUMAN"]);
        assert_eq!(
            entries[1].protein_ids,
            vec!["sp|Q9NP61|ARFG3_HUMAN", "sp|Q04726|TLE3_HUMAN"],
            "an untouched entry must keep its loader's accession order"
        );
    }

    #[test]
    fn test_normalize_keeps_a_real_accession_containing_bare_pep() {
        // Pinned at the normalize level, not only inside strip_pep_suffixes: an accession
        // whose real name contains "_pep" with no digits after it must survive untouched,
        // and must not by itself make the library look dirty.
        let mut entries = vec![make_entry(0, &["sp|P00761_peptidase|TRYP_PIG"])];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 0);
        assert_eq!(entries[0].protein_ids, vec!["sp|P00761_peptidase|TRYP_PIG"]);
    }

    #[test]
    fn test_strips_every_token_not_only_the_first() {
        // Covers the scanner's re-synchronization after a match. With `copied` advanced to the
        // wrong offset, a second token's digits leak back into the result.
        assert_eq!(
            strip_pep_suffixes("sp|A_pep00019|B_pep00020|C").as_deref(),
            Some("sp|A|B|C")
        );
        // Overlapping token: the leading bare "_pep" has no digits and must survive.
        assert_eq!(
            strip_pep_suffixes("sp|A_pep_pep00019|B").as_deref(),
            Some("sp|A_pep|B")
        );
        // A token at the very end of the string - the loop's last position.
        assert_eq!(
            strip_pep_suffixes("sp|A|B_pep00019").as_deref(),
            Some("sp|A|B")
        );
    }

    #[test]
    fn test_normalize_collapses_distinct_peptides_of_one_protein() {
        // The case the whole change exists for, and the one the other tests cannot model:
        // N DIFFERENT peptides of the same protein. Same-precursor entries would be merged by
        // deduplicate_library, so only distinct (sequence, charge) entries show that the
        // protein reaches the >= 2 detected-peptide stratum after stripping.
        let mut entries = vec![
            make_precursor(0, "PEPTIDEK", 2, &["sp|O95139_pep00019|NDUB6_HUMAN"]),
            make_precursor(1, "SAMPLERK", 3, &["sp|O95139_pep00020|NDUB6_HUMAN"]),
        ];
        assert_eq!(normalize_carafe_protein_ids(&mut entries), 2);
        let distinct: HashSet<&str> = entries
            .iter()
            .flat_map(|e| e.protein_ids.iter().map(String::as_str))
            .collect();
        assert_eq!(distinct.len(), 1, "two peptides must name ONE protein");
        assert!(distinct.contains("sp|O95139|NDUB6_HUMAN"));
    }
}
