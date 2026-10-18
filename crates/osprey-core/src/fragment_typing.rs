//! A library spectrum's peaks as a mass spectrometer could have produced them, typed against
//! the peptide's own ions, and the one ranking every top-N selection uses.
//!
//! These are the rules every library reader applies before a search, whichever file the
//! library came from, and they mirror the C# Osprey to the bit (`pwiz_tools/Osprey/Osprey.Core`:
//! `FragmentPeaks`, `FragmentTyping`, `FragmentLadder`, `PeptideFragmentMass`,
//! `FragmentMath::TopIndicesByIntensity`); the cross-implementation gate holds the two to a
//! 1e-9 tolerance on every output.
//!
//! 1. **One peak per m/z, scaled to the base peak** ([`merge_equal_mz`]). A predicted library
//!    can list two ions of one peptide at exactly one m/z, each with its own intensity (Carafe
//!    lists `b2` and `b4^2` of a peptide whose first four residues repeat a composition as two
//!    rows at the same written m/z). No instrument can separate them: the spectrum has one
//!    peak there, with the summed intensity, and a search that extracted and scored two
//!    chromatograms at one m/z would count one measurement twice. Two peaks at different m/z,
//!    however close, are what the library states and stay separate.
//! 2. **Typing from m/z** ([`type_fragments`]). Most intense first, each peak takes the
//!    preferred primary b or y ion within the search's fragment tolerance that no peak has
//!    claimed: the nearest, counting ions within [`TIE_TOLERANCE`] of each other as equally
//!    near, then the lower fragment charge, then y before b, then the shorter ion. The types a
//!    library states are a report only; a peak no ion reaches is `Unknown`.
//! 3. **Residue masses** ([`residue_mass`]) are the values Skyline calculates from the residue
//!    formulas and `pwiz.Common.Chemistry`'s isotope table (its traditional C, H, N, O and S,
//!    NIST's Se), to the last digit; the proton is the one Skyline's m/z use.
//! 4. **Ties on intensity** ([`top_fragment_indices`]) are broken by the fragments themselves -
//!    lower m/z first, then ion type, ordinal and charge - never by their position in the
//!    library, so a selection does not depend on the order a library lists its peaks in.

use crate::{
    FragmentAnnotation, FragmentToleranceConfig, IonType, LibraryEntry, LibraryFragment,
    Modification, ToleranceUnit,
};
use std::cmp::Ordering;
use std::collections::HashMap;

/// The proton mass Skyline's peptide and fragment m/z use (`BioMassCalc.MassProton`).
pub const PROTON_MASS: f64 = 1.007276;
/// Monoisotopic mass of water, from the shared chemistry.
pub const H2O_MASS: f64 = 18.0105647;
/// Monoisotopic mass of ammonia, from the shared chemistry.
pub const NH3_MASS: f64 = 17.026549105;
/// Monoisotopic mass of phosphoric acid, from the shared chemistry.
pub const H3PO4_MASS: f64 = 97.976895625;

/// Ions whose distances to a peak differ by less than this (Th) are equally near. Isobaric
/// ions (`b2` and `b4^2` of one peptide) differ by at most a few 1e-6 Th of floating-point
/// rounding; distinct ions of tryptic peptides come within 1e-4 Th of each other.
pub const TIE_TOLERANCE: f64 = 1e-5;

/// Ion-type slots per cleavage position of the ladder.
pub const TYPES_PER_POSITION: usize = 4;
/// The highest fragment charge the ladder carries.
pub const MAX_FRAGMENT_CHARGE: u8 = 2;

/// Monoisotopic residue mass of one of the 20 standard amino acids or selenocysteine, as the
/// C# computes it from `AminoAcidFormulas.Default` (the one-residue peptide less water).
/// `None` for the ambiguity codes (B, J, O, X, Z).
pub fn residue_mass(aa: char) -> Option<f64> {
    Some(match aa {
        'A' => 71.03711380499999,
        'C' => 103.009184505,
        'D' => 115.02694306499998,
        'E' => 129.04259313499998,
        'F' => 147.068413945,
        'G' => 57.021463735,
        'H' => 137.05891187499998,
        'I' => 113.084064015,
        'K' => 128.09496305,
        'L' => 113.084064015,
        'M' => 131.040484645,
        'N' => 114.04292747,
        'P' => 97.05276387499998,
        'Q' => 128.05857754000002,
        'R' => 156.10111105,
        'S' => 87.032028435,
        'T' => 101.04767850500001,
        'V' => 99.068413945,
        'W' => 186.07931298,
        'Y' => 163.063328575,
        'U' => 150.953635605,
        _ => return None,
    })
}

/// Modification mass by 0-based residue position. Modifications at one position add: both
/// library loaders put an N-terminal modification on residue 0, where a modification of the
/// first residue also sits, and every ion spanning that residue carries both.
pub fn mod_masses_by_position(modifications: &[Modification]) -> HashMap<usize, f64> {
    let mut masses: HashMap<usize, f64> = HashMap::new();
    for m in modifications {
        let mass = masses.entry(m.position).or_insert(0.0);
        *mass += m.mass_delta;
    }
    masses
}

/// m/z of the b or y ion of `ordinal` residues at `charge`, or `None` for another ion type, an
/// ordinal below 1 or past the end of the sequence, a charge of 0, or an ion spanning a residue
/// with no standard mass. The summation order is the C#'s, so the two agree to the bit.
pub fn calculate_fragment_mz(
    ion_type: IonType,
    ordinal: usize,
    charge: u8,
    sequence: &[char],
    mod_masses: &HashMap<usize, f64>,
    neutral_loss: Option<f64>,
) -> Option<f64> {
    if charge == 0 || ordinal < 1 {
        return None;
    }
    let seq_len = sequence.len();
    let (start, end) = match ion_type {
        IonType::B => (0, ordinal),
        IonType::Y => (seq_len.checked_sub(ordinal)?, seq_len),
        _ => return None,
    };
    if end > seq_len {
        return None;
    }

    let mut mass = 0.0;
    for (i, aa) in sequence.iter().enumerate().take(end).skip(start) {
        mass += residue_mass(*aa)?;
        if let Some(&mod_mass) = mod_masses.get(&i) {
            mass += mod_mass;
        }
    }

    match ion_type {
        IonType::B => mass += PROTON_MASS,
        IonType::Y => mass += H2O_MASS + PROTON_MASS,
        _ => {}
    }

    if let Some(loss) = neutral_loss {
        mass -= loss;
    }

    Some((mass + (charge as f64 - 1.0) * PROTON_MASS) / charge as f64)
}

/// Number of ladder slots for a peptide of `length` residues: four per cleavage position
/// (b and y at charge 1 and 2).
pub fn slot_count(length: usize) -> usize {
    if length < 2 {
        0
    } else {
        TYPES_PER_POSITION * (length - 1)
    }
}

/// Ion type of a ladder slot: slot `p*4+t` holds `t` in b1, b2, y1, y2.
pub fn ion_type_of(slot: usize) -> IonType {
    if slot % TYPES_PER_POSITION < 2 {
        IonType::B
    } else {
        IonType::Y
    }
}

/// Fragment charge of a ladder slot.
pub fn charge_of(slot: usize) -> u8 {
    if slot % 2 == 0 {
        1
    } else {
        2
    }
}

/// Ordinal (residues in the ion) of a ladder slot in a peptide of `length` residues:
/// position `p` holds b(p+1) and y(L-1-p).
pub fn ordinal_of(slot: usize, length: usize) -> usize {
    let position = slot / TYPES_PER_POSITION;
    if ion_type_of(slot) == IonType::B {
        position + 1
    } else {
        length - 1 - position
    }
}

/// m/z of every ladder slot for `sequence` carrying `modifications` at `precursor_charge`,
/// NaN where a slot is not applicable (a fragment charge above the precursor's). Empty for a
/// sequence shorter than two residues.
pub fn build_ladder(
    sequence: &[char],
    modifications: &[Modification],
    precursor_charge: u8,
) -> Vec<f64> {
    let length = sequence.len();
    let count = slot_count(length);
    if count == 0 {
        return Vec::new();
    }
    let mod_masses = mod_masses_by_position(modifications);
    let max_charge = precursor_charge.clamp(1, MAX_FRAGMENT_CHARGE);
    (0..count)
        .map(|slot| {
            let charge = charge_of(slot);
            if charge > max_charge {
                return f64::NAN;
            }
            calculate_fragment_mz(
                ion_type_of(slot),
                ordinal_of(slot, length),
                charge,
                sequence,
                &mod_masses,
                None,
            )
            .unwrap_or(f64::NAN)
        })
        .collect()
}

/// `|lib_mz - obs_mz| <= tolerance at lib_mz`, the C#'s `WithinTolerance` (a ppm tolerance is
/// taken at the library m/z, as extraction takes it).
fn within_tolerance(tolerance: &FragmentToleranceConfig, lib_mz: f64, obs_mz: f64) -> bool {
    (lib_mz - obs_mz).abs() <= tolerance.tolerance_da(lib_mz)
}

/// The primary b and y ions of one peptide precursor that [`type_fragments`] types peaks as:
/// the applicable ladder slots, sorted by m/z so each peak is placed by binary search.
pub struct FragmentCandidates {
    length: usize,
    slot_count: usize,
    sorted: Vec<(f64, usize)>,
}

impl FragmentCandidates {
    pub fn new(sequence: &str, modifications: &[Modification], precursor_charge: u8) -> Self {
        let chars: Vec<char> = sequence.chars().collect();
        let ladder = build_ladder(&chars, modifications, precursor_charge);
        let mut sorted: Vec<(f64, usize)> = ladder
            .iter()
            .enumerate()
            .filter(|(_, mz)| !mz.is_nan())
            .map(|(slot, mz)| (*mz, slot))
            .collect();
        // Slots are unique, so the order is total.
        sorted.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        Self {
            length: chars.len(),
            slot_count: ladder.len(),
            sorted,
        }
    }

    /// Slots of the ladder, applicable or not.
    pub fn slot_count(&self) -> usize {
        self.slot_count
    }

    /// The slot of the preferred unclaimed candidate within tolerance of the peak, or `None`:
    /// the nearest, taking those within [`TIE_TOLERANCE`] of it as equally near; then the lower
    /// fragment charge; then y before b; then the shorter ion.
    pub fn preferred_slot(
        &self,
        peak_mz: f64,
        tolerance: &FragmentToleranceConfig,
        claimed: &[bool],
    ) -> Option<usize> {
        if self.sorted.is_empty() || !peak_mz.is_finite() {
            return None;
        }
        let reach = tolerance.tolerance_da(peak_mz) * 1.01;
        let mut in_reach: Vec<(f64, usize)> = Vec::new();
        let mut i = self.sorted.partition_point(|c| c.0 < peak_mz - reach);
        while i < self.sorted.len() && self.sorted[i].0 <= peak_mz + reach {
            let (mz, slot) = self.sorted[i];
            if !claimed[slot] && within_tolerance(tolerance, mz, peak_mz) {
                in_reach.push((mz, slot));
            }
            i += 1;
        }
        let nearest = in_reach
            .iter()
            .map(|(mz, _)| (mz - peak_mz).abs())
            .fold(f64::INFINITY, f64::min);
        in_reach
            .iter()
            .filter(|(mz, _)| (mz - peak_mz).abs() - nearest < TIE_TOLERANCE)
            .map(|(_, slot)| *slot)
            .min_by_key(|&slot| {
                (
                    charge_of(slot),
                    if ion_type_of(slot) == IonType::Y {
                        0
                    } else {
                        1
                    },
                    ordinal_of(slot, self.length),
                )
            })
    }

    /// The annotation a ladder slot names.
    pub fn annotation_of(&self, slot: usize) -> FragmentAnnotation {
        FragmentAnnotation {
            ion_type: ion_type_of(slot),
            ordinal: ordinal_of(slot, self.length) as u8,
            charge: charge_of(slot),
            neutral_loss: None,
        }
    }
}

/// Types `fragments` of the peptide `sequence` carrying `modifications` at `precursor_charge`,
/// in place, by the rules in the module summary. A typed peak gets the ion's type, ordinal and
/// charge with no loss; any other is `Unknown`, whatever the library stated. Returns the number
/// of peaks typed.
pub fn type_fragments(
    sequence: &str,
    modifications: &[Modification],
    precursor_charge: u8,
    fragments: &mut [LibraryFragment],
    tolerance: &FragmentToleranceConfig,
) -> usize {
    if fragments.is_empty() {
        return 0;
    }
    let candidates = FragmentCandidates::new(sequence, modifications, precursor_charge);
    let mut claimed = vec![false; candidates.slot_count()];
    // Most intense first, so the ion a peak claims goes to the strongest peak it explains;
    // ties go to the lower m/z rather than the earlier row.
    let mut order: Vec<usize> = (0..fragments.len()).collect();
    order.sort_by(|&a, &b| {
        cmp_f32_dotnet(
            fragments[b].relative_intensity,
            fragments[a].relative_intensity,
        )
        .then(cmp_f64_dotnet(fragments[a].mz, fragments[b].mz))
        .then(a.cmp(&b))
    });
    let mut typed = 0;
    for i in order {
        match candidates.preferred_slot(fragments[i].mz, tolerance, &claimed) {
            Some(slot) => {
                fragments[i].annotation = candidates.annotation_of(slot);
                claimed[slot] = true;
                typed += 1;
            }
            None => {
                fragments[i].annotation = FragmentAnnotation {
                    ion_type: IonType::Unknown,
                    ordinal: 0,
                    charge: 0,
                    neutral_loss: None,
                };
            }
        }
    }
    typed
}

/// Merges peaks of `fragments_by_mz` (in m/z order) that share one m/z into one peak with their
/// summed intensity, and scales the result so the base peak is 1. The merged peak keeps the
/// first peak's annotation. Intensities are summed and scaled in single precision, as the C#
/// does.
pub fn merge_equal_mz(fragments_by_mz: Vec<LibraryFragment>) -> Vec<LibraryFragment> {
    let mut merged: Vec<LibraryFragment> = Vec::with_capacity(fragments_by_mz.len());
    for fragment in fragments_by_mz {
        match merged.last_mut() {
            Some(last) if last.mz == fragment.mz => {
                last.relative_intensity += fragment.relative_intensity;
            }
            _ => merged.push(fragment),
        }
    }
    let mut base_peak = 0.0f32;
    for fragment in &merged {
        if fragment.relative_intensity > base_peak {
            base_peak = fragment.relative_intensity;
        }
    }
    if base_peak > 0.0 && base_peak != 1.0 {
        for fragment in &mut merged {
            fragment.relative_intensity /= base_peak;
        }
    }
    merged
}

/// Puts one entry's fragments in m/z order, one peak per m/z scaled to the base peak, and
/// types them within `tolerance`. Returns the number of peaks typed.
pub fn prepare_library_fragments(
    entry: &mut LibraryEntry,
    tolerance: &FragmentToleranceConfig,
) -> usize {
    let mut fragments = std::mem::take(&mut entry.fragments);
    fragments.sort_by(|a, b| cmp_f64_dotnet(a.mz, b.mz));
    let mut fragments = merge_equal_mz(fragments);
    let typed = type_fragments(
        &entry.sequence,
        &entry.modifications,
        entry.charge,
        &mut fragments,
        tolerance,
    );
    entry.fragments = fragments;
    typed
}

/// [`prepare_library_fragments`] over a whole library. Returns (peaks typed, peaks in all).
pub fn prepare_library(
    entries: &mut [LibraryEntry],
    tolerance: &FragmentToleranceConfig,
) -> (usize, usize) {
    let mut typed = 0;
    let mut total = 0;
    for entry in entries.iter_mut() {
        typed += prepare_library_fragments(entry, tolerance);
        total += entry.fragments.len();
    }
    (typed, total)
}

/// Indices of the `count` most intense of `fragments`, most intense first, or every index in
/// library order when there are at most `count`. Ties on intensity are broken by the
/// fragments themselves - lower m/z first, then ion type, ordinal and charge - never by their
/// position in the library. The one ranking every top-N selection uses.
pub fn top_fragment_indices(fragments: &[LibraryFragment], count: usize) -> Vec<usize> {
    if fragments.len() <= count {
        return (0..fragments.len()).collect();
    }
    let mut order: Vec<usize> = (0..fragments.len()).collect();
    order.sort_by(|&a, &b| {
        let fa = &fragments[a];
        let fb = &fragments[b];
        cmp_f32_dotnet(fb.relative_intensity, fa.relative_intensity)
            .then(cmp_f64_dotnet(fa.mz, fb.mz))
            .then(ion_type_rank(fa.annotation.ion_type).cmp(&ion_type_rank(fb.annotation.ion_type)))
            .then(fa.annotation.ordinal.cmp(&fb.annotation.ordinal))
            .then(fa.annotation.charge.cmp(&fb.annotation.charge))
    });
    order.truncate(count);
    order
}

/// The ion types in the C# enum's order, which its tie-break compares.
fn ion_type_rank(ion_type: IonType) -> u8 {
    match ion_type {
        IonType::B => 0,
        IonType::Y => 1,
        IonType::A => 2,
        IonType::C => 3,
        IonType::X => 4,
        IonType::Z => 5,
        IonType::Precursor => 6,
        IonType::Internal => 7,
        IonType::Immonium => 8,
        IonType::Unknown => 9,
    }
}

/// .NET's `Comparer<double>`: NaN is less than every number and equal to itself.
pub fn cmp_f64_dotnet(a: f64, b: f64) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        (false, false) => a.partial_cmp(&b).unwrap_or(Ordering::Equal),
    }
}

/// .NET's `Comparer<float>`: NaN is less than every number and equal to itself.
pub fn cmp_f32_dotnet(a: f32, b: f32) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Less,
        (false, true) => Ordering::Greater,
        (false, false) => a.partial_cmp(&b).unwrap_or(Ordering::Equal),
    }
}

/// The unit of a tolerance as one byte, for the library cache header.
pub fn tolerance_unit_code(unit: ToleranceUnit) -> u8 {
    match unit {
        ToleranceUnit::Ppm => 0,
        ToleranceUnit::Mz => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQUENCE: &str = "PEPTIDEK";
    // The first twelve residues repeat the first six, so b6 and b12^2 share one m/z.
    const REPEAT: &str = "AAAAAAAAAAAAAAK";

    fn mz(ion_type: IonType, ordinal: usize, charge: u8, sequence: &str) -> f64 {
        let chars: Vec<char> = sequence.chars().collect();
        calculate_fragment_mz(ion_type, ordinal, charge, &chars, &HashMap::new(), None).unwrap()
    }

    fn fragment(mz: f64, intensity: f32) -> LibraryFragment {
        LibraryFragment {
            mz,
            relative_intensity: intensity,
            annotation: FragmentAnnotation::default(),
        }
    }

    fn entry(sequence: &str, fragments: Vec<LibraryFragment>) -> LibraryEntry {
        let mut entry = LibraryEntry::new(
            1,
            sequence.to_string(),
            sequence.to_string(),
            2,
            500.0,
            10.0,
        );
        entry.fragments = fragments;
        entry
    }

    #[test]
    fn residue_masses_are_the_shared_chemistry_values() {
        assert_eq!(residue_mass('C'), Some(103.009184505));
        assert_eq!(residue_mass('T'), Some(101.04767850500001));
        assert_eq!(residue_mass('U'), Some(150.953635605));
        assert_eq!(residue_mass('X'), None);
        assert_eq!(H2O_MASS, 18.0105647);
        assert_eq!(NH3_MASS, 17.026549105);
    }

    #[test]
    fn ladder_slots_hold_b_and_y_at_charge_1_and_2() {
        let chars: Vec<char> = SEQUENCE.chars().collect();
        let ladder = build_ladder(&chars, &[], 2);
        assert_eq!(ladder.len(), 28);
        assert_eq!(ladder[0], mz(IonType::B, 1, 1, SEQUENCE));
        assert_eq!(ladder[1], mz(IonType::B, 1, 2, SEQUENCE));
        assert_eq!(ladder[2], mz(IonType::Y, 7, 1, SEQUENCE));
        assert_eq!(ladder[3], mz(IonType::Y, 7, 2, SEQUENCE));
        // At precursor charge 1 the doubly charged slots are not applicable.
        let singly = build_ladder(&chars, &[], 1);
        assert!(singly[1].is_nan());
        assert!(!singly[0].is_nan());
    }

    #[test]
    fn peaks_take_the_nearest_primary_ion_within_tolerance() {
        let y3 = mz(IonType::Y, 3, 1, SEQUENCE);
        let b4 = mz(IonType::B, 4, 1, SEQUENCE);
        let tolerance = FragmentToleranceConfig::default(); // 10 ppm
        let mut fragments = vec![
            fragment(y3, 1.0),
            fragment(y3 + y3 * 9.9e-6, 0.5),
            fragment(y3 + y3 * 10.1e-6, 0.4),
            fragment(b4, 0.3),
            fragment(700.0, 0.2),
        ];
        let typed = type_fragments(SEQUENCE, &[], 2, &mut fragments, &tolerance);
        assert_eq!(typed, 2);
        assert_eq!(fragments[0].annotation.ion_type, IonType::Y);
        assert_eq!(fragments[0].annotation.ordinal, 3);
        // Each ion types one peak: the most intense peak within reach of y3 took it, and the
        // second, with no other candidate in reach, is Unknown; the third is outside 10 ppm.
        assert_eq!(fragments[1].annotation.ion_type, IonType::Unknown);
        assert_eq!(fragments[2].annotation.ion_type, IonType::Unknown);
        assert_eq!(fragments[3].annotation.ion_type, IonType::B);
        assert_eq!(fragments[3].annotation.ordinal, 4);
        assert_eq!(fragments[4].annotation.ion_type, IonType::Unknown);
        assert_eq!(fragments[4].annotation.ordinal, 0);
    }

    #[test]
    fn isobaric_ions_prefer_the_lower_charge() {
        let b6 = mz(IonType::B, 6, 1, REPEAT);
        assert!((b6 - mz(IonType::B, 12, 2, REPEAT)).abs() < TIE_TOLERANCE);
        let mut fragments = vec![fragment(b6, 1.0)];
        type_fragments(
            REPEAT,
            &[],
            2,
            &mut fragments,
            &FragmentToleranceConfig::default(),
        );
        assert_eq!(fragments[0].annotation.ion_type, IonType::B);
        assert_eq!(fragments[0].annotation.ordinal, 6);
        assert_eq!(fragments[0].annotation.charge, 1);
    }

    #[test]
    fn equal_mz_peaks_merge_and_the_spectrum_rescales() {
        let merged = merge_equal_mz(vec![
            fragment(100.0, 0.5),
            fragment(200.0, 0.2),
            fragment(200.0, 1.0),
            fragment(300.0, 0.4),
        ]);
        assert_eq!(merged.len(), 3);
        let expected: Vec<f32> = vec![0.5f32 / 1.2f32, 1.0, 0.4f32 / 1.2f32];
        let actual: Vec<f32> = merged.iter().map(|f| f.relative_intensity).collect();
        assert_eq!(actual, expected);

        // Distinct m/z, however close, stay separate, and a base peak of 1 leaves the rest.
        let separate = merge_equal_mz(vec![fragment(100.0, 0.5), fragment(100.0 + 1e-9, 1.0)]);
        assert_eq!(separate.len(), 2);
        assert_eq!(separate[0].relative_intensity, 0.5);
    }

    #[test]
    fn preparation_sorts_merges_and_types() {
        let b6 = mz(IonType::B, 6, 1, REPEAT);
        let y3 = mz(IonType::Y, 3, 1, REPEAT);
        let mut entry = entry(
            REPEAT,
            vec![fragment(y3, 1.0), fragment(b6, 1.0), fragment(b6, 1.0)],
        );
        let typed = prepare_library_fragments(&mut entry, &FragmentToleranceConfig::default());
        assert_eq!(typed, 2);
        assert_eq!(entry.fragments.len(), 2);
        assert_eq!(entry.fragments[0].mz, b6.min(y3));
        let b6_peak = entry.fragments.iter().find(|f| f.mz == b6).unwrap();
        assert_eq!(b6_peak.relative_intensity, 1.0);
        assert_eq!(b6_peak.annotation.ordinal, 6);
        let y3_peak = entry.fragments.iter().find(|f| f.mz == y3).unwrap();
        assert_eq!(y3_peak.relative_intensity, 0.5);
    }

    #[test]
    fn top_indices_break_ties_by_mz_then_ion() {
        let fragments: Vec<LibraryFragment> = [
            (300.0, 1.0),
            (350.0, 0.9),
            (400.0, 0.8),
            (450.0, 0.5),
            (500.0, 0.5),
            (550.0, 0.5),
            (600.0, 0.5),
            (650.0, 0.1),
        ]
        .iter()
        .map(|&(mz, intensity)| fragment(mz, intensity))
        .collect();
        assert_eq!(top_fragment_indices(&fragments, 6), vec![0, 1, 2, 3, 4, 5]);
        let mut reversed = fragments.clone();
        reversed.reverse();
        let mut chosen: Vec<f64> = top_fragment_indices(&reversed, 6)
            .iter()
            .map(|&i| reversed[i].mz)
            .collect();
        chosen.sort_by(f64::total_cmp);
        assert_eq!(chosen, vec![300.0, 350.0, 400.0, 450.0, 500.0, 550.0]);
        // At most `count` fragments come back in library order.
        assert_eq!(top_fragment_indices(&fragments[..3], 6), vec![0, 1, 2]);

        // Two peaks at one m/z and intensity: the ion decides, in either order.
        let b = LibraryFragment {
            mz: 225.6266,
            relative_intensity: 0.5,
            annotation: FragmentAnnotation {
                ion_type: IonType::B,
                ordinal: 4,
                charge: 2,
                neutral_loss: None,
            },
        };
        let y = LibraryFragment {
            mz: 225.6266,
            relative_intensity: 0.5,
            annotation: FragmentAnnotation {
                ion_type: IonType::Y,
                ordinal: 3,
                charge: 2,
                neutral_loss: None,
            },
        };
        for pair in [vec![b.clone(), y.clone()], vec![y, b]] {
            let top = top_fragment_indices(&pair, 1);
            assert_eq!(pair[top[0]].annotation.ion_type, IonType::B);
        }
    }
}
