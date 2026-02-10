//! Isotope distribution calculations from exact elemental composition
//!
//! This module calculates theoretical isotope distributions based on the exact
//! elemental composition of peptides (C, H, N, O, S), rather than using the
//! averagine approximation.
//!
//! The isotope distribution is calculated using the polynomial expansion method,
//! accounting for natural isotope abundances of each element.

/// Elemental composition (atom counts)
#[derive(Debug, Clone, Default)]
pub struct ElementalComposition {
    /// Carbon atoms
    pub c: usize,
    /// Hydrogen atoms
    pub h: usize,
    /// Nitrogen atoms
    pub n: usize,
    /// Oxygen atoms
    pub o: usize,
    /// Sulfur atoms
    pub s: usize,
}

impl ElementalComposition {
    /// Create new elemental composition
    pub fn new(c: usize, h: usize, n: usize, o: usize, s: usize) -> Self {
        Self { c, h, n, o, s }
    }

    /// Add two compositions
    pub fn add(&self, other: &ElementalComposition) -> ElementalComposition {
        ElementalComposition {
            c: self.c + other.c,
            h: self.h + other.h,
            n: self.n + other.n,
            o: self.o + other.o,
            s: self.s + other.s,
        }
    }

    /// Subtract water (for peptide bond formation)
    pub fn subtract_water(&mut self) {
        self.h = self.h.saturating_sub(2);
        self.o = self.o.saturating_sub(1);
    }

    /// Add water (for N and C terminus)
    pub fn add_water(&mut self) {
        self.h += 2;
        self.o += 1;
    }

    /// Calculate monoisotopic mass
    pub fn monoisotopic_mass(&self) -> f64 {
        // Exact monoisotopic masses from NIST
        const C12: f64 = 12.000000;
        const H1: f64 = 1.007825;
        const N14: f64 = 14.003074;
        const O16: f64 = 15.994915;
        const S32: f64 = 31.972071;

        self.c as f64 * C12
            + self.h as f64 * H1
            + self.n as f64 * N14
            + self.o as f64 * O16
            + self.s as f64 * S32
    }
}

/// Get amino acid elemental composition (residue form, without water)
///
/// These are the compositions for amino acid residues in a peptide chain,
/// NOT the free amino acid form (which includes an extra H2O).
pub fn amino_acid_composition(aa: char) -> Option<ElementalComposition> {
    // Amino acid residue compositions (formula - H2O for peptide bond)
    // Reference: NIST Chemistry WebBook
    let comp = match aa.to_ascii_uppercase() {
        'A' => ElementalComposition::new(3, 5, 1, 1, 0),   // Alanine: C3H5NO
        'C' => ElementalComposition::new(3, 5, 1, 1, 1),   // Cysteine: C3H5NOS
        'D' => ElementalComposition::new(4, 5, 1, 3, 0),   // Aspartic acid: C4H5NO3
        'E' => ElementalComposition::new(5, 7, 1, 3, 0),   // Glutamic acid: C5H7NO3
        'F' => ElementalComposition::new(9, 9, 1, 1, 0),   // Phenylalanine: C9H9NO
        'G' => ElementalComposition::new(2, 3, 1, 1, 0),   // Glycine: C2H3NO
        'H' => ElementalComposition::new(6, 7, 3, 1, 0),   // Histidine: C6H7N3O
        'I' => ElementalComposition::new(6, 11, 1, 1, 0),  // Isoleucine: C6H11NO
        'K' => ElementalComposition::new(6, 12, 2, 1, 0),  // Lysine: C6H12N2O
        'L' => ElementalComposition::new(6, 11, 1, 1, 0),  // Leucine: C6H11NO
        'M' => ElementalComposition::new(5, 9, 1, 1, 1),   // Methionine: C5H9NOS
        'N' => ElementalComposition::new(4, 6, 2, 2, 0),   // Asparagine: C4H6N2O2
        'P' => ElementalComposition::new(5, 7, 1, 1, 0),   // Proline: C5H7NO
        'Q' => ElementalComposition::new(5, 8, 2, 2, 0),   // Glutamine: C5H8N2O2
        'R' => ElementalComposition::new(6, 12, 4, 1, 0),  // Arginine: C6H12N4O
        'S' => ElementalComposition::new(3, 5, 1, 2, 0),   // Serine: C3H5NO2
        'T' => ElementalComposition::new(4, 7, 1, 2, 0),   // Threonine: C4H7NO2
        'V' => ElementalComposition::new(5, 9, 1, 1, 0),   // Valine: C5H9NO
        'W' => ElementalComposition::new(11, 10, 2, 1, 0), // Tryptophan: C11H10N2O
        'Y' => ElementalComposition::new(9, 9, 1, 2, 0),   // Tyrosine: C9H9NO2
        'U' => ElementalComposition::new(3, 5, 1, 1, 0),   // Selenocysteine (approx as Cys-S)
        _ => return None,
    };
    Some(comp)
}

/// Calculate elemental composition for a peptide sequence
///
/// Accounts for peptide bond formation (loss of H2O per bond)
/// and adds terminal H2O (H at N-term, OH at C-term).
pub fn peptide_composition(sequence: &str) -> Option<ElementalComposition> {
    let mut total = ElementalComposition::default();

    // Sum up all amino acid residue compositions
    for aa in sequence.chars() {
        if aa == '[' || aa == ']' || aa == '(' || aa == ')' || aa.is_ascii_digit()
            || aa == '+' || aa == '-' || aa == '.' {
            // Skip modification notation characters
            continue;
        }
        if let Some(comp) = amino_acid_composition(aa) {
            total = total.add(&comp);
        } else if aa.is_ascii_alphabetic() {
            // Unknown amino acid
            return None;
        }
    }

    // Add terminal H2O (H at N-term, OH at C-term)
    total.add_water();

    Some(total)
}

/// Natural isotope abundances
/// Reference: IUPAC 2016 atomic weights and isotopic compositions,
/// except 13C which uses D.E. Matthews value from IDCalc
mod abundances {
    // Carbon: 13C abundance from D.E. Matthews (IDCalc)
    // IDCalc uses 12C = 100, 13C = 1.0958793 (relative)
    // Converted: 13C / (100 + 1.0958793) = 0.01084
    pub const C13: f64 = 0.01084;

    // Hydrogen: 2H (deuterium) abundance (0.0115%)
    pub const H2: f64 = 0.000115;

    // Nitrogen: 15N abundance (0.364%)
    pub const N15: f64 = 0.00364;

    // Oxygen: 16O (99.757%), 17O (0.038%), 18O (0.205%)
    pub const O16: f64 = 0.99757;
    pub const O17: f64 = 0.00038;
    pub const O18: f64 = 0.00205;

    // Sulfur: 32S (94.99%), 33S (0.75%), 34S (4.25%), 36S (0.01%)
    pub const S32: f64 = 0.9499;
    pub const S33: f64 = 0.0075;
    pub const S34: f64 = 0.0425;
    pub const S36: f64 = 0.0001;
}

/// Calculate isotope distribution for an elemental composition
///
/// Uses the polynomial expansion method, considering up to M+4 peaks.
/// Returns normalized relative intensities for [M+0, M+1, M+2, M+3, M+4].
///
/// For efficiency, uses the binomial approximation for C and H (which dominate),
/// and direct calculation for N, O, S.
pub fn calculate_isotope_distribution(composition: &ElementalComposition) -> [f64; 5] {
    // For very large molecules, the isotope distribution approaches a normal distribution
    // centered near the average mass. For peptides (typically <5000 Da), we use
    // direct polynomial expansion.

    // Start with M+0 = 1.0, then convolve with each element's contribution
    let mut dist = [1.0, 0.0, 0.0, 0.0, 0.0];

    // Carbon contribution (13C = 1.07%)
    if composition.c > 0 {
        dist = convolve_binomial(&dist, composition.c, abundances::C13);
    }

    // Hydrogen contribution (2H = 0.0115%)
    // Usually negligible but included for accuracy
    if composition.h > 0 {
        dist = convolve_binomial(&dist, composition.h, abundances::H2);
    }

    // Nitrogen contribution (15N = 0.364%)
    if composition.n > 0 {
        dist = convolve_binomial(&dist, composition.n, abundances::N15);
    }

    // Oxygen contribution (17O = 0.038%, 18O = 0.205%)
    // More complex due to three isotopes
    if composition.o > 0 {
        dist = convolve_oxygen(&dist, composition.o);
    }

    // Sulfur contribution (33S = 0.75%, 34S = 4.25%)
    // Significant for Cys/Met-containing peptides
    if composition.s > 0 {
        dist = convolve_sulfur(&dist, composition.s);
    }

    // Normalize to sum = 1
    let sum: f64 = dist.iter().sum();
    if sum > 0.0 {
        for v in &mut dist {
            *v /= sum;
        }
    }

    dist
}

/// Convolve distribution with binomial (two-isotope element like C, H, N)
fn convolve_binomial(dist: &[f64; 5], n: usize, heavy_prob: f64) -> [f64; 5] {
    let light_prob = 1.0 - heavy_prob;
    let mut new_dist = [0.0; 5];

    // Calculate binomial coefficients and probabilities for 0..=4 heavy atoms
    for k in 0..=4.min(n) {
        let binom_coeff = binomial_coefficient(n, k);
        let prob = binom_coeff * heavy_prob.powi(k as i32) * light_prob.powi((n - k) as i32);

        // Shift and add
        for i in 0..5 {
            if i + k < 5 {
                new_dist[i + k] += dist[i] * prob;
            }
        }
    }

    new_dist
}

/// Convolve with oxygen isotope distribution (three isotopes)
fn convolve_oxygen(dist: &[f64; 5], n: usize) -> [f64; 5] {
    let mut new_dist = [0.0; 5];

    // For oxygen: 16O (99.757%), 17O (+1, 0.038%), 18O (+2, 0.205%)
    // Iterate over possible counts of 17O and 18O
    for n17 in 0..=n.min(4) {
        for n18 in 0..=(n - n17).min(4) {
            let n16 = n - n17 - n18;
            let mass_shift = n17 + 2 * n18;
            if mass_shift > 4 {
                continue;
            }

            // Multinomial probability
            let prob = multinomial_prob(
                n,
                &[n16, n17, n18],
                &[abundances::O16, abundances::O17, abundances::O18],
            );

            // Convolve
            for i in 0..5 {
                if i + mass_shift < 5 {
                    new_dist[i + mass_shift] += dist[i] * prob;
                }
            }
        }
    }

    new_dist
}

/// Convolve with sulfur isotope distribution (four isotopes)
fn convolve_sulfur(dist: &[f64; 5], n: usize) -> [f64; 5] {
    let mut new_dist = [0.0; 5];

    // For sulfur: 32S, 33S (+1), 34S (+2), 36S (+4)
    // Iterate over possible counts
    for n33 in 0..=n.min(4) {
        for n34 in 0..=(n - n33).min(4) {
            for n36 in 0..=(n - n33 - n34).min(1) {
                // 36S is very rare
                let n32 = n - n33 - n34 - n36;
                let mass_shift = n33 + 2 * n34 + 4 * n36;
                if mass_shift > 4 {
                    continue;
                }

                // Multinomial probability
                let prob = multinomial_prob(
                    n,
                    &[n32, n33, n34, n36],
                    &[abundances::S32, abundances::S33, abundances::S34, abundances::S36],
                );

                // Convolve
                for i in 0..5 {
                    if i + mass_shift < 5 {
                        new_dist[i + mass_shift] += dist[i] * prob;
                    }
                }
            }
        }
    }

    new_dist
}

/// Calculate binomial coefficient n choose k
fn binomial_coefficient(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    if k == 0 || k == n {
        return 1.0;
    }

    // Use the multiplicative formula to avoid overflow
    let k = k.min(n - k); // Use smaller k for efficiency
    let mut result = 1.0;
    for i in 0..k {
        result *= (n - i) as f64 / (i + 1) as f64;
    }
    result
}

/// Calculate multinomial probability
fn multinomial_prob(n: usize, counts: &[usize], probs: &[f64]) -> f64 {
    if counts.iter().sum::<usize>() != n {
        return 0.0;
    }

    // Multinomial coefficient
    let mut coeff = 1.0;
    let mut remaining = n;
    for &k in counts {
        coeff *= binomial_coefficient(remaining, k);
        remaining -= k;
    }

    // Probability product
    let mut prob = coeff;
    for (&count, &p) in counts.iter().zip(probs.iter()) {
        prob *= p.powi(count as i32);
    }

    prob
}

/// Calculate cosine angle score between observed and theoretical isotope distributions
///
/// # Arguments
/// * `observed` - Observed isotope intensities [M-1, M+0, M+1, M+2, M+3]
/// * `theoretical` - Theoretical relative intensities [M+0, M+1, M+2, M+3, M+4]
///
/// # Returns
/// Cosine similarity (0.0 to 1.0), or None if vectors are invalid
///
/// Note: The observed array has M-1 first (from pyXcorrDIA convention),
/// while theoretical starts at M+0. This function aligns them appropriately.
pub fn isotope_cosine_score(observed: &[f64; 5], theoretical: &[f64; 5]) -> Option<f64> {
    // Align: observed[1..5] corresponds to M+0 to M+3
    // theoretical[0..4] corresponds to M+0 to M+3
    let obs = &observed[1..5]; // Skip M-1
    let theo = &theoretical[0..4]; // M+0 to M+3

    // Calculate dot product and norms
    let mut dot = 0.0;
    let mut obs_norm_sq = 0.0;
    let mut theo_norm_sq = 0.0;

    for i in 0..4 {
        dot += obs[i] * theo[i];
        obs_norm_sq += obs[i] * obs[i];
        theo_norm_sq += theo[i] * theo[i];
    }

    let obs_norm = obs_norm_sq.sqrt();
    let theo_norm = theo_norm_sq.sqrt();

    if obs_norm < 1e-10 || theo_norm < 1e-10 {
        return None;
    }

    let cosine = dot / (obs_norm * theo_norm);
    Some(cosine.clamp(0.0, 1.0))
}

/// Calculate isotope cosine score for a peptide sequence
///
/// # Arguments
/// * `sequence` - Peptide sequence (unmodified amino acids only)
/// * `observed` - Observed isotope intensities [M-1, M+0, M+1, M+2, M+3]
///
/// # Returns
/// Cosine similarity score (0.0 to 1.0), or None if calculation fails
pub fn peptide_isotope_cosine(sequence: &str, observed: &[f64; 5]) -> Option<f64> {
    let composition = peptide_composition(sequence)?;
    let theoretical = calculate_isotope_distribution(&composition);
    isotope_cosine_score(observed, &theoretical)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies elemental compositions for Ala, Cys, Met and that unknown amino acids return None.
    #[test]
    fn test_amino_acid_compositions() {
        // Test a few amino acids
        let ala = amino_acid_composition('A').unwrap();
        assert_eq!(ala.c, 3);
        assert_eq!(ala.h, 5);
        assert_eq!(ala.n, 1);
        assert_eq!(ala.o, 1);
        assert_eq!(ala.s, 0);

        let cys = amino_acid_composition('C').unwrap();
        assert_eq!(cys.s, 1);

        let met = amino_acid_composition('M').unwrap();
        assert_eq!(met.s, 1);

        // Unknown amino acid
        assert!(amino_acid_composition('X').is_none());
    }

    /// Verifies peptide composition includes terminal H2O and sums residue atom counts correctly.
    #[test]
    fn test_peptide_composition() {
        // Simple peptide: AA (Ala-Ala)
        let comp = peptide_composition("AA").unwrap();
        // Two Ala residues: 2 * C3H5NO = C6H10N2O2
        // Plus terminal H2O: C6H12N2O3
        assert_eq!(comp.c, 6);
        assert_eq!(comp.h, 12);
        assert_eq!(comp.n, 2);
        assert_eq!(comp.o, 3);

        // Peptide with sulfur: CM (Cys-Met)
        let comp = peptide_composition("CM").unwrap();
        assert_eq!(comp.s, 2); // One from Cys, one from Met
    }

    /// Verifies binomial coefficient calculation for known values of C(5,k) and C(10,3).
    #[test]
    fn test_binomial_coefficient() {
        assert!((binomial_coefficient(5, 0) - 1.0).abs() < 1e-10);
        assert!((binomial_coefficient(5, 1) - 5.0).abs() < 1e-10);
        assert!((binomial_coefficient(5, 2) - 10.0).abs() < 1e-10);
        assert!((binomial_coefficient(5, 3) - 10.0).abs() < 1e-10);
        assert!((binomial_coefficient(5, 5) - 1.0).abs() < 1e-10);
        assert!((binomial_coefficient(10, 3) - 120.0).abs() < 1e-10);
    }

    /// Verifies that a small peptide (GGG) has dominant M+0 peak and distribution sums to 1.
    #[test]
    fn test_isotope_distribution_small_peptide() {
        // Small peptide: GGG (Gly-Gly-Gly)
        // C6H11N3O4 (3 Gly residues + terminal H2O)
        let comp = peptide_composition("GGG").unwrap();
        let dist = calculate_isotope_distribution(&comp);

        // M+0 should be highest
        assert!(dist[0] > dist[1]);
        assert!(dist[0] > 0.9); // Small molecule, mostly monoisotopic

        // Sum should be 1
        let sum: f64 = dist.iter().sum();
        assert!((sum - 1.0).abs() < 0.001);
    }

    /// Verifies that a larger peptide (ACDEFGHIK) has significant M+1 and distribution sums to 1.
    #[test]
    fn test_isotope_distribution_larger_peptide() {
        // Larger peptide: ACDEFGHIK (9 residues)
        let comp = peptide_composition("ACDEFGHIK").unwrap();
        let dist = calculate_isotope_distribution(&comp);

        // For larger peptides, M+1 becomes more significant
        assert!(dist[0] > 0.4); // M+0 still significant
        assert!(dist[1] > 0.2); // M+1 significant

        // Sum should be 1
        let sum: f64 = dist.iter().sum();
        assert!((sum - 1.0).abs() < 0.001);
    }

    /// Verifies that identical observed and theoretical distributions yield a cosine score of 1.0.
    #[test]
    fn test_isotope_cosine_perfect_match() {
        // Perfect match: theoretical distribution matches observed
        let theoretical = [0.6, 0.3, 0.08, 0.02, 0.0];
        // observed has M-1 prefix: [M-1, M+0, M+1, M+2, M+3]
        let observed = [0.0, 0.6, 0.3, 0.08, 0.02];

        let score = isotope_cosine_score(&observed, &theoretical).unwrap();
        assert!((score - 1.0).abs() < 0.001);
    }

    /// Verifies that orthogonal observed and theoretical distributions yield a near-zero score.
    #[test]
    fn test_isotope_cosine_orthogonal() {
        // Completely different distributions
        let theoretical = [1.0, 0.0, 0.0, 0.0, 0.0];
        let observed = [0.0, 0.0, 0.0, 1.0, 0.0]; // All in M+2 position

        let score = isotope_cosine_score(&observed, &theoretical).unwrap();
        assert!(score < 0.1); // Should be near zero
    }

    /// Verifies end-to-end peptide isotope cosine scoring with a self-matching distribution.
    #[test]
    fn test_peptide_isotope_cosine() {
        // Test end-to-end calculation
        let sequence = "PEPTIDE";
        let comp = peptide_composition(sequence).unwrap();
        let theoretical = calculate_isotope_distribution(&comp);

        // Create observed that matches theoretical
        let observed = [
            0.0,              // M-1
            theoretical[0],   // M+0
            theoretical[1],   // M+1
            theoretical[2],   // M+2
            theoretical[3],   // M+3
        ];

        let score = peptide_isotope_cosine(sequence, &observed).unwrap();
        assert!((score - 1.0).abs() < 0.001);
    }

    /// Verifies that sulfur-containing peptides have higher M+2 ratios due to 34S abundance.
    #[test]
    fn test_sulfur_effect() {
        // Compare peptides with and without sulfur
        let no_sulfur = peptide_composition("AAAA").unwrap();
        let with_sulfur = peptide_composition("AAMC").unwrap(); // Met and Cys

        let dist_no_s = calculate_isotope_distribution(&no_sulfur);
        let dist_with_s = calculate_isotope_distribution(&with_sulfur);

        // Sulfur-containing peptide should have relatively higher M+2
        // due to 34S (4.25%)
        let m2_ratio_no_s = dist_no_s[2] / dist_no_s[0];
        let m2_ratio_with_s = dist_with_s[2] / dist_with_s[0];

        assert!(m2_ratio_with_s > m2_ratio_no_s);
    }
}
