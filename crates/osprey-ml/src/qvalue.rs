//! Originally from Sage (https://github.com/lazear/sage)
//! Copyright (c) 2022 Michael Lazear
//! Licensed under the MIT License

//! Q-value calculation using target-decoy competition
//!
//! # Invariants
//! * Input must be sorted in descending order by score (best matches first)

/// Calculate q-values for a sorted list of target and decoy matches
///
/// # Arguments
/// * `is_decoy` - Boolean slice where `true` indicates a decoy match
/// * `q_values` - Output slice to store calculated q-values
///
/// # Returns
/// Number of matches passing 1% FDR threshold
///
/// # Panics
/// Panics if `is_decoy.len() != q_values.len()`
pub fn calculate_q_values(is_decoy: &[bool], q_values: &mut [f64]) -> usize {
    assert_eq!(
        is_decoy.len(),
        q_values.len(),
        "is_decoy and q_values must have same length"
    );

    let mut decoy = 0; // Count actual decoys (not conservative +1)
    let mut target = 0;

    // Calculate FDR for each position
    for (i, &is_dec) in is_decoy.iter().enumerate() {
        match is_dec {
            true => decoy += 1,
            false => target += 1,
        }
        q_values[i] = decoy as f64 / target.max(1) as f64;
    }

    // Reverse pass to calculate cumulative minimum (q-value)
    let mut q_min = 1.0f64;
    let mut passing = 0;
    for q in q_values.iter_mut().rev() {
        q_min = q_min.min(*q);
        *q = q_min;
        if q_min <= 0.01 {
            passing += 1;
        }
    }

    passing
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_qvalue_calculation() {
        // 5 targets, 2 decoys (all targets before decoys)
        let is_decoy = vec![false, false, false, false, false, true, true];
        let mut q_values = vec![0.0; 7];

        let passing = calculate_q_values(&is_decoy, &mut q_values);

        // Expected FDR progression: 0/1, 0/2, 0/3, 0/4, 0/5, 1/5, 2/5
        // Q-values (cumulative min from end): 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.4
        assert_eq!(q_values[0], 0.0); // 0/1 = 0.0
        assert_eq!(q_values[4], 0.0); // 0/5 = 0.0
        assert_eq!(q_values[5], 0.2); // 1/5 = 0.2
        assert_eq!(q_values[6], 0.4); // 2/5 = 0.4

        // All 5 targets pass 1% FDR (q = 0.0 <= 0.01)
        assert_eq!(passing, 5);
    }

    #[test]
    fn test_qvalue_mixed() {
        // Interleaved targets and decoys
        let is_decoy = vec![false, false, true, false, false, true, false];
        let mut q_values = vec![0.0; 7];

        calculate_q_values(&is_decoy, &mut q_values);

        // Verify q-values are monotonically non-increasing (cumulative minimum property)
        for i in 1..q_values.len() {
            assert!(
                q_values[i] >= q_values[i - 1] || (q_values[i] - q_values[i - 1]).abs() < 1e-10,
                "Q-values should be monotonically non-increasing: q[{}]={} > q[{}]={}",
                i - 1,
                q_values[i - 1],
                i,
                q_values[i]
            );
        }

        // First q-value should be <= 1.0 (worst case: all decoys)
        assert!(q_values[0] <= 1.0);
    }
}
