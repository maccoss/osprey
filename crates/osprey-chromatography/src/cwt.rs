//! Continuous Wavelet Transform (CWT) consensus peak detection for DIA proteomics.
//!
//! Uses the Mexican Hat (Ricker) wavelet as a matched filter for Gaussian-like
//! chromatographic peaks. Multi-transition consensus — the pointwise median of
//! CWT coefficients across all fragment XICs — is only high where the majority
//! of transitions simultaneously exhibit peak-like shapes, naturally rejecting
//! single-fragment interference.
//!
//! Reference: DIA-PeakDetectionPlan.md

use log;
use osprey_core::XICPeakBounds;
use std::f64::consts::PI;

/// Generate a discrete Mexican Hat (Ricker) wavelet kernel.
///
/// The Mexican Hat wavelet is the negative normalized second derivative of a
/// Gaussian, acting as a matched filter for Gaussian-like chromatographic peaks:
///
/// `ψ(t) = (2 / sqrt(3σ) π^(1/4)) × (1 - (t/σ)²) × exp(-t²/(2σ²))`
///
/// The kernel is zero-mean corrected to ensure pure wavelet behavior (zero
/// response to constant signals).
///
/// # Arguments
/// * `sigma` - Scale parameter in scan-index units. Controls matched peak width.
/// * `kernel_radius` - Points on each side of center. Total size = 2*radius + 1.
///   Recommended: `ceil(5 * sigma)` to capture >99.99% of wavelet energy.
pub fn mexican_hat_kernel(sigma: f64, kernel_radius: usize) -> Vec<f64> {
    let len = 2 * kernel_radius + 1;
    let mut kernel = vec![0.0; len];
    let center = kernel_radius as f64;

    let norm = 2.0 / ((3.0 * sigma).sqrt() * PI.powf(0.25));

    for (i, val) in kernel.iter_mut().enumerate() {
        let t = (i as f64 - center) / sigma;
        *val = norm * (1.0 - t * t) * (-0.5 * t * t).exp();
    }

    // Zero-mean correction: numerical discretization can leave a tiny DC offset
    let mean = kernel.iter().sum::<f64>() / len as f64;
    for v in &mut kernel {
        *v -= mean;
    }

    kernel
}

/// Convolve a signal with a kernel using "same" output size.
///
/// Output length equals input length. Edges are zero-padded: points beyond
/// the signal boundary contribute 0. Direct convolution is used because
/// typical XIC lengths (30-200 scans) and kernel sizes (11-51 points) make
/// O(N×K) trivially fast.
fn convolve_same(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let k = kernel.len();
    if n == 0 || k == 0 {
        return vec![0.0; n];
    }

    let half_k = k / 2;
    let mut output = vec![0.0; n];

    for (i, out) in output.iter_mut().enumerate() {
        let mut sum = 0.0;
        for (j, &kval) in kernel.iter().enumerate() {
            let signal_idx = i as isize + j as isize - half_k as isize;
            if signal_idx >= 0 && (signal_idx as usize) < n {
                sum += signal[signal_idx as usize] * kval;
            }
        }
        *out = sum;
    }

    output
}

/// Estimate the CWT scale parameter (sigma) from fragment XICs.
///
/// Estimates sigma from the median FWHM across fragment XICs that have
/// detectable signal. The FWHM is measured in scan-index units and converted
/// via `sigma = fwhm / 2.355` (Gaussian FWHM = 2.355σ).
///
/// Falls back to 4.0 scans if no FWHM can be estimated.
pub fn estimate_cwt_scale(xics: &[(usize, Vec<(f64, f64)>)]) -> f64 {
    let mut fwhm_values: Vec<f64> = Vec::new();

    for (_, xic) in xics {
        if xic.len() < 5 {
            continue;
        }

        let intensities: Vec<f64> = xic.iter().map(|(_, v)| *v).collect();

        // Find apex
        let (apex_idx, apex_val) = intensities
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, &v)| (i, v))
            .unwrap_or((0, 0.0));

        if apex_val <= 0.0 {
            continue;
        }

        let half_max = apex_val / 2.0;

        // Walk left to find half-height crossing
        let mut left_idx_f: Option<f64> = None;
        for i in (1..=apex_idx).rev() {
            if intensities[i] >= half_max && intensities[i - 1] < half_max {
                let denom = intensities[i] - intensities[i - 1];
                if denom > 1e-30 {
                    let frac = (intensities[i] - half_max) / denom;
                    left_idx_f = Some(i as f64 - frac);
                } else {
                    left_idx_f = Some(i as f64);
                }
                break;
            }
        }

        // Walk right to find half-height crossing
        let mut right_idx_f: Option<f64> = None;
        for i in apex_idx..intensities.len() - 1 {
            if intensities[i] >= half_max && intensities[i + 1] < half_max {
                let denom = intensities[i] - intensities[i + 1];
                if denom > 1e-30 {
                    let frac = (intensities[i] - half_max) / denom;
                    right_idx_f = Some(i as f64 + frac);
                } else {
                    right_idx_f = Some(i as f64);
                }
                break;
            }
        }

        if let (Some(l), Some(r)) = (left_idx_f, right_idx_f) {
            let fwhm = r - l;
            if fwhm > 1.0 {
                fwhm_values.push(fwhm);
            }
        }
    }

    if fwhm_values.is_empty() {
        return 4.0; // Default fallback
    }

    // Median FWHM
    fwhm_values.sort_by(|a, b| a.total_cmp(b));
    let median_fwhm = if fwhm_values.len() % 2 == 0 {
        let mid = fwhm_values.len() / 2;
        (fwhm_values[mid - 1] + fwhm_values[mid]) / 2.0
    } else {
        fwhm_values[fwhm_values.len() / 2]
    };

    // sigma = FWHM / 2.355 (Gaussian relationship)
    let sigma = median_fwhm / 2.355;
    sigma.clamp(2.0, 20.0)
}

/// Compute the median of a small slice of f64 values.
///
/// Sorts the values and returns the middle element (or average of two middle
/// elements for even-length slices). NaN values are treated as negative infinity
/// for sorting purposes.
fn small_median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.total_cmp(b));
    let n = values.len();
    if n % 2 == 0 {
        (values[n / 2 - 1] + values[n / 2]) / 2.0
    } else {
        values[n / 2]
    }
}

/// Detect peaks using multi-transition CWT consensus.
///
/// Drop-in replacement for `detect_all_xic_peaks`. Integrates multi-transition
/// information directly into peak detection by:
///
/// 1. Convolving each fragment XIC with the Mexican Hat wavelet (1D, along time)
/// 2. Computing the pointwise **median** of CWT coefficients across all transitions
/// 3. Finding local maxima in the consensus CWT signal
/// 4. Defining boundaries by walking from each apex to **zero-crossings**
/// 5. Computing area/SNR from a reference signal (sum of fragment XICs)
///
/// The consensus is only high where the majority of transitions simultaneously
/// exhibit peak-like shapes. Single-fragment interference is naturally rejected
/// by the median.
///
/// # Arguments
/// * `xics` - Fragment XICs: `Vec<(fragment_index, Vec<(RT, intensity)>)>`.
///   All XICs must share the same time axis (same length, same RT values).
/// * `min_consensus_height` - Minimum consensus CWT coefficient to consider
///   a peak candidate. Use 0.0 to accept any positive consensus peak.
///
/// # Returns
/// `Vec<XICPeakBounds>` sorted by consensus CWT coefficient descending.
/// The `apex_intensity`, `area`, and `signal_to_noise` come from the summed
/// raw fragment signal, not from CWT coefficients.
pub fn detect_cwt_consensus_peaks(
    xics: &[(usize, Vec<(f64, f64)>)],
    min_consensus_height: f64,
) -> Vec<XICPeakBounds> {
    // Validate input
    if xics.len() < 2 {
        return Vec::new();
    }

    let n_scans = xics[0].1.len();
    if n_scans < 5 {
        return Vec::new();
    }

    // Verify all XICs have the same length
    if xics.iter().any(|(_, xic)| xic.len() != n_scans) {
        return Vec::new();
    }

    let n_frags = xics.len();

    // Extract RT axis and intensity matrix
    let rts: Vec<f64> = xics[0].1.iter().map(|(rt, _)| *rt).collect();
    let intensities: Vec<Vec<f64>> = xics
        .iter()
        .map(|(_, xic)| xic.iter().map(|(_, v)| *v).collect())
        .collect();

    // Check for any signal
    let has_signal = intensities.iter().any(|row| row.iter().any(|&v| v > 0.0));
    if !has_signal {
        return Vec::new();
    }

    // Estimate scale and generate kernel
    let sigma = estimate_cwt_scale(xics);
    let kernel_radius = ((5.0 * sigma).ceil() as usize).min(n_scans / 2);
    let kernel = mexican_hat_kernel(sigma, kernel_radius);

    // Convolve each fragment XIC with the Mexican Hat kernel
    let cwt_coeffs: Vec<Vec<f64>> = intensities
        .iter()
        .map(|signal| convolve_same(signal, &kernel))
        .collect();

    // Compute pointwise median consensus
    let mut consensus = vec![0.0; n_scans];
    let mut buf = vec![0.0; n_frags];
    for s in 0..n_scans {
        for f in 0..n_frags {
            buf[f] = cwt_coeffs[f][s];
        }
        consensus[s] = small_median(&mut buf);
    }

    // Build reference signal (sum of raw fragment intensities) for quantitation
    let mut ref_signal = vec![0.0; n_scans];
    for s in 0..n_scans {
        for row in &intensities {
            ref_signal[s] += row[s];
        }
    }

    // Find local maxima in consensus signal
    let mut apex_indices: Vec<(usize, f64)> = Vec::new();

    for i in 1..n_scans - 1 {
        if consensus[i] > min_consensus_height
            && consensus[i] > consensus[i - 1]
            && consensus[i] > consensus[i + 1]
        {
            apex_indices.push((i, consensus[i]));
        }
    }

    // Check endpoints
    if n_scans >= 2 && consensus[0] > min_consensus_height && consensus[0] > consensus[1] {
        apex_indices.push((0, consensus[0]));
    }
    if n_scans >= 2
        && consensus[n_scans - 1] > min_consensus_height
        && consensus[n_scans - 1] > consensus[n_scans - 2]
    {
        apex_indices.push((n_scans - 1, consensus[n_scans - 1]));
    }

    // Sort by consensus coefficient descending
    apex_indices.sort_by(|a, b| b.1.total_cmp(&a.1));

    // Coverage factor: zero-crossings are at ±σ (~68% coverage).
    // Multiply by this factor to get ~95% coverage (±2σ for Gaussian peaks).
    let coverage_factor = 2.0;

    // Build peaks with extended boundaries (±2σ with valley guard)
    let mut peaks: Vec<XICPeakBounds> = Vec::new();

    for &(apex_idx, _consensus_coeff) in &apex_indices {
        // Step 1: Find zero-crossings (at ±σ from apex in CWT space)
        let mut left_zc = apex_idx;
        while left_zc > 0 && consensus[left_zc - 1] > 0.0 {
            left_zc -= 1;
        }

        let mut right_zc = apex_idx;
        while right_zc < n_scans - 1 && consensus[right_zc + 1] > 0.0 {
            right_zc += 1;
        }

        // Step 2: Compute asymmetric σ estimates from zero-crossing distances
        let left_sigma = apex_idx - left_zc;
        let right_sigma = right_zc - apex_idx;

        // Step 3: Extend to ±2σ for ~95% coverage
        let target_start = apex_idx.saturating_sub(coverage_factor as usize * left_sigma.max(1));
        let target_end =
            (apex_idx + coverage_factor as usize * right_sigma.max(1)).min(n_scans - 1);

        // Step 4: Valley guard — walk outward from zero-crossing on the raw
        // reference signal. If the signal rises >5% of apex above a running
        // minimum, we've entered a neighboring peak; stop at the valley.
        let ref_apex_val_approx = ref_signal[apex_idx];
        let rise_threshold = ref_apex_val_approx * 0.05;

        // Extend left with valley guard
        let mut start_idx = left_zc;
        let mut running_min = ref_signal[left_zc];
        let mut min_pos = left_zc;
        for i in (target_start..left_zc).rev() {
            if ref_signal[i] < running_min {
                running_min = ref_signal[i];
                min_pos = i;
            } else if ref_signal[i] - running_min > rise_threshold {
                // Rising into another peak — stop at the valley
                start_idx = min_pos;
                break;
            }
            start_idx = i;
        }

        // Extend right with valley guard
        let mut end_idx = right_zc;
        let mut running_min = ref_signal[right_zc];
        let mut min_pos = right_zc;
        for (i, &val) in ref_signal
            .iter()
            .enumerate()
            .take(target_end + 1)
            .skip(right_zc + 1)
        {
            if val < running_min {
                running_min = val;
                min_pos = i;
            } else if val - running_min > rise_threshold {
                // Rising into another peak — stop at the valley
                end_idx = min_pos;
                break;
            }
            end_idx = i;
        }

        // Require at least 3 scans
        if end_idx - start_idx + 1 < 3 {
            continue;
        }

        // Find the actual apex in the reference signal within boundaries
        let (ref_apex_idx, ref_apex_val) = ref_signal[start_idx..=end_idx]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, &v)| (start_idx + i, v))
            .unwrap_or((apex_idx, 0.0));

        // Compute area from reference signal
        let ref_series: Vec<(f64, f64)> = (start_idx..=end_idx)
            .map(|i| (rts[i], ref_signal[i]))
            .collect();
        let area = super::trapezoidal_area(&ref_series);

        // Compute SNR from reference signal
        let snr = super::compute_snr(&ref_signal, ref_apex_idx, start_idx, end_idx);

        peaks.push(XICPeakBounds {
            apex_rt: rts[ref_apex_idx],
            apex_intensity: ref_apex_val,
            apex_index: ref_apex_idx,
            start_rt: rts[start_idx],
            end_rt: rts[end_idx],
            start_index: start_idx,
            end_index: end_idx,
            area,
            signal_to_noise: snr,
        });
    }

    peaks
}

/// Diagnostic version of `detect_cwt_consensus_peaks` that logs detailed
/// information about the peak detection process. Call this for specific
/// precursors you want to investigate.
pub fn detect_cwt_consensus_peaks_diagnostic(
    xics: &[(usize, Vec<(f64, f64)>)],
    min_consensus_height: f64,
    label: &str,
) -> Vec<XICPeakBounds> {
    // Validate input
    if xics.len() < 2 {
        log::info!(
            "[CWT-DIAG {}] Skipped: only {} XICs (need ≥2)",
            label,
            xics.len()
        );
        return Vec::new();
    }

    let n_scans = xics[0].1.len();
    if n_scans < 5 {
        log::info!(
            "[CWT-DIAG {}] Skipped: only {} scans (need ≥5)",
            label,
            n_scans
        );
        return Vec::new();
    }

    if xics.iter().any(|(_, xic)| xic.len() != n_scans) {
        log::info!("[CWT-DIAG {}] Skipped: XICs have different lengths", label);
        return Vec::new();
    }

    let n_frags = xics.len();

    let rts: Vec<f64> = xics[0].1.iter().map(|(rt, _)| *rt).collect();
    let intensities: Vec<Vec<f64>> = xics
        .iter()
        .map(|(_, xic)| xic.iter().map(|(_, v)| *v).collect())
        .collect();

    let has_signal = intensities.iter().any(|row| row.iter().any(|&v| v > 0.0));
    if !has_signal {
        log::info!("[CWT-DIAG {}] Skipped: no signal", label);
        return Vec::new();
    }

    // Scale estimation
    let sigma = estimate_cwt_scale(xics);
    let kernel_radius = ((5.0 * sigma).ceil() as usize).min(n_scans / 2);
    let kernel = mexican_hat_kernel(sigma, kernel_radius);

    log::info!(
        "[CWT-DIAG {}] n_scans={}, n_frags={}, sigma={:.2} scans, kernel_radius={}, RT range=[{:.3}, {:.3}]",
        label, n_scans, n_frags, sigma, kernel_radius, rts[0], rts[n_scans - 1]
    );

    // Log per-fragment max intensities
    for (i, (frag_idx, _)) in xics.iter().enumerate() {
        let max_int: f64 = intensities[i].iter().cloned().fold(0.0f64, f64::max);
        let max_pos = intensities[i]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(j, _)| j)
            .unwrap_or(0);
        log::info!(
            "[CWT-DIAG {}]   frag[{}] idx={}: max_intensity={:.0} at scan {} (RT={:.3})",
            label,
            i,
            frag_idx,
            max_int,
            max_pos,
            rts[max_pos]
        );
    }

    // CWT convolution
    let cwt_coeffs: Vec<Vec<f64>> = intensities
        .iter()
        .map(|signal| convolve_same(signal, &kernel))
        .collect();

    // Consensus
    let mut consensus = vec![0.0; n_scans];
    let mut buf = vec![0.0; n_frags];
    for s in 0..n_scans {
        for f in 0..n_frags {
            buf[f] = cwt_coeffs[f][s];
        }
        consensus[s] = small_median(&mut buf);
    }

    // Reference signal
    let mut ref_signal = vec![0.0; n_scans];
    for s in 0..n_scans {
        for row in &intensities {
            ref_signal[s] += row[s];
        }
    }

    // Log consensus signal around the peak region (where consensus > 0)
    let positive_scans: Vec<usize> = (0..n_scans).filter(|&s| consensus[s] > 0.0).collect();
    if !positive_scans.is_empty() {
        let first = positive_scans[0].saturating_sub(2);
        let last = (positive_scans[positive_scans.len() - 1] + 3).min(n_scans);
        log::info!(
            "[CWT-DIAG {}] Consensus signal (scans {}..{}):",
            label,
            first,
            last
        );
        for s in first..last {
            log::info!(
                "[CWT-DIAG {}]   scan={:3} RT={:.3} consensus={:12.1} ref_signal={:12.0}",
                label,
                s,
                rts[s],
                consensus[s],
                ref_signal[s]
            );
        }
    }

    // Find local maxima
    let mut apex_indices: Vec<(usize, f64)> = Vec::new();
    for i in 1..n_scans - 1 {
        if consensus[i] > min_consensus_height
            && consensus[i] > consensus[i - 1]
            && consensus[i] > consensus[i + 1]
        {
            apex_indices.push((i, consensus[i]));
        }
    }
    if n_scans >= 2 && consensus[0] > min_consensus_height && consensus[0] > consensus[1] {
        apex_indices.push((0, consensus[0]));
    }
    if n_scans >= 2
        && consensus[n_scans - 1] > min_consensus_height
        && consensus[n_scans - 1] > consensus[n_scans - 2]
    {
        apex_indices.push((n_scans - 1, consensus[n_scans - 1]));
    }
    apex_indices.sort_by(|a, b| b.1.total_cmp(&a.1));

    log::info!(
        "[CWT-DIAG {}] Found {} consensus apex candidates",
        label,
        apex_indices.len()
    );
    for (rank, &(idx, coeff)) in apex_indices.iter().enumerate() {
        log::info!(
            "[CWT-DIAG {}]   candidate {}: scan={} RT={:.3} consensus_coeff={:.1}",
            label,
            rank,
            idx,
            rts[idx],
            coeff
        );
    }

    let coverage_factor = 2.0;
    let mut peaks: Vec<XICPeakBounds> = Vec::new();

    for (rank, &(apex_idx, consensus_coeff)) in apex_indices.iter().enumerate() {
        // Zero-crossings
        let mut left_zc = apex_idx;
        while left_zc > 0 && consensus[left_zc - 1] > 0.0 {
            left_zc -= 1;
        }
        let mut right_zc = apex_idx;
        while right_zc < n_scans - 1 && consensus[right_zc + 1] > 0.0 {
            right_zc += 1;
        }

        let left_sigma = apex_idx - left_zc;
        let right_sigma = right_zc - apex_idx;

        let target_start = apex_idx.saturating_sub(coverage_factor as usize * left_sigma.max(1));
        let target_end =
            (apex_idx + coverage_factor as usize * right_sigma.max(1)).min(n_scans - 1);

        log::info!(
            "[CWT-DIAG {}] Peak {} (scan={}, RT={:.3}, coeff={:.1}):",
            label,
            rank,
            apex_idx,
            rts[apex_idx],
            consensus_coeff
        );
        log::info!(
            "[CWT-DIAG {}]   zero-crossings: left_zc={} (RT={:.3}), right_zc={} (RT={:.3})",
            label,
            left_zc,
            rts[left_zc],
            right_zc,
            rts[right_zc]
        );
        log::info!(
            "[CWT-DIAG {}]   sigma estimates: left_sigma={} scans, right_sigma={} scans",
            label,
            left_sigma,
            right_sigma
        );
        log::info!(
            "[CWT-DIAG {}]   2σ targets: target_start={} (RT={:.3}), target_end={} (RT={:.3})",
            label,
            target_start,
            rts[target_start],
            target_end,
            rts[target_end]
        );

        // Valley guard
        let ref_apex_val_approx = ref_signal[apex_idx];
        let rise_threshold = ref_apex_val_approx * 0.05;

        log::info!(
            "[CWT-DIAG {}]   valley guard: ref_apex={:.0}, rise_threshold={:.0} (5%)",
            label,
            ref_apex_val_approx,
            rise_threshold
        );

        // Left extension
        let mut start_idx = left_zc;
        let mut running_min = ref_signal[left_zc];
        let mut min_pos = left_zc;
        let mut left_guard_triggered = false;
        for i in (target_start..left_zc).rev() {
            if ref_signal[i] < running_min {
                running_min = ref_signal[i];
                min_pos = i;
            } else if ref_signal[i] - running_min > rise_threshold {
                start_idx = min_pos;
                left_guard_triggered = true;
                log::info!(
                    "[CWT-DIAG {}]   LEFT valley guard triggered at scan={} (RT={:.3}): \
                     ref_signal={:.0}, running_min={:.0}, rise={:.0} > threshold={:.0}. \
                     Stopped at valley scan={} (RT={:.3})",
                    label,
                    i,
                    rts[i],
                    ref_signal[i],
                    running_min,
                    ref_signal[i] - running_min,
                    rise_threshold,
                    min_pos,
                    rts[min_pos]
                );
                break;
            }
            start_idx = i;
        }
        if !left_guard_triggered {
            log::info!(
                "[CWT-DIAG {}]   LEFT extension: reached target_start={} (RT={:.3}), no valley guard",
                label, start_idx, rts[start_idx]
            );
        }

        // Right extension
        let mut end_idx = right_zc;
        let mut running_min = ref_signal[right_zc];
        let mut min_pos = right_zc;
        let mut right_guard_triggered = false;
        for (i, &val) in ref_signal
            .iter()
            .enumerate()
            .take(target_end + 1)
            .skip(right_zc + 1)
        {
            if val < running_min {
                running_min = val;
                min_pos = i;
            } else if val - running_min > rise_threshold {
                end_idx = min_pos;
                right_guard_triggered = true;
                log::info!(
                    "[CWT-DIAG {}]   RIGHT valley guard triggered at scan={} (RT={:.3}): \
                     ref_signal={:.0}, running_min={:.0}, rise={:.0} > threshold={:.0}. \
                     Stopped at valley scan={} (RT={:.3})",
                    label,
                    i,
                    rts[i],
                    val,
                    running_min,
                    val - running_min,
                    rise_threshold,
                    min_pos,
                    rts[min_pos]
                );
                break;
            }
            end_idx = i;
        }
        if !right_guard_triggered {
            log::info!(
                "[CWT-DIAG {}]   RIGHT extension: reached target_end={} (RT={:.3}), no valley guard",
                label, end_idx, rts[end_idx]
            );
        }

        if end_idx - start_idx + 1 < 3 {
            log::info!(
                "[CWT-DIAG {}]   REJECTED: only {} scans (need ≥3)",
                label,
                end_idx - start_idx + 1
            );
            continue;
        }

        let (ref_apex_idx, ref_apex_val) = ref_signal[start_idx..=end_idx]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, &v)| (start_idx + i, v))
            .unwrap_or((apex_idx, 0.0));

        let ref_series: Vec<(f64, f64)> = (start_idx..=end_idx)
            .map(|i| (rts[i], ref_signal[i]))
            .collect();
        let area = super::trapezoidal_area(&ref_series);
        let snr = super::compute_snr(&ref_signal, ref_apex_idx, start_idx, end_idx);

        log::info!(
            "[CWT-DIAG {}]   FINAL: apex_scan={} (RT={:.3}), boundaries=[{}, {}] (RT=[{:.3}, {:.3}]), \
             width={:.3} min, area={:.0}, S/N={:.1}",
            label, ref_apex_idx, rts[ref_apex_idx],
            start_idx, end_idx, rts[start_idx], rts[end_idx],
            rts[end_idx] - rts[start_idx], area, snr
        );

        peaks.push(XICPeakBounds {
            apex_rt: rts[ref_apex_idx],
            apex_intensity: ref_apex_val,
            apex_index: ref_apex_idx,
            start_rt: rts[start_idx],
            end_rt: rts[end_idx],
            start_index: start_idx,
            end_index: end_idx,
            area,
            signal_to_noise: snr,
        });
    }

    log::info!("[CWT-DIAG {}] Returning {} peaks", label, peaks.len());

    peaks
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: create a Gaussian XIC centered at `center` with given sigma and amplitude.
    fn make_gaussian_xic(
        center: f64,
        sigma: f64,
        amplitude: f64,
        n: usize,
        spacing: f64,
    ) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let rt = i as f64 * spacing;
                let intensity =
                    amplitude * (-((rt - center).powi(2)) / (2.0 * sigma * sigma)).exp();
                (rt, intensity)
            })
            .collect()
    }

    // =========================================================================
    // Kernel tests
    // =========================================================================

    #[test]
    fn test_kernel_zero_mean() {
        let kernel = mexican_hat_kernel(5.0, 25);
        let sum: f64 = kernel.iter().sum();
        assert!(
            sum.abs() < 1e-10,
            "Kernel should be zero-mean, got sum = {:.2e}",
            sum
        );
    }

    #[test]
    fn test_kernel_symmetric() {
        let kernel = mexican_hat_kernel(5.0, 25);
        let len = kernel.len();
        for i in 0..len / 2 {
            assert!(
                (kernel[i] - kernel[len - 1 - i]).abs() < 1e-12,
                "Kernel should be symmetric: kernel[{}]={:.6} != kernel[{}]={:.6}",
                i,
                kernel[i],
                len - 1 - i,
                kernel[len - 1 - i]
            );
        }
    }

    #[test]
    fn test_kernel_positive_center_negative_tails() {
        let kernel = mexican_hat_kernel(5.0, 25);
        let center = kernel.len() / 2;
        assert!(
            kernel[center] > 0.0,
            "Kernel center should be positive: {}",
            kernel[center]
        );
        assert!(
            kernel[0] < 0.0,
            "Kernel tail should be negative: {}",
            kernel[0]
        );
        assert!(
            kernel[kernel.len() - 1] < 0.0,
            "Kernel tail should be negative: {}",
            kernel[kernel.len() - 1]
        );
    }

    #[test]
    fn test_kernel_size() {
        let kernel = mexican_hat_kernel(3.0, 15);
        assert_eq!(kernel.len(), 31); // 2 * 15 + 1
    }

    // =========================================================================
    // Convolution tests
    // =========================================================================

    #[test]
    fn test_convolve_same_length() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let kernel = vec![0.25, 0.5, 0.25];
        let result = convolve_same(&signal, &kernel);
        assert_eq!(result.len(), signal.len());
    }

    #[test]
    fn test_convolve_delta_function() {
        // Delta function convolved with kernel should produce the kernel
        let n = 21;
        let center = n / 2;
        let mut signal = vec![0.0; n];
        signal[center] = 1.0;
        let kernel = mexican_hat_kernel(3.0, 5);

        let result = convolve_same(&signal, &kernel);

        // Result at center should match kernel center
        let k_center = kernel.len() / 2;
        assert!(
            (result[center] - kernel[k_center]).abs() < 1e-10,
            "Delta response at center: got {:.6}, expected {:.6}",
            result[center],
            kernel[k_center]
        );
    }

    #[test]
    fn test_convolve_gaussian_response() {
        // Mexican Hat convolved with Gaussian should produce positive response at center
        let sigma = 5.0;
        let signal: Vec<f64> = (0..100)
            .map(|i| {
                let t = i as f64 - 50.0;
                1000.0 * (-t * t / (2.0 * sigma * sigma)).exp()
            })
            .collect();
        let kernel = mexican_hat_kernel(sigma, (5.0 * sigma) as usize);

        let result = convolve_same(&signal, &kernel);

        // Response should be positive at peak center
        assert!(
            result[50] > 0.0,
            "CWT of Gaussian should be positive at center: {}",
            result[50]
        );
        // Response should be maximal near the center
        let max_idx = result
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(i, _)| i)
            .unwrap();
        assert!(
            (max_idx as i32 - 50).unsigned_abs() <= 1,
            "Max CWT response should be near center, got index {}",
            max_idx
        );
    }

    // =========================================================================
    // Scale estimation tests
    // =========================================================================

    #[test]
    fn test_estimate_scale_known_peak() {
        // Create XICs with known FWHM. Gaussian FWHM = 2.355 * sigma.
        // Use sigma=5 scans → FWHM ≈ 11.775 scans → estimated CWT sigma ≈ 5
        let sigma: f64 = 5.0;
        let sigma_rt = sigma * 0.01; // sigma in RT units (0.01 min spacing)
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let amplitude = 1000.0 * (i + 1) as f64;
                let xic: Vec<(f64, f64)> = (0..100)
                    .map(|s| {
                        let rt = s as f64 * 0.01; // 0.01 min spacing
                        let intensity =
                            amplitude * (-(rt - 0.5).powi(2) / (2.0 * sigma_rt * sigma_rt)).exp();
                        (rt, intensity)
                    })
                    .collect();
                (i, xic)
            })
            .collect();

        let est_sigma = estimate_cwt_scale(&xics);
        // Should be close to 5.0 scans (sigma in scan units)
        assert!(
            est_sigma > 3.0 && est_sigma < 8.0,
            "Estimated sigma should be near 5.0, got {:.2}",
            est_sigma
        );
    }

    #[test]
    fn test_estimate_scale_fallback() {
        // All-zero XICs → should return fallback
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let xic: Vec<(f64, f64)> = (0..50).map(|s| (s as f64 * 0.1, 0.0)).collect();
                (i, xic)
            })
            .collect();

        let est = estimate_cwt_scale(&xics);
        assert!(
            (est - 4.0).abs() < 1e-10,
            "Should fall back to 4.0, got {}",
            est
        );
    }

    // =========================================================================
    // Consensus peak detection tests
    // =========================================================================

    #[test]
    fn test_consensus_single_gaussian_peak() {
        // 6 coeluting Gaussian peaks at RT=5.0, different amplitudes
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let amplitude = 1000.0 * (i + 1) as f64;
                (i, make_gaussian_xic(5.0, 0.3, amplitude, 100, 0.1))
            })
            .collect();

        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);

        assert!(!peaks.is_empty(), "Should find at least 1 peak");
        assert!(
            (peaks[0].apex_rt - 5.0).abs() < 0.3,
            "Peak apex should be near 5.0, got {:.2}",
            peaks[0].apex_rt
        );
        assert!(
            peaks[0].start_rt < 5.0 && peaks[0].end_rt > 5.0,
            "Peak boundaries should bracket apex: [{:.2}, {:.2}]",
            peaks[0].start_rt,
            peaks[0].end_rt
        );
    }

    #[test]
    fn test_consensus_interference_rejection() {
        // 5 fragments with real peak at RT=5.0, 1 fragment with interference at RT=2.0
        let mut xics: Vec<(usize, Vec<(f64, f64)>)> = (0..5)
            .map(|i| (i, make_gaussian_xic(5.0, 0.3, 1000.0, 100, 0.1)))
            .collect();
        // Interference fragment: peak at RT=2.0, no signal at RT=5.0
        xics.push((5, make_gaussian_xic(2.0, 0.3, 5000.0, 100, 0.1)));

        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);

        assert!(!peaks.is_empty(), "Should find at least 1 peak");
        // The first (strongest) peak should be the real one at RT=5.0
        assert!(
            (peaks[0].apex_rt - 5.0).abs() < 0.5,
            "Best peak should be near real peak at 5.0, got {:.2}",
            peaks[0].apex_rt
        );
        // Interference at RT=2.0 should NOT be the strongest peak
        // (only 1/6 fragments have signal there, median suppresses it)
    }

    #[test]
    fn test_consensus_two_separated_peaks() {
        // 6 fragments each with two peaks at RT=3.0 and RT=7.0
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let amplitude = 1000.0 * (i + 1) as f64;
                let xic: Vec<(f64, f64)> = (0..100)
                    .map(|s| {
                        let rt = s as f64 * 0.1;
                        let p1 = amplitude * (-((rt - 3.0).powi(2)) / (2.0 * 0.3 * 0.3)).exp();
                        let p2 =
                            amplitude * 0.5 * (-((rt - 7.0).powi(2)) / (2.0 * 0.3 * 0.3)).exp();
                        (rt, p1 + p2)
                    })
                    .collect();
                (i, xic)
            })
            .collect();

        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);

        assert!(
            peaks.len() >= 2,
            "Should find at least 2 peaks, found {}",
            peaks.len()
        );

        // Peaks should be near RT=3.0 and RT=7.0
        let apex_rts: Vec<f64> = peaks.iter().map(|p| p.apex_rt).collect();
        assert!(
            apex_rts.iter().any(|&rt| (rt - 3.0).abs() < 0.5),
            "Should find peak near 3.0, got {:?}",
            apex_rts
        );
        assert!(
            apex_rts.iter().any(|&rt| (rt - 7.0).abs() < 0.5),
            "Should find peak near 7.0, got {:?}",
            apex_rts
        );

        // First peak should have highest consensus coefficient (strongest)
        assert!(
            peaks[0].apex_intensity >= peaks[1].apex_intensity,
            "Peaks should be sorted by intensity descending"
        );
    }

    #[test]
    fn test_consensus_degenerate_single_xic() {
        // Only 1 XIC → should return empty (need ≥2 for consensus)
        let xics = vec![(0, make_gaussian_xic(5.0, 0.3, 1000.0, 100, 0.1))];
        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);
        assert!(peaks.is_empty(), "Should return empty for single XIC");
    }

    #[test]
    fn test_consensus_degenerate_short_xic() {
        // Only 3 scans → should return empty (need ≥5)
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let xic = vec![(0.0, 100.0), (0.1, 500.0), (0.2, 100.0)];
                (i, xic)
            })
            .collect();
        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);
        assert!(peaks.is_empty(), "Should return empty for short XICs");
    }

    #[test]
    fn test_consensus_all_zero() {
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let xic: Vec<(f64, f64)> = (0..50).map(|s| (s as f64 * 0.1, 0.0)).collect();
                (i, xic)
            })
            .collect();
        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);
        assert!(peaks.is_empty(), "Should return empty for all-zero XICs");
    }

    #[test]
    fn test_consensus_peak_bounds_valid() {
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| (i, make_gaussian_xic(5.0, 0.3, 1000.0, 100, 0.1)))
            .collect();

        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);
        let n_scans = xics[0].1.len();

        for peak in &peaks {
            assert!(
                peak.start_index <= peak.apex_index,
                "start_index {} > apex_index {}",
                peak.start_index,
                peak.apex_index
            );
            assert!(
                peak.apex_index <= peak.end_index,
                "apex_index {} > end_index {}",
                peak.apex_index,
                peak.end_index
            );
            assert!(
                peak.end_index < n_scans,
                "end_index {} >= n_scans {}",
                peak.end_index,
                n_scans
            );
            assert!(peak.start_rt <= peak.apex_rt);
            assert!(peak.apex_rt <= peak.end_rt);
            assert!(peak.area >= 0.0);
        }
    }

    #[test]
    fn test_consensus_noise_robustness() {
        // Peak + deterministic noise pattern
        let xics: Vec<(usize, Vec<(f64, f64)>)> = (0..6)
            .map(|i| {
                let xic: Vec<(f64, f64)> = (0..100)
                    .map(|s| {
                        let rt = s as f64 * 0.1;
                        let signal = 1000.0 * (-((rt - 5.0).powi(2)) / (2.0 * 0.3 * 0.3)).exp();
                        // Deterministic "noise": small sinusoid unique per fragment
                        let noise = 50.0 * ((s as f64 * 0.7 + i as f64 * 1.3).sin());
                        (rt, (signal + noise).max(0.0))
                    })
                    .collect();
                (i, xic)
            })
            .collect();

        let peaks = detect_cwt_consensus_peaks(&xics, 0.0);
        assert!(!peaks.is_empty(), "Should find peak despite noise");
        assert!(
            (peaks[0].apex_rt - 5.0).abs() < 0.5,
            "Peak should be near 5.0 despite noise, got {:.2}",
            peaks[0].apex_rt
        );
    }
}
