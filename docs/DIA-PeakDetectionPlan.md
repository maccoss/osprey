# Technical Specification: Multi-Transition CWT Peak Detection for DIA Proteomics
**Context:** Rust Implementation for Data Independent Acquisition (DIA)

---

## 1. Executive Summary
This document outlines a strategy for detecting chromatographic peaks in noisy Data Independent Acquisition (DIA) data. The core approach utilizes a **Continuous Wavelet Transform (CWT)** to suppress noise and a **Multi-Transition Consensus** mechanism to "borrow information" across co-varying product ions. This allows for robust detection of peak boundaries (start/end times) even when individual transitions are compromised by interference.

## 2. Mathematical Theory
The peak detection relies on the **Mexican Hat Wavelet** (Ricker Wavelet), which serves as a matched filter for Gaussian-like chromatographic peaks.

### 2.1 The Mother Wavelet
The continuous form of the Mexican Hat wavelet is the negative normalized second derivative of a Gaussian:

$$\psi(t) = \frac{2}{\sqrt{3\sigma}\pi^{1/4}} \left( 1 - \left(\frac{t}{\sigma}\right)^2 \right) e^{-\frac{t^2}{2\sigma^2}}$$

Where:
* $t$: Time (or index in the discrete domain).
* $\sigma$: Scale parameter, related to the peak width.

### 2.2 Consensus CWT (The "Borrowing" Step)
Standard peak pickers operate on single signals. To leverage the correlated nature of DIA transitions ($y_1, y_2, ... y_n$), we compute a **Consensus Scalogram**:

$$S_{consensus}(a, b) = \text{Median}_{i=1..n} \left( CWT(y_i, a, b) \right)$$

Using the **median** (or a weighted sum) effectively filters out non-Gaussian noise spikes that occur in only one transition (interference), as the consensus signal will only rise if the majority of transitions exhibit a peak-like shape at the same time $b$ and scale $a$.

---

## 3. Algorithm Lifecycle

### Phase 1: Preprocessing & Transformation
1.  **Input:** A matrix of intensities `[n_transitions, n_timepoints]`.
2.  **Smoothing (Optional):** Apply a Savitzky-Golay filter if data is extremely jagged, though CWT handles noise naturally.
3.  **Wavelet Convolution:** For each transition, convolve the intensity vector with the Mexican Hat kernel at a specific scale $a$ (or a range of scales).

### Phase 2: Peak Candidate Detection
1.  **Consensus Aggregation:** Collapse the $n$ convolution results into a single 1D array (the consensus signal).
2.  **Ridge Detection:** Identify local maxima (peaks) in the consensus signal.
3.  **Thresholding:** Discard candidates with a consensus score below a noise floor (e.g., Signal-to-Noise Ratio < 3).

### Phase 3: Boundary Definition (The "Zero-Crossing" Logic)
1.  **Apex:** The index of the maximum in the consensus signal.
2.  **Boundaries:** Traverse left and right from the apex until the consensus signal drops below zero (or a strict percentage, e.g., 5% of max).
    * *Note:* Using the wavelet signal for boundaries is more robust than raw intensity because the wavelet naturally decays to zero at the inflection points of the peak.

### Phase 4: Quantitation
1.  **Integration:** Switch back to **raw intensity data**.
2.  **Background Subtraction:** Estimate a linear baseline connecting the start and end intensity points.
3.  **Calculation:** Use the Trapezoidal Rule to calculate the area under the curve (AUC) minus the trapezoidal background area.

---

## 4. Rust Implementation Guide

The following Rust code illustrates the core structures and logic. It utilizes `ndarray` for efficient matrix operations.

### 4.1 Data Structures

```rust
use ndarray::{Array1, Array2, Axis};

/// Container for the raw chromatographic data of a precursor
pub struct PrecursorData {
    pub retention_times: Array1<f64>,
    pub intensities: Array2<f64>, // Shape: [n_transitions, n_timepoints]
}

/// Output metrics for a detected peak
pub struct PeakMetrics {
    pub apex_rt: f64,
    pub centroid_rt: f64,
    pub left_rt: f64,
    pub right_rt: f64,
    pub total_area: f64,
    pub transition_areas: Vec<f64>,
    pub transition_backgrounds: Vec<f64>,
}
```

### 4.2 The Mexican Hat Kernel Generator

```rust
use std::f64::consts::PI;
/// Generates a discrete Mexican Hat wavelet kernel (normalized).
pub fn generate_mexican_hat(scale: f64, num_points: usize) -> Array1<f64> {
    let mut kernel = Array1::<f64>::zeros(num_points);
    let half_width = (num_points as f64 - 1.0) / 2.0;
    
    // Normalization factor for energy preservation
    let norm = 2.0 / ((3.0 * scale).sqrt() * PI.powf(0.25));

    for i in 0..num_points {
        let t = (i as f64 - half_width) / scale;
        // The Mexican Hat function
        let val = norm * (1.0 - t.powi(2)) * (-0.5 * t.powi(2)).exp();
        kernel[i] = val;
    }

    // Zero-mean correction (crucial for wavelet behavior)
    let mean = kernel.sum() / num_points as f64;
    kernel.mapv_inplace(|x| x - mean);
    
    kernel
}
```

### 4.3 Consensus Logic & Peak Picking
```rust
/// Main function to detect peaks and return boundaries
pub fn detect_peak_boundaries(
    data: &PrecursorData,
    scale: f64,
) -> Option<(usize, usize, usize)> { // (Left_Idx, Apex_Idx, Right_Idx)
    
    let n_points = data.intensities.ncols();
    let n_transitions = data.intensities.nrows();
    
    // 1. Generate Kernel (Window size usually 10x scale)
    let kernel_size = (scale * 10.0) as usize;
    let kernel = generate_mexican_hat(scale, kernel_size);
    
    // 2. Compute Consensus CWT
    let mut consensus = Array1::<f64>::zeros(n_points);
    
    // Simplified: Summing CWTs (Real implementation: use Median for robustness)
    for i in 0..n_transitions {
        let signal = data.intensities.row(i);
        // Note: You will need a convolution function (e.g. from `scipy` equivalent or FFT)
        let cwt = perform_convolution(&signal, &kernel); 
        consensus = consensus + cwt;
    }

    // 3. Find Apex (Global Max for this example)
    let (apex_idx, &max_val) = consensus.iter().enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())?;

    if max_val <= 0.0 { return None; } // No peak found

    // 4. Walk Left/Right to Zero Crossings
    let mut left_idx = apex_idx;
    while left_idx > 0 && consensus[left_idx] > 0.0 {
        left_idx -= 1;
    }

    let mut right_idx = apex_idx;
    while right_idx < n_points - 1 && consensus[right_idx] > 0.0 {
        right_idx += 1;
    }

    Some((left_idx, apex_idx, right_idx))
}
```


### 4.4 Quantitation (Trapezoidal + Background Subtraction)
```rust
pub fn quantify_peak(
    data: &PrecursorData,
    bounds: (usize, usize, usize),
) -> PeakMetrics {
    let (left, apex, right) = bounds;
    let rts = &data.retention_times;
    
    let mut areas = Vec::new();
    let mut backgrounds = Vec::new();
    
    // Iterate over each transition to integrate
    for i in 0..data.intensities.nrows() {
        let signal = data.intensities.row(i);
        
        // Define Linear Background
        let y_start = signal[left];
        let y_end = signal[right];
        let t_start = rts[left];
        let t_end = rts[right];
        let width = t_end - t_start;
        
        // Area of trapezoid under the baseline
        let bg_area = 0.5 * (y_start + y_end) * width;
        
        // Trapezoidal Integration of Signal
        let mut raw_area = 0.0;
        for k in left..right {
            let dt = rts[k+1] - rts[k];
            let height = 0.5 * (signal[k] + signal[k+1]);
            raw_area += height * dt;
        }
        
        areas.push(raw_area - bg_area);
        backgrounds.push(bg_area);
    }
    
    // Calculate Centroid (Intensity Weighted RT)
    // Implementation omitted for brevity (Sum(RT * I) / Sum(I))

    PeakMetrics {
        apex_rt: rts[apex],
        centroid_rt: 0.0, // Calculate as above
        left_rt: rts[left],
        right_rt: rts[right],
        total_area: areas.iter().sum(),
        transition_areas: areas,
        transition_backgrounds: backgrounds,
    }
}
```

5. Advanced Considerations
5.1 Handling Interference
In DIA, it is common for 1 out of 6 transitions to be contaminated. The Consensus step in 4.3 uses a simple sum. Changing this to a Tukey median polished chromatograms to improve robustness:

Calculate the correlation of each transition to the "median" shape.


5.2 Scale Space Search
A robust detector should not rely on a single scale. Wrap the logic in 4.3 in a loop over scales = [2.0, 4.0, 8.0, 16.0] (seconds). The "true" peak is the one that appears consistently across the most scales (Ridge Tracking).
