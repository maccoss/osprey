//! Window-aware spectrum accumulator for pipelined scoring
//!
//! This module collects preprocessed spectra by isolation window,
//! detects when windows are complete, and triggers batch scoring.

use std::collections::HashMap;

use super::preprocessor::PreprocessedSpectrum;

/// Key for identifying isolation windows
///
/// Windows are identified by their discretized bounds to handle
/// floating point comparison issues.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct WindowKey {
    /// Lower bound * 1000 (3 decimal places)
    pub lower_millimz: i64,
    /// Upper bound * 1000 (3 decimal places)
    pub upper_millimz: i64,
}

impl WindowKey {
    /// Create a window key from isolation window bounds
    pub fn from_bounds(lower: f64, upper: f64) -> Self {
        Self {
            lower_millimz: (lower * 1000.0).round() as i64,
            upper_millimz: (upper * 1000.0).round() as i64,
        }
    }

    /// Get the center m/z of this window
    pub fn center(&self) -> f64 {
        (self.lower_millimz + self.upper_millimz) as f64 / 2000.0
    }

    /// Get the lower bound
    pub fn lower(&self) -> f64 {
        self.lower_millimz as f64 / 1000.0
    }

    /// Get the upper bound
    pub fn upper(&self) -> f64 {
        self.upper_millimz as f64 / 1000.0
    }
}

/// Buffer for spectra within a single isolation window
#[derive(Debug)]
struct WindowBuffer {
    /// Collected spectra for this window
    spectra: Vec<PreprocessedSpectrum>,
    /// Last scan number seen in this window
    last_scan: u32,
    /// Number of times this window has cycled (for future cycle detection)
    _cycle_count: u32,
    /// Whether this window is marked complete
    is_complete: bool,
}

impl WindowBuffer {
    fn new() -> Self {
        Self {
            spectra: Vec::new(),
            last_scan: 0,
            _cycle_count: 0,
            is_complete: false,
        }
    }
}

/// Accumulator that groups preprocessed spectra by isolation window
///
/// Detects when windows are complete based on:
/// 1. DIA cycle detection (seeing window again after full cycle)
/// 2. Explicit completion markers
/// 3. End of file
#[derive(Debug)]
pub struct WindowAccumulator {
    /// Buffers for each isolation window
    windows: HashMap<WindowKey, WindowBuffer>,
    /// Order in which windows were first seen (for cycle detection)
    window_order: Vec<WindowKey>,
    /// Current position in window cycle (for future cycle detection)
    _current_cycle_index: usize,
    /// Number of complete DIA cycles seen (for future cycle detection)
    _total_cycles: u32,
    /// Minimum spectra required to consider a window valid
    min_spectra_per_window: usize,
}

impl Default for WindowAccumulator {
    fn default() -> Self {
        Self::new()
    }
}

impl WindowAccumulator {
    /// Create a new window accumulator
    pub fn new() -> Self {
        Self {
            windows: HashMap::new(),
            window_order: Vec::new(),
            _current_cycle_index: 0,
            _total_cycles: 0,
            min_spectra_per_window: 1,
        }
    }

    /// Set minimum spectra required per window
    pub fn with_min_spectra(mut self, min: usize) -> Self {
        self.min_spectra_per_window = min;
        self
    }

    /// Add a preprocessed spectrum to the accumulator
    ///
    /// Returns a list of windows that became complete after this addition.
    pub fn add_spectrum(&mut self, spectrum: PreprocessedSpectrum) -> Vec<WindowKey> {
        let key = WindowKey::from_bounds(
            spectrum.isolation_window.lower_bound(),
            spectrum.isolation_window.upper_bound(),
        );

        // Track window order for cycle detection
        if !self.windows.contains_key(&key) {
            self.window_order.push(key);
        }

        // Add spectrum to window buffer
        let buffer = self.windows.entry(key).or_insert_with(WindowBuffer::new);
        buffer.last_scan = spectrum.scan_number;
        buffer.spectra.push(spectrum);

        // Check for cycle completion
        self.detect_completed_windows()
    }

    /// Detect windows that have completed based on cycle detection
    fn detect_completed_windows(&mut self) -> Vec<WindowKey> {
        // Simple heuristic: if we've seen all windows at least once,
        // consider earlier windows complete when we see them again
        // This is based on DIA cycling through windows: 1,2,3,...,N,1,2,3,...

        let completed = Vec::new();

        // If we have enough windows to form a cycle
        if self.window_order.len() >= 2 {
            // TODO: Implement cycle detection logic
            // Check if we've cycled back to an earlier window
            // For simplicity, we mark windows complete when all windows
            // have at least some spectra and we've completed a full cycle
        }

        completed
    }

    /// Mark all windows as complete (call at end of file)
    pub fn mark_all_complete(&mut self) {
        for buffer in self.windows.values_mut() {
            buffer.is_complete = true;
        }
    }

    /// Get all complete windows and their spectra
    ///
    /// Returns windows that are marked complete and have enough spectra.
    /// Removes the returned windows from the accumulator.
    pub fn drain_complete_windows(&mut self) -> Vec<(WindowKey, Vec<PreprocessedSpectrum>)> {
        let mut result = Vec::new();

        let complete_keys: Vec<WindowKey> = self
            .windows
            .iter()
            .filter(|(_, buf)| buf.is_complete && buf.spectra.len() >= self.min_spectra_per_window)
            .map(|(k, _)| *k)
            .collect();

        for key in complete_keys {
            if let Some(buffer) = self.windows.remove(&key) {
                result.push((key, buffer.spectra));
            }
        }

        result
    }

    /// Get all windows (complete or not) and drain the accumulator
    ///
    /// Call this at the end of processing to get remaining windows.
    pub fn drain_all(&mut self) -> Vec<(WindowKey, Vec<PreprocessedSpectrum>)> {
        self.mark_all_complete();
        self.drain_complete_windows()
    }

    /// Get the number of windows currently being accumulated
    pub fn num_windows(&self) -> usize {
        self.windows.len()
    }

    /// Get the total number of spectra accumulated
    pub fn num_spectra(&self) -> usize {
        self.windows.values().map(|b| b.spectra.len()).sum()
    }

    /// Get statistics about accumulated windows
    pub fn stats(&self) -> AccumulatorStats {
        let window_sizes: Vec<usize> = self.windows.values().map(|b| b.spectra.len()).collect();
        let total_spectra: usize = window_sizes.iter().sum();
        let avg_per_window = if window_sizes.is_empty() {
            0.0
        } else {
            total_spectra as f64 / window_sizes.len() as f64
        };

        AccumulatorStats {
            num_windows: self.windows.len(),
            total_spectra,
            avg_spectra_per_window: avg_per_window,
            min_spectra_per_window: window_sizes.iter().cloned().min().unwrap_or(0),
            max_spectra_per_window: window_sizes.iter().cloned().max().unwrap_or(0),
        }
    }
}

/// Statistics about accumulated windows
#[derive(Debug, Clone)]
pub struct AccumulatorStats {
    pub num_windows: usize,
    pub total_spectra: usize,
    pub avg_spectra_per_window: f64,
    pub min_spectra_per_window: usize,
    pub max_spectra_per_window: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array1;
    use osprey_core::IsolationWindow;

    fn make_preprocessed(scan: u32, rt: f64, window: (f64, f64)) -> PreprocessedSpectrum {
        let center = (window.0 + window.1) / 2.0;
        let half_width = (window.1 - window.0) / 2.0;
        PreprocessedSpectrum {
            scan_number: scan,
            retention_time: rt,
            isolation_window: IsolationWindow::new(center, half_width, half_width),
            libcosine_vector: Array1::zeros(100),
            xcorr_vector: vec![0.0; 100],
        }
    }

    #[test]
    fn test_accumulator_grouping() {
        let mut acc = WindowAccumulator::new();

        // Add spectra from different windows
        acc.add_spectrum(make_preprocessed(1, 1.0, (400.0, 425.0)));
        acc.add_spectrum(make_preprocessed(2, 1.1, (425.0, 450.0)));
        acc.add_spectrum(make_preprocessed(3, 1.2, (400.0, 425.0)));
        acc.add_spectrum(make_preprocessed(4, 1.3, (425.0, 450.0)));

        assert_eq!(acc.num_windows(), 2);
        assert_eq!(acc.num_spectra(), 4);

        let stats = acc.stats();
        assert_eq!(stats.avg_spectra_per_window, 2.0);
    }

    #[test]
    fn test_drain_all() {
        let mut acc = WindowAccumulator::new();

        acc.add_spectrum(make_preprocessed(1, 1.0, (400.0, 425.0)));
        acc.add_spectrum(make_preprocessed(2, 1.1, (425.0, 450.0)));

        let drained = acc.drain_all();
        assert_eq!(drained.len(), 2);
        assert_eq!(acc.num_windows(), 0);
    }

    #[test]
    fn test_window_key() {
        let key1 = WindowKey::from_bounds(400.0, 425.0);
        let key2 = WindowKey::from_bounds(400.0, 425.0);
        let key3 = WindowKey::from_bounds(400.0004, 425.0004); // Within rounding tolerance

        assert_eq!(key1, key2);
        assert_eq!(key1, key3); // Close enough to round to same value
        assert!((key1.center() - 412.5).abs() < 0.01);
    }
}
