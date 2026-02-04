//! Pipeline coordinator for streaming mzML processing
//!
//! This module orchestrates the pipelined architecture:
//! 1. Async mzML parsing (I/O bound)
//! 2. Parallel preprocessing (CPU bound)
//! 3. Window-based accumulation
//! 4. Batch scoring

use std::path::Path;
use std::sync::Arc;
use std::thread;

use crossbeam::channel::{bounded, Receiver, Sender};
use rayon::prelude::*;

use osprey_core::{LibraryEntry, Spectrum};

use super::accumulator::{WindowAccumulator, WindowKey};
use super::preprocessor::{PreprocessedSpectrum, PreprocessingConfig, PreprocessingWorker};
use crate::batch::{CalibrationMatch, PreprocessedLibrary};
use crate::SpectralScorer;

/// Configuration for the streaming pipeline
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Number of preprocessing threads (default: num_cpus)
    pub num_preprocessing_threads: usize,
    /// Channel buffer size for backpressure (default: 2 * num_cpus)
    pub channel_buffer_size: usize,
    /// RT tolerance for matching (minutes)
    pub rt_tolerance: f64,
    /// Preprocessing configuration
    pub preprocessing: PreprocessingConfig,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        let num_cpus = rayon::current_num_threads();
        Self {
            num_preprocessing_threads: num_cpus,
            channel_buffer_size: num_cpus * 2,
            rt_tolerance: 5.0,
            preprocessing: PreprocessingConfig::default(),
        }
    }
}

/// Result of pipeline processing for a single window
#[derive(Debug)]
pub struct WindowResult {
    /// Window key
    pub window: WindowKey,
    /// Number of spectra in window
    pub num_spectra: usize,
    /// Number of library entries scored
    pub num_entries: usize,
    /// Best matches found
    pub matches: Vec<CalibrationMatch>,
}

/// Pipeline coordinator that orchestrates streaming processing
pub struct PipelineCoordinator {
    config: PipelineConfig,
}

impl Default for PipelineCoordinator {
    fn default() -> Self {
        Self::new()
    }
}

impl PipelineCoordinator {
    /// Create a new pipeline coordinator with default config
    pub fn new() -> Self {
        Self {
            config: PipelineConfig::default(),
        }
    }

    /// Create with custom config
    pub fn with_config(config: PipelineConfig) -> Self {
        Self { config }
    }

    /// Run the streaming pipeline on a collection of spectra
    ///
    /// This is a simpler interface that takes pre-loaded spectra and
    /// processes them through the pipeline. Useful for testing and
    /// when spectra are already in memory.
    pub fn run_on_spectra(
        &self,
        spectra: &[Spectrum],
        library: &[LibraryEntry],
    ) -> Vec<CalibrationMatch> {
        log::info!(
            "Starting streaming pipeline: {} spectra, {} library entries",
            spectra.len(),
            library.len()
        );

        // Create channels
        let (raw_tx, raw_rx): (Sender<Spectrum>, Receiver<Spectrum>) =
            bounded(self.config.channel_buffer_size);
        let (preprocessed_tx, preprocessed_rx): (
            Sender<PreprocessedSpectrum>,
            Receiver<PreprocessedSpectrum>,
        ) = bounded(self.config.channel_buffer_size);

        // Clone config for threads
        let preprocessing_config = self.config.preprocessing.clone();
        let rt_tolerance = self.config.rt_tolerance;

        // Preprocess library upfront (this is fast compared to spectra)
        let preprocessed_library = PreprocessedLibrary::from_entries(library);
        let library_arc = Arc::new(library.to_vec());
        let library_arc_clone = library_arc.clone();

        // Spawn preprocessing threads
        let num_workers = self.config.num_preprocessing_threads;
        let preprocessing_handles: Vec<_> = (0..num_workers)
            .map(|_| {
                let rx = raw_rx.clone();
                let tx = preprocessed_tx.clone();
                let config = preprocessing_config.clone();

                thread::spawn(move || {
                    let worker = PreprocessingWorker::with_config(config);
                    while let Ok(spectrum) = rx.recv() {
                        let preprocessed = worker.process(&spectrum);
                        if tx.send(preprocessed).is_err() {
                            break;
                        }
                    }
                })
            })
            .collect();

        // Drop extra senders so receivers know when to stop
        drop(raw_rx);
        drop(preprocessed_tx);

        // Send all spectra to preprocessing
        for spectrum in spectra {
            if raw_tx.send(spectrum.clone()).is_err() {
                break;
            }
        }
        drop(raw_tx);

        // Accumulate preprocessed spectra by window
        let mut accumulator = WindowAccumulator::new();
        while let Ok(preprocessed) = preprocessed_rx.recv() {
            accumulator.add_spectrum(preprocessed);
        }

        // Wait for preprocessing to complete
        for handle in preprocessing_handles {
            let _ = handle.join();
        }

        log::debug!(
            "Accumulator stats: {:?}",
            accumulator.stats()
        );

        // Score all windows
        let windows = accumulator.drain_all();
        let scorer = SpectralScorer::new();

        let all_matches: Vec<CalibrationMatch> = windows
            .into_par_iter()
            .flat_map(|(window_key, window_spectra)| {
                self.score_window(
                    window_key,
                    &window_spectra,
                    &library_arc_clone,
                    &preprocessed_library,
                    &scorer,
                    rt_tolerance,
                )
            })
            .collect();

        // Deduplicate (keep best score per entry)
        let mut best_matches: std::collections::HashMap<u32, CalibrationMatch> =
            std::collections::HashMap::new();
        for m in all_matches {
            best_matches
                .entry(m.entry_id)
                .and_modify(|existing| {
                    if m.score > existing.score {
                        *existing = m.clone();
                    }
                })
                .or_insert(m);
        }

        let mut results: Vec<CalibrationMatch> = best_matches.into_values().collect();
        results.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        log::info!(
            "Streaming pipeline complete: {} unique matches",
            results.len()
        );

        results
    }

    /// Score a single window's spectra against the library
    fn score_window(
        &self,
        window_key: WindowKey,
        window_spectra: &[PreprocessedSpectrum],
        library: &[LibraryEntry],
        _preprocessed_library: &PreprocessedLibrary,
        scorer: &SpectralScorer,
        rt_tolerance: f64,
    ) -> Vec<CalibrationMatch> {
        // Filter library entries to this window
        let lower = window_key.lower();
        let upper = window_key.upper();

        let window_entries: Vec<&LibraryEntry> = library
            .iter()
            .filter(|e| e.precursor_mz >= lower && e.precursor_mz <= upper)
            .collect();

        if window_entries.is_empty() || window_spectra.is_empty() {
            return Vec::new();
        }

        // Precompute XCorr vectors for entries
        let entry_xcorr: std::collections::HashMap<u32, Vec<f64>> = window_entries
            .iter()
            .map(|e| (e.id, scorer.preprocess_library_for_xcorr(e)))
            .collect();

        // Score each entry against window spectra
        let mut matches = Vec::new();

        for entry in &window_entries {
            let mut best_score = 0.0f64;
            let mut best_rt = 0.0f64;
            let mut best_xcorr = 0.0f64;
            let mut best_scan = 0u32;

            for preprocessed in window_spectra {
                // RT tolerance check
                let rt_diff = (preprocessed.retention_time - entry.retention_time).abs();
                if rt_diff > rt_tolerance {
                    continue;
                }

                // LibCosine score from preprocessed vectors
                // This is a simplified version - in production we'd use
                // the full BLAS-accelerated batch scoring
                let libcosine_score = self.compute_libcosine_from_preprocessed(
                    &preprocessed.libcosine_vector,
                    entry,
                    &self.config.preprocessing,
                );

                // XCorr score
                let xcorr_score = if let Some(entry_vec) = entry_xcorr.get(&entry.id) {
                    SpectralScorer::xcorr_from_preprocessed(&preprocessed.xcorr_vector, entry_vec)
                } else {
                    0.0
                };

                if libcosine_score > best_score {
                    best_score = libcosine_score;
                    best_rt = preprocessed.retention_time;
                    best_xcorr = xcorr_score;
                    best_scan = preprocessed.scan_number;
                }
            }

            if best_score > 0.0 {
                matches.push(CalibrationMatch {
                    entry_id: entry.id,
                    is_decoy: entry.is_decoy,
                    library_rt: entry.retention_time,
                    measured_rt: best_rt,
                    score: best_score,
                    ms1_ppm_error: None,
                    library_precursor_mz: entry.precursor_mz,
                    observed_precursor_mz: Some(window_key.center()),
                    ms2_mass_errors: Vec::new(),
                    avg_ms2_error: None,
                    n_matched_fragments: 0,
                    n_library_fragments: entry.fragments.len(),
                    xcorr_score: best_xcorr,
                    isotope_cosine_score: None,
                    sequence: entry.sequence.clone(),
                    charge: entry.charge,
                    scan_number: best_scan,
                });
            }
        }

        matches
    }

    /// Compute LibCosine score from preprocessed spectrum vector and library entry
    fn compute_libcosine_from_preprocessed(
        &self,
        spectrum_vec: &ndarray::Array1<f64>,
        entry: &LibraryEntry,
        config: &PreprocessingConfig,
    ) -> f64 {
        // Build library vector
        let mut lib_vec = ndarray::Array1::zeros(config.num_bins);
        for frag in &entry.fragments {
            let bin = ((frag.mz - config.min_mz) / config.bin_width) as usize;
            if bin < config.num_bins {
                let smz = (frag.relative_intensity as f64).sqrt() * frag.mz.powi(2);
                lib_vec[bin] += smz;
            }
        }

        // L2 normalize
        let lib_norm: f64 = lib_vec.iter().map(|x| x * x).sum::<f64>().sqrt();
        if lib_norm < 1e-10 {
            return 0.0;
        }
        lib_vec /= lib_norm;

        // Dot product (spectrum_vec is already normalized)
        let dot: f64 = spectrum_vec
            .iter()
            .zip(lib_vec.iter())
            .map(|(s, l)| s * l)
            .sum();

        dot.max(0.0).min(1.0)
    }
}

/// Run streaming pipeline from an mzML file
///
/// This is the main entry point for streaming processing.
/// Requires the `streaming` feature.
#[cfg(feature = "streaming")]
pub async fn run_streaming_pipeline(
    mzml_path: impl AsRef<Path>,
    library: &[LibraryEntry],
    config: PipelineConfig,
) -> Result<Vec<CalibrationMatch>, StreamingPipelineError> {
    use osprey_io::mzml::streaming::StreamingMzmlReader;

    let path = mzml_path.as_ref();
    log::info!("Starting streaming pipeline for: {}", path.display());

    // Create channels for spectrum streaming
    let (spectrum_tx, spectrum_rx) = bounded(config.channel_buffer_size);
    let (preprocessed_tx, preprocessed_rx) = bounded(config.channel_buffer_size);

    // Clone config for tasks
    let preprocessing_config = config.preprocessing.clone();
    let rt_tolerance = config.rt_tolerance;

    // Preprocess library upfront
    let preprocessed_library = PreprocessedLibrary::from_entries(library);
    let library_arc = Arc::new(library.to_vec());

    // Spawn preprocessing threads
    let num_workers = config.num_preprocessing_threads;
    let preprocessing_handles: Vec<_> = (0..num_workers)
        .map(|_| {
            let rx = spectrum_rx.clone();
            let tx = preprocessed_tx.clone();
            let cfg = preprocessing_config.clone();

            thread::spawn(move || {
                let worker = PreprocessingWorker::with_config(cfg);
                while let Ok(result) = rx.recv() {
                    match result {
                        Ok(spectrum) => {
                            let preprocessed = worker.process(&spectrum);
                            if tx.send(preprocessed).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            log::warn!("Error parsing spectrum: {:?}", e);
                        }
                    }
                }
            })
        })
        .collect();

    // Drop extra receivers/senders
    drop(spectrum_rx);
    drop(preprocessed_tx);

    // Start async mzML parsing
    let reader = StreamingMzmlReader::ms2_only();
    let stats = reader.parse_file_to_channel(path, spectrum_tx).await
        .map_err(|e| StreamingPipelineError::ParseError(e.to_string()))?;

    log::info!(
        "Parsed {} spectra (MS2: {}, skipped: {})",
        stats.total_spectra,
        stats.ms2_spectra,
        stats.skipped_spectra
    );

    // Accumulate preprocessed spectra
    let mut accumulator = WindowAccumulator::new();
    while let Ok(preprocessed) = preprocessed_rx.recv() {
        accumulator.add_spectrum(preprocessed);
    }

    // Wait for preprocessing threads
    for handle in preprocessing_handles {
        let _ = handle.join();
    }

    log::debug!("Accumulator stats: {:?}", accumulator.stats());

    // Score all windows
    let windows = accumulator.drain_all();
    let scorer = SpectralScorer::new();
    let coordinator = PipelineCoordinator::with_config(config);

    let all_matches: Vec<CalibrationMatch> = windows
        .into_par_iter()
        .flat_map(|(window_key, window_spectra)| {
            coordinator.score_window(
                window_key,
                &window_spectra,
                &library_arc,
                &preprocessed_library,
                &scorer,
                rt_tolerance,
            )
        })
        .collect();

    // Deduplicate
    let mut best_matches: std::collections::HashMap<u32, CalibrationMatch> =
        std::collections::HashMap::new();
    for m in all_matches {
        best_matches
            .entry(m.entry_id)
            .and_modify(|existing| {
                if m.score > existing.score {
                    *existing = m.clone();
                }
            })
            .or_insert(m);
    }

    let mut results: Vec<CalibrationMatch> = best_matches.into_values().collect();
    results.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    log::info!("Streaming pipeline complete: {} matches", results.len());

    Ok(results)
}

/// Error type for streaming pipeline
#[cfg(feature = "streaming")]
#[derive(thiserror::Error, Debug)]
pub enum StreamingPipelineError {
    #[error("mzML parsing error: {0}")]
    ParseError(String),

    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
}

#[cfg(test)]
mod tests {
    use super::*;
    use osprey_core::{IsolationWindow, LibraryFragment};

    fn make_test_spectrum(scan: u32, rt: f64, precursor_mz: f64, peaks: Vec<(f64, f32)>) -> Spectrum {
        Spectrum {
            scan_number: scan,
            retention_time: rt,
            precursor_mz,
            isolation_window: IsolationWindow::symmetric(precursor_mz, 12.5),
            mzs: peaks.iter().map(|(mz, _)| *mz).collect(),
            intensities: peaks.iter().map(|(_, int)| *int).collect(),
        }
    }

    fn make_test_entry(id: u32, precursor_mz: f64, rt: f64, fragments: Vec<(f64, f32)>) -> LibraryEntry {
        let mut entry = LibraryEntry::new(id, "TEST".into(), "TEST".into(), 2, precursor_mz, rt);
        entry.fragments = fragments
            .into_iter()
            .map(|(mz, int)| LibraryFragment {
                mz,
                relative_intensity: int,
                annotation: osprey_core::FragmentAnnotation::default(),
            })
            .collect();
        entry
    }

    #[test]
    fn test_pipeline_basic() {
        let spectra = vec![
            make_test_spectrum(1, 10.0, 500.0, vec![(300.0, 100.0), (400.0, 50.0)]),
            make_test_spectrum(2, 10.5, 500.0, vec![(300.0, 80.0), (400.0, 60.0)]),
        ];

        let library = vec![make_test_entry(
            1,
            500.0,
            10.0,
            vec![(300.0, 100.0), (400.0, 50.0)],
        )];

        let coordinator = PipelineCoordinator::new();
        let matches = coordinator.run_on_spectra(&spectra, &library);

        assert!(!matches.is_empty());
        assert_eq!(matches[0].entry_id, 1);
        assert!(matches[0].score > 0.0);
    }
}
