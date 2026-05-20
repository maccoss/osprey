//! mzML file parser using mzdata crate
//!
//! Converts mzdata spectrum representations to Osprey's internal types.
//!
//! ## MS1 Spectra for Mass Calibration
//!
//! Extracts the M+0 isotope peak from MS1 spectra for accurate mass calibration.
//! This module provides:
//! - `load_all_spectra()` - Load both MS1 and MS2 spectra in a single pass (most efficient)
//! - `load_ms1_spectra()` - Load only MS1 spectra for isotope extraction
//! - `MS1Index` - Index MS1 spectra by retention time for efficient lookup

use mzdata::io::mzml::MzMLReader;
use mzdata::prelude::*;
use mzdata::spectrum::RefPeakDataLevel;
use osprey_core::{IsolationWindow, MS1Spectrum, OspreyError, Result, Spectrum, SpectrumSource};
use quick_xml::events::Event;
use quick_xml::reader::Reader;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

/// Raw f64 isolation window cvParam values lifted directly from the
/// mzML XML, indexed by spectrum index (matching mzdata's
/// `MultiLayerSpectrum.description().index`). Used to override
/// mzdata 0.63's f32-quantized `IsolationWindow.lower_bound` /
/// `upper_bound` values; see [`read_isolation_cvparams_f64`] for the
/// rationale.
#[derive(Default, Clone, Debug)]
struct IsolationCvParams {
    /// MS:1000827 isolation window target m/z (f64 from XML text).
    target_mz: Option<f64>,
    /// MS:1000828 isolation window lower offset (positive magnitude).
    lower_offset: Option<f64>,
    /// MS:1000829 isolation window upper offset (positive magnitude).
    upper_offset: Option<f64>,
}

/// One-pass streaming scan of an mzML file to extract every
/// `<isolationWindow>` block's cvParam values at f64 precision,
/// directly from the XML text. mzdata 0.63's reader pipes these
/// values through `param.to_f32()` (see mzdata-0.63
/// `src/io/mzml/reader.rs:281`), which quantizes through f32 and
/// produces ~3e-5 m/z drift on isolation-window edges that land
/// between two f32-representable values. The C# port
/// (OspreySharp) parses the same cvParams at f64, so the f32
/// quantization on the Rust side is the sole source of a 1,732-row
/// cross-impl `iso_upper` mismatch on Stellar 3-file. This
/// function reads the same XML text mzdata sees but preserves full
/// f64 precision.
///
/// The override is wired in via [`make_isolation_window`] at the two
/// MS2 parsing sites in this module. When upstream mzdata moves to
/// f64 storage for `IsolationWindow` bounds, this pre-pass + override
/// can be deleted in one commit. See workspace `Cargo.toml`
/// `quick-xml = "0.30"` line for the lifetime expectation.
fn read_isolation_cvparams_f64(path: &Path) -> Result<HashMap<u32, IsolationCvParams>> {
    let mut reader = Reader::from_file(path).map_err(|e| {
        OspreyError::MzmlParseError(format!("quick-xml open '{}': {}", path.display(), e))
    })?;
    let mut buf = Vec::new();
    let mut result: HashMap<u32, IsolationCvParams> = HashMap::new();
    let mut current_index: Option<u32> = None;
    let mut depth_in_isolation_window: i32 = 0;
    let mut current_params = IsolationCvParams::default();

    loop {
        let event = reader.read_event_into(&mut buf).map_err(|e| {
            OspreyError::MzmlParseError(format!(
                "quick-xml read at pos {}: {}",
                reader.buffer_position(),
                e
            ))
        })?;
        match event {
            Event::Start(ref e) if e.name().as_ref() == b"spectrum" => {
                current_index = None;
                for attr in e.attributes().flatten() {
                    if attr.key.as_ref() == b"index" {
                        if let Ok(s) = std::str::from_utf8(attr.value.as_ref()) {
                            current_index = s.parse().ok();
                        }
                    }
                }
                current_params = IsolationCvParams::default();
            }
            Event::Start(ref e) if e.name().as_ref() == b"isolationWindow" => {
                depth_in_isolation_window += 1;
            }
            Event::End(ref e) if e.name().as_ref() == b"isolationWindow" => {
                depth_in_isolation_window -= 1;
                if depth_in_isolation_window == 0 {
                    if let Some(idx) = current_index {
                        result.insert(idx, current_params.clone());
                    }
                }
            }
            // cvParam elements inside isolationWindow are typically
            // self-closing (`Event::Empty`); handle the start variant
            // too for defensive parsing.
            Event::Empty(ref e) | Event::Start(ref e)
                if depth_in_isolation_window > 0 && e.name().as_ref() == b"cvParam" =>
            {
                let mut accession: Option<&[u8]> = None;
                let mut value: Option<&[u8]> = None;
                let attrs: Vec<_> = e.attributes().flatten().collect();
                for attr in &attrs {
                    match attr.key.as_ref() {
                        b"accession" => accession = Some(attr.value.as_ref()),
                        b"value" => value = Some(attr.value.as_ref()),
                        _ => {}
                    }
                }
                if let (Some(acc), Some(v)) = (accession, value) {
                    if let (Ok(acc_str), Ok(v_str)) =
                        (std::str::from_utf8(acc), std::str::from_utf8(v))
                    {
                        if let Ok(f) = v_str.parse::<f64>() {
                            match acc_str {
                                "MS:1000827" => current_params.target_mz = Some(f),
                                "MS:1000828" => current_params.lower_offset = Some(f),
                                "MS:1000829" => current_params.upper_offset = Some(f),
                                _ => {}
                            }
                        }
                    }
                }
            }
            Event::Eof => break,
            _ => {}
        }
        buf.clear();
    }
    Ok(result)
}

/// Build an [`IsolationWindow`] for a spectrum, preferring f64 cvParam
/// values from the pre-parsed map if available. Falls back to
/// mzdata's f32-quantized bounds when (a) no map is provided, (b) the
/// map has no entry for this scan index, or (c) the entry is missing
/// the requested cvParams (e.g., older mzML converters that omit
/// MS:1000827 / 828 / 829).
fn make_isolation_window(
    precursor: &mzdata::spectrum::Precursor,
    scan_number: u32,
    iso_cv_params: Option<&HashMap<u32, IsolationCvParams>>,
) -> Option<IsolationWindow> {
    let ion = precursor.ion()?;
    let center_fallback = ion.mz;

    if let Some(map) = iso_cv_params {
        if let Some(cv) = map.get(&scan_number) {
            if cv.target_mz.is_some() || cv.lower_offset.is_some() || cv.upper_offset.is_some() {
                let center = cv.target_mz.unwrap_or(center_fallback);
                let lower_offset = cv.lower_offset.unwrap_or(12.5);
                let upper_offset = cv.upper_offset.unwrap_or(12.5);
                return Some(IsolationWindow::new(center, lower_offset, upper_offset));
            }
        }
    }

    // Fallback: mzdata-quantized values (f32 → f64 cast).
    let isolation = &precursor.isolation_window;
    let center = center_fallback;
    let lower_offset = if isolation.lower_bound > 0.0 {
        center - isolation.lower_bound as f64
    } else {
        12.5
    };
    let upper_offset = if isolation.upper_bound > 0.0 {
        isolation.upper_bound as f64 - center
    } else {
        12.5
    };
    Some(IsolationWindow::new(center, lower_offset, upper_offset))
}

/// Some mzML producers emit peaks that are not strictly ascending in m/z
/// (observed in a HeLa Astral 3 mz DIA file: ~1 row in 1.7M had a single
/// inverted pair of consecutive centroids). Downstream fragment matching
/// binary-searches the spectrum; the Rust `partition_point` and the C# port's
/// `BinarySearchLowerBound` use procedurally different step patterns, so an
/// unsorted region produces UB-style divergence between the two impls. This
/// helper sorts when an inversion is detected so downstream consumers see a
/// well-defined ordering. The leading O(n) sortedness check is the
/// common-case fast path; the actual sort only runs on inversions.
///
/// Returns `(mzs, intensities, did_sort)`. `did_sort = true` when an
/// inversion was found and the arrays were reordered; the bulk loaders
/// aggregate this into a single per-file end-of-load info log so the
/// default-verbosity log stays uncluttered while `--verbose` still gets
/// the per-scan trace via the `log::debug!` below.
fn ensure_sorted(
    mzs: Vec<f64>,
    intensities: Vec<f32>,
    scan_number: u32,
) -> (Vec<f64>, Vec<f32>, bool) {
    if mzs.len() < 2 || mzs.windows(2).all(|w| w[0] <= w[1]) {
        return (mzs, intensities, false);
    }
    // Defensive guard: a malformed mzML where the m/z and intensity arrays
    // are not the same length would panic on `intensities[i]` during the
    // permutation step. Skip sorting in that case and let downstream length
    // checks decide what to do with the spectrum.
    if mzs.len() != intensities.len() {
        log::warn!(
            "[unsorted-spectrum] scan_number={} mz/intensity length mismatch ({} vs {}); skipping sort",
            scan_number,
            mzs.len(),
            intensities.len()
        );
        return (mzs, intensities, false);
    }
    // Per-scan trace is `debug` (visible with `--verbose`) so a typical
    // run with ~10 inversions across 200K spectra doesn't spam the default
    // log. The bulk loaders emit a one-line aggregate summary at info level.
    log::debug!(
        "[unsorted-spectrum] scan_number={} n_peaks={}",
        scan_number,
        mzs.len()
    );
    let mut idx: Vec<usize> = (0..mzs.len()).collect();
    idx.sort_by(|&a, &b| mzs[a].total_cmp(&mzs[b]));
    let sorted_mzs: Vec<f64> = idx.iter().map(|&i| mzs[i]).collect();
    let sorted_int: Vec<f32> = idx.iter().map(|&i| intensities[i]).collect();
    (sorted_mzs, sorted_int, true)
}

/// Reader for mzML files
pub struct MzmlReader {
    path: PathBuf,
    reader: MzMLReader<BufReader<File>>,
    total_spectra: Option<usize>,
    current_index: usize,
    /// Pre-parsed f64 isolation window cvParams, keyed by spectrum
    /// index. Overrides mzdata 0.63's f32-quantized bounds in
    /// `convert_spectrum`. See [`read_isolation_cvparams_f64`].
    iso_cv_params: HashMap<u32, IsolationCvParams>,
}

impl MzmlReader {
    /// Open an mzML file for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        // First pass: scan the XML for isolation-window cvParams as
        // f64. Cheap streaming pass (~5% of bytes parsed) that lets us
        // override mzdata 0.63's f32 quantization downstream.
        let iso_cv_params = read_isolation_cvparams_f64(&path)?;
        let file = File::open(&path).map_err(|e| {
            OspreyError::MzmlParseError(format!("Failed to open file '{}': {}", path.display(), e))
        })?;
        let reader = BufReader::new(file);

        let mzml_reader = MzMLReader::new(reader);

        Ok(Self {
            path,
            reader: mzml_reader,
            total_spectra: None, // Will be determined on first iteration
            current_index: 0,
            iso_cv_params,
        })
    }

    /// Convert an mzdata spectrum to an Osprey spectrum
    fn convert_spectrum(
        &self,
        mz_spectrum: mzdata::spectrum::MultiLayerSpectrum,
    ) -> Result<Option<Spectrum>> {
        // Only process MS2 spectra
        if mz_spectrum.description().ms_level != 2 {
            return Ok(None);
        }

        let desc = mz_spectrum.description();

        // Get scan number from index
        let scan_number = desc.index as u32;

        // Get retention time in minutes
        let retention_time = desc.acquisition.first_scan().map_or(0.0, |scan| {
            scan.start_time // mzdata returns time in minutes
        });

        // Get isolation window from precursor (now a Vec in mzdata 0.63).
        // Routes through `make_isolation_window` so we use f64 cvParam
        // values from the pre-parsed XML, bypassing mzdata 0.63's f32
        // quantization. See `read_isolation_cvparams_f64`.
        let isolation_window = match desc.precursor.first() {
            Some(precursor) => {
                match make_isolation_window(precursor, scan_number, Some(&self.iso_cv_params)) {
                    Some(w) => w,
                    None => return Ok(None), // No selected ion
                }
            }
            None => return Ok(None), // No precursor info
        };

        // Get peaks - use the peaks() method which returns RefPeakDataLevel
        let peaks = mz_spectrum.peaks();
        let (mzs, intensities) = match peaks {
            RefPeakDataLevel::Missing => {
                // No peak data available
                return Ok(None);
            }
            RefPeakDataLevel::RawData(arrays) => {
                // mzs() returns a Cow, we need to convert it properly
                let mz_data = arrays.mzs();
                let int_data = arrays.intensities();
                let mzs: Vec<f64> = match mz_data {
                    Ok(cow) => cow.iter().copied().collect(),
                    Err(_) => return Ok(None),
                };
                let intensities: Vec<f32> = match int_data {
                    Ok(cow) => cow.iter().copied().collect(),
                    Err(_) => return Ok(None),
                };
                (mzs, intensities)
            }
            RefPeakDataLevel::Centroid(peaks) => {
                let mzs: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
                let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                (mzs, intensities)
            }
            RefPeakDataLevel::Deconvoluted(peaks) => {
                let mzs: Vec<f64> = peaks.iter().map(|p| p.neutral_mass).collect();
                let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                (mzs, intensities)
            }
        };

        // Iterator path: per-spectrum sort, no aggregation. The `did_sort`
        // bool is discarded here; the bulk loaders below aggregate it
        // across the file instead.
        let (mzs, intensities, _did_sort) = ensure_sorted(mzs, intensities, scan_number);

        Ok(Some(Spectrum {
            scan_number,
            retention_time,
            precursor_mz: isolation_window.center,
            isolation_window,
            mzs,
            intensities,
        }))
    }
}

impl Iterator for MzmlReader {
    type Item = Result<Spectrum>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // mzdata 0.63 iterator returns MultiLayerSpectrum directly, not Result
            match self.reader.next() {
                Some(mz_spectrum) => {
                    self.current_index += 1;
                    match self.convert_spectrum(mz_spectrum) {
                        Ok(Some(spectrum)) => return Some(Ok(spectrum)),
                        Ok(None) => continue, // Skip non-MS2 spectra
                        Err(e) => return Some(Err(e)),
                    }
                }
                None => return None,
            }
        }
    }
}

impl SpectrumSource for MzmlReader {
    fn total_spectra(&self) -> Option<usize> {
        self.total_spectra
    }

    fn reset(&mut self) -> Result<()> {
        // Re-open the file
        *self = Self::open(&self.path)?;
        Ok(())
    }

    fn file_path(&self) -> &Path {
        &self.path
    }
}

// ============================================================================
// MS1 Spectra Loading (for MS1 mass calibration)
// ============================================================================

/// Index of MS1 spectra sorted by retention time for efficient lookup
///
/// Used for finding the nearest MS1 scan to extract isotope envelopes
/// for MS1 mass calibration.
#[derive(Debug, Clone)]
pub struct MS1Index {
    /// MS1 spectra sorted by retention time
    spectra: Vec<MS1Spectrum>,
}

impl MS1Index {
    /// Create an MS1 index from a vector of MS1 spectra
    ///
    /// The spectra will be sorted by retention time.
    pub fn new(mut spectra: Vec<MS1Spectrum>) -> Self {
        spectra.sort_by(|a, b| {
            a.retention_time
                .partial_cmp(&b.retention_time)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        Self { spectra }
    }

    /// Get the number of MS1 spectra
    pub fn len(&self) -> usize {
        self.spectra.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.spectra.is_empty()
    }

    /// Find the nearest MS1 spectrum to a given retention time
    ///
    /// Returns None if no MS1 spectra are available.
    pub fn find_nearest(&self, retention_time: f64) -> Option<&MS1Spectrum> {
        if self.spectra.is_empty() {
            return None;
        }

        // Binary search for the insertion point
        let idx = self
            .spectra
            .partition_point(|s| s.retention_time < retention_time);

        // Check adjacent spectra to find the nearest
        let candidates: Vec<usize> = match idx {
            0 => vec![0],
            n if n >= self.spectra.len() => vec![self.spectra.len() - 1],
            n => vec![n - 1, n],
        };

        candidates
            .into_iter()
            .min_by(|&a, &b| {
                let diff_a = (self.spectra[a].retention_time - retention_time).abs();
                let diff_b = (self.spectra[b].retention_time - retention_time).abs();
                diff_a.total_cmp(&diff_b)
            })
            .map(|i| &self.spectra[i])
    }

    /// Find the nearest MS1 spectrum within a maximum RT tolerance
    ///
    /// Returns None if no MS1 spectrum is within the tolerance.
    pub fn find_nearest_within(
        &self,
        retention_time: f64,
        max_delta_rt: f64,
    ) -> Option<&MS1Spectrum> {
        self.find_nearest(retention_time)
            .filter(|s| (s.retention_time - retention_time).abs() <= max_delta_rt)
    }

    /// Get RT range of MS1 spectra
    pub fn rt_range(&self) -> Option<(f64, f64)> {
        if self.spectra.is_empty() {
            return None;
        }
        Some((
            self.spectra.first().unwrap().retention_time,
            self.spectra.last().unwrap().retention_time,
        ))
    }

    /// Get a reference to the underlying MS1 spectra
    pub fn spectra(&self) -> &[MS1Spectrum] {
        &self.spectra
    }
}

/// Load all spectra (both MS1 and MS2) from an mzML file in a single pass
///
/// This is more efficient than calling MzmlReader and load_ms1_spectra separately,
/// as it only parses the file once.
///
/// Returns: (MS2 spectra, MS1 index)
pub fn load_all_spectra<P: AsRef<Path>>(path: P) -> Result<(Vec<Spectrum>, MS1Index)> {
    let path = path.as_ref();
    // First pass: f64 isolation-window cvParams from raw XML.
    let iso_cv_params = read_isolation_cvparams_f64(path)?;
    let file = File::open(path).map_err(|e| {
        OspreyError::MzmlParseError(format!("Failed to open file '{}': {}", path.display(), e))
    })?;
    let reader = BufReader::new(file);
    let mzml_reader = MzMLReader::new(reader);

    let mut ms2_spectra = Vec::new();
    let mut ms1_spectra = Vec::new();
    // Per-file count of spectra that had non-monotonic centroid pairs and
    // were sorted by `ensure_sorted`. Reported once at end-of-load (info)
    // so the default-verbosity log shows "Sorted N spectra" instead of one
    // line per scan. Per-scan detail is available with `--verbose` via the
    // `log::debug!` inside `ensure_sorted`.
    let mut n_unsorted = 0usize;

    for mz_spectrum in mzml_reader {
        let desc = mz_spectrum.description();
        let ms_level = desc.ms_level;
        let scan_number = desc.index as u32;
        let retention_time = desc.acquisition.first_scan().map_or(0.0, |scan| {
            scan.start_time // mzdata returns time in minutes
        });

        match ms_level {
            1 => {
                // Process MS1 spectrum
                let peaks = mz_spectrum.peaks();
                let (mzs, intensities) = match peaks {
                    RefPeakDataLevel::Missing => continue,
                    RefPeakDataLevel::RawData(arrays) => {
                        let mz_data = arrays.mzs();
                        let int_data = arrays.intensities();
                        let mzs: Vec<f64> = match mz_data {
                            Ok(cow) => cow.iter().copied().collect(),
                            Err(_) => continue,
                        };
                        let intensities: Vec<f32> = match int_data {
                            Ok(cow) => cow.iter().copied().collect(),
                            Err(_) => continue,
                        };
                        (mzs, intensities)
                    }
                    RefPeakDataLevel::Centroid(peaks) => {
                        let mzs: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
                        let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                        (mzs, intensities)
                    }
                    RefPeakDataLevel::Deconvoluted(peaks) => {
                        let mzs: Vec<f64> = peaks.iter().map(|p| p.neutral_mass).collect();
                        let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                        (mzs, intensities)
                    }
                };

                let (mzs, intensities, did_sort) = ensure_sorted(mzs, intensities, scan_number);
                if did_sort {
                    n_unsorted += 1;
                }

                ms1_spectra.push(MS1Spectrum {
                    scan_number,
                    retention_time,
                    mzs,
                    intensities,
                });
            }
            2 => {
                // Process MS2 spectrum.
                // Route through `make_isolation_window` for the f64
                // cvParam override path; see `convert_spectrum` and
                // `read_isolation_cvparams_f64` for the rationale.
                let isolation_window = match desc.precursor.first() {
                    Some(precursor) => {
                        match make_isolation_window(precursor, scan_number, Some(&iso_cv_params)) {
                            Some(w) => w,
                            None => continue, // No selected ion
                        }
                    }
                    None => continue,
                };

                // Get peaks
                let peaks = mz_spectrum.peaks();
                let (mzs, intensities) = match peaks {
                    RefPeakDataLevel::Missing => continue,
                    RefPeakDataLevel::RawData(arrays) => {
                        let mz_data = arrays.mzs();
                        let int_data = arrays.intensities();
                        let mzs: Vec<f64> = match mz_data {
                            Ok(cow) => cow.iter().copied().collect(),
                            Err(_) => continue,
                        };
                        let intensities: Vec<f32> = match int_data {
                            Ok(cow) => cow.iter().copied().collect(),
                            Err(_) => continue,
                        };
                        (mzs, intensities)
                    }
                    RefPeakDataLevel::Centroid(peaks) => {
                        let mzs: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
                        let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                        (mzs, intensities)
                    }
                    RefPeakDataLevel::Deconvoluted(peaks) => {
                        let mzs: Vec<f64> = peaks.iter().map(|p| p.neutral_mass).collect();
                        let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                        (mzs, intensities)
                    }
                };

                // See `convert_spectrum` for the rationale; some mzML
                // producers emit peaks not strictly ascending in m/z, and
                // downstream binary-search consumers need a sorted spectrum.
                let (mzs, intensities, did_sort) = ensure_sorted(mzs, intensities, scan_number);
                if did_sort {
                    n_unsorted += 1;
                }

                ms2_spectra.push(Spectrum {
                    scan_number,
                    retention_time,
                    precursor_mz: isolation_window.center,
                    isolation_window,
                    mzs,
                    intensities,
                });
            }
            _ => {
                // Skip other MS levels (MS3, etc.)
                continue;
            }
        }
    }

    log::info!(
        "Loaded {} MS2 spectra and {} MS1 spectra from '{}'",
        ms2_spectra.len(),
        ms1_spectra.len(),
        path.display()
    );
    if n_unsorted > 0 {
        // Aggregate summary; per-scan detail is at `debug` level.
        log::info!(
            "Sorted {} spectra with non-monotonic centroids (run with --verbose to list scan_numbers)",
            n_unsorted
        );
    }

    Ok((ms2_spectra, MS1Index::new(ms1_spectra)))
}

/// Load all MS1 spectra from an mzML file
///
/// This is used for MS1 mass calibration (extracting M+0 isotope peaks).
/// Returns an MS1Index for efficient nearest-neighbor lookup.
///
/// Note: If you need both MS1 and MS2, use `load_all_spectra()` instead
/// to avoid parsing the file twice.
pub fn load_ms1_spectra<P: AsRef<Path>>(path: P) -> Result<MS1Index> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        OspreyError::MzmlParseError(format!("Failed to open file '{}': {}", path.display(), e))
    })?;
    let reader = BufReader::new(file);
    let mzml_reader = MzMLReader::new(reader);

    let mut ms1_spectra = Vec::new();
    let mut n_unsorted = 0usize;

    for mz_spectrum in mzml_reader {
        // Only process MS1 spectra
        if mz_spectrum.description().ms_level != 1 {
            continue;
        }

        let desc = mz_spectrum.description();
        let scan_number = desc.index as u32;
        let retention_time = desc.acquisition.first_scan().map_or(0.0, |scan| {
            scan.start_time // mzdata returns time in minutes
        });

        // Get peaks
        let peaks = mz_spectrum.peaks();
        let (mzs, intensities) = match peaks {
            RefPeakDataLevel::Missing => continue,
            RefPeakDataLevel::RawData(arrays) => {
                let mz_data = arrays.mzs();
                let int_data = arrays.intensities();
                let mzs: Vec<f64> = match mz_data {
                    Ok(cow) => cow.iter().copied().collect(),
                    Err(_) => continue,
                };
                let intensities: Vec<f32> = match int_data {
                    Ok(cow) => cow.iter().copied().collect(),
                    Err(_) => continue,
                };
                (mzs, intensities)
            }
            RefPeakDataLevel::Centroid(peaks) => {
                let mzs: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
                let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                (mzs, intensities)
            }
            RefPeakDataLevel::Deconvoluted(peaks) => {
                let mzs: Vec<f64> = peaks.iter().map(|p| p.neutral_mass).collect();
                let intensities: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
                (mzs, intensities)
            }
        };

        let (mzs, intensities, did_sort) = ensure_sorted(mzs, intensities, scan_number);
        if did_sort {
            n_unsorted += 1;
        }

        ms1_spectra.push(MS1Spectrum {
            scan_number,
            retention_time,
            mzs,
            intensities,
        });
    }

    log::info!(
        "Loaded {} MS1 spectra from '{}'",
        ms1_spectra.len(),
        path.display()
    );
    if n_unsorted > 0 {
        log::info!(
            "Sorted {} spectra with non-monotonic centroids (run with --verbose to list scan_numbers)",
            n_unsorted
        );
    }

    Ok(MS1Index::new(ms1_spectra))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies that a symmetric isolation window has the correct center m/z and total width.
    #[test]
    fn test_isolation_window() {
        let window = IsolationWindow::symmetric(500.0, 12.5);
        assert!((window.center - 500.0).abs() < 1e-6);
        assert!((window.width() - 25.0).abs() < 1e-6);
    }

    /// Already-sorted input: ensure_sorted returns the inputs unchanged
    /// (fast path; no permutation, no extra allocations). `did_sort` must
    /// be false so the bulk loaders don't credit this spectrum to the
    /// "sorted N spectra" summary log line.
    #[test]
    fn ensure_sorted_already_sorted_fast_path() {
        let mzs = vec![100.0, 200.0, 300.0];
        let intensities = vec![1.0_f32, 2.0, 3.0];
        let (sorted_mzs, sorted_int, did_sort) = ensure_sorted(mzs.clone(), intensities.clone(), 1);
        assert_eq!(sorted_mzs, mzs);
        assert_eq!(sorted_int, intensities);
        assert!(!did_sort, "fast path must report did_sort = false");
    }

    /// Single inversion: both arrays reorder consistently so each
    /// intensity stays paired with its original m/z value. `did_sort`
    /// must be true so the bulk loaders increment their summary counter.
    #[test]
    fn ensure_sorted_inversion_reorders_both_arrays() {
        let mzs = vec![100.0, 300.0, 200.0];
        let intensities = vec![1.0_f32, 3.0, 2.0];
        let (sorted_mzs, sorted_int, did_sort) = ensure_sorted(mzs, intensities, 1);
        assert_eq!(sorted_mzs, vec![100.0, 200.0, 300.0]);
        assert_eq!(sorted_int, vec![1.0_f32, 2.0, 3.0]);
        assert!(did_sort, "inversion-path must report did_sort = true");
    }

    /// Stable sort: equal m/z values keep their original relative order
    /// (matches `slice::sort_by` and the C# port's stable LINQ OrderBy).
    /// The deliberate inversion at the end (200.0 before 100.0) forces
    /// the sort path; the duplicate-m/z pair at the front must retain
    /// input order in the result.
    #[test]
    fn ensure_sorted_equal_mz_preserves_order() {
        let mzs = vec![100.0, 100.0, 200.0, 100.0];
        let intensities = vec![1.0_f32, 2.0, 3.0, 4.0];
        let (sorted_mzs, sorted_int, did_sort) = ensure_sorted(mzs, intensities, 1);
        assert_eq!(sorted_mzs, vec![100.0, 100.0, 100.0, 200.0]);
        assert_eq!(sorted_int, vec![1.0_f32, 2.0, 4.0, 3.0]);
        assert!(did_sort);
    }

    /// Length mismatch (malformed mzML): skip sorting and return the
    /// inputs untouched so downstream length checks can act on the
    /// malformed spectrum without panicking on out-of-bounds indexing.
    /// `did_sort` is false because we did not actually sort (we bailed).
    #[test]
    fn ensure_sorted_length_mismatch_skips() {
        let mzs = vec![300.0, 100.0, 200.0];
        let intensities = vec![1.0_f32, 2.0]; // shorter than mzs
        let (out_mzs, out_int, did_sort) = ensure_sorted(mzs.clone(), intensities.clone(), 1);
        assert_eq!(out_mzs, mzs);
        assert_eq!(out_int, intensities);
        assert!(
            !did_sort,
            "length mismatch skips the sort, so did_sort = false"
        );
    }
}
