//! mzML file parser using mzdata crate
//!
//! Converts mzdata spectrum representations to Osprey's internal types.
//!
//! ## MS1 Spectra for Mass Calibration
//!
//! pyXcorrDIA extracts the M+0 isotope peak from MS1 spectra for accurate mass calibration.
//! This module provides:
//! - `load_all_spectra()` - Load both MS1 and MS2 spectra in a single pass (most efficient)
//! - `load_ms1_spectra()` - Load only MS1 spectra for isotope extraction
//! - `MS1Index` - Index MS1 spectra by retention time for efficient lookup

use mzdata::io::mzml::MzMLReader;
use mzdata::prelude::*;
use mzdata::spectrum::RefPeakDataLevel;
use osprey_core::{IsolationWindow, MS1Spectrum, OspreyError, Result, Spectrum, SpectrumSource};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

/// Reader for mzML files
pub struct MzmlReader {
    path: PathBuf,
    reader: MzMLReader<BufReader<File>>,
    total_spectra: Option<usize>,
    current_index: usize,
}

impl MzmlReader {
    /// Open an mzML file for reading
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
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

        // Get isolation window from precursor (now a Vec in mzdata 0.63)
        let isolation_window = if let Some(precursor) = desc.precursor.first() {
            let ion = match precursor.ion() {
                Some(i) => i,
                None => return Ok(None),
            };
            let isolation = &precursor.isolation_window;

            let center = ion.mz;
            // In mzdata 0.63, bounds are f32 directly
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

            IsolationWindow::new(center, lower_offset, upper_offset)
        } else {
            // No precursor info - skip this spectrum
            return Ok(None);
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
// MS1 Spectra Loading (for MS1 mass calibration - pyXcorrDIA compatible)
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
    pub fn find_nearest_within(&self, retention_time: f64, max_delta_rt: f64) -> Option<&MS1Spectrum> {
        self.find_nearest(retention_time).filter(|s| {
            (s.retention_time - retention_time).abs() <= max_delta_rt
        })
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
}

/// Load all spectra (both MS1 and MS2) from an mzML file in a single pass
///
/// This is more efficient than calling MzmlReader and load_ms1_spectra separately,
/// as it only parses the file once.
///
/// Returns: (MS2 spectra, MS1 index)
pub fn load_all_spectra<P: AsRef<Path>>(path: P) -> Result<(Vec<Spectrum>, MS1Index)> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        OspreyError::MzmlParseError(format!("Failed to open file '{}': {}", path.display(), e))
    })?;
    let reader = BufReader::new(file);
    let mzml_reader = MzMLReader::new(reader);

    let mut ms2_spectra = Vec::new();
    let mut ms1_spectra = Vec::new();

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

                ms1_spectra.push(MS1Spectrum {
                    scan_number,
                    retention_time,
                    mzs,
                    intensities,
                });
            }
            2 => {
                // Process MS2 spectrum
                // Get isolation window from precursor
                let isolation_window = if let Some(precursor) = desc.precursor.first() {
                    let ion = match precursor.ion() {
                        Some(i) => i,
                        None => continue,
                    };
                    let isolation = &precursor.isolation_window;

                    let center = ion.mz;
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

                    IsolationWindow::new(center, lower_offset, upper_offset)
                } else {
                    continue;
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
}
