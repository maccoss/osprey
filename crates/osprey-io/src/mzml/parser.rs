//! mzML file parser using mzdata crate
//!
//! Converts mzdata spectrum representations to Osprey's internal types.

use mzdata::io::mzml::MzMLReader;
use mzdata::prelude::*;
use mzdata::spectrum::RefPeakDataLevel;
use osprey_core::{IsolationWindow, OspreyError, Result, Spectrum, SpectrumSource};
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isolation_window() {
        let window = IsolationWindow::symmetric(500.0, 12.5);
        assert!((window.center - 500.0).abs() < 1e-6);
        assert!((window.width() - 25.0).abs() < 1e-6);
    }
}
