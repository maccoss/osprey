//! Binary spectra cache for fast re-loading.
//!
//! After parsing mzML (which takes ~5 minutes for large Astral files),
//! this module writes a compact `.spectra.bin` file that can be reloaded
//! in seconds. Used by the post-FDR re-scoring phase to avoid re-parsing mzML.
//!
//! Format (little-endian):
//! ```text
//! [magic: 8 bytes "OSPRSPC\0"]
//! [version: u32]
//! [n_ms2: u32]
//! [n_ms1: u32]
//! For each MS2 spectrum:
//!   [scan_number: u32] [retention_time: f64] [precursor_mz: f64]
//!   [iso_center: f64] [iso_lower: f64] [iso_upper: f64]
//!   [n_peaks: u32] [mzs: f64 × n_peaks] [intensities: f32 × n_peaks]
//! For each MS1 spectrum:
//!   [scan_number: u32] [retention_time: f64]
//!   [n_peaks: u32] [mzs: f64 × n_peaks] [intensities: f32 × n_peaks]
//! ```

use osprey_core::{IsolationWindow, MS1Spectrum, OspreyError, Result, Spectrum};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use super::MS1Index;

const MAGIC: &[u8; 8] = b"OSPRSPC\0";
const VERSION: u32 = 1;

/// Save spectra to a binary cache file for fast reload.
pub fn save_spectra_cache(path: &Path, spectra: &[Spectrum], ms1_index: &MS1Index) -> Result<()> {
    // Write to a LOCAL temp file first, then move to final destination.
    // This avoids corrupt/0-byte files on network filesystems (NAS).
    let tmp_path = std::env::temp_dir().join(format!("osprey_spectra_{}.bin", std::process::id()));
    let file = File::create(&tmp_path).map_err(|e| {
        OspreyError::IoError(std::io::Error::other(format!(
            "Failed to create spectra cache '{}': {}",
            tmp_path.display(),
            e
        )))
    })?;
    let mut w = BufWriter::new(file);

    let ms1_spectra = ms1_index.spectra();

    // Header
    w.write_all(MAGIC).map_err(write_err)?;
    w.write_all(&VERSION.to_le_bytes()).map_err(write_err)?;
    w.write_all(&(spectra.len() as u32).to_le_bytes())
        .map_err(write_err)?;
    w.write_all(&(ms1_spectra.len() as u32).to_le_bytes())
        .map_err(write_err)?;

    // MS2 spectra
    for s in spectra {
        w.write_all(&s.scan_number.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.retention_time.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.precursor_mz.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.isolation_window.center.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.isolation_window.lower_offset.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.isolation_window.upper_offset.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&(s.mzs.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        // Write mzs as raw bytes
        let mz_bytes: &[u8] =
            unsafe { std::slice::from_raw_parts(s.mzs.as_ptr() as *const u8, s.mzs.len() * 8) };
        w.write_all(mz_bytes).map_err(write_err)?;
        // Write intensities as raw bytes
        let int_bytes: &[u8] = unsafe {
            std::slice::from_raw_parts(s.intensities.as_ptr() as *const u8, s.intensities.len() * 4)
        };
        w.write_all(int_bytes).map_err(write_err)?;
    }

    // MS1 spectra
    for s in ms1_spectra {
        w.write_all(&s.scan_number.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&s.retention_time.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&(s.mzs.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        let mz_bytes: &[u8] =
            unsafe { std::slice::from_raw_parts(s.mzs.as_ptr() as *const u8, s.mzs.len() * 8) };
        w.write_all(mz_bytes).map_err(write_err)?;
        let int_bytes: &[u8] = unsafe {
            std::slice::from_raw_parts(s.intensities.as_ptr() as *const u8, s.intensities.len() * 4)
        };
        w.write_all(int_bytes).map_err(write_err)?;
    }

    w.flush().map_err(write_err)?;
    drop(w); // close file before move

    osprey_core::move_file_safe(&tmp_path, path)?;
    Ok(())
}

/// Load spectra from a binary cache file.
pub fn load_spectra_cache(path: &Path) -> Result<(Vec<Spectrum>, MS1Index)> {
    let file = File::open(path).map_err(|e| {
        OspreyError::IoError(std::io::Error::other(format!(
            "Failed to open spectra cache '{}': {}",
            path.display(),
            e
        )))
    })?;
    let mut r = BufReader::new(file);

    // Validate header
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic).map_err(read_err)?;
    if &magic != MAGIC {
        return Err(OspreyError::IoError(std::io::Error::other(
            "Invalid spectra cache magic bytes",
        )));
    }

    let version = read_u32(&mut r)?;
    if version != VERSION {
        return Err(OspreyError::IoError(std::io::Error::other(format!(
            "Unsupported spectra cache version: {} (expected {})",
            version, VERSION
        ))));
    }

    let n_ms2 = read_u32(&mut r)? as usize;
    let n_ms1 = read_u32(&mut r)? as usize;

    // Read MS2 spectra
    let mut spectra = Vec::with_capacity(n_ms2);
    for _ in 0..n_ms2 {
        let scan_number = read_u32(&mut r)?;
        let retention_time = read_f64(&mut r)?;
        let precursor_mz = read_f64(&mut r)?;
        let iso_center = read_f64(&mut r)?;
        let iso_lower = read_f64(&mut r)?;
        let iso_upper = read_f64(&mut r)?;
        let n_peaks = read_u32(&mut r)? as usize;

        let mut mzs = vec![0f64; n_peaks];
        let mz_bytes: &mut [u8] =
            unsafe { std::slice::from_raw_parts_mut(mzs.as_mut_ptr() as *mut u8, n_peaks * 8) };
        r.read_exact(mz_bytes).map_err(read_err)?;

        let mut intensities = vec![0f32; n_peaks];
        let int_bytes: &mut [u8] = unsafe {
            std::slice::from_raw_parts_mut(intensities.as_mut_ptr() as *mut u8, n_peaks * 4)
        };
        r.read_exact(int_bytes).map_err(read_err)?;

        spectra.push(Spectrum {
            scan_number,
            retention_time,
            precursor_mz,
            isolation_window: IsolationWindow::new(iso_center, iso_lower, iso_upper),
            mzs,
            intensities,
        });
    }

    // Read MS1 spectra
    let mut ms1_spectra = Vec::with_capacity(n_ms1);
    for _ in 0..n_ms1 {
        let scan_number = read_u32(&mut r)?;
        let retention_time = read_f64(&mut r)?;
        let n_peaks = read_u32(&mut r)? as usize;

        let mut mzs = vec![0f64; n_peaks];
        let mz_bytes: &mut [u8] =
            unsafe { std::slice::from_raw_parts_mut(mzs.as_mut_ptr() as *mut u8, n_peaks * 8) };
        r.read_exact(mz_bytes).map_err(read_err)?;

        let mut intensities = vec![0f32; n_peaks];
        let int_bytes: &mut [u8] = unsafe {
            std::slice::from_raw_parts_mut(intensities.as_mut_ptr() as *mut u8, n_peaks * 4)
        };
        r.read_exact(int_bytes).map_err(read_err)?;

        ms1_spectra.push(MS1Spectrum {
            scan_number,
            retention_time,
            mzs,
            intensities,
        });
    }

    Ok((spectra, MS1Index::new(ms1_spectra)))
}

/// Get the spectra cache path for a given input file.
pub fn spectra_cache_path(input_file: &Path) -> std::path::PathBuf {
    input_file.with_extension("spectra.bin")
}

fn write_err(e: std::io::Error) -> OspreyError {
    OspreyError::IoError(e)
}

fn read_err(e: std::io::Error) -> OspreyError {
    OspreyError::IoError(e)
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_f64<R: Read>(r: &mut R) -> Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(f64::from_le_bytes(buf))
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_spectra_cache_round_trip() {
        let spectra = vec![
            Spectrum {
                scan_number: 1,
                retention_time: 10.5,
                precursor_mz: 500.0,
                isolation_window: IsolationWindow::new(500.0, 1.5, 1.5),
                mzs: vec![100.0, 200.0, 300.0],
                intensities: vec![1000.0, 2000.0, 500.0],
            },
            Spectrum {
                scan_number: 2,
                retention_time: 11.0,
                precursor_mz: 600.0,
                isolation_window: IsolationWindow::new(600.0, 2.0, 3.0),
                mzs: vec![150.0, 250.0],
                intensities: vec![800.0, 1200.0],
            },
        ];

        let ms1_spectra = vec![MS1Spectrum {
            scan_number: 0,
            retention_time: 10.0,
            mzs: vec![400.0, 500.0, 600.0],
            intensities: vec![5000.0, 3000.0, 1000.0],
        }];
        let ms1_index = MS1Index::new(ms1_spectra);

        let tmp = NamedTempFile::new().unwrap();
        save_spectra_cache(tmp.path(), &spectra, &ms1_index).unwrap();
        let (loaded_spectra, loaded_ms1) = load_spectra_cache(tmp.path()).unwrap();

        assert_eq!(loaded_spectra.len(), 2);
        assert_eq!(loaded_spectra[0].scan_number, 1);
        assert!((loaded_spectra[0].retention_time - 10.5).abs() < 1e-10);
        assert_eq!(loaded_spectra[0].mzs, vec![100.0, 200.0, 300.0]);
        assert_eq!(loaded_spectra[0].intensities, vec![1000.0, 2000.0, 500.0]);
        assert!((loaded_spectra[1].isolation_window.lower_offset - 2.0).abs() < 1e-10);
        assert!((loaded_spectra[1].isolation_window.upper_offset - 3.0).abs() < 1e-10);

        assert_eq!(loaded_ms1.len(), 1);
    }

    #[test]
    fn test_empty_spectra_cache() {
        let tmp = NamedTempFile::new().unwrap();
        let ms1_index = MS1Index::new(vec![]);
        save_spectra_cache(tmp.path(), &[], &ms1_index).unwrap();
        let (spectra, ms1) = load_spectra_cache(tmp.path()).unwrap();
        assert!(spectra.is_empty());
        assert!(ms1.is_empty());
    }
}
