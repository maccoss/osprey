//! Binary library cache for fast re-loading.
//!
//! Saves/loads `Vec<LibraryEntry>` in a compact binary format, avoiding the
//! overhead of TSV/blib/elib parsing on repeated runs.

use osprey_core::{
    FragmentAnnotation, IonType, LibraryEntry, LibraryFragment, Modification, NeutralLoss,
    OspreyError, Result,
};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Magic bytes at the start of every cache file.
const MAGIC: &[u8; 8] = b"OSPRLBR\0";
/// Current cache format version.
const VERSION: u32 = 1;

/// Save library entries to a binary cache file.
pub fn save_library_cache(entries: &[LibraryEntry], path: &Path) -> Result<()> {
    let file = std::fs::File::create(path).map_err(|e| {
        OspreyError::LibraryLoadError(format!("Failed to create cache file: {}", e))
    })?;
    let mut w = BufWriter::new(file);

    w.write_all(MAGIC).map_err(write_err)?;
    w.write_all(&VERSION.to_le_bytes()).map_err(write_err)?;
    w.write_all(&(entries.len() as u64).to_le_bytes())
        .map_err(write_err)?;

    for entry in entries {
        w.write_all(&entry.id.to_le_bytes()).map_err(write_err)?;
        write_string(&mut w, &entry.sequence)?;
        write_string(&mut w, &entry.modified_sequence)?;
        w.write_all(&[entry.charge]).map_err(write_err)?;
        w.write_all(&entry.precursor_mz.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&entry.retention_time.to_le_bytes())
            .map_err(write_err)?;
        w.write_all(&[entry.rt_calibrated as u8])
            .map_err(write_err)?;
        w.write_all(&[entry.is_decoy as u8]).map_err(write_err)?;

        // Modifications
        w.write_all(&(entry.modifications.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        for m in &entry.modifications {
            w.write_all(&(m.position as u32).to_le_bytes())
                .map_err(write_err)?;
            match m.unimod_id {
                Some(id) => {
                    w.write_all(&[1]).map_err(write_err)?;
                    w.write_all(&id.to_le_bytes()).map_err(write_err)?;
                }
                None => {
                    w.write_all(&[0]).map_err(write_err)?;
                }
            }
            w.write_all(&m.mass_delta.to_le_bytes())
                .map_err(write_err)?;
            match &m.name {
                Some(name) => {
                    w.write_all(&[1]).map_err(write_err)?;
                    write_string(&mut w, name)?;
                }
                None => {
                    w.write_all(&[0]).map_err(write_err)?;
                }
            }
        }

        // Fragments
        w.write_all(&(entry.fragments.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        for frag in &entry.fragments {
            w.write_all(&frag.mz.to_le_bytes()).map_err(write_err)?;
            w.write_all(&frag.relative_intensity.to_le_bytes())
                .map_err(write_err)?;
            w.write_all(&[ion_type_to_u8(frag.annotation.ion_type)])
                .map_err(write_err)?;
            w.write_all(&[frag.annotation.ordinal]).map_err(write_err)?;
            w.write_all(&[frag.annotation.charge]).map_err(write_err)?;
            write_neutral_loss(&mut w, &frag.annotation.neutral_loss)?;
        }

        // Protein IDs
        w.write_all(&(entry.protein_ids.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        for pid in &entry.protein_ids {
            write_string(&mut w, pid)?;
        }

        // Gene names
        w.write_all(&(entry.gene_names.len() as u32).to_le_bytes())
            .map_err(write_err)?;
        for gn in &entry.gene_names {
            write_string(&mut w, gn)?;
        }
    }

    w.flush().map_err(write_err)?;
    Ok(())
}

/// Load library entries from a binary cache file.
pub fn load_library_cache(path: &Path) -> Result<Vec<LibraryEntry>> {
    let file = std::fs::File::open(path)
        .map_err(|e| OspreyError::LibraryLoadError(format!("Failed to open cache file: {}", e)))?;
    let mut r = BufReader::new(file);

    // Validate magic
    let mut magic = [0u8; 8];
    r.read_exact(&mut magic).map_err(read_err)?;
    if &magic != MAGIC {
        return Err(OspreyError::LibraryLoadError(
            "Invalid library cache file (bad magic)".into(),
        ));
    }

    let version = read_u32(&mut r)?;
    if version != VERSION {
        return Err(OspreyError::LibraryLoadError(format!(
            "Unsupported cache version {} (expected {})",
            version, VERSION
        )));
    }

    let count = read_u64(&mut r)? as usize;
    let mut entries = Vec::with_capacity(count);

    for _ in 0..count {
        let id = read_u32(&mut r)?;
        let sequence = read_string(&mut r)?;
        let modified_sequence = read_string(&mut r)?;
        let charge = read_u8(&mut r)?;
        let precursor_mz = read_f64(&mut r)?;
        let retention_time = read_f64(&mut r)?;
        let rt_calibrated = read_u8(&mut r)? != 0;
        let is_decoy = read_u8(&mut r)? != 0;

        // Modifications
        let n_mods = read_u32(&mut r)? as usize;
        let mut modifications = Vec::with_capacity(n_mods);
        for _ in 0..n_mods {
            let position = read_u32(&mut r)? as usize;
            let has_unimod = read_u8(&mut r)? != 0;
            let unimod_id = if has_unimod {
                Some(read_u32(&mut r)?)
            } else {
                None
            };
            let mass_delta = read_f64(&mut r)?;
            let has_name = read_u8(&mut r)? != 0;
            let name = if has_name {
                Some(read_string(&mut r)?)
            } else {
                None
            };
            modifications.push(Modification {
                position,
                unimod_id,
                mass_delta,
                name,
            });
        }

        // Fragments
        let n_frags = read_u32(&mut r)? as usize;
        let mut fragments = Vec::with_capacity(n_frags);
        for _ in 0..n_frags {
            let mz = read_f64(&mut r)?;
            let relative_intensity = read_f32(&mut r)?;
            let ion_type = u8_to_ion_type(read_u8(&mut r)?);
            let ordinal = read_u8(&mut r)?;
            let charge = read_u8(&mut r)?;
            let neutral_loss = read_neutral_loss(&mut r)?;
            fragments.push(LibraryFragment {
                mz,
                relative_intensity,
                annotation: FragmentAnnotation {
                    ion_type,
                    ordinal,
                    charge,
                    neutral_loss,
                },
            });
        }

        // Protein IDs
        let n_proteins = read_u32(&mut r)? as usize;
        let mut protein_ids = Vec::with_capacity(n_proteins);
        for _ in 0..n_proteins {
            protein_ids.push(read_string(&mut r)?);
        }

        // Gene names
        let n_genes = read_u32(&mut r)? as usize;
        let mut gene_names = Vec::with_capacity(n_genes);
        for _ in 0..n_genes {
            gene_names.push(read_string(&mut r)?);
        }

        entries.push(LibraryEntry {
            id,
            sequence,
            modified_sequence,
            modifications,
            charge,
            precursor_mz,
            retention_time,
            rt_calibrated,
            is_decoy,
            fragments,
            protein_ids,
            gene_names,
        });
    }

    Ok(entries)
}

// ---------- helpers ----------

fn write_err(e: std::io::Error) -> OspreyError {
    OspreyError::LibraryLoadError(format!("Cache write error: {}", e))
}

fn read_err(e: std::io::Error) -> OspreyError {
    OspreyError::LibraryLoadError(format!("Cache read error: {}", e))
}

fn write_string(w: &mut impl Write, s: &str) -> Result<()> {
    let bytes = s.as_bytes();
    w.write_all(&(bytes.len() as u32).to_le_bytes())
        .map_err(write_err)?;
    w.write_all(bytes).map_err(write_err)?;
    Ok(())
}

fn read_string(r: &mut impl Read) -> Result<String> {
    let len = read_u32(r)? as usize;
    let mut buf = vec![0u8; len];
    r.read_exact(&mut buf).map_err(read_err)?;
    String::from_utf8(buf)
        .map_err(|e| OspreyError::LibraryLoadError(format!("Invalid UTF-8 in cache: {}", e)))
}

fn read_u8(r: &mut impl Read) -> Result<u8> {
    let mut buf = [0u8; 1];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(buf[0])
}

fn read_u32(r: &mut impl Read) -> Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64(r: &mut impl Read) -> Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_f32(r: &mut impl Read) -> Result<f32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(f32::from_le_bytes(buf))
}

fn read_f64(r: &mut impl Read) -> Result<f64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf).map_err(read_err)?;
    Ok(f64::from_le_bytes(buf))
}

fn ion_type_to_u8(ion: IonType) -> u8 {
    match ion {
        IonType::A => 0,
        IonType::B => 1,
        IonType::C => 2,
        IonType::X => 3,
        IonType::Y => 4,
        IonType::Z => 5,
        IonType::Precursor => 6,
        IonType::Immonium => 7,
        IonType::Internal => 8,
        IonType::Unknown => 9,
    }
}

fn u8_to_ion_type(v: u8) -> IonType {
    match v {
        0 => IonType::A,
        1 => IonType::B,
        2 => IonType::C,
        3 => IonType::X,
        4 => IonType::Y,
        5 => IonType::Z,
        6 => IonType::Precursor,
        7 => IonType::Immonium,
        8 => IonType::Internal,
        _ => IonType::Unknown,
    }
}

fn write_neutral_loss(w: &mut impl Write, nl: &Option<NeutralLoss>) -> Result<()> {
    match nl {
        None => w.write_all(&[0]).map_err(write_err)?,
        Some(NeutralLoss::H2O) => w.write_all(&[1]).map_err(write_err)?,
        Some(NeutralLoss::NH3) => w.write_all(&[2]).map_err(write_err)?,
        Some(NeutralLoss::H3PO4) => w.write_all(&[3]).map_err(write_err)?,
        Some(NeutralLoss::Custom(mass)) => {
            w.write_all(&[4]).map_err(write_err)?;
            w.write_all(&mass.to_le_bytes()).map_err(write_err)?;
        }
    }
    Ok(())
}

fn read_neutral_loss(r: &mut impl Read) -> Result<Option<NeutralLoss>> {
    match read_u8(r)? {
        0 => Ok(None),
        1 => Ok(Some(NeutralLoss::H2O)),
        2 => Ok(Some(NeutralLoss::NH3)),
        3 => Ok(Some(NeutralLoss::H3PO4)),
        4 => {
            let mass = read_f64(r)?;
            Ok(Some(NeutralLoss::Custom(mass)))
        }
        tag => Err(OspreyError::LibraryLoadError(format!(
            "Unknown neutral loss tag: {}",
            tag
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_entry(id: u32) -> LibraryEntry {
        LibraryEntry {
            id,
            sequence: "PEPTIDEK".into(),
            modified_sequence: "PEPTIDEK".into(),
            modifications: vec![Modification {
                position: 0,
                unimod_id: Some(1),
                mass_delta: 42.010565,
                name: Some("Acetyl".into()),
            }],
            charge: 2,
            precursor_mz: 471.2567,
            retention_time: 25.3,
            rt_calibrated: false,
            is_decoy: false,
            fragments: vec![
                LibraryFragment {
                    mz: 147.1128,
                    relative_intensity: 1.0,
                    annotation: FragmentAnnotation {
                        ion_type: IonType::Y,
                        ordinal: 1,
                        charge: 1,
                        neutral_loss: None,
                    },
                },
                LibraryFragment {
                    mz: 262.1398,
                    relative_intensity: 0.5,
                    annotation: FragmentAnnotation {
                        ion_type: IonType::B,
                        ordinal: 3,
                        charge: 1,
                        neutral_loss: Some(NeutralLoss::H2O),
                    },
                },
                LibraryFragment {
                    mz: 500.0,
                    relative_intensity: 0.3,
                    annotation: FragmentAnnotation {
                        ion_type: IonType::Y,
                        ordinal: 5,
                        charge: 2,
                        neutral_loss: Some(NeutralLoss::Custom(98.0)),
                    },
                },
            ],
            protein_ids: vec!["P12345".into(), "Q67890".into()],
            gene_names: vec!["GENE1".into()],
        }
    }

    #[test]
    fn test_library_cache_round_trip() {
        let entries = vec![make_test_entry(0), make_test_entry(1)];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.libcache");

        save_library_cache(&entries, &path).unwrap();
        let loaded = load_library_cache(&path).unwrap();

        assert_eq!(loaded.len(), 2);
        for (orig, loaded) in entries.iter().zip(loaded.iter()) {
            assert_eq!(orig.id, loaded.id);
            assert_eq!(orig.sequence, loaded.sequence);
            assert_eq!(orig.modified_sequence, loaded.modified_sequence);
            assert_eq!(orig.charge, loaded.charge);
            assert!((orig.precursor_mz - loaded.precursor_mz).abs() < 1e-10);
            assert!((orig.retention_time - loaded.retention_time).abs() < 1e-10);
            assert_eq!(orig.rt_calibrated, loaded.rt_calibrated);
            assert_eq!(orig.is_decoy, loaded.is_decoy);
            assert_eq!(orig.fragments.len(), loaded.fragments.len());
            assert_eq!(orig.modifications.len(), loaded.modifications.len());
            assert_eq!(orig.protein_ids, loaded.protein_ids);
            assert_eq!(orig.gene_names, loaded.gene_names);

            // Check fragment details
            for (of, lf) in orig.fragments.iter().zip(loaded.fragments.iter()) {
                assert!((of.mz - lf.mz).abs() < 1e-10);
                assert!((of.relative_intensity - lf.relative_intensity).abs() < 1e-6);
                assert_eq!(of.annotation.ion_type, lf.annotation.ion_type);
                assert_eq!(of.annotation.ordinal, lf.annotation.ordinal);
                assert_eq!(of.annotation.charge, lf.annotation.charge);
                assert_eq!(of.annotation.neutral_loss, lf.annotation.neutral_loss);
            }

            // Check modification details
            for (om, lm) in orig.modifications.iter().zip(loaded.modifications.iter()) {
                assert_eq!(om.position, lm.position);
                assert_eq!(om.unimod_id, lm.unimod_id);
                assert!((om.mass_delta - lm.mass_delta).abs() < 1e-10);
                assert_eq!(om.name, lm.name);
            }
        }
    }

    #[test]
    fn test_library_cache_empty() {
        let entries: Vec<LibraryEntry> = vec![];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.libcache");

        save_library_cache(&entries, &path).unwrap();
        let loaded = load_library_cache(&path).unwrap();
        assert!(loaded.is_empty());
    }

    #[test]
    fn test_library_cache_invalid_magic() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.libcache");
        std::fs::write(&path, b"NOTVALID").unwrap();

        let result = load_library_cache(&path);
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(err_msg.contains("bad magic"));
    }
}
