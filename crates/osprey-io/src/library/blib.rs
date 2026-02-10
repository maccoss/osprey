//! BiblioSpec blib format loader
//!
//! Parses spectral libraries in BiblioSpec's SQLite-based blib format.
//! This format is commonly used by Skyline and other proteomics tools.

use osprey_core::{
    FragmentAnnotation, IonType, LibraryEntry, LibraryFragment, Modification, OspreyError, Result,
};
use rusqlite::{Connection, OpenFlags};
use std::path::Path;

/// BiblioSpec blib format loader
pub struct BlibLoader;

impl BlibLoader {
    pub fn new() -> Self {
        Self
    }

    /// Load library entries from a blib file
    pub fn load<P: AsRef<Path>>(&self, path: P) -> Result<Vec<LibraryEntry>> {
        let path = path.as_ref();

        // Open in read-only mode
        let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)
            .map_err(|e| OspreyError::LibraryLoadError(format!("Failed to open blib: {}", e)))?;

        self.load_from_connection(&conn)
    }

    /// Load entries from an open database connection
    fn load_from_connection(&self, conn: &Connection) -> Result<Vec<LibraryEntry>> {
        // Check if this is a valid blib file
        if !self.table_exists(conn, "RefSpectra")? {
            return Err(OspreyError::LibraryLoadError(
                "Invalid blib file: RefSpectra table not found".to_string(),
            ));
        }

        // Query all spectra with their peaks
        let mut stmt = conn
            .prepare(
                r#"
                SELECT
                    r.id,
                    r.peptideSeq,
                    r.peptideModSeq,
                    r.precursorMZ,
                    r.precursorCharge,
                    r.retentionTime,
                    r.numPeaks,
                    p.peakMZ,
                    p.peakIntensity
                FROM RefSpectra r
                LEFT JOIN RefSpectraPeaks p ON r.id = p.RefSpectraID
                ORDER BY r.id
                "#,
            )
            .map_err(|e| OspreyError::LibraryLoadError(format!("SQL prepare error: {}", e)))?;

        let mut entries: Vec<LibraryEntry> = Vec::new();
        let mut current_id: Option<i64> = None;
        let mut current_entry: Option<LibraryEntry> = None;

        let rows = stmt
            .query_map([], |row| {
                let id: i64 = row.get(0)?;
                let peptide_seq: String = row.get(1)?;
                let peptide_mod_seq: String = row.get(2)?;
                let precursor_mz: f64 = row.get(3)?;
                let precursor_charge: i32 = row.get(4)?;
                let retention_time: Option<f64> = row.get(5)?;
                let _num_peaks: Option<i32> = row.get(6)?;
                let peak_mz_blob: Option<Vec<u8>> = row.get(7)?;
                let peak_int_blob: Option<Vec<u8>> = row.get(8)?;

                Ok((
                    id,
                    peptide_seq,
                    peptide_mod_seq,
                    precursor_mz,
                    precursor_charge,
                    retention_time,
                    peak_mz_blob,
                    peak_int_blob,
                ))
            })
            .map_err(|e| OspreyError::LibraryLoadError(format!("Query error: {}", e)))?;

        for row_result in rows {
            let (
                id,
                peptide_seq,
                peptide_mod_seq,
                precursor_mz,
                precursor_charge,
                retention_time,
                peak_mz_blob,
                peak_int_blob,
            ) = row_result
                .map_err(|e| OspreyError::LibraryLoadError(format!("Row read error: {}", e)))?;

            // Check if this is a new spectrum or continuation of previous
            if current_id != Some(id) {
                // Save previous entry if exists
                if let Some(entry) = current_entry.take() {
                    entries.push(entry);
                }

                // Parse modifications from modified sequence
                let modifications = parse_blib_modifications(&peptide_mod_seq);

                // Decode peaks
                let fragments =
                    if let (Some(mz_blob), Some(int_blob)) = (&peak_mz_blob, &peak_int_blob) {
                        decode_blib_peaks(mz_blob, int_blob)?
                    } else {
                        Vec::new()
                    };

                current_entry = Some(LibraryEntry {
                    id: id as u32,
                    sequence: peptide_seq,
                    modified_sequence: peptide_mod_seq,
                    modifications,
                    charge: precursor_charge as u8,
                    precursor_mz,
                    retention_time: retention_time.unwrap_or(0.0),
                    rt_calibrated: false,
                    fragments,
                    protein_ids: Vec::new(),
                    gene_names: Vec::new(),
                    is_decoy: false,
                });
                current_id = Some(id);
            }
        }

        // Don't forget the last entry
        if let Some(entry) = current_entry {
            entries.push(entry);
        }

        // Try to load protein mappings
        self.load_protein_mappings(conn, &mut entries)?;

        log::info!("Loaded {} library entries from blib", entries.len());
        Ok(entries)
    }

    /// Load protein mappings from RefSpectraProteins table
    fn load_protein_mappings(
        &self,
        conn: &Connection,
        entries: &mut [LibraryEntry],
    ) -> Result<()> {
        if !self.table_exists(conn, "RefSpectraProteins")? {
            return Ok(());
        }

        // Build a map of RefSpectraID -> protein accessions
        let mut stmt = conn
            .prepare(
                r#"
                SELECT rsp.RefSpectraID, p.accession
                FROM RefSpectraProteins rsp
                JOIN Proteins p ON rsp.ProteinID = p.id
                ORDER BY rsp.RefSpectraID
                "#,
            )
            .map_err(|e| {
                OspreyError::LibraryLoadError(format!("Failed to query proteins: {}", e))
            })?;

        let rows = stmt
            .query_map([], |row| {
                let ref_id: i64 = row.get(0)?;
                let accession: String = row.get(1)?;
                Ok((ref_id, accession))
            })
            .map_err(|e| OspreyError::LibraryLoadError(format!("Query error: {}", e)))?;

        // Create a lookup map
        let mut protein_map: std::collections::HashMap<u32, Vec<String>> =
            std::collections::HashMap::new();

        for row_result in rows {
            if let Ok((ref_id, accession)) = row_result {
                protein_map
                    .entry(ref_id as u32)
                    .or_default()
                    .push(accession);
            }
        }

        // Apply to entries
        for entry in entries.iter_mut() {
            if let Some(proteins) = protein_map.remove(&entry.id) {
                entry.protein_ids = proteins;
            }
        }

        Ok(())
    }

    /// Check if a table exists in the database
    fn table_exists(&self, conn: &Connection, table_name: &str) -> Result<bool> {
        let count: i32 = conn
            .query_row(
                "SELECT COUNT(*) FROM sqlite_master WHERE type='table' AND name=?",
                [table_name],
                |row| row.get(0),
            )
            .map_err(|e| OspreyError::LibraryLoadError(format!("Failed to check table: {}", e)))?;
        Ok(count > 0)
    }
}

impl Default for BlibLoader {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse modifications from BiblioSpec modified sequence format
/// BiblioSpec uses formats like:
/// - "PEPTC[+57.0]IDE" (mass shift in brackets)
/// - "PEPTC[160.0]IDE" (absolute mass)
fn parse_blib_modifications(mod_seq: &str) -> Vec<Modification> {
    let mut modifications = Vec::new();
    let mut chars = mod_seq.chars().peekable();
    let mut position: usize = 0;

    while let Some(c) = chars.next() {
        if c.is_ascii_alphabetic() {
            position += 1;
        } else if c == '[' {
            // Parse modification mass
            let mut mass_str = String::new();
            while let Some(&next) = chars.peek() {
                if next == ']' {
                    chars.next();
                    break;
                }
                mass_str.push(chars.next().unwrap());
            }

            // Parse the mass value
            if let Ok(mass) = mass_str.trim_start_matches('+').parse::<f64>() {
                let mod_position = if position > 0 { position - 1 } else { 0 };

                // Try to identify common modifications
                let (mass_delta, unimod_id, name) = identify_modification(mass, mod_position == 0);

                modifications.push(Modification {
                    position: mod_position,
                    unimod_id,
                    mass_delta,
                    name,
                });
            }
        }
    }

    modifications
}

/// Try to identify a modification by its mass
fn identify_modification(mass: f64, is_nterm: bool) -> (f64, Option<u32>, Option<String>) {
    // Common modifications with their UniMod IDs
    const CARBAMIDOMETHYL: f64 = 57.021464;
    const OXIDATION: f64 = 15.994915;
    const ACETYL: f64 = 42.010565;
    const PHOSPHO: f64 = 79.966331;
    const DEAMIDATION: f64 = 0.984016;
    const TMT6PLEX: f64 = 229.162932;

    let tolerance = 0.01;

    // Check if this is a mass shift (starts with + or -) or absolute mass
    // For absolute mass, we need to subtract the amino acid mass
    // For simplicity, we'll treat values > 100 as likely absolute masses for Cys (103.0)
    let mass_delta = if mass > 100.0 && mass < 200.0 {
        // Likely absolute mass on Cys
        mass - 103.009185 // Cys residue mass
    } else {
        mass
    };

    if (mass_delta - CARBAMIDOMETHYL).abs() < tolerance {
        (CARBAMIDOMETHYL, Some(4), Some("Carbamidomethyl".to_string()))
    } else if (mass_delta - OXIDATION).abs() < tolerance {
        (OXIDATION, Some(35), Some("Oxidation".to_string()))
    } else if (mass_delta - ACETYL).abs() < tolerance && is_nterm {
        (ACETYL, Some(1), Some("Acetyl".to_string()))
    } else if (mass_delta - PHOSPHO).abs() < tolerance {
        (PHOSPHO, Some(21), Some("Phospho".to_string()))
    } else if (mass_delta - DEAMIDATION).abs() < tolerance {
        (DEAMIDATION, Some(7), Some("Deamidated".to_string()))
    } else if (mass_delta - TMT6PLEX).abs() < tolerance {
        (TMT6PLEX, Some(737), Some("TMT6plex".to_string()))
    } else {
        (mass_delta, None, None)
    }
}

/// Decode peak blobs from blib format
/// BiblioSpec stores peaks as little-endian doubles
fn decode_blib_peaks(mz_blob: &[u8], intensity_blob: &[u8]) -> Result<Vec<LibraryFragment>> {
    // BiblioSpec typically stores as doubles (f64)
    if mz_blob.len() % 8 != 0 {
        return Err(OspreyError::LibraryLoadError(
            "Invalid peak m/z blob size".to_string(),
        ));
    }

    let n_peaks = mz_blob.len() / 8;

    // Intensity can be float (4 bytes) or double (8 bytes)
    let intensity_size = if intensity_blob.len() == n_peaks * 4 {
        4
    } else if intensity_blob.len() == n_peaks * 8 {
        8
    } else {
        return Err(OspreyError::LibraryLoadError(format!(
            "Invalid peak intensity blob size: expected {} or {} bytes, got {}",
            n_peaks * 4,
            n_peaks * 8,
            intensity_blob.len()
        )));
    };

    let mut fragments = Vec::with_capacity(n_peaks);

    for i in 0..n_peaks {
        let mz_bytes: [u8; 8] = mz_blob[i * 8..(i + 1) * 8].try_into().unwrap();
        let mz = f64::from_le_bytes(mz_bytes);

        let intensity = if intensity_size == 4 {
            let int_bytes: [u8; 4] = intensity_blob[i * 4..(i + 1) * 4].try_into().unwrap();
            f32::from_le_bytes(int_bytes)
        } else {
            let int_bytes: [u8; 8] = intensity_blob[i * 8..(i + 1) * 8].try_into().unwrap();
            f64::from_le_bytes(int_bytes) as f32
        };

        fragments.push(LibraryFragment {
            mz,
            relative_intensity: intensity,
            annotation: FragmentAnnotation {
                ion_type: IonType::Unknown,
                ordinal: (i + 1) as u8,
                charge: 1,
                neutral_loss: None,
            },
        });
    }

    // Normalize intensities to max = 1.0
    let max_intensity = fragments
        .iter()
        .map(|f| f.relative_intensity)
        .fold(0.0f32, |a, b| a.max(b));

    if max_intensity > 0.0 {
        for f in &mut fragments {
            f.relative_intensity /= max_intensity;
        }
    }

    Ok(fragments)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies that an unmodified peptide sequence produces an empty modifications list.
    #[test]
    fn test_parse_blib_modifications_simple() {
        let mods = parse_blib_modifications("PEPTIDE");
        assert!(mods.is_empty());
    }

    /// Verifies that a carbamidomethyl modification is parsed with correct position, mass, and UniMod ID.
    #[test]
    fn test_parse_blib_modifications_carbamidomethyl() {
        let mods = parse_blib_modifications("PEPTC[+57.021]IDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 4);
        assert!((mods[0].mass_delta - 57.021464).abs() < 0.01);
        assert_eq!(mods[0].unimod_id, Some(4));
    }

    /// Verifies that an oxidation modification is parsed with the correct residue position and UniMod ID.
    #[test]
    fn test_parse_blib_modifications_oxidation() {
        let mods = parse_blib_modifications("PEPTM[+15.995]IDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 4);
        assert_eq!(mods[0].unimod_id, Some(35));
    }

    /// Verifies that a mass delta is correctly identified as carbamidomethyl with its UniMod ID and name.
    #[test]
    fn test_identify_modification() {
        let (mass, unimod, name) = identify_modification(57.021, false);
        assert!((mass - 57.021464).abs() < 0.01);
        assert_eq!(unimod, Some(4));
        assert_eq!(name, Some("Carbamidomethyl".to_string()));
    }
}
