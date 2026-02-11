//! EncyclopeDIA elib format loader
//!
//! Parses spectral libraries in EncyclopeDIA's SQLite-based elib format.
//! The elib format stores peptide precursors with their fragment spectra
//! in a SQLite database.

use osprey_core::{
    FragmentAnnotation, IonType, LibraryEntry, LibraryFragment, Modification, OspreyError, Result,
};
use rusqlite::{Connection, OpenFlags};
use std::path::Path;

/// EncyclopeDIA elib format loader
pub struct ElibLoader;

impl ElibLoader {
    pub fn new() -> Self {
        Self
    }

    /// Load library entries from an elib file
    pub fn load<P: AsRef<Path>>(&self, path: P) -> Result<Vec<LibraryEntry>> {
        let path = path.as_ref();

        // Open in read-only mode
        let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY)
            .map_err(|e| OspreyError::LibraryLoadError(format!("Failed to open elib: {}", e)))?;

        // First, check which schema version we have
        let has_entries_table = self.table_exists(&conn, "entries")?;
        let has_peptide_to_protein = self.table_exists(&conn, "peptidetoprotein")?;

        if has_entries_table {
            self.load_standard_elib(&conn)
        } else if has_peptide_to_protein {
            self.load_chromatogram_elib(&conn)
        } else {
            Err(OspreyError::LibraryLoadError(
                "Unknown elib schema - no recognized tables found".to_string(),
            ))
        }
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

    /// Load from standard elib format (EncyclopeDIA library format)
    fn load_standard_elib(&self, conn: &Connection) -> Result<Vec<LibraryEntry>> {
        // Query all entries with their peaks
        let mut stmt = conn
            .prepare(
                r#"
                SELECT
                    e.PeptideModSeq,
                    e.PrecursorCharge,
                    e.PrecursorMz,
                    e.RTInSeconds,
                    e.RTInSecondsStart,
                    e.RTInSecondsStop,
                    e.MedianChromatogramEncodedLength,
                    e.ProteinAccession,
                    p.MassEncodedLength,
                    p.MassArray,
                    p.IntensityEncodedLength,
                    p.IntensityArray
                FROM entries e
                LEFT JOIN peaks p ON e.PeptideModSeq = p.PeptideModSeq
                    AND e.PrecursorCharge = p.PrecursorCharge
                ORDER BY e.PeptideModSeq, e.PrecursorCharge
                "#,
            )
            .map_err(|e| OspreyError::LibraryLoadError(format!("SQL prepare error: {}", e)))?;

        let mut entries: Vec<LibraryEntry> = Vec::new();
        let mut lib_id = 0u32;

        let rows = stmt
            .query_map([], |row| {
                let peptide_mod_seq: String = row.get(0)?;
                let precursor_charge: i32 = row.get(1)?;
                let precursor_mz: f64 = row.get(2)?;
                let rt_seconds: Option<f64> = row.get(3)?;
                let protein: Option<String> = row.get(7)?;
                let mass_array: Option<Vec<u8>> = row.get(9)?;
                let intensity_array: Option<Vec<u8>> = row.get(11)?;

                Ok((
                    peptide_mod_seq,
                    precursor_charge,
                    precursor_mz,
                    rt_seconds,
                    protein,
                    mass_array,
                    intensity_array,
                ))
            })
            .map_err(|e| OspreyError::LibraryLoadError(format!("Query error: {}", e)))?;

        for row_result in rows {
            let (
                peptide_mod_seq,
                precursor_charge,
                precursor_mz,
                rt_seconds,
                protein,
                mass_array,
                intensity_array,
            ) = row_result
                .map_err(|e| OspreyError::LibraryLoadError(format!("Row read error: {}", e)))?;

            // Parse sequence and modifications
            let (sequence, modifications) = parse_modified_sequence(&peptide_mod_seq)?;

            // Decode fragment peaks
            let fragments = if let (Some(mass_blob), Some(int_blob)) = (mass_array, intensity_array)
            {
                decode_peaks(&mass_blob, &int_blob)?
            } else {
                Vec::new()
            };

            // Convert RT from seconds to minutes
            let rt_minutes = rt_seconds.map(|rt| rt / 60.0).unwrap_or(0.0);

            let entry = LibraryEntry {
                id: lib_id,
                sequence: sequence.clone(),
                modified_sequence: peptide_mod_seq,
                modifications,
                charge: precursor_charge as u8,
                precursor_mz,
                retention_time: rt_minutes,
                rt_calibrated: false,
                fragments,
                protein_ids: protein.map(|p| vec![p]).unwrap_or_default(),
                gene_names: Vec::new(),
                is_decoy: false,
            };

            entries.push(entry);
            lib_id += 1;
        }

        log::info!(
            "Loaded {} library entries from elib (standard format)",
            entries.len()
        );
        Ok(entries)
    }

    /// Load from chromatogram library elib format (Skyline-exported or similar)
    fn load_chromatogram_elib(&self, conn: &Connection) -> Result<Vec<LibraryEntry>> {
        // This format has peptidetoprotein and separate spectrum tables
        let mut stmt = conn
            .prepare(
                r#"
                SELECT DISTINCT
                    PeptideModSeq,
                    PrecursorCharge,
                    PrecursorMz
                FROM peptidetoprotein
                ORDER BY PeptideModSeq, PrecursorCharge
                "#,
            )
            .map_err(|e| OspreyError::LibraryLoadError(format!("SQL prepare error: {}", e)))?;

        let mut entries: Vec<LibraryEntry> = Vec::new();
        let mut lib_id = 0u32;

        let rows = stmt
            .query_map([], |row| {
                let peptide_mod_seq: String = row.get(0)?;
                let precursor_charge: i32 = row.get(1)?;
                let precursor_mz: f64 = row.get(2)?;
                Ok((peptide_mod_seq, precursor_charge, precursor_mz))
            })
            .map_err(|e| OspreyError::LibraryLoadError(format!("Query error: {}", e)))?;

        for row_result in rows {
            let (peptide_mod_seq, precursor_charge, precursor_mz) = row_result
                .map_err(|e| OspreyError::LibraryLoadError(format!("Row read error: {}", e)))?;

            let (sequence, modifications) = parse_modified_sequence(&peptide_mod_seq)?;

            let entry = LibraryEntry {
                id: lib_id,
                sequence,
                modified_sequence: peptide_mod_seq,
                modifications,
                charge: precursor_charge as u8,
                precursor_mz,
                retention_time: 0.0,
                rt_calibrated: false,
                fragments: Vec::new(), // Will need separate query for peaks
                protein_ids: Vec::new(),
                gene_names: Vec::new(),
                is_decoy: false,
            };

            entries.push(entry);
            lib_id += 1;
        }

        log::info!(
            "Loaded {} library entries from elib (chromatogram format)",
            entries.len()
        );
        Ok(entries)
    }
}

impl Default for ElibLoader {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a modified sequence string (e.g., "PEPTC[+57.021]IDE") into
/// sequence and modifications
fn parse_modified_sequence(mod_seq: &str) -> Result<(String, Vec<Modification>)> {
    let mut sequence = String::new();
    let mut modifications = Vec::new();
    let mut chars = mod_seq.chars().peekable();
    let mut position: usize = 0;

    while let Some(c) = chars.next() {
        if c.is_ascii_alphabetic() {
            sequence.push(c.to_ascii_uppercase());
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

            // Parse the mass value (may have + or - prefix)
            if let Ok(mass) = mass_str.trim_start_matches('+').parse::<f64>() {
                // Modification on previous residue (position is 0-indexed)
                let mod_position = if position > 0 { position - 1 } else { 0 };
                modifications.push(Modification {
                    position: mod_position,
                    unimod_id: None,
                    mass_delta: mass,
                    name: None,
                });
            }
        } else if c == '(' {
            // Some formats use parentheses for modifications
            let mut mass_str = String::new();
            while let Some(&next) = chars.peek() {
                if next == ')' {
                    chars.next();
                    break;
                }
                mass_str.push(chars.next().unwrap());
            }

            if let Ok(mass) = mass_str.trim_start_matches('+').parse::<f64>() {
                let mod_position = if position > 0 { position - 1 } else { 0 };
                modifications.push(Modification {
                    position: mod_position,
                    unimod_id: None,
                    mass_delta: mass,
                    name: None,
                });
            }
        }
        // Ignore other characters (like underscores, periods for flanking residues)
    }

    Ok((sequence, modifications))
}

/// Decode peaks from elib blob format
/// Elib stores peaks as compressed double arrays
fn decode_peaks(mass_blob: &[u8], intensity_blob: &[u8]) -> Result<Vec<LibraryFragment>> {
    // Elib stores masses and intensities as arrays of doubles (f64)
    // Check if the blobs are the right size for doubles
    if mass_blob.len() % 8 != 0 || intensity_blob.len() % 8 != 0 {
        // Try as floats (f32)
        if mass_blob.len() % 4 == 0 && intensity_blob.len() % 4 == 0 {
            return decode_peaks_f32(mass_blob, intensity_blob);
        }
        return Err(OspreyError::LibraryLoadError(
            "Invalid peak blob size".to_string(),
        ));
    }

    let n_masses = mass_blob.len() / 8;
    let n_intensities = intensity_blob.len() / 8;

    if n_masses != n_intensities {
        return Err(OspreyError::LibraryLoadError(format!(
            "Mismatched peak arrays: {} masses, {} intensities",
            n_masses, n_intensities
        )));
    }

    let mut fragments = Vec::with_capacity(n_masses);

    for i in 0..n_masses {
        let mass_bytes: [u8; 8] = mass_blob[i * 8..(i + 1) * 8].try_into().unwrap();
        let int_bytes: [u8; 8] = intensity_blob[i * 8..(i + 1) * 8].try_into().unwrap();

        let mz = f64::from_le_bytes(mass_bytes);
        let intensity = f64::from_le_bytes(int_bytes) as f32;

        // Create a generic fragment (we don't have ion type info in basic elib)
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

    // Normalize intensities
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

/// Decode peaks stored as f32 arrays
fn decode_peaks_f32(mass_blob: &[u8], intensity_blob: &[u8]) -> Result<Vec<LibraryFragment>> {
    let n_masses = mass_blob.len() / 4;
    let n_intensities = intensity_blob.len() / 4;

    if n_masses != n_intensities {
        return Err(OspreyError::LibraryLoadError(format!(
            "Mismatched peak arrays: {} masses, {} intensities",
            n_masses, n_intensities
        )));
    }

    let mut fragments = Vec::with_capacity(n_masses);

    for i in 0..n_masses {
        let mass_bytes: [u8; 4] = mass_blob[i * 4..(i + 1) * 4].try_into().unwrap();
        let int_bytes: [u8; 4] = intensity_blob[i * 4..(i + 1) * 4].try_into().unwrap();

        let mz = f32::from_le_bytes(mass_bytes) as f64;
        let intensity = f32::from_le_bytes(int_bytes);

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

    // Normalize intensities
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

    /// Verifies that an unmodified sequence is parsed into the correct plain sequence with no modifications.
    #[test]
    fn test_parse_modified_sequence_simple() {
        let (seq, mods) = parse_modified_sequence("PEPTIDE").unwrap();
        assert_eq!(seq, "PEPTIDE");
        assert!(mods.is_empty());
    }

    /// Verifies that a bracketed mass shift modification is extracted with the correct residue position and mass.
    #[test]
    fn test_parse_modified_sequence_with_mod() {
        let (seq, mods) = parse_modified_sequence("PEPTC[+57.021]IDE").unwrap();
        assert_eq!(seq, "PEPTCIDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 4); // 0-indexed, C is at position 4
        assert!((mods[0].mass_delta - 57.021).abs() < 0.001);
    }

    /// Verifies that an N-terminal modification is parsed and assigned to position 0.
    #[test]
    fn test_parse_modified_sequence_nterm() {
        let (seq, mods) = parse_modified_sequence("[+42.011]PEPTIDE").unwrap();
        assert_eq!(seq, "PEPTIDE");
        assert_eq!(mods.len(), 1);
        // N-terminal mod is stored at position 0 (the first residue)
        assert_eq!(mods[0].position, 0);
    }

    /// Verifies that multiple modifications on different residues are all correctly parsed.
    #[test]
    fn test_parse_modified_sequence_multiple() {
        let (seq, mods) = parse_modified_sequence("PEP[+15.995]TC[+57.021]IDE").unwrap();
        assert_eq!(seq, "PEPTCIDE");
        assert_eq!(mods.len(), 2);
    }
}
