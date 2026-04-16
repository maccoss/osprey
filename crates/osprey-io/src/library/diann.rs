//! DIA-NN TSV library format parser
//!
//! Parses spectral libraries in the DIA-NN tab-separated format.
//! The format includes predicted fragment intensities and retention times.

use csv::ReaderBuilder;
use osprey_core::{
    FragmentAnnotation, IonType, LibraryEntry, LibraryFragment, LibraryLoader, Modification,
    NeutralLoss, OspreyError, Result,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Loader for DIA-NN TSV format libraries
pub struct DiannTsvLoader {
    /// Minimum number of fragments required per precursor
    min_fragments: usize,
}

impl DiannTsvLoader {
    /// Create a new DIA-NN TSV loader with default settings
    pub fn new() -> Self {
        Self { min_fragments: 3 }
    }

    /// Set minimum number of fragments required per precursor
    pub fn with_min_fragments(mut self, min: usize) -> Self {
        self.min_fragments = min;
        self
    }

    /// Load library from a path
    pub fn load<P: AsRef<Path>>(&self, path: P) -> Result<Vec<LibraryEntry>> {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| {
            OspreyError::LibraryLoadError(format!(
                "Failed to open library file '{}': {}",
                path.display(),
                e
            ))
        })?;

        let reader = BufReader::new(file);
        self.parse_reader(reader, path)
    }

    /// Parse library from a reader
    fn parse_reader<R: BufRead>(&self, reader: R, path: &Path) -> Result<Vec<LibraryEntry>> {
        let mut csv_reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_reader(reader);

        // Get headers and find column indices
        let headers = csv_reader.headers().map_err(|e| {
            OspreyError::LibraryLoadError(format!("Failed to read TSV headers: {}", e))
        })?;

        let col_indices = ColumnIndices::from_headers(headers)?;

        // Group fragments by precursor
        let mut precursor_map: HashMap<String, PrecursorData> = HashMap::new();

        for (row_num, result) in csv_reader.records().enumerate() {
            let record = result.map_err(|e| {
                OspreyError::LibraryLoadError(format!("Error at row {}: {}", row_num + 2, e))
            })?;

            // Parse row
            let row_data = self.parse_row(&record, &col_indices, row_num + 2)?;

            // Group by precursor key (modified sequence + charge)
            let key = format!("{}_{}", row_data.modified_sequence, row_data.charge);

            let precursor = precursor_map.entry(key).or_insert_with(|| PrecursorData {
                sequence: row_data.stripped_sequence.clone(),
                modified_sequence: row_data.modified_sequence.clone(),
                charge: row_data.charge,
                precursor_mz: row_data.precursor_mz,
                retention_time: row_data.retention_time,
                protein_ids: row_data.protein_ids.clone(),
                gene_names: row_data.gene_names.clone(),
                fragments: Vec::new(),
            });

            // Add fragment
            precursor.fragments.push(LibraryFragment {
                mz: row_data.fragment_mz,
                relative_intensity: row_data.relative_intensity,
                annotation: row_data.annotation,
            });
        }

        // Convert to LibraryEntry list
        let mut entries: Vec<LibraryEntry> = Vec::with_capacity(precursor_map.len());
        let mut id = 0u32;

        for (_key, data) in precursor_map {
            // Skip if too few fragments
            if data.fragments.len() < self.min_fragments {
                continue;
            }

            // Parse modifications from modified sequence
            let modifications = parse_modifications(&data.modified_sequence);

            let mut entry = LibraryEntry::new(
                id,
                data.sequence,
                data.modified_sequence,
                data.charge,
                data.precursor_mz,
                data.retention_time,
            );

            entry.modifications = modifications;
            entry.fragments = data.fragments;
            entry.protein_ids = data.protein_ids;
            entry.gene_names = data.gene_names;

            entries.push(entry);
            id += 1;
        }

        log::info!(
            "Loaded {} precursors from {}",
            entries.len(),
            path.display()
        );

        Ok(entries)
    }

    /// Parse a single row from the TSV
    fn parse_row(
        &self,
        record: &csv::StringRecord,
        cols: &ColumnIndices,
        row_num: usize,
    ) -> Result<RowData> {
        let get_field = |idx: Option<usize>, name: &str| -> Result<&str> {
            idx.and_then(|i| record.get(i)).ok_or_else(|| {
                OspreyError::LibraryLoadError(format!("Missing {} at row {}", name, row_num))
            })
        };

        let parse_f64 = |s: &str, name: &str| -> Result<f64> {
            s.parse().map_err(|_| {
                OspreyError::LibraryLoadError(format!(
                    "Invalid {} '{}' at row {}",
                    name, s, row_num
                ))
            })
        };

        let parse_f32 = |s: &str, name: &str| -> Result<f32> {
            s.parse().map_err(|_| {
                OspreyError::LibraryLoadError(format!(
                    "Invalid {} '{}' at row {}",
                    name, s, row_num
                ))
            })
        };

        let parse_u8 = |s: &str, name: &str| -> Result<u8> {
            s.parse().map_err(|_| {
                OspreyError::LibraryLoadError(format!(
                    "Invalid {} '{}' at row {}",
                    name, s, row_num
                ))
            })
        };

        // Required fields
        let precursor_mz = parse_f64(get_field(cols.precursor_mz, "PrecursorMz")?, "PrecursorMz")?;
        let charge = parse_u8(
            get_field(cols.precursor_charge, "PrecursorCharge")?,
            "PrecursorCharge",
        )?;
        let modified_sequence =
            strip_flanking_chars(get_field(cols.modified_peptide, "ModifiedPeptide")?);
        let fragment_mz = parse_f64(get_field(cols.fragment_mz, "FragmentMz")?, "FragmentMz")?;
        let relative_intensity = parse_f32(
            get_field(cols.relative_intensity, "RelativeIntensity")?,
            "RelativeIntensity",
        )?;

        // Get retention time (try multiple column names)
        let retention_time = if let Some(idx) = cols.irt {
            parse_f64(record.get(idx).unwrap_or("0"), "iRT")?
        } else if let Some(idx) = cols.normalized_rt {
            parse_f64(record.get(idx).unwrap_or("0"), "NormalizedRetentionTime")?
        } else {
            0.0
        };

        // Optional stripped sequence
        let stripped_sequence = cols
            .stripped_peptide
            .and_then(|i| record.get(i))
            .map(|s| s.to_string())
            .unwrap_or_else(|| strip_modifications(&modified_sequence));

        // Fragment annotation
        let ion_type = cols
            .fragment_type
            .and_then(|i| record.get(i))
            .and_then(|s| s.chars().next())
            .map(IonType::from_char)
            .unwrap_or(IonType::Unknown);

        let ordinal = cols
            .fragment_series_number
            .and_then(|i| record.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let fragment_charge = cols
            .fragment_charge
            .and_then(|i| record.get(i))
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);

        let neutral_loss = cols
            .fragment_loss_type
            .and_then(|i| record.get(i))
            .and_then(NeutralLoss::parse);

        let annotation = FragmentAnnotation {
            ion_type,
            ordinal,
            charge: fragment_charge,
            neutral_loss,
        };

        // Protein and gene info (may be semicolon or comma separated)
        let protein_ids = cols
            .protein_id
            .and_then(|i| record.get(i))
            .map(split_list)
            .unwrap_or_default();

        let gene_names = cols
            .gene_name
            .and_then(|i| record.get(i))
            .map(split_list)
            .unwrap_or_default();

        Ok(RowData {
            precursor_mz,
            charge,
            modified_sequence,
            stripped_sequence,
            retention_time,
            fragment_mz,
            relative_intensity,
            annotation,
            protein_ids,
            gene_names,
        })
    }
}

impl Default for DiannTsvLoader {
    fn default() -> Self {
        Self::new()
    }
}

impl LibraryLoader for DiannTsvLoader {
    fn load(&self, path: &Path) -> Result<Vec<LibraryEntry>> {
        DiannTsvLoader::load(self, path)
    }

    fn supports_format(&self, path: &Path) -> bool {
        path.extension()
            .is_some_and(|ext| ext == "tsv" || ext == "txt" || ext == "csv")
    }

    fn format_name(&self) -> &'static str {
        "DIA-NN TSV"
    }
}

/// Column indices for parsing
struct ColumnIndices {
    precursor_mz: Option<usize>,
    precursor_charge: Option<usize>,
    modified_peptide: Option<usize>,
    stripped_peptide: Option<usize>,
    fragment_mz: Option<usize>,
    relative_intensity: Option<usize>,
    fragment_type: Option<usize>,
    fragment_series_number: Option<usize>,
    fragment_charge: Option<usize>,
    fragment_loss_type: Option<usize>,
    irt: Option<usize>,
    normalized_rt: Option<usize>,
    protein_id: Option<usize>,
    gene_name: Option<usize>,
}

impl ColumnIndices {
    fn from_headers(headers: &csv::StringRecord) -> Result<Self> {
        let find = |names: &[&str]| -> Option<usize> {
            for name in names {
                for (i, h) in headers.iter().enumerate() {
                    if h.eq_ignore_ascii_case(name) {
                        return Some(i);
                    }
                }
            }
            None
        };

        let indices = Self {
            precursor_mz: find(&["PrecursorMz", "Precursor.Mz", "Q1"]),
            precursor_charge: find(&["PrecursorCharge", "Precursor.Charge"]),
            modified_peptide: find(&["ModifiedPeptide", "Modified.Peptide", "FullPeptideName"]),
            stripped_peptide: find(&["StrippedPeptide", "Stripped.Peptide", "PeptideSequence"]),
            fragment_mz: find(&["FragmentMz", "Fragment.Mz", "ProductMz", "Q3"]),
            relative_intensity: find(&[
                "RelativeIntensity",
                "Relative.Intensity",
                "LibraryIntensity",
            ]),
            fragment_type: find(&["FragmentType", "Fragment.Type", "IonType"]),
            fragment_series_number: find(&["FragmentSeriesNumber", "FragmentNumber", "IonNumber"]),
            fragment_charge: find(&["FragmentCharge", "Fragment.Charge", "ProductCharge"]),
            fragment_loss_type: find(&["FragmentLossType", "LossType", "NeutralLoss"]),
            irt: find(&["iRT", "iRt"]),
            normalized_rt: find(&["NormalizedRetentionTime", "Tr_recalibrated", "RT"]),
            protein_id: find(&["ProteinId", "Protein.Id", "ProteinName", "Protein"]),
            gene_name: find(&["GeneName", "Gene.Name", "Genes"]),
        };

        // Validate required columns exist
        if indices.precursor_mz.is_none() {
            return Err(OspreyError::LibraryLoadError(
                "Missing required column: PrecursorMz".to_string(),
            ));
        }
        if indices.precursor_charge.is_none() {
            return Err(OspreyError::LibraryLoadError(
                "Missing required column: PrecursorCharge".to_string(),
            ));
        }
        if indices.modified_peptide.is_none() {
            return Err(OspreyError::LibraryLoadError(
                "Missing required column: ModifiedPeptide".to_string(),
            ));
        }
        if indices.fragment_mz.is_none() {
            return Err(OspreyError::LibraryLoadError(
                "Missing required column: FragmentMz".to_string(),
            ));
        }
        if indices.relative_intensity.is_none() {
            return Err(OspreyError::LibraryLoadError(
                "Missing required column: RelativeIntensity".to_string(),
            ));
        }

        Ok(indices)
    }
}

/// Intermediate row data
struct RowData {
    precursor_mz: f64,
    charge: u8,
    modified_sequence: String,
    stripped_sequence: String,
    retention_time: f64,
    fragment_mz: f64,
    relative_intensity: f32,
    annotation: FragmentAnnotation,
    protein_ids: Vec<String>,
    gene_names: Vec<String>,
}

/// Intermediate precursor data for grouping fragments
struct PrecursorData {
    sequence: String,
    modified_sequence: String,
    charge: u8,
    precursor_mz: f64,
    retention_time: f64,
    protein_ids: Vec<String>,
    gene_names: Vec<String>,
    fragments: Vec<LibraryFragment>,
}

/// Strip flanking characters from peptide sequences
///
/// Handles various formats from different library sources:
/// - `_PEPTIDE_` → `PEPTIDE` (DIA-NN format)
/// - `K.PEPTIDE.R` → `PEPTIDE` (Percolator format)
/// - `-PEPTIDE-` → `PEPTIDE`
fn strip_flanking_chars(seq: &str) -> String {
    let trimmed = seq.trim_matches(|c| c == '_' || c == '.' || c == '-');

    // Also handle internal patterns like "K.PEPTIDE.R" -> "PEPTIDE"
    if let Some(start) = trimmed.find('.') {
        if let Some(end) = trimmed.rfind('.') {
            if start < end {
                // Extract the middle part between the first and last dots
                return trimmed[start + 1..end].to_string();
            }
        }
    }

    trimmed.to_string()
}

/// Strip modifications from a modified peptide sequence
fn strip_modifications(modified: &str) -> String {
    let mut result = String::with_capacity(modified.len());
    let mut in_bracket = false;
    let mut in_paren = false;

    for c in modified.chars() {
        match c {
            '[' => in_bracket = true,
            ']' => in_bracket = false,
            '(' => in_paren = true,
            ')' => in_paren = false,
            _ if !in_bracket && !in_paren && c.is_ascii_alphabetic() => {
                result.push(c);
            }
            _ => {}
        }
    }

    result
}

/// Parse modifications from modified sequence
fn parse_modifications(modified: &str) -> Vec<Modification> {
    let mut modifications = Vec::new();
    let mut position = 0usize;
    let mut chars = modified.chars().peekable();

    while let Some(c) = chars.next() {
        if c.is_ascii_alphabetic() {
            // Check for modification after this residue
            if chars.peek() == Some(&'[') {
                chars.next(); // consume '['
                let mut mod_str = String::new();
                while let Some(&next) = chars.peek() {
                    if next == ']' {
                        chars.next();
                        break;
                    }
                    mod_str.push(chars.next().unwrap());
                }

                // Parse modification
                if let Some(mass) = parse_mod_mass(&mod_str) {
                    modifications.push(Modification {
                        position,
                        unimod_id: parse_unimod_id(&mod_str),
                        mass_delta: mass,
                        name: Some(mod_str),
                    });
                }
            }
            position += 1;
        } else if c == '(' {
            // Handle parenthetical modifications (e.g., "(UniMod:35)")
            let mut mod_str = String::new();
            while let Some(&next) = chars.peek() {
                if next == ')' {
                    chars.next();
                    break;
                }
                mod_str.push(chars.next().unwrap());
            }

            if let Some(mass) = parse_mod_mass(&mod_str) {
                modifications.push(Modification {
                    position: position.saturating_sub(1),
                    unimod_id: parse_unimod_id(&mod_str),
                    mass_delta: mass,
                    name: Some(mod_str),
                });
            }
        }
    }

    modifications
}

/// Parse mass from modification string
fn parse_mod_mass(s: &str) -> Option<f64> {
    // Try parsing as a number directly
    if let Ok(mass) = s.parse::<f64>() {
        return Some(mass);
    }

    // Handle +/- prefix
    if s.starts_with('+') || s.starts_with('-') {
        if let Ok(mass) = s[1..].parse::<f64>() {
            return Some(if s.starts_with('-') { -mass } else { mass });
        }
    }

    // Try UniMod notation (e.g., "UniMod:4" or "UNIMOD:4")
    let unimod_prefix = s
        .strip_prefix("UniMod:")
        .or_else(|| s.strip_prefix("UNIMOD:"));
    if let Some(id_str) = unimod_prefix {
        if let Ok(id) = id_str.parse::<u32>() {
            if let Some(mass) = crate::output::unimod_id_to_mass(id) {
                return Some(mass);
            }
        }
    }

    // Known modifications by name
    match s.to_uppercase().as_str() {
        "OXIDATION" => Some(15.9949),
        "CARBAMIDOMETHYL" | "CAM" => Some(57.0215),
        "PHOSPHO" => Some(79.9663),
        "ACETYL" => Some(42.0106),
        "DEAMIDATED" | "DEAMIDATION" => Some(0.9840),
        _ => None,
    }
}

/// Parse UniMod ID from modification string
fn parse_unimod_id(s: &str) -> Option<u32> {
    if let Some(start) = s.find("UniMod:") {
        let rest = &s[start + 7..];
        let end = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        rest[..end].parse().ok()
    } else {
        None
    }
}

/// Split a semicolon or comma separated list
fn split_list(s: &str) -> Vec<String> {
    s.split([';', ','])
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies that modification notations in brackets and parentheses are correctly stripped from peptide sequences.
    #[test]
    fn test_strip_modifications() {
        assert_eq!(strip_modifications("PEPTIDE"), "PEPTIDE");
        assert_eq!(strip_modifications("PEP[+15.99]TIDE"), "PEPTIDE");
        assert_eq!(strip_modifications("(UniMod:1)PEPTIDE"), "PEPTIDE");
        assert_eq!(strip_modifications("PEPTM[Oxidation]IDE"), "PEPTMIDE");
    }

    /// Verifies that bracketed mass shift modifications are parsed with correct position and mass delta.
    #[test]
    fn test_parse_modifications() {
        let mods = parse_modifications("PEPTM[+15.9949]IDE");
        assert_eq!(mods.len(), 1);
        assert_eq!(mods[0].position, 4);
        assert!((mods[0].mass_delta - 15.9949).abs() < 1e-4);
    }

    /// Verifies that semicolon-delimited and comma-delimited protein/gene lists are correctly split and trimmed.
    #[test]
    fn test_split_list() {
        assert_eq!(split_list("A;B;C"), vec!["A", "B", "C"]);
        assert_eq!(split_list("A,B,C"), vec!["A", "B", "C"]);
        assert_eq!(split_list("A; B; C"), vec!["A", "B", "C"]);
    }

    /// Verifies that modification mass parsing handles signed numeric strings and named modifications like Oxidation.
    #[test]
    fn test_parse_mod_mass() {
        assert!((parse_mod_mass("+15.9949").unwrap() - 15.9949).abs() < 1e-4);
        assert!((parse_mod_mass("-17.026").unwrap() + 17.026).abs() < 1e-3);
        assert!((parse_mod_mass("Oxidation").unwrap() - 15.9949).abs() < 1e-4);
    }
}
