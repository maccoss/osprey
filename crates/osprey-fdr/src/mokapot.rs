//! Mokapot integration for semi-supervised FDR control
//!
//! Mokapot is a Python tool for semi-supervised learning in proteomics.
//! This module provides:
//! - PIN file generation
//! - Mokapot execution
//! - Result parsing

use osprey_core::{FeatureSet, OspreyError, Result};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::process::{Command, Stdio};

/// PSM (Peptide-Spectrum Match) with features for mokapot
#[derive(Debug, Clone)]
pub struct PsmFeatures {
    /// Unique PSM identifier
    pub psm_id: String,
    /// Peptide sequence (modified)
    pub peptide: String,
    /// Protein accessions
    pub proteins: Vec<String>,
    /// Scan number
    pub scan_number: u32,
    /// Run/file name
    pub file_name: String,
    /// Precursor charge
    pub charge: u8,
    /// Is this a decoy?
    pub is_decoy: bool,
    /// Feature set
    pub features: FeatureSet,
    /// Optional pre-computed score (for initial ranking)
    pub initial_score: Option<f64>,
}

/// Mokapot results for a PSM
#[derive(Debug, Clone)]
pub struct MokapotResult {
    /// PSM identifier
    pub psm_id: String,
    /// Mokapot score
    pub score: f64,
    /// Q-value
    pub q_value: f64,
    /// Posterior error probability
    pub pep: f64,
}

/// Mokapot integration
pub struct MokapotRunner {
    /// Path to mokapot executable (or "mokapot" if in PATH)
    mokapot_path: String,
    /// Training FDR threshold
    train_fdr: f64,
    /// Test FDR threshold
    test_fdr: f64,
    /// Maximum iterations for training
    max_iter: u32,
    /// Random seed
    seed: u64,
}

impl MokapotRunner {
    /// Create a new mokapot runner with default settings
    pub fn new() -> Self {
        Self {
            mokapot_path: "mokapot".to_string(),
            train_fdr: 0.01,
            test_fdr: 0.01,
            max_iter: 10,
            seed: 42,
        }
    }

    /// Set the mokapot executable path
    pub fn with_path(mut self, path: &str) -> Self {
        self.mokapot_path = path.to_string();
        self
    }

    /// Set training FDR threshold
    pub fn with_train_fdr(mut self, fdr: f64) -> Self {
        self.train_fdr = fdr;
        self
    }

    /// Set test FDR threshold
    pub fn with_test_fdr(mut self, fdr: f64) -> Self {
        self.test_fdr = fdr;
        self
    }

    /// Check if mokapot is available
    pub fn is_available(&self) -> bool {
        Command::new(&self.mokapot_path)
            .arg("--help")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok()
    }

    /// Write PSMs to a PIN file
    pub fn write_pin<P: AsRef<Path>>(&self, psms: &[PsmFeatures], path: P) -> Result<()> {
        let mut file = std::fs::File::create(path.as_ref()).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create PIN file: {}", e))
        })?;

        // Write header
        writeln!(
            file,
            "SpecId\tLabel\tScanNr\tChargeState\t{}\tPeptide\tProteins",
            get_feature_header()
        )
        .map_err(|e| OspreyError::OutputError(format!("Failed to write PIN header: {}", e)))?;

        // Write PSMs
        for psm in psms {
            let label = if psm.is_decoy { -1 } else { 1 };
            let proteins = if psm.proteins.is_empty() {
                "UNKNOWN".to_string()
            } else {
                psm.proteins.join(",")
            };

            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                psm.psm_id,
                label,
                psm.scan_number,
                psm.charge,
                format_features(&psm.features),
                format_peptide(&psm.peptide),
                proteins
            )
            .map_err(|e| OspreyError::OutputError(format!("Failed to write PSM: {}", e)))?;
        }

        Ok(())
    }

    /// Run mokapot on a PIN file
    pub fn run<P: AsRef<Path>>(
        &self,
        pin_file: P,
        output_dir: P,
    ) -> Result<Vec<MokapotResult>> {
        let pin_path = pin_file.as_ref();
        let out_path = output_dir.as_ref();

        // Create output directory if needed
        std::fs::create_dir_all(out_path).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create output directory: {}", e))
        })?;

        // Build mokapot command
        let output = Command::new(&self.mokapot_path)
            .arg(pin_path)
            .arg("--dest_dir")
            .arg(out_path)
            .arg("--train_fdr")
            .arg(self.train_fdr.to_string())
            .arg("--test_fdr")
            .arg(self.test_fdr.to_string())
            .arg("--max_iter")
            .arg(self.max_iter.to_string())
            .arg("--seed")
            .arg(self.seed.to_string())
            .output()
            .map_err(|e| {
                OspreyError::ExternalToolError(format!(
                    "Failed to run mokapot: {}. Is mokapot installed? Try: pip install mokapot",
                    e
                ))
            })?;

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(OspreyError::ExternalToolError(format!(
                "Mokapot failed: {}",
                stderr
            )));
        }

        // Parse results
        let results_file = out_path.join("mokapot.psms.txt");
        self.parse_results(&results_file)
    }

    /// Parse mokapot results file
    pub fn parse_results<P: AsRef<Path>>(&self, path: P) -> Result<Vec<MokapotResult>> {
        let file = std::fs::File::open(path.as_ref()).map_err(|e| {
            OspreyError::OutputError(format!("Failed to open mokapot results: {}", e))
        })?;

        let reader = BufReader::new(file);
        let mut results = Vec::new();
        let mut header_parsed = false;
        let mut psm_id_idx = 0;
        let mut score_idx = 0;
        let mut qvalue_idx = 0;
        let mut pep_idx = 0;

        for line_result in reader.lines() {
            let line = line_result.map_err(|e| {
                OspreyError::OutputError(format!("Failed to read line: {}", e))
            })?;

            let fields: Vec<&str> = line.split('\t').collect();

            if !header_parsed {
                // Parse header to find column indices
                for (i, field) in fields.iter().enumerate() {
                    match field.to_lowercase().as_str() {
                        "specid" | "psm_id" | "psmid" => psm_id_idx = i,
                        "score" | "mokapot_score" | "mokapot score" => score_idx = i,
                        "q-value" | "qvalue" | "mokapot_qvalue" | "mokapot q-value" => qvalue_idx = i,
                        "posterior_error_prob" | "pep" | "mokapot_pep" => pep_idx = i,
                        _ => {}
                    }
                }
                header_parsed = true;
                continue;
            }

            if fields.len() <= psm_id_idx.max(score_idx).max(qvalue_idx).max(pep_idx) {
                continue;
            }

            let psm_id = fields[psm_id_idx].to_string();
            let score = fields[score_idx].parse::<f64>().unwrap_or(0.0);
            let q_value = fields[qvalue_idx].parse::<f64>().unwrap_or(1.0);
            let pep = fields.get(pep_idx).and_then(|s| s.parse::<f64>().ok()).unwrap_or(1.0);

            results.push(MokapotResult {
                psm_id,
                score,
                q_value,
                pep,
            });
        }

        log::info!("Parsed {} mokapot results", results.len());
        Ok(results)
    }

    /// Run mokapot using the Python API directly (alternative method)
    /// This is useful when mokapot CLI is not available
    pub fn run_python_api<P: AsRef<Path>>(
        &self,
        pin_file: P,
        output_dir: P,
    ) -> Result<Vec<MokapotResult>> {
        let pin_path = pin_file.as_ref();
        let out_path = output_dir.as_ref();

        std::fs::create_dir_all(out_path).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create output directory: {}", e))
        })?;

        // Python script to run mokapot
        let python_script = format!(
            r#"
import mokapot
import pandas as pd

# Read PIN file
psms = mokapot.read_pin("{pin_file}")

# Configure and run mokapot
results, models = mokapot.brew(
    psms,
    train_fdr={train_fdr},
    test_fdr={test_fdr},
    max_iter={max_iter},
    rng={seed},
)

# Write results
results.to_txt("{output_dir}")
print("SUCCESS")
"#,
            pin_file = pin_path.display(),
            train_fdr = self.train_fdr,
            test_fdr = self.test_fdr,
            max_iter = self.max_iter,
            seed = self.seed,
            output_dir = out_path.display(),
        );

        let output = Command::new("python3")
            .arg("-c")
            .arg(&python_script)
            .output()
            .map_err(|e| {
                OspreyError::ExternalToolError(format!("Failed to run Python: {}", e))
            })?;

        let stdout = String::from_utf8_lossy(&output.stdout);
        if !stdout.contains("SUCCESS") {
            let stderr = String::from_utf8_lossy(&output.stderr);
            return Err(OspreyError::ExternalToolError(format!(
                "Mokapot Python API failed: {}",
                stderr
            )));
        }

        // Parse results
        let results_file = out_path.join("mokapot.psms.txt");
        self.parse_results(&results_file)
    }
}

impl Default for MokapotRunner {
    fn default() -> Self {
        Self::new()
    }
}

/// Get the feature header for PIN format
fn get_feature_header() -> String {
    [
        "peak_apex",
        "peak_area",
        "emg_fit_quality",
        "peak_width",
        "peak_symmetry",
        "rt_deviation",
        "rt_deviation_normalized",
        "n_contributing_scans",
        "coefficient_stability",
        "peak_sharpness",
        "peak_prominence",
        "hyperscore",
        "spectral_contrast_angle",
        "dot_product",
        "pearson_correlation",
        "spearman_correlation",
        "fragment_coverage",
        "sequence_coverage",
        "consecutive_ions",
        "base_peak_rank",
        "top3_matches",
        "explained_intensity",
        "n_competitors",
        "relative_coefficient",
        "local_peptide_density",
        "spectral_complexity",
        "regression_residual",
        "modification_count",
    ]
    .join("\t")
}

/// Format features for PIN output
fn format_features(features: &FeatureSet) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        features.peak_apex,
        features.peak_area,
        features.emg_fit_quality,
        features.peak_width,
        features.peak_symmetry,
        features.rt_deviation,
        features.rt_deviation_normalized,
        features.n_contributing_scans,
        features.coefficient_stability,
        features.peak_sharpness,
        features.peak_prominence,
        features.hyperscore,
        features.spectral_contrast_angle,
        features.dot_product,
        features.pearson_correlation,
        features.spearman_correlation,
        features.fragment_coverage,
        features.sequence_coverage,
        features.consecutive_ions,
        features.base_peak_rank,
        features.top3_matches,
        features.explained_intensity,
        features.n_competitors,
        features.relative_coefficient,
        features.local_peptide_density,
        features.spectral_complexity,
        features.regression_residual,
        features.modification_count,
    )
}

/// Format peptide for PIN output
/// Percolator/mokapot expects format: FLANKING.SEQUENCE.FLANKING
fn format_peptide(peptide: &str) -> String {
    // If peptide already has flanking residues, use as-is
    if peptide.contains('.') {
        peptide.to_string()
    } else {
        // Add placeholder flanking residues
        format!("-.{}.-", peptide)
    }
}

/// Write a simple TSV report of mokapot results
pub fn write_mokapot_report<P: AsRef<Path>>(
    results: &[MokapotResult],
    path: P,
) -> Result<()> {
    let mut file = std::fs::File::create(path.as_ref()).map_err(|e| {
        OspreyError::OutputError(format!("Failed to create report file: {}", e))
    })?;

    writeln!(file, "PSM_ID\tMokapot_Score\tQ_Value\tPEP").map_err(|e| {
        OspreyError::OutputError(format!("Failed to write header: {}", e))
    })?;

    for r in results {
        writeln!(
            file,
            "{}\t{:.6}\t{:.6}\t{:.6}",
            r.psm_id, r.score, r.q_value, r.pep
        )
        .map_err(|e| OspreyError::OutputError(format!("Failed to write result: {}", e)))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_peptide() {
        assert_eq!(format_peptide("PEPTIDE"), "-.PEPTIDE.-");
        assert_eq!(format_peptide("K.PEPTIDE.R"), "K.PEPTIDE.R");
    }

    #[test]
    fn test_feature_header() {
        let header = get_feature_header();
        assert!(header.contains("peak_apex"));
        assert!(header.contains("dot_product"));
        assert!(!header.contains("precursor_intensity")); // Optional field not included
    }

    #[test]
    fn test_mokapot_runner_availability() {
        let runner = MokapotRunner::new();
        // This will be false on systems without mokapot, which is expected
        let _available = runner.is_available();
    }
}
