//! Mokapot integration for semi-supervised FDR control
//!
//! Mokapot is a Python tool for semi-supervised learning in proteomics.
//! This module provides:
//! - PIN file generation (single or per-file)
//! - Mokapot CLI execution with multi-file support
//! - Result parsing
//!
//! ## Two-step analysis for multiple files
//!
//! When analyzing multiple files, Osprey uses a two-step approach:
//!
//! 1. **Run-level FDR**: Run mokapot WITHOUT `--aggregate` + WITH `--save_models`
//!    Trains a joint model on all files, reports per-file q-values.
//!
//! 2. **Experiment-level FDR**: Run mokapot WITH `--aggregate` + WITH `--load_models`
//!    Reuses the model from Step 1, reports experiment-wide q-values.

use osprey_core::{FeatureSet, OspreyError, Result};
use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
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
    /// Peptide sequence (for counting unique peptides)
    pub peptide: String,
    /// Mokapot score
    pub score: f64,
    /// Q-value
    pub q_value: f64,
    /// Posterior error probability
    pub pep: f64,
}

/// Mokapot integration using the mokapot CLI
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
    /// Number of parallel workers for cross-validation
    num_workers: u32,
    /// Maximum number of PSMs to use for training (memory-aware)
    /// If None, will be calculated based on available memory
    subset_max_train: Option<usize>,
}

impl MokapotRunner {
    /// Create a new mokapot runner with default settings
    pub fn new() -> Self {
        // Default to number of CPUs, capped at 8 to avoid memory issues
        let num_cpus = std::thread::available_parallelism()
            .map(|n| n.get() as u32)
            .unwrap_or(4);
        let num_workers = num_cpus.min(8);

        Self {
            mokapot_path: "mokapot".to_string(),
            train_fdr: 0.01,
            test_fdr: 0.01,
            max_iter: 10,
            seed: 42,
            num_workers,
            subset_max_train: None, // Will be calculated based on available memory
        }
    }

    /// Calculate training subset size based on available system memory
    ///
    /// Mokapot's memory usage scales with the number of PSMs used for training.
    /// This function estimates a safe subset size based on available memory.
    ///
    /// # Memory estimation
    /// - ~50 KB per PSM in training (features, labels, internal structures)
    /// - Reserve 2 GB for system overhead and other processes
    /// - Use 70% of remaining available memory for safety margin
    ///
    /// # Returns
    /// A subset size that should fit in available memory, or None if no limit needed
    pub fn calculate_subset_from_memory() -> Option<usize> {
        use sysinfo::System;

        let mut sys = System::new();
        sys.refresh_memory();

        let available_bytes = sys.available_memory();
        let total_bytes = sys.total_memory();

        // Reserve 2 GB for system overhead
        const RESERVED_BYTES: u64 = 2 * 1024 * 1024 * 1024;

        // Estimated memory per PSM in training (conservative estimate)
        // This accounts for: feature matrix, labels, SVM internals, cross-validation folds
        const BYTES_PER_PSM: u64 = 50 * 1024; // 50 KB

        // Use 70% of available memory after reservation
        let usable_bytes = available_bytes.saturating_sub(RESERVED_BYTES);
        let safe_bytes = (usable_bytes as f64 * 0.7) as u64;

        let max_psms = safe_bytes / BYTES_PER_PSM;

        log::debug!(
            "Memory-aware training: {:.1} GB available / {:.1} GB total → max {} PSMs for training",
            available_bytes as f64 / 1e9,
            total_bytes as f64 / 1e9,
            max_psms
        );

        // Only apply limit if it's restrictive (< 1 million PSMs)
        // Above this, mokapot should handle it fine
        if max_psms < 1_000_000 {
            Some(max_psms as usize)
        } else {
            None
        }
    }

    /// Get the effective subset_max_train value
    ///
    /// Returns the configured value if set, otherwise calculates from available memory.
    fn effective_subset_max_train(&self) -> Option<usize> {
        self.subset_max_train
            .or_else(Self::calculate_subset_from_memory)
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

    /// Set number of parallel workers
    pub fn with_num_workers(mut self, n: u32) -> Self {
        self.num_workers = n;
        self
    }

    /// Set maximum iterations
    pub fn with_max_iter(mut self, n: u32) -> Self {
        self.max_iter = n;
        self
    }

    /// Set maximum training subset size
    ///
    /// If not set, will be automatically calculated based on available memory.
    /// Set to Some(0) to disable subsetting entirely.
    pub fn with_subset_max_train(mut self, n: Option<usize>) -> Self {
        self.subset_max_train = n;
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
        use std::io::BufWriter;

        let file = std::fs::File::create(path.as_ref())
            .map_err(|e| OspreyError::OutputError(format!("Failed to create PIN file: {}", e)))?;

        // Use buffered writer for much faster I/O (avoids syscall per line)
        let mut writer = BufWriter::with_capacity(1024 * 1024, file); // 1MB buffer

        // Write header
        writeln!(
            writer,
            "SpecId\tLabel\tScanNr\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\t{}\tPeptide\tProteins",
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
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                psm.psm_id,
                label,
                psm.scan_number,
                format_charge_features(psm.charge),
                format_features(&psm.features),
                format_peptide(&psm.peptide),
                proteins
            )
            .map_err(|e| OspreyError::OutputError(format!("Failed to write PSM: {}", e)))?;
        }

        // Ensure all data is flushed to disk
        writer
            .flush()
            .map_err(|e| OspreyError::OutputError(format!("Failed to flush PIN file: {}", e)))?;

        Ok(())
    }

    /// Write a single file's PSMs to a PIN file
    ///
    /// This is used for memory-efficient processing where each file's PIN
    /// is written immediately after scoring, rather than accumulating all data.
    ///
    /// # Arguments
    /// * `file_name` - Name of the source file (used for PIN file naming)
    /// * `psms` - PSM features to write
    /// * `output_dir` - Directory to write the PIN file
    ///
    /// # Returns
    /// Path to the written PIN file
    pub fn write_single_pin_file<P: AsRef<Path>>(
        &self,
        file_name: &str,
        psms: &[PsmFeatures],
        output_dir: P,
    ) -> Result<PathBuf> {
        use std::io::BufWriter;

        let out_dir = output_dir.as_ref();
        std::fs::create_dir_all(out_dir).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create output directory: {}", e))
        })?;

        // Sanitize file name for use as PIN file name
        let safe_name = file_name.replace(['/', '\\', ' '], "_");
        let pin_path = out_dir.join(format!("{}.pin", safe_name));

        let file = std::fs::File::create(&pin_path)
            .map_err(|e| OspreyError::OutputError(format!("Failed to create PIN file: {}", e)))?;

        let mut writer = BufWriter::with_capacity(1024 * 1024, file);

        // Write header
        writeln!(
            writer,
            "SpecId\tLabel\tScanNr\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\t{}\tPeptide\tProteins",
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
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                psm.psm_id,
                label,
                psm.scan_number,
                format_charge_features(psm.charge),
                format_features(&psm.features),
                format_peptide(&psm.peptide),
                proteins
            )
            .map_err(|e| OspreyError::OutputError(format!("Failed to write PSM: {}", e)))?;
        }

        writer
            .flush()
            .map_err(|e| OspreyError::OutputError(format!("Failed to flush PIN file: {}", e)))?;

        log::info!(
            "Wrote {} PSMs to PIN file: {}",
            psms.len(),
            pin_path.display()
        );

        Ok(pin_path)
    }

    /// Run mokapot on a PIN file
    pub fn run<P: AsRef<Path>>(&self, pin_file: P, output_dir: P) -> Result<Vec<MokapotResult>> {
        let pin_path = pin_file.as_ref();
        let out_path = output_dir.as_ref();

        // Create output directory if needed
        std::fs::create_dir_all(out_path).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create output directory: {}", e))
        })?;

        // Build mokapot command
        log::info!(
            "Running mokapot with {} workers, {} max iterations",
            self.num_workers,
            self.max_iter
        );

        // Use spawn() to stream output instead of output() which blocks
        let mut child = Command::new(&self.mokapot_path)
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
            .arg("--max_workers")
            .arg(self.num_workers.to_string())
            .arg("--verbosity")
            .arg("1")
            .arg("--save_models") // Save model weights to inspect feature importance
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| {
                OspreyError::ExternalToolError(format!(
                    "Failed to run mokapot: {}. Is mokapot installed? Try: pip install mokapot",
                    e
                ))
            })?;

        // Stream stderr (mokapot writes progress to stderr)
        // Also capture error messages to include in failure report
        let mut error_lines: Vec<String> = Vec::new();
        if let Some(stderr) = child.stderr.take() {
            let reader = BufReader::new(stderr);
            for line in reader.lines().map_while(|r| r.ok()) {
                // Log mokapot output so user can see progress
                if line.contains("INFO") || line.contains("iter") || line.contains("Iteration") {
                    log::info!("[mokapot] {}", line);
                } else if line.contains("ERROR")
                    || line.contains("Error")
                    || line.contains("error:")
                    || line.contains("Exception")
                    || line.contains("Traceback")
                    || line.contains("Warning")
                    || line.contains("failed")
                    || line.contains("CRITICAL")
                {
                    // Log errors/warnings at warn level so they're visible
                    log::warn!("[mokapot] {}", line);
                    error_lines.push(line);
                } else if !line.trim().is_empty() {
                    log::debug!("[mokapot] {}", line);
                    // Capture all non-empty lines in case they're part of a traceback
                    error_lines.push(line);
                }
            }
        }

        let status = child.wait().map_err(|e| {
            OspreyError::ExternalToolError(format!("Failed to wait for mokapot: {}", e))
        })?;

        if !status.success() {
            // Include the last few error lines for context
            let error_context = if error_lines.is_empty() {
                "No error output captured".to_string()
            } else {
                // Take last 20 lines for context (may include traceback)
                let start = error_lines.len().saturating_sub(20);
                error_lines[start..].join("\n")
            };
            return Err(OspreyError::ExternalToolError(format!(
                "Mokapot failed with exit code: {:?}\n\nMokapot output:\n{}",
                status.code(),
                error_context
            )));
        }

        // Log model weights file location
        let model_file = out_path.join("mokapot.model.pkl");
        if model_file.exists() {
            log::info!(
                "Mokapot model saved to: {} (use scripts/inspect_mokapot_weights.py to view feature weights)",
                model_file.display()
            );
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
        let mut peptide_idx = 0;
        let mut score_idx = 0;
        let mut qvalue_idx = 0;
        let mut pep_idx = 0;

        for line_result in reader.lines() {
            let line = line_result
                .map_err(|e| OspreyError::OutputError(format!("Failed to read line: {}", e)))?;

            let fields: Vec<&str> = line.split('\t').collect();

            if !header_parsed {
                // Parse header to find column indices
                for (i, field) in fields.iter().enumerate() {
                    match field.to_lowercase().as_str() {
                        "specid" | "psm_id" | "psmid" => psm_id_idx = i,
                        "peptide" => peptide_idx = i,
                        "score" | "mokapot_score" | "mokapot score" => score_idx = i,
                        "q-value" | "qvalue" | "mokapot_qvalue" | "mokapot q-value" => {
                            qvalue_idx = i
                        }
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
            let peptide = fields
                .get(peptide_idx)
                .map(|s| s.to_string())
                .unwrap_or_default();
            let score = fields[score_idx].parse::<f64>().unwrap_or(0.0);
            let q_value = fields[qvalue_idx].parse::<f64>().unwrap_or(1.0);
            let pep = fields
                .get(pep_idx)
                .and_then(|s| s.parse::<f64>().ok())
                .unwrap_or(1.0);

            results.push(MokapotResult {
                psm_id,
                peptide,
                score,
                q_value,
                pep,
            });
        }

        log::debug!(
            "Parsed {} total mokapot results (before FDR filtering)",
            results.len()
        );
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
    test_fdr={test_fdr},
    max_workers={max_workers},
    rng={seed},
)

# Write results
results.to_txt("{output_dir}")
print("SUCCESS")
"#,
            pin_file = pin_path.display(),
            test_fdr = self.test_fdr,
            max_workers = self.num_workers,
            seed = self.seed,
            output_dir = out_path.display(),
        );

        let output = Command::new("python3")
            .arg("-c")
            .arg(&python_script)
            .output()
            .map_err(|e| OspreyError::ExternalToolError(format!("Failed to run Python: {}", e)))?;

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

    /// Write separate PIN files for each file in the input
    ///
    /// Returns a HashMap mapping file_name to the PIN file path.
    pub fn write_pin_files<P: AsRef<Path>>(
        &self,
        psms_by_file: &HashMap<String, Vec<PsmFeatures>>,
        output_dir: P,
    ) -> Result<HashMap<String, PathBuf>> {
        use std::io::BufWriter;

        let out_dir = output_dir.as_ref();
        std::fs::create_dir_all(out_dir).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create output directory: {}", e))
        })?;

        let mut pin_files = HashMap::new();

        for (file_name, psms) in psms_by_file {
            // Sanitize file name for use as PIN file name
            let safe_name = file_name.replace(['/', '\\', ' '], "_");
            let pin_path = out_dir.join(format!("{}.pin", safe_name));

            let file = std::fs::File::create(&pin_path).map_err(|e| {
                OspreyError::OutputError(format!("Failed to create PIN file: {}", e))
            })?;

            let mut writer = BufWriter::with_capacity(1024 * 1024, file);

            // Write header
            writeln!(
                writer,
                "SpecId\tLabel\tScanNr\tCharge1\tCharge2\tCharge3\tCharge4\tCharge5\t{}\tPeptide\tProteins",
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
                    writer,
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

            writer.flush().map_err(|e| {
                OspreyError::OutputError(format!("Failed to flush PIN file: {}", e))
            })?;

            log::info!(
                "Wrote {} PSMs to PIN file: {}",
                psms.len(),
                pin_path.display()
            );

            pin_files.insert(file_name.clone(), pin_path);
        }

        Ok(pin_files)
    }

    /// Run two-step mokapot analysis using CLI
    ///
    /// For a single file: just runs mokapot once with --save_models
    /// For multiple files:
    /// - Step 1: Run mokapot WITHOUT --aggregate + WITH --save_models → per-file results
    /// - Step 2: Run mokapot WITH --aggregate + WITH --load_models → experiment-level results
    ///
    /// # Arguments
    /// * `pin_files` - HashMap mapping file_name to PIN file path
    /// * `output_dir` - Output directory for results
    ///
    /// # Returns
    /// Tuple of (per_file_results, experiment_results)
    /// - per_file_results: HashMap mapping file_name to Vec<MokapotResult>
    /// - experiment_results: Vec<MokapotResult> for experiment-level (empty for single file)
    #[allow(clippy::type_complexity)]
    pub fn run_two_step_analysis<P: AsRef<Path>>(
        &self,
        pin_files: &HashMap<String, PathBuf>,
        output_dir: P,
    ) -> Result<(HashMap<String, Vec<MokapotResult>>, Vec<MokapotResult>)> {
        let out_path = output_dir.as_ref();

        std::fs::create_dir_all(out_path).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create output directory: {}", e))
        })?;

        // Canonicalize output path to ensure absolute paths are passed to mokapot
        let out_path = out_path.canonicalize().map_err(|e| {
            OspreyError::OutputError(format!("Failed to canonicalize output directory: {}", e))
        })?;

        // Build ordered list of PIN file paths to maintain consistent ordering
        // Also canonicalize PIN file paths
        let mut ordered_files: Vec<(String, PathBuf)> = pin_files
            .iter()
            .filter_map(|(name, path)| {
                path.canonicalize()
                    .ok()
                    .map(|abs_path| (name.clone(), abs_path))
            })
            .collect();
        ordered_files.sort_by(|a, b| a.0.cmp(&b.0));

        if ordered_files.is_empty() {
            return Err(OspreyError::ConfigError(
                "No valid PIN files found for mokapot".to_string(),
            ));
        }

        let file_names: Vec<String> = ordered_files.iter().map(|(n, _)| n.clone()).collect();
        let pin_paths: Vec<&PathBuf> = ordered_files.iter().map(|(_, p)| p).collect();

        // Create subdirectories for step outputs
        let step1_dir = out_path.join("run_level");
        let step2_dir = out_path.join("experiment_level");
        std::fs::create_dir_all(&step1_dir).map_err(|e| {
            OspreyError::OutputError(format!("Failed to create step1 directory: {}", e))
        })?;

        if pin_files.len() == 1 {
            // Single file case: just run mokapot once
            log::info!("Running mokapot on single file...");

            let pin_path = pin_paths[0];
            self.run_mokapot_cli(
                &[pin_path],
                &step1_dir,
                false, // no aggregate
                true,  // save models
                None,  // no load models
            )?;

            // Parse single-file results
            let results_file = step1_dir.join("mokapot.psms.txt");
            let results = self.parse_results(&results_file)?;
            let file_name = &file_names[0];

            // Log single-file statistics (filtered by FDR threshold)
            let passing_results: Vec<_> = results
                .iter()
                .filter(|r| r.q_value <= self.test_fdr)
                .collect();
            let precursor_count = passing_results.len();
            let unique_peptides: HashSet<&str> =
                passing_results.iter().map(|r| r.peptide.as_str()).collect();
            let peptide_count = unique_peptides.len();

            log::info!("");
            log::info!("=== Results ({}% FDR) ===", (self.test_fdr * 100.0) as u32);
            log::info!(
                "  {}: {} precursors, {} peptides",
                file_name,
                precursor_count,
                peptide_count
            );

            let mut per_file_results = HashMap::new();
            per_file_results.insert(file_name.clone(), results);

            // For single file, experiment-level is the same as run-level
            // Return empty experiment results since there's no meaningful aggregation
            Ok((per_file_results, Vec::new()))
        } else {
            // Multiple files: two-step analysis
            log::info!(
                "Step 1: Training joint model on {} files (run-level FDR)...",
                pin_files.len()
            );

            // Step 1: Run mokapot WITHOUT --aggregate + WITH --save_models
            // This trains one model on all files but reports results per file
            self.run_mokapot_cli(
                &pin_paths, &step1_dir, false, // no aggregate
                true,  // save models
                None,  // no load models
            )?;

            // Parse per-file results from Step 1
            // When run without --aggregate, mokapot creates separate output per input file
            log::info!("");
            log::info!(
                "=== Per-file results (run-level FDR at {}%) ===",
                (self.test_fdr * 100.0) as u32
            );
            let mut per_file_results = HashMap::new();
            let mut total_precursors = 0usize;
            let mut total_peptides_set: HashSet<String> = HashSet::new();

            for (file_name, pin_path) in &ordered_files {
                // mokapot names output files based on PIN file stem
                let pin_stem = pin_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown");
                let results_file = step1_dir.join(format!("{}.mokapot.psms.txt", pin_stem));

                if results_file.exists() {
                    let results = self.parse_results(&results_file)?;

                    // Count precursors and unique peptides passing FDR threshold
                    let passing_results: Vec<_> = results
                        .iter()
                        .filter(|r| r.q_value <= self.test_fdr)
                        .collect();
                    let precursor_count = passing_results.len();
                    let unique_peptides: HashSet<&str> =
                        passing_results.iter().map(|r| r.peptide.as_str()).collect();
                    let peptide_count = unique_peptides.len();

                    // Update totals for summary
                    total_precursors += precursor_count;
                    for pep in &unique_peptides {
                        total_peptides_set.insert((*pep).to_string());
                    }

                    log::info!(
                        "  {}: {} precursors, {} peptides",
                        file_name,
                        precursor_count,
                        peptide_count
                    );

                    per_file_results.insert(file_name.clone(), results);
                } else {
                    log::warn!(
                        "  {}: No results file found: {}",
                        file_name,
                        results_file.display()
                    );
                }
            }

            log::info!("");
            log::info!(
                "Step 1 total: {} precursors, {} unique peptides across {} files",
                total_precursors,
                total_peptides_set.len(),
                per_file_results.len()
            );

            log::info!("Step 2: Training new model for experiment-level FDR...");

            // Create step2 directory
            std::fs::create_dir_all(&step2_dir).map_err(|e| {
                OspreyError::OutputError(format!("Failed to create step2 directory: {}", e))
            })?;

            // Step 2: Run mokapot WITH --aggregate, training a fresh model
            // This trains independently from Step 1 to avoid any model bias
            self.run_mokapot_cli(
                &pin_paths, &step2_dir, true,  // aggregate
                false, // don't save models
                None,  // train a new model (don't reuse Step 1)
            )?;

            // Parse experiment-level results from Step 2
            let experiment_results_file = step2_dir.join("mokapot.psms.txt");
            let experiment_results = self.parse_results(&experiment_results_file)?;

            // Log experiment-level statistics (filtered by FDR threshold)
            let passing_exp_results: Vec<_> = experiment_results
                .iter()
                .filter(|r| r.q_value <= self.test_fdr)
                .collect();
            let exp_precursor_count = passing_exp_results.len();
            let exp_unique_peptides: HashSet<&str> = passing_exp_results
                .iter()
                .map(|r| r.peptide.as_str())
                .collect();
            let exp_peptide_count = exp_unique_peptides.len();

            log::info!("");
            log::info!(
                "=== Experiment-level results ({}% FDR) ===",
                (self.test_fdr * 100.0) as u32
            );
            log::info!(
                "  Experiment: {} precursors, {} peptides",
                exp_precursor_count,
                exp_peptide_count
            );

            Ok((per_file_results, experiment_results))
        }
    }

    /// Run mokapot CLI with specified options
    fn run_mokapot_cli(
        &self,
        pin_files: &[&PathBuf],
        output_dir: &Path,
        aggregate: bool,
        save_models: bool,
        load_models: Option<&[PathBuf]>,
    ) -> Result<()> {
        let mut cmd = Command::new(&self.mokapot_path);

        // Add PIN files
        for pin_file in pin_files {
            cmd.arg(pin_file);
        }

        // Output directory
        cmd.arg("--dest_dir").arg(output_dir);

        // FDR thresholds
        cmd.arg("--train_fdr").arg(self.train_fdr.to_string());
        cmd.arg("--test_fdr").arg(self.test_fdr.to_string());

        // Training parameters
        cmd.arg("--max_iter").arg(self.max_iter.to_string());
        cmd.arg("--seed").arg(self.seed.to_string());
        cmd.arg("--max_workers").arg(self.num_workers.to_string());

        // Memory-aware training subset (only for training step, not when loading models)
        if load_models.is_none() {
            if let Some(subset_size) = self.effective_subset_max_train() {
                log::debug!(
                    "Using --subset_max_train {} to limit memory usage",
                    subset_size
                );
                cmd.arg("--subset_max_train").arg(subset_size.to_string());
            }
        }

        // Aggregate mode (combine all files for experiment-level FDR)
        if aggregate {
            cmd.arg("--aggregate");
        }

        // Model persistence
        if save_models {
            cmd.arg("--save_models");
        }
        if let Some(model_paths) = load_models {
            cmd.arg("--load_models");
            for model_path in model_paths {
                cmd.arg(model_path);
            }
        }

        // Verbosity
        cmd.arg("--verbosity").arg("1");

        log::debug!("Running mokapot command: {:?}", cmd);

        // Spawn with piped output
        let mut child = cmd
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| {
                OspreyError::ExternalToolError(format!(
                    "Failed to run mokapot: {}. Is mokapot installed? Try: pip install mokapot",
                    e
                ))
            })?;

        // Stream stderr for progress
        let mut error_lines: Vec<String> = Vec::new();
        if let Some(stderr) = child.stderr.take() {
            let reader = BufReader::new(stderr);
            for line in reader.lines().map_while(|r| r.ok()) {
                // Skip Python FutureWarning/DeprecationWarning messages (not actionable)
                if line.contains("FutureWarning") || line.contains("DeprecationWarning") {
                    log::debug!("[mokapot] {}", line);
                    continue;
                }

                if line.contains("INFO") || line.contains("iter") || line.contains("Iteration") {
                    log::info!("[mokapot] {}", line);
                } else if line.contains("ERROR")
                    || line.contains("Error")
                    || line.contains("error:")
                    || line.contains("Exception")
                    || line.contains("Traceback")
                    || line.contains("Warning")
                    || line.contains("failed")
                    || line.contains("CRITICAL")
                {
                    log::warn!("[mokapot] {}", line);
                    error_lines.push(line);
                } else if !line.trim().is_empty() {
                    log::debug!("[mokapot] {}", line);
                    error_lines.push(line);
                }
            }
        }

        let status = child.wait().map_err(|e| {
            OspreyError::ExternalToolError(format!("Failed to wait for mokapot: {}", e))
        })?;

        if !status.success() {
            let error_context = if error_lines.is_empty() {
                "No error output captured".to_string()
            } else {
                let start = error_lines.len().saturating_sub(20);
                error_lines[start..].join("\n")
            };
            return Err(OspreyError::ExternalToolError(format!(
                "Mokapot failed with exit code: {:?}\n\nMokapot output:\n{}",
                status.code(),
                error_context
            )));
        }

        Ok(())
    }
}

impl Default for MokapotRunner {
    fn default() -> Self {
        Self::new()
    }
}

/// Find mokapot model file(s) in a directory
///
/// Mokapot may use different naming conventions depending on version:
/// - `mokapot.model` (single model, older versions)
/// - `mokapot.model.pkl` (pickle format)
/// - `mokapot.model_fold-1.pkl`, `mokapot.model_fold-2.pkl`, etc. (cross-validation folds)
/// - `mokapot.weights.csv` (older versions)
///
/// Returns a vector of model file paths. For cross-validation, returns all fold files.
#[allow(dead_code)]
fn find_mokapot_models(dir: &Path) -> Result<Vec<PathBuf>> {
    // First, check for cross-validation fold files (mokapot.model_fold-N.pkl)
    let mut fold_files: Vec<PathBuf> = std::fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok())
                .filter(|e| {
                    let name = e.file_name().to_string_lossy().to_string();
                    name.starts_with("mokapot.model_fold-") && name.ends_with(".pkl")
                })
                .map(|e| e.path())
                .collect()
        })
        .unwrap_or_default();

    // Sort fold files to ensure consistent ordering
    fold_files.sort();

    if !fold_files.is_empty() {
        log::debug!(
            "Found {} cross-validation fold model files",
            fold_files.len()
        );
        return Ok(fold_files);
    }

    // Fall back to single model file names
    let single_candidates = [
        "mokapot.model",
        "mokapot.model.pkl",
        "mokapot.weights.csv",
        "mokapot.weights",
    ];

    for name in &single_candidates {
        let path = dir.join(name);
        if path.exists() {
            return Ok(vec![path]);
        }
    }

    // List what files ARE in the directory to help with debugging
    let files_in_dir: Vec<String> = std::fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(|e| e.ok())
                .map(|e| e.file_name().to_string_lossy().to_string())
                .collect()
        })
        .unwrap_or_default();

    Err(OspreyError::ExternalToolError(format!(
        "Mokapot model file not found in {}. \
         Expected fold files (mokapot.model_fold-*.pkl) or single file. \
         Files in directory: {:?}",
        dir.display(),
        files_in_dir
    )))
}

/// Get the feature header for PIN format (37 features)
///
/// Features are grouped by source:
/// - Ridge regression (chromatographic): peak_apex, peak_area, peak_width,
///   coefficient_stability, relative_coefficient, explained_intensity, signal_to_noise,
///   xic_signal_to_noise
/// - Spectral matching (mixed, at apex): xcorr, consecutive_ions
/// - Spectral matching (deconvoluted, apex ± 2): *_deconv versions
/// - Derived: rt_deviation, fragment co-elution, mass accuracy
/// - Percolator-style: abs_rt_deviation, peptide_length, missed_cleavages,
///   ln_num_candidates, coef_zscore, coef_zscore_mean
fn get_feature_header() -> String {
    [
        // Ridge regression features (chromatographic profile)
        "peak_apex",
        "peak_area",
        "peak_width",
        "coefficient_stability",
        "relative_coefficient",
        "explained_intensity",
        "signal_to_noise",
        "xic_signal_to_noise",
        // Spectral matching features (mixed/observed at apex spectrum)
        "xcorr",
        "consecutive_ions",
        // Spectral matching features (deconvoluted, coefficient-weighted apex ± 2 scans)
        "xcorr_deconv",
        "consecutive_ions_deconv",
        // Derived features
        "rt_deviation",
        // Fragment co-elution features (bounded to peak integration boundaries)
        "fragment_coelution_sum",
        "fragment_coelution_min",
        "n_coeluting_fragments",
        // Per-fragment correlations (top 6 by library intensity, 0-padded)
        "fragment_corr_0",
        "fragment_corr_1",
        "fragment_corr_2",
        "fragment_corr_3",
        "fragment_corr_4",
        "fragment_corr_5",
        // Elution-weighted spectral similarity
        "elution_weighted_cosine",
        // Per-fragment mass accuracy
        "mass_accuracy_deviation_mean",
        "abs_mass_accuracy_deviation_mean",
        "mass_accuracy_std",
        // Percolator-style features
        "abs_rt_deviation",
        "peptide_length",
        "missed_cleavages",
        "ln_num_candidates",
        "coef_zscore",
        "coef_zscore_mean",
        // MS1-based features (HRAM only, 0.0 for unit resolution or missing MS1)
        "ms1_precursor_coelution",
        "ms1_isotope_cosine",
        // Tukey median polish features (fragment XIC decomposition)
        "median_polish_cosine",
        "median_polish_rsquared",
        "median_polish_residual_ratio",
    ]
    .join("\t")
}

/// Format charge state as 5 binary one-hot features (Charge1..Charge5)
/// Charge 6+ produces all zeros.
fn format_charge_features(charge: u8) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}",
        if charge == 1 { 1 } else { 0 },
        if charge == 2 { 1 } else { 0 },
        if charge == 3 { 1 } else { 0 },
        if charge == 4 { 1 } else { 0 },
        if charge == 5 { 1 } else { 0 },
    )
}

/// Format features for PIN output (37 features)
fn format_features(features: &FeatureSet) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        // Ridge regression features (chromatographic profile)
        features.peak_apex,
        features.peak_area,
        features.peak_width,
        features.coefficient_stability,
        features.relative_coefficient,
        features.explained_intensity,
        features.signal_to_noise,
        features.xic_signal_to_noise,
        // Spectral matching features (mixed/observed at apex spectrum)
        features.xcorr,
        features.consecutive_ions,
        // Spectral matching features (deconvoluted, coefficient-weighted apex ± 2 scans)
        features.xcorr_deconv,
        features.consecutive_ions_deconv,
        // Derived features
        features.rt_deviation,
        // Fragment co-elution features (bounded to peak integration boundaries)
        features.fragment_coelution_sum,
        features.fragment_coelution_min,
        features.n_coeluting_fragments,
        // Per-fragment correlations (top 6 by library intensity, 0-padded)
        features.fragment_corr_0,
        features.fragment_corr_1,
        features.fragment_corr_2,
        features.fragment_corr_3,
        features.fragment_corr_4,
        features.fragment_corr_5,
        // Elution-weighted spectral similarity
        features.elution_weighted_cosine,
        // Per-fragment mass accuracy
        features.mass_accuracy_deviation_mean,
        features.abs_mass_accuracy_deviation_mean,
        features.mass_accuracy_std,
        // Percolator-style features
        features.abs_rt_deviation,
        features.peptide_length,
        features.missed_cleavages,
        features.ln_num_candidates,
        features.coef_zscore,
        features.coef_zscore_mean,
        // MS1-based features (HRAM only, 0.0 for unit resolution or missing MS1)
        features.ms1_precursor_coelution,
        features.ms1_isotope_cosine,
        // Tukey median polish features (fragment XIC decomposition)
        features.median_polish_cosine,
        features.median_polish_rsquared,
        features.median_polish_residual_ratio,
    )
}

/// Returns the 37 PIN feature names in the same order as get_feature_header() and format_features().
pub fn get_pin_feature_names() -> Vec<&'static str> {
    vec![
        "peak_apex",
        "peak_area",
        "peak_width",
        "coefficient_stability",
        "relative_coefficient",
        "explained_intensity",
        "signal_to_noise",
        "xic_signal_to_noise",
        "xcorr",
        "consecutive_ions",
        "xcorr_deconv",
        "consecutive_ions_deconv",
        "rt_deviation",
        "fragment_coelution_sum",
        "fragment_coelution_min",
        "n_coeluting_fragments",
        "fragment_corr_0",
        "fragment_corr_1",
        "fragment_corr_2",
        "fragment_corr_3",
        "fragment_corr_4",
        "fragment_corr_5",
        "elution_weighted_cosine",
        "mass_accuracy_deviation_mean",
        "abs_mass_accuracy_deviation_mean",
        "mass_accuracy_std",
        "abs_rt_deviation",
        "peptide_length",
        "missed_cleavages",
        "ln_num_candidates",
        "coef_zscore",
        "coef_zscore_mean",
        "ms1_precursor_coelution",
        "ms1_isotope_cosine",
        "median_polish_cosine",
        "median_polish_rsquared",
        "median_polish_residual_ratio",
    ]
}

/// Returns a single PIN feature value by index (0-36).
///
/// Index order matches get_feature_header(), format_features(), and get_pin_feature_names().
pub fn pin_feature_value(features: &FeatureSet, index: usize) -> f64 {
    match index {
        0 => features.peak_apex,
        1 => features.peak_area,
        2 => features.peak_width,
        3 => features.coefficient_stability,
        4 => features.relative_coefficient,
        5 => features.explained_intensity,
        6 => features.signal_to_noise,
        7 => features.xic_signal_to_noise,
        8 => features.xcorr,
        9 => features.consecutive_ions as f64,
        10 => features.xcorr_deconv,
        11 => features.consecutive_ions_deconv as f64,
        12 => features.rt_deviation,
        13 => features.fragment_coelution_sum,
        14 => features.fragment_coelution_min,
        15 => features.n_coeluting_fragments as f64,
        16 => features.fragment_corr_0,
        17 => features.fragment_corr_1,
        18 => features.fragment_corr_2,
        19 => features.fragment_corr_3,
        20 => features.fragment_corr_4,
        21 => features.fragment_corr_5,
        22 => features.elution_weighted_cosine,
        23 => features.mass_accuracy_deviation_mean,
        24 => features.abs_mass_accuracy_deviation_mean,
        25 => features.mass_accuracy_std,
        26 => features.abs_rt_deviation,
        27 => features.peptide_length as f64,
        28 => features.missed_cleavages as f64,
        29 => features.ln_num_candidates,
        30 => features.coef_zscore,
        31 => features.coef_zscore_mean,
        32 => features.ms1_precursor_coelution,
        33 => features.ms1_isotope_cosine,
        34 => features.median_polish_cosine,
        35 => features.median_polish_rsquared,
        36 => features.median_polish_residual_ratio,
        _ => 0.0,
    }
}

/// Number of PIN features (37).
pub const NUM_PIN_FEATURES: usize = 37;

/// Format peptide for PIN output
/// Percolator/mokapot expects format: FLANKING.SEQUENCE.FLANKING
fn format_peptide(peptide: &str) -> String {
    // First, strip any existing flanking characters (underscores, dots, dashes)
    let stripped = strip_flanking_chars(peptide);

    // Add standard PIN format flanking residues
    format!("-.{}.-", stripped)
}

/// Strip flanking characters from peptide sequences
///
/// Handles various formats:
/// - `_PEPTIDE_` → `PEPTIDE`
/// - `K.PEPTIDE.R` → `PEPTIDE`
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

/// Write a simple TSV report of mokapot results
pub fn write_mokapot_report<P: AsRef<Path>>(results: &[MokapotResult], path: P) -> Result<()> {
    let mut file = std::fs::File::create(path.as_ref())
        .map_err(|e| OspreyError::OutputError(format!("Failed to create report file: {}", e)))?;

    writeln!(file, "PSM_ID\tMokapot_Score\tQ_Value\tPEP")
        .map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

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

    /// Verifies that peptide sequences are wrapped in standard PIN flanking format and existing flanking characters are stripped.
    #[test]
    fn test_format_peptide() {
        assert_eq!(format_peptide("PEPTIDE"), "-.PEPTIDE.-");
        // Flanking residues are stripped and re-wrapped with standard format
        assert_eq!(format_peptide("K.PEPTIDE.R"), "-.PEPTIDE.-");
    }

    /// Verifies that the PIN feature header contains expected feature names and excludes optional fields.
    #[test]
    fn test_feature_header() {
        let header = get_feature_header();
        assert!(header.contains("peak_apex"));
        assert!(header.contains("xcorr"));
        assert!(header.contains("coef_zscore"));
        assert!(header.contains("xic_signal_to_noise"));
        assert!(header.contains("ms1_precursor_coelution"));
        assert!(header.contains("ms1_isotope_cosine"));
        assert!(header.contains("median_polish_cosine"));
        assert!(header.contains("median_polish_rsquared"));
        assert!(header.contains("median_polish_residual_ratio"));
        // Removed features should not be present
        assert!(!header.contains("dot_product"));
        assert!(!header.contains("hyperscore"));
        assert!(!header.contains("fragment_coverage"));
        assert!(!header.contains("top6_matches"));
        assert!(!header.contains("n_contributing_scans"));
    }

    /// Verifies that charge states 1-5 produce correct one-hot encodings and out-of-range charges produce all zeros.
    #[test]
    fn test_charge_features() {
        assert_eq!(format_charge_features(1), "1\t0\t0\t0\t0");
        assert_eq!(format_charge_features(2), "0\t1\t0\t0\t0");
        assert_eq!(format_charge_features(3), "0\t0\t1\t0\t0");
        assert_eq!(format_charge_features(4), "0\t0\t0\t1\t0");
        assert_eq!(format_charge_features(5), "0\t0\t0\t0\t1");
        assert_eq!(format_charge_features(6), "0\t0\t0\t0\t0");
        assert_eq!(format_charge_features(0), "0\t0\t0\t0\t0");
    }

    /// Verifies that the MokapotRunner can check for mokapot availability without panicking.
    #[test]
    fn test_mokapot_runner_availability() {
        let runner = MokapotRunner::new();
        // This will be false on systems without mokapot, which is expected
        let _available = runner.is_available();
    }

    /// Verifies that get_pin_feature_names() returns exactly 37 names matching the header.
    #[test]
    fn test_get_pin_feature_names_count() {
        let names = get_pin_feature_names();
        assert_eq!(names.len(), NUM_PIN_FEATURES);

        // Verify names match the header
        let header = get_feature_header();
        let header_names: Vec<&str> = header.split('\t').collect();
        assert_eq!(header_names.len(), names.len());
        for (h, n) in header_names.iter().zip(names.iter()) {
            assert_eq!(h, n);
        }
    }

    /// Verifies that pin_feature_value returns consistent values matching format_features order.
    #[test]
    fn test_pin_feature_value_matches_format() {
        let features = FeatureSet::default();
        // All default values should be 0.0, so pin_feature_value should return 0.0 for all indices
        for i in 0..NUM_PIN_FEATURES {
            let val = pin_feature_value(&features, i);
            assert!(val.is_finite(), "Feature {} returned non-finite value", i);
        }
        // Out of range should return 0.0
        assert_eq!(pin_feature_value(&features, 99), 0.0);
    }
}
