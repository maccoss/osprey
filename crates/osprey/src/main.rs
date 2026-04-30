//! Osprey CLI - Command-line interface for peptide-centric DIA analysis

use anyhow::Result;
use clap::Parser;
use log::LevelFilter;
use osprey::{
    run_analysis, ConfigOverrides, FdrLevel, FdrMethod, OspreyConfig, ParquetCompression,
    ResolutionMode, SharedPeptideMode, ToleranceUnit,
};
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

/// Format duration into human-readable string
fn format_duration(duration: std::time::Duration) -> String {
    let total_seconds = duration.as_secs();
    let days = total_seconds / 86400;
    let hours = (total_seconds % 86400) / 3600;
    let minutes = (total_seconds % 3600) / 60;
    let seconds = total_seconds % 60;
    let millis = duration.subsec_millis();

    if days > 0 {
        if hours > 0 {
            format!("{} days {} hours", days, hours)
        } else {
            format!("{} days", days)
        }
    } else if hours > 0 {
        if minutes > 0 {
            format!("{} hours {} minutes", hours, minutes)
        } else {
            format!("{} hours", hours)
        }
    } else if minutes > 0 {
        if seconds > 0 {
            format!("{} minutes {} seconds", minutes, seconds)
        } else {
            format!("{} minutes", minutes)
        }
    } else if seconds > 0 {
        format!("{}.{:03} seconds", seconds, millis)
    } else {
        format!("{} ms", millis)
    }
}

/// Osprey: Peptide-centric DIA analysis tool
#[derive(Parser, Debug)]
#[command(name = "osprey")]
#[command(
    version,
    about = "Peptide-centric DIA analysis with Skyline integration"
)]
#[command(long_about = r#"
Osprey is an open-source tool for peptide detection and quantification
in data-independent acquisition (DIA) mass spectrometry data.

It uses fragment XIC co-elution analysis to detect peptides in DIA data,
with machine learning scoring and rigorous FDR control.

EXAMPLES:
    # Basic analysis with DIA-NN library
    osprey -i sample.mzML -l library.tsv -o results.blib

    # Multiple input files
    osprey -i *.mzML -l library.tsv -o results.blib

    # Write additional TSV report
    osprey -i sample.mzML -l library.tsv -o results.blib --report results.tsv
"#)]
struct Args {
    /// Configuration file (YAML format)
    #[arg(short, long)]
    config: Option<PathBuf>,

    /// Generate a template configuration file
    #[arg(long)]
    generate_config: Option<PathBuf>,

    /// Input mzML file(s)
    #[arg(short, long, num_args = 1..)]
    input: Option<Vec<PathBuf>>,

    /// Spectral library file (.tsv for DIA-NN, .blib, or .elib)
    #[arg(short, long)]
    library: Option<PathBuf>,

    /// Output results file (.blib format for Skyline)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Resolution mode: unit, hram, auto
    #[arg(long, default_value = "auto")]
    resolution: String,

    /// Fragment m/z tolerance (e.g., 10 for 10 ppm, or 0.3 for 0.3 Th)
    #[arg(long)]
    fragment_tolerance: Option<f64>,

    /// Fragment tolerance unit: ppm, mz (Thompson)
    #[arg(long)]
    fragment_unit: Option<String>,

    /// Precursor m/z tolerance (e.g., 10 for 10 ppm, or 1.0 for 1.0 Th)
    #[arg(long)]
    precursor_tolerance: Option<f64>,

    /// Precursor tolerance unit: ppm, mz (Thompson)
    #[arg(long)]
    precursor_unit: Option<String>,

    /// RT tolerance in minutes (fallback when calibration disabled)
    #[arg(long, default_value_t = 2.0)]
    rt_tolerance: f64,

    /// Disable RT calibration (use fixed rt_tolerance instead)
    #[arg(long)]
    no_rt_calibration: bool,

    /// Run-level FDR threshold
    #[arg(long, default_value_t = 0.01)]
    run_fdr: f64,

    /// Experiment-level FDR threshold (for multi-file analyses)
    #[arg(long, default_value_t = 0.01)]
    experiment_fdr: f64,

    /// Number of threads (default: all available)
    #[arg(long)]
    threads: Option<usize>,

    /// Write TSV report to this file
    #[arg(long)]
    report: Option<PathBuf>,

    /// FDR method: percolator (native SVM, default), mokapot (external Python), or simple (no ML)
    #[arg(long, default_value = "percolator")]
    fdr_method: String,

    /// FDR filtering level: precursor (default), peptide, protein, or both.
    /// Controls which peptide identities are eligible for blib output;
    /// precursor-level FDR is always enforced within each eligible peptide.
    /// 'protein' requires --protein-fdr.
    #[arg(long)]
    fdr_level: Option<String>,

    /// Protein-level FDR threshold [default: 0.01]. Protein parsimony and
    /// picked-protein FDR always run; this sets the cutoff for reporting and
    /// for the compaction/reconciliation rescue rule.
    #[arg(long)]
    protein_fdr: Option<f64>,

    /// How to handle shared peptides for protein FDR: all (default), razor, or unique
    #[arg(long)]
    shared_peptides: Option<String>,

    /// Write PIN files for external tools
    #[arg(long)]
    write_pin: bool,

    /// Disable the coelution signal pre-filter (3-of-4 consecutive scans with ≥2 top-6
    /// fragments). The pre-filter speeds up HRAM searches ~30% with minimal sensitivity loss.
    #[arg(long)]
    no_prefilter: bool,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,

    /// HPC modifier: run only the per-file fan-out from the entry
    /// point, skipping the next join. With `-i ...` (Stage 1 entry):
    /// runs Stages 1-4 only — per-file `.scores.parquet` is written
    /// next to each input mzML, no FDR, no blib. With
    /// `--join-at-pass=<N>`: skip the join at <N> and run only the
    /// per-file phase that follows it. Mutually exclusive with
    /// --join-only.
    #[arg(long)]
    no_join: bool,

    /// HPC modifier: run only the join phase at the entry point,
    /// stopping before the per-file fan-out that follows. Used when an
    /// HPC coordinator finishes a join (e.g. Stage 5 first-pass
    /// Percolator + reconciliation planning) and wants to ship the
    /// resulting plan files to N worker nodes for the per-file phase.
    /// Requires `--join-at-pass=<N>`; not valid alone. Mutually
    /// exclusive with --no-join.
    #[arg(long)]
    join_only: bool,

    /// HPC: enter the pipeline at a specific join checkpoint, consuming
    /// the per-file Parquet score files for that pass via
    /// --input-scores. `1` = consume Stage 4 outputs (raw scoring) and
    /// run Stages 5-8. `2` = consume Stage 6 outputs (reconciled
    /// scoring) and run Stages 7-8. Combine with `--join-only` to stop
    /// after the next join, or `--no-join` to skip the next join and run
    /// only the per-file phase. Mutually exclusive with --input.
    #[arg(long, value_name = "N")]
    join_at_pass: Option<u8>,

    /// HPC: one or more `.scores.parquet` files (or a single directory
    /// scanned non-recursively for `*.scores.parquet`). Required when
    /// --join-at-pass is set; ignored otherwise.
    #[arg(long, num_args = 1..)]
    input_scores: Option<Vec<PathBuf>>,

    /// Compression codec for `.scores.parquet` writes. Default `zstd`
    /// keeps the historical Osprey production behavior. `snappy` is
    /// offered for cross-impl interop with OspreySharp (Parquet.Net 3.x
    /// supports Snappy only). Reading is always auto-dispatched on
    /// per-column-chunk metadata, so this flag only governs writes.
    #[arg(long)]
    parquet_compression: Option<String>,
}

/// Normalize the HPC entry-point + modifier flags before validation.
/// Resolves the (`--join-at-pass=<N>`, `--join-only`, `--no-join`) triple
/// into the single `args.join_only` boolean that downstream code already
/// reads (and leaves `args.no_join` alone since its semantics are
/// unchanged for the Stage 1 entry point).
///
/// Modifier semantics:
/// - `--join-only` runs only the next join from the entry point.
/// - `--no-join` runs only the per-file fan-out from the entry point.
///
/// Entry points:
/// - `-i ...` (no `--join-at-pass`): Stage 1 (raw mzML).
/// - `--join-at-pass=1 --input-scores ...`: post-Stage-4 (raw parquets).
/// - `--join-at-pass=2 --input-scores ...`: post-Stage-6 (reconciled
///   parquets).
///
/// PR 1 (this commit) wires the rename only — combinations that need
/// the Stage 5 → Stage 6 boundary persistence (`--join-at-pass=1`
/// combined with either modifier) error out as "not yet implemented",
/// and `--join-at-pass=2` errors the same way until the Stage 6 →
/// Stage 7 path lands.
fn normalize_hpc_args(args: &mut Args) -> Result<()> {
    // Modifiers are mutually exclusive: can't be both per-file-only and
    // join-only simultaneously.
    if args.no_join && args.join_only {
        anyhow::bail!("--no-join and --join-only are mutually exclusive modifiers.");
    }

    // `--join-only` is a modifier of `--join-at-pass=<N>`; standalone use
    // has no entry point to modify. The old standalone spelling that meant
    // "run Stages 5-8 from Stage 4 parquets" is now `--join-at-pass=1`.
    if args.join_only && args.join_at_pass.is_none() {
        anyhow::bail!(
            "--join-only is a modifier and requires --join-at-pass=<N>. \
             To run Stages 5-8 from Stage 4 parquets, use --join-at-pass=1 --input-scores ..."
        );
    }

    match args.join_at_pass {
        Some(1) => {
            // PR 2 will implement these modifier combinations against
            // persisted Stage 5 → Stage 6 boundary files.
            if args.join_only {
                anyhow::bail!(
                    "--join-at-pass=1 --join-only (run only Stage 5) is not yet implemented."
                );
            }
            if args.no_join {
                anyhow::bail!(
                    "--join-at-pass=1 --no-join (run only Stage 6 from persisted Stage 5 outputs) is not yet implemented."
                );
            }
            // Plain `--join-at-pass=1`: route through the existing
            // Stage 5+ entry path that reads `args.join_only`.
            args.join_only = true;
        }
        Some(2) => {
            anyhow::bail!(
                "--join-at-pass=2 (Stage 6 reconciled-parquet input) is not yet implemented."
            );
        }
        Some(n) => {
            anyhow::bail!("--join-at-pass must be 1 or 2 (got {}).", n);
        }
        None => {
            // Stage 1 entry point (with `-i ...`). `--no-join` keeps its
            // existing meaning: do per-file work only, which is exactly
            // Stages 1-4. No remapping needed.
        }
    }
    Ok(())
}

/// Validate the HPC mode flags (`--no-join`, `--join-only`,
/// `--input-scores`) for mutual exclusion and required-companion errors.
/// Does not warn (warnings stay in `main`).
fn validate_hpc_args(args: &Args) -> Result<()> {
    if args.no_join && args.join_only {
        anyhow::bail!("--no-join and --join-only are mutually exclusive.");
    }
    if args.join_only {
        // Reachable only via `--join-at-pass=1` after `normalize_hpc_args`
        // (standalone `--join-only` errors there). Error text references the
        // canonical flag the user typed.
        if args.input.is_some() {
            anyhow::bail!(
                "--join-at-pass=1 cannot be combined with --input. Use --input-scores instead."
            );
        }
        if args.input_scores.is_none() {
            anyhow::bail!("--join-at-pass=1 requires --input-scores <path...>.");
        }
        if args.library.is_none() || args.output.is_none() {
            anyhow::bail!("--join-at-pass=1 requires --library and --output.");
        }
    }
    if args.no_join {
        if args.input_scores.is_some() {
            anyhow::bail!("--no-join cannot be combined with --input-scores.");
        }
        if args.input.is_none() {
            anyhow::bail!("--no-join requires --input <mzML...>.");
        }
    }
    Ok(())
}

/// Expand `--input-scores` arguments: a single directory becomes the
/// non-recursive list of `*.scores.parquet` files in it; explicit file
/// paths are passed through unchanged. Returns an error if the directory
/// is empty or any explicit path doesn't exist.
fn resolve_input_scores(paths: Vec<PathBuf>) -> Result<Vec<PathBuf>> {
    if paths.len() == 1 && paths[0].is_dir() {
        let dir = &paths[0];
        let mut found: Vec<PathBuf> = std::fs::read_dir(dir)
            .map_err(|e| anyhow::anyhow!("Failed to read --input-scores dir {:?}: {}", dir, e))?
            .filter_map(|entry| entry.ok().map(|e| e.path()))
            .filter(|p| p.is_file() && p.to_string_lossy().ends_with(".scores.parquet"))
            .collect();
        if found.is_empty() {
            anyhow::bail!(
                "No *.scores.parquet files found in --input-scores directory {:?}",
                dir
            );
        }
        found.sort();
        return Ok(found);
    }
    for p in &paths {
        if !p.exists() {
            anyhow::bail!("--input-scores path not found: {:?}", p);
        }
    }
    Ok(paths)
}

fn main() -> Result<()> {
    let mut args = Args::parse();

    // Initialize two-tier logging: clean terminal output + verbose log file.
    // Terminal shows info-level with elapsed time prefixes by default.
    // Log file always captures debug-level with full timestamps.
    // --verbose makes the terminal show everything too.
    let terminal_level = if args.verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    let log_filename = format!(
        "osprey_{}.log",
        chrono::Local::now().format("%Y-%m-%d_%H%M%S")
    );
    let log_file = std::fs::File::create(&log_filename)?;
    osprey::logging::init(terminal_level, log_file);

    // Handle --generate-config
    if let Some(config_path) = args.generate_config {
        log::info!(
            "Generating template configuration file: {}",
            config_path.display()
        );
        OspreyConfig::create_template(&config_path)?;
        log::info!(
            "Template created. Edit the file and run with --config {}",
            config_path.display()
        );
        return Ok(());
    }

    normalize_hpc_args(&mut args)?;
    validate_hpc_args(&args)?;
    if args.no_join && args.output.is_some() {
        log::warn!(
            "--no-join: --output is ignored (no blib is written). \
             Per-file `.scores.parquet` files will be written next to each input mzML."
        );
    }

    // Load config from file or build from CLI args
    let mut config = if let Some(config_path) = &args.config {
        log::info!("Loading configuration from: {}", config_path.display());
        OspreyConfig::from_yaml(config_path)?
    } else {
        // Require input, library, output if no config file
        // (skipped in HPC modes — handled above)
        if !args.no_join
            && !args.join_only
            && (args.input.is_none() || args.library.is_none() || args.output.is_none())
        {
            anyhow::bail!(
                "Either --config <file> or --input, --library, and --output are required.\n\
                 Use --generate-config <file> to create a template configuration file."
            );
        }
        OspreyConfig::default()
    };

    // Parse tolerance units if provided
    let fragment_unit = args
        .fragment_unit
        .as_ref()
        .map(|s| match s.to_lowercase().as_str() {
            "ppm" => ToleranceUnit::Ppm,
            "mz" | "th" | "da" => ToleranceUnit::Mz,
            _ => {
                log::warn!("Unknown fragment unit '{}', defaulting to ppm", s);
                ToleranceUnit::Ppm
            }
        });

    let precursor_unit = args
        .precursor_unit
        .as_ref()
        .map(|s| match s.to_lowercase().as_str() {
            "ppm" => ToleranceUnit::Ppm,
            "mz" | "th" | "da" => ToleranceUnit::Mz,
            _ => {
                log::warn!("Unknown precursor unit '{}', defaulting to ppm", s);
                ToleranceUnit::Ppm
            }
        });

    // Parse FDR method
    let fdr_method = match args.fdr_method.to_lowercase().as_str() {
        "percolator" => Some(FdrMethod::Percolator),
        "mokapot" => Some(FdrMethod::Mokapot),
        "simple" => Some(FdrMethod::Simple),
        other => {
            log::warn!("Unknown FDR method '{}', defaulting to percolator", other);
            Some(FdrMethod::Percolator)
        }
    };

    // Parse FDR level
    let fdr_level = args
        .fdr_level
        .as_ref()
        .map(|s| match s.to_lowercase().as_str() {
            "precursor" => FdrLevel::Precursor,
            "peptide" => FdrLevel::Peptide,
            "protein" => FdrLevel::Protein,
            "both" => FdrLevel::Both,
            other => {
                log::warn!("Unknown FDR level '{}', defaulting to peptide", other);
                FdrLevel::Peptide
            }
        });

    // Parse shared peptide mode
    let shared_peptides = args
        .shared_peptides
        .as_ref()
        .map(|s| match s.to_lowercase().as_str() {
            "all" => SharedPeptideMode::All,
            "razor" => SharedPeptideMode::Razor,
            "unique" => SharedPeptideMode::Unique,
            other => {
                log::warn!(
                    "Unknown shared peptides mode '{}', defaulting to all",
                    other
                );
                SharedPeptideMode::All
            }
        });

    // Create overrides from CLI args (these take precedence over config file)
    let overrides = ConfigOverrides {
        input_files: args.input,
        library: args.library,
        output: args.output,
        report: args.report,
        rt_tolerance: Some(args.rt_tolerance),
        run_fdr: Some(args.run_fdr),
        experiment_fdr: Some(args.experiment_fdr),
        n_threads: args.threads,
        verbose: args.verbose,
        disable_rt_calibration: args.no_rt_calibration,
        fragment_tolerance: args.fragment_tolerance,
        fragment_unit,
        precursor_tolerance: args.precursor_tolerance,
        precursor_unit,
        fdr_method,
        fdr_level,
        protein_fdr: args.protein_fdr,
        shared_peptides,
        write_pin: args.write_pin,
    };

    // Apply CLI overrides
    config.merge_with_args(&overrides);

    // HPC scoring split: set the new config fields after merge_with_args.
    // These don't go through ConfigOverrides because they're CLI-only and
    // not commonly set in YAML (HPC orchestration is the use case).
    config.no_join = args.no_join;
    if let Some(ref s) = args.parquet_compression {
        config.parquet_compression = match s.to_lowercase().as_str() {
            "zstd" => ParquetCompression::Zstd,
            "snappy" => ParquetCompression::Snappy,
            other => {
                anyhow::bail!(
                    "Unknown --parquet-compression value '{}'. Expected 'zstd' or 'snappy'.",
                    other
                );
            }
        };
    }
    if args.join_only {
        let resolved = resolve_input_scores(args.input_scores.unwrap())?;
        log::info!(
            "--join-only: skipping Stages 1-4, loading {} `.scores.parquet` file(s)",
            resolved.len()
        );
        config.input_scores = Some(resolved);
    }

    // Override resolution mode if specified
    let resolution_mode = match args.resolution.to_lowercase().as_str() {
        "unit" => ResolutionMode::UnitResolution,
        "hram" => ResolutionMode::HRAM,
        _ => config.resolution_mode,
    };
    config.resolution_mode = resolution_mode;

    // Pre-filter: on by default for all modes; --no-prefilter disables it.
    if args.no_prefilter {
        config.prefilter_enabled = false;
    }

    // Apply unit resolution defaults (Th units, 1.0 Th precursor, 0.5 Th fragment)
    // Explicit CLI args override these defaults
    if resolution_mode == ResolutionMode::UnitResolution {
        // Precursor tolerance defaults for unit resolution
        if args.precursor_unit.is_none() {
            config.precursor_tolerance.unit = ToleranceUnit::Mz;
        }
        if args.precursor_tolerance.is_none() {
            config.precursor_tolerance.tolerance = 1.0;
        }

        // Fragment tolerance defaults for unit resolution
        if args.fragment_unit.is_none() {
            config.fragment_tolerance.unit = ToleranceUnit::Mz;
            if args.fragment_tolerance.is_none() {
                config.fragment_tolerance.tolerance = 0.5;
            }
        }

        log::info!(
            "Unit resolution mode: MS1 {:.2} {:?}, MS2 {:.2} {:?}",
            config.precursor_tolerance.tolerance,
            config.precursor_tolerance.unit,
            config.fragment_tolerance.tolerance,
            config.fragment_tolerance.unit
        );
    }

    // Set thread count
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    // Log reproducibility header
    log::info!("Osprey v{}", env!("CARGO_PKG_VERSION"));
    log::info!(
        "Command: {}",
        std::env::args().collect::<Vec<_>>().join(" ")
    );
    log::info!("Log file: {}", log_filename);
    match config.to_yaml_string() {
        Ok(yaml) => log::info!("Resolved configuration:\n{}", yaml),
        Err(e) => log::warn!("Could not serialize config: {}", e),
    }

    // Run analysis
    let start_time = Instant::now();

    run_analysis(config)?;

    let elapsed = start_time.elapsed();
    let elapsed_str = format_duration(elapsed);
    log::info!("Analysis complete in {}", elapsed_str);

    // Flush to ensure all log lines are written to the log file and terminal
    log::logger().flush();
    let _ = std::io::stderr().flush();

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::fs::File;
    use tempfile::TempDir;

    fn parse(extra: &[&str]) -> Args {
        let mut argv: Vec<&str> = vec!["osprey"];
        argv.extend_from_slice(extra);
        Args::try_parse_from(argv).unwrap()
    }

    fn assert_err_contains(result: Result<()>, needle: &str) {
        let err = result.expect_err("expected error").to_string();
        assert!(
            err.contains(needle),
            "expected error to contain {:?}, got: {}",
            needle,
            err
        );
    }

    // --- validate_hpc_args ----------------------------------------------

    #[test]
    fn validate_no_join_and_join_only_is_mutex() {
        let args = parse(&["--no-join", "--join-only", "-i", "x.mzML", "-l", "x.blib"]);
        assert_err_contains(validate_hpc_args(&args), "mutually exclusive");
    }

    #[test]
    fn validate_join_only_requires_input_scores() {
        let args = parse(&["--join-only", "-l", "x.blib", "-o", "y.blib"]);
        let err = validate_hpc_args(&args)
            .expect_err("expected error")
            .to_string();
        // Error refers to the canonical flag a user types to reach this branch.
        assert!(err.contains("--join-at-pass=1"), "got: {}", err);
        assert!(err.contains("--input-scores"), "got: {}", err);
    }

    #[test]
    fn validate_join_only_rejects_input_mzml() {
        let args = parse(&[
            "--join-only",
            "-i",
            "a.mzML",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "x.blib",
            "-o",
            "y.blib",
        ]);
        let err = validate_hpc_args(&args)
            .expect_err("expected error")
            .to_string();
        assert!(err.contains("--join-at-pass=1"), "got: {}", err);
        assert!(err.contains("--input"), "got: {}", err);
    }

    #[test]
    fn validate_join_only_requires_library_and_output() {
        let args = parse(&["--join-only", "--input-scores", "a.scores.parquet"]);
        let err = validate_hpc_args(&args)
            .expect_err("expected error")
            .to_string();
        assert!(err.contains("--join-at-pass=1"), "got: {}", err);
        assert!(err.contains("--library and --output"), "got: {}", err);
    }

    #[test]
    fn validate_no_join_requires_input() {
        let args = parse(&["--no-join", "-l", "x.blib"]);
        assert_err_contains(validate_hpc_args(&args), "--input");
    }

    #[test]
    fn validate_no_join_rejects_input_scores() {
        let args = parse(&[
            "--no-join",
            "-i",
            "a.mzML",
            "--input-scores",
            "a.scores.parquet",
        ]);
        assert_err_contains(validate_hpc_args(&args), "--input-scores");
    }

    #[test]
    fn validate_no_join_happy_path() {
        let args = parse(&["--no-join", "-i", "a.mzML", "-l", "ref.blib"]);
        validate_hpc_args(&args).unwrap();
    }

    #[test]
    fn validate_join_only_happy_path() {
        let args = parse(&[
            "--join-only",
            "--input-scores",
            "a.scores.parquet",
            "b.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        validate_hpc_args(&args).unwrap();
    }

    #[test]
    fn validate_default_mode_is_unaffected() {
        let args = parse(&["-i", "a.mzML", "-l", "ref.blib", "-o", "out.blib"]);
        validate_hpc_args(&args).unwrap();
    }

    // --- normalize_hpc_args (--join-at-pass) ----------------------------

    #[test]
    fn normalize_join_at_pass_1_sets_join_only() {
        let mut args = parse(&[
            "--join-at-pass=1",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert!(!args.join_only);
        normalize_hpc_args(&mut args).unwrap();
        assert!(args.join_only);
    }

    #[test]
    fn normalize_join_at_pass_2_errors_until_implemented() {
        let mut args = parse(&[
            "--join-at-pass=2",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert_err_contains(normalize_hpc_args(&mut args), "not yet implemented");
    }

    #[test]
    fn normalize_join_at_pass_invalid_value_errors() {
        let mut args = parse(&[
            "--join-at-pass=3",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert_err_contains(normalize_hpc_args(&mut args), "must be 1 or 2");
    }

    #[test]
    fn normalize_join_at_pass_1_with_join_only_modifier_errors_until_implemented() {
        // PR 2 will implement "run only Stage 5" via persisted plan files.
        let mut args = parse(&[
            "--join-at-pass=1",
            "--join-only",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert_err_contains(normalize_hpc_args(&mut args), "not yet implemented");
    }

    #[test]
    fn normalize_join_at_pass_1_with_no_join_modifier_errors_until_implemented() {
        // PR 2 will implement "run only Stage 6 from persisted Stage 5
        // outputs" once the boundary files exist.
        let mut args = parse(&[
            "--join-at-pass=1",
            "--no-join",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert_err_contains(normalize_hpc_args(&mut args), "not yet implemented");
    }

    #[test]
    fn normalize_join_only_alone_errors_no_entry_point() {
        let mut args = parse(&[
            "--join-only",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert_err_contains(normalize_hpc_args(&mut args), "modifier");
    }

    #[test]
    fn normalize_no_join_and_join_only_modifiers_are_mutex() {
        let mut args = parse(&[
            "--join-at-pass=1",
            "--no-join",
            "--join-only",
            "--input-scores",
            "a.scores.parquet",
            "-l",
            "ref.blib",
            "-o",
            "out.blib",
        ]);
        assert_err_contains(normalize_hpc_args(&mut args), "mutually exclusive");
    }

    #[test]
    fn normalize_no_join_alone_unchanged() {
        // The Stage 1 entry path with `-i ...` + `--no-join` keeps its
        // existing meaning: do per-file work only = Stages 1-4.
        let mut args = parse(&["--no-join", "-i", "a.mzML", "-l", "ref.blib"]);
        normalize_hpc_args(&mut args).unwrap();
        assert!(args.no_join);
        assert!(!args.join_only);
    }

    #[test]
    fn normalize_default_mode_is_a_noop() {
        let mut args = parse(&["-i", "a.mzML", "-l", "ref.blib", "-o", "out.blib"]);
        normalize_hpc_args(&mut args).unwrap();
        assert!(!args.join_only);
        assert!(args.join_at_pass.is_none());
    }

    // --- resolve_input_scores -------------------------------------------

    fn touch(path: &std::path::Path) {
        File::create(path).unwrap();
    }

    #[test]
    fn resolve_explicit_files_pass_through() {
        let dir = TempDir::new().unwrap();
        let a = dir.path().join("a.scores.parquet");
        let b = dir.path().join("b.scores.parquet");
        touch(&a);
        touch(&b);
        let resolved = resolve_input_scores(vec![a.clone(), b.clone()]).unwrap();
        assert_eq!(resolved, vec![a, b]);
    }

    #[test]
    fn resolve_explicit_missing_file_errors() {
        let dir = TempDir::new().unwrap();
        let missing = dir.path().join("does-not-exist.scores.parquet");
        let err = resolve_input_scores(vec![missing]).unwrap_err().to_string();
        assert!(err.contains("not found"), "got: {}", err);
    }

    #[test]
    fn resolve_directory_scans_and_sorts() {
        let dir = TempDir::new().unwrap();
        // Mix in a non-matching file to confirm the filter.
        touch(&dir.path().join("z.scores.parquet"));
        touch(&dir.path().join("a.scores.parquet"));
        touch(&dir.path().join("m.scores.parquet"));
        touch(&dir.path().join("readme.txt"));
        let resolved = resolve_input_scores(vec![dir.path().to_path_buf()]).unwrap();
        let names: Vec<String> = resolved
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
            .collect();
        assert_eq!(
            names,
            vec![
                "a.scores.parquet".to_string(),
                "m.scores.parquet".to_string(),
                "z.scores.parquet".to_string()
            ]
        );
    }

    #[test]
    fn resolve_empty_directory_errors() {
        let dir = TempDir::new().unwrap();
        touch(&dir.path().join("not-a-match.txt"));
        let err = resolve_input_scores(vec![dir.path().to_path_buf()])
            .unwrap_err()
            .to_string();
        assert!(err.contains("No *.scores.parquet"), "got: {}", err);
    }

    // --- OspreyConfig defaults ------------------------------------------

    #[test]
    fn config_defaults_disable_hpc_mode() {
        let cfg = OspreyConfig::default();
        assert!(!cfg.no_join, "no_join should default to false");
        assert!(
            cfg.input_scores.is_none(),
            "input_scores should default to None"
        );
    }
}
