//! Osprey CLI - Command-line interface for peptide-centric DIA analysis

use anyhow::Result;
use clap::Parser;
use log::LevelFilter;
use osprey::{
    run_analysis, ConfigOverrides, FdrLevel, FdrMethod, OspreyConfig, ResolutionMode,
    SharedPeptideMode, ToleranceUnit,
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

    /// FDR filtering level: precursor, peptide, or both (default)
    #[arg(long)]
    fdr_level: Option<String>,

    /// Protein-level FDR threshold (enables protein parsimony and picked-protein FDR)
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
}

fn main() -> Result<()> {
    let args = Args::parse();

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

    // Load config from file or build from CLI args
    let mut config = if let Some(config_path) = &args.config {
        log::info!("Loading configuration from: {}", config_path.display());
        OspreyConfig::from_yaml(config_path)?
    } else {
        // Require input, library, output if no config file
        if args.input.is_none() || args.library.is_none() || args.output.is_none() {
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
            "both" => FdrLevel::Both,
            other => {
                log::warn!("Unknown FDR level '{}', defaulting to both", other);
                FdrLevel::Both
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

    // Flush to ensure the final log line appears on the terminal before the prompt
    let _ = std::io::stderr().flush();

    Ok(())
}
