//! Osprey CLI - Command-line interface for peptide-centric DIA analysis

use anyhow::Result;
use clap::Parser;
use osprey::{run_analysis, ConfigOverrides, OspreyConfig, ResolutionMode};
use std::path::PathBuf;

/// Osprey: Peptide-centric DIA analysis tool
#[derive(Parser, Debug)]
#[command(name = "osprey")]
#[command(version, about = "Peptide-centric DIA analysis with Skyline integration")]
#[command(long_about = r#"
Osprey is an open-source tool for peptide detection and quantification
in data-independent acquisition (DIA) mass spectrometry data.

It uses ridge regression to deconvolute mixed MS/MS spectra, aggregates
evidence across the chromatographic dimension, and uses machine learning
to score peptide detections with rigorous FDR control.

EXAMPLES:
    # Basic analysis with DIA-NN library
    osprey -i sample.mzML -l library.tsv -o results.blib

    # Multiple input files
    osprey -i *.mzML -l library.tsv -o results.blib

    # Unit resolution mode with specific lambda
    osprey -i sample.mzML -l library.tsv -o results.blib --resolution unit --lambda 0.5

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

    /// RT tolerance in minutes (fallback when calibration disabled)
    #[arg(long, default_value_t = 2.0)]
    rt_tolerance: f64,

    /// Disable RT calibration (use fixed rt_tolerance instead)
    #[arg(long)]
    no_rt_calibration: bool,

    /// Fixed regularization parameter lambda
    #[arg(long)]
    lambda: Option<f64>,

    /// Maximum candidates per spectrum
    #[arg(long, default_value_t = 200)]
    max_candidates: usize,

    /// Run-level FDR threshold
    #[arg(long, default_value_t = 0.01)]
    run_fdr: f64,

    /// Number of threads (default: all available)
    #[arg(long)]
    threads: Option<usize>,

    /// Write TSV report to this file
    #[arg(long)]
    report: Option<PathBuf>,

    /// Use streaming mode for memory-efficient processing (requires streaming feature)
    #[arg(long)]
    streaming: bool,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Initialize logging
    let log_level = if args.verbose { "debug" } else { "info" };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level)).init();

    // Handle --generate-config
    if let Some(config_path) = args.generate_config {
        log::info!("Generating template configuration file: {}", config_path.display());
        OspreyConfig::create_template(&config_path)?;
        log::info!("Template created. Edit the file and run with --config {}", config_path.display());
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

    // Create overrides from CLI args (these take precedence over config file)
    let overrides = ConfigOverrides {
        input_files: args.input,
        library: args.library,
        output: args.output,
        report: args.report,
        rt_tolerance: Some(args.rt_tolerance),
        run_fdr: Some(args.run_fdr),
        n_threads: args.threads,
        lambda: args.lambda,
        verbose: args.verbose,
        disable_rt_calibration: args.no_rt_calibration,
        streaming: args.streaming,
    };

    // Apply CLI overrides
    config.merge_with_args(&overrides);

    // Override resolution mode if specified
    let resolution_mode = match args.resolution.to_lowercase().as_str() {
        "unit" => ResolutionMode::UnitResolution,
        "hram" => ResolutionMode::HRAM,
        "auto" | _ => config.resolution_mode,
    };
    config.resolution_mode = resolution_mode;

    // Set max candidates
    config.max_candidates_per_spectrum = args.max_candidates;

    // Set thread count
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()?;
    }

    // Check streaming mode
    #[cfg(feature = "streaming")]
    if config.streaming {
        log::info!("Streaming mode enabled");
    }

    #[cfg(not(feature = "streaming"))]
    if args.streaming {
        log::warn!("Streaming mode requested but 'streaming' feature not compiled. Ignoring.");
        // Don't set config.streaming since we already merged args
    }

    // Run analysis
    log::info!("Osprey v{}", env!("CARGO_PKG_VERSION"));
    log::info!("Input files: {:?}", config.input_files);
    log::info!("Library: {:?}", config.library_source.path());
    log::info!("Output: {:?}", config.output_blib);

    let results = run_analysis(config)?;

    log::info!("Analysis complete.");
    log::info!("Total regression results: {}", results.len());

    // Summary statistics
    let total_coefficients: usize = results.iter().map(|r| r.coefficients.len()).sum();
    log::info!("Total non-zero coefficients: {}", total_coefficients);

    Ok(())
}
