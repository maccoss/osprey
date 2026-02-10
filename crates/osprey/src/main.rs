//! Osprey CLI - Command-line interface for peptide-centric DIA analysis

use anyhow::Result;
use clap::Parser;
use osprey::{run_analysis, ConfigOverrides, OspreyConfig, ResolutionMode, ToleranceUnit};
use std::io::Write;
use std::path::PathBuf;
use std::sync::Mutex;

/// Writer that tees output to both stderr and a log file
struct TeeWriter {
    stderr: std::io::Stderr,
    file: Mutex<std::io::BufWriter<std::fs::File>>,
}

impl Write for TeeWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.stderr.write_all(buf)?;
        self.file.lock().unwrap().write_all(buf)?;
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.stderr.flush()?;
        self.file.lock().unwrap().flush()?;
        Ok(())
    }
}

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

    /// Fixed regularization parameter lambda
    #[arg(long)]
    lambda: Option<f64>,

    /// Maximum candidates per spectrum (use 0 for unlimited)
    #[arg(long, default_value_t = 5250)]
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

    /// Export coefficient time series to parquet file(s)
    #[arg(long)]
    export_coefficients: bool,

    /// Use streaming mode for memory-efficient processing (requires streaming feature)
    #[arg(long)]
    streaming: bool,

    /// Verbose output
    #[arg(short, long)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Initialize logging — tee to both stderr and a log file in the current directory
    let log_level = if args.verbose { "debug" } else { "info" };
    let log_filename = format!("osprey_{}.log", chrono::Local::now().format("%Y-%m-%d_%H%M%S"));
    let log_file = std::fs::File::create(&log_filename)?;
    let tee = TeeWriter {
        stderr: std::io::stderr(),
        file: Mutex::new(std::io::BufWriter::new(log_file)),
    };
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(log_level))
        .target(env_logger::Target::Pipe(Box::new(tee)))
        .init();

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

    // Parse tolerance units if provided
    let fragment_unit = args.fragment_unit.as_ref().map(|s| {
        match s.to_lowercase().as_str() {
            "ppm" => ToleranceUnit::Ppm,
            "mz" | "th" | "da" => ToleranceUnit::Mz,
            _ => {
                log::warn!("Unknown fragment unit '{}', defaulting to ppm", s);
                ToleranceUnit::Ppm
            }
        }
    });

    let precursor_unit = args.precursor_unit.as_ref().map(|s| {
        match s.to_lowercase().as_str() {
            "ppm" => ToleranceUnit::Ppm,
            "mz" | "th" | "da" => ToleranceUnit::Mz,
            _ => {
                log::warn!("Unknown precursor unit '{}', defaulting to ppm", s);
                ToleranceUnit::Ppm
            }
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
        n_threads: args.threads,
        lambda: args.lambda,
        verbose: args.verbose,
        disable_rt_calibration: args.no_rt_calibration,
        streaming: args.streaming,
        fragment_tolerance: args.fragment_tolerance,
        fragment_unit,
        precursor_tolerance: args.precursor_tolerance,
        precursor_unit,
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

        log::info!("Unit resolution mode: MS1 {:.2} {:?}, MS2 {:.2} {:?}",
            config.precursor_tolerance.tolerance, config.precursor_tolerance.unit,
            config.fragment_tolerance.tolerance, config.fragment_tolerance.unit);
    }

    // Set max candidates
    config.max_candidates_per_spectrum = args.max_candidates;

    // Set export coefficients
    if args.export_coefficients {
        config.export_coefficients = true;
    }

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

    // Log reproducibility header
    log::info!("Osprey v{}", env!("CARGO_PKG_VERSION"));
    log::info!("Command: {}", std::env::args().collect::<Vec<_>>().join(" "));
    log::info!("Log file: {}", log_filename);
    match config.to_yaml_string() {
        Ok(yaml) => log::info!("Resolved configuration:\n{}", yaml),
        Err(e) => log::warn!("Could not serialize config: {}", e),
    }

    // Run analysis

    run_analysis(config)?;

    log::info!("Analysis complete.");

    Ok(())
}
