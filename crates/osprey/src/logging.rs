//! Two-tier logging: clean terminal output + verbose log file.
//!
//! By default, the terminal shows concise DIA-NN-style progress with elapsed
//! time prefixes. The log file always captures full detail (debug level) with
//! timestamps and module paths. The `--verbose` flag makes the terminal show
//! everything too.

use std::io::Write;
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

use indicatif::ProgressBar;
use log::{Level, LevelFilter, Log, Metadata, Record};

/// Global reference to the progress bar for coordinating with log output.
static PROGRESS_BAR: OnceLock<Mutex<Option<ProgressBar>>> = OnceLock::new();

fn progress_bar_slot() -> &'static Mutex<Option<ProgressBar>> {
    PROGRESS_BAR.get_or_init(|| Mutex::new(None))
}

/// Set the active progress bar for the logger.
///
/// While a progress bar is active, terminal log lines are routed through
/// `pb.println()` to prevent interleaving. Call `clear_progress_bar()`
/// when the progress bar is finished.
pub fn set_progress_bar(pb: &ProgressBar) {
    if let Ok(mut guard) = progress_bar_slot().lock() {
        *guard = Some(pb.clone());
    }
}

/// Clear the active progress bar from the logger.
pub fn clear_progress_bar() {
    if let Ok(mut guard) = progress_bar_slot().lock() {
        *guard = None;
    }
}

/// Logger that routes to two destinations with independent level filters.
///
/// - **Terminal** (stderr): shows only `terminal_level` and above, formatted
///   as `[M:SS] message` with no module path. Warn/Error get a level prefix.
/// - **File**: always captures `file_level` and above with full timestamps
///   and module paths, matching the old `env_logger` format.
pub struct TwoTierLogger {
    terminal_level: LevelFilter,
    file_level: LevelFilter,
    file: Mutex<std::io::BufWriter<std::fs::File>>,
    start_time: Instant,
}

impl TwoTierLogger {
    /// Create a new two-tier logger.
    ///
    /// `terminal_level`: minimum level for stderr output (Info or Debug).
    /// `file_level`: minimum level for log file output (typically Debug).
    pub fn new(terminal_level: LevelFilter, file_level: LevelFilter, file: std::fs::File) -> Self {
        Self {
            terminal_level,
            file_level,
            file: Mutex::new(std::io::BufWriter::new(file)),
            start_time: Instant::now(),
        }
    }

    /// Format elapsed time as `[M:SS]` or `[H:MM:SS]` for runs over an hour.
    fn format_elapsed(&self) -> String {
        let secs = self.start_time.elapsed().as_secs();
        let m = secs / 60;
        let s = secs % 60;
        if m >= 60 {
            let h = m / 60;
            format!("[{}:{:02}:{:02}]", h, m % 60, s)
        } else {
            format!("[{}:{:02}]", m, s)
        }
    }

    /// Format a log record for the terminal: `[M:SS] message`
    /// Warn/Error records get a level prefix: `[M:SS] WARN: message`
    fn format_terminal(&self, record: &Record) -> String {
        let elapsed = self.format_elapsed();
        match record.level() {
            Level::Error => format!("{} ERROR: {}", elapsed, record.args()),
            Level::Warn => format!("{} WARN: {}", elapsed, record.args()),
            _ => format!("{} {}", elapsed, record.args()),
        }
    }

    /// Format a log record for the log file with full detail:
    /// `[2026-04-09T10:30:15Z INFO  module::path] message`
    fn format_file(&self, record: &Record) -> String {
        let now = chrono::Utc::now().format("%Y-%m-%dT%H:%M:%SZ");
        let module = record.module_path().unwrap_or("unknown");
        format!(
            "[{} {:5} {}] {}",
            now,
            record.level(),
            module,
            record.args()
        )
    }
}

impl Log for TwoTierLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= self.terminal_level || metadata.level() <= self.file_level
    }

    fn log(&self, record: &Record) {
        if !self.enabled(record.metadata()) {
            return;
        }

        // Always write to file if level passes file filter
        if record.level() <= self.file_level {
            let line = self.format_file(record);
            if let Ok(mut f) = self.file.lock() {
                let _ = writeln!(f, "{}", line);
                let _ = f.flush();
            }
        }

        // Write to terminal if level passes terminal filter
        if record.level() <= self.terminal_level {
            let line = self.format_terminal(record);
            // Route through progress bar if active to avoid interleaving
            let pb = progress_bar_slot().lock().ok().and_then(|g| g.clone());
            if let Some(pb) = pb {
                pb.println(&line);
            } else {
                let mut stderr = std::io::stderr().lock();
                let _ = writeln!(stderr, "{}", line);
            }
        }
    }

    fn flush(&self) {
        if let Ok(mut f) = self.file.lock() {
            let _ = f.flush();
        }
        let _ = std::io::stderr().flush();
    }
}

/// Install the two-tier logger as the global logger.
///
/// `terminal_level`: `LevelFilter::Info` for clean output, `LevelFilter::Debug` for verbose.
/// `file`: the log file to write full detail to.
pub fn init(terminal_level: LevelFilter, file: std::fs::File) {
    let logger = TwoTierLogger::new(terminal_level, LevelFilter::Debug, file);
    log::set_boxed_logger(Box::new(logger)).expect("Failed to set logger");
    log::set_max_level(LevelFilter::Debug);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_elapsed_under_one_hour() {
        let logger = TwoTierLogger::new(
            LevelFilter::Info,
            LevelFilter::Debug,
            tempfile::tempfile().unwrap(),
        );
        let elapsed = logger.format_elapsed();
        assert!(elapsed.starts_with('['));
        assert!(elapsed.ends_with(']'));
        assert!(elapsed.contains(':'));
    }

    #[test]
    fn test_format_terminal_info_no_level_prefix() {
        let logger = TwoTierLogger::new(
            LevelFilter::Info,
            LevelFilter::Debug,
            tempfile::tempfile().unwrap(),
        );
        let record = log::Record::builder()
            .args(format_args!("Loading library"))
            .level(Level::Info)
            .target("test")
            .build();
        let formatted = logger.format_terminal(&record);
        // Should have elapsed prefix but no "INFO:" tag
        assert!(formatted.contains("Loading library"));
        assert!(!formatted.contains("INFO"));
    }

    #[test]
    fn test_format_terminal_warn_has_prefix() {
        let logger = TwoTierLogger::new(
            LevelFilter::Info,
            LevelFilter::Debug,
            tempfile::tempfile().unwrap(),
        );
        let record = log::Record::builder()
            .args(format_args!("Cache stale"))
            .level(Level::Warn)
            .target("test")
            .build();
        let formatted = logger.format_terminal(&record);
        assert!(formatted.contains("WARN: Cache stale"));
    }

    #[test]
    fn test_format_file_has_full_detail() {
        let logger = TwoTierLogger::new(
            LevelFilter::Info,
            LevelFilter::Debug,
            tempfile::tempfile().unwrap(),
        );
        let record = log::Record::builder()
            .args(format_args!("test message"))
            .level(Level::Debug)
            .target("test")
            .module_path(Some("osprey::pipeline"))
            .build();
        let formatted = logger.format_file(&record);
        assert!(formatted.contains("DEBUG"));
        assert!(formatted.contains("osprey::pipeline"));
        assert!(formatted.contains("test message"));
    }
}
