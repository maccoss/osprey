//! Cross-implementation bisection diagnostic primitives.
//!
//! Shared helpers for the diagnostic dump infrastructure used to compare
//! Osprey (this crate) and OspreySharp (the C# port at
//! `pwiz_tools/OspreySharp/`). Each pipeline checkpoint that participates
//! in cross-impl parity testing dumps state via env-var-gated functions in
//! `osprey-scoring::diagnostics` and `osprey::diagnostics`. This module
//! holds the primitives those functions share.
//!
//! ## Filename convention
//!
//! Rust dumps begin with `rust_`; C# OspreySharp dumps begin with `cs_`.
//! Diff `cs_X.txt` against `rust_X.txt` to locate divergence at the
//! corresponding stage.
//!
//! ## Numeric formatting
//!
//! Most dumped floating-point values use F10 (1e-10) fixed precision via
//! [`format_f10`] or an inline `{:.10}` format spec. A few dumps that need
//! full-bit fidelity (e.g., LOESS input pairs) intentionally use `{:.17}`
//! or full-precision default formatting. Rust's `{:.N}` format spec rounds
//! half-to-even (banker's). The C# side mirrors this in
//! `pwiz.OspreySharp.OspreyDiagnostics.F10()` (which forces half-to-even
//! via `Math.Round(MidpointRounding.ToEven)`) so the formatted strings
//! diff cleanly.

/// Returns true if the named env var is set (any value).
///
/// Used to gate dump emission. Most diagnostic functions take this form:
///
/// ```ignore
/// pub fn dump_X(...) {
///     if !diagnostics::is_dump_enabled("OSPREY_DUMP_X") { return; }
///     // ... write rust_X.txt ...
///     diagnostics::exit_if_only("OSPREY_X_ONLY", "X dump");
/// }
/// ```
pub fn is_dump_enabled(name: &str) -> bool {
    std::env::var(name).is_ok()
}

/// If the named env var is set, log an aborting message and exit the
/// process. Used for `_ONLY` variants of dump env vars to short-circuit
/// the pipeline after a checkpoint dump for fast bisection cycles.
pub fn exit_if_only(name: &str, label: &str) {
    if std::env::var(name).is_ok() {
        log::info!("{} set - aborting after {}", name, label);
        std::process::exit(0);
    }
}

/// Format a `f64` at F10 (1e-10) fixed precision matching the C#
/// `OspreyDiagnostics.F10()` output. Rust's `{:.10}` already rounds
/// half-to-even (banker's), so this is a thin wrapper for symmetry with
/// the C# API and explicit intent at call sites.
pub fn format_f10(x: f64) -> String {
    format!("{:.10}", x)
}

/// If `OSPREY_EXIT_AFTER_CALIBRATION=1` is set, log a benchmark message
/// and return true so the caller can short-circuit the pipeline after
/// Stage 3. Used to time and diff calibration in isolation.
pub fn should_exit_after_calibration() -> bool {
    if std::env::var("OSPREY_EXIT_AFTER_CALIBRATION").is_ok() {
        log::info!(
            "[BENCH] OSPREY_EXIT_AFTER_CALIBRATION set - exiting after Stage 3 (calibration done)"
        );
        true
    } else {
        false
    }
}

// Note: the `OSPREY_EXIT_AFTER_SCORING` env-var gate that used to live
// here was retired in favor of the `--no-join` CLI flag. See the HPC
// scoring split work in `crates/osprey/src/pipeline.rs::run_analysis`.
// `OSPREY_EXIT_AFTER_CALIBRATION` (Stage 3) stays because it has no
// production CLI analog -- its purpose is purely diagnostic.
