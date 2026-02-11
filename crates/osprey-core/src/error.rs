//! Error types for Osprey
//!
//! This module provides a unified error type for all Osprey operations,
//! with specific variants for different failure modes.

use thiserror::Error;

/// Unified error type for Osprey operations
#[derive(Error, Debug)]
pub enum OspreyError {
    /// Failed to parse mzML file
    #[error("Failed to parse mzML file: {0}")]
    MzmlParseError(String),

    /// Failed to load spectral library
    #[error("Failed to load library: {0}")]
    LibraryLoadError(String),

    /// Invalid library format
    #[error("Invalid library format: {0}")]
    InvalidLibraryFormat(String),

    /// Library entry is invalid
    #[error("Invalid library entry: {0}")]
    InvalidLibraryEntry(String),

    /// Regression failed
    #[error("Regression failed: {0}")]
    RegressionError(String),

    /// Matrix operation failed
    #[error("Matrix operation failed: {0}")]
    MatrixError(String),

    /// Peak detection failed
    #[error("Peak detection failed: {0}")]
    PeakDetectionError(String),

    /// FDR calculation failed
    #[error("FDR calculation failed: {0}")]
    FdrError(String),

    /// Configuration error
    #[error("Configuration error: {0}")]
    ConfigError(String),

    /// Invalid input file
    #[error("Invalid input file '{path}': {reason}")]
    InvalidInputFile { path: String, reason: String },

    /// File not found
    #[error("File not found: {0}")]
    FileNotFound(String),

    /// I/O error
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),

    /// CSV parsing error
    #[error("CSV parsing error: {0}")]
    CsvError(String),

    /// SQLite error (for blib/elib files)
    #[error("SQLite error: {0}")]
    SqliteError(String),

    /// XML parsing error
    #[error("XML parsing error: {0}")]
    XmlError(String),

    /// Numeric conversion error
    #[error("Numeric error: {0}")]
    NumericError(String),

    /// Thread pool error
    #[error("Thread pool error: {0}")]
    ThreadPoolError(String),

    /// External tool error (mokapot, etc.)
    #[error("External tool error: {0}")]
    ExternalToolError(String),

    /// Feature extraction error
    #[error("Feature extraction error: {0}")]
    FeatureError(String),

    /// Output writing error
    #[error("Output writing error: {0}")]
    OutputError(String),

    /// Generic internal error
    #[error("Internal error: {0}")]
    InternalError(String),
}

/// Result type alias for Osprey operations
pub type Result<T> = std::result::Result<T, OspreyError>;

impl OspreyError {
    /// Create a library load error
    pub fn library_load(msg: impl Into<String>) -> Self {
        OspreyError::LibraryLoadError(msg.into())
    }

    /// Create a regression error
    pub fn regression(msg: impl Into<String>) -> Self {
        OspreyError::RegressionError(msg.into())
    }

    /// Create a configuration error
    pub fn config(msg: impl Into<String>) -> Self {
        OspreyError::ConfigError(msg.into())
    }

    /// Create an input file error
    pub fn invalid_input(path: impl Into<String>, reason: impl Into<String>) -> Self {
        OspreyError::InvalidInputFile {
            path: path.into(),
            reason: reason.into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Verifies that OspreyError Display impl formats the error message correctly.
    #[test]
    fn test_error_display() {
        let err = OspreyError::LibraryLoadError("test error".to_string());
        assert_eq!(format!("{}", err), "Failed to load library: test error");
    }

    /// Verifies that helper constructors produce the correct OspreyError variants.
    #[test]
    fn test_error_helpers() {
        let err = OspreyError::library_load("test");
        assert!(matches!(err, OspreyError::LibraryLoadError(_)));

        let err = OspreyError::config("bad config");
        assert!(matches!(err, OspreyError::ConfigError(_)));
    }
}
