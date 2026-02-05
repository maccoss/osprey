//! TSV report writer for Osprey results
//!
//! Writes human-readable tab-separated reports of analysis results.

use osprey_core::{OspreyError, RegressionResult, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Writer for TSV coefficient reports
pub struct ReportWriter {
    writer: BufWriter<File>,
    header_written: bool,
}

impl ReportWriter {
    /// Create a new report writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path.as_ref()).map_err(|e| {
            OspreyError::OutputError(format!(
                "Failed to create output file '{}': {}",
                path.as_ref().display(),
                e
            ))
        })?;

        Ok(Self {
            writer: BufWriter::new(file),
            header_written: false,
        })
    }

    /// Write header row
    fn write_header(&mut self) -> Result<()> {
        if !self.header_written {
            writeln!(
                self.writer,
                "scan_number\tretention_time\tlibrary_id\tcoefficient"
            )
            .map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;
            self.header_written = true;
        }
        Ok(())
    }

    /// Write a regression result
    pub fn write_result(&mut self, result: &RegressionResult) -> Result<()> {
        self.write_header()?;

        for (library_id, coefficient) in result.library_ids.iter().zip(result.coefficients.iter()) {
            if *coefficient > 0.0 {
                writeln!(
                    self.writer,
                    "{}\t{:.4}\t{}\t{:.6}",
                    result.scan_number, result.retention_time, library_id, coefficient
                )
                .map_err(|e| OspreyError::OutputError(format!("Failed to write row: {}", e)))?;
            }
        }

        Ok(())
    }

    /// Flush and close the writer
    pub fn finish(mut self) -> Result<()> {
        self.writer
            .flush()
            .map_err(|e| OspreyError::OutputError(format!("Failed to flush output: {}", e)))?;
        Ok(())
    }
}

/// Write a simple coefficient summary to a file
pub fn write_coefficient_summary<P: AsRef<Path>>(
    path: P,
    results: &[RegressionResult],
    library_sequences: &[String],
) -> Result<()> {
    let file = File::create(path.as_ref()).map_err(|e| {
        OspreyError::OutputError(format!(
            "Failed to create output file '{}': {}",
            path.as_ref().display(),
            e
        ))
    })?;

    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(
        writer,
        "scan_number\tretention_time\tlibrary_id\tmodified_sequence\tcoefficient"
    )
    .map_err(|e| OspreyError::OutputError(format!("Failed to write header: {}", e)))?;

    // Write results
    for result in results {
        for (library_id, coefficient) in result.library_ids.iter().zip(result.coefficients.iter()) {
            if *coefficient > 0.0 {
                let seq = library_sequences
                    .get(*library_id as usize)
                    .map(|s| s.as_str())
                    .unwrap_or("");

                writeln!(
                    writer,
                    "{}\t{:.4}\t{}\t{}\t{:.6}",
                    result.scan_number, result.retention_time, library_id, seq, coefficient
                )
                .map_err(|e| OspreyError::OutputError(format!("Failed to write row: {}", e)))?;
            }
        }
    }

    writer
        .flush()
        .map_err(|e| OspreyError::OutputError(format!("Failed to flush output: {}", e)))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_write_result() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.tsv");

        let mut writer = ReportWriter::new(&path).unwrap();

        let result = RegressionResult {
            scan_number: 100,
            retention_time: 10.5,
            library_ids: vec![1, 2, 3],
            coefficients: vec![0.5, 0.0, 0.3],
            residual: 0.1,
            n_candidates: 5,
            coefficient_sum: 0.8,
            observed_norm: 1.0,
        };

        writer.write_result(&result).unwrap();
        writer.finish().unwrap();

        // Verify file was created
        assert!(path.exists());
    }
}
