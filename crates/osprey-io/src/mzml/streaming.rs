//! Streaming mzML parser using async I/O
//!
//! This module provides an async streaming mzML parser adapted from Sage
//! (https://github.com/lazear/sage). Sage is Copyright (c) 2022 Michael Lazear,
//! licensed under the MIT License. See NOTICE.md in the repository root for details.
//!
//! It uses event-based XML parsing for memory efficiency and can send
//! spectra through channels as they are parsed.
//!
//! ## Features
//!
//! - Async I/O with tokio
//! - Event-based streaming (O(1) memory for parsing)
//! - On-the-fly zlib decompression
//! - Keeps f64 precision for m/z values (for calibration accuracy)
//!
//! ## Usage
//!
//! ```ignore
//! use osprey_io::mzml::streaming::StreamingMzmlReader;
//! use crossbeam::channel;
//!
//! let (tx, rx) = channel::bounded(1000);
//! let reader = StreamingMzmlReader::new();
//! reader.parse_file_to_channel("file.mzML", tx).await?;
//! ```

use async_compression::tokio::bufread::ZlibDecoder;
use crossbeam::channel::Sender;
use osprey_core::{IsolationWindow, Spectrum};
use quick_xml::events::Event;
use quick_xml::Reader;
use std::path::Path;
use tokio::fs::File;
use tokio::io::{AsyncBufRead, AsyncReadExt, BufReader};

// CV parameter constants (MS ontology accession numbers)
// Compression types
const ZLIB_COMPRESSION: &[u8] = b"MS:1000574";
const NO_COMPRESSION: &[u8] = b"MS:1000576";

// Array types
const INTENSITY_ARRAY: &[u8] = b"MS:1000515";
const MZ_ARRAY: &[u8] = b"MS:1000514";

// Data types
const FLOAT_64: &[u8] = b"MS:1000523";
const FLOAT_32: &[u8] = b"MS:1000521";

// Spectrum metadata
const MS_LEVEL: &[u8] = b"MS:1000511";
const TOTAL_ION_CURRENT: &[u8] = b"MS:1000285";

// Scan metadata
const SCAN_START_TIME: &[u8] = b"MS:1000016";
const UNIT_SECONDS: &[u8] = b"UO:0000010";
const UNIT_MINUTES: &[u8] = b"UO:0000031";

// Precursor/isolation window
const SELECTED_ION_MZ: &[u8] = b"MS:1000744";
const ISO_WINDOW_TARGET: &[u8] = b"MS:1000827";
const ISO_WINDOW_LOWER: &[u8] = b"MS:1000828";
const ISO_WINDOW_UPPER: &[u8] = b"MS:1000829";

/// State machine for tracking position in XML hierarchy
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum State {
    Spectrum,
    Scan,
    BinaryDataArray,
    Binary,
    Precursor,
    SelectedIon,
}

/// Type of binary data array
#[derive(Copy, Clone, Debug)]
enum BinaryKind {
    Intensity,
    Mz,
}

/// Data type of binary array
#[derive(Copy, Clone, Debug)]
enum Dtype {
    F32,
    F64,
}

/// Statistics from parsing
#[derive(Debug, Clone, Default)]
pub struct ParseStats {
    pub total_spectra: usize,
    pub ms1_spectra: usize,
    pub ms2_spectra: usize,
    pub skipped_spectra: usize,
}

/// Streaming mzML reader using async I/O
///
/// Adapted from Sage's mzML parser with modifications for Osprey:
/// - Keeps f64 precision for m/z values (calibration accuracy)
/// - Outputs Osprey's Spectrum type
/// - Sends spectra through channels instead of collecting
pub struct StreamingMzmlReader {
    /// Filter to specific MS level (None = all levels)
    ms_level_filter: Option<u8>,
}

impl Default for StreamingMzmlReader {
    fn default() -> Self {
        Self::new()
    }
}

impl StreamingMzmlReader {
    /// Create a new streaming reader (no filtering)
    pub fn new() -> Self {
        Self {
            ms_level_filter: None,
        }
    }

    /// Create a reader that filters to a specific MS level
    pub fn with_ms_level_filter(ms_level: u8) -> Self {
        Self {
            ms_level_filter: Some(ms_level),
        }
    }

    /// Create a reader for MS2 spectra only
    pub fn ms2_only() -> Self {
        Self::with_ms_level_filter(2)
    }

    /// Parse an mzML file and send spectra to a channel
    pub async fn parse_file_to_channel(
        &self,
        path: impl AsRef<Path>,
        sender: Sender<Result<Spectrum, StreamingMzmlError>>,
    ) -> Result<ParseStats, StreamingMzmlError> {
        let file = File::open(path).await?;
        let reader = BufReader::new(file);
        self.parse_to_channel(reader, sender).await
    }

    /// Parse from any async buffered reader and send spectra to a channel
    pub async fn parse_to_channel<R: AsyncBufRead + Unpin>(
        &self,
        reader: R,
        sender: Sender<Result<Spectrum, StreamingMzmlError>>,
    ) -> Result<ParseStats, StreamingMzmlError> {
        let mut xml_reader = Reader::from_reader(reader);
        let mut buf = Vec::new();
        let mut output_buffer = Vec::with_capacity(4096);

        // Parsing state
        let mut state: Option<State> = None;
        let mut compression = false;
        let mut binary_dtype = Dtype::F64;
        let mut binary_array: Option<BinaryKind> = None;

        // Current spectrum being built
        let mut scan_number: u32 = 0;
        let mut ms_level: u8 = 0;
        let mut retention_time: f64 = 0.0;
        let mut precursor_mz: f64 = 0.0;
        let mut iso_window_lo: Option<f64> = None;
        let mut iso_window_hi: Option<f64> = None;
        let mut mzs: Vec<f64> = Vec::new();
        let mut intensities: Vec<f32> = Vec::new();

        // Statistics
        let mut stats = ParseStats::default();

        // Macro to extract attribute value
        macro_rules! extract {
            ($ev:expr, $key:expr) => {
                $ev.try_get_attribute($key)?
                    .ok_or(StreamingMzmlError::MalformedXml)?
                    .value
            };
        }

        macro_rules! extract_value {
            ($ev:expr) => {{
                let s = $ev
                    .try_get_attribute(b"value")?
                    .ok_or(StreamingMzmlError::MalformedXml)?
                    .value;
                std::str::from_utf8(&s)?.parse()?
            }};
        }

        loop {
            match xml_reader.read_event_into_async(&mut buf).await {
                Ok(Event::Start(ref ev)) => {
                    // State transition into child tag
                    state = match (ev.name().into_inner(), state) {
                        (b"spectrum", _) => Some(State::Spectrum),
                        (b"scan", Some(State::Spectrum)) => Some(State::Scan),
                        (b"binaryDataArray", Some(State::Spectrum)) => Some(State::BinaryDataArray),
                        (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
                        (b"precursor", Some(State::Spectrum)) => Some(State::Precursor),
                        (b"selectedIon", Some(State::Precursor)) => Some(State::SelectedIon),
                        _ => state,
                    };

                    // Extract spectrum ID for scan number
                    if ev.name().into_inner() == b"spectrum" {
                        let id = extract!(ev, b"id");
                        let id_str = std::str::from_utf8(&id)?;
                        scan_number = parse_scan_number(id_str);
                    }
                }
                Ok(Event::Empty(ref ev)) => match (state, ev.name().into_inner()) {
                    (Some(State::BinaryDataArray), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ZLIB_COMPRESSION => compression = true,
                            NO_COMPRESSION => compression = false,
                            FLOAT_64 => binary_dtype = Dtype::F64,
                            FLOAT_32 => binary_dtype = Dtype::F32,
                            INTENSITY_ARRAY => binary_array = Some(BinaryKind::Intensity),
                            MZ_ARRAY => binary_array = Some(BinaryKind::Mz),
                            _ => {}
                        }
                    }
                    (Some(State::Spectrum), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            MS_LEVEL => {
                                let level: u8 = extract_value!(ev);
                                ms_level = level;
                                // Early skip if filtering by MS level
                                if let Some(filter) = self.ms_level_filter {
                                    if level != filter {
                                        state = None;
                                    }
                                }
                            }
                            TOTAL_ION_CURRENT => {
                                let tic: f64 = extract_value!(ev);
                                if tic == 0.0 {
                                    // No ion current, skip this spectrum
                                    state = None;
                                }
                            }
                            _ => {}
                        }
                    }
                    (Some(State::Precursor), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        match accession.as_ref() {
                            ISO_WINDOW_TARGET => {
                                if precursor_mz == 0.0 {
                                    precursor_mz = extract_value!(ev);
                                }
                            }
                            ISO_WINDOW_LOWER => iso_window_lo = Some(extract_value!(ev)),
                            ISO_WINDOW_UPPER => iso_window_hi = Some(extract_value!(ev)),
                            _ => {}
                        }
                    }
                    (Some(State::SelectedIon), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        if accession.as_ref() == SELECTED_ION_MZ && precursor_mz == 0.0 {
                            precursor_mz = extract_value!(ev);
                        }
                    }
                    (Some(State::Scan), b"cvParam") => {
                        let accession = extract!(ev, b"accession");
                        if accession.as_ref() == SCAN_START_TIME {
                            let time: f64 = extract_value!(ev);
                            let unit = extract!(ev, b"unitAccession");
                            retention_time = match unit.as_ref() {
                                UNIT_SECONDS => time / 60.0,
                                UNIT_MINUTES => time,
                                _ => time, // Assume minutes if unknown
                            };
                        }
                    }
                    _ => {}
                },
                Ok(Event::Text(text)) => {
                    if let Some(State::Binary) = state {
                        // Skip if we filtered out this spectrum
                        if let Some(filter) = self.ms_level_filter {
                            if ms_level != filter {
                                continue;
                            }
                        }

                        let raw = text.unescape()?;
                        if raw.is_empty() || binary_array.is_none() {
                            continue;
                        }

                        // Decode base64
                        use base64::Engine;
                        let decoded = base64::engine::general_purpose::STANDARD
                            .decode(raw.as_bytes())?;

                        // Decompress if needed
                        let bytes = match compression {
                            false => &decoded[..],
                            true => {
                                let mut decoder = ZlibDecoder::new(decoded.as_slice());
                                let n = decoder.read_to_end(&mut output_buffer).await?;
                                &output_buffer[..n]
                            }
                        };

                        // Parse binary data - KEEP f64 PRECISION for m/z
                        match binary_array {
                            Some(BinaryKind::Mz) => {
                                mzs = match binary_dtype {
                                    Dtype::F64 => {
                                        let mut buf: [u8; 8] = [0; 8];
                                        bytes
                                            .chunks(8)
                                            .filter(|chunk| chunk.len() == 8)
                                            .map(|chunk| {
                                                buf.copy_from_slice(chunk);
                                                f64::from_le_bytes(buf)
                                            })
                                            .collect()
                                    }
                                    Dtype::F32 => {
                                        let mut buf: [u8; 4] = [0; 4];
                                        bytes
                                            .chunks(4)
                                            .filter(|chunk| chunk.len() == 4)
                                            .map(|chunk| {
                                                buf.copy_from_slice(chunk);
                                                f32::from_le_bytes(buf) as f64
                                            })
                                            .collect()
                                    }
                                };
                            }
                            Some(BinaryKind::Intensity) => {
                                intensities = match binary_dtype {
                                    Dtype::F32 => {
                                        let mut buf: [u8; 4] = [0; 4];
                                        bytes
                                            .chunks(4)
                                            .filter(|chunk| chunk.len() == 4)
                                            .map(|chunk| {
                                                buf.copy_from_slice(chunk);
                                                f32::from_le_bytes(buf)
                                            })
                                            .collect()
                                    }
                                    Dtype::F64 => {
                                        let mut buf: [u8; 8] = [0; 8];
                                        bytes
                                            .chunks(8)
                                            .filter(|chunk| chunk.len() == 8)
                                            .map(|chunk| {
                                                buf.copy_from_slice(chunk);
                                                f64::from_le_bytes(buf) as f32
                                            })
                                            .collect()
                                    }
                                };
                            }
                            None => {}
                        }

                        output_buffer.clear();
                        binary_array = None;
                    }
                }
                Ok(Event::End(ev)) => {
                    state = match (state, ev.name().into_inner()) {
                        (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
                        (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
                        (Some(State::SelectedIon), b"selectedIon") => Some(State::Precursor),
                        (Some(State::Precursor), b"precursor") => Some(State::Spectrum),
                        (Some(State::Scan), b"scan") => Some(State::Spectrum),
                        (_, b"spectrum") => {
                            stats.total_spectra += 1;

                            // Check if we should emit this spectrum
                            let should_emit = self
                                .ms_level_filter
                                .map(|filter| filter == ms_level)
                                .unwrap_or(true);

                            if should_emit && !mzs.is_empty() && precursor_mz > 0.0 {
                                // Build isolation window
                                let isolation_window = match (iso_window_lo, iso_window_hi) {
                                    (Some(lo), Some(hi)) => {
                                        IsolationWindow::new(precursor_mz, lo, hi)
                                    }
                                    _ => IsolationWindow::symmetric(precursor_mz, 12.5),
                                };

                                let spectrum = Spectrum {
                                    scan_number,
                                    retention_time,
                                    precursor_mz,
                                    isolation_window,
                                    mzs: std::mem::take(&mut mzs),
                                    intensities: std::mem::take(&mut intensities),
                                };

                                // Track MS level stats
                                match ms_level {
                                    1 => stats.ms1_spectra += 1,
                                    2 => stats.ms2_spectra += 1,
                                    _ => {}
                                }

                                // Send to channel (blocks if full - backpressure)
                                if sender.send(Ok(spectrum)).is_err() {
                                    // Receiver dropped, stop parsing
                                    break;
                                }
                            } else {
                                stats.skipped_spectra += 1;
                            }

                            // Reset for next spectrum
                            scan_number = 0;
                            ms_level = 0;
                            retention_time = 0.0;
                            precursor_mz = 0.0;
                            iso_window_lo = None;
                            iso_window_hi = None;
                            mzs.clear();
                            intensities.clear();

                            None
                        }
                        _ => state,
                    };
                }
                Ok(Event::Eof) => break,
                Ok(_) => {}
                Err(err) => {
                    log::error!("XML parsing error: {}", err);
                    let _ = sender.send(Err(StreamingMzmlError::XmlError(err)));
                }
            }
            buf.clear();
        }

        Ok(stats)
    }

    /// Parse mzML file and collect all spectra (convenience method for testing)
    pub async fn parse_file(
        &self,
        path: impl AsRef<Path>,
    ) -> Result<Vec<Spectrum>, StreamingMzmlError> {
        let (tx, rx) = crossbeam::channel::unbounded();
        let stats = self.parse_file_to_channel(path, tx).await?;
        log::debug!(
            "Parsed {} spectra (MS1: {}, MS2: {}, skipped: {})",
            stats.total_spectra,
            stats.ms1_spectra,
            stats.ms2_spectra,
            stats.skipped_spectra
        );
        Ok(rx.into_iter().filter_map(|r| r.ok()).collect())
    }
}

/// Parse scan number from spectrum ID string
///
/// Handles common formats:
/// - "scan=123" -> 123
/// - "spectrum=123" -> 123
/// - "index=123" -> 123
/// - Just extract any number found
fn parse_scan_number(id: &str) -> u32 {
    // Try common patterns
    for prefix in &["scan=", "spectrum=", "index="] {
        if let Some(rest) = id.strip_prefix(prefix) {
            if let Ok(n) = rest.split_whitespace().next().unwrap_or("0").parse() {
                return n;
            }
        }
    }

    // Try to find "=N" pattern anywhere
    for part in id.split('=') {
        if let Ok(n) = part.trim().parse() {
            return n;
        }
    }

    // Last resort: extract any digits
    let digits: String = id.chars().filter(|c| c.is_ascii_digit()).collect();
    digits.parse().unwrap_or(0)
}

/// Error type for streaming mzML parsing
#[derive(thiserror::Error, Debug)]
pub enum StreamingMzmlError {
    #[error("malformed mzML XML")]
    MalformedXml,

    #[error("XML parsing error: {0}")]
    XmlError(#[from] quick_xml::Error),

    #[error("XML attribute error: {0}")]
    AttrError(#[from] quick_xml::events::attributes::AttrError),

    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("UTF-8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),

    #[error("float parsing error: {0}")]
    FloatError(#[from] std::num::ParseFloatError),

    #[error("integer parsing error: {0}")]
    IntError(#[from] std::num::ParseIntError),

    #[error("base64 decoding error: {0}")]
    Base64Error(#[from] base64::DecodeError),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_scan_number() {
        assert_eq!(parse_scan_number("scan=123"), 123);
        assert_eq!(parse_scan_number("spectrum=456"), 456);
        assert_eq!(parse_scan_number("index=789"), 789);
        assert_eq!(parse_scan_number("scan=123 foo=bar"), 123);
        assert_eq!(parse_scan_number("some_prefix_scan=42"), 42);
        assert_eq!(parse_scan_number("no_number"), 0);
    }

    #[tokio::test]
    async fn test_parse_minimal_spectrum() {
        // Minimal mzML spectrum for testing
        let mzml = r#"<?xml version="1.0" encoding="UTF-8"?>
<mzML>
  <run>
    <spectrumList count="1">
      <spectrum id="scan=1" index="0" defaultArrayLength="3">
        <cvParam accession="MS:1000511" name="ms level" value="2"/>
        <cvParam accession="MS:1000285" name="total ion current" value="1000"/>
        <scanList count="1">
          <scan>
            <cvParam accession="MS:1000016" name="scan start time" value="1.5" unitAccession="UO:0000031"/>
          </scan>
        </scanList>
        <precursorList count="1">
          <precursor>
            <isolationWindow>
              <cvParam accession="MS:1000827" name="isolation window target m/z" value="500.0"/>
              <cvParam accession="MS:1000828" name="isolation window lower offset" value="12.5"/>
              <cvParam accession="MS:1000829" name="isolation window upper offset" value="12.5"/>
            </isolationWindow>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam accession="MS:1000744" name="selected ion m/z" value="500.0"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
        </precursorList>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="32">
            <cvParam accession="MS:1000514" name="m/z array"/>
            <cvParam accession="MS:1000523" name="64-bit float"/>
            <cvParam accession="MS:1000576" name="no compression"/>
            <binary>AAAAAAAAcEAAAAAAAAB5QAAAAAAAAIFAAAAAAACAgkA=</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="16">
            <cvParam accession="MS:1000515" name="intensity array"/>
            <cvParam accession="MS:1000521" name="32-bit float"/>
            <cvParam accession="MS:1000576" name="no compression"/>
            <binary>AADIQgAAyEIAAMhC</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"#;

        let reader = StreamingMzmlReader::ms2_only();
        let (tx, rx) = crossbeam::channel::unbounded();

        let stats = reader
            .parse_to_channel(mzml.as_bytes(), tx)
            .await
            .expect("parsing failed");

        assert_eq!(stats.total_spectra, 1);
        assert_eq!(stats.ms2_spectra, 1);

        let spectra: Vec<Spectrum> = rx.into_iter().filter_map(|r| r.ok()).collect();
        assert_eq!(spectra.len(), 1);

        let spec = &spectra[0];
        assert_eq!(spec.scan_number, 1);
        assert!((spec.retention_time - 1.5).abs() < 0.001);
        assert!((spec.precursor_mz - 500.0).abs() < 0.001);
        assert_eq!(spec.mzs.len(), 4);
        assert_eq!(spec.intensities.len(), 3);
    }
}
