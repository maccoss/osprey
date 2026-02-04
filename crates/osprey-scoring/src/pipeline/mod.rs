//! Pipelined processing for streaming mzML analysis
//!
//! This module provides a streaming architecture for processing large mzML files
//! with overlapped I/O, preprocessing, and scoring.
//!
//! ## Architecture
//!
//! ```text
//! mzML File
//!     │
//!     ▼ (async streaming)
//! StreamingMzmlReader
//!     │
//!     ▼ (crossbeam channel)
//! PreprocessingWorker(s)
//!     │
//!     ▼
//! WindowAccumulator
//!     │
//!     ▼ (when window complete)
//! BatchScorer
//!     │
//!     ▼
//! CalibrationMatch Results
//! ```
//!
//! ## Usage
//!
//! ```ignore
//! use osprey_scoring::pipeline::{PipelineCoordinator, PipelineConfig};
//!
//! let coordinator = PipelineCoordinator::new();
//! let matches = coordinator.run_on_spectra(&spectra, &library);
//! ```
//!
//! For streaming from files (requires `streaming` feature):
//!
//! ```ignore
//! use osprey_scoring::pipeline::run_streaming_pipeline;
//!
//! let matches = run_streaming_pipeline("file.mzML", &library, config).await?;
//! ```

mod accumulator;
mod preprocessor;

#[cfg(feature = "streaming")]
mod coordinator;

pub use accumulator::{AccumulatorStats, WindowAccumulator, WindowKey};
pub use preprocessor::{PreprocessedSpectrum, PreprocessingConfig, PreprocessingWorker};

#[cfg(feature = "streaming")]
pub use coordinator::{PipelineConfig, PipelineCoordinator, WindowResult, run_streaming_pipeline, StreamingPipelineError};
