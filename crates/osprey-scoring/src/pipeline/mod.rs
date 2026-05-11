//! Pipeline components for XCorr preprocessing
//!
//! This module provides preprocessing workers for XCorr scoring.
//! XCorr uses binning (Comet-style) and can be preprocessed ahead of time.
//!
//! **Note:** LibCosine scoring uses PPM matching (not binning) and is computed
//! on-demand using `LibCosineScorer`. It does not use preprocessed vectors.

mod accumulator;
mod preprocessor;

pub use accumulator::{AccumulatorStats, WindowAccumulator, WindowKey};
pub use preprocessor::{PreprocessedSpectrum, PreprocessingConfig, PreprocessingWorker};
