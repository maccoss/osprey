//! Buffer pool for XCorr preprocessing scratch space.
//!
//! XCorr preprocessing allocates three `NBins`-sized `f32` buffers per call
//! (binned accumulator, windowed result, sliding-window prefix sum). On
//! HRAM binning (NBins ~= 100K) each buffer is ~400 KB, so per-window
//! scoring across ~1000 spectra churns through ~1.2 GB of transient
//! `Vec<f32>` allocations plus the ~400 MB per-window preprocessed cache
//! itself. All of this memory is dropped at the end of each window and
//! re-allocated by the next window -- pure churn.
//!
//! This pool hands out reusable scratch bundles and single output buffers
//! so that after the first few windows reach the high-water mark the hot
//! path never allocates. Mirrors the C# OspreySharp `XcorrScratchPool` in
//! `pwiz_tools/OspreySharp/OspreySharp.Scoring/XcorrScratchPool.cs`.
//!
//! Thread-safety: the pool uses `Mutex<Vec<_>>` bags. Contention is low
//! because each Rayon worker rents at most a few items per window, and the
//! critical section is a single push/pop.

use std::sync::Mutex;

/// A bundle of scratch buffers reused across XCorr preprocess calls.
///
/// Rented by value, returned by value via [`XcorrScratchPool::recycle`].
pub struct XcorrScratch {
    /// `NBins` accumulator for binned-and-sqrt-transformed spectrum
    /// intensities. Written via `+=` so must start zeroed; zeroed on
    /// recycle.
    pub binned: Vec<f32>,
    /// `NBins` windowing-normalization output. Filled by
    /// `apply_windowing_normalization_into`; zeroed on recycle so
    /// below-threshold positions retain 0.
    pub windowed: Vec<f32>,
    /// `NBins + 1` sliding-window prefix sum. `prefix[0]` is always 0;
    /// positions `1..=NBins` are fully overwritten each call.
    pub prefix: Vec<f32>,
}

impl XcorrScratch {
    /// Allocate a fresh scratch bundle not backed by a pool. Useful for
    /// one-off preprocessing calls where setting up a pool is overkill
    /// (e.g. unit tests, `SpectralScorer::preprocess_spectrum_for_xcorr`).
    pub fn new_unpooled(n_bins: usize) -> Self {
        Self {
            binned: vec![0.0; n_bins],
            windowed: vec![0.0; n_bins],
            prefix: vec![0.0; n_bins + 1],
        }
    }
}

/// Thread-safe pool of XCorr scratch bundles and per-window output buffers.
///
/// Two independent bags:
///
/// * `scratch` -- [`XcorrScratch`] bundles, one per in-flight preprocess
///   call. High-water mark is roughly the Rayon worker count.
/// * `bins` -- single `Vec<f32>` of length `NBins`, used for the per-window
///   `preprocessed_xcorr` cache. High-water mark is roughly
///   `n_workers * spectra_per_window`.
pub struct XcorrScratchPool {
    n_bins: usize,
    scratch: Mutex<Vec<XcorrScratch>>,
    bins: Mutex<Vec<Vec<f32>>>,
}

impl XcorrScratchPool {
    pub fn new(n_bins: usize) -> Self {
        Self {
            n_bins,
            scratch: Mutex::new(Vec::new()),
            bins: Mutex::new(Vec::new()),
        }
    }

    pub fn n_bins(&self) -> usize {
        self.n_bins
    }

    /// Number of scratch bundles currently idle in the pool. Diagnostic
    /// only; the value may change before the caller sees it.
    pub fn scratch_idle(&self) -> usize {
        self.scratch.lock().unwrap().len()
    }

    /// Number of bins buffers currently idle in the pool. Diagnostic only.
    pub fn bins_idle(&self) -> usize {
        self.bins.lock().unwrap().len()
    }

    /// Rent a scratch bundle. Allocates a fresh one if the bag is empty.
    pub fn rent(&self) -> XcorrScratch {
        if let Some(s) = self.scratch.lock().unwrap().pop() {
            return s;
        }
        XcorrScratch::new_unpooled(self.n_bins)
    }

    /// Return a scratch bundle to the pool. Does not zero the buffers --
    /// the preprocess functions zero the accumulator fields on each call,
    /// which is both the correct behaviour when callers reuse a single
    /// rented scratch across many preprocess calls and the cheapest
    /// policy overall (the memset happens once per use regardless).
    ///
    /// Scratches with mismatched buffer sizes are dropped defensively.
    pub fn recycle(&self, s: XcorrScratch) {
        if s.binned.len() != self.n_bins
            || s.windowed.len() != self.n_bins
            || s.prefix.len() != self.n_bins + 1
        {
            return;
        }
        self.scratch.lock().unwrap().push(s);
    }

    /// Rent a single `Vec<f32>` of length `NBins`. The buffer is always
    /// fully overwritten by the preprocessing pipeline (sliding-window
    /// output), so no zero-on-rent is required.
    pub fn rent_bins(&self) -> Vec<f32> {
        if let Some(b) = self.bins.lock().unwrap().pop() {
            return b;
        }
        vec![0.0; self.n_bins]
    }

    /// Return a single `Vec<f32>` to the pool. Buffers of the wrong length
    /// are dropped (defensive -- should not happen in practice).
    pub fn recycle_bins(&self, b: Vec<f32>) {
        if b.len() != self.n_bins {
            return;
        }
        self.bins.lock().unwrap().push(b);
    }

    /// Return every buffer in a per-window cache.
    pub fn recycle_bins_all(&self, cache: Vec<Vec<f32>>) {
        if cache.is_empty() {
            return;
        }
        let mut guard = self.bins.lock().unwrap();
        for b in cache {
            if b.len() == self.n_bins {
                guard.push(b);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rent_and_recycle_reuses_allocation() {
        let pool = XcorrScratchPool::new(100);
        let s = pool.rent();
        assert_eq!(s.binned.len(), 100);
        assert_eq!(s.windowed.len(), 100);
        assert_eq!(s.prefix.len(), 101);
        let binned_ptr = s.binned.as_ptr();
        pool.recycle(s);
        let s2 = pool.rent();
        assert_eq!(
            s2.binned.as_ptr(),
            binned_ptr,
            "allocation should be reused"
        );
    }

    #[test]
    fn bins_rent_and_recycle() {
        let pool = XcorrScratchPool::new(50);
        let b = pool.rent_bins();
        assert_eq!(b.len(), 50);
        let ptr = b.as_ptr();
        pool.recycle_bins(b);
        let b2 = pool.rent_bins();
        assert_eq!(b2.as_ptr(), ptr, "bins allocation should be reused");
    }

    #[test]
    fn wrong_length_bins_dropped() {
        let pool = XcorrScratchPool::new(50);
        pool.recycle_bins(vec![0.0; 49]);
        assert_eq!(pool.bins_idle(), 0, "wrong-length buffer should be dropped");
    }

    #[test]
    fn recycle_bins_all_returns_every_buffer() {
        let pool = XcorrScratchPool::new(10);
        let cache: Vec<Vec<f32>> = (0..5).map(|_| pool.rent_bins()).collect();
        pool.recycle_bins_all(cache);
        assert_eq!(pool.bins_idle(), 5);
    }

    #[test]
    fn pool_is_send_and_sync() {
        fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<XcorrScratchPool>();
    }
}
