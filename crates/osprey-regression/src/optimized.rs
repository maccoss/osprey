//! Optimized regression with pre-binned library and f32 operations
//!
//! Performance optimizations:
//! 1. Pre-bin library entries once at startup (avoids redundant binning)
//! 2. Use f32 for matrix operations (faster BLAS, less memory)
//! 3. Cache binned observed spectra (reuse between calibration and main search)
//! 4. Thread-local buffer reuse (avoid allocation churn)

use crate::cd_nnls::{solve_cd_nnls_f32, CdNnlsParams};
use crate::Binner;
use ndarray::{Array1, Array2};
use osprey_core::{LibraryEntry, Result, Spectrum};
use std::collections::HashMap;

/// Pre-binned library for efficient regression
///
/// Bins all library entries once at startup, storing them in a contiguous
/// matrix for fast submatrix extraction during regression.
#[derive(Debug, Clone)]
pub struct BinnedLibrary {
    /// Pre-computed normalized dense vectors for each library entry
    /// Shape: (n_entries, n_bins) - row-major for efficient row access
    dense_matrix: Array2<f32>,
    /// Mapping from library entry ID to row index in dense_matrix
    id_to_row: HashMap<u32, usize>,
    /// Mapping from row index to library entry ID
    row_to_id: Vec<u32>,
    /// Number of bins
    n_bins: usize,
}

impl BinnedLibrary {
    /// Create a new pre-binned library
    ///
    /// Bins all library entries once and stores them in a contiguous matrix.
    /// This is done once at startup and amortizes the binning cost.
    pub fn new(library: &[LibraryEntry], binner: &Binner) -> Self {
        let n_bins = binner.n_bins();
        let n_entries = library.len();

        log::info!(
            "Pre-binning {} library entries into {} bins...",
            n_entries,
            n_bins
        );

        let mut dense_matrix = Array2::<f32>::zeros((n_entries, n_bins));
        let mut id_to_row = HashMap::with_capacity(n_entries);
        let mut row_to_id = Vec::with_capacity(n_entries);

        for (row, entry) in library.iter().enumerate() {
            // Bin the library entry
            let binned = binner.bin_library_entry(entry);

            // Fill the row
            for (bin, intensity) in binned.bin_indices.iter().zip(binned.intensities.iter()) {
                dense_matrix[[row, *bin as usize]] = *intensity;
            }

            // Normalize the row to unit sum
            let sum: f32 = dense_matrix.row(row).sum();
            if sum > 0.0 {
                dense_matrix.row_mut(row).mapv_inplace(|x| x / sum);
            }

            id_to_row.insert(entry.id, row);
            row_to_id.push(entry.id);
        }

        log::info!(
            "Pre-binned library: {} entries × {} bins = {:.1} MB",
            n_entries,
            n_bins,
            (n_entries * n_bins * 4) as f64 / 1_000_000.0
        );

        Self {
            dense_matrix,
            id_to_row,
            row_to_id,
            n_bins,
        }
    }

    /// Get the number of library entries
    pub fn len(&self) -> usize {
        self.dense_matrix.nrows()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.dense_matrix.nrows() == 0
    }

    /// Get the number of bins
    pub fn n_bins(&self) -> usize {
        self.n_bins
    }

    /// Get the row index for a library ID
    pub fn get_row(&self, id: u32) -> Option<usize> {
        self.id_to_row.get(&id).copied()
    }

    /// Get the library ID for a row index
    pub fn get_id(&self, row: usize) -> Option<u32> {
        self.row_to_id.get(row).copied()
    }

    /// Extract a design matrix for the given candidate indices
    ///
    /// Returns a matrix of shape (n_bins, n_candidates) where each column
    /// is the pre-binned spectrum for that candidate.
    ///
    /// This is the key optimization: no re-binning, just copying pre-computed rows.
    pub fn extract_design_matrix(&self, candidate_indices: &[usize]) -> (Array2<f32>, Vec<u32>) {
        let n_candidates = candidate_indices.len();
        let mut matrix = Array2::<f32>::zeros((self.n_bins, n_candidates));
        let mut ids = Vec::with_capacity(n_candidates);

        for (col, &row_idx) in candidate_indices.iter().enumerate() {
            // Copy pre-binned row into column (transpose)
            let row = self.dense_matrix.row(row_idx);
            for (bin_idx, &val) in row.iter().enumerate() {
                matrix[[bin_idx, col]] = val;
            }
            ids.push(self.row_to_id[row_idx]);
        }

        (matrix, ids)
    }

    /// Extract design matrix for candidates specified by library ID
    pub fn extract_design_matrix_by_id(&self, candidate_ids: &[u32]) -> (Array2<f32>, Vec<u32>) {
        let indices: Vec<usize> = candidate_ids
            .iter()
            .filter_map(|id| self.id_to_row.get(id).copied())
            .collect();
        self.extract_design_matrix(&indices)
    }
}

/// Cache for binned observed spectra
///
/// Stores binned spectra to avoid re-binning the same spectrum multiple times
/// (e.g., once during calibration discovery, again during main search).
#[derive(Debug, Clone)]
pub struct BinnedSpectraCache {
    /// scan_number → binned dense vector (f32)
    cache: HashMap<u32, Array1<f32>>,
    /// Number of bins expected
    n_bins: usize,
}

impl BinnedSpectraCache {
    /// Create a new cache with expected bin count
    pub fn new(n_bins: usize) -> Self {
        Self {
            cache: HashMap::new(),
            n_bins,
        }
    }

    /// Create a cache with pre-allocated capacity
    pub fn with_capacity(n_bins: usize, capacity: usize) -> Self {
        Self {
            cache: HashMap::with_capacity(capacity),
            n_bins,
        }
    }

    /// Get a cached binned spectrum, or bin it and cache
    pub fn get_or_bin(&mut self, spectrum: &Spectrum, binner: &Binner) -> Array1<f32> {
        if let Some(cached) = self.cache.get(&spectrum.scan_number) {
            return cached.clone();
        }

        // Bin the spectrum
        let binned = binner.bin_spectrum(spectrum);
        let mut dense = Array1::<f32>::zeros(self.n_bins);
        for (bin, intensity) in binned.bin_indices.iter().zip(binned.intensities.iter()) {
            dense[*bin as usize] = *intensity;
        }

        // Cache and return
        self.cache.insert(spectrum.scan_number, dense.clone());
        dense
    }

    /// Pre-populate cache with all spectra (for batch processing)
    pub fn populate(&mut self, spectra: &[Spectrum], binner: &Binner) {
        log::info!("Pre-binning {} observed spectra...", spectra.len());

        for spectrum in spectra {
            if !self.cache.contains_key(&spectrum.scan_number) {
                let binned = binner.bin_spectrum(spectrum);
                let mut dense = Array1::<f32>::zeros(self.n_bins);
                for (bin, intensity) in binned.bin_indices.iter().zip(binned.intensities.iter()) {
                    dense[*bin as usize] = *intensity;
                }
                self.cache.insert(spectrum.scan_number, dense);
            }
        }

        log::info!(
            "Cached {} binned spectra ({:.1} MB)",
            self.cache.len(),
            (self.cache.len() * self.n_bins * 4) as f64 / 1_000_000.0
        );
    }

    /// Get a cached spectrum by scan number
    pub fn get(&self, scan_number: u32) -> Option<&Array1<f32>> {
        self.cache.get(&scan_number)
    }

    /// Check if a spectrum is cached
    pub fn contains(&self, scan_number: u32) -> bool {
        self.cache.contains_key(&scan_number)
    }

    /// Get the number of cached spectra
    pub fn len(&self) -> usize {
        self.cache.len()
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.cache.is_empty()
    }

    /// Clear the cache
    pub fn clear(&mut self) {
        self.cache.clear();
    }
}

// Note: Thread-local buffer reuse can be added here for further optimization
// Currently not implemented as the main gains come from pre-binning and f32.

/// Optimized ridge regression solver using f32
///
/// Key optimizations:
/// - Uses f32 instead of f64 (faster BLAS, less memory)
/// - Coordinate descent NNLS solver (O(m*k) per sweep vs O(k³) for Cholesky)
/// - Active set acceleration (skip variables at zero)
/// - Works with pre-binned library (BinnedLibrary)
#[derive(Debug, Clone)]
pub struct OptimizedSolver {
    /// Default regularization parameter
    default_lambda: f32,
    /// Coordinate descent parameters
    cd_params: CdNnlsParams,
}

impl OptimizedSolver {
    /// Create a new optimized solver with default CD parameters
    pub fn new(default_lambda: f32) -> Self {
        Self {
            default_lambda,
            cd_params: CdNnlsParams::default(),
        }
    }

    /// Create a new optimized solver with custom CD parameters
    pub fn with_params(default_lambda: f32, cd_params: CdNnlsParams) -> Self {
        Self {
            default_lambda,
            cd_params,
        }
    }

    /// Solve the NNLS ridge regression problem using coordinate descent
    ///
    /// Given design matrix A (m × k) and observed vector b (m × 1),
    /// returns coefficient vector x (k × 1) that minimizes:
    ///   ‖Ax - b‖² + λ‖x‖² subject to x ≥ 0
    ///
    /// Uses coordinate descent with active set acceleration for O(m*k) per sweep
    /// instead of O(k³) for Cholesky decomposition.
    pub fn solve_nonnegative(
        &self,
        a: &Array2<f32>,
        b: &Array1<f32>,
        lambda: Option<f32>,
    ) -> Result<Array1<f32>> {
        let lambda = lambda.unwrap_or(self.default_lambda);
        solve_cd_nnls_f32(a, b, lambda, &self.cd_params)
    }

    /// Compute residual norm ‖Ax - b‖²
    pub fn residual(&self, a: &Array2<f32>, x: &Array1<f32>, b: &Array1<f32>) -> f32 {
        let ax = a.dot(x);
        let diff = &ax - b;
        diff.dot(&diff)
    }
}

impl Default for OptimizedSolver {
    fn default() -> Self {
        Self {
            default_lambda: 1.0,
            cd_params: CdNnlsParams::default(),
        }
    }
}

/// Result from optimized regression
#[derive(Debug, Clone)]
pub struct OptimizedRegressionResult {
    /// Scan number
    pub scan_number: u32,
    /// Retention time
    pub retention_time: f64,
    /// Library IDs with non-zero coefficients
    pub library_ids: Vec<u32>,
    /// Corresponding coefficients (f32 for memory efficiency)
    pub coefficients: Vec<f32>,
    /// Residual norm ‖Ax - b‖²
    pub residual: f32,
    /// Number of candidates considered
    pub n_candidates: u32,
    /// Sum of all coefficients (before filtering)
    pub coefficient_sum: f32,
    /// Observed spectrum norm ‖b‖²
    pub observed_norm: f32,
}

impl OptimizedRegressionResult {
    /// Create an empty result
    pub fn empty(scan_number: u32, retention_time: f64) -> Self {
        Self {
            scan_number,
            retention_time,
            library_ids: Vec::new(),
            coefficients: Vec::new(),
            residual: 0.0,
            n_candidates: 0,
            coefficient_sum: 0.0,
            observed_norm: 0.0,
        }
    }

    /// Compute explained variance (1 - residual/||b||²)
    pub fn explained_variance(&self) -> f32 {
        if self.observed_norm > 1e-10 {
            1.0 - self.residual / self.observed_norm
        } else {
            0.0
        }
    }

    /// Get relative coefficient for a library ID
    pub fn relative_coefficient(&self, lib_id: u32) -> f32 {
        if self.coefficient_sum > 1e-10 {
            if let Some(pos) = self.library_ids.iter().position(|&id| id == lib_id) {
                return self.coefficients[pos] / self.coefficient_sum;
            }
        }
        0.0
    }

    /// Convert to f64 coefficients for compatibility
    pub fn coefficients_f64(&self) -> Vec<f64> {
        self.coefficients.iter().map(|&c| c as f64).collect()
    }
}

/// Optimized regression processor that combines all optimizations
pub struct OptimizedProcessor {
    /// Pre-binned library
    pub binned_library: BinnedLibrary,
    /// Cached binned spectra
    pub spectra_cache: BinnedSpectraCache,
    /// The binner (for new spectra)
    pub binner: Binner,
    /// The solver
    pub solver: OptimizedSolver,
}

impl OptimizedProcessor {
    /// Create a new optimized processor
    pub fn new(library: &[LibraryEntry], binner: Binner, lambda: f32) -> Self {
        let binned_library = BinnedLibrary::new(library, &binner);
        let spectra_cache = BinnedSpectraCache::new(binner.n_bins());
        let solver = OptimizedSolver::new(lambda);

        Self {
            binned_library,
            spectra_cache,
            binner,
            solver,
        }
    }

    /// Pre-cache all spectra for batch processing
    pub fn cache_spectra(&mut self, spectra: &[Spectrum]) {
        self.spectra_cache.populate(spectra, &self.binner);
    }

    /// Process a single spectrum with candidate indices
    ///
    /// candidate_indices: indices into the original library (not row indices)
    pub fn process_spectrum(
        &mut self,
        spectrum: &Spectrum,
        candidate_indices: &[usize],
    ) -> Result<OptimizedRegressionResult> {
        if candidate_indices.is_empty() {
            return Ok(OptimizedRegressionResult::empty(
                spectrum.scan_number,
                spectrum.retention_time,
            ));
        }

        // Get or compute binned observed spectrum
        let observed = self.spectra_cache.get_or_bin(spectrum, &self.binner);
        let observed_norm = observed.dot(&observed);

        // Extract design matrix from pre-binned library
        let (design_matrix, library_ids) = self.binned_library.extract_design_matrix(candidate_indices);

        // Solve regression
        let coefficients = self.solver.solve_nonnegative(&design_matrix, &observed, None)?;
        let coefficient_sum: f32 = coefficients.iter().sum();

        // Compute residual
        let predicted = design_matrix.dot(&coefficients);
        let residual = (&observed - &predicted).mapv(|x| x * x).sum();

        // Filter to non-zero coefficients
        let mut result = OptimizedRegressionResult::empty(spectrum.scan_number, spectrum.retention_time);
        result.n_candidates = candidate_indices.len() as u32;
        result.coefficient_sum = coefficient_sum;
        result.observed_norm = observed_norm;
        result.residual = residual;

        for (id, coef) in library_ids.iter().zip(coefficients.iter()) {
            if *coef > 1e-6 {
                result.library_ids.push(*id);
                result.coefficients.push(*coef);
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;
    use approx::assert_abs_diff_eq;
    use osprey_core::{FragmentAnnotation, LibraryFragment, IsolationWindow};

    fn make_test_entry(id: u32, fragments: Vec<(f64, f32)>) -> LibraryEntry {
        let mut entry = LibraryEntry::new(
            id,
            "PEPTIDE".to_string(),
            "PEPTIDE".to_string(),
            2,
            500.0,
            10.0,
        );

        entry.fragments = fragments
            .into_iter()
            .map(|(mz, intensity)| LibraryFragment {
                mz,
                relative_intensity: intensity,
                annotation: FragmentAnnotation::default(),
            })
            .collect();

        entry
    }

    #[test]
    fn test_binned_library_creation() {
        let binner = Binner::unit_resolution();
        let entries = vec![
            make_test_entry(1, vec![(300.0, 100.0), (400.0, 200.0)]),
            make_test_entry(2, vec![(350.0, 150.0), (450.0, 250.0)]),
        ];

        let lib = BinnedLibrary::new(&entries, &binner);

        assert_eq!(lib.len(), 2);
        assert_eq!(lib.get_row(1), Some(0));
        assert_eq!(lib.get_row(2), Some(1));
        assert_eq!(lib.get_id(0), Some(1));
    }

    #[test]
    fn test_extract_design_matrix() {
        let binner = Binner::unit_resolution();
        let entries = vec![
            make_test_entry(1, vec![(300.0, 100.0), (400.0, 200.0)]),
            make_test_entry(2, vec![(350.0, 150.0), (450.0, 250.0)]),
            make_test_entry(3, vec![(300.0, 50.0), (500.0, 100.0)]),
        ];

        let lib = BinnedLibrary::new(&entries, &binner);

        // Extract subset
        let (matrix, ids) = lib.extract_design_matrix(&[0, 2]);

        assert_eq!(matrix.ncols(), 2);
        assert_eq!(ids, vec![1, 3]);

        // Columns should be normalized
        let col0_sum: f32 = matrix.column(0).sum();
        let col1_sum: f32 = matrix.column(1).sum();
        assert!((col0_sum - 1.0).abs() < 1e-5);
        assert!((col1_sum - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_optimized_solver_nonnegative() {
        let solver = OptimizedSolver::new(0.1);

        let a = array![
            [1.0_f32, 0.5, 0.2],
            [0.5, 1.0, 0.3],
            [0.2, 0.3, 1.0],
        ];
        let b = array![1.0_f32, 0.5, 0.8];

        let x = solver.solve_nonnegative(&a, &b, None).unwrap();

        // All non-negative
        assert!(x.iter().all(|&v| v >= 0.0));

        // Reasonable residual
        let residual = solver.residual(&a, &x, &b);
        assert!(residual < 1.0);
    }

    #[test]
    fn test_spectra_cache() {
        let binner = Binner::unit_resolution();
        let mut cache = BinnedSpectraCache::new(binner.n_bins());

        let spectrum = Spectrum {
            scan_number: 42,
            retention_time: 10.0,
            precursor_mz: 500.0,
            isolation_window: IsolationWindow::symmetric(500.0, 12.5),
            mzs: vec![300.0, 400.0, 500.0],
            intensities: vec![100.0, 200.0, 300.0],
        };

        // First call should bin
        let binned1 = cache.get_or_bin(&spectrum, &binner);
        assert!(cache.contains(42));

        // Second call should return cached
        let binned2 = cache.get_or_bin(&spectrum, &binner);
        assert_eq!(binned1, binned2);
    }

    #[test]
    fn test_optimized_result_explained_variance() {
        let mut result = OptimizedRegressionResult::empty(1, 10.0);
        result.observed_norm = 100.0;
        result.residual = 20.0;

        assert_abs_diff_eq!(result.explained_variance(), 0.8, epsilon = 1e-6);
    }
}
