# BLAS Vectorization in Osprey

## What is BLAS?

BLAS (Basic Linear Algebra Subprograms) is a specification for a set of low-level linear algebra routines -- dot products, matrix-vector multiplication, and matrix-matrix multiplication. The specification dates to the 1970s, but modern implementations (OpenBLAS, Intel MKL, Apple Accelerate) are among the most highly optimized numerical code in existence.

What makes BLAS fast is not just vectorization (SIMD instructions like AVX2/AVX-512 that process 8-16 floats per clock cycle) -- it is the combination of:

- **SIMD vectorization**: Processes multiple floats per instruction (e.g., 8 floats at once with AVX2, 16 with AVX-512)
- **Cache blocking**: Decomposes large matrix operations into sub-blocks that fit in L1/L2 cache, minimizing memory stalls
- **Register tiling**: Keeps partial results in CPU registers across inner loop iterations, avoiding store/reload overhead
- **Loop unrolling and software pipelining**: Overlaps memory loads with computation so the CPU is never waiting
- **Platform-specific tuning**: OpenBLAS includes hand-written assembly kernels for specific CPU microarchitectures (Haswell, Zen, Skylake, etc.)

The result is that a BLAS matrix multiplication typically achieves 80-95% of a CPU's theoretical peak FLOPS -- far beyond what a naive loop implementation can reach. A naive dot product in Rust might achieve 2-4 GFLOPS; OpenBLAS on the same hardware can sustain 50-100+ GFLOPS for large matrices.

## How Osprey Uses BLAS

Osprey uses BLAS for **batch spectral scoring** during the calibration discovery phase (Phase 2 of the pipeline). The key insight is that scoring thousands of library entries against hundreds of spectra can be reformulated as a single matrix multiplication.

### The Calibration Scoring Problem

During calibration, Osprey needs to find the best-matching spectrum for each library peptide across a retention time window. This means computing spectral similarity (LibCosine) between every library entry and every spectrum in the window:

- ~5,000 library entries per isolation window
- ~500 spectra per RT window
- = 2.5 million pairwise similarity scores

### Matrix Formulation

Each library entry and each spectrum is preprocessed into a fixed-length vector of binned intensities (2,001 bins spanning 200-2000 Da at 1.0005 Da width):

1. **Bin** fragment m/z values into the bin array
2. **sqrt(intensity)** transformation to down-weight dominant peaks
3. **L2 normalize** so the dot product equals cosine similarity

This converts spectral matching into a matrix multiplication:

```
scores = library_matrix @ spectra_matrix.T

  (n_library × n_bins) · (n_bins × n_spectra) → (n_library × n_spectra)
```

A single BLAS call (`sgemm` under the hood) computes all 2.5 million cosine similarities at once. This is approximately 20x faster than computing each pair individually, because BLAS can exploit cache locality across the shared bin dimension and amortize memory access costs.

### Implementation

The BLAS-accelerated scoring lives in `crates/osprey-scoring/src/batch.rs`:

- **`PreprocessedLibrary`**: Stores library fragments as an f32 matrix (n_entries x n_bins), built once per isolation window
- **`PreprocessedSpectra`**: Stores preprocessed spectra as an f32 matrix (n_spectra x n_bins), built once per RT window
- **`BatchScorer::score_all()`**: The single BLAS call -- `library.matrix.dot(&spectra.matrix.t())` -- via ndarray's BLAS-backed `dot()` method
- **`BatchScorer::score_entry_vs_all()`**: Scores a single library entry against all spectra (matrix-vector multiply, `sgemv`)

The f32 data type is deliberate: it halves memory relative to f64, doubles SIMD throughput (16 floats/AVX-512 vs 8 doubles), and calibration scoring does not need f64 precision. Results are converted to f64 only on output.

### Dependency Chain

```
ndarray (with "blas" feature)
  └── blas-src (with "openblas" feature)
       └── openblas-src (with "system" and "cblas" features)
            └── system-installed OpenBLAS (libopenblas-dev on Linux)
```

The `system` feature means Osprey links against the platform's installed OpenBLAS rather than compiling it from source, which avoids a lengthy build step.

## Where BLAS is Not Used

### Main Coelution Search (Phase 3)

The main search does not use batch matrix multiplication. Instead, each precursor's XCorr is computed individually using a precomputed bin vector (`preprocess_spectrum_for_xcorr()`) and an O(n_fragments) bin lookup. The per-window preprocessing is done once and reused across all candidates in that window, but the scoring itself is a simple dot product over the matched bins -- not a full matrix multiply.

The reason is architectural: the main search processes precursors in parallel across isolation windows via Rayon, and each precursor only needs to be scored against a handful of spectra within its detected peak boundaries. The all-pairs matrix formulation would waste computation on library-spectrum pairs outside the RT window.

### ML Crate (osprey-ml)

The `osprey-ml` crate implements its own matrix operations (`Matrix::dot()`, `Matrix::dotv()`) using manual loops with Rayon parallelization, without BLAS. These are used for LDA scatter matrices and SVM kernel computations where the matrices are small (e.g., 4x4 for LDA, n_features x n_features for SVM) and the overhead of BLAS dispatch would negate any benefit.

## Why OpenBLAS?

Osprey uses OpenBLAS specifically because:

- **Open source**: BSD-licensed, no commercial restrictions (unlike Intel MKL)
- **Cross-platform**: Works on Linux, macOS, and Windows
- **Good performance**: Within 5-15% of Intel MKL on Intel CPUs, and often faster on AMD
- **System packages available**: `apt install libopenblas-dev` on Ubuntu, `brew install openblas` on macOS
- **Pre-built binaries**: Available for Windows from the OpenBLAS GitHub releases

Alternative BLAS implementations could be swapped in by changing the `blas-src` feature flags (e.g., `features = ["intel-mkl"]` for MKL), but OpenBLAS provides the best balance of performance, portability, and ease of installation.
