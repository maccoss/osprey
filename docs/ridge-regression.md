# Ridge Regression for Peptide Detection

Osprey uses ridge regression (L2-regularized least squares) to deconvolute DIA spectra and detect peptides. This is the core algorithm that distinguishes peptide-centric from spectrum-centric approaches.

## The DIA Problem

In Data-Independent Acquisition (DIA), each MS2 spectrum contains fragments from **multiple precursors**:

```
Isolation window: 500-525 m/z

Precursors in window:
  - Peptide A: 502.3 m/z
  - Peptide B: 510.7 m/z
  - Peptide C: 518.2 m/z

Observed spectrum = fragments from A + B + C (mixed)
```

## The Regression Model

We model the observed spectrum as a linear combination of library spectra:

```
observed = c₁ × library₁ + c₂ × library₂ + ... + cₙ × libraryₙ + noise
```

Or in matrix form:

```
b = Ax + ε

Where:
  b = observed spectrum (m bins)
  A = design matrix (m × n), each column is a library spectrum
  x = coefficients (n × 1), one per candidate peptide
  ε = noise/unexplained signal
```

## Why Ridge Regression?

### The Problem with OLS

Ordinary Least Squares (OLS) minimizes `‖Ax - b‖²` but has problems:

1. **Ill-conditioned**: Similar library spectra cause numerical instability
2. **Overfitting**: May assign signal to wrong peptides
3. **Negative coefficients**: Physically meaningless

### Ridge Regression Solution

Ridge regression adds L2 regularization:

```
minimize: ‖Ax - b‖² + λ‖x‖²
```

This:
- **Stabilizes** the solution when columns are correlated
- **Shrinks** coefficients toward zero (reduces overfitting)
- **Always has a unique solution** (even when A is singular)

### Closed-Form Solution

```
x = (AᵀA + λI)⁻¹ Aᵀb
```

## Design Matrix Construction

### Unit Resolution Mode

For unit-resolution instruments (~1 Th bins):

```
Bin width: 1.0005 Th
Range: 200-2000 m/z
Bins: ~1800

For each library spectrum:
  1. Take fragment (m/z, intensity) pairs
  2. Assign each to bin: bin_idx = floor((mz - 200) / 1.0005)
  3. Sum intensities in each bin
  4. Normalize to sum = 1
```

### HRAM Mode

For high-resolution accurate mass (HRAM) instruments:

```
Tolerance: 20 ppm

Instead of fixed bins:
  1. Build sparse matrix (CSC format)
  2. Match observed peaks to library peaks within ppm tolerance
  3. Use conjugate gradient solver for large problems
```

## Candidate Selection

Before building the design matrix, we filter candidates:

```python
def select_candidates(spectrum, library, rt_tolerance, calibration):
    candidates = []

    for entry in library:
        # 1. Check isolation window (from mzML)
        if not spectrum.isolation_window.contains(entry.precursor_mz):
            continue

        # 2. Check RT tolerance (with calibration if available)
        if calibration:
            expected_rt = calibration.predict(entry.retention_time)
        else:
            expected_rt = entry.retention_time

        if abs(expected_rt - spectrum.retention_time) > rt_tolerance:
            continue

        candidates.append(entry)

    # 3. Limit to max candidates (keep closest in RT)
    if len(candidates) > max_candidates:
        candidates.sort(key=lambda e: abs(e.rt - spectrum.rt))
        candidates = candidates[:max_candidates]

    return candidates
```

## Regularization Parameter (λ)

### Fixed λ

```yaml
regularization_lambda:
  Fixed: 0.5
```

### Cross-Validated λ

Automatically select λ using leave-one-out cross-validation:

```python
def cross_validate_lambda(A, b, lambdas=[0.01, 0.1, 1.0, 10.0]):
    best_lambda = None
    best_error = float('inf')

    for lam in lambdas:
        # Leave-one-out CV using the "hat matrix" trick
        H = A @ inv(A.T @ A + lam * I) @ A.T
        residuals = b - H @ b
        cv_error = sum((r / (1 - h_ii))**2 for r, h_ii in zip(residuals, diag(H)))

        if cv_error < best_error:
            best_error = cv_error
            best_lambda = lam

    return best_lambda
```

### Adaptive λ

Scale λ based on matrix condition number or spectrum complexity.

## Non-Negative Constraints

Peptide abundances can't be negative. Osprey enforces this with:

### Active Set Method

```python
def solve_nonnegative(A, b, lambda_):
    x = ridge_solve(A, b, lambda_)  # Initial solution

    while any(x < 0):
        # Remove candidates with negative coefficients
        positive_mask = x >= 0
        A_reduced = A[:, positive_mask]
        x_reduced = ridge_solve(A_reduced, b, lambda_)

        x = zeros(n)
        x[positive_mask] = x_reduced

    return x
```

## Solving the System

### Small Problems (< 200 candidates)

Use Cholesky decomposition:

```
(AᵀA + λI) x = Aᵀb

1. Compute C = AᵀA + λI
2. Cholesky: C = LLᵀ
3. Forward solve: Ly = Aᵀb
4. Backward solve: Lᵀx = y
```

### Large Problems (≥ 200 candidates)

Use Conjugate Gradient (CG) solver:

```python
def cg_solve(A, b, lambda_, tol=1e-6, max_iter=1000):
    # Solve (AᵀA + λI)x = Aᵀb iteratively
    x = zeros(n)
    r = A.T @ b - (A.T @ A + lambda_ * I) @ x
    p = r.copy()

    for i in range(max_iter):
        Ap = A.T @ (A @ p) + lambda_ * p
        alpha = (r @ r) / (p @ Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap

        if norm(r_new) < tol:
            break

        beta = (r_new @ r_new) / (r @ r)
        p = r_new + beta * p
        r = r_new

    return x
```

## Interpreting Coefficients

### Coefficient Meaning

- **Coefficient > 0**: Peptide contributes to the spectrum
- **Coefficient ≈ 0**: Peptide not detected in this spectrum
- **Larger coefficient**: Higher relative abundance (in this spectrum)

### Important Notes

1. **Coefficients are NOT abundances** - They're relative contributions to one spectrum
2. **Aggregate across RT** - Sum coefficients across the peak for quantification
3. **Compare within-spectrum** - Coefficients are comparable within a spectrum, not across spectra

## Residual Analysis

The residual quantifies unexplained signal:

```
residual = ‖b - Ax‖² = sum of squared unexplained intensities
```

High residual may indicate:
- Missing peptides in the library
- Contaminants or noise
- Poor spectral prediction quality

## Example Walkthrough

```
Spectrum at RT=25.3 min, isolation window 500-525 m/z

Candidates (after filtering):
  1. PEPTIDEK (m/z 501.2, library RT 25.0)
  2. ANOTHERK (m/z 515.3, library RT 25.5)
  3. THIRDPEP (m/z 520.1, library RT 24.8)

Design matrix A (1800 bins × 3 candidates):
  Column 1: binned PEPTIDEK spectrum
  Column 2: binned ANOTHERK spectrum
  Column 3: binned THIRDPEP spectrum

Observed spectrum b (1800 bins)

Solve: (AᵀA + λI)x = Aᵀb

Result:
  x = [0.85, 0.42, 0.00]

Interpretation:
  - PEPTIDEK: major contributor (coef = 0.85)
  - ANOTHERK: minor contributor (coef = 0.42)
  - THIRDPEP: not detected (coef = 0.00)
```

## Configuration

```yaml
regularization_lambda: CrossValidated
# Options:
#   CrossValidated - automatic selection
#   Adaptive - scale by matrix properties
#   Fixed: 0.5 - use specific value

max_candidates_per_spectrum: 200
max_iterations: 1000
convergence_threshold: 1e-6
```

## Implementation

Key files:
- `crates/osprey-regression/src/ridge.rs` - RidgeSolver
- `crates/osprey-regression/src/matrix.rs` - DesignMatrixBuilder
- `crates/osprey-regression/src/binning.rs` - Binner for unit resolution
- `crates/osprey-regression/src/sparse.rs` - Sparse solver for HRAM
- `crates/osprey/src/pipeline.rs` - Integration in analysis pipeline

## Performance Considerations

| Factor | Impact | Mitigation |
|--------|--------|------------|
| Many candidates | Slow matrix operations | Limit max_candidates |
| High resolution | Sparse matrices | Use CSC format, CG solver |
| Many spectra | Linear scaling | Parallel processing (rayon) |
| Large design matrix | Memory | Process spectra in batches |
