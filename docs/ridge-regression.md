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

### Binning: Comet BIN Macro

Both unit resolution and HRAM modes use XCorr-style binning via Comet's BIN macro:

```
BIN(mass) = (int)(mass / bin_width + (1 - offset))
```

### Unit Resolution Mode

For unit-resolution instruments (~1 Da bins):

```
Bin width: 1.0005079 Da
Offset: 0.4
Formula: bin = (int)(mz / 1.0005079 + 0.6)
Bins: 2001

For each library spectrum:
  1. Take fragment (m/z, intensity) pairs
  2. Assign each to bin using BIN macro
  3. Sum intensities in each bin
  4. Normalize to sum = 1
```

### HRAM Mode

For high-resolution accurate mass (HRAM) instruments, the same binning approach is used with narrower bins:

```
Bin width: 0.02 Da
Offset: 0.0
Formula: bin = (int)(mz / 0.02 + 1.0)
Bins: 100,001

Same binning flow as unit resolution, just with finer bins.
Library candidates are binned on-the-fly per spectrum (after filtering)
to avoid pre-binning the entire library into a 100K-bin dense matrix.
```

**On-the-fly binning for HRAM**: Pre-binning the entire library (100K bins × 500K entries = 190 GB) is infeasible. Instead, for each spectrum, candidates are filtered first (isolation window + calibrated RT + top-3 fragment match), then only the ~200 passing candidates are binned into a column-major design matrix.

## M/Z Calibration Before Regression

Before candidate selection and ridge regression, observed spectra are calibrated using MS2 calibration data from the discovery phase.

### Calibration Application

```python
def apply_spectrum_calibration(spectrum, ms2_calibration):
    """
    Apply m/z calibration to shift observed peaks.

    The calibration corrects systematic mass measurement errors by
    shifting each centroid so the mean error becomes 0.
    """
    if not ms2_calibration.calibrated:
        return spectrum  # No correction needed

    corrected_mzs = []
    for mz in spectrum.mzs:
        if ms2_calibration.unit == "ppm":
            # PPM correction: shift = mz × mean_error / 1e6
            correction = mz * ms2_calibration.mean / 1e6
        else:
            # Absolute correction (Th)
            correction = ms2_calibration.mean

        corrected_mzs.append(mz - correction)

    return Spectrum(mzs=corrected_mzs, intensities=spectrum.intensities, ...)
```

### Why Calibrate Before Regression?

1. **Top-3 pre-filter accuracy**: Matching top-3 peaks uses calibrated m/z and 3×SD tolerance
2. **Better binning**: Calibrated peaks land in correct bins for regression
3. **Consistent scoring**: Both regression and post-regression scoring use calibrated spectra

### Calibration Flow

```
1. Discovery phase:
   - Find high-confidence matches
   - Compute MS2 error distribution (mean, SD)

2. Full search phase:
   - Apply mean shift to all spectra (center error at 0)
   - Use 3×SD as matching tolerance
   - Run candidate selection + ridge regression on calibrated spectra

3. Scoring phase:
   - Spectra already calibrated from step 2
   - Spectral scores computed on calibrated data
```

## Candidate Selection

Before building the design matrix, we filter candidates using a multi-stage approach:

```python
def select_candidates(spectrum, library, rt_tolerance, calibration, ms2_calibration):
    candidates = []

    for entry in library:
        # 1. Check isolation window (from mzML)
        if not spectrum.isolation_window.contains(entry.precursor_mz):
            continue

        # 2. Check RT tolerance (with calibration if available)
        if calibration:
            expected_rt = calibration.predict(entry.retention_time)
            # Use local tolerance based on residuals at this RT
            effective_tolerance = calibration.local_tolerance(entry.retention_time)
        else:
            expected_rt = entry.retention_time
            effective_tolerance = rt_tolerance

        if abs(expected_rt - spectrum.retention_time) > effective_tolerance:
            continue

        # 3. Top-3 fragment pre-filter (NEW)
        # Require at least 1 of the top 3 most intense library peaks
        # to be present in the observed spectrum
        if not has_top3_fragment_match(entry, spectrum, ms2_calibration):
            continue

        candidates.append(entry)

    # 4. Limit to max candidates (keep closest in RT)
    if len(candidates) > max_candidates:
        candidates.sort(key=lambda e: abs(e.rt - spectrum.rt))
        candidates = candidates[:max_candidates]

    return candidates
```

### Top-3 Fragment Pre-Filter

The top-3 fragment pre-filter is a fast check to eliminate candidates that have no signal overlap with the observed spectrum before expensive ridge regression.

**Rationale**: If a peptide is truly present, its most intense predicted fragments should be visible. Candidates where none of the top 3 peaks match are very unlikely to contribute signal.

```python
def has_top3_fragment_match(library_entry, spectrum, ms2_calibration):
    """
    Check if at least 1 of the top 3 most intense library peaks
    is present in the observed spectrum.

    Args:
        library_entry: Library spectrum with fragments
        spectrum: Observed spectrum (m/z calibrated)
        ms2_calibration: MS2 calibration for tolerance

    Returns:
        True if at least 1 of top 3 library peaks has a match
    """
    # Get top 3 fragments by intensity
    fragments = sorted(library_entry.fragments,
                       key=lambda f: f.intensity,
                       reverse=True)[:3]

    # Calculate tolerance from calibration (3×SD) or use default
    if ms2_calibration and ms2_calibration.calibrated:
        tolerance = 3.0 * ms2_calibration.sd  # ppm or Th
    else:
        tolerance = default_tolerance

    # Check if any top-3 peak matches an observed peak
    for frag in fragments:
        lib_mz = frag.mz

        # Calculate m/z window
        if unit == "ppm":
            tol_da = lib_mz * tolerance / 1e6
        else:
            tol_da = tolerance

        lower = lib_mz - tol_da
        upper = lib_mz + tol_da

        # Binary search in sorted observed m/z array
        idx = binary_search(spectrum.mzs, lower)
        if idx < len(spectrum.mzs) and spectrum.mzs[idx] <= upper:
            return True  # Found a match

    return False  # No matches
```

**Key points:**

1. **Uses calibrated spectrum**: The observed spectrum m/z values are already shifted by the mean MS2 error
2. **Uses calibrated tolerance**: 3×SD from MS2 calibration ensures tight matching
3. **Fast binary search**: O(log n) per peak, only 3 peaks checked
4. **Conservative filter**: Only rejects candidates with zero top-3 overlap

**Impact**: This filter typically removes 30-50% of candidates that passed RT filtering but have no spectral evidence, significantly reducing the design matrix size for ridge regression.

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

## Non-Negative Constraints and Iterative Solver

Peptide abundances can't be negative. Osprey solves this using **Coordinate Descent Non-Negative Least Squares (CD-NNLS)**, an iterative algorithm that:

1. Enforces non-negativity at each step (no post-hoc clipping)
2. Avoids forming the expensive normal equations (AᵀA)
3. Uses active set acceleration to skip zero coefficients

### The CD-NNLS Algorithm

The algorithm solves:
```
minimize: ‖Ax - b‖² + λ‖x‖²  subject to x ≥ 0
```

**Key insight**: For a single variable xⱼ, holding all other variables fixed, the optimal update has a closed-form solution:

```
xⱼ* = max(0, (aⱼ · residual + xⱼ · ‖aⱼ‖²) / (‖aⱼ‖² + λ))

where:
  aⱼ = column j of design matrix A (the j-th library spectrum)
  residual = b - Ax (what's left unexplained)
  ‖aⱼ‖² = sum of squared intensities in library spectrum j
```

### Detailed Pseudocode

```python
def solve_cd_nnls(A, b, lambda_, max_sweeps=50, tol=1e-6):
    """
    Coordinate Descent NNLS for ridge regression with non-negativity.

    Args:
        A: Design matrix (m bins × k candidates), each column is a binned library spectrum
        b: Observed spectrum (m bins)
        lambda_: Regularization parameter
        max_sweeps: Maximum full passes over all variables
        tol: Convergence tolerance (max coefficient change)

    Returns:
        x: Coefficient vector (k × 1) with all x ≥ 0
    """
    m, k = A.shape  # m bins, k candidates

    # Precompute column norms: ‖aⱼ‖² + λ for each candidate
    col_norm_sq = [sum(A[:, j]**2) + lambda_ for j in range(k)]

    # Initialize: all coefficients zero, residual = observed spectrum
    x = zeros(k)
    residual = b.copy()  # residual = b - Ax = b (since x=0)
    active = [True] * k  # All candidates active initially

    for sweep in range(max_sweeps):
        max_change = 0

        # Use active set after warmup (first 3 sweeps)
        use_active_set = (sweep >= 3)

        for j in range(k):
            # Skip inactive variables (those at zero) after warmup
            if use_active_set and not active[j]:
                continue

            old_xj = x[j]

            # Step 1: Compute correlation between library spectrum j and residual
            # This tells us "how much does candidate j explain the unexplained signal?"
            corr = dot(A[:, j], residual)

            # Step 2: Compute unconstrained optimal update
            # The closed-form solution for coordinate j, holding others fixed:
            #   xⱼ* = (aⱼ·(b - A₋ⱼx₋ⱼ)) / (‖aⱼ‖² + λ)
            #       = (aⱼ·residual + xⱼ·‖aⱼ‖²) / (‖aⱼ‖² + λ)
            col_norm_only = col_norm_sq[j] - lambda_  # ‖aⱼ‖²
            unconstrained = (corr + old_xj * col_norm_only) / col_norm_sq[j]

            # Step 3: Apply non-negativity constraint
            new_xj = max(0, unconstrained)

            # Step 4: Update residual incrementally
            # residual_new = residual - (new_xj - old_xj) × aⱼ
            delta = new_xj - old_xj
            if abs(delta) > 1e-12:
                residual -= delta * A[:, j]  # O(m) operation
                x[j] = new_xj

            max_change = max(max_change, abs(delta))

        # Update active set: candidates at zero become inactive
        if use_active_set:
            for j in range(k):
                active[j] = (x[j] > 0)

            # Every 5 sweeps: check if inactive variables want to become positive
            if sweep % 5 == 0:
                for j in range(k):
                    if not active[j]:
                        # Gradient at xⱼ=0: if negative, variable wants to increase
                        grad = -dot(A[:, j], residual)
                        if grad < -tol:
                            active[j] = True  # Re-activate

        # Check convergence
        if max_change < tol:
            break

    # Final validation: check all inactive variables one more time
    for j in range(k):
        if not active[j]:
            corr = dot(A[:, j], residual)
            if corr / col_norm_sq[j] > tol:
                # This variable should be positive
                new_xj = max(0, corr / col_norm_sq[j])
                residual -= new_xj * A[:, j]
                x[j] = new_xj

    return x
```

### Why Coordinate Descent?

| Approach | Complexity | Notes |
|----------|------------|-------|
| Direct solve (AᵀA + λI)⁻¹Aᵀb | O(mk² + k³) | Must form AᵀA, then factor |
| Coordinate Descent | O(mk) per sweep | Never forms AᵀA, ~10-20 sweeps typical |

For 500 candidates and 1800 bins:
- Direct: 500² × 1800 = 450M ops to form AᵀA, plus O(500³) = 125M to factor
- CD: 500 × 1800 × 15 sweeps = 13.5M ops total

### Active Set Acceleration

The key optimization is **skipping zero coefficients**:

1. **Warmup phase** (sweeps 0-2): Update all candidates to get initial solution
2. **Active set phase** (sweeps 3+): Only update candidates with x > 0
3. **Reactivation checks**: Every 5 sweeps, check if any inactive candidate should become positive (gradient test)
4. **Final validation**: One last check of all inactive candidates

In practice, only 10-50 candidates have non-zero coefficients per spectrum, so active set acceleration provides 10-50× speedup after warmup.

### Convergence

The algorithm converges when the maximum coefficient change in a sweep is below the tolerance:

```
Typical convergence: 10-20 sweeps for tol=1e-6
Fast cases (few active): 5-10 sweeps
Slow cases (many correlated candidates): 30-50 sweeps
```

### Residual Maintenance

The residual `r = b - Ax` is maintained incrementally:
```
When xⱼ changes by Δ:
  r_new = r - Δ × aⱼ
```

This is O(m) per update instead of O(mk) to recompute from scratch.

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
regularization_lambda:
  Fixed: 1.0   # Default lambda value

max_candidates_per_spectrum: 500  # Limit candidates per spectrum
```

### CD-NNLS Solver Parameters

The solver uses these defaults (in `CdNnlsParams`):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_sweeps` | 50 | Maximum full passes over variables |
| `tolerance` | 1e-6 | Convergence threshold (max coefficient change) |
| `relative_tolerance` | 1e-4 | Relative convergence threshold |
| `warmup_sweeps` | 3 | Sweeps before activating active set |
| `reactivation_interval` | 5 | How often to check inactive variables |

## Implementation

Key files:
- `crates/osprey-regression/src/cd_nnls.rs` - **Coordinate Descent NNLS solver** (main algorithm)
- `crates/osprey-regression/src/solver.rs` - Unified f32 solver (unit res + HRAM)
- `crates/osprey-regression/src/binning.rs` - Binner (Comet BIN macro, configurable bin width)
- `crates/osprey-regression/src/matrix.rs` - DesignMatrixBuilder (column-major for BLAS)
- `crates/osprey/src/pipeline.rs` - Integration in analysis pipeline

## Performance Considerations

| Factor | Impact | Mitigation |
|--------|--------|------------|
| Many candidates | More iterations, memory | Limit max_candidates to 500 |
| Active set efficiency | 10-50× speedup | Most candidates converge to zero |
| HRAM bins (100K) | Larger design matrix | On-the-fly binning per spectrum |
| Many spectra | Linear scaling | Parallel processing (rayon) |
| f32 vs f64 | 2× memory, ~1.5× speed | Use f32 for production |
