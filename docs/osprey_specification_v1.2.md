# Software Specification: Peptide-Centric DIA Analysis Tool

## Project Name: **Osprey** (Open-Source Peptide Recognition and Elution Yield)

**Version**: 1.0 Specification  
**Date**: February 2026  
**Authors**: MacCoss Lab, University of Washington  

---

## 1. Executive Summary

Osprey is an open-source software tool for peptide-centric detection and quantification in data-independent acquisition (DIA) mass spectrometry data. It employs regularized regression to deconvolute mixed MS/MS spectra, aggregates evidence across the chromatographic dimension, and uses machine learning to score peptide detections with rigorous FDR control. The tool supports both unit resolution (nominal mass) and high-resolution accurate mass (HRAM) instruments, filling a critical gap in the proteomics software ecosystem where no fully open-source DIA analysis tool currently exists.

### 1.1 Key Objectives

1. Enable peptide-centric DIA analysis on low-cost nominal mass instruments
2. Provide competitive performance with existing closed-source tools (Spectronaut, DIA-NN, CHIMERYS)
3. Deliver a fully open-source, community-extensible codebase
4. Integrate seamlessly with the Skyline targeted proteomics ecosystem
5. Process data efficiently (<10 minutes per 30-minute DIA run)

### 1.2 Scope

- **In Scope**: DIA data analysis, peptide detection, peptide quantification, FDR control, Skyline integration, HRAM and unit resolution support
- **Out of Scope**: Spectral library generation, de novo sequencing, protein inference (handled by downstream tools), instrument control

---

## 2. System Architecture

### 2.1 High-Level Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              Osprey Pipeline                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌───────────┐ │
│  │   Input      │    │  Spectrum    │    │  Chromato-   │    │  Output   │ │
│  │   Module     │───▶│  Regression  │───▶│  graphic     │───▶│  Module   │ │
│  │              │    │  Engine      │    │  Scoring     │    │           │ │
│  └──────────────┘    └──────────────┘    └──────────────┘    └───────────┘ │
│         │                   │                   │                   │       │
│         ▼                   ▼                   ▼                   ▼       │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌───────────┐ │
│  │   Library    │    │   Matrix     │    │   Feature    │    │  Skyline  │ │
│  │   Loader     │    │   Algebra    │    │   Extractor  │    │  Bridge   │ │
│  │ (TSV, blib)  │    │   (BLAS)     │    │              │    │           │ │
│  └──────────────┘    └──────────────┘    └──────────────┘    └───────────┘ │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                        FDR Control Module                            │   │
│  │                    (Target-Decoy / Mokapot)                          │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Library-First Design Philosophy:**

Osprey adopts a library-first design where spectral predictions and retention times are provided via pre-built libraries rather than generated on-the-fly. This approach offers several advantages:

1. **Flexibility**: Users can choose their preferred predictor (Carafe, Prosit, DIA-NN, etc.) or use empirical libraries
2. **Reproducibility**: The exact library used is recorded, enabling precise replication
3. **Efficiency**: Expensive predictions are computed once and reused across analyses
4. **Filtering**: Libraries can be pre-filtered (HRAM-constrained, first-pass, etc.) before analysis
5. **Modularity**: Osprey focuses on detection/quantification; library generation is handled by specialized tools

### 2.2 Module Overview

| Module | Responsibility | Language |
|--------|---------------|----------|
| Input Module | Parse mzML/raw files, manage file I/O | Rust |
| Library Loader | Parse DIA-NN TSV, blib, elib formats | Rust |
| Decoy Generator | Generate reversed/shuffled decoys from library entries | Rust |
| Spectrum Regression Engine | Ridge regression with non-negativity constraints | Rust |
| Matrix Algebra Layer | BLAS operations, sparse matrix support | Rust (ndarray + sprs) |
| Chromatographic Scoring | Peak detection, EMG fitting, feature extraction | Rust |
| Feature Extractor | Compute all scoring features per peptide | Rust |
| FDR Control Module | Target-decoy analysis, mokapot integration | Rust + Python FFI |
| Output Module | Generate results, reports, Skyline-compatible files | Rust |
| Skyline Bridge | Direct integration with Skyline application | C# (.NET) |
| CLI Interface | Command-line interface for standalone use | Rust |

### 2.3 Technology Stack

| Component | Technology | Rationale |
|-----------|------------|-----------|
| Core Engine | Rust 1.75+ | Performance, memory safety, fearless concurrency |
| Linear Algebra | ndarray + ndarray-linalg | Native Rust with BLAS/LAPACK backend |
| Sparse Matrices | sprs | Efficient CSR/CSC for HRAM mode |
| BLAS Backend | OpenBLAS or Intel MKL | Vectorized matrix operations |
| Parallelism | Rayon | Data parallelism across spectra |
| File I/O | mzml-rs, quick-xml | mzML parsing |
| Python Interface | PyO3 | Carafe and mokapot integration |
| Skyline Integration | C# / .NET 6+ | Native Skyline plugin |
| Build System | Cargo | Standard Rust toolchain |
| Testing | cargo test + criterion | Unit tests and benchmarks |
| Documentation | rustdoc + mdBook | API docs and user guide |

---

## 3. Functional Requirements

### 3.1 Input Processing

#### FR-1.1: File Format Support
- **FR-1.1.1**: Parse mzML files (version 1.1+) with indexed and non-indexed variants
- **FR-1.1.2**: Support gzip-compressed mzML (.mzML.gz)
- **FR-1.1.3**: Read Thermo RAW files via ProteoWizard bridge (optional)
- **FR-1.1.4**: Parse centroided and profile-mode spectra (centroid on-the-fly if needed)

#### FR-1.2: Spectral Library Input

Osprey requires a pre-built spectral library containing predicted or empirical spectra and retention times. This design provides flexibility: users can supply a comprehensive library (all peptides from a FASTA prediction), a constrained library (from HRAM-filtered results), or anything in between. Library generation is handled by upstream tools (Skyline, DIA-NN, Carafe, etc.).

**Supported Formats:**

- **FR-1.2.1**: DIA-NN TSV library format
  - Tab-delimited with columns: PrecursorMz, ProductMz, FragmentType, FragmentSeriesNumber, FragmentCharge, FragmentLossType, RelativeIntensity, PeptideSequence, ModifiedPeptide, PrecursorCharge, NormalizedRetentionTime, ProteinId, GeneName, etc.
  - Support both predicted libraries (from DIA-NN --gen-spec-lib) and empirical libraries
  - Handle iRT and normalized RT formats

- **FR-1.2.2**: Skyline blib format (BiblioSpec)
  - SQLite-based spectral library format
  - Contains: peptide sequences, modifications, precursor m/z, fragment m/z and intensities, retention times
  - Support redundant and non-redundant blib files
  - Read via SQLite interface

- **FR-1.2.3**: Skyline elib format (EncyclopeDIA chromatogram libraries)
  - SQLite-based format with chromatographic peak shapes
  - Contains calibrated RT and fragment intensities
  - Preferred format when available (highest quality RT calibration)

**Library Content Requirements:**

- **FR-1.2.4**: Each library entry must contain:
  - Peptide sequence with modification positions and masses
  - Precursor m/z and charge state
  - Fragment m/z values and relative intensities (minimum 3 fragments)
  - Predicted or measured retention time

- **FR-1.2.5**: Support standard modification notation:
  - UniMod identifiers (preferred)
  - Mass delta notation (e.g., M[+15.9949])
  - DIA-NN modification format

**Library Scope Options:**

- **FR-1.2.6**: Comprehensive library: All tryptic peptides from proteome FASTA (predicted)
  - Typical size: 500,000-2,000,000 precursors
  - Use case: Discovery experiments, maximum sensitivity

- **FR-1.2.7**: HRAM-constrained library: Only peptides detected in HRAM DIA of same sample type
  - Typical size: 30,000-100,000 precursors
  - Use case: Unit resolution analysis with improved FDR control

- **FR-1.2.8**: First-pass library: Peptides from initial search (DIA-NN, MSFragger-DIA, etc.)
  - Typical size: Variable based on first-pass results
  - Use case: Iterative refinement, targeted reanalysis

- **FR-1.2.9**: Targeted library: Small curated list of peptides of interest
  - Typical size: 100-10,000 precursors
  - Use case: Biomarker panels, pathway-focused analysis

#### FR-1.3: Decoy Generation
- **FR-1.3.1**: Generate decoys from library entries by sequence reversal (default)
- **FR-1.3.2**: Generate decoys by sequence shuffling (alternative method)
- **FR-1.3.3**: For each target library entry, create corresponding decoy with:
  - Reversed/shuffled sequence
  - Recalculated precursor m/z
  - Fragment m/z values from reversed/shuffled sequence
  - Same relative intensity pattern as target (reordered appropriately)
  - Adjusted retention time based on reversed sequence hydrophobicity (optional)
- **FR-1.3.4**: Maintain 1:1 target-decoy ratio by default (configurable)
- **FR-1.3.5**: Option to use pre-computed decoys if provided in library

### 3.2 Spectrum Processing

#### FR-2.0: Library Format Specifications

**DIA-NN TSV Format:**

The DIA-NN library format is a tab-separated file with the following required columns:

| Column | Type | Description |
|--------|------|-------------|
| PrecursorMz | float | Precursor m/z |
| PrecursorCharge | int | Precursor charge state |
| ModifiedPeptide | string | Peptide sequence with modifications |
| StrippedPeptide | string | Unmodified peptide sequence |
| FragmentMz | float | Fragment m/z |
| RelativeIntensity | float | Normalized fragment intensity (0-1) |
| FragmentType | char | Ion type (b, y, etc.) |
| FragmentSeriesNumber | int | Ion ordinal |
| FragmentCharge | int | Fragment charge |
| FragmentLossType | string | Neutral loss (empty, H2O, NH3, etc.) |
| iRT / NormalizedRetentionTime | float | Retention time |
| ProteinId | string | Protein accession(s) |
| GeneName | string | Gene name(s) |

- **FR-2.0.1**: Parse standard DIA-NN output from `--gen-spec-lib`
- **FR-2.0.2**: Handle both comma and semicolon delimiters for protein/gene lists
- **FR-2.0.3**: Support iRT scale and raw minute RT values
- **FR-2.0.4**: Group fragments by precursor for efficient loading

**Skyline blib Format (BiblioSpec):**

SQLite database with the following key tables:

| Table | Key Columns | Description |
|-------|-------------|-------------|
| RefSpectra | peptideSeq, precursorMZ, precursorCharge, retentionTime | Spectrum metadata |
| RefSpectraPeaks | RefSpectraId, peakMZ, peakIntensity | Fragment peaks |
| Modifications | RefSpectraId, position, mass | Modification positions |
| RefSpectraProteins | RefSpectraId, accession | Protein mappings |

- **FR-2.0.5**: Read blib files via SQLite (rusqlite crate)
- **FR-2.0.6**: Handle both redundant and non-redundant blib formats
- **FR-2.0.7**: Parse modification masses and map to UniMod where possible
- **FR-2.0.8**: Support compressed blib files

**EncyclopeDIA elib Format:**

SQLite database extending blib with chromatographic information:

| Additional Table | Key Columns | Description |
|------------------|-------------|-------------|
| PeptideQuants | peptideSeq, fileName, retentionTime, totalAreaFragment | Quantification |
| PrecursorQuants | precursorId, area, apexRT | Precursor-level quants |

- **FR-2.0.9**: Read elib files with chromatogram library extensions
- **FR-2.0.10**: Use calibrated RT when available (preferred over predicted)
- **FR-2.0.11**: Extract peak shape information for validation

#### FR-2.1: Resolution-Adaptive Binning
- **FR-2.1.1**: Implement Comet-style binning with configurable bin width
- **FR-2.1.2**: Unit resolution mode: 1.0005 Th bins (~1,000 bins for typical fragment range)
- **FR-2.1.3**: HRAM mode: 0.02-0.05 Th bins (~10,000-50,000 bins)
- **FR-2.1.4**: Auto-detect resolution from data if not specified
- **FR-2.1.5**: Apply bin offset for optimal peak alignment

#### FR-2.2: Candidate Selection
- **FR-2.2.1**: Filter peptides by precursor m/z within isolation window (± tolerance)
- **FR-2.2.2**: Filter peptides by predicted RT within configurable window (default ±2 min)
- **FR-2.2.3**: Consider isotope envelope overlap with isolation window
- **FR-2.2.4**: Limit maximum candidates per spectrum (configurable, default 200)

#### FR-2.3: Design Matrix Construction
- **FR-2.3.1**: Build matrix A where columns are predicted spectra for candidates
- **FR-2.3.2**: Normalize each column to unit total intensity
- **FR-2.3.3**: Use dense matrix for unit resolution (ndarray)
- **FR-2.3.4**: Use sparse matrix (CSR) for HRAM mode (sprs)
- **FR-2.3.5**: Bin predicted spectra according to resolution mode

### 3.3 Regression Engine

#### FR-3.1: Ridge Regression with Non-Negativity
- **FR-3.1.1**: Solve: minimize ‖Ax - b‖² + λ‖x‖² subject to x ≥ 0
- **FR-3.1.2**: Compute (AᵀA + λI) via Cholesky factorization
- **FR-3.1.3**: Apply non-negativity via active-set method or coordinate descent
- **FR-3.1.4**: Warm-start from previous scan's solution
- **FR-3.1.5**: Support configurable regularization parameter λ

#### FR-3.2: Regularization Selection
- **FR-3.2.1**: Default λ selection via cross-validation on subset of spectra
- **FR-3.2.2**: Allow user-specified fixed λ
- **FR-3.2.3**: Support adaptive λ per spectrum based on candidate count

#### FR-3.3: Coefficient Output
- **FR-3.3.1**: Store coefficient vector for each spectrum
- **FR-3.3.2**: Associate coefficients with peptide identifiers and scan numbers
- **FR-3.3.3**: Store in memory-efficient columnar format
- **FR-3.3.4**: Support streaming to disk for very large files

### 3.4 Chromatographic Analysis

#### FR-4.1: Coefficient Time Series Extraction
- **FR-4.1.1**: For each peptide, extract coefficient values across all relevant scans
- **FR-4.1.2**: Associate with retention time for each scan
- **FR-4.1.3**: Handle gaps in scan coverage (DIA window cycling)

#### FR-4.2: Peak Detection
- **FR-4.2.1**: Identify contiguous regions with non-zero coefficients
- **FR-4.2.2**: Find local maxima (apex candidates)
- **FR-4.2.3**: Determine peak boundaries (start/end RT)
- **FR-4.2.4**: Handle multiple peaks per peptide (report best or all)

#### FR-4.3: Peak Shape Modeling
- **FR-4.3.1**: Fit exponentially modified Gaussian (EMG) to coefficient time series
- **FR-4.3.2**: Extract EMG parameters: μ (center), σ (width), τ (tailing)
- **FR-4.3.3**: Calculate fit quality (residual sum of squares)
- **FR-4.3.4**: Calculate peak metrics: FWHM, symmetry, prominence

#### FR-4.4: Background Estimation
- **FR-4.4.1**: Compute features at off-target RT windows (e.g., predicted RT ± 5, ±8 min)
- **FR-4.4.2**: Calculate median/mean background for each feature
- **FR-4.4.3**: Compute background-subtracted feature values
- **FR-4.4.4**: Compute contrast ratios and z-scores

### 3.5 Feature Extraction

#### FR-5.1: Chromatographic Features
- **FR-5.1.1**: Peak apex coefficient (maximum value)
- **FR-5.1.2**: Integrated peak area (AUC of coefficients)
- **FR-5.1.3**: EMG fit quality (RSS)
- **FR-5.1.4**: Peak width (FWHM)
- **FR-5.1.5**: Peak symmetry (leading/trailing ratio)
- **FR-5.1.6**: RT deviation from prediction
- **FR-5.1.7**: Normalized RT deviation (deviation / uncertainty)
- **FR-5.1.8**: Number of contributing scans
- **FR-5.1.9**: Coefficient stability (variance near apex)
- **FR-5.1.10**: Peak boundary sharpness
- **FR-5.1.11**: Peak prominence (apex / baseline)

#### FR-5.2: Spectral Features
- **FR-5.2.1**: Aggregate observed spectrum across peak apex region
- **FR-5.2.2**: Hyperscore (X!Tandem style)
- **FR-5.2.3**: Normalized spectral contrast angle
- **FR-5.2.4**: Dot product
- **FR-5.2.5**: Pearson intensity correlation
- **FR-5.2.6**: Spearman rank correlation
- **FR-5.2.7**: Fragment coverage (fraction detected)
- **FR-5.2.8**: Sequence coverage (backbone coverage)
- **FR-5.2.9**: Consecutive ion count (longest b/y run)
- **FR-5.2.10**: Base peak match rank
- **FR-5.2.11**: Top-3 match count
- **FR-5.2.12**: Explained intensity fraction

#### FR-5.3: Contextual Features
- **FR-5.3.1**: Number of competing candidates
- **FR-5.3.2**: Relative coefficient magnitude
- **FR-5.3.3**: Local peptide density
- **FR-5.3.4**: Spectral complexity estimate
- **FR-5.3.5**: Regression residual
- **FR-5.3.6**: Precursor intensity (if MS1 available)
- **FR-5.3.7**: Modification count

#### FR-5.4: Background-Corrected Features
- **FR-5.4.1**: For each feature F: F_background, F_subtracted, F_contrast, F_zscore
- **FR-5.4.2**: Configurable off-target RT offsets

### 3.6 FDR Control

Osprey implements a hierarchical FDR control strategy with both run-level and experiment-level q-values, enabling appropriate confidence assessment at different scopes.

#### FR-6.1: Target-Decoy Framework
- **FR-6.1.1**: Compute all features for both targets and decoys
- **FR-6.1.2**: Apply identical pipeline to targets and decoys
- **FR-6.1.3**: Support peptide-level FDR (not spectrum-level)
- **FR-6.1.4**: Decoys compete fairly with targets at all stages

#### FR-6.2: Machine Learning Scoring
- **FR-6.2.1**: Interface with mokapot for semi-supervised learning
- **FR-6.2.2**: Train SVM to discriminate targets from decoys
- **FR-6.2.3**: Output discriminant score per peptide
- **FR-6.2.4**: Support Percolator as alternative backend
- **FR-6.2.5**: Allow pre-trained models for consistency across experiments

#### FR-6.3: Run-Level FDR Control
- **FR-6.3.1**: Compute q-values independently within each run
- **FR-6.3.2**: Target-decoy competition uses only peptides searched in that run
- **FR-6.3.3**: Mokapot training can be per-run or shared across runs
- **FR-6.3.4**: Report run-level q-value and posterior error probability (PEP)
- **FR-6.3.5**: Run-level FDR answers: "Is this peptide confidently detected in this specific file?"

#### FR-6.4: Experiment-Level FDR Control
- **FR-6.4.1**: Compute q-values across all runs in the experiment
- **FR-6.4.2**: Aggregate evidence: peptides detected in multiple runs have increased confidence
- **FR-6.4.3**: Account for multiple testing burden across runs
- **FR-6.4.4**: Experiment-level FDR answers: "Is this peptide confidently detected somewhere in this experiment?"

**Experiment-Level Q-Value Calculation:**

Two approaches are supported:

*Approach A: Best-run q-value with correction*
- Take the best (lowest) run-level q-value for each peptide
- Apply Benjamini-Hochberg correction for number of runs searched
- Simple but may be conservative

*Approach B: Combined evidence scoring*
- Combine features across runs (e.g., max score, sum of scores, # runs detected)
- Train experiment-level model on combined features
- More powerful but requires more computation

- **FR-6.4.5**: Default to Approach B when multiple runs available
- **FR-6.4.6**: Fall back to Approach A for single-run experiments

#### FR-6.5: Two-Step Search FDR Strategy

The two-step search (FR-8.6) interacts with FDR control as follows:

**Step 1 FDR:**
- Run-level FDR at liberal threshold (e.g., 1%) to maximize discovery
- Union of detected peptides forms constrained library

**Step 2 FDR:**
- Constrained search space dramatically reduces multiple testing burden
- Run-level and experiment-level q-values computed on Step 2 results
- Final reported q-values come from Step 2

- **FR-6.5.1**: Step 1 q-values are intermediate, not reported to user
- **FR-6.5.2**: Step 2 q-values are final and reported in output
- **FR-6.5.3**: Document which peptides were detected only in Step 2 (not Step 1)

#### FR-6.6: FDR Thresholding and Reporting
- **FR-6.6.1**: Filter peptides at user-specified q-value threshold (default 0.01)
- **FR-6.6.2**: Report number of detections at multiple FDR levels (0.1%, 1%, 5%, 10%)
- **FR-6.6.3**: Separate reporting for run-level and experiment-level
- **FR-6.6.4**: Option to require detection in minimum N runs at experiment level

### 3.7 Peak Detection (for Skyline Quantification)

Osprey does **not** perform quantification—that is Skyline's role. Osprey's responsibility is to determine **where** peptides elute so Skyline can quantify them accurately.

#### FR-7.1: Peak Boundary Determination
- **FR-7.1.1**: Identify peak boundaries from regression coefficient time series
- **FR-7.1.2**: Apply Savitzky-Golay smoothing to reduce noise effects on boundaries
- **FR-7.1.3**: Define start boundary as first point where coefficient exceeds 10% of apex
- **FR-7.1.4**: Define end boundary as last point where coefficient exceeds 10% of apex
- **FR-7.1.5**: Ensure boundaries span at least 3 scans (configurable minimum)
- **FR-7.1.6**: Handle asymmetric peaks (tailing) appropriately

#### FR-7.2: Peak Quality Assessment
- **FR-7.2.1**: Flag peaks with poor EMG fit (potential interference)
- **FR-7.2.2**: Flag peaks with abnormal width (too narrow or too wide)
- **FR-7.2.3**: Flag split peaks or shoulders
- **FR-7.2.4**: Include quality flags in output for Skyline to display

#### FR-7.3: Apex Determination
- **FR-7.3.1**: Report apex as time of maximum coefficient
- **FR-7.3.2**: Use parabolic interpolation for sub-scan apex precision
- **FR-7.3.3**: Compare apex to predicted RT as quality metric

#### FR-7.4: Coefficient Export for Visualization
- **FR-7.4.1**: Optionally export full coefficient time series
- **FR-7.4.2**: Skyline can display as "Osprey coefficient chromatogram"
- **FR-7.4.3**: Useful for understanding interference and peak quality
- **FR-7.4.4**: Compress time series to reduce file size (store only non-zero regions)

### 3.8 Output Generation

Osprey's primary output is a blib-format SQLite file containing detection results, peak boundaries, and q-values that Skyline uses to perform quantification. Osprey focuses on **detection** (is this peptide present?) and **localization** (where are the peak boundaries?), while Skyline handles **quantification** (how much is present?).

#### FR-8.1: Primary Output - Results blib File

The output blib extends the standard BiblioSpec format with Osprey-specific tables:

**Standard blib Tables (populated):**

| Table | Contents |
|-------|----------|
| RefSpectra | Detected peptide precursors with consensus spectrum |
| RefSpectraPeaks | Fragment m/z and intensities (aggregated across peak) |
| Modifications | Modification positions and masses |
| RefSpectraProteins | Protein mappings |

**Osprey Extension Tables:**

| Table | Columns | Description |
|-------|---------|-------------|
| OspreyPeakBoundaries | RefSpectraId, FileName, StartRT, EndRT, ApexRT | Peak boundaries per run |
| OspreyRunScores | RefSpectraId, FileName, RunQValue, DiscriminantScore, FeatureVector | Run-level detection scores |
| OspreyExperimentScores | RefSpectraId, ExperimentQValue, NRunsDetected, NRunsSearched | Experiment-level scores |
| OspreyCoefficients | RefSpectraId, FileName, ScanNumber, RT, Coefficient | Optional: full coefficient time series |
| OspreyMetadata | Key, Value | Analysis parameters, version info |

- **FR-8.1.1**: Output SQLite blib file readable by Skyline
- **FR-8.1.2**: Include peak boundaries (StartRT, EndRT, ApexRT) for each peptide in each run
- **FR-8.1.3**: Include run-level q-values for per-file FDR control
- **FR-8.1.4**: Include experiment-level q-values for cross-run FDR control
- **FR-8.1.5**: Include discriminant scores and feature vectors for transparency
- **FR-8.1.6**: Optionally include full coefficient time series for visualization

#### FR-8.2: Peak Boundary Determination

- **FR-8.2.1**: Determine peak boundaries from regression coefficient time series
- **FR-8.2.2**: StartRT = first scan where coefficient exceeds baseline threshold
- **FR-8.2.3**: EndRT = last scan where coefficient exceeds baseline threshold
- **FR-8.2.4**: ApexRT = scan with maximum coefficient value
- **FR-8.2.5**: Apply smoothing to avoid boundary jitter from noise
- **FR-8.2.6**: Enforce minimum peak width (configurable, default 3 scans)
- **FR-8.2.7**: Handle split peaks (report primary peak, flag secondary)

#### FR-8.3: Run-Level Q-Values

- **FR-8.3.1**: Compute q-values independently for each run
- **FR-8.3.2**: Use target-decoy competition within each file
- **FR-8.3.3**: Apply mokapot/Percolator scoring per run
- **FR-8.3.4**: Report both q-value and posterior error probability (PEP)

#### FR-8.4: Experiment-Level Q-Values

- **FR-8.4.1**: Compute experiment-level q-values across all runs in batch
- **FR-8.4.2**: Peptides detected in multiple runs receive boosted confidence
- **FR-8.4.3**: Account for multiple testing across runs
- **FR-8.4.4**: Use two-step search strategy for improved experiment-level FDR (see FR-8.6)

#### FR-8.5: Additional Output Files

- **FR-8.5.1**: Tab-delimited peptide report (human-readable summary)
- **FR-8.5.2**: Tab-delimited run-level report with all features
- **FR-8.5.3**: QC metrics report (# detections, score distributions, etc.)
- **FR-8.5.4**: Log file with analysis parameters and warnings

#### FR-8.6: Two-Step Search Strategy

Osprey implements a two-step search strategy similar to Spectronaut and DIA-NN to improve detection consistency and experiment-level FDR control:

**Step 1: Discovery Search (per-file)**
- **FR-8.6.1**: Search each file independently against the full input library
- **FR-8.6.2**: Apply run-level FDR threshold (default 1%)
- **FR-8.6.3**: Collect union of all peptides detected across any run
- **FR-8.6.4**: This produces a "detected peptide set" for the experiment

**Step 2: Constrained Search (experiment-wide)**
- **FR-8.6.5**: Create constrained library containing only peptides from Step 1
- **FR-8.6.6**: Re-search all files against constrained library
- **FR-8.6.7**: Reduced search space improves statistical power
- **FR-8.6.8**: Enables detection of peptides missed in Step 1 due to stringent FDR
- **FR-8.6.9**: Compute experiment-level q-values on Step 2 results

**Two-Step Configuration:**
- **FR-8.6.10**: Enable/disable two-step search (default: enabled for multi-file experiments)
- **FR-8.6.11**: Configurable Step 1 FDR threshold (default: 1%)
- **FR-8.6.12**: Configurable minimum runs for peptide inclusion in Step 2 (default: 1)
- **FR-8.6.13**: Option to include "high-confidence" peptides from library even if not detected in Step 1

### 3.9 Skyline Integration

Osprey and Skyline have complementary roles: Osprey performs peptide **detection** and **peak boundary determination**, while Skyline performs **quantification** and provides the user interface for review and refinement.

#### FR-9.1: Workflow Integration

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│   Library   │────▶│   Osprey    │────▶│  Results    │────▶│   Skyline   │
│  (blib/tsv) │     │  Detection  │     │   (blib)    │     │   Quant     │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
                          │                    │
                          ▼                    ▼
                    ┌───────────┐      ┌─────────────────┐
                    │  mzML     │      │ • Peak bounds   │
                    │  files    │      │ • Run q-values  │
                    └───────────┘      │ • Exp q-values  │
                                       │ • Coefficients  │
                                       └─────────────────┘
```

- **FR-9.1.1**: Skyline imports Osprey results blib directly
- **FR-9.1.2**: Peak boundaries from Osprey guide Skyline's integration windows
- **FR-9.1.3**: Q-values from Osprey inform Skyline's display and filtering
- **FR-9.1.4**: Skyline performs fragment-level quantification using its existing algorithms
- **FR-9.1.5**: Users can manually adjust boundaries in Skyline; Osprey's boundaries are starting points

#### FR-9.2: Skyline Plugin (Optional)

- **FR-9.2.1**: Implement as Skyline external tool for seamless invocation
- **FR-9.2.2**: Accept Skyline document as input (extracts library from document)
- **FR-9.2.3**: Automatically import results blib back into Skyline document
- **FR-9.2.4**: Display Osprey coefficient chromatograms alongside standard XICs
- **FR-9.2.5**: Show run-level and experiment-level q-values in Skyline UI
- **FR-9.2.6**: Color-code peptides by detection confidence

#### FR-9.3: Results blib Compatibility

- **FR-9.3.1**: Output blib must be importable by Skyline without modification
- **FR-9.3.2**: Extension tables (OspreyPeakBoundaries, etc.) are ignored by Skyline if not recognized
- **FR-9.3.3**: Core blib tables follow BiblioSpec schema exactly
- **FR-9.3.4**: Test compatibility with Skyline versions 21.1+

#### FR-9.4: Peak Boundary Communication

The key interface between Osprey and Skyline is the peak boundary information:

| Column | Type | Description |
|--------|------|-------------|
| RefSpectraId | int | Links to peptide in RefSpectra table |
| FileName | string | Run identifier (matches mzML filename) |
| StartRT | float | Peak start time (minutes) |
| EndRT | float | Peak end time (minutes) |
| ApexRT | float | Peak apex time (minutes) |
| RunQValue | float | Detection confidence for this run |
| ExperimentQValue | float | Detection confidence across experiment |

- **FR-9.4.1**: Skyline reads OspreyPeakBoundaries table to set initial integration windows
- **FR-9.4.2**: Boundaries are suggestions; Skyline's peak detection can refine them
- **FR-9.4.3**: Q-values can be used for filtering in Skyline's Document Grid

---

## 4. Non-Functional Requirements

### 4.1 Performance

| Requirement | Target | Measurement |
|-------------|--------|-------------|
| NFR-1.1 | Process 30-min DIA run in <10 min | Wall clock time on reference hardware |
| NFR-1.2 | Memory usage <16 GB for typical run | Peak RSS |
| NFR-1.3 | Linear scaling with file size | Time vs spectra count |
| NFR-1.4 | Efficient parallelization | >80% CPU utilization on 16 cores |
| NFR-1.5 | HRAM mode <2x slower than unit resolution | Relative processing time |

**Reference Hardware**: 16-core CPU (AMD Ryzen 9 or Intel i9), 32 GB RAM, NVMe SSD

### 4.2 Scalability

- **NFR-2.1**: Handle files up to 50 GB (extended DIA runs)
- **NFR-2.2**: Support peptide lists up to 500,000 entries
- **NFR-2.3**: Process batch of 100+ files without restart
- **NFR-2.4**: Streaming mode for memory-constrained environments

### 4.3 Reliability

- **NFR-3.1**: Graceful handling of malformed input files
- **NFR-3.2**: Detailed error messages with actionable guidance
- **NFR-3.3**: Checkpoint/resume for long-running analyses
- **NFR-3.4**: Deterministic results (same input → same output)

### 4.4 Portability

- **NFR-4.1**: Run on Linux (primary), macOS, Windows
- **NFR-4.2**: No proprietary dependencies in core engine
- **NFR-4.3**: Containerized deployment option (Docker)
- **NFR-4.4**: Conda package for easy installation

### 4.5 Maintainability

- **NFR-5.1**: Modular architecture with clear interfaces
- **NFR-5.2**: >80% code coverage in unit tests
- **NFR-5.3**: Comprehensive API documentation
- **NFR-5.4**: Semantic versioning for releases

### 4.6 Usability

- **NFR-6.1**: Sensible defaults requiring minimal configuration
- **NFR-6.2**: Clear progress indication during processing
- **NFR-6.3**: Example datasets and tutorials
- **NFR-6.4**: Active user support channels (GitHub issues, forums)

---

## 5. Data Structures

### 5.1 Core Types

```rust
/// Library entry representing a peptide precursor with spectral information
pub struct LibraryEntry {
    pub id: u32,                          // Unique identifier
    pub sequence: String,                  // Unmodified amino acid sequence
    pub modified_sequence: String,         // Sequence with modification notation
    pub modifications: Vec<Modification>,
    pub charge: u8,
    pub precursor_mz: f64,
    pub retention_time: f64,               // Normalized or calibrated RT
    pub rt_calibrated: bool,               // Whether RT is calibrated to this run
    pub fragments: Vec<LibraryFragment>,
    pub protein_ids: Vec<String>,
    pub gene_names: Vec<String>,
    pub is_decoy: bool,
}

/// Fragment ion from library
pub struct LibraryFragment {
    pub mz: f64,
    pub relative_intensity: f32,           // Normalized 0-1 or 0-100
    pub annotation: FragmentAnnotation,
}

/// Modification on a peptide
pub struct Modification {
    pub position: usize,          // 0-indexed position in sequence
    pub unimod_id: Option<u32>,   // UniMod accession if known
    pub mass_delta: f64,          // Monoisotopic mass change
    pub name: Option<String>,     // Human-readable name (e.g., "Oxidation")
}

/// Fragment ion annotation
pub struct FragmentAnnotation {
    pub ion_type: IonType,        // b, y, etc.
    pub ordinal: u8,              // ion number
    pub charge: u8,
    pub neutral_loss: Option<NeutralLoss>,
}

/// Ion types
pub enum IonType {
    B, Y, A, C, X, Z,
    Precursor,
    Internal,
    Immonium,
    Unknown,
}

/// Neutral losses
pub enum NeutralLoss {
    H2O,      // -18.0106
    NH3,      // -17.0265
    H3PO4,    // -97.9769 (phospho)
    Custom(f64),
}

/// MS/MS spectrum from data file
pub struct Spectrum {
    pub scan_number: u32,
    pub retention_time: f64,
    pub precursor_mz: f64,
    pub isolation_window: IsolationWindow,
    pub mzs: Vec<f64>,
    pub intensities: Vec<f32>,
}

/// DIA isolation window
pub struct IsolationWindow {
    pub center: f64,
    pub lower_offset: f64,
    pub upper_offset: f64,
}

/// Binned spectrum for regression
pub struct BinnedSpectrum {
    pub bin_indices: Vec<u32>,    // Which bins have signal
    pub intensities: Vec<f32>,    // Intensity in each bin
}

/// Regression result for one spectrum
pub struct RegressionResult {
    pub scan_number: u32,
    pub retention_time: f64,
    pub library_ids: Vec<u32>,    // Indices into library
    pub coefficients: Vec<f64>,   // Corresponding coefficients
    pub residual: f64,            // Unexplained intensity
}

/// Peptide detection result
pub struct PeptideDetection {
    pub library_entry: u32,       // Reference to library entry
    pub run_results: Vec<RunDetection>,  // Results per run
    pub experiment_q_value: f64,  // Experiment-level q-value
    pub n_runs_detected: u32,     // Number of runs with detection
    pub n_runs_searched: u32,     // Total runs searched
    pub detected_step1: bool,     // Was this detected in Step 1?
}

/// Per-run detection result
pub struct RunDetection {
    pub file_name: String,
    pub detected: bool,
    pub peak_boundaries: Option<PeakBoundaries>,
    pub run_q_value: f64,
    pub discriminant_score: f64,
    pub features: FeatureSet,
    pub background_features: FeatureSet,
}

/// Peak boundary information for Skyline
pub struct PeakBoundaries {
    pub start_rt: f64,            // Peak start time (minutes)
    pub end_rt: f64,              // Peak end time (minutes)  
    pub apex_rt: f64,             // Peak apex time (minutes)
    pub apex_coefficient: f64,    // Coefficient at apex
    pub integrated_area: f64,     // Sum of coefficients (for QC, not quant)
    pub peak_quality: PeakQuality,
}

/// Peak quality flags
pub struct PeakQuality {
    pub emg_fit_r2: f64,          // R² of EMG fit
    pub is_split: bool,           // Peak appears split
    pub is_truncated: bool,       // Peak truncated at gradient edge
    pub has_shoulder: bool,       // Shoulder detected
    pub width_percentile: f64,    // Width relative to other peaks (0-100)
}

/// Complete feature set for scoring
pub struct FeatureSet {
    // Chromatographic features
    pub peak_apex: f64,
    pub peak_area: f64,
    pub emg_fit_quality: f64,
    pub peak_width: f64,
    pub peak_symmetry: f64,
    pub rt_deviation: f64,
    pub rt_deviation_normalized: f64,
    pub n_contributing_scans: u32,
    pub coefficient_stability: f64,
    pub peak_sharpness: f64,
    pub peak_prominence: f64,
    
    // Spectral features
    pub hyperscore: f64,
    pub spectral_contrast_angle: f64,
    pub dot_product: f64,
    pub pearson_correlation: f64,
    pub spearman_correlation: f64,
    pub fragment_coverage: f64,
    pub sequence_coverage: f64,
    pub consecutive_ions: u32,
    pub base_peak_rank: u32,
    pub top3_matches: u32,
    pub explained_intensity: f64,
    
    // Contextual features
    pub n_competitors: u32,
    pub relative_coefficient: f64,
    pub local_peptide_density: f64,
    pub spectral_complexity: f64,
    pub regression_residual: f64,
    pub precursor_intensity: Option<f64>,
    pub modification_count: u32,
}
```

### 5.2 Matrix Types

```rust
/// Dense design matrix (unit resolution)
pub type DenseDesignMatrix = ndarray::Array2<f64>;

/// Sparse design matrix (HRAM mode)
pub type SparseDesignMatrix = sprs::CsMat<f64>;

/// Resolution-agnostic matrix wrapper
pub enum DesignMatrix {
    Dense(DenseDesignMatrix),
    Sparse(SparseDesignMatrix),
}

/// Binning configuration
pub struct BinConfig {
    pub bin_width: f64,           // Th
    pub bin_offset: f64,          // Th
    pub min_mz: f64,
    pub max_mz: f64,
    pub n_bins: usize,
}

impl BinConfig {
    pub fn unit_resolution() -> Self {
        Self {
            bin_width: 1.0005079,
            bin_offset: 0.4,
            min_mz: 100.0,
            max_mz: 2000.0,
            n_bins: 1899,
        }
    }
    
    pub fn hram(tolerance_th: f64) -> Self {
        Self {
            bin_width: tolerance_th,
            bin_offset: 0.0,
            min_mz: 100.0,
            max_mz: 2000.0,
            n_bins: ((2000.0 - 100.0) / tolerance_th).ceil() as usize,
        }
    }
    
    pub fn mz_to_bin(&self, mz: f64) -> usize {
        ((mz - self.min_mz + self.bin_offset) / self.bin_width).floor() as usize
    }
}
```

### 5.3 Configuration

```rust
/// Main configuration structure
pub struct OspreyConfig {
    // Input/Output
    pub input_files: Vec<PathBuf>,
    pub library_source: LibrarySource,
    pub output_blib: PathBuf,         // Primary output: blib for Skyline
    pub output_report: Option<PathBuf>, // Optional: TSV report
    
    // Resolution settings
    pub resolution_mode: ResolutionMode,
    pub custom_bin_width: Option<f64>,
    
    // Candidate selection
    pub rt_tolerance: f64,            // minutes, default 2.0
    pub precursor_tolerance: f64,     // Th, default 0.5
    pub max_candidates_per_spectrum: usize,  // default 200
    
    // Regression
    pub regularization_lambda: RegularizationSetting,
    pub max_iterations: usize,        // default 1000
    pub convergence_threshold: f64,   // default 1e-6
    
    // Background correction
    pub background_rt_offsets: Vec<f64>,  // default [5.0, 8.0]
    
    // Two-step search
    pub two_step_search: TwoStepConfig,
    
    // FDR control
    pub run_fdr: f64,                 // default 0.01
    pub experiment_fdr: f64,          // default 0.01
    pub decoy_method: DecoyMethod,
    pub decoys_in_library: bool,      // true if library already contains decoys
    
    // Performance
    pub n_threads: usize,             // default: all cores
    pub memory_limit_gb: Option<f64>,
    
    // Output options
    pub export_coefficients: bool,    // Include coefficient time series in blib
    pub export_features: bool,        // Write features to separate TSV
}

/// Two-step search configuration
pub struct TwoStepConfig {
    pub enabled: bool,                // default: true for multi-file
    pub step1_fdr: f64,               // default: 0.01
    pub min_runs_for_step2: usize,    // default: 1
    pub include_high_confidence: bool, // Include library peptides even if not detected
}

impl TwoStepConfig {
    pub fn enabled() -> Self {
        Self {
            enabled: true,
            step1_fdr: 0.01,
            min_runs_for_step2: 1,
            include_high_confidence: false,
        }
    }
    
    pub fn disabled() -> Self {
        Self {
            enabled: false,
            ..Self::enabled()
        }
    }
}

pub enum ResolutionMode {
    UnitResolution,
    HRAM { tolerance_ppm: f64 },
    Auto,
}

pub enum RegularizationSetting {
    Fixed(f64),
    CrossValidated,
    Adaptive,
}

pub enum DecoyMethod {
    Reverse,
    Shuffle,
    FromLibrary,  // Use decoys already in library
}

/// Spectral library source
pub enum LibrarySource {
    /// DIA-NN TSV format library
    DiannTsv(PathBuf),
    /// Skyline BiblioSpec library (.blib)
    Blib(PathBuf),
    /// EncyclopeDIA chromatogram library (.elib)
    Elib(PathBuf),
    /// Skyline document (extracts library from document)
    SkylineDocument(PathBuf),
}

/// Library metadata
pub struct LibraryInfo {
    pub source: LibrarySource,
    pub n_targets: usize,
    pub n_decoys: usize,
    pub rt_type: RtType,
    pub modifications: Vec<String>,    // Unique modifications in library
    pub proteome_coverage: Option<f64>, // If known
}

pub enum RtType {
    Predicted,           // From predictor, not calibrated
    Normalized,          // iRT or similar normalized scale
    CalibratedMinutes,   // Calibrated to specific LC conditions
}
```

---

## 6. Algorithm Details

### 6.1 Ridge Regression with Non-Negativity

**Input**: Design matrix A (m bins × k candidates), observed spectrum b (m × 1), regularization λ

**Output**: Coefficient vector x (k × 1) with x ≥ 0

```
Algorithm: Non-Negative Ridge Regression (Active Set Method)

1. Initialize:
   - Compute G = AᵀA + λI  (k × k matrix)
   - Compute c = Aᵀb       (k × 1 vector)
   - Perform Cholesky factorization: G = LLᵀ
   - Initialize x = 0, active_set = {}

2. Iterate until convergence:
   a. Compute gradient: g = Gx - c
   
   b. If all(g[i] ≥ 0 for i not in active_set) and all(x[i] ≥ 0):
      - Converged, return x
   
   c. Find most violated constraint:
      - j = argmin(g[i]) for i not in active_set with g[i] < 0
      - Add j to active_set
   
   d. Solve unconstrained subproblem on active set:
      - x_active = solve(G[active, active], c[active]) using Cholesky
      - x[not active] = 0
   
   e. If any x[active] < 0:
      - Find limiting constraint via line search
      - Remove limiting index from active_set
      - Go to step 2d
   
3. Return x
```

**Warm-Start Optimization**:
- Initialize active_set from previous scan's solution
- Skip Cholesky factorization if G is unchanged (same candidates)
- Typically converges in 1-3 iterations with warm start

### 6.2 EMG Peak Fitting

**Input**: Coefficient time series {(t_i, x_i)} for a peptide

**Output**: EMG parameters (μ, σ, τ, A) and fit quality

```
Algorithm: EMG Fitting via Levenberg-Marquardt

1. Initialize parameters:
   - μ = weighted mean of t_i by x_i
   - σ = weighted std of t_i
   - τ = σ / 2 (typical tailing)
   - A = max(x_i)

2. Define EMG function:
   EMG(t; μ, σ, τ, A) = (A / 2τ) * exp((σ²/2τ²) + (μ-t)/τ) * erfc((μ + σ²/τ - t) / (√2 * σ))

3. Minimize via Levenberg-Marquardt:
   minimize Σ(x_i - EMG(t_i; μ, σ, τ, A))²

4. Compute fit quality:
   RSS = Σ(x_i - EMG(t_i))²
   R² = 1 - RSS / Σ(x_i - mean(x))²

5. Return (μ, σ, τ, A, RSS, R²)
```

### 6.3 Hyperscore Calculation

**Input**: Predicted spectrum P, observed (deconvoluted) spectrum O

**Output**: Hyperscore

```
Algorithm: Hyperscore (X!Tandem style)

1. Match predicted fragments to observed:
   - For each predicted fragment f in P:
     - Find closest observed peak within tolerance
     - Record matched intensity I_f

2. Count matched ions by type:
   - n_b = count of matched b-ions
   - n_y = count of matched y-ions

3. Calculate score:
   hyperscore = log(n_b!) + log(n_y!) + Σ log(I_f + 1)
   
   where sum is over all matched fragments

4. Return hyperscore
```

### 6.4 Two-Step Search Algorithm

The two-step search improves detection consistency across runs and provides better experiment-level FDR control by reducing the search space in the second step.

```
Algorithm: Two-Step Search

INPUT:
  - input_files: List of mzML files [F1, F2, ..., Fn]
  - library: Full spectral library L with targets and decoys
  - config: Search parameters

OUTPUT:
  - results: PeptideDetection for each library entry with run and experiment q-values

STEP 1: DISCOVERY SEARCH
========================
detected_peptides = empty set

FOR each file Fi in input_files:
    # Search against full library
    run_results[Fi] = search_file(Fi, library, config)
    
    # Apply run-level FDR
    run_detections = filter_by_qvalue(run_results[Fi], config.step1_fdr)
    
    # Add to detected set
    FOR each peptide P in run_detections:
        detected_peptides.add(P)
    
    LOG: "Step 1: File {Fi} - {count} peptides at {step1_fdr} FDR"

LOG: "Step 1 complete: {detected_peptides.size()} unique peptides detected across {n} files"

STEP 2: CONSTRAINED SEARCH  
==========================
# Build constrained library
constrained_library = filter_library(library, detected_peptides)
constrained_library += generate_decoys(constrained_library)  # Fresh decoys for Step 2

LOG: "Step 2: Searching against {constrained_library.size()} peptides (reduced from {library.size()})"

FOR each file Fi in input_files:
    # Search against constrained library
    run_results[Fi] = search_file(Fi, constrained_library, config)
    
    # Compute run-level q-values (more power due to smaller search space)
    run_qvalues[Fi] = compute_run_fdr(run_results[Fi])

# Compute experiment-level q-values
experiment_qvalues = compute_experiment_fdr(run_results, run_qvalues)

# Assemble final results
FOR each peptide P in constrained_library.targets:
    detection = PeptideDetection {
        library_entry: P.id,
        experiment_q_value: experiment_qvalues[P],
        n_runs_detected: count_detections(P, run_qvalues, config.run_fdr),
        detected_step1: P in detected_peptides,
        run_results: [
            RunDetection {
                file_name: Fi,
                run_q_value: run_qvalues[Fi][P],
                peak_boundaries: extract_boundaries(run_results[Fi][P]),
                ...
            } for Fi in input_files
        ]
    }
    results.add(detection)

RETURN results
```

**Key Benefits of Two-Step Search:**

1. **Reduced multiple testing burden**: Step 2 searches ~10-50x fewer peptides, dramatically improving statistical power

2. **Improved detection consistency**: Peptides missed in Step 1 due to stringent FDR can be recovered in Step 2 when competing against fewer decoys

3. **Better experiment-level FDR**: The constrained search space makes cross-run evidence aggregation more meaningful

4. **Computational efficiency**: Step 2 is faster because fewer candidates are evaluated per spectrum

**Comparison to DIA-NN/Spectronaut:**

| Aspect | DIA-NN | Spectronaut | Osprey |
|--------|--------|-------------|--------|
| Step 1 scope | Per-file | Per-file | Per-file |
| Step 2 trigger | Automatic | Automatic | Configurable |
| Decoy handling | Reuse | Regenerate | Regenerate |
| MBR-like feature | Yes | Yes | Via Step 2 |

### 6.5 Feature Computation Pipeline

```
Algorithm: Peptide Feature Extraction

Input: Peptide P, coefficient time series X, predicted spectrum S

1. Peak Detection:
   - Smooth X with Savitzky-Golay filter
   - Find peaks using prominence threshold
   - Select best peak (highest apex near predicted RT)

2. Chromatographic Features:
   - apex_coefficient = max(X) in peak region
   - peak_area = trapezoid integration of X
   - Fit EMG → emg_fit_quality, peak_width, peak_symmetry
   - rt_deviation = |apex_rt - predicted_rt|
   - n_scans = count(X > threshold) in peak region
   - coefficient_stability = var(X) near apex
   - peak_prominence = (apex - baseline) / baseline

3. Spectral Features:
   - Aggregate observed spectra across peak apex (±N scans)
   - Compute consensus spectrum O
   - hyperscore = hyperscore(S, O)
   - spectral_angle = arccos(dot(S, O) / (|S| * |O|))
   - dot_product = dot(S, O) / (|S| * |O|)
   - fragment_coverage = |matched| / |predicted|
   - ... (remaining spectral features)

4. Background Features:
   - For each offset in background_offsets:
     - Repeat steps 1-3 at predicted_rt ± offset
     - Record background feature values
   - Compute median background for each feature
   - Compute subtracted/contrast/zscore versions

5. Contextual Features:
   - n_competitors = count candidates in same RT/mz window
   - relative_coefficient = coefficient / sum(all coefficients)
   - ... (remaining contextual features)

6. Return FeatureSet
```

---

## 7. API Specification

### 7.1 Rust Public API

```rust
// Main entry point
pub fn analyze(config: OspreyConfig) -> Result<AnalysisResults, OspreyError>;

// Streaming API for large files
pub fn analyze_streaming(
    config: OspreyConfig,
    callback: impl FnMut(PeptideDetection) -> bool,
) -> Result<AnalysisSummary, OspreyError>;

// Individual components for advanced usage
pub mod components {
    // Library loading
    pub fn load_library(source: &LibrarySource) -> Result<Vec<LibraryEntry>, OspreyError>;
    pub fn load_diann_tsv(path: &Path) -> Result<Vec<LibraryEntry>, OspreyError>;
    pub fn load_blib(path: &Path) -> Result<Vec<LibraryEntry>, OspreyError>;
    pub fn load_elib(path: &Path) -> Result<Vec<LibraryEntry>, OspreyError>;
    
    // Library utilities
    pub fn library_info(source: &LibrarySource) -> Result<LibraryInfo, OspreyError>;
    pub fn filter_library_by_proteins(library: &[LibraryEntry], proteins: &[String]) -> Vec<LibraryEntry>;
    pub fn subset_library_by_rt_range(library: &[LibraryEntry], rt_min: f64, rt_max: f64) -> Vec<LibraryEntry>;
    
    // Decoy generation
    pub fn generate_decoys(library: &[LibraryEntry], method: DecoyMethod) -> Vec<LibraryEntry>;
    
    // Core analysis
    pub fn build_design_matrix(spectrum: &BinnedSpectrum, candidates: &[&LibraryEntry], config: &BinConfig) -> DesignMatrix;
    pub fn solve_ridge_regression(matrix: &DesignMatrix, observed: &[f64], lambda: f64) -> Vec<f64>;
    pub fn extract_features(entry: &LibraryEntry, coefficients: &[(f64, f64)], config: &OspreyConfig) -> FeatureSet;
    pub fn compute_fdr(targets: &[FeatureSet], decoys: &[FeatureSet]) -> Vec<f64>;
}

// Skyline integration
pub mod skyline {
    pub fn export_to_skyline(results: &AnalysisResults, path: &Path) -> Result<(), OspreyError>;
    pub fn load_library_from_document(sky_path: &Path) -> Result<Vec<LibraryEntry>, OspreyError>;
}

// Library format conversion utilities
pub mod convert {
    pub fn diann_to_blib(tsv_path: &Path, blib_path: &Path) -> Result<(), OspreyError>;
    pub fn blib_to_diann(blib_path: &Path, tsv_path: &Path) -> Result<(), OspreyError>;
}
```

### 7.2 Command-Line Interface

```
osprey 1.0.0
Peptide-centric DIA analysis with Skyline integration

USAGE:
    osprey [OPTIONS] --input <FILES>... --library <LIBRARY> --output <BLIB>

OPTIONS:
    -i, --input <FILES>...           Input mzML file(s)
    -l, --library <LIBRARY>          Spectral library file (.tsv for DIA-NN, .blib, or .elib)
    -o, --output <BLIB>              Output results file (.blib format for Skyline)
    
    --library-format <FORMAT>        Library format: diann, blib, elib [default: auto-detect]
    
    --resolution <MODE>              Resolution mode: unit, hram, auto [default: auto]
    --hram-tolerance <PPM>           HRAM tolerance in ppm [default: 20]
    
    --rt-tolerance <MIN>             RT window for candidates [default: 2.0]
    --precursor-tolerance <TH>       Precursor m/z tolerance [default: 0.5]
    
    --lambda <VALUE>                 Fixed regularization parameter
    --lambda-cv                      Select lambda via cross-validation

TWO-STEP SEARCH OPTIONS:
    --two-step                       Enable two-step search [default: enabled for multi-file]
    --no-two-step                    Disable two-step search
    --step1-fdr <VALUE>              Step 1 FDR threshold [default: 0.01]
    --step2-min-runs <N>             Minimum runs for peptide inclusion in Step 2 [default: 1]

FDR OPTIONS:
    --run-fdr <VALUE>                Run-level FDR threshold [default: 0.01]
    --experiment-fdr <VALUE>         Experiment-level FDR threshold [default: 0.01]
    --decoy-method <METHOD>          Decoy generation: reverse, shuffle, library [default: reverse]
    --decoys-in-library              Library already contains decoys (skip generation)

PERFORMANCE OPTIONS:
    --threads <N>                    Number of threads [default: all]
    --memory-limit <GB>              Memory limit in GB

OUTPUT OPTIONS:
    --export-coefficients            Include coefficient time series in output blib
    --export-features                Export all features to separate TSV file
    --report <TSV>                   Also write human-readable TSV report
    
    -v, --verbose                    Verbose output
    -h, --help                       Print help
    -V, --version                    Print version

OUTPUT:
    The primary output is a .blib file containing:
    - Detected peptides with peak boundaries (StartRT, EndRT, ApexRT)
    - Run-level q-values for each peptide in each file
    - Experiment-level q-values across all files
    - Consensus spectra for detected peptides
    
    This blib is designed for direct import into Skyline for quantification.

EXAMPLES:
    # Basic analysis with DIA-NN library (two-step enabled by default)
    osprey -i *.mzML -l predicted_library.tsv -o results.blib
    
    # Single file analysis (two-step disabled automatically)
    osprey -i sample.mzML -l library.blib -o results.blib
    
    # HRAM-constrained library for unit resolution data
    osprey -i *.mzML -l hram_filtered.blib -o results.blib --resolution unit
    
    # Strict FDR control
    osprey -i *.mzML -l library.tsv -o results.blib --run-fdr 0.005 --experiment-fdr 0.01
    
    # Disable two-step search
    osprey -i *.mzML -l library.tsv -o results.blib --no-two-step
    
    # Export additional files for debugging
    osprey -i *.mzML -l library.tsv -o results.blib --export-coefficients --report report.tsv
```

### 7.3 Python Bindings (via PyO3)

```python
import osprey

# Simple analysis with DIA-NN library
results = osprey.analyze(
    input_files=["sample.mzML"],
    library="predicted_library.tsv",  # DIA-NN format
    fdr=0.01,
)

# Using blib library (HRAM-constrained)
results = osprey.analyze(
    input_files=["sample.mzML"],
    library="hram_filtered.blib",
    library_format="blib",
    resolution_mode="unit",
    fdr=0.01,
)

# Access results
for detection in results.detections:
    print(f"{detection.sequence}: q={detection.q_value:.4f}, area={detection.area:.2e}")

# Advanced configuration
config = osprey.Config(
    resolution_mode="unit",
    rt_tolerance=2.0,
    lambda_setting="cross_validated",
    decoys_in_library=False,
    decoy_method="reverse",
)
results = osprey.analyze(
    input_files=["sample.mzML"],
    library="library.tsv",
    config=config,
)

# Export to pandas
df = results.to_dataframe()

# Get library statistics
info = osprey.library_info("library.tsv")
print(f"Library contains {info.n_targets} targets, {info.n_decoys} decoys")
print(f"RT type: {info.rt_type}")
```

---

## 8. Testing Strategy

### 8.1 Unit Tests

| Component | Test Focus | Coverage Target |
|-----------|------------|-----------------|
| DIA-NN TSV parsing | Column mapping, edge cases, large files | 100% |
| blib parsing | SQLite queries, modification handling | 100% |
| elib parsing | Chromatogram data, RT calibration | 100% |
| Decoy generation | Reversal, shuffling, fragment recalculation | 100% |
| Binning | Correct bin assignment, edge cases | 100% |
| Matrix construction | Sparse/dense equivalence, normalization | 95% |
| Ridge regression | Known solutions, convergence, non-negativity | 100% |
| EMG fitting | Synthetic peaks, edge cases | 95% |
| Feature extraction | Each feature independently | 100% |
| FDR calculation | Known target/decoy mixtures | 100% |
| mzML parsing | Valid/invalid inputs, large files | 90% |

### 8.2 Integration Tests

| Test | Description |
|------|-------------|
| End-to-end unit resolution | Full pipeline on unit resolution test data |
| End-to-end HRAM | Full pipeline on HRAM test data |
| Skyline round-trip | Export → Skyline import → verify |
| Batch processing | Multiple files in sequence |
| Memory limits | Processing under constrained memory |
| Determinism | Same input → identical output |

### 8.3 Benchmark Tests

| Benchmark | Metric | Target |
|-----------|--------|--------|
| Regression speed | Spectra/second | >10,000 |
| Feature extraction | Peptides/second | >50,000 |
| File parsing | MB/second | >100 |
| Memory efficiency | Peak RSS / file size | <3x |
| Parallel scaling | Speedup vs cores | >0.8 × cores |

### 8.4 Validation Tests

| Test | Ground Truth | Success Criterion |
|------|--------------|-------------------|
| Spike-in detection | Known peptide concentrations | Sensitivity >80% at 1% FDR |
| Quantitative accuracy | Known ratios | Pearson r >0.95 |
| FDR calibration | Entrapment peptides | Observed FDR ≤ reported FDR |
| HRAM comparison | DIA-NN/Spectronaut results | Overlap >80%, comparable sensitivity |

---

## 9. Deployment

### 9.1 Distribution Channels

| Channel | Format | Priority |
|---------|--------|----------|
| GitHub Releases | Pre-built binaries (Linux, macOS, Windows) | High |
| Conda | conda-forge package | High |
| Cargo | crates.io | Medium |
| Docker | Container image | Medium |
| Bioconda | Bioinformatics package | Medium |
| Skyline Tool Store | Skyline plugin | High |

### 9.2 Build Requirements

```toml
# Cargo.toml dependencies (key items)
[dependencies]
ndarray = { version = "0.15", features = ["blas"] }
ndarray-linalg = { version = "0.16", features = ["openblas-static"] }
sprs = "0.11"
rayon = "1.8"
quick-xml = "0.31"
pyo3 = { version = "0.20", features = ["extension-module"], optional = true }
clap = { version = "4.4", features = ["derive"] }
serde = { version = "1.0", features = ["derive"] }
thiserror = "1.0"
log = "0.4"
env_logger = "0.10"

[build-dependencies]
cbindgen = "0.26"  # For C header generation (Skyline FFI)

[features]
default = ["python"]
python = ["pyo3"]
```

### 9.3 CI/CD Pipeline

```yaml
# GitHub Actions workflow (simplified)
name: CI

on: [push, pull_request]

jobs:
  test:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - run: cargo test --all-features
      - run: cargo bench --no-run

  lint:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          components: clippy, rustfmt
      - run: cargo fmt --check
      - run: cargo clippy -- -D warnings

  release:
    needs: [test, lint]
    if: startsWith(github.ref, 'refs/tags/')
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        include:
          - os: ubuntu-latest
            target: x86_64-unknown-linux-gnu
          - os: macos-latest
            target: x86_64-apple-darwin
          - os: macos-latest
            target: aarch64-apple-darwin
          - os: windows-latest
            target: x86_64-pc-windows-msvc
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
        with:
          targets: ${{ matrix.target }}
      - run: cargo build --release --target ${{ matrix.target }}
      - uses: actions/upload-artifact@v4
        with:
          name: osprey-${{ matrix.target }}
          path: target/${{ matrix.target }}/release/osprey*
```

---

## 10. Documentation Plan

### 10.1 Documentation Types

| Type | Location | Format |
|------|----------|--------|
| API Reference | docs.rs/osprey | rustdoc |
| User Guide | osprey.maccosslab.org | mdBook |
| Tutorials | GitHub wiki | Markdown |
| Method Paper | Journal submission | PDF/LaTeX |

### 10.2 User Guide Outline

1. **Getting Started**
   - Installation
   - Quick start example
   - Input file requirements

2. **Spectral Libraries**
   - Supported formats (DIA-NN TSV, blib, elib)
   - Creating libraries with DIA-NN
   - Creating libraries with Skyline/EncyclopeDIA
   - HRAM-constrained libraries for improved FDR
   - Library size considerations and filtering strategies

3. **Core Concepts**
   - Peptide-centric analysis
   - Ridge regression deconvolution
   - Feature scoring and FDR

4. **Configuration**
   - Resolution modes
   - Candidate selection parameters
   - Regularization settings
   - FDR thresholds

5. **Skyline Integration**
   - Setting up the plugin
   - Workflow walkthrough
   - Exporting libraries from Skyline
   - Interpreting results in Skyline

6. **Advanced Usage**
   - Batch processing multiple files
   - Programmatic API (Rust and Python)
   - Custom feature extraction
   - Extending Osprey

7. **Best Practices**
   - Library selection for different applications
   - Unit resolution vs HRAM workflows
   - Optimizing FDR control

8. **Troubleshooting**
   - Common errors
   - Performance tuning
   - FAQ

### 10.3 Tutorial Datasets

| Dataset | Description | Size | Purpose |
|---------|-------------|------|---------|
| HeLa-unit | HeLa digest, unit resolution DIA | 500 MB | Basic tutorial |
| HeLa-HRAM | Same sample, HRAM DIA | 2 GB | Resolution comparison |
| Spike-in | Known concentrations | 1 GB | Quantitative validation |
| Phospho | Phosphopeptide-enriched | 1 GB | Modification handling |

---

## 11. Project Timeline

### Phase 1: Core Engine (Months 1-6)

| Month | Milestone |
|-------|-----------|
| 1 | Project setup, mzML parsing, data structures |
| 2 | Binning implementation, design matrix construction |
| 3 | Ridge regression solver, warm-start optimization |
| 4 | Carafe integration, spectral prediction interface |
| 5 | Basic feature extraction (chromatographic) |
| 6 | Validation on HRAM data, initial benchmarks |

**Deliverable**: Working prototype processing HRAM DIA data

### Phase 2: Full Feature Set (Months 7-12)

| Month | Milestone |
|-------|-----------|
| 7 | Unit resolution binning, sparse matrix support |
| 8 | Complete spectral feature extraction |
| 9 | Background correction implementation |
| 10 | Mokapot/Percolator integration |
| 11 | FDR calibration and validation |
| 12 | Performance optimization, benchmarking |

**Deliverable**: Feature-complete tool with validated FDR control

### Phase 3: Integration and Release (Months 13-18)

| Month | Milestone |
|-------|-----------|
| 13 | Skyline plugin development |
| 14 | Skyline integration testing |
| 15 | User documentation, tutorials |
| 16 | Beta release, community testing |
| 17 | Bug fixes, performance tuning |
| 18 | v1.0 release, manuscript submission |

**Deliverable**: Production release with Skyline integration

### Phase 4: Refinement (Months 19-24)

| Month | Milestone |
|-------|-----------|
| 19-20 | Community feedback incorporation |
| 21-22 | Advanced features (user requests) |
| 23-24 | Long-term maintenance planning, v1.1 |

**Deliverable**: Mature, community-tested tool

---

## 12. Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Carafe prediction accuracy insufficient | Low | High | Fall back to theoretical spectra; retrain Carafe |
| Performance targets not met | Medium | Medium | Profile early; consider GPU acceleration |
| FDR not well-calibrated at unit resolution | Medium | High | Extensive validation; conservative defaults |
| Skyline integration complexity | Medium | Medium | Early engagement with Skyline team |
| Limited community adoption | Low | Medium | Comprehensive documentation; active support |
| Dependency on external libraries | Low | Low | Pin versions; consider vendoring critical deps |

---

## 13. Success Metrics

### 13.1 Technical Metrics

| Metric | Target | Measurement |
|--------|--------|-------------|
| Detection sensitivity (HRAM) | Within 10% of DIA-NN | Benchmark datasets |
| Detection sensitivity (unit) | >70% of HRAM | Same sample comparison |
| Quantitative accuracy | Pearson r >0.95 | Spike-in experiments |
| FDR calibration | Observed ≤ reported | Entrapment validation |
| Processing speed | <10 min / 30 min run | Wall clock time |

### 13.2 Adoption Metrics

| Metric | Target (Year 1) | Target (Year 2) |
|--------|-----------------|-----------------|
| GitHub stars | 100 | 500 |
| Downloads (all channels) | 1,000 | 10,000 |
| Citations | 5 | 50 |
| Community contributions | 5 PRs | 20 PRs |
| Active users (Skyline) | 500 | 2,000 |

---

## 14. License and Governance

### 14.1 License

- **Core software**: Apache License 2.0
- **Documentation**: CC BY 4.0
- **Example data**: CC0 (public domain)

### 14.2 Governance

- **Maintainers**: MacCoss Lab, University of Washington
- **Contribution process**: GitHub pull requests with review
- **Decision making**: Maintainer consensus; community RFC for major changes
- **Code of conduct**: Contributor Covenant

### 14.3 Long-Term Sustainability

- Integration with Skyline ensures ongoing maintenance
- Apache 2.0 license allows commercial support options
- Modular architecture enables component reuse
- Community contributions reduce single-point-of-failure risk

---

## Appendix A: Glossary

| Term | Definition |
|------|------------|
| DIA | Data-Independent Acquisition |
| HRAM | High-Resolution Accurate Mass |
| XIC | Extracted Ion Chromatogram |
| EMG | Exponentially Modified Gaussian |
| FDR | False Discovery Rate |
| RT | Retention Time |
| PSM | Peptide-Spectrum Match |
| Ridge regression | L2-regularized least squares |
| LASSO | L1-regularized least squares (Least Absolute Shrinkage and Selection Operator) |
| Comet binning | Fragment m/z binning scheme from Comet search engine |
| mokapot | Python tool for semi-supervised PSM rescoring |
| Percolator | SVM-based PSM rescoring tool |

## Appendix B: References

1. Ting YS et al. (2015) Peptide-Centric Proteome Analysis. MCP. DOI: 10.1074/mcp.O114.047035
2. Eng JK et al. (2015) A deeper look into Comet. JASMS. DOI: 10.1007/s13361-015-1179-x
3. Pino LK et al. (2020) Acquiring and Analyzing DIA Experiments without Spectrum Libraries. MCP. DOI: 10.1074/mcp.P119.001913
4. Searle BC et al. (2018) Chromatogram libraries improve peptide detection. Nat Commun. DOI: 10.1038/s41467-018-07454-w
5. Demichev V et al. (2020) DIA-NN. Nat Methods. DOI: 10.1038/s41592-019-0638-x
6. Bekker-Jensen DB et al. (2020) Rapid and site-specific deep phosphoproteome profiling. Nat Commun. DOI: 10.1038/s41467-020-14609-1
7. The M et al. (2016) Fast and Accurate Protein FDR with Picked Protein FDR. J Proteome Res. DOI: 10.1021/acs.jproteome.6b00144

---

*Document version: 1.0*  
*Last updated: February 2026*
