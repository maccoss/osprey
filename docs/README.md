# Osprey Algorithm Documentation

This folder contains detailed documentation of the algorithms used in Osprey, a peptide-centric DIA analysis tool. Osprey uses a coelution search approach (DIA-NN-style fragment XIC correlation) to detect peptides in DIA data.

## Current Status (Working Prototype)

Osprey v0.1.0 features:
- Parse mzML files with DIA data
- Load spectral libraries (DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib)
- Generate enzyme-aware decoys
- Auto-calibration with target-decoy FDR control (RT, MS1, MS2)
- Coelution search with 21 features extracted per precursor
- FDR control via native Percolator (default), external Mokapot, or simple TDC
- Two-level FDR (run + experiment) at precursor and peptide levels
- Protein-level FDR with parsimony and picked-protein approach
- Cross-run peak reconciliation with consensus RT alignment
- Disk-backed memory architecture for 1000+ file experiments
- Automatic cache invalidation via SHA-256 parameter hashing
- BiblioSpec .blib output for Skyline (library theoretical fragments)
- Safe file writes for NAS/CIFS mounts (copy_and_verify pattern)

## Algorithm Documentation

| Document | Description |
|----------|-------------|
| [01 - Decoy Generation](01-decoy-generation.md) | Enzyme-aware sequence reversal for FDR control |
| [02 - Calibration](02-calibration.md) | RT, MS1, and MS2 calibration with LDA scoring and LOESS fitting |
| [03 - Spectral Scoring](03-spectral-scoring.md) | Spectral scoring features for calibration (XCorr, E-value, isotope cosine) |
| [04 - XCorr Scoring](04-xcorr-scoring.md) | Comet-style XCorr implementation |
| [05 - Peak Detection](05-peak-detection.md) | CWT consensus peak detection and selection in fragment XIC time series |
| [06 - Multi-Charge Consensus](06-multi-charge-consensus.md) | Cross-charge-state peak boundary sharing |
| [07 - FDR Control](07-fdr-control.md) | Two-level FDR: native Percolator, Mokapot, and simple TDC |
| [08 - BiblioSpec Output Schema](08-blib-output-schema.md) | blib output format and Skyline integration |
| [09 - Determinism](09-determinism.md) | Deterministic analysis: patterns, invariants, and maintenance |
| [10 - Cross-Run Reconciliation](10-cross-run-reconciliation.md) | Consensus RT computation and peak alignment across replicates |
| [12 - Intermediate File Formats](12-intermediate-files.md) | Calibration JSON, spectra cache, Parquet scores, blib output |
| [13 - Boundary Overrides & Forced Integration](13-boundary-overrides.md) | How `run_search()` handles override boundaries for consensus and reconciliation |
| [14 - RT Alignment & Consensus Library RT](14-rt-alignment.md) | LOWESS calibration, inverse prediction, weighted median consensus RT, peak imputation |
| [15 - Protein FDR](07-fdr-control.md#protein-level-fdr) | Protein parsimony, picked-protein FDR, shared peptide modes, protein PEP |

---

## Pipeline Overview

Osprey's pipeline has five major phases.

```
INPUT
  mzML files (DIA data) + Spectral Library (predicted spectra)
  |
  v
PHASE 1: INITIALIZATION
  +- Load spectral library (DIA-NN TSV, elib, or blib)
  +- Deduplicate by (modified_sequence, charge)
  +- Generate decoys (enzyme-aware reversal)
  +- Build m/z index for candidate lookup
  |
  v
PHASE 2: CALIBRATION DISCOVERY (per file)
  +- Calculate library-to-measured RT mapping
  +- Set wide RT tolerance (20-50% of gradient)
  +- Co-elution search with CWT peak detection
  |     +- Fragment XIC correlation, spectral similarity (libcosine), XCorr
  |     +- Top-3 fragment matching (binary search) -> MS2 mass errors
  +- LDA scoring (4 features: correlation, libcosine, top6_matched, xcorr)
  +- Paired target-decoy competition on LDA score -> 1% FDR filtering
  +- Fit LOESS RT calibration from confident targets
  +- Calculate MS1/MS2 mass error statistics
  +- Set tight RT tolerance for main search (3x residual SD)
  |
  v
PHASE 3: COELUTION SEARCH (per file)
  +- Pre-filter: 3-of-4 consecutive scans must have ≥2 top-6 fragment matches
  |    (~30% speedup, minimal sensitivity loss; disable with --no-prefilter)
  +- Extract fragment XICs per precursor within RT window
  +- CWT consensus peak detection (Mexican Hat wavelet, median across transitions)
  +- Fallback: median polish profile or reference XIC peak detection
  +- Peak boundaries: zero-crossings ±2σ (valley guard)
  +- Feature extraction (21 PIN features used for scoring, ~47 computed total)
  +- Write full scored entries to per-file Parquet cache (ZSTD compressed)
  +- Convert to FdrEntry stubs (~80 bytes inline each) — drop features, fragments, CWT
  |
  v
PHASE 4: FIRST-PASS FDR + POST-FDR CONSENSUS
  +- Best-per-precursor subsampling: find best observation per base_id across files
  |    (maximizes peptide diversity for SVM training on multi-file experiments)
  +- Load PIN features on-demand from per-file Parquet caches
  +- FDR control (Percolator SVM or Mokapot) -> run-level q-values
  +- Multi-charge consensus: charge states of same peptide share peak boundaries
  |    (post-FDR — uses SVM scores to select best charge state as consensus leader)
  |    Re-scored via run_search() with boundary_overrides
  |
  v
PHASE 5: CROSS-RUN RECONCILIATION (multi-file only)
  +- Compact FDR stubs: drop non-passing entries (saves ~21 GB for 240-file experiments)
  |    (parquet_index field preserves original Parquet row reference after compaction)
  +- Collect peptides passing run-level FDR in each replicate (only FDR-passing detections)
  +- Compute consensus library RTs: weighted median of per-run detections
  |    (weight = coelution_sum; mapped back to library RT space via inverse calibration)
  |    (peak widths also use weighted median to resist noisy wide detections)
  +- Refit per-run LOESS calibration using consensus peptides
  +- Plan reconciliation for each entry:
  |    Keep | UseCwtPeak (alternate stored candidate) | ForcedIntegration
  |    (CWT candidates loaded selectively from Parquet — no full entry reload)
  +- Merge reconciliation targets with multi-charge consensus targets
  +- Re-score via run_search() with boundary_overrides per file
  |    (parallel window processing + shared XCorr preprocessing — same as first pass)
  |    (files processed sequentially to limit memory — each loads ~3 GB)
  +- Second-pass FDR -> final q-values
  |
  v
PHASE 6: PROTEIN FDR (optional, --protein-fdr)
  +- Build protein parsimony from library (bipartite graph, identical-set grouping,
  |    subset elimination, shared peptide assignment: All/Razor/Unique)
  +- Picked-protein FDR: target group vs decoy shadow (best peptide SVM score)
  |    (peptide scores collected BEFORE stub compaction to include decoy winners)
  +- Compute protein PEP from picked competition winners
  +- Propagate protein q-values to FdrEntry stubs
  +- Write protein CSV report (accessions, gene names, q-value, PEP, peptides)
  |
  v
OUTPUT
  +- Load lightweight BlibPlanEntry (~96 bytes) from 5-column Parquet projection
  +- Look up fragments, modifications, proteins from in-memory library (no full entry reload)
  +- BiblioSpec (.blib) for Skyline (library theoretical fragments)
  +- Protein CSV report (if --protein-fdr enabled)
  +- SVM/Mokapot model weights for feature analysis
```

---

## Feature Set (21 PIN Features)

Osprey computes ~47 features per precursor in the `CoelutionFeatureSet` struct, but only **21 are written to the PIN file** for semi-supervised FDR control via native Percolator (default) or external Mokapot. The remaining features were removed during feature weight optimization (they provided little discriminative value or were redundant). All intensity-based spectral similarity scores include ALL library fragments within the spectrum's mass range, using 0 intensity for unmatched peaks.

### Pairwise Coelution (3)

| Feature | Description |
|---------|-------------|
| `fragment_coelution_sum` | Sum of all pairwise Pearson correlations between fragment XICs |
| `fragment_coelution_max` | Maximum pairwise correlation (best pair) |
| `n_coeluting_fragments` | Number of fragments with positive mean correlation |

### Peak Shape (3)

| Feature | Description |
|---------|-------------|
| `peak_apex` | Peak apex intensity (background-subtracted) |
| `peak_area` | Integrated peak area within boundaries |
| `peak_sharpness` | Steepness of peak edges |

### Spectral at Apex (3)

| Feature | Description |
|---------|-------------|
| `xcorr` | Comet-style cross-correlation |
| `consecutive_ions` | Longest consecutive b or y ion series |
| `explained_intensity` | Fraction of observed intensity explained by matches |

### Mass Accuracy (2)

| Feature | Description |
|---------|-------------|
| `mass_accuracy_deviation_mean` | Mean signed mass error (ppm) across matched fragments |
| `abs_mass_accuracy_deviation_mean` | Mean absolute mass error (ppm) |

### RT Deviation (2)

| Feature | Description |
|---------|-------------|
| `rt_deviation` | Observed apex RT - predicted RT from LOESS calibration (minutes) |
| `abs_rt_deviation` | Absolute value of RT deviation |

### MS1 Features (2)

HRAM only; 0.0 for unit resolution or when MS1 data is unavailable:

| Feature | Description |
|---------|-------------|
| `ms1_precursor_coelution` | Pearson correlation between reference XIC and MS1 monoisotopic precursor XIC |
| `ms1_isotope_cosine` | Cosine similarity between observed MS1 isotope envelope and theoretical distribution |

### Tukey Median Polish (2)

| Feature | Description |
|---------|-------------|
| `median_polish_cosine` | Cosine similarity between data-derived fragment intensities (row effects) and library intensities, sqrt preprocessing |
| `median_polish_residual_ratio` | Fraction of total signal unexplained by the model: sum|obs-pred| / sum(obs) in linear space. Lower = cleaner co-elution. |

### SG-Weighted Multi-Scan (4)

Spectral scores computed on raw DIA spectra at apex ±2 scans, weighted by Savitzky-Golay smoothed reference XIC:

| Feature | Description |
|---------|-------------|
| `sg_weighted_xcorr` | XCorr on raw DIA spectra, Savitzky-Golay weighted |
| `sg_weighted_cosine` | Cosine similarity on raw DIA spectra, SG weighted |
| `median_polish_min_fragment_r2` | Minimum per-fragment R² with the shared elution profile (weakest link) |
| `median_polish_residual_correlation` | Mean pairwise correlation of median polish residuals across fragments |

### Removed Features (26, computed but not in PIN)

These features are still in the `CoelutionFeatureSet` struct but are **not written to the PIN file**:

| Category | Removed Features |
|----------|-----------------|
| Pairwise coelution | `fragment_corr_0..5`, `fragment_coelution_min`, `n_fragment_pairs` |
| Peak shape | `peak_width`, `peak_symmetry`, `signal_to_noise`, `n_scans` |
| Spectral at apex | `hyperscore`, `dot_product`, `dot_product_smz`, `dot_product_top4..6`, `dot_product_smz_top4..6`, `fragment_coverage`, `sequence_coverage`, `elution_weighted_cosine` |
| Mass accuracy | `mass_accuracy_std` |
| Peptide properties | `peptide_length`, `missed_cleavages` |
| Median polish | `median_polish_rsquared` |

---

## Tukey Median Polish

The Tukey median polish decomposes the fragment XIC matrix into an additive model:

```
ln(Observed[f,s]) = mu + alpha_f + beta_s + epsilon_fs

  f = fragment index (1..6, top 6 by library intensity)
  s = scan index (1..N)
  mu = overall effect (grand median)
  alpha_f = row effect = data-derived fragment relative intensity
  beta_s = column effect = shared elution profile shape
  epsilon_fs = residual = interference + noise
```

The algorithm iterates (max 20 iterations, convergence tolerance 1e-4):
1. **Row sweep**: subtract nanmedian of each row from residuals
2. **Column sweep**: subtract nanmedian of each column from residuals
3. Check convergence: max(|change|) < tolerance

### Uses

- **Peak boundaries**: Column effects give a robust elution profile by borrowing
  strength across all 6 transitions. The median suppresses interference from any
  single fragment. FWHM computed on this profile determines peak boundaries:
  sigma = FWHM / 2.355, boundaries at +/-1.96*sigma from apex.

- **Scoring features**: Row effects give interference-free fragment intensities
  that are scored against the library (median_polish_cosine). R^2 and residual
  ratio measure how well the additive model explains the data.

---

## FDR Control

Osprey uses **two-level FDR control** with three available methods (`--fdr-method`):

- **Percolator** (default): Native linear SVM with cross-validation. Both targets and decoys enter; Percolator does internal paired competition within each fold.
- **Mokapot**: External Python tool. Requires pre-competed PIN files (only competition winners). Uses ROC AUC to select the best feature for upstream competition.
- **Simple**: Direct target-decoy competition on the best single feature.

### Run-Level FDR
- Per-file target-decoy competition and q-value computation
- PIN features loaded on-demand from per-file Parquet caches (not held in memory)
- Q-values assigned at both precursor and peptide level
- Effective q-value: `max(precursor_qvalue, peptide_qvalue)` — must pass both

### Experiment-Level FDR
- Takes the **single best-scoring observation** per precursor (modified_sequence + charge) across all files
- Target-decoy competition on this deduplicated set at both precursor and peptide level
- Effective q-value: `max(precursor_qvalue, peptide_qvalue)`
- Single replicate skips experiment-level (uses run-level directly)

### Multi-File Observation Propagation
- After experiment-level FDR determines passing precursors, **all per-file target observations** for those precursors are included in the blib output
- Each file gets its own RT boundaries; best experiment_qvalue is propagated to all observations
- This enables Skyline to use per-file peak boundaries for quantification across replicates/GPF files

### Memory Architecture: FdrEntry Stubs + Parquet Caching

For multi-file experiments (100s-1000s of files), Osprey uses a disk-backed memory architecture to keep RAM usage manageable:

1. **Per-file Parquet caching**: After scoring each file, the full `CoelutionScoredEntry` data (features, fragments, CWT candidates) is written to a ZSTD-compressed Parquet cache file (`.scores.parquet`). All cache writers use a safe write pattern: write to local temp file, then `copy_and_verify` to the final destination. This prevents corrupt files on NAS/CIFS mounts if the process is interrupted.

2. **FdrEntry conversion**: Full entries are converted to lightweight `FdrEntry` stubs (~128 bytes each vs ~940 bytes for the full entry). Stubs retain only the fields needed for FDR control: `entry_id`, `parquet_index`, `is_decoy`, `charge`, `scan_number`, RTs, `coelution_sum`, `score`, 6 q-value fields (run/experiment x precursor/peptide/protein), `pep`, `modified_sequence`. The `modified_sequence` field uses `Arc<str>` for string interning, deduplicating identical peptide sequences across files (e.g., 240M entries -> ~3.5M unique strings).

3. **FDR stub compaction**: After first-pass FDR, stubs for precursors that didn't pass run-level FDR in any replicate are dropped. The `parquet_index` field preserves the original Parquet row reference for CWT/feature lookup after compaction. For 240-file experiments this reduces ~271M entries to ~106M, freeing ~21 GB of RAM for reconciliation.

4. **On-demand loading**: Downstream phases load selectively from the Parquet cache:
   - **FDR**: Loads PIN feature columns using `parquet_index` for lookup
   - **Reconciliation**: Loads CWT candidates only via `load_cwt_candidates_from_parquet()`
   - **Re-scoring**: Loads full entries per file sequentially (not parallel) to limit memory
   - **Output**: Loads a 5-column Parquet projection (entry_id, apex_rt, start_rt, end_rt, bounds_area) into `BlibPlanEntry` (~96 bytes each), with fragments and metadata looked up from the in-memory library

5. **LightFdr extraction**: Before blib output, the remaining FDR stubs (~128 bytes each) are converted to `LightFdr` structs (~48 bytes each) containing only q-values, score, pep, and grouping fields. The heavy stubs are then dropped, freeing memory for the blib loading phase.

6. **Cache invalidation**: Each Parquet file stores SHA-256 hashes in its footer metadata (`osprey.version`, `osprey.search_hash`, `osprey.library_hash`, `osprey.reconciled`, `osprey.reconciliation_hash`). Calibration JSON files also store a `search_hash`. On rerun, stale caches are automatically detected, deleted, and re-scored.

This architecture enables 240-file experiments (~275M total entries) on machines with 54 GB RAM. The blib output phase uses ~5 GB (BlibPlanEntry + library) instead of ~22 GB (full entries).

See [Intermediate File Formats](12-intermediate-files.md) for cache file formats, metadata hashing, and memory tier details.

### Protein-Level FDR

When `--protein-fdr <threshold>` is specified, Osprey runs native protein inference after precursor/peptide FDR:

1. **Protein parsimony**: Builds a bipartite peptide-protein graph from the library (using only detected peptides). Proteins with identical peptide sets are merged into groups. Groups whose peptides are a strict subset of another group are eliminated.

2. **Shared peptide handling** (`--shared-peptides`): Controls how peptides mapping to multiple protein groups are assigned:
   - **All** (default): Include shared peptides for all their groups
   - **Razor**: Assign shared peptides to the group with the most unique peptides
   - **Unique**: Use only proteotypic (unique) peptides

3. **Picked-protein FDR**: Each target protein group is paired with a decoy "shadow" group (same peptide set with `DECOY_` prefix). The best peptide SVM score determines the group score. Target and decoy groups compete; only the winner enters the ranked list for TDC q-value computation. Peptide scores are collected before FDR stub compaction to include decoy competition winners.

4. **Protein PEP**: Posterior error probabilities are estimated from the picked-protein competition winners using the same kernel density approach as peptide PEP.

5. **Output**: A CSV protein report (`{output}.proteins.csv`) with columns: protein_group, protein_accessions, gene_names, protein_qvalue, protein_pep, best_score, n_unique_peptides, n_shared_peptides, unique_peptides, shared_peptides.

See [FDR Control](07-fdr-control.md) for full algorithm details.

---

## Phase 2: Calibration Scoring (Quick Pass)

The calibration phase is a **fast scoring pass** designed to establish RT and mass calibration before the main search. Calibration always uses **unit resolution XCorr** (2001 bins, 1.0005 Da) regardless of data type, providing ~50x faster scoring for HRAM data.

### How Calibration Scoring Works

1. **Pre-filter** (top-3 fragment match): Binary search check -- at least 1 of the top-3 most intense library fragments must be present in the observed spectrum. Eliminates candidates with no spectral evidence before expensive XCorr preprocessing. O(3 x log n) per candidate.

2. **XCorr** (all passing spectra): Preprocess experimental spectrum with sqrt, windowing normalization, and flanking bin subtraction. Score against all spectra in RT window via BLAS sdot. Used to find the best-matching RT for each peptide.

3. **Top-3 fragment matching** (at best XCorr spectrum): Binary search the top-3 most intense library fragments in the best XCorr spectrum. Returns signed mass errors (ppm) for MS2 mass calibration.

4. **E-value** (significance): Calculated from the XCorr survival function (Comet-style) -- fit log-linear regression to survival counts. Computed and stored but not used for FDR -- LDA scoring is used instead.

5. **LDA scoring** (FDR): Linear Discriminant Analysis on 4 features (correlation, libcosine, top6_matched, xcorr). Produces a discriminant score used for paired target-decoy competition. See [Calibration](02-calibration.md).

5. **MS1 isotope envelope** (at best XCorr spectrum): Extract monoisotopic peak from nearest MS1 spectrum. Returns MS1 mass error (ppm) for precursor mass calibration.

See: [Spectral Scoring](03-spectral-scoring.md), [XCorr Scoring](04-xcorr-scoring.md)

### What Calibration Produces

- **RT calibration curve** (LOESS): Maps library RT -> expected measured RT
- **RT tolerance**: Residual SD x 3 (typically 0.5-2 min)
- **Mass calibration**: Mean/median PPM error for MS1 and MS2
- **Isolation scheme**: DIA window widths extracted from mzML

See: [Calibration](02-calibration.md)

### Multi-File Strategy

- **Each file**: Independent calibration discovery (each file may have different LC conditions)
- Calibration saved to JSON per file for reuse on reruns

---

## Phase 3: Main Search

The coelution search extracts fragment XICs directly from the DIA data and scores precursors based on how well their fragment ion chromatograms co-elute.

Before per-precursor scoring, a **pre-filter** eliminates candidates with no spectral evidence: at least 3 of 4 consecutive scans must contain ≥2 of the top-6 library fragment matches. This reduces scoring work by ~30% with negligible sensitivity loss. Disable with `--no-prefilter`.

For each precursor passing the pre-filter:

1. **Fragment XIC extraction**: Extract chromatograms for each library fragment ion at the expected m/z (ppm tolerance)
2. **CWT consensus peak detection**: Convolve each fragment XIC with a Mexican Hat wavelet, compute pointwise median across transitions, find peaks in consensus signal. Boundaries via zero-crossings extended to ±2σ with valley guard. See [Peak Detection](05-peak-detection.md).
3. **Spectral scoring**: At the apex scan, score the observed spectrum against the library (hyperscore, xcorr, dot products)
4. **Tukey median polish**: Decompose fragment XIC matrix for robust scoring features and blib peak boundaries

---

## Phase 4: Post-Search Scoring and First-Pass FDR

After the main search, Osprey computes scoring features per precursor from fragment XIC profiles.

### Step 1: Peak Detection

CWT consensus peak detection uses Mexican Hat wavelet convolution of each fragment XIC, pointwise median across transitions, and ±2σ boundary extension with valley guard. Falls back to SG-smoothed peak detection if CWT finds no peaks.

See: [Peak Detection](05-peak-detection.md)

### Step 2: Tukey Median Polish and Peak Boundaries

The top 6 library fragments' XICs are decomposed via Tukey median polish. The column
effects give a robust shared elution profile (resistant to single-fragment interference).
FWHM on this profile determines peak integration boundaries at +/-1.96*sigma from apex.

### Step 3: Score Computation

All 21 PIN features are computed within or at the peak boundaries:

- **Pairwise coelution**: Fragment XIC pairwise correlations
- **Peak shape**: Apex, area, width, symmetry, S/N, sharpness
- **Spectral at apex**: Hyperscore, xcorr, dot products, coverage, consecutive ions
- **Elution-weighted cosine**: LibCosine at each scan, weighted by reference XIC intensity^2
- **Mass accuracy**: Mean/std of mass errors across matched fragments
- **Tukey median polish features**: Cosine, R^2, residual ratio from the XIC decomposition
- **MS1 features**: Precursor co-elution and isotope cosine (HRAM only)
- **Peptide properties**: Peptide length, missed cleavages, RT deviation

### Step 4: First-Pass FDR and Multi-Charge Consensus

All 21 PIN features are combined into an optimal discriminant score via semi-supervised learning (linear SVM). The default native Percolator receives both targets and decoys and performs internal paired competition within each fold. For external Mokapot, targets and decoys are pre-competed using the best feature (selected by ROC AUC) and only winners are written to PIN files.

This produces **run-level q-values** for each file. Then **multi-charge consensus** is applied: all charge states of the same peptide are forced to share the same peak RT and integration boundaries, with the best SVM-scoring FDR-passing charge state as the consensus leader. Disagreeing charge states are re-scored at the consensus boundaries via `run_search()` with boundary overrides. See [Multi-Charge Consensus](06-multi-charge-consensus.md).

The run-level q-values are used as input to cross-run reconciliation.

See: [FDR Control](07-fdr-control.md)

### Step 5: Cross-Run Reconciliation (multi-file only)

For multi-file experiments, cross-run reconciliation aligns peak boundaries across replicates before the final FDR pass:

1. **Consensus RTs**: Peptides passing experiment-level FDR contribute (apex_rt → library_rt) pairs. The consensus library RT is computed as the weighted median across runs (weight = coelution_sum). See [RT Alignment & Consensus Library RT](14-rt-alignment.md).
2. **Calibration refit**: Each run's LOESS calibration is refit using consensus peptides, improving RT prediction accuracy for the second pass.
3. **Reconciliation plan**: For each entry, the expected measured RT is predicted from the refined calibration. Three outcomes: **Keep** (current peak contains expected RT), **UseCwtPeak** (alternate stored CWT candidate overlaps), or **ForcedIntegration** (no overlapping candidate — integrate at expected RT ± half of median peak width). See [Boundary Overrides & Forced Integration](13-boundary-overrides.md).
4. **Re-scoring via `run_search()`**: Reconciliation targets (merged with any remaining multi-charge consensus targets) are re-scored by calling `run_search()` with `boundary_overrides`. This reuses the same parallel window processing and per-window XCorr preprocessing as the first-pass search.
5. **Second-pass FDR**: Reconciled entries are re-scored and FDR is recomputed, producing final experiment-level q-values.

**Why this matters**: In low-signal runs or windows with interference, the peak selector may choose a completely wrong chromatographic peak for a peptide. Reconciliation uses the high-confidence runs to establish where the peptide truly elutes, then corrects the other runs — either by switching to an alternate stored CWT candidate or by imputing boundaries at the expected RT.

See: [Cross-Run Reconciliation](10-cross-run-reconciliation.md)

**Final two-level FDR with dual precursor+peptide control**:
1. **Run-level**: Per-file FDR at both precursor and peptide level
2. **Experiment-level**: Best observation per precursor across all files, FDR at both levels
3. **Effective q-value**: `max(precursor_qvalue, peptide_qvalue)` at each level
4. **Observation propagation**: All per-file observations for passing precursors included in output

See: [FDR Control](07-fdr-control.md)

---

## Key Distinctions

### Calibration vs Main Search

| Aspect | Calibration (Phase 2) | Main Search (Phase 3-4) |
|--------|----------------------|---------------------------|
| **Purpose** | Establish RT/mass calibration | Determine peptide detections |
| **Peak detection** | CWT consensus (with fallback) | CWT consensus (with fallback) |
| **Method** | Co-elution search + LDA | Coelution search |
| **Features** | 4 LDA features | 21 PIN features |
| **FDR method** | Target-decoy on LDA score | Semi-supervised SVM (Percolator/Mokapot) |
| **Scope** | Per file, wide tolerance | All files, tight calibrated tolerance |

---

## Why Peptide-Centric Analysis?

Unlike spectrum-centric approaches (e.g., Comet, MSFragger) that identify each spectrum independently, Osprey:

1. **Aggregates evidence across the chromatographic dimension** before scoring, using fragment XIC profiles rather than a single spectrum
2. **Scores at the best RT** for each peptide, where the co-elution signal is strongest
3. **Uses pairwise fragment co-elution** to distinguish true signals from interference

### Why Coelution Search?

DIA spectra are mixed spectra containing fragments from multiple precursors. The coelution approach:
- Extracts fragment XICs and measures how well they co-elute (pairwise correlations)
- Fragments from the same peptide should rise and fall together
- Robust to interference because co-elution metrics are computed pairwise across fragments
- Scores spectral similarity at the detected apex for additional discrimination

### Why RT Calibration?

Library retention times may be predicted by deep learning models (Prosit, DeepLC), measured on different LC systems, or normalized to iRT scale. RT calibration converts library RTs to expected measured RTs, enabling tight RT filtering that reduces the candidate set and improves search quality.

---

## Crate Architecture

```
osprey/
+-- osprey-core/              # Core types, config, errors, traits
|   +-- config.rs             # YAML configuration structures (incl. ReconciliationConfig)
|   +-- isotope.rs            # Isotope distribution calculations
|   +-- types.rs              # LibraryEntry, Spectrum, CoelutionFeatureSet, FdrEntry
|   +-- traits.rs             # MS1SpectrumLookup, etc.
|
+-- osprey-io/                # File I/O
|   +-- library/              # DIA-NN TSV, EncyclopeDIA elib, BiblioSpec blib
|   +-- mzml/                 # mzML parsing (via mzdata crate)
|   +-- output/               # blib and report output
|
+-- osprey-scoring/           # Spectral scoring and feature extraction
|   +-- lib.rs                # SpectralScorer, FeatureExtractor, DecoyGenerator,
|   |                         # Tukey median polish, fragment co-elution
|   +-- batch.rs              # BLAS-accelerated batch scoring (calibration)
|   +-- calibration_ml.rs     # LDA calibration scoring (classical LDA, 3-fold CV)
|   +-- pipeline/             # Streaming pipeline components
|
+-- osprey-chromatography/    # Chromatographic analysis
|   +-- lib.rs                # PeakDetector, detect_xic_peak (legacy fallback)
|   +-- cwt.rs                # CWT consensus peak detection (Mexican Hat wavelet)
|   +-- calibration/          # LOESS RT calibration, mass calibration, I/O
|
+-- osprey-ml/                # Machine learning primitives
|   +-- svm.rs                # Linear SVM (dual coordinate descent, Percolator training)
|   +-- linear_discriminant.rs# Classical LDA (scatter matrices + eigenvector, for calibration)
|   +-- pep.rs                # PEP estimation via KDE + isotonic regression
|   +-- kde.rs                # Kernel density estimation
|   +-- qvalue.rs             # Q-value computation
|   +-- matrix.rs             # Linear algebra operations
|
+-- osprey-fdr/               # FDR control
|   +-- lib.rs                # Target-decoy competition, q-values
|   +-- percolator.rs         # Native Percolator (SVM, fold assignment, score calibration)
|   +-- mokapot.rs            # PIN file output, Mokapot integration
|
+-- osprey/                   # Main library + CLI binary
    +-- main.rs               # CLI entry point
    +-- pipeline.rs           # Analysis pipeline orchestration
    +-- reconciliation.rs     # Cross-run consensus RT and peak alignment
```

## Utility Scripts

```
scripts/
+-- evaluate_calibration.py      # Generate HTML report from calibration JSON
+-- inspect_mokapot_weights.py   # Display Mokapot feature weights/importance
```

### evaluate_calibration.py

Generates an interactive HTML report visualizing calibration quality:
- RT calibration curve with residuals
- MS1/MS2 mass error histograms
- Candidate density heatmap (RT x m/z)

```bash
python scripts/evaluate_calibration.py calibration.json --output report.html
```

### inspect_mokapot_weights.py

Displays feature weights from trained Mokapot models to identify important features:

```bash
python scripts/inspect_mokapot_weights.py mokapot.model.pkl
```

---

## Configuration

Osprey supports YAML configuration files:

```bash
# Generate template
osprey --generate-config config.yaml

# Run with config
osprey --config config.yaml
```

Example configuration:

```yaml
# Input/Output
input_files:
  - sample1.mzML
  - sample2.mzML
spectral_library: library.tsv
output_blib: results.blib

# Resolution mode
resolution_mode:
  HRAM:
    tolerance_ppm: 10.0

# RT Calibration
rt_calibration:
  enabled: true
  loess_bandwidth: 0.3
  min_calibration_points: 50
  rt_tolerance_factor: 3.0

# FDR Control
run_fdr: 0.01

# Decoy Generation
decoy_method: Reverse

# Search pre-filter (speeds up ~30%, disable only if investigating sensitivity)
prefilter_enabled: true

# Cross-run peak reconciliation (multi-file only)
reconciliation:
  enabled: true
  consensus_fdr: 0.01   # FDR threshold for selecting consensus peptides
```

## References

- pyXcorrDIA: https://github.com/maccoss/pyXcorrDIA
- Comet: https://comet-ms.sourceforge.io/
- Mokapot: https://github.com/wfondrie/mokapot
- LOESS: Cleveland, W.S. (1979). "Robust locally weighted regression and smoothing scatterplots"
