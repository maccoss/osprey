# Third-Party Notices

Osprey includes code, libraries, and algorithmic ideas from the following projects.

---

## Percolator

**Algorithm:** Semi-supervised learning for peptide detection
**Citation:** Käll L, Canterbury JD, Weston J, Noble WS, MacCoss MJ. Semi-supervised
learning for peptide identification from shotgun proteomics datasets. *Nature Methods*
4(11), 923-925 (2007). <https://doi.org/10.1038/nmeth1113>

Osprey implements a native Percolator-style FDR control pipeline inspired by the
original Percolator algorithm:

- **Iterative SVM training**: Train on high-confidence targets (passing FDR threshold)
  vs all decoys, rescore, repeat (`crates/osprey-fdr/src/percolator.rs`)
- **3-fold cross-validation**: With peptide-grouped fold assignment to prevent
  target-decoy pair leakage between folds
- **Conservative q-values**: `(n_decoy + 1) / n_target` formula at both precursor
  and peptide levels

No source code was copied from the C++ Percolator implementation. This is an
independent Rust implementation of the published algorithm.

---

## Mokapot

**Project:** <https://github.com/wfondrie/mokapot>
**Authors:** William E. Fondrie, William S. Noble
**License:** Apache-2.0
**Citation:** Fondrie WE, Noble WS. mokapot: Fast and Flexible Semisupervised
Learning for Peptide Detection. *Journal of Proteome Research* 20(4), 1966-1971
(2021). <https://doi.org/10.1021/acs.jproteome.0c01010>

Osprey supports mokapot as an external FDR engine via PIN file output and CLI
invocation (`crates/osprey-fdr/src/mokapot.rs`). The native Percolator
implementation was also informed by mokapot's refinements to the original
algorithm, including fold assignment strategy and training set selection.

---

## Score Calibration Between Cross-Validation Folds

**Citation:** Granholm V, Noble WS, Käll L. A cross-validation scheme for machine
learning algorithms in shotgun proteomics. *BMC Bioinformatics* 13(Suppl 16), S3
(2012). <https://doi.org/10.1186/1471-2105-13-S16-S3>

Osprey uses the Granholm et al. score calibration method to normalize SVM scores
across cross-validation folds (`crates/osprey-fdr/src/percolator.rs`). Each fold's
scores are linearly transformed so that the FDR threshold maps to 0 and the median
decoy score maps to -1, making scores comparable across folds.

---

## Posterior Error Probability (PEP) Estimation

**Citation:** Käll L, Storey JD, MacCoss MJ, Noble WS. Posterior error probabilities
and false discovery rates: two sides of the same coin. *Journal of Proteome Research*
7(1), 40-44 (2008). <https://doi.org/10.1021/pr700739d>

PEP estimation uses kernel density estimation (KDE) to model target and decoy score
distributions, Bayes' rule to compute P(incorrect | score), and isotonic regression
(pool adjacent violators) to enforce monotonicity (`crates/osprey-ml/src/pep.rs`).

---

## Linear SVM (Dual Coordinate Descent)

**Citation:** Hsieh CJ, Chang KW, Lin CJ, Keerthi SS, Sundararajan S. A Dual
Coordinate Descent Method for Large-scale Linear SVM. *Proceedings of the 25th
International Conference on Machine Learning (ICML)* (2008).
<https://doi.org/10.1145/1390156.1390208>

Osprey's linear SVM solver uses the dual coordinate descent algorithm for
L2-regularized L2-loss support vector classification (`crates/osprey-ml/src/svm.rs`).
This is the same algorithm used by LIBLINEAR and scikit-learn's LinearSVC. The
implementation includes grid search for the regularization parameter C and feature
standardization.

No source code was copied from LIBLINEAR. This is an independent Rust implementation
of the published algorithm.

---

## Comet

**Project:** <https://comet-ms.sourceforge.io/>
**Authors:** Jimmy Eng et al.
**Citation:** Eng JK, Jahan TA, Hoopmann MR. Comet: an open-source MS/MS sequence
database search tool. *Proteomics* 13(1), 22-24 (2013).
<https://doi.org/10.1002/pmic.201200439>

Osprey's XCorr implementation follows Comet's approach
(`crates/osprey-scoring/src/lib.rs`):

- **Binning**: 1.0005079 Da bins with 0.4 offset (unit resolution); 0.02 Da (HRAM)
- **Windowing normalization**: Divide spectrum into 10 windows, normalize peak
  intensities to 50.0 within each window
- **Flanking bin subtraction**: Sliding window of 75 bins subtracted from each bin
  to create the preprocessed spectrum
- **E-value**: Calculated from XCorr survival function histogram (log-linear fit)

No source code was copied from Comet. These are independent Rust implementations
of the published spectral scoring concepts.

---

## X!Tandem

**Project:** <https://www.thegpm.org/tandem/>
**Authors:** Ronald C. Beavis, Global Proteome Machine Organization
**Citation:** Craig R, Beavis RC. TANDEM: matching proteins with tandem mass spectra.
*Bioinformatics* 20(9), 1466-1467 (2004).
<https://doi.org/10.1093/bioinformatics/bth092>

Osprey implements the X!Tandem hyperscore for spectral matching
(`crates/osprey-scoring/src/lib.rs`):

```text
hyperscore = log(n_b!) + log(n_y!) + sum(log(matched_intensity + 1))
```

where n_b and n_y are the counts of matched b and y ions. This score rewards both
the number and intensity of matched fragment ions.

---

## DIA-NN

**Project:** <https://github.com/vdemichev/DiaNN>
**Authors:** Vadim Demichev et al.
**License:** BSD 2-Clause / Creative Commons Attribution 4.0 International (CC-BY 4.0)
**Citation:** Demichev V, Messner CB, Vernardis SI et al. DIA-NN: neural networks and
interference correction enable deep proteome coverage in high throughput. *Nature
Methods* 17, 41-44 (2020). <https://doi.org/10.1038/s41592-019-0638-x>

Osprey's coelution search mode and related algorithms are inspired by DIA-NN's
approach. Specifically:

- **Fragment co-elution correlation** (`compute_fragment_coelution`,
  `run_coelution_calibration_scoring`): The idea that co-eluting fragment XICs
  from the same peptide correlate temporally while interference does not,
  inspired by DIA-NN's pTimeCorr score.
- **Fragment XIC extraction** (`extract_fragment_xics`): Binary search-based
  extraction of fragment ion chromatograms from DIA spectra.
- **FWHM computation** (`compute_fwhm_interpolated`): Linear interpolation
  approach for chromatographic peak width estimation.
- **Valley-based boundary detection** (`walk_boundary_left/right`): Peak boundary
  detection using intensity threshold with valley detection for adjacent peaks.

No source code was directly copied from DIA-NN. These are independent Rust
implementations inspired by the published algorithmic concepts.

---

## Tukey Median Polish

**Reference:** Tukey JW. *Exploratory Data Analysis*. Addison-Wesley (1977).

**Implementation reference:** PRISM `rollup.py` (lines 640-731) for the application
of Tukey median polish to fragment XIC matrices.

Osprey applies Tukey's iterative median decomposition to fragment XIC matrices in log
space (`crates/osprey-scoring/src/lib.rs`):

```text
ln(Observed[f,s]) = mu + alpha_f + beta_s + epsilon_fs
```

Row effects give data-derived fragment intensities; column effects give a robust
shared elution profile. The median operation suppresses interference on individual
transitions, providing better peak boundaries than any single fragment XIC.

---

## LOESS

**Citation:** Cleveland WS. Robust locally weighted regression and smoothing
scatterplots. *Journal of the American Statistical Association* 74(368), 829-836
(1979). <https://doi.org/10.1080/01621459.1979.10481038>

Osprey uses LOESS (Locally Estimated Scatterplot Smoothing) for RT calibration,
fitting local linear regressions with tricube weighting and bisquare robustness
iterations (`crates/osprey-chromatography/src/calibration/rt.rs`).

---

## Target-Decoy Approach

**Citation:** Elias JE, Gygi SP. Target-decoy search strategy for increased
confidence in large-scale protein identifications by mass spectrometry. *Nature
Methods* 4(3), 207-214 (2007). <https://doi.org/10.1038/nmeth1019>

**Conservative FDR:** Levitsky LI, Ivanov MV, Lobas AA, Gorshkov MV. Unbiased
False Discovery Rate Estimation for Shotgun Proteomics Based on the Target-Decoy
Approach. *Journal of Proteome Research* 16(2), 689-694 (2017).
<https://doi.org/10.1021/acs.jproteome.6b00144>

Osprey uses enzyme-aware sequence reversal for decoy generation and the conservative
`(n_decoy + 1) / n_target` formula for FDR estimation
(`crates/osprey-scoring/src/lib.rs`, `crates/osprey-fdr/src/lib.rs`).

---

## Sage

**Project:** <https://github.com/lazear/sage>
**Author:** Michael Lazear
**License:** MIT
**Copyright:** Copyright (c) 2022 Michael Lazear

Osprey incorporates machine learning modules from Sage for linear discriminant
analysis, kernel density estimation, and q-value calculation. These modules
enable data-driven scoring and FDR control for peptide identification.

The following modules were adapted from Sage:

- `osprey-ml/src/linear_discriminant.rs` - Linear Discriminant Analysis for target-decoy separation
- `osprey-ml/src/matrix.rs` - Matrix operations for scatter matrix computation
- `osprey-ml/src/gauss.rs` - Gauss-Jordan elimination for solving linear systems
- `osprey-ml/src/kde.rs` - Kernel Density Estimation for posterior error probabilities
- `osprey-ml/src/qvalue.rs` - Q-value calculation via target-decoy competition

The MIT License permits use, modification, and distribution with proper attribution.
Full license text is included in the source files.

---

## mzdata

**Project:** <https://github.com/mobiusklein/mzdata>
**Author:** Joshua Klein
**License:** Apache-2.0

Osprey uses the mzdata crate for parsing mzML files and reading mass spectrometry
data formats.

---

## Rust Dependencies

Osprey uses the following notable open-source Rust crates:

| Crate | License | Purpose |
|-------|---------|---------|
| `arrow` / `parquet` | Apache-2.0 | Apache Arrow columnar format and Parquet file writing |
| `ndarray` | BSD-2-Clause / Apache-2.0 | Matrix operations |
| `rayon` | Apache-2.0 / MIT | Parallel processing |
| `rusqlite` | MIT | SQLite for BiblioSpec blib output |
| `clap` | Apache-2.0 / MIT | CLI argument parsing |
| `serde` / `serde_yaml` | Apache-2.0 / MIT | YAML configuration |
| `chrono` | Apache-2.0 / MIT | Timestamps |
| `blas-src` / `openblas-src` | BSD-3-Clause | BLAS linear algebra |
