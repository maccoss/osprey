# Third-Party Notices

Osprey includes code, libraries, and algorithmic ideas from the following projects.

---

## mzdata

**Project:** https://github.com/mobiusklein/mzdata
**Author:** Joshua Klein
**License:** Apache-2.0

Osprey uses the mzdata crate for parsing mzML files and reading mass spectrometry
data formats.

---

## DIA-NN

**Project:** https://github.com/vdemichev/DiaNN
**Authors:** Vadim Demichev et al.
**License:** BSD 2-Clause / Creative Commons Attribution 4.0 International (CC-BY 4.0)
**Citation:** Demichev, V., Messner, C.B., Vernardis, S.I. et al. DIA-NN: neural
networks and interference correction enable deep proteome coverage in high
throughput. *Nature Methods* 17, 41-44 (2020). https://doi.org/10.1038/s41592-019-0638-x

Osprey's fragment co-elution scoring and related algorithms are inspired by
DIA-NN's approach. Specifically:

- **Fragment co-elution correlation** (`compute_fragment_coelution`,
  `run_coelution_calibration_scoring`): The idea that co-eluting fragment XICs
  from the same peptide correlate temporally while interference does not,
  inspired by DIA-NN's pTimeCorr score.
- **Fragment XIC extraction** (`extract_fragment_xics`): Binary search-based
  extraction of fragment ion chromatograms from DIA spectra.
- **FWHM computation** (`compute_fwhm_interpolated`): Linear interpolation
  approach for chromatographic peak width estimation.
- **Consensus XIC for peak boundaries** (`compute_fragment_fwhm`): Building
  consensus XICs from co-eluting fragments for peak boundary determination.

No source code was directly copied from DIA-NN. These are independent Rust
implementations inspired by the published algorithmic concepts.

---

## Sage

**Project:** https://github.com/lazear/sage
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
