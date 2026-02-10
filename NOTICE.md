# Third-Party Notices

Osprey includes code and algorithmic ideas from the following projects.

---

## Sage

**Project:** https://github.com/lazear/sage
**Copyright:** (c) 2022 Michael Lazear
**License:** MIT

The streaming mzML parser (`crates/osprey-io/src/mzml/streaming.rs`) is adapted
from Sage's mzML parser with modifications for Osprey (f64 precision, Osprey
Spectrum types, channel-based output).

### MIT License

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

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
