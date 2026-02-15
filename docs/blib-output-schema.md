# BiblioSpec Output Schema

Osprey outputs results in the **BiblioSpec (.blib)** format, a SQLite database schema used by Skyline for spectral library management. This document describes the tables written by Osprey, including standard BiblioSpec tables and Osprey-specific extension tables.

## Overview

The blib output contains:
- **Standard BiblioSpec tables**: Compatible with Skyline and other tools
- **Osprey extension tables**: Additional scoring and peak boundary information

```
┌─────────────────────────────────────────────────────────────────┐
│                    BiblioSpec Standard Tables                    │
├─────────────────────────────────────────────────────────────────┤
│  LibInfo              │  Library metadata (version, timestamp)  │
│  ScoreTypes           │  Score type definitions                 │
│  SpectrumSourceFiles  │  Source mzML file paths                 │
│  RefSpectra           │  Peptide spectrum entries               │
│  RefSpectraPeaks      │  Fragment m/z and intensities (blobs)   │
│  RefSpectraPeakAnnotations │  Peak annotations (empty)          │
│  Modifications        │  Modification positions and masses      │
│  Proteins             │  Protein accessions                     │
│  RefSpectraProteins   │  Spectrum-to-protein mappings           │
│  IonMobilityTypes     │  Ion mobility type definitions          │
└─────────────────────────────────────────────────────────────────┘
┌─────────────────────────────────────────────────────────────────┐
│                    Osprey Extension Tables                       │
├─────────────────────────────────────────────────────────────────┤
│  OspreyPeakBoundaries │  Peak boundaries per run (StartRT, etc) │
│  OspreyRunScores      │  Run-level q-values and scores          │
│  OspreyExperimentScores │  Experiment-level q-values            │
│  OspreyCoefficients   │  Coefficient time series (optional)     │
│  OspreyMetadata       │  Analysis metadata (version, settings)  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Standard BiblioSpec Tables

### LibInfo

Library metadata.

| Column | Type | Description |
|--------|------|-------------|
| libLSID | TEXT | Unique identifier (LSID format) |
| createTime | TEXT | Creation timestamp |
| numSpecs | INTEGER | Number of spectra (updated at finalize) |
| majorVersion | INTEGER | BiblioSpec major version (1) |
| minorVersion | INTEGER | BiblioSpec minor version (10) |

### ScoreTypes

Score type definitions. Osprey uses PERCOLATOR QVALUE (ID 14) which Skyline recognizes as a q-value.

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Score type ID (14 = PERCOLATOR QVALUE) |
| scoreType | TEXT | Score type name |

### SpectrumSourceFiles

Source mzML files. **Contains relative paths** from the blib file location for cross-platform compatibility.

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | File ID (referenced by RefSpectra.fileID) |
| fileName | TEXT | Relative path to mzML file (from blib location) |
| idFileName | TEXT | Same as fileName |
| cutoffScore | REAL | Score cutoff (0.0) |

**Path handling**: Osprey computes relative paths from the blib output directory to the mzML files (e.g., `../data/sample.mzML` or just `sample.mzML` if in the same directory). This ensures compatibility when running Osprey in WSL2/Linux and opening results in Skyline on Windows. If the mzML files are in the same directory as the blib, Skyline will find them automatically.

### RefSpectra

Main spectrum table. Each row represents a detected peptide precursor.

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Spectrum ID (primary key) |
| peptideSeq | TEXT | Unmodified amino acid sequence |
| precursorMZ | REAL | Precursor m/z |
| precursorCharge | INTEGER | Charge state |
| peptideModSeq | TEXT | Modified sequence with mass annotations |
| prevAA | TEXT | Previous amino acid ('-' for N-term) |
| nextAA | TEXT | Next amino acid ('-' for C-term) |
| copies | INTEGER | Copy count (1) |
| numPeaks | INTEGER | Number of fragment peaks |
| ionMobility | REAL | Ion mobility value (0.0) |
| collisionalCrossSectionSqA | REAL | CCS (0.0) |
| ionMobilityHighEnergyOffset | REAL | High energy offset (0.0) |
| ionMobilityType | INTEGER | Ion mobility type (0 = none) |
| retentionTime | REAL | Apex retention time (minutes) |
| startTime | REAL | Peak start time (minutes) |
| endTime | REAL | Peak end time (minutes) |
| totalIonCurrent | REAL | Total ion current (0.0) |
| moleculeName | TEXT | Molecule name (NULL) |
| chemicalFormula | TEXT | Chemical formula (NULL) |
| precursorAdduct | TEXT | Precursor adduct (NULL) |
| inchiKey | TEXT | InChI key (NULL) |
| otherKeys | TEXT | Other identifiers (NULL) |
| fileID | INTEGER | Source file ID (FK to SpectrumSourceFiles) |
| SpecIDinFile | TEXT | Spectrum ID in file (NULL) |
| score | REAL | Score (1 - q_value, higher is better) |
| scoreType | INTEGER | Score type ID (14 = PERCOLATOR QVALUE) |

**Note on sequences**: Osprey automatically:
- Strips flanking characters (underscores, periods, dashes) from sequences. Input formats like `_PEPTIDE_` or `K.PEPTIDE.R` are cleaned to `PEPTIDE`.
- Converts UniMod notation to mass notation for Skyline compatibility. `PEPTC[UniMod:4]IDE` becomes `PEPTC[+57.0215]IDE`.

### RefSpectraPeaks

Fragment peak data stored as binary blobs.

| Column | Type | Description |
|--------|------|-------------|
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| peakMZ | BLOB | m/z values (little-endian f64 array, zlib-compressed) |
| peakIntensity | BLOB | Intensity values (little-endian f32 array, zlib-compressed) |

**Fragment data source**: The m/z and intensity values stored here are the **library theoretical fragments** (b/y ions from the spectral library), not observed peaks from the DIA spectrum. This applies to both regression and coelution search modes. Skyline uses these theoretical fragment m/z values to build extracted-ion chromatograms for the correct transitions.

**Compression**: Blobs are zlib-compressed. If compression does not reduce size (e.g., for very small arrays), the raw bytes are stored instead. To read, attempt zlib decompression; if it fails or the decompressed size doesn't match `numPeaks * sizeof(type)`, treat the blob as raw bytes.

### RefSpectraPeakAnnotations

Peak annotations (required by Skyline but left empty by Osprey).

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Annotation ID |
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| peakIndex | INTEGER | Peak index |
| name | TEXT | Annotation name |
| formula | TEXT | Chemical formula |
| inchiKey | TEXT | InChI key |
| otherKeys | TEXT | Other identifiers |
| charge | INTEGER | Charge state |
| adduct | TEXT | Adduct |
| comment | TEXT | Comment |
| mzTheoretical | REAL | Theoretical m/z |
| mzObserved | REAL | Observed m/z |

### Modifications

Modification positions and masses for each spectrum.

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Modification ID |
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| position | INTEGER | 0-indexed position in sequence |
| mass | REAL | Modification mass delta (Da) |

### Proteins

Protein accessions (de-duplicated).

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Protein ID |
| accession | TEXT | Protein accession/identifier |

### RefSpectraProteins

Spectrum-to-protein mappings (many-to-many).

| Column | Type | Description |
|--------|------|-------------|
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| ProteinID | INTEGER | FK to Proteins.id |

### IonMobilityTypes

Ion mobility type definitions (standard BiblioSpec).

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Type ID |
| ionMobilityType | TEXT | Type name |

Types: 0=none, 1=driftTime(msec), 2=inverseK0(Vsec/cm^2), 3=compensation(V)

---

## Osprey Extension Tables

These tables provide additional information not part of the standard BiblioSpec schema.

### OspreyPeakBoundaries

Peak boundaries per detection. Useful for quantification and visualization.

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Boundary ID |
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| FileName | TEXT | Source file name |
| StartRT | REAL | Peak start time (minutes) |
| EndRT | REAL | Peak end time (minutes) |
| ApexRT | REAL | Peak apex time (minutes) |
| ApexIntensity | REAL | Apex coefficient value |
| IntegratedArea | REAL | Integrated peak area |

### OspreyRunScores

Run-level (per-file) scores and q-values.

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Score ID |
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| FileName | TEXT | Source file name |
| RunQValue | REAL | Run-level q-value (FDR) |
| DiscriminantScore | REAL | Mokapot discriminant score |
| PosteriorErrorProb | REAL | Posterior error probability (PEP) |

### OspreyExperimentScores

Experiment-level scores (aggregated across replicates).

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Score ID |
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| ExperimentQValue | REAL | Experiment-level q-value |
| NRunsDetected | INTEGER | Number of runs where detected |
| NRunsSearched | INTEGER | Total number of runs searched |

### OspreyCoefficients

Coefficient time series (optional, enabled with `--export-coefficients`).

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Coefficient ID |
| RefSpectraID | INTEGER | FK to RefSpectra.id |
| FileName | TEXT | Source file name |
| ScanNumber | INTEGER | Scan number |
| RT | REAL | Retention time (minutes) |
| Coefficient | REAL | Regression coefficient value |

### OspreyMetadata

Analysis metadata as key-value pairs.

| Column | Type | Description |
|--------|------|-------------|
| Key | TEXT | Metadata key (primary key) |
| Value | TEXT | Metadata value |

Common keys:
- `osprey_version`: Osprey version
- `rt_calibration_enabled`: Whether RT calibration was used
- `run_fdr`: Run-level FDR threshold
- `fdr_method`: FDR method used

---

## Indices

The following indices are created for query performance:

```sql
CREATE INDEX idx_refspectra_peptide ON RefSpectra(peptideSeq);
CREATE INDEX idx_refspectra_modseq ON RefSpectra(peptideModSeq);
CREATE INDEX idx_refspectra_mz ON RefSpectra(precursorMZ);
CREATE INDEX idx_peaks_refid ON RefSpectraPeaks(RefSpectraID);
CREATE INDEX idx_mods_refid ON Modifications(RefSpectraID);
CREATE INDEX idx_boundaries_refid ON OspreyPeakBoundaries(RefSpectraID);
CREATE INDEX idx_runscores_refid ON OspreyRunScores(RefSpectraID);
```

---

## Querying the blib

### Example: List all detected peptides with q-values

```sql
SELECT
    r.peptideSeq,
    r.peptideModSeq,
    r.precursorCharge,
    r.precursorMZ,
    r.retentionTime,
    e.ExperimentQValue
FROM RefSpectra r
JOIN OspreyExperimentScores e ON r.id = e.RefSpectraID
WHERE e.ExperimentQValue <= 0.01
ORDER BY e.ExperimentQValue;
```

### Example: Get peak boundaries for a peptide

```sql
SELECT
    r.peptideModSeq,
    b.FileName,
    b.StartRT,
    b.ApexRT,
    b.EndRT,
    b.IntegratedArea
FROM RefSpectra r
JOIN OspreyPeakBoundaries b ON r.id = b.RefSpectraID
WHERE r.peptideSeq = 'PEPTIDE';
```

### Example: Get run-level scores

```sql
SELECT
    r.peptideModSeq,
    s.FileName,
    s.RunQValue,
    s.DiscriminantScore,
    s.PosteriorErrorProb
FROM RefSpectra r
JOIN OspreyRunScores s ON r.id = s.RefSpectraID
WHERE s.RunQValue <= 0.01;
```

### Example: Get protein mappings

```sql
SELECT
    r.peptideSeq,
    p.accession as ProteinAccession
FROM RefSpectra r
JOIN RefSpectraProteins rp ON r.id = rp.RefSpectraID
JOIN Proteins p ON rp.ProteinID = p.id;
```

### Example: Get modifications for a spectrum

```sql
SELECT
    r.peptideSeq,
    r.peptideModSeq,
    m.position,
    m.mass
FROM RefSpectra r
JOIN Modifications m ON r.id = m.RefSpectraID
WHERE r.id = 1;
```

---

## Skyline Compatibility

The blib output is designed for seamless import into Skyline:

1. **Score type**: Uses `PERCOLATOR QVALUE` (ID 14) which Skyline recognizes as a q-value
2. **File paths**: Full paths in `SpectrumSourceFiles` populate the Extract Chromatograms dialog
3. **Peak boundaries**: `startTime` and `endTime` in RefSpectra provide integration boundaries
4. **Modifications**: Standard BiblioSpec Modifications table for PTM display
5. **Proteins**: Standard BiblioSpec Proteins/RefSpectraProteins tables for grouping

### Import into Skyline

1. File → Import → Peptide Search
2. Select the `.blib` file
3. Skyline will recognize the PERCOLATOR QVALUE scores as q-values
4. The Extract Chromatograms dialog will be pre-populated with mzML file paths

---

## Implementation

The blib writer is implemented in:
- `crates/osprey-io/src/output/blib.rs` - `BlibWriter` struct

Key methods:
- `create()` - Create new blib file with schema
- `add_source_file()` - Add mzML source file
- `add_spectrum()` - Add detected peptide spectrum
- `add_modifications()` - Add modification positions
- `add_protein_mapping()` - Add protein accessions
- `add_peak_boundaries()` - Add peak boundaries
- `add_run_scores()` - Add run-level scores
- `add_experiment_scores()` - Add experiment-level scores
- `finalize()` - Update spectrum count and close

## Data Flow

How data flows from library through search to blib output:

1. **Library loading**: Spectral library loaded from blib/elib/DIA-NN TSV format
2. **Deduplication**: Entries grouped by (modified_sequence, charge); best spectrum kept per precursor
3. **Decoy generation**: Reversed-sequence decoys generated for each target
4. **Search**: Each library precursor scored against DIA data (coelution or regression mode)
5. **FDR control**: Mokapot semi-supervised learning assigns q-values
6. **Blib output**: Passing precursors written with:
   - **RefSpectra**: Peptide sequence, precursor m/z, charge, apex RT, peak boundaries
   - **RefSpectraPeaks**: Library theoretical fragment m/z and intensities (not observed DIA peaks)
   - **Modifications**: From library entry (UniMod converted to mass notation)
   - **Proteins**: From library entry protein accessions
   - **RetentionTimes**: Per-run peak boundaries for multi-file experiments
   - **Osprey tables**: Run/experiment q-values, peak boundaries, coefficient series
