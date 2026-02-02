//! BiblioSpec (blib) output writer
//!
//! Writes Osprey detection results to SQLite blib format for Skyline integration.
//! Follows the BiblioSpec schema with Osprey-specific extension tables.

use osprey_core::{OspreyError, PeakBoundaries, Result};
use rusqlite::{params, Connection};
use std::path::Path;

/// BiblioSpec library version
const BLIB_MAJOR_VERSION: i32 = 1;
const BLIB_MINOR_VERSION: i32 = 10;

/// Score type for Osprey detections
const SCORE_TYPE_OSPREY: i32 = 100; // Custom score type for Osprey

/// blib file writer
pub struct BlibWriter {
    conn: Connection,
}

impl BlibWriter {
    /// Create a new blib file at the given path
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let conn = Connection::open(path.as_ref()).map_err(|e| {
            OspreyError::SqliteError(format!("Failed to create blib file: {}", e))
        })?;

        let writer = Self { conn };
        writer.create_schema()?;

        Ok(writer)
    }

    /// Create the blib schema tables
    fn create_schema(&self) -> Result<()> {
        // Create standard BiblioSpec tables
        self.conn.execute_batch(
            r#"
            -- Library metadata
            CREATE TABLE LibInfo (
                libLSID TEXT PRIMARY KEY,
                createTime TEXT,
                numSpecs INTEGER,
                majorVersion INTEGER,
                minorVersion INTEGER
            );

            -- Score type definitions
            CREATE TABLE ScoreTypes (
                id INTEGER PRIMARY KEY,
                scoreType TEXT
            );

            -- Ion mobility type definitions
            CREATE TABLE IonMobilityTypes (
                id INTEGER PRIMARY KEY,
                ionMobilityType TEXT
            );

            -- Source files
            CREATE TABLE SpectrumSourceFiles (
                id INTEGER PRIMARY KEY,
                fileName TEXT,
                idFileName TEXT,
                cutoffScore REAL
            );

            -- Main spectrum table
            CREATE TABLE RefSpectra (
                id INTEGER PRIMARY KEY,
                peptideSeq TEXT,
                precursorMZ REAL,
                precursorCharge INTEGER,
                peptideModSeq TEXT,
                prevAA TEXT,
                nextAA TEXT,
                copies INTEGER,
                numPeaks INTEGER,
                ionMobility REAL,
                collisionalCrossSectionSqA REAL,
                ionMobilityHighEnergyOffset REAL,
                ionMobilityType INTEGER,
                retentionTime REAL,
                startTime REAL,
                endTime REAL,
                totalIonCurrent REAL,
                moleculeName TEXT,
                chemicalFormula TEXT,
                precursorAdduct TEXT,
                inchiKey TEXT,
                otherKeys TEXT,
                fileID INTEGER,
                SpecIDinFile TEXT,
                score REAL,
                scoreType INTEGER,
                FOREIGN KEY (fileID) REFERENCES SpectrumSourceFiles(id),
                FOREIGN KEY (scoreType) REFERENCES ScoreTypes(id),
                FOREIGN KEY (ionMobilityType) REFERENCES IonMobilityTypes(id)
            );

            -- Peak data (stored as blobs)
            CREATE TABLE RefSpectraPeaks (
                RefSpectraID INTEGER,
                peakMZ BLOB,
                peakIntensity BLOB,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Modifications
            CREATE TABLE Modifications (
                id INTEGER PRIMARY KEY,
                RefSpectraID INTEGER,
                position INTEGER,
                mass REAL,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Protein mappings
            CREATE TABLE Proteins (
                id INTEGER PRIMARY KEY,
                accession TEXT
            );

            CREATE TABLE RefSpectraProteins (
                RefSpectraID INTEGER,
                ProteinID INTEGER,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id),
                FOREIGN KEY (ProteinID) REFERENCES Proteins(id)
            );

            -- Osprey extension tables

            -- Peak boundaries per run
            CREATE TABLE OspreyPeakBoundaries (
                id INTEGER PRIMARY KEY,
                RefSpectraID INTEGER,
                FileName TEXT,
                StartRT REAL,
                EndRT REAL,
                ApexRT REAL,
                ApexIntensity REAL,
                IntegratedArea REAL,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Run-level scores
            CREATE TABLE OspreyRunScores (
                id INTEGER PRIMARY KEY,
                RefSpectraID INTEGER,
                FileName TEXT,
                RunQValue REAL,
                DiscriminantScore REAL,
                PosteriorErrorProb REAL,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Experiment-level scores
            CREATE TABLE OspreyExperimentScores (
                id INTEGER PRIMARY KEY,
                RefSpectraID INTEGER,
                ExperimentQValue REAL,
                NRunsDetected INTEGER,
                NRunsSearched INTEGER,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Coefficient time series (optional)
            CREATE TABLE OspreyCoefficients (
                id INTEGER PRIMARY KEY,
                RefSpectraID INTEGER,
                FileName TEXT,
                ScanNumber INTEGER,
                RT REAL,
                Coefficient REAL,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Analysis metadata
            CREATE TABLE OspreyMetadata (
                Key TEXT PRIMARY KEY,
                Value TEXT
            );

            -- Create indices for performance
            CREATE INDEX idx_refspectra_peptide ON RefSpectra(peptideSeq);
            CREATE INDEX idx_refspectra_modseq ON RefSpectra(peptideModSeq);
            CREATE INDEX idx_refspectra_mz ON RefSpectra(precursorMZ);
            CREATE INDEX idx_peaks_refid ON RefSpectraPeaks(RefSpectraID);
            CREATE INDEX idx_mods_refid ON Modifications(RefSpectraID);
            CREATE INDEX idx_boundaries_refid ON OspreyPeakBoundaries(RefSpectraID);
            CREATE INDEX idx_runscores_refid ON OspreyRunScores(RefSpectraID);
            "#,
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to create schema: {}", e)))?;

        // Insert initial metadata
        self.conn.execute(
            "INSERT INTO LibInfo (libLSID, createTime, numSpecs, majorVersion, minorVersion)
             VALUES (?, datetime('now'), 0, ?, ?)",
            params![
                format!("urn:lsid:osprey:blib:{}", uuid_simple()),
                BLIB_MAJOR_VERSION,
                BLIB_MINOR_VERSION
            ],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to insert LibInfo: {}", e)))?;

        // Insert score types
        self.conn.execute(
            "INSERT INTO ScoreTypes (id, scoreType) VALUES (?, ?)",
            params![SCORE_TYPE_OSPREY, "OSPREY q-value"],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to insert ScoreTypes: {}", e)))?;

        // Insert ion mobility types
        for (id, name) in [
            (0, "none"),
            (1, "driftTime(msec)"),
            (2, "inverseK0(Vsec/cm^2)"),
            (3, "compensation(V)"),
        ] {
            self.conn.execute(
                "INSERT INTO IonMobilityTypes (id, ionMobilityType) VALUES (?, ?)",
                params![id, name],
            ).map_err(|e| OspreyError::SqliteError(format!("Failed to insert IonMobilityTypes: {}", e)))?;
        }

        Ok(())
    }

    /// Add a source file and return its ID
    pub fn add_source_file(&self, file_name: &str) -> Result<i64> {
        self.conn.execute(
            "INSERT INTO SpectrumSourceFiles (fileName, idFileName, cutoffScore) VALUES (?, ?, ?)",
            params![file_name, file_name, 0.0],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add source file: {}", e)))?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Add a detected peptide spectrum
    pub fn add_spectrum(
        &self,
        peptide_seq: &str,
        peptide_mod_seq: &str,
        precursor_mz: f64,
        precursor_charge: i32,
        retention_time: f64,
        start_time: f64,
        end_time: f64,
        mzs: &[f64],
        intensities: &[f32],
        score: f64,
        file_id: i64,
    ) -> Result<i64> {
        // Convert peaks to blobs
        let mz_blob = f64_vec_to_blob(mzs);
        let int_blob = f32_vec_to_blob(intensities);

        self.conn.execute(
            r#"INSERT INTO RefSpectra (
                peptideSeq, precursorMZ, precursorCharge, peptideModSeq,
                prevAA, nextAA, copies, numPeaks, ionMobility,
                collisionalCrossSectionSqA, ionMobilityHighEnergyOffset, ionMobilityType,
                retentionTime, startTime, endTime, totalIonCurrent,
                moleculeName, chemicalFormula, precursorAdduct, inchiKey, otherKeys,
                fileID, SpecIDinFile, score, scoreType
            ) VALUES (
                ?, ?, ?, ?,
                '-', '-', 1, ?, 0.0,
                0.0, 0.0, 0,
                ?, ?, ?, 0.0,
                NULL, NULL, NULL, NULL, NULL,
                ?, NULL, ?, ?
            )"#,
            params![
                peptide_seq,
                precursor_mz,
                precursor_charge,
                peptide_mod_seq,
                mzs.len() as i32,
                retention_time,
                start_time,
                end_time,
                file_id,
                score,
                SCORE_TYPE_OSPREY
            ],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add spectrum: {}", e)))?;

        let ref_id = self.conn.last_insert_rowid();

        // Insert peaks
        self.conn.execute(
            "INSERT INTO RefSpectraPeaks (RefSpectraID, peakMZ, peakIntensity) VALUES (?, ?, ?)",
            params![ref_id, mz_blob, int_blob],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add peaks: {}", e)))?;

        Ok(ref_id)
    }

    /// Add peak boundaries for a detection
    pub fn add_peak_boundaries(
        &self,
        ref_id: i64,
        file_name: &str,
        boundaries: &PeakBoundaries,
    ) -> Result<()> {
        self.conn.execute(
            r#"INSERT INTO OspreyPeakBoundaries (
                RefSpectraID, FileName, StartRT, EndRT, ApexRT, ApexIntensity, IntegratedArea
            ) VALUES (?, ?, ?, ?, ?, ?, ?)"#,
            params![
                ref_id,
                file_name,
                boundaries.start_rt,
                boundaries.end_rt,
                boundaries.apex_rt,
                boundaries.apex_coefficient,
                boundaries.integrated_area
            ],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add peak boundaries: {}", e)))?;

        Ok(())
    }

    /// Add run-level scores
    pub fn add_run_scores(
        &self,
        ref_id: i64,
        file_name: &str,
        q_value: f64,
        score: f64,
        pep: f64,
    ) -> Result<()> {
        self.conn.execute(
            r#"INSERT INTO OspreyRunScores (
                RefSpectraID, FileName, RunQValue, DiscriminantScore, PosteriorErrorProb
            ) VALUES (?, ?, ?, ?, ?)"#,
            params![ref_id, file_name, q_value, score, pep],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add run scores: {}", e)))?;

        Ok(())
    }

    /// Add coefficient time series point
    pub fn add_coefficient(
        &self,
        ref_id: i64,
        file_name: &str,
        scan_number: u32,
        rt: f64,
        coefficient: f64,
    ) -> Result<()> {
        self.conn.execute(
            r#"INSERT INTO OspreyCoefficients (
                RefSpectraID, FileName, ScanNumber, RT, Coefficient
            ) VALUES (?, ?, ?, ?, ?)"#,
            params![ref_id, file_name, scan_number, rt, coefficient],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add coefficient: {}", e)))?;

        Ok(())
    }

    /// Add metadata key-value pair
    pub fn add_metadata(&self, key: &str, value: &str) -> Result<()> {
        self.conn.execute(
            "INSERT OR REPLACE INTO OspreyMetadata (Key, Value) VALUES (?, ?)",
            params![key, value],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to add metadata: {}", e)))?;

        Ok(())
    }

    /// Update the spectrum count in LibInfo
    pub fn finalize(&self) -> Result<()> {
        self.conn.execute(
            "UPDATE LibInfo SET numSpecs = (SELECT COUNT(*) FROM RefSpectra)",
            [],
        ).map_err(|e| OspreyError::SqliteError(format!("Failed to update LibInfo: {}", e)))?;

        Ok(())
    }
}

/// Convert f64 vector to byte blob (little-endian)
fn f64_vec_to_blob(values: &[f64]) -> Vec<u8> {
    let mut blob = Vec::with_capacity(values.len() * 8);
    for v in values {
        blob.extend_from_slice(&v.to_le_bytes());
    }
    blob
}

/// Convert f32 vector to byte blob (little-endian)
fn f32_vec_to_blob(values: &[f32]) -> Vec<u8> {
    let mut blob = Vec::with_capacity(values.len() * 4);
    for v in values {
        blob.extend_from_slice(&v.to_le_bytes());
    }
    blob
}

/// Generate a simple UUID-like string (not cryptographically secure)
fn uuid_simple() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    format!("{:032x}", now)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn test_create_blib() {
        let temp = NamedTempFile::new().unwrap();
        let writer = BlibWriter::create(temp.path()).unwrap();

        // Add a source file
        let file_id = writer.add_source_file("test.mzML").unwrap();
        assert!(file_id > 0);

        // Add a spectrum
        let ref_id = writer.add_spectrum(
            "PEPTIDE",
            "PEPTIDE",
            500.0,
            2,
            10.0,
            9.0,
            11.0,
            &[300.0, 400.0, 500.0],
            &[100.0, 200.0, 300.0],
            0.01,
            file_id,
        ).unwrap();
        assert!(ref_id > 0);

        // Add peak boundaries
        let boundaries = PeakBoundaries {
            start_rt: 9.0,
            end_rt: 11.0,
            apex_rt: 10.0,
            apex_coefficient: 0.5,
            integrated_area: 1.5,
            peak_quality: Default::default(),
        };
        writer.add_peak_boundaries(ref_id, "test.mzML", &boundaries).unwrap();

        // Add metadata
        writer.add_metadata("osprey_version", "0.1.0").unwrap();

        // Finalize
        writer.finalize().unwrap();
    }

    #[test]
    fn test_blob_conversion() {
        let values = vec![1.0f64, 2.0, 3.0];
        let blob = f64_vec_to_blob(&values);
        assert_eq!(blob.len(), 24); // 3 * 8 bytes

        let float_values = vec![1.0f32, 2.0, 3.0];
        let float_blob = f32_vec_to_blob(&float_values);
        assert_eq!(float_blob.len(), 12); // 3 * 4 bytes
    }
}
