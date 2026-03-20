//! BiblioSpec (blib) output writer
//!
//! Writes Osprey detection results to SQLite blib format for Skyline integration.
//! Follows the BiblioSpec schema with Osprey-specific extension tables.

use flate2::write::ZlibEncoder;
use flate2::Compression;
use osprey_core::{Modification, OspreyError, PeakBoundaries, Result};
use rusqlite::{params, Connection};
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// BiblioSpec library version
const BLIB_MAJOR_VERSION: i32 = 1;
const BLIB_MINOR_VERSION: i32 = 11;

/// Score type for Osprey detections
/// Using GENERIC Q-VALUE (19) which Skyline recognizes as a q-value
const SCORE_TYPE_GENERIC_QVALUE: i32 = 19;

/// blib file writer
pub struct BlibWriter {
    conn: Connection,
    in_transaction: bool,
    next_spec_id: i64,
    /// In-memory cache of protein accession → SQLite row ID
    protein_cache: HashMap<String, i64>,
}

impl BlibWriter {
    /// Create a new blib file at the given path
    ///
    /// If the file already exists, it will be removed and recreated.
    pub fn create<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();

        // Remove existing file if it exists to ensure a fresh database
        if path.exists() {
            std::fs::remove_file(path).map_err(|e| {
                OspreyError::SqliteError(format!(
                    "Failed to remove existing blib file '{}': {}",
                    path.display(),
                    e
                ))
            })?;
        }

        let conn = Connection::open(path)
            .map_err(|e| OspreyError::SqliteError(format!("Failed to create blib file: {}", e)))?;

        // Enable WAL mode for better write performance
        conn.pragma_update(None, "journal_mode", "WAL").ok();
        // Reduce sync frequency (still safe with WAL)
        conn.pragma_update(None, "synchronous", "NORMAL").ok();

        let writer = Self {
            conn,
            in_transaction: false,
            next_spec_id: 0,
            protein_cache: HashMap::new(),
        };
        writer.create_schema()?;

        Ok(writer)
    }

    /// Begin a batch transaction for faster writes
    ///
    /// Call this before writing multiple spectra, then call `commit()` when done.
    /// All writes between `begin_batch()` and `commit()` will be in a single transaction.
    pub fn begin_batch(&mut self) -> Result<()> {
        if self.in_transaction {
            return Ok(()); // Already in a transaction
        }
        self.conn
            .execute("BEGIN TRANSACTION", [])
            .map_err(|e| OspreyError::SqliteError(format!("Failed to begin transaction: {}", e)))?;
        self.in_transaction = true;
        Ok(())
    }

    /// Commit the batch transaction
    pub fn commit(&mut self) -> Result<()> {
        if !self.in_transaction {
            return Ok(()); // Nothing to commit
        }
        self.conn.execute("COMMIT", []).map_err(|e| {
            OspreyError::SqliteError(format!("Failed to commit transaction: {}", e))
        })?;
        self.in_transaction = false;
        Ok(())
    }

    /// Create the blib schema tables
    fn create_schema(&self) -> Result<()> {
        // Create standard BiblioSpec tables matching Skyline's expected schema
        self.conn
            .execute_batch(
                r#"
            -- Library metadata
            CREATE TABLE LibInfo (
                libLSID TEXT PRIMARY KEY,
                createTime TEXT,
                numSpecs INTEGER,
                majorVersion INTEGER,
                minorVersion INTEGER
            );

            -- Score type definitions (3 columns to match Skyline)
            CREATE TABLE ScoreTypes (
                id INTEGER PRIMARY KEY,
                scoreType VARCHAR(128),
                probabilityType VARCHAR(128)
            );

            -- Ion mobility type definitions
            CREATE TABLE IonMobilityTypes (
                id INTEGER PRIMARY KEY,
                ionMobilityType TEXT
            );

            -- Source files
            CREATE TABLE SpectrumSourceFiles (
                id INTEGER PRIMARY KEY autoincrement not null,
                fileName VARCHAR(512),
                idFileName VARCHAR(512),
                cutoffScore REAL,
                workflowType TINYINT
            );

            -- Main spectrum table
            CREATE TABLE RefSpectra (
                id INTEGER PRIMARY KEY autoincrement not null,
                peptideSeq VARCHAR(150),
                precursorMZ REAL,
                precursorCharge INTEGER,
                peptideModSeq VARCHAR(200),
                prevAA CHAR(1),
                nextAA CHAR(1),
                copies INTEGER,
                numPeaks INTEGER,
                ionMobility REAL,
                collisionalCrossSectionSqA REAL,
                ionMobilityHighEnergyOffset REAL,
                ionMobilityType TINYINT,
                retentionTime REAL,
                startTime REAL,
                endTime REAL,
                totalIonCurrent REAL,
                moleculeName VARCHAR(128),
                chemicalFormula VARCHAR(128),
                precursorAdduct VARCHAR(128),
                inchiKey VARCHAR(128),
                otherKeys VARCHAR(128),
                fileID INTEGER,
                SpecIDinFile VARCHAR(256),
                score REAL,
                scoreType TINYINT
            );

            -- Peak data (stored as blobs)
            CREATE TABLE RefSpectraPeaks (
                RefSpectraID INTEGER,
                peakMZ BLOB,
                peakIntensity BLOB,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Peak annotations (required by Skyline BlibBuild)
            CREATE TABLE RefSpectraPeakAnnotations (
                id INTEGER PRIMARY KEY,
                RefSpectraID INTEGER,
                peakIndex INTEGER,
                name TEXT,
                formula TEXT,
                inchiKey TEXT,
                otherKeys TEXT,
                charge INTEGER,
                adduct TEXT,
                comment TEXT,
                mzTheoretical REAL,
                mzObserved REAL,
                FOREIGN KEY (RefSpectraID) REFERENCES RefSpectra(id)
            );

            -- Modifications (1-based positions for Skyline)
            CREATE TABLE Modifications (
                id INTEGER PRIMARY KEY autoincrement not null,
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

            -- Per-run retention times and peak boundaries (Skyline standard table)
            CREATE TABLE RetentionTimes (
                RefSpectraID INTEGER,
                RedundantRefSpectraID INTEGER,
                SpectrumSourceID INTEGER,
                ionMobility REAL,
                collisionalCrossSectionSqA REAL,
                ionMobilityHighEnergyOffset REAL,
                ionMobilityType TINYINT,
                retentionTime REAL,
                startTime REAL,
                endTime REAL,
                score REAL,
                bestSpectrum INTEGER,
                FOREIGN KEY(RefSpectraID) REFERENCES RefSpectra(id)
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

            "#,
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to create schema: {}", e)))?;

        // Insert initial metadata
        self.conn
            .execute(
                "INSERT INTO LibInfo (libLSID, createTime, numSpecs, majorVersion, minorVersion)
             VALUES (?, datetime('now'), 0, ?, ?)",
                params![
                    format!("urn:lsid:osprey:blib:{}", uuid_simple()),
                    BLIB_MAJOR_VERSION,
                    BLIB_MINOR_VERSION
                ],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to insert LibInfo: {}", e)))?;

        // Insert all standard BiblioSpec score types (matching Skyline's schema)
        let score_types = [
            (0, "UNKNOWN", "NOT_A_PROBABILITY_VALUE"),
            (
                1,
                "PERCOLATOR QVALUE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                2,
                "PEPTIDE PROPHET SOMETHING",
                "PROBABILITY_THAT_IDENTIFICATION_IS_CORRECT",
            ),
            (3, "SPECTRUM MILL", "NOT_A_PROBABILITY_VALUE"),
            (
                4,
                "IDPICKER FDR",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                5,
                "MASCOT IONS SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                6,
                "TANDEM EXPECTATION VALUE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                7,
                "PROTEIN PILOT CONFIDENCE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_CORRECT",
            ),
            (
                8,
                "SCAFFOLD SOMETHING",
                "PROBABILITY_THAT_IDENTIFICATION_IS_CORRECT",
            ),
            (9, "WATERS MSE PEPTIDE SCORE", "NOT_A_PROBABILITY_VALUE"),
            (
                10,
                "OMSSA EXPECTATION SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                11,
                "PROTEIN PROSPECTOR EXPECTATION SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                12,
                "SEQUEST XCORR",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                13,
                "MAXQUANT SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                14,
                "MORPHEUS SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                15,
                "MSGF+ SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                16,
                "PEAKS CONFIDENCE SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                17,
                "BYONIC SCORE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                18,
                "PEPTIDE SHAKER CONFIDENCE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_CORRECT",
            ),
            (
                19,
                "GENERIC Q-VALUE",
                "PROBABILITY_THAT_IDENTIFICATION_IS_INCORRECT",
            ),
            (
                20,
                "HARDKLOR IDOTP",
                "PROBABILITY_THAT_IDENTIFICATION_IS_CORRECT",
            ),
        ];
        for (id, name, prob_type) in score_types {
            self.conn
                .execute(
                    "INSERT INTO ScoreTypes (id, scoreType, probabilityType) VALUES (?, ?, ?)",
                    params![id, name, prob_type],
                )
                .map_err(|e| {
                    OspreyError::SqliteError(format!("Failed to insert ScoreTypes: {}", e))
                })?;
        }

        // Insert ion mobility types
        for (id, name) in [
            (0, "none"),
            (1, "driftTime(msec)"),
            (2, "inverseK0(Vsec/cm^2)"),
            (3, "compensation(V)"),
        ] {
            self.conn
                .execute(
                    "INSERT INTO IonMobilityTypes (id, ionMobilityType) VALUES (?, ?)",
                    params![id, name],
                )
                .map_err(|e| {
                    OspreyError::SqliteError(format!("Failed to insert IonMobilityTypes: {}", e))
                })?;
        }

        Ok(())
    }

    /// Add a source file and return its ID
    ///
    /// The `id_file_name` is the identification file (e.g., library or search output).
    /// The `cutoff_score` is the FDR cutoff used for filtering.
    pub fn add_source_file(
        &self,
        file_name: &str,
        id_file_name: &str,
        cutoff_score: f64,
    ) -> Result<i64> {
        self.conn
            .execute(
                "INSERT INTO SpectrumSourceFiles (fileName, idFileName, cutoffScore, workflowType) VALUES (?, ?, ?, 1)",
                params![file_name, id_file_name, cutoff_score],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to add source file: {}", e)))?;

        Ok(self.conn.last_insert_rowid())
    }

    /// Add a detected peptide spectrum
    ///
    /// The `score` should be the raw q-value (lower is better).
    /// The `copies` indicates how many source files this precursor was detected in.
    #[allow(clippy::too_many_arguments)]
    pub fn add_spectrum(
        &mut self,
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
        copies: i32,
        total_ion_current: f64,
    ) -> Result<i64> {
        // Clean sequences - strip flanking characters (underscores, periods, dashes)
        // that are used in some library formats for flanking amino acids
        let clean_seq = strip_flanking_chars(peptide_seq);
        // Also convert UniMod notation (e.g., [UniMod:4]) to mass notation (e.g., [+57.0215])
        // because Skyline expects mass values, not UniMod references
        let clean_mod_seq = convert_unimod_to_mass(&strip_flanking_chars(peptide_mod_seq));

        // Convert peaks to blobs - zlib-compressed when it reduces size
        let mz_blob = compress_bytes(f64_vec_to_bytes(mzs));
        let int_blob = compress_bytes(f32_vec_to_bytes(intensities));

        let spec_id_in_file = self.next_spec_id.to_string();
        self.next_spec_id += 1;

        self.conn
            .execute(
                r#"INSERT INTO RefSpectra (
                peptideSeq, precursorMZ, precursorCharge, peptideModSeq,
                prevAA, nextAA, copies, numPeaks, ionMobility,
                collisionalCrossSectionSqA, ionMobilityHighEnergyOffset, ionMobilityType,
                retentionTime, startTime, endTime, totalIonCurrent,
                moleculeName, chemicalFormula, precursorAdduct, inchiKey, otherKeys,
                fileID, SpecIDinFile, score, scoreType
            ) VALUES (
                ?, ?, ?, ?,
                '-', '-', ?, ?, 0.0,
                0.0, 0.0, 0,
                ?, ?, ?, ?,
                '', '', '', '', '',
                ?, ?, ?, ?
            )"#,
                params![
                    clean_seq,
                    precursor_mz,
                    precursor_charge,
                    clean_mod_seq,
                    copies,
                    mzs.len() as i32,
                    retention_time,
                    start_time,
                    end_time,
                    total_ion_current,
                    file_id,
                    spec_id_in_file,
                    score,
                    SCORE_TYPE_GENERIC_QVALUE
                ],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to add spectrum: {}", e)))?;

        let ref_id = self.conn.last_insert_rowid();

        // Insert peaks
        self.conn
            .execute(
                "INSERT INTO RefSpectraPeaks (RefSpectraID, peakMZ, peakIntensity) VALUES (?, ?, ?)",
                params![ref_id, mz_blob, int_blob],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to add peaks: {}", e)))?;

        Ok(ref_id)
    }

    /// Add modifications for a spectrum
    ///
    /// Writes modification positions and masses to the Modifications table.
    /// Positions are converted to 1-based indexing as expected by Skyline.
    pub fn add_modifications(&self, ref_id: i64, modifications: &[Modification]) -> Result<()> {
        for modif in modifications {
            // Skyline expects 1-based modification positions
            let position_1based = modif.position as i32 + 1;
            self.conn
                .execute(
                    "INSERT INTO Modifications (RefSpectraID, position, mass) VALUES (?, ?, ?)",
                    params![ref_id, position_1based, modif.mass_delta],
                )
                .map_err(|e| {
                    OspreyError::SqliteError(format!("Failed to add modification: {}", e))
                })?;
        }
        Ok(())
    }

    /// Add protein mapping for a spectrum
    ///
    /// Writes protein accessions to the Proteins table and creates the mapping
    /// in RefSpectraProteins. Uses an in-memory cache to avoid repeated SQL
    /// lookups for the same protein accession.
    pub fn add_protein_mapping(&mut self, ref_id: i64, protein_ids: &[String]) -> Result<()> {
        for accession in protein_ids {
            if accession.is_empty() {
                continue;
            }

            let protein_id = if let Some(&id) = self.protein_cache.get(accession) {
                id
            } else {
                // Insert new protein
                self.conn
                    .execute(
                        "INSERT INTO Proteins (accession) VALUES (?)",
                        params![accession],
                    )
                    .map_err(|e| {
                        OspreyError::SqliteError(format!("Failed to add protein: {}", e))
                    })?;
                let id = self.conn.last_insert_rowid();
                self.protein_cache.insert(accession.clone(), id);
                id
            };

            // Create mapping
            self.conn
                .execute(
                    "INSERT INTO RefSpectraProteins (RefSpectraID, ProteinID) VALUES (?, ?)",
                    params![ref_id, protein_id],
                )
                .map_err(|e| {
                    OspreyError::SqliteError(format!("Failed to add protein mapping: {}", e))
                })?;
        }
        Ok(())
    }

    /// Add a retention time entry for per-run peak boundaries (Skyline standard table)
    ///
    /// Each precursor gets one RetentionTimes row per source file where it was detected.
    /// The `best_spectrum` flag marks which run's data was used for the RefSpectra entry.
    ///
    /// `retention_time` should be `Some(apex_rt)` for runs where the precursor passed
    /// run-level FDR (Skyline shows an ID line), or `None` for runs where it did not
    /// (Skyline uses startTime/endTime for quantification without showing an ID line).
    #[allow(clippy::too_many_arguments)]
    pub fn add_retention_time(
        &self,
        ref_id: i64,
        source_file_id: i64,
        retention_time: Option<f64>,
        start_time: f64,
        end_time: f64,
        score: f64,
        best_spectrum: bool,
    ) -> Result<()> {
        self.conn
            .execute(
                r#"INSERT INTO RetentionTimes (
                RefSpectraID, RedundantRefSpectraID, SpectrumSourceID,
                ionMobility, collisionalCrossSectionSqA,
                ionMobilityHighEnergyOffset, ionMobilityType,
                retentionTime, startTime, endTime, score, bestSpectrum
            ) VALUES (?, 0, ?, 0.0, 0.0, 0.0, 0, ?, ?, ?, ?, ?)"#,
                params![
                    ref_id,
                    source_file_id,
                    retention_time,
                    start_time,
                    end_time,
                    score,
                    if best_spectrum { 1 } else { 0 }
                ],
            )
            .map_err(|e| {
                OspreyError::SqliteError(format!("Failed to add retention time: {}", e))
            })?;

        Ok(())
    }

    /// Add peak boundaries for a detection
    pub fn add_peak_boundaries(
        &self,
        ref_id: i64,
        file_name: &str,
        boundaries: &PeakBoundaries,
    ) -> Result<()> {
        self.conn
            .execute(
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
            )
            .map_err(|e| {
                OspreyError::SqliteError(format!("Failed to add peak boundaries: {}", e))
            })?;

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
        self.conn
            .execute(
                r#"INSERT INTO OspreyRunScores (
                RefSpectraID, FileName, RunQValue, DiscriminantScore, PosteriorErrorProb
            ) VALUES (?, ?, ?, ?, ?)"#,
                params![ref_id, file_name, q_value, score, pep],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to add run scores: {}", e)))?;

        Ok(())
    }

    /// Add experiment-level scores (one per RefSpectra entry)
    pub fn add_experiment_scores(
        &self,
        ref_id: i64,
        experiment_q_value: f64,
        n_runs_detected: i32,
        n_runs_searched: i32,
    ) -> Result<()> {
        self.conn
            .execute(
                r#"INSERT INTO OspreyExperimentScores (
                RefSpectraID, ExperimentQValue, NRunsDetected, NRunsSearched
            ) VALUES (?, ?, ?, ?)"#,
                params![ref_id, experiment_q_value, n_runs_detected, n_runs_searched],
            )
            .map_err(|e| {
                OspreyError::SqliteError(format!("Failed to add experiment scores: {}", e))
            })?;

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
        self.conn
            .execute(
                r#"INSERT INTO OspreyCoefficients (
                RefSpectraID, FileName, ScanNumber, RT, Coefficient
            ) VALUES (?, ?, ?, ?, ?)"#,
                params![ref_id, file_name, scan_number, rt, coefficient],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to add coefficient: {}", e)))?;

        Ok(())
    }

    /// Add metadata key-value pair
    pub fn add_metadata(&self, key: &str, value: &str) -> Result<()> {
        self.conn
            .execute(
                "INSERT OR REPLACE INTO OspreyMetadata (Key, Value) VALUES (?, ?)",
                params![key, value],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to add metadata: {}", e)))?;

        Ok(())
    }

    /// Update the spectrum count in LibInfo and checkpoint the WAL
    pub fn finalize(&self) -> Result<()> {
        self.conn
            .execute(
                "UPDATE LibInfo SET numSpecs = (SELECT COUNT(*) FROM RefSpectra)",
                [],
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to update LibInfo: {}", e)))?;

        // Create indices after all data is inserted (faster than indexing during inserts)
        self.conn
            .execute_batch(
                r#"
            CREATE INDEX idx_refspectra_peptide ON RefSpectra(peptideSeq);
            CREATE INDEX idx_refspectra_modseq ON RefSpectra(peptideModSeq);
            CREATE INDEX idx_refspectra_mz ON RefSpectra(precursorMZ);
            CREATE INDEX idx_peaks_refid ON RefSpectraPeaks(RefSpectraID);
            CREATE INDEX idx_mods_refid ON Modifications(RefSpectraID);
            CREATE INDEX idx_boundaries_refid ON OspreyPeakBoundaries(RefSpectraID);
            CREATE INDEX idx_runscores_refid ON OspreyRunScores(RefSpectraID);
            CREATE INDEX idx_rettimes_refid ON RetentionTimes(RefSpectraID);
            "#,
            )
            .map_err(|e| OspreyError::SqliteError(format!("Failed to create indices: {}", e)))?;

        // Checkpoint WAL into the main database file and switch to DELETE journal mode
        // so the -wal and -shm files are removed, leaving a single clean .blib file
        self.conn
            .pragma_update(None, "wal_checkpoint", "TRUNCATE")
            .ok();
        self.conn.pragma_update(None, "journal_mode", "DELETE").ok();

        Ok(())
    }
}

impl Drop for BlibWriter {
    fn drop(&mut self) {
        // Rollback any uncommitted transaction on drop
        if self.in_transaction {
            let _ = self.conn.execute("ROLLBACK", []);
        }
    }
}

/// Zlib-compress a byte buffer, returning raw bytes if compression doesn't help.
///
/// BiblioSpec readers determine whether a blob is compressed by comparing its
/// length to the expected uncompressed size (numPeaks * sizeof(type)). If the
/// sizes match, the blob is treated as raw bytes. This means we must store the
/// uncompressed bytes whenever compression does not reduce the size.
fn compress_bytes(raw: Vec<u8>) -> Vec<u8> {
    let mut encoder = ZlibEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(&raw).expect("zlib compression failed");
    let compressed = encoder.finish().expect("zlib finish failed");
    if compressed.len() >= raw.len() {
        raw
    } else {
        compressed
    }
}

/// Convert f64 vector to little-endian bytes
fn f64_vec_to_bytes(values: &[f64]) -> Vec<u8> {
    let mut blob = Vec::with_capacity(values.len() * 8);
    for v in values {
        blob.extend_from_slice(&v.to_le_bytes());
    }
    blob
}

/// Convert f32 vector to little-endian bytes
fn f32_vec_to_bytes(values: &[f32]) -> Vec<u8> {
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

/// Strip flanking characters from peptide sequences
///
/// Some library formats (DIA-NN, Prosit) use formats like "_.PEPTIDE._" or "_PEPTIDE_"
/// where underscores/periods represent flanking amino acids or N/C-terminal markers.
/// BiblioSpec expects clean peptide sequences without these flanking characters.
fn strip_flanking_chars(seq: &str) -> String {
    let trimmed = seq.trim_matches(|c| c == '_' || c == '.' || c == '-');

    // Also handle internal patterns like "K.PEPTIDE.R" -> "PEPTIDE"
    // But ignore dots inside modification brackets like [+57.02146]
    let first_dot = find_dot_outside_brackets(trimmed, false);
    let last_dot = find_dot_outside_brackets(trimmed, true);
    if let (Some(start), Some(end)) = (first_dot, last_dot) {
        if start < end {
            return trimmed[start + 1..end].to_string();
        }
    }

    trimmed.to_string()
}

/// Find the first (or last) dot that is NOT inside square brackets.
fn find_dot_outside_brackets(s: &str, reverse: bool) -> Option<usize> {
    let bytes = s.as_bytes();
    if reverse {
        let mut depth: i32 = 0;
        for i in (0..bytes.len()).rev() {
            if bytes[i] == b']' {
                depth += 1;
            } else if bytes[i] == b'[' {
                depth -= 1;
            } else if bytes[i] == b'.' && depth == 0 {
                return Some(i);
            }
        }
        return None;
    }
    let mut depth: i32 = 0;
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'[' {
            depth += 1;
        } else if b == b']' {
            depth -= 1;
        } else if b == b'.' && depth == 0 {
            return Some(i);
        }
    }
    None
}

/// Convert UniMod notation to mass notation in modified sequences
///
/// Skyline expects mass notation like `PEPTC[+57.021]IDE` but some libraries
/// use UniMod notation like `PEPTC[UniMod:4]IDE`. This function converts
/// UniMod references to their corresponding mass shifts.
fn convert_unimod_to_mass(seq: &str) -> String {
    let mut result = String::with_capacity(seq.len());
    let mut chars = seq.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '[' {
            // Look for UniMod pattern
            let mut bracket_content = String::new();
            while let Some(&next) = chars.peek() {
                if next == ']' {
                    chars.next();
                    break;
                }
                bracket_content.push(chars.next().unwrap());
            }

            // Check if this is a UniMod reference
            if bracket_content.starts_with("UniMod:") || bracket_content.starts_with("UNIMOD:") {
                let id_str = bracket_content
                    .trim_start_matches("UniMod:")
                    .trim_start_matches("UNIMOD:");
                if let Ok(id) = id_str.parse::<u32>() {
                    if let Some(mass) = unimod_id_to_mass(id) {
                        // Convert to mass notation
                        result.push('[');
                        if mass >= 0.0 {
                            result.push_str(&format!("+{:.4}", mass));
                        } else {
                            result.push_str(&format!("{:.4}", mass));
                        }
                        result.push(']');
                        continue;
                    }
                }
            }

            // Not a recognized UniMod, keep as-is
            result.push('[');
            result.push_str(&bracket_content);
            result.push(']');
        } else {
            result.push(c);
        }
    }

    result
}

/// Look up monoisotopic mass for a UniMod accession ID
///
/// Returns the mass delta for common UniMod modifications.
/// Used by both the blib writer (sequence conversion) and library loaders (modification parsing).
pub fn unimod_id_to_mass(id: u32) -> Option<f64> {
    match id {
        1 => Some(42.010565),    // Acetyl
        4 => Some(57.021464),    // Carbamidomethyl
        5 => Some(43.005814),    // Carbamyl
        7 => Some(0.984016),     // Deamidated
        21 => Some(79.966331),   // Phospho
        28 => Some(-18.010565),  // Glu->pyro-Glu (loss of water)
        34 => Some(14.015650),   // Methyl
        35 => Some(15.994915),   // Oxidation
        36 => Some(28.031300),   // Dimethyl
        37 => Some(42.046950),   // Trimethyl
        121 => Some(114.042927), // Ubiquitin (GlyGly)
        122 => Some(383.228102), // SUMO
        214 => Some(44.985078),  // Nitro
        312 => Some(-17.026549), // Ammonia loss (Gln->pyro-Glu)
        385 => Some(229.162932), // TMT6plex
        737 => Some(229.162932), // TMT6plex (alternate ID)
        747 => Some(304.207146), // TMTpro
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    /// Verifies that a blib file can be created with a source file, spectrum, peak boundaries, and metadata.
    #[test]
    fn test_create_blib() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = BlibWriter::create(temp.path()).unwrap();

        // Add a source file
        let file_id = writer
            .add_source_file("test.mzML", "library.tsv", 0.01)
            .unwrap();
        assert!(file_id > 0);

        // Add a spectrum (score is raw q-value, lower is better)
        let ref_id = writer
            .add_spectrum(
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
                1,
                0.0,
            )
            .unwrap();
        assert!(ref_id > 0);

        // Add retention time entry
        writer
            .add_retention_time(ref_id, file_id, Some(10.0), 9.0, 11.0, 0.01, true)
            .unwrap();

        // Add peak boundaries
        let boundaries = PeakBoundaries {
            start_rt: 9.0,
            end_rt: 11.0,
            apex_rt: 10.0,
            apex_coefficient: 0.5,
            integrated_area: 1.5,
            peak_quality: Default::default(),
        };
        writer
            .add_peak_boundaries(ref_id, "test.mzML", &boundaries)
            .unwrap();

        // Add metadata
        writer.add_metadata("osprey_version", "0.1.0").unwrap();

        // Finalize
        writer.finalize().unwrap();
    }

    /// Verifies that compress_bytes produces valid zlib or falls back to raw.
    #[test]
    fn test_blob_compression() {
        use flate2::read::ZlibDecoder;
        use std::io::Read;

        // f64 blob: with enough values, compression should reduce size
        let values: Vec<f64> = (0..100).map(|i| i as f64 * 1.5).collect();
        let raw = f64_vec_to_bytes(&values);
        let blob = compress_bytes(f64_vec_to_bytes(&values));
        if blob.len() < raw.len() {
            // Compressed — verify we can decompress back
            let mut decoder = ZlibDecoder::new(&blob[..]);
            let mut decompressed = Vec::new();
            decoder.read_to_end(&mut decompressed).unwrap();
            assert_eq!(decompressed.len(), raw.len());
        } else {
            // Fell back to raw — sizes must match
            assert_eq!(blob.len(), raw.len());
        }

        // f32 blob: same logic
        let float_values: Vec<f32> = (0..100).map(|i| i as f32 * 0.5).collect();
        let raw_f32 = f32_vec_to_bytes(&float_values);
        let blob_f32 = compress_bytes(f32_vec_to_bytes(&float_values));
        if blob_f32.len() < raw_f32.len() {
            let mut decoder = ZlibDecoder::new(&blob_f32[..]);
            let mut decompressed = Vec::new();
            decoder.read_to_end(&mut decompressed).unwrap();
            assert_eq!(decompressed.len(), raw_f32.len());
        } else {
            assert_eq!(blob_f32.len(), raw_f32.len());
        }
    }

    /// Verifies that modifications are written with 1-based positions.
    #[test]
    fn test_modifications_1based() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = BlibWriter::create(temp.path()).unwrap();
        let file_id = writer
            .add_source_file("test.mzML", "library.tsv", 0.01)
            .unwrap();
        let ref_id = writer
            .add_spectrum(
                "PEPTCIDE",
                "PEPTC[+57.021]IDE",
                500.0,
                2,
                10.0,
                9.0,
                11.0,
                &[300.0],
                &[100.0],
                0.01,
                file_id,
                1,
                0.0,
            )
            .unwrap();

        let mods = vec![Modification {
            position: 4, // 0-based: C at index 4
            unimod_id: Some(4),
            mass_delta: 57.021464,
            name: Some("Carbamidomethyl".to_string()),
        }];
        writer.add_modifications(ref_id, &mods).unwrap();

        // Verify stored as 1-based
        let stored_pos: i32 = writer
            .conn
            .query_row(
                "SELECT position FROM Modifications WHERE RefSpectraID = ?",
                params![ref_id],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(stored_pos, 5); // 0-based 4 → 1-based 5
    }

    /// Verifies that flanking amino acid characters, underscores, periods, and dashes are stripped from peptide sequences.
    #[test]
    fn test_strip_flanking_chars() {
        // Simple case - no flanking chars
        assert_eq!(strip_flanking_chars("PEPTIDE"), "PEPTIDE");

        // Underscores at start/end (DIA-NN format)
        assert_eq!(strip_flanking_chars("_PEPTIDE_"), "PEPTIDE");

        // Periods/underscores (Prosit format)
        assert_eq!(strip_flanking_chars("_.PEPTIDE._"), "PEPTIDE");

        // Flanking amino acids with periods (standard notation)
        assert_eq!(strip_flanking_chars("K.PEPTIDE.R"), "PEPTIDE");

        // Just leading/trailing periods
        assert_eq!(strip_flanking_chars(".PEPTIDE."), "PEPTIDE");

        // Modified sequence with flanking chars
        assert_eq!(
            strip_flanking_chars("_PEPTC[+57.021]IDE_"),
            "PEPTC[+57.021]IDE"
        );

        // Dashes at ends
        assert_eq!(strip_flanking_chars("-PEPTIDE-"), "PEPTIDE");

        // Dots inside modification brackets must NOT be treated as flanking separators
        // This was the bug: C[+57.02146]GPSPC[+57.02146]GLLK was truncated to 02146]GPSPC[+57
        assert_eq!(
            strip_flanking_chars("C[+57.02146]GPSPC[+57.02146]GLLK"),
            "C[+57.02146]GPSPC[+57.02146]GLLK"
        );

        // Multiple modifications with decimal dots
        assert_eq!(
            strip_flanking_chars("AC[+57.02146]AC[+57.02146]AHGMLAEDGASC[+57.02146]R"),
            "AC[+57.02146]AC[+57.02146]AHGMLAEDGASC[+57.02146]R"
        );

        // Flanking chars with modification dots — flanking stripped, mod dots preserved
        assert_eq!(
            strip_flanking_chars("K.C[+57.02146]PEPTIDE.R"),
            "C[+57.02146]PEPTIDE"
        );

        // N-terminal modification with dot
        assert_eq!(
            strip_flanking_chars("[+42.011]PEPTC[+57.02146]IDE"),
            "[+42.011]PEPTC[+57.02146]IDE"
        );

        // Single modification with decimal
        assert_eq!(
            strip_flanking_chars("PEPTM[+15.995]IDE"),
            "PEPTM[+15.995]IDE"
        );
    }

    /// Verifies that UniMod notation in modified sequences is converted to mass shift notation for Skyline compatibility.
    #[test]
    fn test_convert_unimod_to_mass() {
        // No modifications
        assert_eq!(convert_unimod_to_mass("PEPTIDE"), "PEPTIDE");

        // Already mass notation - should be unchanged
        assert_eq!(
            convert_unimod_to_mass("PEPTC[+57.021]IDE"),
            "PEPTC[+57.021]IDE"
        );

        // UniMod:4 (Carbamidomethyl)
        let result = convert_unimod_to_mass("PEPTC[UniMod:4]IDE");
        assert!(result.starts_with("PEPTC[+57."));
        assert!(result.ends_with("]IDE"));

        // UniMod:35 (Oxidation)
        let result = convert_unimod_to_mass("PEPTM[UniMod:35]IDE");
        assert!(result.starts_with("PEPTM[+15."));
        assert!(result.ends_with("]IDE"));

        // Multiple modifications
        let result = convert_unimod_to_mass("PEPTC[UniMod:4]M[UniMod:35]IDE");
        assert!(result.contains("[+57."));
        assert!(result.contains("[+15."));

        // Unknown UniMod - should be kept as-is
        assert_eq!(
            convert_unimod_to_mass("PEPTX[UniMod:99999]IDE"),
            "PEPTX[UniMod:99999]IDE"
        );

        // Case insensitive UNIMOD
        let result = convert_unimod_to_mass("PEPTC[UNIMOD:4]IDE");
        assert!(result.starts_with("PEPTC[+57."));
    }

    /// Helper: decompress a zlib-or-raw blob back to bytes
    fn decompress_blob(blob: &[u8], expected_len: usize) -> Vec<u8> {
        if blob.len() == expected_len {
            // Raw (uncompressed) — compression didn't help
            blob.to_vec()
        } else {
            use flate2::read::ZlibDecoder;
            use std::io::Read;
            let mut decoder = ZlibDecoder::new(blob);
            let mut decompressed = Vec::new();
            decoder
                .read_to_end(&mut decompressed)
                .expect("zlib decompression failed");
            decompressed
        }
    }

    /// Helper: convert bytes back to Vec<f64>
    fn bytes_to_f64_vec(bytes: &[u8]) -> Vec<f64> {
        bytes
            .chunks_exact(8)
            .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
            .collect()
    }

    /// Helper: convert bytes back to Vec<f32>
    fn bytes_to_f32_vec(bytes: &[u8]) -> Vec<f32> {
        bytes
            .chunks_exact(4)
            .map(|chunk| f32::from_le_bytes(chunk.try_into().unwrap()))
            .collect()
    }

    /// Verifies that fragment m/z and intensity values survive the write-read round trip.
    #[test]
    fn test_fragment_roundtrip() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = BlibWriter::create(temp.path()).unwrap();
        let file_id = writer
            .add_source_file("test.mzML", "library.tsv", 0.01)
            .unwrap();

        // Known b/y ion m/z values (y1, y2, y3 for a real peptide)
        let mzs = [175.1190, 288.2030, 375.2351, 488.3191, 601.4032];
        let intensities: [f32; 5] = [100.0, 80.5, 60.0, 45.2, 30.0];

        let ref_id = writer
            .add_spectrum(
                "PEPTIDE",
                "PEPTIDE",
                400.0,
                2,
                10.0,
                9.0,
                11.0,
                &mzs,
                &intensities,
                0.005,
                file_id,
                1,
                0.0,
            )
            .unwrap();

        // Read back from SQLite
        let (mz_blob, int_blob, num_peaks): (Vec<u8>, Vec<u8>, i32) = writer
            .conn
            .query_row(
                "SELECT peakMZ, peakIntensity, r.numPeaks FROM RefSpectraPeaks p JOIN RefSpectra r ON r.id = p.RefSpectraID WHERE p.RefSpectraID = ?",
                params![ref_id],
                |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)),
            )
            .unwrap();

        assert_eq!(num_peaks, 5);

        // Decompress and verify m/z values
        let mz_bytes = decompress_blob(&mz_blob, mzs.len() * 8);
        let read_mzs = bytes_to_f64_vec(&mz_bytes);
        assert_eq!(read_mzs.len(), mzs.len());
        for (got, expected) in read_mzs.iter().zip(mzs.iter()) {
            assert!(
                (got - expected).abs() < 1e-10,
                "m/z mismatch: got {}, expected {}",
                got,
                expected
            );
        }

        // Decompress and verify intensities
        let int_bytes = decompress_blob(&int_blob, intensities.len() * 4);
        let read_intensities = bytes_to_f32_vec(&int_bytes);
        assert_eq!(read_intensities.len(), intensities.len());
        for (got, expected) in read_intensities.iter().zip(intensities.iter()) {
            assert!(
                (got - expected).abs() < 1e-5,
                "intensity mismatch: got {}, expected {}",
                got,
                expected
            );
        }
    }

    /// Verifies that UniMod notation in modified sequences is converted to mass notation in the database.
    #[test]
    fn test_modseq_roundtrip() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = BlibWriter::create(temp.path()).unwrap();
        let file_id = writer
            .add_source_file("test.mzML", "library.tsv", 0.01)
            .unwrap();

        let ref_id = writer
            .add_spectrum(
                "PEPTCIDE",
                "PEPTC[UniMod:4]IDE",
                500.0,
                2,
                10.0,
                9.0,
                11.0,
                &[300.0],
                &[100.0],
                0.01,
                file_id,
                1,
                0.0,
            )
            .unwrap();

        // Read back peptideModSeq — should be mass notation, not UniMod
        let stored_modseq: String = writer
            .conn
            .query_row(
                "SELECT peptideModSeq FROM RefSpectra WHERE id = ?",
                params![ref_id],
                |row| row.get(0),
            )
            .unwrap();

        assert!(
            stored_modseq.contains("[+57."),
            "Expected mass notation, got: {}",
            stored_modseq
        );
        assert!(
            !stored_modseq.contains("UniMod"),
            "UniMod notation should be converted: {}",
            stored_modseq
        );
    }

    /// Verifies that protein mappings are correctly stored and retrievable.
    #[test]
    fn test_protein_mapping_roundtrip() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = BlibWriter::create(temp.path()).unwrap();
        let file_id = writer
            .add_source_file("test.mzML", "library.tsv", 0.01)
            .unwrap();

        let ref_id = writer
            .add_spectrum(
                "PEPTIDE",
                "PEPTIDE",
                400.0,
                2,
                10.0,
                9.0,
                11.0,
                &[300.0],
                &[100.0],
                0.01,
                file_id,
                1,
                0.0,
            )
            .unwrap();

        let proteins = vec![
            "sp|P12345|PROT_HUMAN".to_string(),
            "sp|P67890|PROT2_HUMAN".to_string(),
        ];
        writer.add_protein_mapping(ref_id, &proteins).unwrap();

        // Read back protein accessions via join
        let mut stmt = writer
            .conn
            .prepare(
                "SELECT p.accession FROM Proteins p \
                 JOIN RefSpectraProteins rp ON p.id = rp.ProteinID \
                 WHERE rp.RefSpectraID = ? ORDER BY p.accession",
            )
            .unwrap();
        let accessions: Vec<String> = stmt
            .query_map(params![ref_id], |row| row.get(0))
            .unwrap()
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(accessions.len(), 2);
        assert_eq!(accessions[0], "sp|P12345|PROT_HUMAN");
        assert_eq!(accessions[1], "sp|P67890|PROT2_HUMAN");
    }

    /// Verifies that multi-run retention times are stored with correct best_spectrum flag.
    #[test]
    fn test_retention_times_multirun() {
        let temp = NamedTempFile::new().unwrap();
        let mut writer = BlibWriter::create(temp.path()).unwrap();

        let file_id_1 = writer
            .add_source_file("run1.mzML", "library.tsv", 0.01)
            .unwrap();
        let file_id_2 = writer
            .add_source_file("run2.mzML", "library.tsv", 0.01)
            .unwrap();
        let file_id_3 = writer
            .add_source_file("run3.mzML", "library.tsv", 0.01)
            .unwrap();

        let ref_id = writer
            .add_spectrum(
                "PEPTIDE",
                "PEPTIDE",
                400.0,
                2,
                10.5,
                9.5,
                11.5,
                &[300.0, 400.0],
                &[100.0, 80.0],
                0.005,
                file_id_1,
                3,
                0.0,
            )
            .unwrap();

        // Add retention times for 3 runs — run 2 is best, run 3 has null RT (didn't pass FDR)
        writer
            .add_retention_time(ref_id, file_id_1, Some(10.5), 9.5, 11.5, 0.008, false)
            .unwrap();
        writer
            .add_retention_time(ref_id, file_id_2, Some(10.3), 9.3, 11.3, 0.005, true)
            .unwrap();
        writer
            .add_retention_time(ref_id, file_id_3, None, 9.7, 11.7, 0.012, false)
            .unwrap();

        // Read back RetentionTimes (retentionTime can be NULL)
        let mut stmt = writer
            .conn
            .prepare(
                "SELECT SpectrumSourceID, retentionTime, startTime, endTime, bestSpectrum \
                 FROM RetentionTimes WHERE RefSpectraID = ? \
                 ORDER BY SpectrumSourceID",
            )
            .unwrap();
        let rows: Vec<(i64, Option<f64>, f64, f64, i32)> = stmt
            .query_map(params![ref_id], |row| {
                Ok((
                    row.get(0)?,
                    row.get(1)?,
                    row.get(2)?,
                    row.get(3)?,
                    row.get(4)?,
                ))
            })
            .unwrap()
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(rows.len(), 3);
        // Run 1: has ID (retentionTime set), not best
        assert_eq!(rows[0].4, 0);
        assert!((rows[0].1.unwrap() - 10.5).abs() < 0.01);
        assert!((rows[0].2 - 9.5).abs() < 0.01);
        // Run 2: has ID, best
        assert_eq!(rows[1].4, 1);
        assert!((rows[1].1.unwrap() - 10.3).abs() < 0.01);
        // Run 3: no ID (retentionTime NULL), but has boundaries
        assert_eq!(rows[2].4, 0);
        assert!(rows[2].1.is_none(), "Run 3 should have NULL retentionTime");
        assert!((rows[2].2 - 9.7).abs() < 0.01);
        assert!((rows[2].3 - 11.7).abs() < 0.01);
    }
}
