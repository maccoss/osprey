# Osprey Utility Scripts

Python scripts for evaluating and visualizing Osprey results.

## Setup

```bash
pip install -r scripts/requirements.txt
```

## Scripts

### evaluate_calibration.py

Generates an interactive HTML report from Osprey calibration JSON files. Useful for inspecting calibration quality before trusting main search results.

**Report contents:**

- MS1 and MS2 mass accuracy histograms
- RT calibration curve and RT shift plot
- Per-file calibration metrics
- Candidate density heatmap (requires `--library`)

```bash
# Single file
python scripts/evaluate_calibration.py sample.calibration.json

# Multiple files (overlay)
python scripts/evaluate_calibration.py *.calibration.json --output report.html

# With candidate density heatmap
python scripts/evaluate_calibration.py sample.calibration.json --library library.tsv

# Specify isolation width (auto-detected if available in JSON)
python scripts/evaluate_calibration.py sample.calibration.json --library library.tsv --isolation-width 8
```

| Argument | Description |
|----------|-------------|
| `calibration_files` | One or more `*.calibration.json` files |
| `--library` | Spectral library in DIA-NN TSV format (for density heatmap) |
| `--isolation-width` | DIA isolation window width in Da (auto-detected if omitted) |
| `--output`, `-o` | Output HTML file (default: `calibration_report.html`) |

---

### visualize_pin_features.py

Generates an HTML report comparing target vs decoy distributions for each feature in a PIN file. Features are sorted by separation score (Cohen's d), making it easy to identify which features have discriminative power.

**Separation score:**  `|mean_targets - mean_decoys| / sqrt((sd_targets² + sd_decoys²) / 2` — Cohen's d (effect size) — the absolute difference in means divided by the pooled standard deviation

```bash
# Default (3 columns, 50 bins)
python scripts/visualize_pin_features.py input.pin

# Custom layout
python scripts/visualize_pin_features.py input.pin --output report.html --cols-per-row 4 --bins 30
```

| Argument | Description |
|----------|-------------|
| `pin_file` | Input PIN file (Percolator/Mokapot format) |
| `--output`, `-o` | Output HTML file (default: `<input>.features.html`) |
| `--cols-per-row`, `-c` | Number of plots per row (default: 3) |
| `--bins`, `-b` | Number of histogram bins (default: 50) |

---

### inspect_mokapot_weights.py

Displays feature weights from a trained Mokapot model, ranked by absolute importance. Requires Osprey to have been run with `--save_models` to produce the model pickle file.

```bash
# Show all feature weights
python scripts/inspect_mokapot_weights.py mokapot/mokapot.model.pkl

# Show top 10 only
python scripts/inspect_mokapot_weights.py mokapot/mokapot.model.pkl --top 10

# Save to TSV
python scripts/inspect_mokapot_weights.py mokapot/mokapot.model.pkl --output weights.tsv
```

| Argument | Description |
|----------|-------------|
| `model_file` | Mokapot model pickle file (`mokapot.model.pkl`) |
| `--top` | Show only top N features by importance |
| `--output` | Save weights to TSV file |

---

### build_fdrbench_input.py

Builds an [FDRBench](https://github.com/Noble-Lab/FDRBench)-compatible input TSV from either a DIA-NN `report.parquet` or an Osprey `.blib` (+ companion `proteins.csv`). Used when you want to evaluate FDR-control quality via entrapment counting *after the fact*, against results you don't have the raw scoring artefacts for.

**Note:** Osprey can emit this format natively via the `--fdrbench` CLI flag — that path is preferred when you control the Osprey run, because it sees every scored target (not just FDR-passing ones) and uses the raw SVM discriminant as `score`. This script is the fallback when you only have the finished blib / report.

```bash
# DIA-NN report -> FDRBench precursor-level input
python scripts/build_fdrbench_input.py diann \
    --input report.parquet --level precursor --output diann_precursor.tsv

# Osprey blib -> FDRBench protein-level input (reads sibling proteins.csv)
python scripts/build_fdrbench_input.py osprey \
    --input osprey_output.blib --level protein --output osprey_protein.tsv

# Add --per-run for one row per (precursor, run) using run-level q-values
python scripts/build_fdrbench_input.py osprey \
    --input osprey_output.blib --level precursor --per-run \
    --output osprey_precursor_per_run.tsv
```

| Argument | Description |
|----------|-------------|
| `{diann,osprey}` | Source format |
| `--input`, `-i` | DIA-NN `report.parquet` or Osprey `.blib` |
| `--output`, `-o` | Output TSV path |
| `--level` | `peptide`, `precursor`, or `protein` |
| `--per-run` | One row per (precursor, run); default is best-across-runs |
| `--proteins-csv` | Override path to Osprey's `proteins.csv` (default: sibling of the blib) |

Invoke FDRBench with `-score 'score:1'` (higher is better — `score` is the raw upstream discriminant: DIA-NN `Evidence`, Osprey `DiscriminantScore` with `1 - q_value` fallback when the blib doesn't carry it, or `best_peptide_score` for protein groups).

---

### build_entrapment_peptide_fasta.py

Generates a peptide-level FASTA + FDRBench-style pairing manifest from a target-only uniprot FASTA. Use as the input to a library predictor (Carafe, etc.) to produce a balanced target / entrapment / decoy / p_decoy library.

For each tryptic peptide:

- **target** is the peptide as digested from the source protein.
- **`p_target`** (entrapment, optional): deterministic shuffle of the target with the C-terminal residue preserved, using a per-peptide RNG seed derived as `SHA-1(entrapment_seed:sequence)`.
- **decoy** (optional): same shuffle rule but with a different master seed.
- **`p_decoy`** (optional, when entrapment also enabled): shuffle of the entrapment peptide using the decoy master seed.

Any synthetic peptide that happens to collide with a real target somewhere in the input drops the whole quartet (matches FDRBench's policy).

Each FASTA entry carries the **source protein accession** (no per-peptide ID in the accession), so multiple peptides from one protein share the same `ProteinID` in the downstream library and Osprey's protein parsimony / picked-protein FDR work as on any normal library.

```bash
# Targets + entrapment + decoys with all defaults (len 7-35, m/z 400-900,
# charges 2 and 3, 1 missed cleavage, decoy_ / _p_target conventions)
python scripts/build_entrapment_peptide_fasta.py \
    --input uniprot_human.fasta \
    --output entrapment_peptides.fasta \
    --manifest pairing_manifest.tsv \
    --add-entrapment --add-decoys

# Tighter precursor m/z window for an Astral run, charge 2-4
python scripts/build_entrapment_peptide_fasta.py \
    --input uniprot_human.fasta --output peptides.fasta --manifest manifest.tsv \
    --add-entrapment --add-decoys \
    --min-mz 400 --max-mz 900 --charges 2 3 4 \
    --min-length 7 --max-length 35 --missed-cleavages 2
```

| Argument | Description |
|----------|-------------|
| `--input` | Source uniprot FASTA (targets only) |
| `--output` | Output peptide-level FASTA |
| `--manifest` | Output FDRBench-style pairing manifest TSV (optional but recommended) |
| `--add-entrapment` | Include `p_target` entrapment peptides |
| `--add-decoys` | Include decoy peptides (and `p_decoys` when `--add-entrapment` is also set) |
| `--missed-cleavages` | Max missed cleavages [default: 1; allowed: 0..5]. A value of N includes peptides with 0..N missed cleavages. |
| `--min-length`, `--max-length` | Peptide length window [default: 7-35] |
| `--min-mz`, `--max-mz` | Precursor m/z window; a peptide is kept iff at least one allowed charge state produces an m/z inside this window [default: 400-900] |
| `--charges` | Allowed precursor charge states [default: 2 3] |
| `--entrapment-seed` | Master RNG seed for entrapment shuffling [default: 42] |
| `--decoy-seed` | Master RNG seed for decoy shuffling [default: 24] |
| `--entrapment-suffix` | Suffix appended to accession/entry-name for entrapment proteins [default: `_p_target`] |
| `--decoy-prefix` | Prefix prepended to db\|accession\|entry for decoy headers [default: `decoy_`] |

The manifest output uses FDRBench's 5-column format (`sequence`, `decoy`, `proteins`, `peptide_type`, `peptide_pair_index`) and can be passed directly to Osprey via `--decoy-pairing-manifest`.

---

### fix_library_decoy_column.py

Stream-processes a DIA-NN-style library TSV to set the `Decoy` column to `1` on rows whose `ProteinID` contains a configured decoy-prefix accession (`decoy_` / `rev_` / `DECOY_` by default, case-insensitive).

Useful when Carafe generates a library that reliably sets the protein-accession prefix but leaves the `Decoy` column at `0` for all rows, which blocks DIA-NN from recognising the decoys.

**Not required for Osprey workflows**: Osprey's DIA-NN TSV loader reads BOTH the `Decoy` column and the protein-prefix scan, so a library with only the prefix set will work directly. Kept around for DIA-NN-only consumers of the library.

```bash
# Use default decoy prefixes (decoy_, rev_, DECOY_)
python scripts/fix_library_decoy_column.py \
    --input  carafe_library.tsv \
    --output carafe_library_fixed.tsv

# Custom prefix list
python scripts/fix_library_decoy_column.py \
    --input lib.tsv --output lib_fixed.tsv \
    --decoy-prefixes decoy_ rev_
```

| Argument | Description |
|----------|-------------|
| `--input` | Input library TSV |
| `--output` | Output library TSV (written via temp file + atomic rename) |
| `--decoy-prefixes` | Protein-accession prefixes that mark decoys, case-insensitive [default: `decoy_ rev_ DECOY_`] |
| `--protein-column` | Name of the protein column [default: `ProteinID`] |
| `--decoy-column` | Name of the Decoy column [default: `Decoy`] |
