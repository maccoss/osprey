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
