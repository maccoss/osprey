#!/usr/bin/env python3
"""
Feature Distribution Visualization for Osprey

Reads an Osprey report parquet file or Percolator/Mokapot PIN file and generates
histograms comparing target vs decoy distributions for each feature. Features are
ranked by ROC AUC to identify discriminative power.

Usage:
  python scripts/visualize_pin_features.py report.parquet
  python scripts/visualize_pin_features.py input.pin
  python scripts/visualize_pin_features.py report.parquet --output report.html
  python scripts/visualize_pin_features.py input.pin --cols-per-row 4

Dependencies:
  pip install plotly numpy pandas
  pip install pyarrow  # only needed for parquet input
"""

import argparse
import sys
from pathlib import Path

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.io as pio
except ImportError:
    print("Error: plotly is required. Install with: pip install plotly", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Error: numpy is required. Install with: pip install numpy", file=sys.stderr)
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("Error: pandas is required. Install with: pip install pandas", file=sys.stderr)
    sys.exit(1)


def load_pin_file(path):
    """Load a PIN file into a pandas DataFrame.

    Returns (df, is_target_mask) where is_target_mask is a boolean Series.
    """
    df = pd.read_csv(path, sep='\t')
    is_target = df['Label'] == 1
    return df, is_target


def load_parquet_file(path):
    """Load an Osprey report parquet file into a pandas DataFrame.

    Returns (df, is_target_mask) where is_target_mask is a boolean Series.
    """
    try:
        import pyarrow  # noqa: F401
    except ImportError:
        print("Error: pyarrow is required for parquet files. Install with: pip install pyarrow",
              file=sys.stderr)
        sys.exit(1)

    df = pd.read_parquet(path)
    is_target = ~df['Is.Decoy']
    return df, is_target


def get_feature_columns(df, is_parquet=False):
    """Get the list of numeric feature columns (excluding metadata)."""
    # Metadata columns to exclude
    pin_metadata = {'SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'ChargeState'}
    parquet_metadata = {
        'Run', 'Modified.Sequence', 'Stripped.Sequence', 'Precursor.Charge',
        'Precursor.Mz', 'Protein.Ids', 'Is.Decoy', 'RT', 'RT.Start', 'RT.Stop',
        'Peak.Width', 'Library.RT', 'Scan.Number', 'Score', 'Q.Value',
        'Global.Q.Value', 'PEP', 'Search.Mode',
    }
    metadata_cols = pin_metadata | parquet_metadata

    feature_cols = []
    for col in df.columns:
        if col in metadata_cols:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            feature_cols.append(col)

    return feature_cols


def compute_roc_auc(targets, decoys):
    """Compute ROC AUC using the Mann-Whitney U statistic.

    AUC > 0.5 means targets tend to have higher values.
    AUC < 0.5 means decoys tend to have higher values.
    AUC = 0.5 means no discrimination.
    """
    n_t = len(targets)
    n_d = len(decoys)
    if n_t == 0 or n_d == 0:
        return 0.5

    combined = np.concatenate([targets, decoys])
    is_target = np.concatenate([np.ones(n_t, dtype=bool), np.zeros(n_d, dtype=bool)])

    # Sort by value (mergesort for stable ordering of ties)
    order = np.argsort(combined, kind='mergesort')
    sorted_vals = combined[order]
    is_target_sorted = is_target[order]

    # Compute ranks with tie handling (average method)
    n = len(sorted_vals)
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i
        while j < n and sorted_vals[j] == sorted_vals[i]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0  # 1-indexed average
        ranks[i:j] = avg_rank
        i = j

    # Mann-Whitney U: sum of target ranks minus minimum possible
    target_rank_sum = ranks[is_target_sorted].sum()
    u = target_rank_sum - n_t * (n_t + 1) / 2.0
    auc = u / (n_t * n_d)
    return auc


def generate_report(df, is_target, output_path, input_path, cols_per_row=4, n_bins=50):
    """Generate an HTML report with feature histograms."""

    targets = df[is_target]
    decoys = df[~is_target]

    n_targets = len(targets)
    n_decoys = len(decoys)

    print(f"Loaded {len(df)} entries: {n_targets} targets, {n_decoys} decoys")

    # Get feature columns
    is_parquet = input_path.suffix.lower() == '.parquet'
    feature_cols = get_feature_columns(df, is_parquet=is_parquet)
    n_features = len(feature_cols)

    print(f"Found {n_features} feature columns")

    # Compute ROC AUC for each feature
    auc_scores = {}
    for col in feature_cols:
        t_vals = targets[col].dropna().values
        d_vals = decoys[col].dropna().values
        auc_scores[col] = compute_roc_auc(t_vals, d_vals)

    # Sort features by discriminative power: |AUC - 0.5| descending
    sorted_features = sorted(feature_cols, key=lambda x: abs(auc_scores[x] - 0.5), reverse=True)

    # Calculate grid dimensions
    n_rows = (n_features + cols_per_row - 1) // cols_per_row

    # Create subplot titles with AUC scores
    subplot_titles = []
    for col in sorted_features:
        auc = auc_scores[col]
        direction = "T>D" if auc >= 0.5 else "D>T"
        subplot_titles.append(f"{col}<br><sup>AUC={auc:.3f} ({direction})</sup>")

    # Pad with empty strings if needed
    while len(subplot_titles) < n_rows * cols_per_row:
        subplot_titles.append("")

    # Create figure with subplots
    v_spacing = max(0.02, 0.4 / n_rows) if n_rows > 1 else 0.1
    fig = make_subplots(
        rows=n_rows,
        cols=cols_per_row,
        subplot_titles=subplot_titles,
        vertical_spacing=v_spacing,
        horizontal_spacing=0.03,
    )

    # Add histograms for each feature
    for idx, col in enumerate(sorted_features):
        row = idx // cols_per_row + 1
        col_idx = idx % cols_per_row + 1

        t_vals = targets[col].dropna().values
        d_vals = decoys[col].dropna().values

        if len(t_vals) == 0 and len(d_vals) == 0:
            continue

        # Compute common bin edges
        all_vals = np.concatenate([t_vals, d_vals])

        val_range = np.nanmax(all_vals) - np.nanmin(all_vals)
        if val_range < 1e-10:
            bin_edges = np.array([np.nanmin(all_vals) - 0.5, np.nanmax(all_vals) + 0.5])
        else:
            p1, p99 = np.nanpercentile(all_vals, [1, 99])
            if p99 - p1 < 1e-10:
                bin_edges = np.linspace(np.nanmin(all_vals), np.nanmax(all_vals), n_bins + 1)
            else:
                bin_edges = np.linspace(p1, p99, n_bins + 1)

        # Compute histograms and normalize to density
        t_counts, _ = np.histogram(t_vals, bins=bin_edges)
        d_counts, _ = np.histogram(d_vals, bins=bin_edges)

        t_density = t_counts / (len(t_vals) + 1e-10) if len(t_vals) > 0 else t_counts
        d_density = d_counts / (len(d_vals) + 1e-10) if len(d_vals) > 0 else d_counts

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        fig.add_trace(
            go.Bar(
                x=bin_centers,
                y=t_density,
                name="Targets",
                marker_color="rgba(55, 128, 191, 0.6)",
                showlegend=(idx == 0),
                legendgroup="targets",
            ),
            row=row, col=col_idx,
        )

        fig.add_trace(
            go.Bar(
                x=bin_centers,
                y=d_density,
                name="Decoys",
                marker_color="rgba(219, 64, 82, 0.6)",
                showlegend=(idx == 0),
                legendgroup="decoys",
            ),
            row=row, col=col_idx,
        )

    # Update layout
    fig.update_layout(
        width=1800,
        height=440 * n_rows,
        title_text=(
            f"Feature Distributions: {n_targets:,} targets vs {n_decoys:,} decoys"
            f"<br><sup>Sorted by ROC AUC (further from 0.5 = better discrimination). "
            f"T>D = targets higher, D>T = decoys higher.</sup>"
        ),
        barmode='overlay',
        template="plotly_white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        ),
        margin=dict(l=40, r=20, t=100, b=40),
    )

    fig.update_xaxes(tickfont=dict(size=9))
    fig.update_yaxes(tickfont=dict(size=9), title_text="")

    # Generate HTML
    plot_html = pio.to_html(fig, full_html=False, include_plotlyjs=True)

    # Create summary table
    summary_rows = []
    for col in sorted_features[:20]:
        t_vals = targets[col].dropna()
        d_vals = decoys[col].dropna()
        auc = auc_scores[col]
        direction = "T>D" if auc >= 0.5 else "D>T"

        summary_rows.append(f"""
        <tr>
            <td>{col}</td>
            <td>{auc:.3f}</td>
            <td>{direction}</td>
            <td>{t_vals.mean():.4g}</td>
            <td>{d_vals.mean():.4g}</td>
            <td>{t_vals.std():.4g}</td>
            <td>{d_vals.std():.4g}</td>
        </tr>
        """)

    summary_table = f"""
    <h2>Top 20 Features by ROC AUC</h2>
    <table class="summary-table">
        <thead>
            <tr>
                <th>Feature</th>
                <th>ROC AUC</th>
                <th>Direction</th>
                <th>Target Mean</th>
                <th>Decoy Mean</th>
                <th>Target SD</th>
                <th>Decoy SD</th>
            </tr>
        </thead>
        <tbody>
            {"".join(summary_rows)}
        </tbody>
    </table>
    """

    input_name = input_path.name

    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Feature Analysis — {input_name}</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    margin: 0;
    padding: 20px;
    background: #f5f5f5;
    color: #333;
  }}
  h1 {{
    margin-bottom: 10px;
    color: #1a1a2e;
  }}
  h2 {{
    color: #1a1a2e;
    margin-top: 30px;
  }}
  .subtitle {{
    color: #666;
    margin-bottom: 20px;
  }}
  .summary-table {{
    border-collapse: collapse;
    margin: 20px 0;
    font-size: 14px;
    min-width: 600px;
    background: white;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
  }}
  .summary-table th, .summary-table td {{
    padding: 10px 15px;
    text-align: left;
    border-bottom: 1px solid #ddd;
  }}
  .summary-table th {{
    background: #f8f9fa;
    font-weight: 600;
  }}
  .summary-table tr:hover {{
    background: #f5f5f5;
  }}
  .plot-container {{
    background: white;
    padding: 20px;
    border-radius: 8px;
    box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    margin-top: 20px;
  }}
  .info-box {{
    background: #e7f3ff;
    border-left: 4px solid #2196F3;
    padding: 15px;
    margin: 20px 0;
  }}
</style>
</head>
<body>
<h1>Feature Analysis</h1>
<p class="subtitle">Feature distributions for target vs decoy entries</p>

<div class="info-box">
    <strong>Summary:</strong> {n_targets:,} targets, {n_decoys:,} decoys, {n_features} features<br>
    <strong>Input:</strong> {input_name}<br>
    <strong>Metric:</strong> ROC AUC (0.5 = no discrimination, 1.0 or 0.0 = perfect separation)<br>
    <strong>Interpretation:</strong> Features further from AUC=0.5 have better discriminative power.
    T>D means targets have higher values; D>T means decoys have higher values.
    Blue = targets, Red = decoys.
</div>

{summary_table}

<div class="plot-container">
{plot_html}
</div>

</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(full_html)

    print(f"\nReport written to: {output_path}")
    print(f"\nTop 5 features by ROC AUC:")
    for i, col in enumerate(sorted_features[:5], 1):
        auc = auc_scores[col]
        direction = "T>D" if auc >= 0.5 else "D>T"
        print(f"  {i}. {col}: AUC={auc:.3f} ({direction})")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize feature distributions from Osprey parquet or PIN files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s report.parquet
  %(prog)s input.pin
  %(prog)s report.parquet --output feature_report.html
  %(prog)s input.pin --cols-per-row 4 --bins 30
        """,
    )
    parser.add_argument(
        "input_file",
        help="Input file: Osprey report parquet or Percolator/Mokapot PIN",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output HTML file (default: <input>.features.html)",
    )
    parser.add_argument(
        "--cols-per-row", "-c",
        type=int,
        default=4,
        help="Number of plots per row (default: 4)",
    )
    parser.add_argument(
        "--bins", "-b",
        type=int,
        default=50,
        help="Number of histogram bins (default: 50)",
    )

    args = parser.parse_args()

    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Set default output path
    if args.output is None:
        output_path = input_path.with_suffix('.features.html')
    else:
        output_path = Path(args.output)

    # Load based on file extension
    ext = input_path.suffix.lower()
    if ext == '.parquet':
        print(f"Loading parquet file: {input_path}")
        df, is_target = load_parquet_file(input_path)
    elif ext == '.pin':
        print(f"Loading PIN file: {input_path}")
        df, is_target = load_pin_file(input_path)
    else:
        print(f"Error: Unsupported file format '{ext}'. Use .parquet or .pin", file=sys.stderr)
        sys.exit(1)

    generate_report(df, is_target, str(output_path), input_path,
                    cols_per_row=args.cols_per_row, n_bins=args.bins)


if __name__ == "__main__":
    main()
