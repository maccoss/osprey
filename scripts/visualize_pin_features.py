#!/usr/bin/env python3
"""
PIN File Feature Visualization

Reads a Percolator/Mokapot PIN file and generates histograms comparing
target vs decoy distributions for each feature. Helps identify which
features have discriminative power.

Usage:
  python scripts/visualize_pin_features.py input.pin
  python scripts/visualize_pin_features.py input.pin --output report.html
  python scripts/visualize_pin_features.py input.pin --cols-per-row 4

Dependencies:
  pip install plotly numpy pandas
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
    """Load a PIN file into a pandas DataFrame."""
    df = pd.read_csv(path, sep='\t')
    return df


def get_feature_columns(df):
    """Get the list of numeric feature columns (excluding metadata)."""
    # Standard PIN metadata columns to exclude
    metadata_cols = {'SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'ChargeState'}

    feature_cols = []
    for col in df.columns:
        if col in metadata_cols:
            continue
        # Check if column is numeric
        if pd.api.types.is_numeric_dtype(df[col]):
            feature_cols.append(col)

    return feature_cols


def compute_separation_score(targets, decoys):
    """Compute a simple separation score between target and decoy distributions.

    Returns the absolute difference in means divided by pooled standard deviation.
    Higher values indicate better separation.
    """
    if len(targets) == 0 or len(decoys) == 0:
        return 0.0

    mean_t = np.nanmean(targets)
    mean_d = np.nanmean(decoys)
    std_t = np.nanstd(targets)
    std_d = np.nanstd(decoys)

    # Pooled standard deviation
    pooled_std = np.sqrt((std_t**2 + std_d**2) / 2)

    if pooled_std < 1e-10:
        return 0.0

    return abs(mean_t - mean_d) / pooled_std


def generate_report(df, output_path, cols_per_row=4, n_bins=50):
    """Generate an HTML report with feature histograms."""

    # Separate targets and decoys
    targets = df[df['Label'] == 1]
    decoys = df[df['Label'] == -1]

    n_targets = len(targets)
    n_decoys = len(decoys)

    print(f"Loaded {len(df)} PSMs: {n_targets} targets, {n_decoys} decoys")

    # Get feature columns
    feature_cols = get_feature_columns(df)
    n_features = len(feature_cols)

    print(f"Found {n_features} feature columns")

    # Compute separation scores for each feature
    separation_scores = {}
    for col in feature_cols:
        t_vals = targets[col].dropna().values
        d_vals = decoys[col].dropna().values
        separation_scores[col] = compute_separation_score(t_vals, d_vals)

    # Sort features by separation score (descending)
    sorted_features = sorted(feature_cols, key=lambda x: separation_scores[x], reverse=True)

    # Calculate grid dimensions
    n_rows = (n_features + cols_per_row - 1) // cols_per_row

    # Create subplot titles with separation scores
    subplot_titles = []
    for col in sorted_features:
        score = separation_scores[col]
        subplot_titles.append(f"{col}<br><sup>sep={score:.2f}</sup>")

    # Pad with empty strings if needed
    while len(subplot_titles) < n_rows * cols_per_row:
        subplot_titles.append("")

    # Create figure with subplots
    # Scale spacing inversely with row count to avoid whitespace dominating
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

        # Skip if no data
        if len(t_vals) == 0 and len(d_vals) == 0:
            continue

        # Compute common bin edges
        all_vals = np.concatenate([t_vals, d_vals])

        # Handle constant values
        val_range = np.nanmax(all_vals) - np.nanmin(all_vals)
        if val_range < 1e-10:
            # All values are the same - use a single bin
            bin_edges = np.array([np.nanmin(all_vals) - 0.5, np.nanmax(all_vals) + 0.5])
        else:
            # Use percentiles to handle outliers
            p1, p99 = np.nanpercentile(all_vals, [1, 99])
            if p99 - p1 < 1e-10:
                # Percentiles are the same, use full range
                bin_edges = np.linspace(np.nanmin(all_vals), np.nanmax(all_vals), n_bins + 1)
            else:
                bin_edges = np.linspace(p1, p99, n_bins + 1)

        # Compute histograms
        t_counts, _ = np.histogram(t_vals, bins=bin_edges)
        d_counts, _ = np.histogram(d_vals, bins=bin_edges)

        # Normalize to density
        t_density = t_counts / (len(t_vals) + 1e-10) if len(t_vals) > 0 else t_counts
        d_density = d_counts / (len(d_vals) + 1e-10) if len(d_vals) > 0 else d_counts

        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Add target histogram
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

        # Add decoy histogram
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
        title_text=f"PIN Feature Distributions: {n_targets:,} targets vs {n_decoys:,} decoys<br><sup>Sorted by separation score (higher = better discrimination)</sup>",
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

    # Update all x-axes to have smaller font
    fig.update_xaxes(tickfont=dict(size=9))
    fig.update_yaxes(tickfont=dict(size=9), title_text="")

    # Generate HTML
    plot_html = pio.to_html(fig, full_html=False, include_plotlyjs=True)

    # Create summary table
    summary_rows = []
    for col in sorted_features[:20]:  # Top 20 features
        t_vals = targets[col].dropna()
        d_vals = decoys[col].dropna()

        summary_rows.append(f"""
        <tr>
            <td>{col}</td>
            <td>{separation_scores[col]:.3f}</td>
            <td>{t_vals.mean():.4g}</td>
            <td>{d_vals.mean():.4g}</td>
            <td>{t_vals.std():.4g}</td>
            <td>{d_vals.std():.4g}</td>
        </tr>
        """)

    summary_table = f"""
    <h2>Top 20 Features by Separation Score</h2>
    <table class="summary-table">
        <thead>
            <tr>
                <th>Feature</th>
                <th>Separation</th>
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

    # Build full HTML
    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>PIN Feature Analysis</title>
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
<h1>PIN Feature Analysis</h1>
<p class="subtitle">Feature distributions for target vs decoy PSMs</p>

<div class="info-box">
    <strong>Summary:</strong> {n_targets:,} targets, {n_decoys:,} decoys, {n_features} features<br>
    <strong>Input:</strong> {Path(output_path).stem}<br>
    <strong>Interpretation:</strong> Features with higher separation scores show better discrimination
    between targets and decoys. Blue = targets, Red = decoys.
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
    print(f"\nTop 5 features by separation score:")
    for i, col in enumerate(sorted_features[:5], 1):
        print(f"  {i}. {col}: {separation_scores[col]:.3f}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize PIN file feature distributions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.pin
  %(prog)s input.pin --output feature_report.html
  %(prog)s input.pin --cols-per-row 4 --bins 30
        """,
    )
    parser.add_argument(
        "pin_file",
        help="Input PIN file (Percolator/Mokapot format)",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output HTML file (default: <input>_features.html)",
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

    # Validate input
    pin_path = Path(args.pin_file)
    if not pin_path.exists():
        print(f"Error: PIN file not found: {pin_path}", file=sys.stderr)
        sys.exit(1)

    # Set default output path
    if args.output is None:
        output_path = pin_path.with_suffix('.features.html')
    else:
        output_path = Path(args.output)

    # Load and process
    print(f"Loading PIN file: {pin_path}")
    df = load_pin_file(pin_path)

    generate_report(df, str(output_path), cols_per_row=args.cols_per_row, n_bins=args.bins)


if __name__ == "__main__":
    main()
