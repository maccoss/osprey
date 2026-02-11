#!/usr/bin/env python3
"""
Osprey Calibration Evaluation Report

Reads one or more Osprey calibration.json files and generates an interactive
HTML report with:
  - MS1 and MS2 mass accuracy histograms (side by side)
  - RT calibration curve and RT shift plot (side by side)
  - Per-file calibration metrics
  - Candidate density heatmap (requires --library)

Usage:
  python scripts/evaluate_calibration.py cal1.json [cal2.json ...]
  python scripts/evaluate_calibration.py cal1.json --library library.tsv --output report.html
  python scripts/evaluate_calibration.py cal1.json --library library.tsv --isolation-width 8

Dependencies:
  pip install plotly numpy
"""

import argparse
import json
import csv
import sys
import os
from pathlib import Path
import multiprocessing  # For cpu_count() if needed in future

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


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_calibration(path):
    """Load an Osprey calibration.json file."""
    with open(path) as f:
        data = json.load(f)
    data["_source_file"] = str(path)
    data["_label"] = Path(path).stem.replace(".calibration", "")
    return data


def load_diann_library(path):
    """Load a DIA-NN TSV spectral library, returning (precursor_mz, library_rt, charge) tuples.

    Only needs two columns: PrecursorMz and iRT (or RetentionTime / NormalizedRetentionTime).
    Returns unique precursors (deduplicated by modified_sequence + charge).
    """
    precursors = {}  # key: (modified_sequence, charge) -> (precursor_mz, rt)

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames
        if headers is None:
            raise ValueError(f"Empty library file: {path}")

        # Find the right column names (DIA-NN uses various naming)
        mz_col = None
        rt_col = None
        seq_col = None
        charge_col = None

        for h in headers:
            hl = h.strip().lower()
            if hl in ("precursormz", "precursor.mz", "q1"):
                mz_col = h
            elif hl in ("irt", "retentiontime", "normalizedretentiontime",
                         "rt", "retention_time", "tr_recalibrated"):
                rt_col = h
            elif hl in ("modifiedpeptide", "modified.sequence",
                         "modifiedpeptidesequence", "modified_sequence",
                         "fullunimodeprecessorname"):
                seq_col = h
            elif hl in ("precursorcharge", "precursor.charge", "charge"):
                charge_col = h

        if mz_col is None:
            raise ValueError(f"Cannot find PrecursorMz column in {path}. Headers: {headers}")
        if rt_col is None:
            raise ValueError(f"Cannot find RT column in {path}. Headers: {headers}")

        for row in reader:
            try:
                mz = float(row[mz_col])
                rt = float(row[rt_col])
            except (ValueError, KeyError):
                continue

            # Deduplicate by (sequence, charge) if columns available
            if seq_col and charge_col:
                key = (row.get(seq_col, ""), row.get(charge_col, ""))
            else:
                key = (mz, rt)  # fallback: keep all rows

            if key not in precursors:
                charge = int(row.get(charge_col, 0)) if charge_col else 0
                precursors[key] = (mz, rt, charge)

    return list(precursors.values())


# ---------------------------------------------------------------------------
# RT calibration helpers
# ---------------------------------------------------------------------------

def predict_rt(model_params, library_rt):
    """Replicate Osprey's RTCalibration::predict() in Python.

    Uses linear interpolation within the fitted curve and linear
    extrapolation at the edges.
    """
    import bisect

    lib_rts = model_params["library_rts"]
    fit_rts = model_params["fitted_rts"]
    n = len(lib_rts)
    if n == 0:
        return library_rt

    idx = bisect.bisect_left(lib_rts, library_rt)

    if idx == 0:
        if n < 2:
            return fit_rts[0]
        slope = (fit_rts[1] - fit_rts[0]) / (lib_rts[1] - lib_rts[0])
        return fit_rts[0] + slope * (library_rt - lib_rts[0])
    elif idx >= n:
        if n < 2:
            return fit_rts[-1]
        slope = (fit_rts[-1] - fit_rts[-2]) / (lib_rts[-1] - lib_rts[-2])
        return fit_rts[-1] + slope * (library_rt - lib_rts[-1])
    else:
        x0, x1 = lib_rts[idx - 1], lib_rts[idx]
        y0, y1 = fit_rts[idx - 1], fit_rts[idx]
        if abs(x1 - x0) < 1e-10:
            return y0
        t = (library_rt - x0) / (x1 - x0)
        return y0 + t * (y1 - y0)


# ---------------------------------------------------------------------------
# Subplot axis name helpers
# ---------------------------------------------------------------------------

def _axis_name(row, col, n_cols):
    """Return the plotly axis index for a given (row, col) in the subplot grid.

    Plotly numbers axes sequentially left-to-right, top-to-bottom:
      row=1,col=1 -> '' (first), row=1,col=2 -> '2', row=2,col=1 -> '3', etc.
    """
    idx = (row - 1) * n_cols + col
    return "" if idx == 1 else str(idx)


# ---------------------------------------------------------------------------
# Plot builders
# ---------------------------------------------------------------------------

def make_mass_accuracy_histogram(cal, level, field, color, fig, row, col, n_cols):
    """Add a single mass accuracy histogram to fig at the given row/col."""
    mz_cal = cal.get(field, {})
    hist = mz_cal.get("histogram")
    if hist is None:
        return

    edges = hist["bin_edges"]
    counts = hist["counts"]
    centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(counts))]

    # Get unit from calibration (ppm or Th)
    unit = mz_cal.get("unit", "ppm")
    is_th = unit == "Th"

    # Adjust format precision based on unit
    if is_th:
        hover_fmt = "%{x:.4f} Th: %{y} observations<extra></extra>"
        val_fmt = ".4f"
    else:
        hover_fmt = "%{x:.2f} ppm: %{y} observations<extra></extra>"
        val_fmt = ".3f"

    fig.add_trace(
        go.Bar(
            x=centers,
            y=counts,
            marker_color=color,
            name=f"{level} Mass Error",
            showlegend=False,
            hovertemplate=hover_fmt,
        ),
        row=row, col=col,
    )

    mean_val = mz_cal.get("mean", 0)
    median_val = mz_cal.get("median", 0)
    sd_val = mz_cal.get("sd", 0)
    count_val = mz_cal.get("count", 0)

    # Vertical lines for mean and median (no annotation text on the line)
    fig.add_vline(
        x=mean_val, line_dash="dash", line_color="red",
        row=row, col=col,
    )
    fig.add_vline(
        x=median_val, line_dash="dot", line_color="green",
        row=row, col=col,
    )

    # Stats annotation box in top-right corner of the subplot
    ax = _axis_name(row, col, n_cols)
    if is_th:
        stats_text = (
            f"<span style='color:red'>--- Mean: {mean_val:+.4f} {unit}</span><br>"
            f"<span style='color:green'>··· Median: {median_val:+.4f} {unit}</span><br>"
            f"SD: {sd_val:.4f} {unit}<br>"
            f"n = {count_val:,}"
        )
    else:
        stats_text = (
            f"<span style='color:red'>--- Mean: {mean_val:+.3f} {unit}</span><br>"
            f"<span style='color:green'>··· Median: {median_val:+.3f} {unit}</span><br>"
            f"SD: {sd_val:.3f} {unit}<br>"
            f"n = {count_val:,}"
        )
    fig.add_annotation(
        text=stats_text,
        xref=f"x{ax} domain", yref=f"y{ax} domain",
        x=0.03, y=0.95,
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=11),
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#ccc",
        borderwidth=1,
        align="left",
    )

    fig.update_xaxes(title_text=f"Mass Error ({unit})", row=row, col=col)
    fig.update_yaxes(title_text="Count", row=row, col=col)


def make_rt_alignment_plot(cal, fig, row, col, n_cols):
    """Add RT alignment plot: library RT vs fitted RT."""
    rt_cal = cal.get("rt_calibration", {})
    model = rt_cal.get("model_params")
    if model is None:
        return

    lib_rts = model["library_rts"]
    fit_rts = model["fitted_rts"]

    if len(lib_rts) == 0:
        return

    # Subsample for performance if very large
    n = len(lib_rts)
    if n > 5000:
        step = n // 5000
        indices = list(range(0, n, step))
        lib_rts_sub = [lib_rts[i] for i in indices]
        fit_rts_sub = [fit_rts[i] for i in indices]
    else:
        lib_rts_sub = lib_rts
        fit_rts_sub = fit_rts

    # Fitted curve
    fig.add_trace(
        go.Scattergl(
            x=lib_rts_sub,
            y=fit_rts_sub,
            mode="markers",
            marker=dict(size=2, color="#1f77b4", opacity=0.5),
            name="LOESS fit",
            showlegend=False,
            hovertemplate="Library RT: %{x:.2f}<br>Fitted RT: %{y:.2f}<extra></extra>",
        ),
        row=row, col=col,
    )

    # Identity line
    rt_min = min(lib_rts)
    rt_max = max(lib_rts)
    fig.add_trace(
        go.Scatter(
            x=[rt_min, rt_max],
            y=[rt_min, rt_max],
            mode="lines",
            line=dict(dash="dash", color="gray"),
            name="Identity",
            showlegend=False,
        ),
        row=row, col=col,
    )

    # Figures of merit
    r_squared = rt_cal.get("r_squared", 0)
    residual_sd = rt_cal.get("residual_sd", 0)
    n_points = rt_cal.get("n_points", n)
    rt_tolerance = max(residual_sd * 3.0, 0.5) if residual_sd > 0 else 0
    fit_min = min(fit_rts)
    fit_max = max(fit_rts)

    ax = _axis_name(row, col, n_cols)
    stats_text = (
        f"R\u00b2 = {r_squared:.6f}<br>"
        f"Residual SD = {residual_sd:.4f} min<br>"
        f"RT tolerance (3\u00d7SD) = {rt_tolerance:.3f} min<br>"
        f"Calibration points: {n_points:,}<br>"
        f"Library RT range: {rt_min:.1f} \u2013 {rt_max:.1f}<br>"
        f"Fitted RT range: {fit_min:.1f} \u2013 {fit_max:.1f}"
    )
    fig.add_annotation(
        text=stats_text,
        xref=f"x{ax} domain", yref=f"y{ax} domain",
        x=0.03, y=0.97,
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=11),
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#ccc",
        borderwidth=1,
        align="left",
    )

    fig.update_xaxes(title_text="Library RT (min)", row=row, col=col)
    fig.update_yaxes(title_text="Fitted Measured RT (min)", row=row, col=col)


def make_rt_residual_plot(cal, fig, row, col, n_cols):
    """Add RT shift plot (library RT vs fitted - library)."""
    rt_cal = cal.get("rt_calibration", {})
    model = rt_cal.get("model_params")
    if model is None:
        return

    lib_rts = model["library_rts"]
    fit_rts = model["fitted_rts"]
    residual_sd = rt_cal.get("residual_sd", 0)

    if len(lib_rts) == 0:
        return

    rt_shift = [f - l for f, l in zip(fit_rts, lib_rts)]

    # Subsample
    n = len(lib_rts)
    if n > 5000:
        step = n // 5000
        indices = list(range(0, n, step))
        lib_sub = [lib_rts[i] for i in indices]
        shift_sub = [rt_shift[i] for i in indices]
    else:
        lib_sub = lib_rts
        shift_sub = rt_shift

    fig.add_trace(
        go.Scattergl(
            x=lib_sub,
            y=shift_sub,
            mode="markers",
            marker=dict(size=2, color="#2ca02c", opacity=0.5),
            name="RT shift",
            showlegend=False,
            hovertemplate="Library RT: %{x:.2f}<br>Shift: %{y:.2f} min<extra></extra>",
        ),
        row=row, col=col,
    )

    # Horizontal line at 0
    fig.add_hline(y=0, line_dash="dash", line_color="gray", row=row, col=col)

    # Compute shift statistics
    shift_arr = np.array(rt_shift)
    shift_mean = float(np.mean(shift_arr))
    shift_median = float(np.median(shift_arr))
    shift_sd = float(np.std(shift_arr))
    shift_min = float(np.min(shift_arr))
    shift_max = float(np.max(shift_arr))
    shift_iqr = float(np.percentile(shift_arr, 75) - np.percentile(shift_arr, 25))

    ax = _axis_name(row, col, n_cols)
    stats_text = (
        f"Mean shift: {shift_mean:+.3f} min<br>"
        f"Median shift: {shift_median:+.3f} min<br>"
        f"SD of shift: {shift_sd:.3f} min<br>"
        f"IQR: {shift_iqr:.3f} min<br>"
        f"Range: [{shift_min:.2f}, {shift_max:+.2f}] min<br>"
        f"Residual SD: {residual_sd:.4f} min"
    )
    fig.add_annotation(
        text=stats_text,
        xref=f"x{ax} domain", yref=f"y{ax} domain",
        x=0.03, y=0.97,
        xanchor="left", yanchor="top",
        showarrow=False,
        font=dict(size=11),
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#ccc",
        borderwidth=1,
        align="left",
    )

    fig.update_xaxes(title_text="Library RT (min)", row=row, col=col)
    fig.update_yaxes(title_text="RT Shift (fitted \u2212 library, min)", row=row, col=col)


def _compute_density_for_mz_bin(args):
    """Worker function to compute candidate density for a single m/z bin.

    This is a module-level function so it can be pickled for multiprocessing.
    """
    mz_idx, mz_center, half_iso, rt_centers, rt_tolerance, prec_mzs, prec_cal_rts = args

    window_lower = mz_center - half_iso
    window_upper = mz_center + half_iso

    # Find precursors in this m/z window (vectorized)
    in_window = (prec_mzs >= window_lower) & (prec_mzs <= window_upper)
    window_rts = prec_cal_rts[in_window]

    # For each RT bin, count candidates
    row_counts = np.zeros(len(rt_centers), dtype=np.int32)
    for rt_idx, spectrum_rt in enumerate(rt_centers):
        row_counts[rt_idx] = np.sum(np.abs(window_rts - spectrum_rt) <= rt_tolerance)

    return mz_idx, row_counts


def _compute_density_vectorized(mz_centers, rt_centers, half_iso, rt_tolerance, prec_mzs, prec_cal_rts, sort_idx):
    """Compute density matrix using sorted m/z values and binary search.

    This is more efficient than the parallel approach for most cases because it
    avoids the multiprocessing overhead and uses O(log n) binary search.
    """
    n_mz_bins = len(mz_centers)
    n_rt_bins = len(rt_centers)
    density = np.zeros((n_mz_bins, n_rt_bins), dtype=np.int32)

    # Sort precursors by m/z for binary search
    sorted_mzs = prec_mzs[sort_idx]
    sorted_rts = prec_cal_rts[sort_idx]

    for mz_idx in range(n_mz_bins):
        mz_center = mz_centers[mz_idx]
        window_lower = mz_center - half_iso
        window_upper = mz_center + half_iso

        # Binary search to find precursors in m/z window
        left = np.searchsorted(sorted_mzs, window_lower, side='left')
        right = np.searchsorted(sorted_mzs, window_upper, side='right')

        if left >= right:
            continue  # No precursors in this window

        window_rts = sorted_rts[left:right]

        # Vectorized RT counting for all RT bins at once
        # Shape: (n_window_precursors, n_rt_bins)
        rt_diffs = np.abs(window_rts[:, np.newaxis] - rt_centers[np.newaxis, :])
        density[mz_idx, :] = np.sum(rt_diffs <= rt_tolerance, axis=0)

        # Progress update
        if (mz_idx + 1) % 10 == 0 or mz_idx == n_mz_bins - 1:
            pct = 100.0 * (mz_idx + 1) / n_mz_bins
            print(f"    Progress: {mz_idx + 1}/{n_mz_bins} m/z bins ({pct:.0f}%) - window {mz_center:.1f} m/z    ", end='\r')
            sys.stdout.flush()

    return density


def make_candidate_density_heatmap(cal, precursors, isolation_width, fig, row, col, n_cols):
    """Create heatmap: RT (x) vs precursor m/z (y), color = candidate count.

    For each grid cell, count how many library precursors would be candidates
    for a spectrum at that (RT, m/z center) given the calibrated RT tolerance
    and isolation window width.
    """
    rt_cal = cal.get("rt_calibration", {})
    model = rt_cal.get("model_params")
    residual_sd = rt_cal.get("residual_sd", 1.0)
    rt_tolerance = max(residual_sd * 3.0, 0.5)

    if model is None:
        return

    # Get m/z and RT ranges from library
    mzs = [p[0] for p in precursors]
    rts = [p[1] for p in precursors]

    if not mzs:
        return

    mz_min, mz_max = min(mzs), max(mzs)
    rt_min, rt_max = min(rts), max(rts)

    # Create grid
    n_rt_bins = 80
    half_iso = isolation_width / 2.0
    n_mz_bins = max(1, int((mz_max - mz_min) / isolation_width) + 1)

    rt_edges = np.linspace(rt_min, rt_max, n_rt_bins + 1)
    rt_centers = (rt_edges[:-1] + rt_edges[1:]) / 2

    mz_edges = np.linspace(mz_min - half_iso, mz_max + half_iso, n_mz_bins + 1)
    mz_centers = (mz_edges[:-1] + mz_edges[1:]) / 2

    # For each precursor, compute calibrated RT
    print(f"  Computing calibrated RTs for {len(precursors):,} precursors...")
    sys.stdout.flush()
    prec_mzs = np.array([p[0] for p in precursors])
    prec_cal_rts = np.array([predict_rt(model, p[1]) for p in precursors])

    # Sort precursors by m/z for efficient binary search
    sort_idx = np.argsort(prec_mzs)

    # Compute density matrix using vectorized approach with binary search
    mz_range_str = f"{mz_centers[0]:.1f}-{mz_centers[-1]:.1f} m/z"
    print(f"  Computing density for {n_mz_bins} m/z bins ({mz_range_str})...")
    sys.stdout.flush()

    density = _compute_density_vectorized(
        mz_centers, rt_centers, half_iso, rt_tolerance,
        prec_mzs, prec_cal_rts, sort_idx
    )

    print(f"    Completed: {n_mz_bins} m/z bins                                        ")

    # Compute summary stats
    median_candidates = int(np.median(density[density > 0])) if np.any(density > 0) else 0
    mean_candidates = float(np.mean(density[density > 0])) if np.any(density > 0) else 0
    max_candidates = int(np.max(density))

    # Apply sqrt scaling for better visualization (gentler than log2)
    density_scaled = np.sqrt(density)

    ax = _axis_name(row, col, n_cols)

    # Create custom hover text showing actual counts
    hover_text = [[f"RT: {rt_centers[j]:.1f} min<br>m/z center: {mz_centers[i]:.1f}<br>Candidates: {density[i, j]}"
                   for j in range(len(rt_centers))] for i in range(len(mz_centers))]

    # Use Turbo colorscale for better color differentiation across the range
    fig.add_trace(
        go.Heatmap(
            x=np.round(rt_centers, 2).tolist(),
            y=np.round(mz_centers, 1).tolist(),
            z=density_scaled.tolist(),
            colorscale="Turbo",
            colorbar=dict(
                title="√count",
                len=0.4,
                thickness=15,
                yanchor="bottom",
                y=0.0,
            ),
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
        ),
        row=row, col=col,
    )

    stats_text = (
        f"Isolation width: {isolation_width} Da<br>"
        f"RT tolerance: \u00b1{rt_tolerance:.2f} min<br>"
        f"Median candidates: {median_candidates}<br>"
        f"Mean candidates: {mean_candidates:.0f}<br>"
        f"Max candidates: {max_candidates}<br>"
        f"<i>(color: √ scale)</i>"
    )
    fig.add_annotation(
        text=stats_text,
        xref=f"x{ax} domain", yref=f"y{ax} domain",
        x=0.97, y=0.97,
        xanchor="right", yanchor="top",
        showarrow=False,
        font=dict(size=11),
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#ccc",
        borderwidth=1,
        align="left",
    )

    fig.update_xaxes(title_text="Retention Time (min)", row=row, col=col)
    fig.update_yaxes(title_text="Isolation Window Center (m/z)", row=row, col=col)


# ---------------------------------------------------------------------------
# Metrics summary
# ---------------------------------------------------------------------------

def format_metrics_table(cal):
    """Generate an HTML table of calibration metrics."""
    meta = cal.get("metadata", {})
    ms1 = cal.get("ms1_calibration", {})
    ms2 = cal.get("ms2_calibration", {})
    rt = cal.get("rt_calibration", {})
    scheme = meta.get("isolation_scheme")

    rt_tolerance = max(rt.get("residual_sd", 0) * 3.0, 0.5) if rt.get("residual_sd", 0) > 0 else "N/A"

    rows = [
        ("Calibration File", cal.get("_source_file", "unknown")),
        ("Timestamp", meta.get("timestamp", "N/A")),
        ("Calibration Successful", str(meta.get("calibration_successful", False))),
        ("Confident Peptides", f"{meta.get('num_confident_peptides', 0):,}"),
        ("Sampled Precursors", f"{meta.get('num_sampled_precursors', 0):,}"),
        ("", ""),
    ]

    # Add isolation scheme if available
    if scheme:
        rows.extend([
            ("<b>DIA Isolation Scheme</b>", ""),
            ("Number of Windows", f"{scheme.get('num_windows', 'N/A')}"),
            ("m/z Range", f"{scheme.get('mz_min', 0):.1f} - {scheme.get('mz_max', 0):.1f}"),
            ("Typical Width (Da)", f"{scheme.get('typical_width', 0):.1f}"),
            ("Uniform Width", "Yes" if scheme.get("uniform_width", False) else "No"),
            ("", ""),
        ])

    # MS1 can be in ppm (HRAM) or Th (unit resolution)
    ms1_unit = ms1.get("unit", "ppm")
    ms1_is_th = ms1_unit == "Th"
    if ms1_is_th:
        # Th units - use 4 decimal places
        rows.extend([
            ("<b>MS1 Mass Calibration (Unit Resolution)</b>", ""),
            (f"Mean Error ({ms1_unit})", f"{ms1.get('mean', 0):.4f}"),
            (f"Median Error ({ms1_unit})", f"{ms1.get('median', 0):.4f}"),
            (f"SD ({ms1_unit})", f"{ms1.get('sd', 0):.4f}"),
            ("Observations", f"{ms1.get('count', 0):,}"),
            (f"Adjusted Tolerance ({ms1_unit})", f"{ms1.get('adjusted_tolerance', 'N/A'):.4f}" if ms1.get('adjusted_tolerance') else "N/A"),
            ("", ""),
        ])
    else:
        # ppm units - use 3 decimal places
        rows.extend([
            ("<b>MS1 Mass Calibration</b>", ""),
            (f"Mean Error ({ms1_unit})", f"{ms1.get('mean', 0):.3f}"),
            (f"Median Error ({ms1_unit})", f"{ms1.get('median', 0):.3f}"),
            (f"SD ({ms1_unit})", f"{ms1.get('sd', 0):.3f}"),
            ("Observations", f"{ms1.get('count', 0):,}"),
            (f"Adjusted Tolerance ({ms1_unit})", f"{ms1.get('adjusted_tolerance', 'N/A')}"),
            ("", ""),
        ])

    # MS2 can be in ppm (HRAM) or Th (unit resolution)
    ms2_unit = ms2.get("unit", "ppm")
    is_th = ms2_unit == "Th"
    if is_th:
        # Th units - use 4 decimal places
        rows.extend([
            ("<b>MS2 Mass Calibration (Unit Resolution)</b>", ""),
            (f"Mean Error ({ms2_unit})", f"{ms2.get('mean', 0):.4f}"),
            (f"Median Error ({ms2_unit})", f"{ms2.get('median', 0):.4f}"),
            (f"SD ({ms2_unit})", f"{ms2.get('sd', 0):.4f}"),
            ("Observations", f"{ms2.get('count', 0):,}"),
            (f"Adjusted Tolerance ({ms2_unit})", f"{ms2.get('adjusted_tolerance', 'N/A'):.4f}" if ms2.get('adjusted_tolerance') else "N/A"),
            ("", ""),
        ])
    else:
        # ppm units - use 3 decimal places
        rows.extend([
            ("<b>MS2 Mass Calibration</b>", ""),
            (f"Mean Error ({ms2_unit})", f"{ms2.get('mean', 0):.3f}"),
            (f"Median Error ({ms2_unit})", f"{ms2.get('median', 0):.3f}"),
            (f"SD ({ms2_unit})", f"{ms2.get('sd', 0):.3f}"),
            ("Observations", f"{ms2.get('count', 0):,}"),
            (f"Adjusted Tolerance ({ms2_unit})", f"{ms2.get('adjusted_tolerance', 'N/A')}"),
            ("", ""),
        ])

    # RT Calibration
    rows.extend([
        ("<b>RT Calibration</b>", ""),
        ("Method", rt.get("method", "N/A")),
        ("Calibration Points", f"{rt.get('n_points', 0):,}"),
        ("R-squared", f"{rt.get('r_squared', 0):.6f}"),
        ("Residual SD (min)", f"{rt.get('residual_sd', 0):.4f}"),
        ("RT Tolerance (3x SD, min)", f"{rt_tolerance:.4f}" if isinstance(rt_tolerance, float) else rt_tolerance),
    ])

    html = '<table class="metrics-table">\n'
    for label, value in rows:
        if label == "" and value == "":
            html += '<tr><td colspan="2" style="height:8px"></td></tr>\n'
        elif value == "":
            html += f'<tr><td colspan="2">{label}</td></tr>\n'
        else:
            html += f"<tr><td>{label}</td><td>{value}</td></tr>\n"
    html += "</table>\n"
    return html


# ---------------------------------------------------------------------------
# HTML report generation
# ---------------------------------------------------------------------------

def generate_report(calibrations, precursors, isolation_width, output_path):
    """Generate the full HTML report."""
    n_files = len(calibrations)

    # Build per-file pages
    pages_html = []

    for file_idx, cal in enumerate(calibrations):
        label = cal["_label"]
        has_rt_model = cal.get("rt_calibration", {}).get("model_params") is not None
        has_histogram = (
            cal.get("ms1_calibration", {}).get("histogram") is not None
            or cal.get("ms2_calibration", {}).get("histogram") is not None
        )
        has_heatmap = precursors is not None and has_rt_model

        # Layout: 2-column grid
        #   Row 1: MS1 histogram | MS2 histogram
        #   Row 2: RT calibration | RT shift
        #   Row 3: Candidate density heatmap (spanning both columns)
        n_rows = 0
        subplot_titles = []
        specs = []
        row_heights = []

        if has_histogram:
            n_rows += 1
            subplot_titles.extend(["MS1 Mass Accuracy", "MS2 Mass Accuracy"])
            specs.append([{}, {}])
            row_heights.append(0.25)

        if has_rt_model:
            n_rows += 1
            subplot_titles.extend(["RT Calibration Curve", "RT Shift (Calibration Effect)"])
            specs.append([{}, {}])
            row_heights.append(0.30)

        if has_heatmap:
            n_rows += 1
            subplot_titles.extend(["Candidate Density per Spectrum", None])
            specs.append([{"colspan": 2}, None])
            row_heights.append(0.45)

        if n_rows == 0:
            n_rows = 1
            subplot_titles = ["No calibration data available", None]
            specs = [[{}, None]]
            row_heights = [1.0]

        # Normalize heights
        total = sum(row_heights)
        row_heights = [h / total for h in row_heights]

        fig = make_subplots(
            rows=n_rows, cols=2,
            subplot_titles=[t for t in subplot_titles],
            row_heights=row_heights,
            specs=specs,
            vertical_spacing=0.10,
            horizontal_spacing=0.08,
        )

        current_row = 1

        # Mass accuracy histograms (side by side)
        if has_histogram:
            make_mass_accuracy_histogram(
                cal, "MS1", "ms1_calibration", "#1f77b4",
                fig, current_row, 1, n_cols=2,
            )
            make_mass_accuracy_histogram(
                cal, "MS2", "ms2_calibration", "#ff7f0e",
                fig, current_row, 2, n_cols=2,
            )
            current_row += 1

        # RT plots (side by side)
        if has_rt_model:
            make_rt_alignment_plot(cal, fig, current_row, 1, n_cols=2)
            make_rt_residual_plot(cal, fig, current_row, 2, n_cols=2)
            current_row += 1

        # Candidate density heatmap (full width)
        if has_heatmap:
            make_candidate_density_heatmap(
                cal, precursors, isolation_width,
                fig, current_row, 1, n_cols=2,
            )
            current_row += 1

        fig.update_layout(
            height=max(500, 400 * n_rows),
            showlegend=False,
            title_text=f"Calibration Report: {label}",
            template="plotly_white",
            margin=dict(l=60, r=40, t=60, b=40),
        )

        plot_html = pio.to_html(fig, full_html=False, include_plotlyjs=(file_idx == 0))
        metrics_html = format_metrics_table(cal)

        page_html = f"""
        <div class="file-page" id="page-{file_idx}">
          <h2>{label}</h2>
          <div class="metrics-container">
            {metrics_html}
          </div>
          <div class="plots-container">
            {plot_html}
          </div>
        </div>
        """
        pages_html.append((label, page_html))

    # Build navigation tabs
    nav_html = '<div class="nav-tabs">\n'
    for i, (label, _) in enumerate(pages_html):
        active = "active" if i == 0 else ""
        nav_html += f'  <button class="tab-btn {active}" onclick="showPage({i})">{label}</button>\n'
    nav_html += "</div>\n"

    # Assemble full HTML
    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Osprey Calibration Report</title>
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
  .subtitle {{
    color: #666;
    margin-bottom: 20px;
    font-size: 14px;
  }}
  .nav-tabs {{
    display: flex;
    gap: 4px;
    margin-bottom: 20px;
    flex-wrap: wrap;
  }}
  .tab-btn {{
    padding: 8px 16px;
    border: 1px solid #ddd;
    background: #fff;
    cursor: pointer;
    border-radius: 4px 4px 0 0;
    font-size: 14px;
    transition: background 0.2s;
  }}
  .tab-btn:hover {{
    background: #e8e8e8;
  }}
  .tab-btn.active {{
    background: #1a1a2e;
    color: white;
    border-color: #1a1a2e;
  }}
  .file-page {{
    display: none;
    background: white;
    border-radius: 8px;
    padding: 20px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
  }}
  .file-page.visible {{
    display: block;
  }}
  .metrics-container {{
    margin-bottom: 20px;
  }}
  .metrics-table {{
    border-collapse: collapse;
    font-size: 14px;
    min-width: 400px;
  }}
  .metrics-table td {{
    padding: 4px 12px;
    border-bottom: 1px solid #eee;
  }}
  .metrics-table td:first-child {{
    font-weight: 500;
    color: #555;
  }}
  .metrics-table td:last-child {{
    font-family: monospace;
  }}
</style>
</head>
<body>
<h1>Osprey Calibration Report</h1>
<p class="subtitle">{n_files} file(s) analyzed</p>
{nav_html}
{"".join(html for _, html in pages_html)}
<script>
function showPage(idx) {{
  document.querySelectorAll('.file-page').forEach((el, i) => {{
    el.classList.toggle('visible', i === idx);
  }});
  document.querySelectorAll('.tab-btn').forEach((el, i) => {{
    el.classList.toggle('active', i === idx);
  }});
}}
// Show first page
showPage(0);
</script>
</body>
</html>"""

    with open(output_path, "w") as f:
        f.write(full_html)

    print(f"Report written to: {output_path}")
    print(f"  {n_files} calibration file(s) processed")
    for cal in calibrations:
        meta = cal.get("metadata", {})
        rt = cal.get("rt_calibration", {})
        print(f"  - {cal['_label']}: {meta.get('num_confident_peptides', 0)} peptides, "
              f"R\u00b2={rt.get('r_squared', 0):.4f}, "
              f"residual_sd={rt.get('residual_sd', 0):.3f} min")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate Osprey calibration evaluation report",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s sample.calibration.json
  %(prog)s *.calibration.json --output report.html
  %(prog)s sample.calibration.json --library library.tsv --isolation-width 8
        """,
    )
    parser.add_argument(
        "calibration_files",
        nargs="+",
        help="One or more Osprey calibration.json files",
    )
    parser.add_argument(
        "--library",
        help="Spectral library in DIA-NN TSV format (for candidate density heatmap)",
    )
    parser.add_argument(
        "--isolation-width",
        type=float,
        default=None,
        help="DIA isolation window width in Da (auto-detected from calibration if available)",
    )
    parser.add_argument(
        "--output", "-o",
        default="calibration_report.html",
        help="Output HTML file (default: calibration_report.html)",
    )

    args = parser.parse_args()

    # Load calibration files
    calibrations = []
    for path in args.calibration_files:
        if not os.path.exists(path):
            print(f"Warning: {path} not found, skipping", file=sys.stderr)
            continue
        try:
            cal = load_calibration(path)
            calibrations.append(cal)
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Warning: Failed to parse {path}: {e}", file=sys.stderr)
            continue

    if not calibrations:
        print("Error: No valid calibration files loaded", file=sys.stderr)
        sys.exit(1)

    # Load library if provided
    precursors = None
    if args.library:
        if not os.path.exists(args.library):
            print(f"Warning: Library file {args.library} not found, skipping heatmap",
                  file=sys.stderr)
        else:
            try:
                precursors = load_diann_library(args.library)
                print(f"Loaded {len(precursors)} unique precursors from library")
            except Exception as e:
                print(f"Warning: Failed to load library: {e}", file=sys.stderr)

    # Determine isolation width: use calibration file if available, else CLI arg, else default
    isolation_width = args.isolation_width
    if isolation_width is None:
        # Try to get from first calibration file
        for cal in calibrations:
            scheme = cal.get("metadata", {}).get("isolation_scheme")
            if scheme:
                isolation_width = scheme.get("typical_width", 8.0)
                print(f"Using isolation width from calibration: {isolation_width:.1f} Da "
                      f"({scheme.get('num_windows', '?')} windows, "
                      f"{scheme.get('mz_min', 0):.0f}-{scheme.get('mz_max', 0):.0f} m/z)")
                break
        if isolation_width is None:
            isolation_width = 8.0
            print(f"Warning: No isolation scheme in calibration file, using default {isolation_width} Da")

    generate_report(calibrations, precursors, isolation_width, args.output)


if __name__ == "__main__":
    # Required for Windows multiprocessing support
    multiprocessing.freeze_support()
    main()
