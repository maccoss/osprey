#!/usr/bin/env python3
"""Compare retention times between Osprey and DIA-NN for matched peptides.

Reads an Osprey parquet report and a DIA-NN report parquet, matches peptides
by stripped sequence, and creates a density plot + residual plot comparing RTs.

Usage:
    python compare_scan_numbers.py osprey_report.parquet diann_report.parquet [--output plot.png]

Example:
    python scripts/compare_scan_numbers.py \
        example_test_data/astral/osprey_report.parquet \
        example_test_data/astral/carafe/diann_train/report.parquet
"""

import argparse
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Compare retention times: Osprey vs DIA-NN"
    )
    parser.add_argument("osprey_file", help="Osprey parquet report file")
    parser.add_argument("diann_file", help="DIA-NN report parquet file")
    parser.add_argument(
        "--output", "-o", default=None, help="Output plot file (default: show interactive)"
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.01,
        help="FDR threshold for filtering (default: 0.01)",
    )
    parser.add_argument(
        "--osprey-run",
        default=None,
        help="Osprey run name to use (default: auto-detect single run or first run)",
    )
    parser.add_argument(
        "--diann-run",
        default=None,
        help="DIA-NN run name to use (default: auto-detect single run or first run)",
    )
    args = parser.parse_args()

    # Load Osprey report
    print(f"Loading Osprey report from {args.osprey_file}")
    osprey = pd.read_parquet(
        args.osprey_file,
        columns=["Run", "Stripped.Sequence", "RT", "Q.Value", "Global.Q.Value", "Is.Decoy"],
    )
    print(f"  {len(osprey)} total entries")

    # Filter to targets only
    osprey = osprey[osprey["Is.Decoy"] == False]
    print(f"  {len(osprey)} target entries")

    # Select Osprey run
    osprey_runs = osprey["Run"].unique()
    if args.osprey_run:
        osprey = osprey[osprey["Run"] == args.osprey_run]
    elif len(osprey_runs) > 1:
        print(f"  Available runs: {osprey_runs.tolist()}")
        print(f"  Using first run: {osprey_runs[0]}")
        osprey = osprey[osprey["Run"] == osprey_runs[0]]

    osprey_run = osprey["Run"].iloc[0] if len(osprey) > 0 else "unknown"
    print(f"  Osprey run: {osprey_run}")

    # Use experiment-level q-value (Global.Q.Value) if available, else run-level
    qval_col = "Global.Q.Value" if osprey["Global.Q.Value"].notna().any() else "Q.Value"
    osprey = osprey[osprey[qval_col] <= args.fdr]
    print(f"  {len(osprey)} targets at {args.fdr:.0%} FDR ({qval_col})")

    # Keep best (lowest q-value) per stripped sequence
    osprey = osprey.sort_values(qval_col).drop_duplicates(
        subset=["Stripped.Sequence"], keep="first"
    )
    print(f"  {len(osprey)} unique peptides")

    # Load DIA-NN report
    print(f"\nLoading DIA-NN report from {args.diann_file}")
    diann = pd.read_parquet(
        args.diann_file,
        columns=["Run", "Stripped.Sequence", "RT", "Q.Value"],
    )
    print(f"  {len(diann)} total precursors")

    # Select DIA-NN run
    diann_runs = diann["Run"].unique()
    if args.diann_run:
        diann = diann[diann["Run"] == args.diann_run]
    elif len(diann_runs) > 1:
        print(f"  Available runs: {diann_runs.tolist()}")
        print(f"  Using first run: {diann_runs[0]}")
        diann = diann[diann["Run"] == diann_runs[0]]

    diann_run = diann["Run"].iloc[0] if len(diann) > 0 else "unknown"
    print(f"  DIA-NN run: {diann_run}")

    # Filter DIA-NN by FDR
    diann = diann[diann["Q.Value"] <= args.fdr]
    print(f"  {len(diann)} precursors at {args.fdr:.0%} FDR")

    # Keep best (lowest q-value) per stripped sequence
    diann = diann.sort_values("Q.Value").drop_duplicates(
        subset=["Stripped.Sequence"], keep="first"
    )
    print(f"  {len(diann)} unique peptides")

    # Merge on stripped sequence
    merged = osprey.merge(
        diann,
        on="Stripped.Sequence",
        how="inner",
        suffixes=("_osprey", "_diann"),
    )
    print(f"\n  Matched peptides: {len(merged)}")
    print(f"  Osprey-only: {len(osprey) - len(merged)}")
    print(f"  DIA-NN-only: {len(diann) - len(merged)}")

    if len(merged) == 0:
        print("No matching peptides found. Check that files are from the same organism.")
        sys.exit(1)

    # Compute statistics (RT in minutes)
    x = merged["RT_diann"].values.astype(float)   # DIA-NN RT
    y = merged["RT_osprey"].values.astype(float)   # Osprey RT
    residuals = y - x
    median_diff = np.median(residuals)
    mad = np.median(np.abs(residuals - median_diff))
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    print(f"\n  RT offset: median={median_diff:.2f} min, MAD={mad:.2f} min, R\u00b2={r2:.4f}")

    # --- Two-panel figure: density plot + residual plot ---
    fig, (ax_density, ax_resid) = plt.subplots(
        2, 1, figsize=(9, 12), height_ratios=[1, 0.6]
    )

    # ===== Panel 1: Hexbin density plot =====
    lo = min(x.min(), y.min())
    hi = max(x.max(), y.max())
    margin = (hi - lo) * 0.02

    hb = ax_density.hexbin(
        x, y,
        gridsize=200,
        cmap="turbo",
        norm=mcolors.LogNorm(),
        mincnt=1,
        rasterized=True,
    )
    cb = fig.colorbar(hb, ax=ax_density, shrink=0.8, pad=0.02)
    cb.set_label("Peptide count")

    ax_density.plot(
        [lo - margin, hi + margin],
        [lo - margin, hi + margin],
        color="white", linewidth=0.8, alpha=0.7, linestyle="--",
    )

    ax_density.set_xlim(lo - margin, hi + margin)
    ax_density.set_ylim(lo - margin, hi + margin)
    ax_density.set_xlabel("DIA-NN RT (min)")
    ax_density.set_ylabel("Osprey RT (min)")
    ax_density.set_title(
        f"RT comparison ({len(merged):,} matched peptides at {args.fdr:.0%} FDR)\n"
        f"median offset={median_diff:.2f} min, MAD={mad:.2f} min, R\u00b2={r2:.4f}"
    )
    ax_density.set_aspect("equal")

    # ===== Panel 2: Residual plot (Osprey - DIA-NN) vs DIA-NN RT =====
    hb2 = ax_resid.hexbin(
        x, residuals,
        gridsize=(200, 100),
        cmap="turbo",
        norm=mcolors.LogNorm(),
        mincnt=1,
        rasterized=True,
    )
    cb2 = fig.colorbar(hb2, ax=ax_resid, shrink=0.8, pad=0.02)
    cb2.set_label("Peptide count")

    # Reference lines
    ax_resid.axhline(median_diff, color="white", linewidth=0.8, alpha=0.7, linestyle="--")
    ax_resid.axhline(0, color="gray", linewidth=0.5, alpha=0.5, linestyle=":")

    ax_resid.set_xlabel("DIA-NN RT (min)")
    ax_resid.set_ylabel("RT difference (min)\n(Osprey \u2212 DIA-NN)")
    ax_resid.set_title("Residuals (RT difference)")

    # Print outlier counts
    n_large_resid = np.sum(np.abs(residuals - median_diff) > 5.0)
    pct_large = 100.0 * n_large_resid / len(residuals)
    print(f"  Outliers (|residual - median| > 5 min): {n_large_resid} ({pct_large:.1f}%)")

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"\nPlot saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
