#!/usr/bin/env python3
"""Compare retention times between Osprey results and a spectral library.

Matches peptides by stripped sequence and creates a density plot + residual plot
comparing observed RT (Osprey) vs predicted RT (Carafe/SkylineAI library). A linear
fit converts between the two scales so residuals reflect genuine disagreements.

Usage:
    python compare_library_rt.py osprey_report.parquet carafe_lib.tsv [--output plot.png]

Example:
    python scripts/compare_library_rt.py \
        report.parquet \
        example_test_data/astral/SkylineAI_spectral_library.tsv
"""

import argparse
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_osprey_report(path: str, qvalue_threshold: float = 0.01) -> pd.DataFrame:
    """Load Osprey scores parquet and deduplicate to one RT per stripped sequence.

    Filters to targets passing the q-value threshold, then keeps the best-scoring
    observation per stripped sequence.
    """
    columns = ["Stripped.Sequence", "RT", "Q.Value", "Global.Q.Value", "Is.Decoy", "Score"]
    df = pd.read_parquet(path, columns=columns)

    # Filter to targets passing FDR threshold
    df = df[~df["Is.Decoy"]]
    qval = df[["Q.Value", "Global.Q.Value"]].max(axis=1)
    df = df[qval <= qvalue_threshold].copy()

    # Keep best score per peptide
    df = df.sort_values("Score", ascending=False).drop_duplicates(
        subset=["Stripped.Sequence"], keep="first"
    )
    df = df.rename(columns={"Stripped.Sequence": "StrippedSequence", "RT": "RT_osprey"})
    return df[["StrippedSequence", "RT_osprey"]]


def load_carafe_library(path: str) -> pd.DataFrame:
    """Load Carafe/SkylineAI TSV library and deduplicate to one RT per stripped sequence."""
    df = pd.read_csv(path, sep="\t", usecols=["StrippedPeptide", "Tr_recalibrated"])
    df = df.drop_duplicates(subset=["StrippedPeptide"], keep="first")
    df = df.rename(
        columns={"StrippedPeptide": "StrippedSequence", "Tr_recalibrated": "RT_library"}
    )
    return df[["StrippedSequence", "RT_library"]]


def main():
    parser = argparse.ArgumentParser(
        description="Compare observed RT (Osprey) vs predicted RT (Carafe/SkylineAI) library"
    )
    parser.add_argument("osprey_report", help="Osprey report parquet (report.parquet)")
    parser.add_argument("carafe_lib", help="Carafe/SkylineAI library TSV")
    parser.add_argument(
        "--output", "-o", default=None, help="Output plot file (default: show interactive)"
    )
    parser.add_argument(
        "--qvalue", "-q", type=float, default=0.01, help="Q-value threshold (default: 0.01)"
    )
    args = parser.parse_args()

    # Load data
    print(f"Loading Osprey report from {args.osprey_report}")
    osprey = load_osprey_report(args.osprey_report, args.qvalue)
    print(f"  {len(osprey)} unique peptides at {args.qvalue:.0%} FDR")

    print(f"Loading Carafe library from {args.carafe_lib}")
    carafe = load_carafe_library(args.carafe_lib)
    print(f"  {len(carafe)} unique peptides")

    # Merge
    merged = osprey.merge(carafe, on="StrippedSequence", how="inner")
    print(f"\n  Matched peptides: {len(merged)}")
    print(f"  Osprey-only: {len(osprey) - len(merged)}")
    print(f"  Library-only: {len(carafe) - len(merged)}")

    if len(merged) == 0:
        print("No matching peptides found.")
        sys.exit(1)

    x = merged["RT_library"].values.astype(float)
    y = merged["RT_osprey"].values.astype(float)

    # Fit linear model: RT_osprey = slope * RT_library + intercept
    slope, intercept = np.polyfit(x, y, 1)
    y_predicted = slope * x + intercept
    residuals = y - y_predicted

    r2 = np.corrcoef(x, y)[0, 1] ** 2
    mad = np.median(np.abs(residuals))
    print(f"\n  Linear fit: RT_osprey = {slope:.4f} \u00d7 RT_library + {intercept:.4f}")
    print(f"  R\u00b2={r2:.6f}, residual MAD={mad:.4f} min")

    # --- Two-panel figure ---
    fig, (ax_density, ax_resid) = plt.subplots(
        2, 1, figsize=(9, 12), height_ratios=[1, 0.6]
    )

    # ===== Panel 1: Hexbin density (library RT vs observed RT) =====
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

    # Plot linear fit line
    x_line = np.array([x.min(), x.max()])
    y_line = slope * x_line + intercept
    ax_density.plot(
        x_line, y_line,
        color="white", linewidth=1.0, alpha=0.8, linestyle="--",
        label=f"RT = {slope:.2f}\u00d7lib + {intercept:.2f}",
    )
    ax_density.legend(loc="upper left", fontsize=9)

    ax_density.set_xlabel("Library predicted RT (min)")
    ax_density.set_ylabel("Osprey observed RT (min)")
    ax_density.set_title(
        f"Osprey vs Library RT ({len(merged):,} matched peptides at {args.qvalue:.0%} FDR)\n"
        f"RT = {slope:.2f}\u00d7lib + {intercept:.2f}, R\u00b2={r2:.6f}"
    )

    # ===== Panel 2: Residuals from linear fit =====
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

    ax_resid.axhline(0, color="white", linewidth=0.8, alpha=0.7, linestyle="--")

    ax_resid.set_xlabel("Library predicted RT (min)")
    ax_resid.set_ylabel("Residual (min)\n(observed \u2212 predicted)")
    ax_resid.set_title(f"Residuals from linear fit (MAD={mad:.4f} min)")

    n_large_resid = np.sum(np.abs(residuals) > 3 * mad)
    pct_large = 100.0 * n_large_resid / len(residuals)
    print(
        f"  Outliers (|residual| > 3\u00d7MAD = {3*mad:.4f} min): "
        f"{n_large_resid} ({pct_large:.1f}%)"
    )

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"\nPlot saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
