#!/usr/bin/env python3
"""Compare retention times between two spectral libraries.

Matches peptides by stripped sequence and creates a density plot + residual plot
comparing predicted iRT (DIA-NN) vs predicted RT (Carafe/SkylineAI). A linear
fit converts between the two scales so residuals reflect genuine disagreements.

Usage:
    python compare_library_rt.py diann_lib.parquet carafe_lib.tsv [--output plot.png]

Example:
    python scripts/compare_library_rt.py \
        example_test_data/astral/carafe/diann_train/report-lib.parquet \
        example_test_data/astral/SkylineAI_spectral_library.tsv
"""

import argparse
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_diann_library(path: str) -> pd.DataFrame:
    """Load DIA-NN library parquet and deduplicate to one iRT per stripped sequence."""
    df = pd.read_parquet(path, columns=["Stripped.Sequence", "RT", "Q.Value"])
    # Keep best q-value per peptide
    df = df.sort_values("Q.Value").drop_duplicates(subset=["Stripped.Sequence"], keep="first")
    df = df.rename(columns={"Stripped.Sequence": "StrippedSequence", "RT": "iRT_diann"})
    return df[["StrippedSequence", "iRT_diann"]]


def load_carafe_library(path: str) -> pd.DataFrame:
    """Load Carafe/SkylineAI TSV library and deduplicate to one RT per stripped sequence."""
    df = pd.read_csv(path, sep="\t", usecols=["StrippedPeptide", "Tr_recalibrated"])
    df = df.drop_duplicates(subset=["StrippedPeptide"], keep="first")
    df = df.rename(columns={"StrippedPeptide": "StrippedSequence", "Tr_recalibrated": "RT_carafe"})
    return df[["StrippedSequence", "RT_carafe"]]


def main():
    parser = argparse.ArgumentParser(
        description="Compare iRT (DIA-NN) vs predicted RT (Carafe/SkylineAI) libraries"
    )
    parser.add_argument("diann_lib", help="DIA-NN library parquet (report-lib.parquet)")
    parser.add_argument("carafe_lib", help="Carafe/SkylineAI library TSV")
    parser.add_argument(
        "--output", "-o", default=None, help="Output plot file (default: show interactive)"
    )
    args = parser.parse_args()

    # Load libraries
    print(f"Loading DIA-NN library from {args.diann_lib}")
    diann = load_diann_library(args.diann_lib)
    print(f"  {len(diann)} unique peptides")

    print(f"Loading Carafe library from {args.carafe_lib}")
    carafe = load_carafe_library(args.carafe_lib)
    print(f"  {len(carafe)} unique peptides")

    # Merge
    merged = diann.merge(carafe, on="StrippedSequence", how="inner")
    print(f"\n  Matched peptides: {len(merged)}")
    print(f"  DIA-NN-only: {len(diann) - len(merged)}")
    print(f"  Carafe-only: {len(carafe) - len(merged)}")

    if len(merged) == 0:
        print("No matching peptides found.")
        sys.exit(1)

    x = merged["RT_carafe"].values.astype(float)
    y = merged["iRT_diann"].values.astype(float)

    # Fit linear model: iRT = slope * RT + intercept
    slope, intercept = np.polyfit(x, y, 1)
    y_predicted = slope * x + intercept
    residuals = y - y_predicted

    r2 = np.corrcoef(x, y)[0, 1] ** 2
    mad = np.median(np.abs(residuals))
    print(f"\n  Linear fit: iRT = {slope:.4f} * RT + {intercept:.4f}")
    print(f"  R\u00b2={r2:.6f}, residual MAD={mad:.2f} iRT units")

    # --- Two-panel figure ---
    fig, (ax_density, ax_resid) = plt.subplots(
        2, 1, figsize=(9, 12), height_ratios=[1, 0.6]
    )

    # ===== Panel 1: Hexbin density (RT vs iRT) =====
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
        label=f"iRT = {slope:.2f}\u00d7RT + {intercept:.2f}",
    )
    ax_density.legend(loc="upper left", fontsize=9)

    ax_density.set_xlabel("Carafe predicted RT (min)")
    ax_density.set_ylabel("DIA-NN predicted iRT")
    ax_density.set_title(
        f"Library RT comparison ({len(merged):,} matched peptides)\n"
        f"iRT = {slope:.2f}\u00d7RT + {intercept:.2f}, R\u00b2={r2:.6f}"
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

    ax_resid.set_xlabel("Carafe predicted RT (min)")
    ax_resid.set_ylabel("Residual (iRT units)\n(DIA-NN iRT \u2212 predicted)")
    ax_resid.set_title(f"Residuals from linear fit (MAD={mad:.2f} iRT units)")

    n_large_resid = np.sum(np.abs(residuals) > 3 * mad)
    pct_large = 100.0 * n_large_resid / len(residuals)
    print(f"  Outliers (|residual| > 3\u00d7MAD = {3*mad:.2f}): {n_large_resid} ({pct_large:.1f}%)")

    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=150)
        print(f"\nPlot saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
