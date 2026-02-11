#!/usr/bin/env python3
"""Compare retention times between Osprey (mokapot) and DIA-NN for matched peptides.

Reads a mokapot peptides output file and a DIA-NN report parquet, matches peptides
by stripped sequence, and creates a density plot + residual plot comparing RTs.

Osprey scan numbers are mapped to RT using the mzML file. DIA-NN RT comes directly
from the report.

Usage:
    python compare_scan_numbers.py mokapot_peptides.txt diann_report.parquet mzml_file [--output plot.png]

Example:
    python scripts/compare_scan_numbers.py \
        example_test_data/astral/mokapot/run_level/file.mokapot.peptides.txt \
        example_test_data/astral/carafe/diann_train/report.parquet \
        example_test_data/astral/file.mzML
"""

import argparse
import re
import sys

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lxml import etree


def strip_mokapot_peptide(peptide: str) -> str:
    """Convert mokapot peptide format to stripped sequence.

    Mokapot format: '-.PEPTC[UniMod:4]IDE.-'
    Returns: 'PEPTCIDE'
    """
    seq = peptide
    if seq.startswith("-."):
        seq = seq[2:]
    if seq.endswith(".-"):
        seq = seq[:-2]
    seq = re.sub(r"\[UniMod:\d+\]", "", seq)
    return seq


def build_scan_rt_map(mzml_path: str) -> dict:
    """Build scan number -> RT (minutes) map from mzML using streaming parse.

    Handles large files efficiently by using iterparse and clearing elements.
    """
    scan_rt = {}
    ns = "http://psi.hupo.org/ms/mzml"
    spectrum_tag = f"{{{ns}}}spectrum"
    cv_tag = f"{{{ns}}}cvParam"
    scan_tag = f"{{{ns}}}scan"

    n_spectra = 0
    for event, elem in etree.iterparse(mzml_path, events=("end",), tag=spectrum_tag):
        # Extract scan number from id attribute
        spec_id = elem.get("id", "")
        scan_match = re.search(r"scan=(\d+)", spec_id)
        if not scan_match:
            elem.clear()
            continue

        scan_num = int(scan_match.group(1))

        # Find RT in scan/cvParam
        rt_min = None
        for scan_elem in elem.iter(scan_tag):
            for cv in scan_elem.iter(cv_tag):
                if cv.get("accession") == "MS:1000016":
                    rt_val = float(cv.get("value"))
                    unit = cv.get("unitAccession", "")
                    if unit == "UO:0000010":  # seconds
                        rt_min = rt_val / 60.0
                    else:  # assume minutes
                        rt_min = rt_val
                    break
            if rt_min is not None:
                break

        if rt_min is not None:
            scan_rt[scan_num] = rt_min

        n_spectra += 1
        if n_spectra % 50000 == 0:
            print(f"    Parsed {n_spectra:,} spectra...", flush=True)

        # Free memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    print(f"    {n_spectra:,} spectra total, {len(scan_rt):,} with RT")
    return scan_rt


def main():
    parser = argparse.ArgumentParser(
        description="Compare retention times: Osprey (mokapot) vs DIA-NN"
    )
    parser.add_argument("mokapot_file", help="Mokapot peptides TSV file")
    parser.add_argument("diann_file", help="DIA-NN report parquet file")
    parser.add_argument("mzml_file", help="mzML file for scan number to RT mapping")
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
        "--diann-run",
        default=None,
        help="DIA-NN run name to use (default: auto-detect single run or first run)",
    )
    args = parser.parse_args()

    # Build scan -> RT map from mzML
    print(f"Building scan-to-RT map from {args.mzml_file}")
    scan_rt = build_scan_rt_map(args.mzml_file)

    # Load mokapot peptides
    print(f"Loading mokapot peptides from {args.mokapot_file}")
    mok = pd.read_csv(args.mokapot_file, sep="\t")
    print(f"  {len(mok)} total peptides")

    # Filter by FDR and target (Label=True)
    mok = mok[(mok["Label"] == True) & (mok["mokapot q-value"] <= args.fdr)]
    print(f"  {len(mok)} target peptides at {args.fdr:.0%} FDR")

    # Strip peptide sequences
    mok["StrippedSequence"] = mok["Peptide"].apply(strip_mokapot_peptide)

    # Map Osprey scan numbers to RT
    mok["RT_osprey"] = mok["ScanNr"].map(scan_rt)
    n_mapped = mok["RT_osprey"].notna().sum()
    n_unmapped = mok["RT_osprey"].isna().sum()
    print(f"  Mapped {n_mapped} scan numbers to RT ({n_unmapped} unmapped)")
    mok = mok.dropna(subset=["RT_osprey"])

    # Load DIA-NN report (with RT column)
    print(f"Loading DIA-NN report from {args.diann_file}")
    diann = pd.read_parquet(
        args.diann_file,
        columns=["Run", "Stripped.Sequence", "RT", "Q.Value"],
    )
    print(f"  {len(diann)} total precursors")

    # Select run
    runs = diann["Run"].unique()
    if args.diann_run:
        diann = diann[diann["Run"] == args.diann_run]
    elif len(runs) == 1:
        pass
    else:
        print(f"  Available runs: {runs.tolist()}")
        print(f"  Using first run: {runs[0]}")
        diann = diann[diann["Run"] == runs[0]]

    diann_run = diann["Run"].iloc[0] if len(diann) > 0 else "unknown"
    print(f"  DIA-NN run: {diann_run}")

    # Filter DIA-NN by FDR
    diann = diann[diann["Q.Value"] <= args.fdr]
    print(f"  {len(diann)} precursors at {args.fdr:.0%} FDR")

    # For DIA-NN, keep best (lowest q-value) per stripped sequence
    diann = diann.sort_values("Q.Value").drop_duplicates(
        subset=["Stripped.Sequence"], keep="first"
    )
    print(f"  {len(diann)} unique peptides")

    # Merge on stripped sequence
    merged = mok.merge(
        diann,
        left_on="StrippedSequence",
        right_on="Stripped.Sequence",
        how="inner",
    )
    print(f"\n  Matched peptides: {len(merged)}")
    print(f"  Osprey-only: {len(mok) - len(merged)}")
    print(f"  DIA-NN-only: {len(diann) - len(merged)}")

    if len(merged) == 0:
        print("No matching peptides found. Check that files are from the same organism.")
        sys.exit(1)

    # Extract mokapot run name from SpecId
    mok_run = mok["SpecId"].iloc[0].rsplit("_", 1)[0] if "SpecId" in mok.columns else "Osprey"

    # Compute statistics (RT in minutes)
    x = merged["RT"].values.astype(float)  # DIA-NN RT
    y = merged["RT_osprey"].values.astype(float)  # Osprey RT from mzML
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
    ax_density.set_xlabel(f"DIA-NN RT (min)")
    ax_density.set_ylabel(f"Osprey RT (min)")
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

    ax_resid.set_xlabel(f"DIA-NN RT (min)")
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
