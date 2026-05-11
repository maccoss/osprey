"""Build FDRBench input TSVs from DIA-NN report.parquet or Osprey .blib output.

FDRBench (https://github.com/Noble-Lab/FDRBench) evaluates FDR control quality by
counting entrapment hits in search results. It consumes a single TSV per
evaluation level (peptide/precursor or protein) and looks up entrapment status
via a separate reference file supplied at invocation time. This script converts
DIA-NN parquet reports and Osprey blib outputs into that TSV format.

Output formats:

Peptide / precursor level (tab-separated):
    peptide  mod_peptide  charge  q_value  score  protein  [run]

Protein level:
    protein  q_value  score

``score`` is emitted as the raw upstream discriminating score that drives the
q-value (DIA-NN ``Evidence`` for peptide/precursor and ``max(Evidence)`` for
protein; Osprey ``DiscriminantScore`` for peptide/precursor and
``best_peptide_score`` for protein). Higher = better; invoke FDRBench with
``-score 'score:1'``. If Osprey's blib was written without
``DiscriminantScore`` populated, the script falls back to ``1 - q_value``
with a warning — use the native Osprey ``--fdrbench`` flag to avoid this.

Usage:
    build_fdrbench_input.py diann  -i report.parquet --level precursor -o out.tsv
    build_fdrbench_input.py osprey -i osprey.blib    --level protein   -o out.tsv [--proteins-csv path]

Add ``--per-run`` to emit one row per (precursor, run) using the per-run
q-value instead of the experiment / best-across-runs aggregation.
"""

from __future__ import annotations

import argparse
import logging
import sqlite3
import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

PEPTIDE_COLUMNS = ["peptide", "mod_peptide", "charge", "q_value", "score", "protein"]
PEPTIDE_COLUMNS_PER_RUN = PEPTIDE_COLUMNS + ["run"]
PROTEIN_COLUMNS = ["protein", "q_value", "score"]


@dataclass
class Args:
    source: str
    input_path: Path
    output_path: Path
    level: str
    per_run: bool
    proteins_csv: Path | None


def parse_args(argv: list[str] | None = None) -> Args:
    """Parse command-line arguments into an Args dataclass."""
    parser = argparse.ArgumentParser(
        description="Build FDRBench input TSVs from DIA-NN or Osprey results.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example FDRBench invocation after this script produces input.tsv:\n"
            "  java -jar fdrbench-0.0.1.jar -i input.tsv -level precursor \\\n"
            "      -pep entrapment_peptides.txt -fold 1 -score 'score:1' -o fdp.csv"
        ),
    )
    subparsers = parser.add_subparsers(dest="source", required=True)

    for name, input_help in [
        ("diann", "Path to DIA-NN report.parquet"),
        ("osprey", "Path to Osprey .blib file"),
    ]:
        sub = subparsers.add_parser(name, help=f"Convert {name} output")
        sub.add_argument("-i", "--input", required=True, type=Path, help=input_help)
        sub.add_argument("-o", "--output", required=True, type=Path, help="Output TSV path")
        sub.add_argument(
            "--level",
            required=True,
            choices=["peptide", "precursor", "protein"],
            help="FDRBench evaluation level",
        )
        sub.add_argument(
            "--per-run",
            action="store_true",
            help="Emit one row per (precursor, run) with per-run q-values "
            "(ignored for --level protein)",
        )
        if name == "osprey":
            sub.add_argument(
                "--proteins-csv",
                type=Path,
                default=None,
                help="Override path to Osprey proteins.csv "
                "(default: <blib-stem>.proteins.csv next to the blib)",
            )

    ns = parser.parse_args(argv)
    return Args(
        source=ns.source,
        input_path=ns.input,
        output_path=ns.output,
        level=ns.level,
        per_run=ns.per_run,
        proteins_csv=getattr(ns, "proteins_csv", None),
    )


def best_by_qvalue(df: pd.DataFrame, key_cols: list[str]) -> pd.DataFrame:
    """Reduce df to one row per key (min q_value, ties broken by max score)."""
    return (
        df.sort_values(["q_value", "score"], ascending=[True, False])
        .drop_duplicates(subset=key_cols, keep="first")
        .reset_index(drop=True)
    )


def build_diann_peptide(parquet_path: Path, level: str, per_run: bool) -> pd.DataFrame:
    """Read a DIA-NN report and return a peptide/precursor-level DataFrame.

    Args:
        parquet_path: Path to DIA-NN report.parquet.
        level: "peptide" or "precursor".
        per_run: If True, keep one row per (precursor, run); otherwise dedup
            across runs to the best q-value.

    Returns:
        DataFrame with FDRBench peptide-level columns in canonical order.
    """
    logger.info("Reading DIA-NN parquet: %s", parquet_path)
    df = pd.read_parquet(parquet_path)
    logger.info("Read %d rows", len(df))

    # Drop only DIA-NN-flagged decoys. Entrapment peptides have Decoy == 0 and
    # carry a `_p_target` suffix on their protein accession; they MUST flow
    # through so FDRBench can use them to compute true-FDR via entrapment counts.
    if "Decoy" in df.columns:
        n_before = len(df)
        df = df[df["Decoy"] == 0]
        logger.debug("Filtered decoys: %d -> %d rows", n_before, len(df))

    df = df.rename(
        columns={
            "Stripped.Sequence": "peptide",
            "Modified.Sequence": "mod_peptide",
            "Precursor.Charge": "charge",
            "Q.Value": "q_value",
            "Protein.Ids": "protein",
            "Run": "run",
            # Evidence is DIA-NN's neural-network discriminant score that drives
            # PEP and ultimately Q.Value. Higher = better. FDRBench needs this
            # raw upstream score (not q_value) to re-rank and count entrapment.
            "Evidence": "score",
        }
    )

    if level == "peptide":
        key = ["mod_peptide"]
    else:
        key = ["mod_peptide", "charge"]
    if per_run:
        key = key + ["run"]

    df = best_by_qvalue(df[["peptide", "mod_peptide", "charge", "q_value", "score", "protein", "run"]], key)
    return df


def build_diann_protein(parquet_path: Path) -> pd.DataFrame:
    """Build a protein-level DataFrame from a DIA-NN parquet report.

    One row per Protein.Group, with the full Protein.Ids list emitted as the
    ``protein`` column so entrapment markers (``_p_target``) survive.
    """
    logger.info("Reading DIA-NN parquet: %s", parquet_path)
    df = pd.read_parquet(
        parquet_path,
        columns=["Decoy", "Protein.Group", "Protein.Ids", "PG.Q.Value", "Evidence"],
    )
    if "Decoy" in df.columns:
        df = df[df["Decoy"] == 0]

    # One row per Protein.Group: take the row with min PG.Q.Value (most confident
    # protein). Score is max Evidence across the group's precursors — the best
    # peptide's NN discriminant ranks proteins for FDRBench's entrapment counting.
    grouped = df.groupby("Protein.Group", as_index=False).agg(
        q_value=("PG.Q.Value", "min"),
        score=("Evidence", "max"),
        protein=("Protein.Ids", "first"),
    )
    return grouped[PROTEIN_COLUMNS].reset_index(drop=True)


def build_osprey_peptide(blib_path: Path, level: str, per_run: bool) -> pd.DataFrame:
    """Read an Osprey blib and return a peptide/precursor-level DataFrame.

    Joins RefSpectra with OspreyRunScores (for --per-run) or OspreyExperimentScores
    (default), and group-concatenates protein accessions per spectrum.
    """
    logger.info("Reading Osprey blib: %s", blib_path)
    con = sqlite3.connect(blib_path)
    try:
        # Always join OspreyRunScores to pull DiscriminantScore (the raw SVM
        # discriminant, which is what FDRBench needs for re-ranking). For
        # experiment-level (default) we still take ExperimentQValue, but the
        # discriminant is per-RefSpectra-row regardless.
        if per_run:
            qvalue_expr = "s.RunQValue"
            run_select = ", s.FileName AS run"
            group_by_extra = ", s.FileName"
        else:
            qvalue_expr = "es.ExperimentQValue"
            run_select = ""
            group_by_extra = ""

        sql = f"""
            SELECT
                r.id                       AS spec_id,
                r.peptideSeq               AS peptide,
                r.peptideModSeq            AS mod_peptide,
                r.precursorCharge          AS charge,
                {qvalue_expr}              AS q_value,
                MAX(s.DiscriminantScore)   AS svm_score,
                GROUP_CONCAT(p.accession, ';') AS protein
                {run_select}
            FROM RefSpectra r
            JOIN OspreyRunScores s ON s.RefSpectraID = r.id
            JOIN OspreyExperimentScores es ON es.RefSpectraID = r.id
            LEFT JOIN RefSpectraProteins rp ON rp.RefSpectraID = r.id
            LEFT JOIN Proteins p ON p.id = rp.ProteinID
            GROUP BY r.id{group_by_extra}
        """
        df = pd.read_sql_query(sql, con)
    finally:
        con.close()

    logger.info("Read %d peptide rows from blib", len(df))

    # If Osprey populated DiscriminantScore (raw SVM score), use it — that's
    # the input score that drove the q-value. If it's all zero (Osprey's blib
    # writer historically did not persist it), fall back to 1 - q_value and
    # warn loudly: 1 - q_value has massive ties and degrades FDRBench's
    # entrapment-FDR estimate. The native Osprey `--fdrbench` flag avoids this.
    if (df["svm_score"].fillna(0.0) != 0.0).any():
        df["score"] = df["svm_score"]
    else:
        logger.warning(
            "Osprey blib has no DiscriminantScore populated (all zeros). "
            "Falling back to score = 1 - q_value, which has heavy ties and "
            "weakens FDRBench's true-FDR estimate. Use Osprey's native "
            "--fdrbench flag for raw SVM scores."
        )
        df["score"] = 1.0 - df["q_value"]
    df = df.drop(columns=["spec_id", "svm_score"])

    # Guard against any decoy contamination (Osprey filters these before blib output,
    # but the convention is DECOY_ prefix on the peptide sequence). Entrapment
    # targets carry `_p_target` on their protein accession, NOT on the peptide
    # sequence, and must pass through unfiltered so FDRBench can score them.
    n_before = len(df)
    df = df[~df["peptide"].astype(str).str.startswith("DECOY_")]
    if len(df) != n_before:
        logger.warning("Dropped %d DECOY_ rows from blib output", n_before - len(df))

    if level == "peptide":
        key = ["mod_peptide"]
    else:
        key = ["mod_peptide", "charge"]
    if per_run:
        key = key + ["run"]
    df = best_by_qvalue(df, key)
    return df


def build_osprey_protein(proteins_csv_path: Path) -> pd.DataFrame:
    """Build a protein-level DataFrame from an Osprey proteins.csv.

    Osprey writes one row per protein group with already-aggregated q-value
    and accessions; the entrapment ``_p_target`` suffix is preserved in
    ``protein_accessions``.
    """
    logger.info("Reading Osprey proteins CSV: %s", proteins_csv_path)
    df = pd.read_csv(proteins_csv_path)
    # best_peptide_score is the raw SVM discriminant of the protein group's
    # best peptide — the upstream score that drives protein_qvalue. Use it
    # directly so FDRBench can re-rank meaningfully.
    df = df.rename(
        columns={
            "protein_accessions": "protein",
            "protein_qvalue": "q_value",
            "best_peptide_score": "score",
        }
    )
    return df[PROTEIN_COLUMNS].reset_index(drop=True)


def write_tsv(df: pd.DataFrame, out_path: Path, level: str, per_run: bool) -> None:
    """Write the DataFrame as a FDRBench-compatible TSV with canonical column order."""
    if level == "protein":
        cols = PROTEIN_COLUMNS
    elif per_run:
        cols = PEPTIDE_COLUMNS_PER_RUN
    else:
        cols = PEPTIDE_COLUMNS

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote %d rows to %s", len(df), out_path)


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args(argv)

    if args.source == "diann":
        if args.level == "protein":
            df = build_diann_protein(args.input_path)
        else:
            df = build_diann_peptide(args.input_path, args.level, args.per_run)
    else:
        if args.level == "protein":
            proteins_csv = args.proteins_csv or args.input_path.with_suffix(".proteins.csv")
            if not proteins_csv.exists():
                # Try removing the .blib suffix and appending .proteins.csv (handles
                # blibs with multi-dotted stems like "foo.bar.blib").
                alt = args.input_path.parent / (args.input_path.stem + ".proteins.csv")
                if alt.exists():
                    proteins_csv = alt
            if not proteins_csv.exists():
                logger.error("Osprey proteins.csv not found: %s", proteins_csv)
                return 1
            df = build_osprey_protein(proteins_csv)
        else:
            df = build_osprey_peptide(args.input_path, args.level, args.per_run)

    write_tsv(df, args.output_path, args.level, args.per_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
