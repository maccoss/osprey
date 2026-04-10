#!/usr/bin/env python3
"""Compare apex_rt for a peptide across .scores.parquet caches and the final .blib.

Helps localize whether divergent RTs in blib output are caused by:
  (1) a blib writing bug (parquet RTs clustered, blib RTs scattered),
  (2) an upstream scoring/reconciliation bug (parquet RTs already scattered),
  (3) expected forced-integration behavior (only one run has non-NULL retentionTime).

Note: q-values live in FDR sidecar files (.fdr_scores.bin), not the Parquet caches.
We pull q-values from the blib's OspreyRunScores / OspreyExperimentScores tables.

Usage:
    python scripts/verify_blib_rt_consistency.py <cache_dir> <blib_path> <peptide> [peptide ...]
"""
import sqlite3
import sys
from pathlib import Path

import pyarrow.parquet as pq


PARQUET_COLS = [
    "modified_sequence",
    "sequence",
    "charge",
    "is_decoy",
    "apex_rt",
    "start_rt",
    "end_rt",
    "file_name",
]


def read_parquet_rows(cache_dir: Path, peptide: str):
    """Read all target rows for the given (unmodified) peptide from every parquet cache."""
    rows = []
    for pq_file in sorted(cache_dir.glob("*.scores.parquet")):
        try:
            tbl = pq.read_table(pq_file, columns=PARQUET_COLS)
        except Exception as e:
            print(f"  ! failed to read {pq_file.name}: {e}")
            continue
        df = tbl.to_pandas()
        # Match either modified_sequence or stripped sequence
        hits = df[
            (~df["is_decoy"])
            & ((df["modified_sequence"] == peptide) | (df["sequence"] == peptide))
        ]
        for _, r in hits.iterrows():
            rows.append(
                {
                    "cache_file": pq_file.stem.replace(".scores", ""),
                    "file_name": str(r["file_name"]),
                    "modseq": str(r["modified_sequence"]),
                    "charge": int(r["charge"]),
                    "apex_rt": float(r["apex_rt"]),
                    "start_rt": float(r["start_rt"]),
                    "end_rt": float(r["end_rt"]),
                }
            )
    return rows


def read_blib_rows(blib_path: Path, peptide: str):
    """Pull RetentionTimes, OspreyPeakBoundaries, and run/experiment q-values for the peptide."""
    conn = sqlite3.connect(str(blib_path))
    cur = conn.cursor()

    # Find matching RefSpectra (match either peptideSeq or peptideModSeq)
    cur.execute(
        """
        SELECT id, peptideSeq, peptideModSeq, precursorCharge
        FROM RefSpectra
        WHERE peptideSeq = ? OR peptideModSeq = ?
        """,
        (peptide, peptide),
    )
    refs = cur.fetchall()
    if not refs:
        conn.close()
        return []

    ref_ids = [r[0] for r in refs]
    placeholders = ",".join("?" * len(ref_ids))

    # RetentionTimes joined with source files
    cur.execute(
        f"""
        SELECT srf.fileName, rt.RefSpectraID, rt.retentionTime, rt.startTime, rt.endTime, rt.score
        FROM RetentionTimes rt
        JOIN SpectrumSourceFiles srf ON rt.SpectrumSourceID = srf.id
        WHERE rt.RefSpectraID IN ({placeholders})
        ORDER BY srf.fileName
        """,
        ref_ids,
    )
    rt_rows = cur.fetchall()

    # OspreyPeakBoundaries (per-file boundaries)
    cur.execute(
        f"""
        SELECT FileName, RefSpectraID, ApexRT, StartRT, EndRT, IntegratedArea
        FROM OspreyPeakBoundaries
        WHERE RefSpectraID IN ({placeholders})
        ORDER BY FileName
        """,
        ref_ids,
    )
    pb_rows = cur.fetchall()

    # OspreyRunScores
    cur.execute(
        f"""
        SELECT FileName, RefSpectraID, RunQValue, DiscriminantScore
        FROM OspreyRunScores
        WHERE RefSpectraID IN ({placeholders})
        """,
        ref_ids,
    )
    run_scores = {(Path(r[0]).stem, r[1]): (r[2], r[3]) for r in cur.fetchall()}

    # OspreyExperimentScores
    cur.execute(
        f"""
        SELECT RefSpectraID, ExperimentQValue, NRunsDetected, NRunsSearched
        FROM OspreyExperimentScores
        WHERE RefSpectraID IN ({placeholders})
        """,
        ref_ids,
    )
    exp_scores = {r[0]: (r[1], r[2], r[3]) for r in cur.fetchall()}

    conn.close()

    ref_meta = {r[0]: {"peptideSeq": r[1], "modseq": r[2], "charge": r[3]} for r in refs}
    return {
        "refs": ref_meta,
        "retention_times": rt_rows,
        "peak_boundaries": pb_rows,
        "run_scores": run_scores,
        "experiment_scores": exp_scores,
    }


def summarize(label: str, vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        print(f"  {label}: no non-NULL values")
        return
    span = max(vals) - min(vals)
    print(
        f"  {label}: n={len(vals)} min={min(vals):.2f} max={max(vals):.2f} span={span:.2f}"
    )


def report_peptide(cache_dir: Path, blib_path: Path, peptide: str):
    print(f"\n========== {peptide} ==========")

    print(f"\n--- Parquet caches ({cache_dir}) ---")
    pq_rows = read_parquet_rows(cache_dir, peptide)
    if not pq_rows:
        print("  (no target hits in any parquet cache)")
    for r in pq_rows:
        print(
            f"  {r['cache_file']:50s} z={r['charge']} apex={r['apex_rt']:7.2f} "
            f"[{r['start_rt']:7.2f}-{r['end_rt']:7.2f}]  modseq={r['modseq']}"
        )
    summarize("parquet apex_rt", [r["apex_rt"] for r in pq_rows])

    print(f"\n--- Blib ({blib_path.name}) ---")
    blib = read_blib_rows(blib_path, peptide)
    if not blib:
        print("  (no matching RefSpectra)")
        return

    print(f"  RefSpectra matches: {len(blib['refs'])}")
    for rid, meta in blib["refs"].items():
        exp = blib["experiment_scores"].get(rid)
        exp_str = (
            f"exp_q={exp[0]:.4f} nDet={exp[1]}/{exp[2]}" if exp else "exp_q=(none)"
        )
        print(
            f"    refid={rid} z={meta['charge']} modseq={meta['modseq']}  {exp_str}"
        )

    print("\n  RetentionTimes table:")
    for fname, rid, rt, start, end, score in blib["retention_times"]:
        rt_str = f"{rt:7.2f}" if rt is not None else "   NULL"
        rs = blib["run_scores"].get((Path(fname).stem, rid))
        rq = f"run_q={rs[0]:.4f}" if rs else "run_q=?"
        print(
            f"    {Path(fname).stem:50s} refid={rid} apex={rt_str} "
            f"[{start:7.2f}-{end:7.2f}] {rq}"
        )

    print("\n  OspreyPeakBoundaries table:")
    for fname, rid, apex, start, end, area in blib["peak_boundaries"]:
        print(
            f"    {Path(fname).stem:50s} refid={rid} apex={apex:7.2f} "
            f"[{start:7.2f}-{end:7.2f}] area={area:.2e}"
        )

    rt_vals = [row[2] for row in blib["retention_times"] if row[2] is not None]
    pb_vals = [row[2] for row in blib["peak_boundaries"]]
    summarize("blib RetentionTimes.retentionTime (non-NULL)", rt_vals)
    summarize("blib OspreyPeakBoundaries.ApexRT", pb_vals)

    n_rt = len(blib["retention_times"])
    n_null = sum(1 for r in blib["retention_times"] if r[2] is None)
    print(f"  RetentionTimes NULL count: {n_null}/{n_rt}")

    print("\n--- Diagnosis ---")
    pq_apex = [r["apex_rt"] for r in pq_rows]
    if pq_apex and pb_vals:
        pq_span = max(pq_apex) - min(pq_apex)
        pb_span = max(pb_vals) - min(pb_vals)
        print(f"  parquet apex span: {pq_span:.2f} min")
        print(f"  blib OspreyPeakBoundaries apex span: {pb_span:.2f} min")
        if pq_span < 1.0 and pb_span > 2.0:
            print("  >>> Parquet RTs clustered, blib RTs scattered: BLIB WRITING BUG")
        elif pq_span > 2.0:
            print("  >>> Parquet RTs already scattered: UPSTREAM SCORING/RECONCILIATION BUG")
        elif n_null == n_rt - 1 and n_rt > 1:
            print("  >>> Only one run has non-NULL retentionTime: expected forced-integration")
        else:
            print("  >>> Inconclusive; inspect rows above")


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    cache_dir = Path(sys.argv[1])
    blib_path = Path(sys.argv[2])
    peptides = sys.argv[3:]
    for pep in peptides:
        report_peptide(cache_dir, blib_path, pep)


if __name__ == "__main__":
    main()
