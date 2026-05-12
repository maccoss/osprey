"""Fix the ``Decoy`` column in a DIA-NN-style library TSV based on protein-prefix detection.

Carafe-generated libraries reliably preserve the ``decoy_`` / ``rev_`` prefix on
decoy protein accessions but do not always populate the TSV's ``Decoy`` column
correctly (it stays ``0`` for all rows). DIA-NN and Osprey both want the
``Decoy`` column set to ``1`` for decoy rows, so this script post-processes the
library TSV: for every row whose ``ProteinID`` column has any
semicolon-separated accession starting (case-insensitively) with one of the
configured ``--decoy-prefixes``, the ``Decoy`` column is set to ``1``; rows
where none of the accessions match are set to ``0``.

Streams the file row-by-row so it works on the multi-GB libraries Carafe
produces from large entrapment FASTAs. Output is written to a temp file
next to the destination, then atomically renamed.

Usage::

    python fix_library_decoy_column.py \\
        --input  human-entrapment-carafe_spectral_library_decoy.tsv \\
        --output human-entrapment-carafe_spectral_library_decoy_fixed.tsv \\
        --decoy-prefixes decoy_ rev_ DECOY_
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from pathlib import Path

logger = logging.getLogger(__name__)


def looks_like_decoy(protein_field: str, prefixes_lower: list[str]) -> bool:
    """True iff any semicolon-separated accession starts (case-insensitively)
    with one of the configured prefixes."""
    for acc in protein_field.split(";"):
        acc_lower = acc.lower()
        for p in prefixes_lower:
            if acc_lower.startswith(p):
                return True
    return False


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Patch the Decoy column of a DIA-NN-style library TSV based on "
            "decoy-prefix detection in ProteinID. Use after Carafe to fix "
            "libraries where Carafe left the Decoy column at 0."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, type=Path, help="Input library TSV.")
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output library TSV (written via temp file + atomic rename).",
    )
    parser.add_argument(
        "--decoy-prefixes",
        nargs="+",
        default=["decoy_", "rev_", "DECOY_"],
        help="Protein-accession prefixes that mark decoys (case-insensitive). Default: decoy_ rev_ DECOY_.",
    )
    parser.add_argument(
        "--protein-column",
        default="ProteinID",
        help="Name of the protein column [default: ProteinID].",
    )
    parser.add_argument(
        "--decoy-column",
        default="Decoy",
        help="Name of the Decoy column [default: Decoy].",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    prefixes_lower = [p.lower() for p in args.decoy_prefixes]
    logger.info("Fixing Decoy column in %s -> %s", args.input, args.output)
    logger.info("Decoy prefixes (case-insensitive): %s", prefixes_lower)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Write to a temp file next to the destination, then rename — this avoids
    # leaving a half-written file at the output path on interrupt.
    tmp = tempfile.NamedTemporaryFile(
        mode="w",
        delete=False,
        dir=args.output.parent,
        prefix=args.output.name + ".",
        suffix=".tmp",
    )
    tmp_path = Path(tmp.name)
    try:
        n_rows = 0
        n_decoy = 0
        n_changed = 0
        with args.input.open() as fi:
            header_line = fi.readline()
            tmp.write(header_line)
            cols = header_line.rstrip("\n").split("\t")
            try:
                prot_i = cols.index(args.protein_column)
            except ValueError:
                logger.error(
                    "Input header does not contain column %r. Found columns: %s",
                    args.protein_column,
                    cols,
                )
                return 1
            try:
                decoy_i = cols.index(args.decoy_column)
            except ValueError:
                logger.error(
                    "Input header does not contain column %r. Found columns: %s",
                    args.decoy_column,
                    cols,
                )
                return 1

            for line in fi:
                # Trailing newline handling: split on '\t' is safe on the raw
                # line (the last field still carries the newline); we re-emit
                # the line as-is when nothing changed.
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= max(prot_i, decoy_i):
                    # Malformed row; pass through unchanged.
                    tmp.write(line)
                    n_rows += 1
                    continue
                is_decoy = looks_like_decoy(fields[prot_i], prefixes_lower)
                new_val = "1" if is_decoy else "0"
                if fields[decoy_i] != new_val:
                    n_changed += 1
                    fields[decoy_i] = new_val
                    tmp.write("\t".join(fields) + "\n")
                else:
                    tmp.write(line)
                if is_decoy:
                    n_decoy += 1
                n_rows += 1
                if n_rows % 10_000_000 == 0:
                    logger.info("  processed %d rows", n_rows)
        tmp.close()
        # Atomic rename (within same filesystem only; that's why we created
        # the temp file alongside the destination).
        os.replace(tmp_path, args.output)
    except Exception:
        tmp.close()
        if tmp_path.exists():
            tmp_path.unlink()
        raise

    logger.info("Wrote %d rows; %d decoy rows; %d rows had Decoy column changed", n_rows, n_decoy, n_changed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
