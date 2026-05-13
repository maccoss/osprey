"""Build a peptide-level FASTA + FDRBench-style pairing manifest from a uniprot FASTA.

Input is a standard uniprot-format FASTA of target proteins. For each tryptic
peptide the script can optionally also emit:

- an entrapment (``p_target``) peptide: deterministic shuffle of the target
  with the C-terminal residue preserved (so the peptide stays tryptic), using
  a per-peptide RNG seed derived from the peptide sequence + a master
  entrapment seed;
- a decoy peptide: same shuffling rule but with a different master seed;
- a ``p_decoy`` peptide: a shuffle of the entrapment peptide using the decoy
  master seed.

If any synthetic peptide (p_target, decoy, p_decoy) happens to match a real
target peptide somewhere in the input, the whole quartet (target, p_target,
decoy, p_decoy) is dropped — matches FDRBench's collision-drop policy.

Output FASTA layout: each peptide gets its own FASTA entry whose accession is
the **source protein accession** (no per-peptide suffix). Multiple peptides
from the same source protein produce multiple FASTA entries that share the
same header. DIA-NN / Carafe digest each entry independently (each one is one
tryptic peptide already), so the resulting library has the source protein
accession in its ProteinID column — Osprey's protein parsimony and
picked-protein FDR work as on any normal library.

Header conventions::

    target:    >sp|<acc>|<entry>
    p_target:  >sp|<acc>_p_target|<entry>_p_target
    decoy:     >rev_sp|<acc>|<entry>
    p_decoy:   >rev_sp|<acc>_p_target|<entry>_p_target

The manifest is the standard FDRBench 5-column TSV
(sequence, decoy, proteins, peptide_type, peptide_pair_index) that Osprey's
``--decoy-pairing-manifest`` reads.

Usage::

    python build_entrapment_peptide_fasta.py \\
        --input uniprot_human.fasta \\
        --output peptide_lib.fasta \\
        --manifest pairing.tsv \\
        --add-entrapment --add-decoys
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import random
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# Standard tryptic cleavage: after K or R, unless the next residue is P.
TRYPSIN_RE = re.compile(r"(?<=[KR])(?!P)")

# Header marker suffix for entrapment proteins (matches FDRBench convention).
DEFAULT_ENTRAPMENT_SUFFIX = "_p_target"
# Prefix on the FASTA header that flags decoy entries. We use `decoy_` (rather
# than the also-common `rev_`) because downstream tooling — DIA-NN /
# Carafe — does not always set the TSV `Decoy` column correctly, and
# `decoy_` is the more recognisable convention for a post-process step that
# patches the column based on protein-accession prefix detection.
DEFAULT_DECOY_PREFIX = "decoy_"

# Monoisotopic residue masses (Da). U/O selenocysteine/pyrrolysine left out; B,
# Z, J, X are ambiguity codes — peptides containing them are skipped silently.
AA_MONO_MASS: dict[str, float] = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406,
    "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
    "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841,
    "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
H2O_MONO = 18.01056
PROTON_MONO = 1.007276


def peptide_neutral_mass(seq: str) -> float | None:
    """Monoisotopic neutral mass of a peptide, or `None` if any residue is not
    in `AA_MONO_MASS` (e.g. contains B, Z, J, X, U, O)."""
    total = H2O_MONO
    for aa in seq:
        m = AA_MONO_MASS.get(aa)
        if m is None:
            return None
        total += m
    return total


def fits_mz_range(neutral_mass: float, charges: list[int], min_mz: float, max_mz: float) -> bool:
    """True iff at least one of the allowed charge states produces an m/z
    inside `[min_mz, max_mz]`."""
    for z in charges:
        mz = (neutral_mass + z * PROTON_MONO) / z
        if min_mz <= mz <= max_mz:
            return True
    return False


@dataclass
class ProteinRecord:
    """A single source protein parsed from the input FASTA."""

    accession: str  # e.g. "P12345"
    entry_name: str  # e.g. "HUMAN_GENE"
    db: str  # "sp" or "tr"
    description: str  # everything after the accession block on the header line
    sequence: str


@dataclass
class Quartet:
    """The four sequences derived from a single target peptide."""

    target: str
    p_target: str | None = None
    decoy: str | None = None
    p_decoy: str | None = None
    # Source proteins this target peptide was tryptic-digested from.
    # If the peptide is shared across multiple proteins it appears once here
    # per source. ``ProteinRecord`` instances are referenced, not copied.
    sources: list[ProteinRecord] = field(default_factory=list)


def parse_fasta(path: Path):
    """Yield ``ProteinRecord`` for each entry in a uniprot-format FASTA.

    Recognises ``>db|accession|entry_name description`` headers and
    accepts headers that don't follow that pattern (passed through as
    accession-only). Strips any trailing ``*`` from the sequence.
    """
    header: str | None = None
    chunks: list[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
            return None
        seq = "".join(chunks).rstrip("*").upper()
        if not seq:
            return None
        # Parse "db|accession|entry_name description"
        parts = header.split(None, 1)
        ident = parts[0]
        description = parts[1] if len(parts) > 1 else ""
        ident_parts = ident.split("|")
        if len(ident_parts) >= 3:
            db, acc, entry_name = ident_parts[0], ident_parts[1], ident_parts[2]
        elif len(ident_parts) == 2:
            db, acc, entry_name = ident_parts[0], ident_parts[1], ident_parts[1]
        else:
            db, acc, entry_name = "sp", ident, ident
        return ProteinRecord(
            accession=acc,
            entry_name=entry_name,
            db=db,
            description=description,
            sequence=seq,
        )

    with path.open() as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                rec = flush()
                if rec is not None:
                    yield rec
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        rec = flush()
        if rec is not None:
            yield rec


def tryptic_digest(seq: str, missed_cleavages: int, min_len: int, max_len: int) -> list[str]:
    """Return all tryptic peptides of length in [min_len, max_len] with up to
    `missed_cleavages` missed cleavages.

    Cleavage rule: after K or R, unless the next residue is P (KP / RP).
    Duplicates (same sequence from different cleavage windows) are dropped.
    """
    if not seq:
        return []
    fragments = TRYPSIN_RE.split(seq)
    # Drop empty fragments at the boundary (can happen if the sequence starts
    # with a cleavage site).
    fragments = [f for f in fragments if f]
    peptides: set[str] = set()
    n = len(fragments)
    for i in range(n):
        for k in range(missed_cleavages + 1):
            j = i + k + 1
            if j > n:
                break
            pep = "".join(fragments[i:j])
            if min_len <= len(pep) <= max_len:
                peptides.add(pep)
    return sorted(peptides)


def shuffle_preserving_cterm(seq: str, master_seed: int) -> str:
    """Deterministically shuffle all but the last residue.

    The per-peptide RNG seed is derived as SHA-1(master_seed:seq) so the same
    (sequence, master_seed) pair always returns the same shuffle, while
    different sequences get independent shuffles. Length-1 and length-2
    inputs are returned unchanged.
    """
    if len(seq) <= 2:
        return seq
    last = seq[-1]
    body = list(seq[:-1])
    digest = hashlib.sha1(f"{master_seed}:{seq}".encode()).digest()
    pep_seed = int.from_bytes(digest[:8], "big", signed=False)
    rng = random.Random(pep_seed)
    rng.shuffle(body)
    return "".join(body) + last


def build_full_protein_label(rec: ProteinRecord, p_target: bool, decoy: bool, entrap_suffix: str, decoy_prefix: str) -> str:
    """Build the ``db|accession|entry_name`` label for a peptide kind."""
    acc = rec.accession + entrap_suffix if p_target else rec.accession
    entry = rec.entry_name + entrap_suffix if p_target else rec.entry_name
    label = f"{rec.db}|{acc}|{entry}"
    if decoy:
        label = decoy_prefix + label
    return label


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a peptide-level entrapment FASTA + pairing manifest from a uniprot FASTA.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, type=Path, help="Source uniprot FASTA (targets only).")
    parser.add_argument("--output", required=True, type=Path, help="Output peptide-level FASTA.")
    parser.add_argument("--manifest", type=Path, default=None, help="Output FDRBench-style pairing manifest TSV (optional).")
    parser.add_argument("--add-entrapment", action="store_true", help="Include p_target entrapment peptides.")
    parser.add_argument("--add-decoys", action="store_true", help="Include decoy peptides (and p_decoys when entrapment is also enabled).")
    parser.add_argument(
        "--missed-cleavages",
        type=int,
        default=1,
        choices=range(0, 6),
        metavar="{0..5}",
        help="Max missed cleavages [default: 1]. A value of N includes peptides with 0..N missed cleavages.",
    )
    parser.add_argument("--min-length", type=int, default=7, help="Minimum peptide length [default: 7].")
    parser.add_argument("--max-length", type=int, default=35, help="Maximum peptide length [default: 35].")
    parser.add_argument(
        "--min-mz",
        type=float,
        default=400.0,
        help="Minimum precursor m/z for at least one allowed charge state [default: 400].",
    )
    parser.add_argument(
        "--max-mz",
        type=float,
        default=900.0,
        help="Maximum precursor m/z for at least one allowed charge state [default: 900].",
    )
    parser.add_argument(
        "--charges",
        type=int,
        nargs="+",
        default=[2, 3],
        help="Allowed precursor charge states; a peptide is kept if any of these charges produces an m/z in [min_mz, max_mz] [default: 2 3].",
    )
    parser.add_argument("--entrapment-seed", type=int, default=42, help="Master RNG seed for entrapment shuffling [default: 42].")
    parser.add_argument("--decoy-seed", type=int, default=24, help="Master RNG seed for decoy shuffling [default: 24].")
    parser.add_argument(
        "--entrapment-suffix",
        type=str,
        default=DEFAULT_ENTRAPMENT_SUFFIX,
        help=f"Suffix appended to accession/entry-name for entrapment proteins [default: {DEFAULT_ENTRAPMENT_SUFFIX!r}].",
    )
    parser.add_argument(
        "--decoy-prefix",
        type=str,
        default=DEFAULT_DECOY_PREFIX,
        help=f"Prefix prepended to db|accession|entry for decoy headers [default: {DEFAULT_DECOY_PREFIX!r}].",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging.")

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    # Step 1: parse the FASTA and tryptic-digest every protein. Collect all
    # unique target peptides, each with a list of source ProteinRecord refs.
    # Apply length, residue (skip unknown AAs), and precursor m/z filters
    # at digest time so we don't carry unusable peptides through the
    # shuffle and collision-check stages.
    logger.info("Reading FASTA: %s", args.input)
    if args.min_length < 1 or args.max_length < args.min_length:
        logger.error("invalid length range: %d-%d", args.min_length, args.max_length)
        return 1
    if args.min_mz >= args.max_mz:
        logger.error("invalid m/z range: %.3f-%.3f", args.min_mz, args.max_mz)
        return 1
    if any(z < 1 or z > 10 for z in args.charges):
        logger.error("invalid charge state in %s (allowed: 1..10)", args.charges)
        return 1

    target_to_sources: dict[str, list[ProteinRecord]] = {}
    n_proteins = 0
    n_dropped_unknown_aa = 0
    n_dropped_out_of_mz = 0
    for rec in parse_fasta(args.input):
        n_proteins += 1
        for pep in tryptic_digest(
            rec.sequence,
            missed_cleavages=args.missed_cleavages,
            min_len=args.min_length,
            max_len=args.max_length,
        ):
            if pep in target_to_sources:
                target_to_sources[pep].append(rec)
                continue
            neutral_mass = peptide_neutral_mass(pep)
            if neutral_mass is None:
                n_dropped_unknown_aa += 1
                continue
            if not fits_mz_range(neutral_mass, args.charges, args.min_mz, args.max_mz):
                n_dropped_out_of_mz += 1
                continue
            target_to_sources[pep] = [rec]
    logger.info(
        "Digested %d proteins (mc=%d, len=%d-%d, charges=%s, m/z=%.1f-%.1f): "
        "%d unique target peptides retained, %d dropped (unknown AA), %d dropped (out of m/z range)",
        n_proteins,
        args.missed_cleavages,
        args.min_length,
        args.max_length,
        args.charges,
        args.min_mz,
        args.max_mz,
        len(target_to_sources),
        n_dropped_unknown_aa,
        n_dropped_out_of_mz,
    )

    # Step 2: build quartets. p_target/decoy/p_decoy are deterministic
    # shuffles of the target (or entrapment) sequence.
    quartets: list[Quartet] = []
    for pep, sources in target_to_sources.items():
        q = Quartet(target=pep, sources=sources)
        if args.add_entrapment:
            q.p_target = shuffle_preserving_cterm(pep, args.entrapment_seed)
        if args.add_decoys:
            q.decoy = shuffle_preserving_cterm(pep, args.decoy_seed)
            if args.add_entrapment and q.p_target is not None:
                q.p_decoy = shuffle_preserving_cterm(q.p_target, args.decoy_seed)
        quartets.append(q)
    logger.info("Built %d quartets", len(quartets))

    # Step 3: collision check. Drop any quartet where a synthetic peptide
    # (p_target / decoy / p_decoy) collides with a real target sequence.
    target_set = set(target_to_sources.keys())
    kept: list[Quartet] = []
    n_dropped_p_target = 0
    n_dropped_decoy = 0
    n_dropped_p_decoy = 0
    for q in quartets:
        drop = False
        if q.p_target is not None and q.p_target in target_set:
            n_dropped_p_target += 1
            drop = True
        if q.decoy is not None and q.decoy in target_set:
            n_dropped_decoy += 1
            drop = True
        if q.p_decoy is not None and q.p_decoy in target_set:
            n_dropped_p_decoy += 1
            drop = True
        if not drop:
            kept.append(q)
    n_dropped_total = len(quartets) - len(kept)
    logger.info(
        "Collision-drop pass: %d quartets dropped (%d via p_target match, %d via decoy, %d via p_decoy)",
        n_dropped_total,
        n_dropped_p_target,
        n_dropped_decoy,
        n_dropped_p_decoy,
    )
    logger.info("%d quartets retained", len(kept))

    # Stable pair_index assignment by sorted target sequence so the manifest
    # is reproducible regardless of dict iteration order.
    kept.sort(key=lambda q: q.target)

    # Step 4: write FASTA. Each peptide gets one entry per source protein
    # (shared peptides produce multiple entries with different accessions
    # but the same peptide sequence). DIA-NN / Carafe will merge by sequence.
    logger.info("Writing FASTA: %s", args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    n_target_entries = 0
    n_ptarget_entries = 0
    n_decoy_entries = 0
    n_pdecoy_entries = 0
    with args.output.open("w") as fo:
        for q in kept:
            for src in q.sources:
                # target
                label = build_full_protein_label(
                    src, p_target=False, decoy=False,
                    entrap_suffix=args.entrapment_suffix,
                    decoy_prefix=args.decoy_prefix,
                )
                fo.write(f">{label}\n{q.target}\n")
                n_target_entries += 1
                # p_target (entrapment)
                if q.p_target is not None:
                    label = build_full_protein_label(
                        src, p_target=True, decoy=False,
                        entrap_suffix=args.entrapment_suffix,
                        decoy_prefix=args.decoy_prefix,
                    )
                    fo.write(f">{label}\n{q.p_target}\n")
                    n_ptarget_entries += 1
                # decoy
                if q.decoy is not None:
                    label = build_full_protein_label(
                        src, p_target=False, decoy=True,
                        entrap_suffix=args.entrapment_suffix,
                        decoy_prefix=args.decoy_prefix,
                    )
                    fo.write(f">{label}\n{q.decoy}\n")
                    n_decoy_entries += 1
                # p_decoy
                if q.p_decoy is not None:
                    label = build_full_protein_label(
                        src, p_target=True, decoy=True,
                        entrap_suffix=args.entrapment_suffix,
                        decoy_prefix=args.decoy_prefix,
                    )
                    fo.write(f">{label}\n{q.p_decoy}\n")
                    n_pdecoy_entries += 1
    logger.info(
        "Wrote FASTA: %d target, %d p_target, %d decoy, %d p_decoy entries (one per peptide-source pair)",
        n_target_entries,
        n_ptarget_entries,
        n_decoy_entries,
        n_pdecoy_entries,
    )

    # Step 5: optional manifest.
    if args.manifest is not None:
        logger.info("Writing manifest: %s", args.manifest)
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        with args.manifest.open("w") as fm:
            fm.write("sequence\tdecoy\tproteins\tpeptide_type\tpeptide_pair_index\n")
            for pair_idx, q in enumerate(kept):
                # The "proteins" column joins all source proteins with ";"
                # following FDRBench/DIA-NN convention. Each peptide_type
                # in the same pair_index lists the SAME set of source
                # proteins (with the appropriate p_target/rev_ wrapping).
                def proteins_for(p_target: bool, decoy: bool) -> str:
                    return ";".join(
                        build_full_protein_label(
                            src, p_target=p_target, decoy=decoy,
                            entrap_suffix=args.entrapment_suffix,
                            decoy_prefix=args.decoy_prefix,
                        )
                        for src in q.sources
                    )

                fm.write(f"{q.target}\tNo\t{proteins_for(False, False)}\ttarget\t{pair_idx}\n")
                if q.p_target is not None:
                    fm.write(f"{q.p_target}\tNo\t{proteins_for(True, False)}\tp_target\t{pair_idx}\n")
                if q.decoy is not None:
                    fm.write(f"{q.decoy}\tYes\t{proteins_for(False, True)}\tdecoy\t{pair_idx}\n")
                if q.p_decoy is not None:
                    fm.write(f"{q.p_decoy}\tYes\t{proteins_for(True, True)}\tp_decoy\t{pair_idx}\n")
        logger.info("Wrote manifest with %d pair_index groups", len(kept))

    return 0


if __name__ == "__main__":
    sys.exit(main())
