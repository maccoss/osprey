# Decoy Generation

Osprey generates decoy peptides for FDR (False Discovery Rate) control using the target-decoy approach. The decoy generation follows the [pyXcorrDIA](https://github.com/maccoss/pyXcorrDIA) methodology.

## Overview

```text
Target peptide: PEPTIDEK
                    ↓
         Enzyme-aware reversal
                    ↓
Decoy peptide:  EDITPEPK
                    ↓
      Fragment m/z recalculation
                    ↓
   Decoy spectrum (same ion types, recomputed m/z)
```

## Why Decoys?

In proteomics, we need to estimate the false discovery rate (FDR) of our identifications. The target-decoy approach:

1. **Generates decoy sequences** that shouldn't exist in the sample
2. **Searches targets and decoys together** with the same scoring
3. **Uses decoy hits to estimate false positives**

```text
FDR ≈ (decoy hits at threshold) / (target hits at threshold)
```

## Enzyme-Aware Sequence Reversal

### The Problem

Simple reversal breaks enzyme cleavage rules:

```text
Trypsin cleaves after K/R

Target:  PEPTIDEK    ← ends with K (valid tryptic peptide)
Reverse: KEDIPTED    ← ends with D (NOT a valid tryptic peptide!)
```

### The Solution

Preserve the C-terminal residue (for C-terminal cleaving enzymes):

```text
Target:  PEPTIDEK
         ├──────┤ reverse this part
Decoy:   EDITPEPK   ← K stays at C-terminus
```

### Enzyme Rules

| Enzyme | Cleavage | Preserved Terminus | Example |
|--------|----------|-------------------|---------|
| Trypsin | After K/R | C-terminus | PEPTIDEK → EDITPEPK |
| Lys-C | After K | C-terminus | PEPTIDEK → EDITPEPK |
| Lys-N | Before K | N-terminus | KPEPTIDE → KEDITPEP |
| Asp-N | Before D | N-terminus | DPEPTIDE → DEDITPEP |

### Algorithm

```python
def reverse_sequence(sequence, enzyme):
    if enzyme.preserves_c_terminus():  # Trypsin, Lys-C
        # Reverse everything except C-terminal residue
        prefix = sequence[:-1]
        suffix = sequence[-1]
        return prefix[::-1] + suffix
    else:  # Lys-N, Asp-N
        # Reverse everything except N-terminal residue
        prefix = sequence[0]
        suffix = sequence[1:]
        return prefix + suffix[::-1]
```

### Position Mapping

Track where each residue moves for modification remapping:

```text
Original:  P  E  P  T  I  D  E  K
Position:  0  1  2  3  4  5  6  7

Reversed:  E  D  I  T  P  E  P  K
Old pos:   6  5  4  3  2  1  0  7

Position mapping: [6, 5, 4, 3, 2, 1, 0, 7]
```

## Fragment m/z Recalculation

### The Problem

When we reverse a sequence, the fragment ion m/z values change:

```text
Target PEPTIDEK:          Decoy EDITPEPK:
  b1 = P                    b1 = E
  b2 = PE                   b2 = ED
  b3 = PEP                  b3 = EDI
  ...                       ...
```

If we don't recalculate, the decoy spectrum won't match any observed peaks!

### The Solution: Recompute the m/z, Keep the Ion

Each decoy fragment keeps its target's ion type and ordinal; only the m/z is recomputed
for the permuted sequence. A target y7 yields a decoy y7, so the copied relative
intensity stays on the same-numbered ion of the same series.

```text
For sequence length N:
  b{i} → b{i}, m/z recomputed from the decoy sequence
  y{i} → y{i}, m/z recomputed from the decoy sequence
```

Osprey used to relabel instead - b{i} → y{N-i} and y{i} → b{N-i} - carrying the
intensity along with the relabel. The residue-coverage reasoning behind that mapping was
sound, but fragment intensity is dominated by ion TYPE, not by which residues an ion
spans. y ions are systematically more intense than b ions, so the swap inverted the decoy
spectrum's intensity structure relative to any real peptide. The decoys were easy to
beat, the target-decoy null was not a null, and the reported q-values were far too
optimistic. Measured with entrapment (pass 1, experiment-wide), true FDP at a claimed
1% q:

```text
                with the swap   without
  Stellar       10.9%           1.5%
  Astral         7.6%           2.0%     (library-decoy reference: 1.9%)
```

Skyline, OpenSWATH, DIA-NN, EncyclopeDIA and SpectraST all map the intensity to the same
ion; none of them swaps.

### Rejecting a candidate that is too similar to its target

A permutation that merely differs from the target is not automatically a usable decoy. If
the candidate's theoretical b/y ladder nearly coincides with the target's, the decoy
cannot lose the target/decoy competition on fragment evidence, so it is not an honest
null. `is_candidate_acceptable` rejects a candidate whose ladder overlaps its target's by
more than **0.4** of the candidate's rungs, within a fixed **0.02 Da** window, and the
cycling fallback then supplies another candidate.

The 0.4 threshold is EncyclopeDIA's (`PeptideUtils.getSmartDecoy` rejects above it and
reshuffles). The window is deliberately NOT the run's fragment tolerance: the decoy set
must be a pure function of the library, and keying it to the search tolerance would make
the same library produce different decoys under unit vs HRAM resolution.

At library scale the effect is nil - on the order of 1e-4 of peptides excluded, entrapment
FDP unchanged within noise. It is kept for robustness at SMALL library scale, where
palindromes, low-complexity runs, and isobaric I/L permutations are a far larger fraction.

### Example: Target and Decoy Library Spectra

**Target: PEPTIDEK (charge 2+, precursor m/z = 464.73)**

```text
Library spectrum for PEPTIDEK:

  m/z        intensity   annotation
  ─────────────────────────────────
  324.16     50          b3  (PEP)
  391.18     30          y3  (DEK)
  425.20     35          b4  (PEPT)
  504.27     40          y4  (IDEK)
  538.29     25          b5  (PEPTI)
  605.31     60          y5  (TIDEK)
  653.31     20          b6  (PEPTID)
  702.37     80          y6  (PTIDEK)
  782.36     15          b7  (PEPTIDE)
  831.41     100         y7  (EPTIDEK)
```

**Decoy: EDITPEPK (charge 2+, precursor m/z = 464.73)**

After enzyme-aware reversal (reverse positions 0-6, keep K at C-terminus):

```text
Library spectrum for EDITPEPK:

  m/z        intensity   annotation
  ─────────────────────────────────
  358.16     50          b3  (EDI)
  373.21     30          y3  (EPK)
  459.21     35          b4  (EDIT)
  470.26     40          y4  (PEPK)
  556.26     25          b5  (EDITP)
  571.31     60          y5  (TPEPK)
  685.30     20          b6  (EDITPE)
  684.39     80          y6  (ITPEPK)
  782.36     15          b7  (EDITPEP)
  799.42     100         y7  (DITPEPK)
```

Note: Same precursor m/z (464.73), same number of fragments, same intensities at corresponding positions, but completely different fragment m/z values.

### Mass Calculation

Fragment masses are calculated from first principles:

```python
def calculate_fragment_mz(ion_type, ordinal, charge, sequence, mod_masses):
    PROTON = 1.007276
    H2O = 18.010565

    if ion_type == 'b':
        # b-ion: N-terminal fragment
        residues = sequence[0:ordinal]
    else:  # y-ion
        # y-ion: C-terminal fragment
        residues = sequence[len(sequence)-ordinal:]

    mass = sum(AA_MASSES[aa] for aa in residues)
    mass += sum(mod_masses.get(i, 0) for i in residue_positions)

    if ion_type == 'b':
        mass += PROTON  # [M+H]+ minus C-terminal OH
    else:  # y-ion
        mass += H2O + PROTON  # [M+H]+ plus H2O

    mz = (mass + (charge - 1) * PROTON) / charge
    return mz
```

## Modification Handling

Modifications must be remapped to new positions:

```text
Target: PEPTIDEK with Oxidation at position 3 (on T)
        P  E  P  T  I  D  E  K
        0  1  2  3* 4  5  6  7

Decoy:  EDITPEPK
        E  D  I  T  P  E  P  K
        0  1  2  3* 4  5  6  7

Position mapping: [6, 5, 4, 3, 2, 1, 0, 7]
Old position 3 → New position 3 (symmetric in this case)
```

### Algorithm

```python
def remap_modifications(modifications, position_mapping):
    # Create reverse map: old_pos → new_pos
    reverse_map = {old: new for new, old in enumerate(position_mapping)}

    new_mods = []
    for mod in modifications:
        if mod.position in reverse_map:
            new_mod = mod.copy()
            new_mod.position = reverse_map[mod.position]
            new_mods.append(new_mod)

    return new_mods
```

## Precursor Mass Conservation

The decoy has the **same precursor m/z** as the target because:

- Same amino acid composition
- Same modifications (just at different positions)
- Same charge state

This ensures decoys compete fairly with targets at the precursor level.

## Decoy Identification

Whether generated by Osprey or supplied by the library, decoys are marked with:

1. **`is_decoy = true`** flag in the library entry
2. **High bit set** on library ID: `id | 0x80000000` (the constant `DECOY_ID_BIT`)
3. **`DECOY_` prefix** on protein accessions (generated decoys; library decoys retain their original prefix)
4. **`DECOY_` prefix** on modified sequence (generated decoys only)

The high bit on `id` lets downstream code recover the paired target via
`base_id = entry_id & 0x7FFFFFFF`. For Osprey-generated decoys this pairing
is exact by construction (the decoy inherits the target's base_id). For
library-supplied decoys Osprey establishes the same pairing in a post-load
step using either a FDRBench manifest or amino-acid composition matching;
see "Library-supplied decoys" below.

## Configuration

```yaml
decoy_method: Reverse       # Options: Reverse, Shuffle, FromLibrary
decoys_in_library: false    # Set true to trust decoys already in the library
decoy_prefixes:             # Protein-accession prefixes that mark library decoys
  - DECOY_                  # (case-insensitive). Default list covers Osprey's
  - rev_                    # own convention plus DIA-NN / EncyclopeDIA / Carafe.
  - decoy_
```

`decoys_in_library: true` and `decoy_method: FromLibrary` are equivalent;
both skip `DecoyGenerator` and instead scan the library for decoys flagged
by prefix.

## Library-supplied decoys

When `decoys_in_library` is enabled, Osprey processes the library in
three post-load steps: decoy marking, target-decoy pairing, and
(when a pairing manifest is supplied) authoritative protein-ID
substitution.

### Step 1: decoy marking

Decoys are detected via two complementary signals, OR'd together:

- **TSV `Decoy` column** (DIA-NN convention, `0` / `1` / case-insensitive
  `true`/`yes`/`y`/`t`). The DIA-NN TSV loader reads it at load time
  and sets `is_decoy = true` accordingly.
- **Protein-accession prefix scan**: each entry's `protein_ids` are
  scanned, and any entry whose protein accession begins
  (case-insensitively) with one of the configured `decoy_prefixes` is
  flagged. Marking is idempotent; entries already flagged by the
  loader just get `DECOY_ID_BIT` canonicalised on their `id`.

The pipeline log breaks the count down: `decoys flagged via TSV Decoy
column = N, via prefix scan = M`. Libraries that set ONLY the column,
ONLY the prefix, or BOTH all work.

### Step 2: target-decoy pairing (hybrid manifest + composition fallback)

Marking alone is not enough — the SVM, LDA calibration, and
cross-validation fold splitting all rely on each decoy sharing a
`base_id` with its target (linked by
`base_id = entry_id & 0x7FFFFFFF`). Without pairing,
`compete_from_indices` treats every entry as an "unpaired winner" and
the SVM trains on the full target population without quality
filtering, producing FDR estimates that are too optimistic. Osprey
runs pairing in two stages, hybrid by design:

- **Stage 2a: manifest-based** (when supplied via
  `decoy_pairing_manifest`): reads a FDRBench-style 5-column TSV
  (`sequence`, `decoy`, `proteins`, `peptide_type`,
  `peptide_pair_index`) and uses each row's `peptide_pair_index` to
  group entries into target/decoy pairs. Two pairs come from each
  group: `(target, decoy)` and `(p_target, p_decoy)`. Lookup is by
  unmodified peptide sequence; charge is honored so charge-2 decoys
  pair with charge-2 targets. **The manifest is authoritative for
  decoy classification too**: any library entry whose sequence is in
  the manifest as `decoy` or `p_decoy` is flagged `is_decoy = true`
  even if the prefix scan missed it (which happens when a library
  predictor strips the `decoy_` / `rev_` prefix from protein
  accessions during processing — the Carafe failure mode).
- **Stage 2b: composition-based fallback** (always runs): for each
  decoy not paired by the manifest, Osprey strips a configured prefix
  from each protein accession to recover the target accession, then
  looks for a target peptide with the same
  `(stripped_accession, charge, sorted_amino_acid_composition)`.
  Within a composition group on a single protein, target and decoy
  entries are sorted by `(sequence, id)` and zipped 1:1, so pairings
  are deterministic regardless of input order. The fallback works for
  any decoy strategy that preserves AA composition (Carafe's
  randomized decoys, sequence-reversal decoys).

Why hybrid? In practice, a manifest generated alongside a library can
still miss most of the library's peptides if the library generator
(e.g. Carafe) uses different digestion rules than FDRBench. The
composition fallback recovers those misses. When both stages run on
real data, the breakdown looks something like:

```text
Library-decoy pairing: paired 10103456/10175597 decoys (99.3%) -
3027526 via manifest, 7075930 via composition; 72141 unpaired decoys,
72141 unpaired targets
```

Both stages write decoy IDs as `target_id | DECOY_ID_BIT`. Osprey
reports the pairing breakdown at INFO level and **bails with a hard
error if fewer than `decoy_pair_min_fraction` (default 80%) of decoys
pair successfully overall**. Below that threshold, FDR estimates would
be unreliable; the user must either supply a more complete manifest,
fix the library's accession conventions, or unset `decoys_in_library`
and let Osprey generate decoys.

### Step 3: protein-ID substitution from the manifest (when supplied)

When the manifest's `proteins` column is non-empty for a matched
library entry, `DecoyPairingManifest::apply_to_library` REPLACES the
library entry's `protein_ids` with the manifest's clean source-protein
list. The pipeline logs `manifest replaced protein_ids on N library
entries`.

This solves two problems at once:

- **Library predictors that deduplicate by FASTA accession**. When
  generating an entrapment FASTA via
  `scripts/build_entrapment_peptide_fasta.py` for a predictor like
  Carafe, the script appends a per-peptide counter
  (`_pep00001`, `_pep00002`, …) to each accession so every FASTA entry
  is unique to the predictor. Carafe then writes that suffixed
  accession into the library's `ProteinID` column. The manifest, by
  contrast, carries the original clean accession; substituting it
  back at load time restores correct protein parsimony / picked-protein
  FDR rollup.
- **Shared peptides**. The FASTA emits the joined source-protein list
  in each header (e.g.
  `>sp|P12345_pep00001|GENE_A;sp|Q67890_pep00001|GENE_B`), matching the
  semicolon convention DIA-NN / Carafe already use in their library
  `ProteinID` column for shared peptides. The library predictor
  typically propagates this joined string verbatim into `ProteinID`,
  so the multi-protein attribution survives even before the manifest
  substitution. When the predictor doesn't preserve the joined header,
  the manifest's `proteins` column carries the same information and
  the substitution restores `entry.protein_ids = ["sp|A|...", "sp|B|..."]`
  at load time.

The substitution is a no-op when the manifest's `proteins` column is
empty or `-`. Library entries whose sequence isn't in the manifest are
left untouched.

Prior to this support, `decoys_in_library: true` skipped DecoyGenerator
but no code set `is_decoy` and no pairing was performed, so real decoys
were silently treated as targets and FDR estimates were corrupted; see
`release-notes/RELEASE_NOTES_next.md` (or current draft).

## Implementation

Key files:

- `crates/osprey-scoring/src/lib.rs` - `DecoyGenerator`, `Enzyme` enum
- `crates/osprey/src/pipeline.rs` - decoy generation / library marking dispatch
- `crates/osprey-core/src/config.rs` - `DecoyMethod` enum, `decoy_prefixes`
- `crates/osprey-core/src/types.rs` - `LibraryEntry::looks_like_library_decoy()`,
  `apply_library_decoy_marking()`, `DECOY_ID_BIT`

## Target-Decoy Competition

After searching both targets and decoys:

```text
Score threshold = 0.5

Targets above threshold: 1000
Decoys above threshold:  10

Estimated FDR = 10 / 1000 = 1%
```

The decoy hits represent false positives, allowing FDR estimation without ground truth.

## Collision Detection (pyXcorrDIA Approach)

### The Problem

When generating decoys by reversal, collisions can occur:

```text
Target database contains:
  PEPTIDEK  →  reverses to  →  EDITPEPK
  EDITPEPK  ←  already exists as a target!
```

If the decoy sequence matches an existing target:

- The decoy will score identically to the real target (same fragments!)
- This causes spurious high-scoring decoy wins
- FDR calculation becomes unreliable (inflated at high scores)

### The Solution: Multi-Strategy Decoy Generation

Following pyXcorrDIA, Osprey uses a collision-aware approach:

```text
┌─────────────────────────────────────────────────────────┐
│  For each target peptide:                               │
│                                                         │
│  1. Try REVERSAL (default)                              │
│     ↓                                                   │
│     Check: Is reversed sequence != target?              │
│            Is reversed sequence NOT in target database? │
│     ↓                                                   │
│     Yes → Use reversed decoy ✓                          │
│     No  → Continue to step 2                            │
│                                                         │
│  2. Try CYCLING (fallback)                              │
│     ↓                                                   │
│     Cycle sequence by 1, 2, 3... positions              │
│     Check for uniqueness at each cycle length           │
│     ↓                                                   │
│     Found unique → Use cycled decoy ✓                   │
│     All collide  → Continue to step 3                   │
│                                                         │
│  3. EXCLUDE peptide                                     │
│     ↓                                                   │
│     No valid decoy possible                             │
│     Remove both target and decoy from analysis          │
└─────────────────────────────────────────────────────────┘
```

### Cycling Method

When reversal collides, cycling shifts the sequence:

```text
Original:  PEPTIDEK
Cycle 1:   EPTIDEPK  (shift by 1, keep C-term K)
Cycle 2:   PTIDEPEK  (shift by 2, keep C-term K)
Cycle 3:   TIDEPEPK  (shift by 3, keep C-term K)
...
```

The first unique cycled sequence is used as the decoy.

### Statistics Tracking

Osprey reports decoy generation statistics:

```text
Decoy generation statistics:
  Reversed: 1,523,456 (97.3%)           ← Primary method succeeded
  Cycling fallback: 41,234 (2.6%)       ← Used cycling
  Excluded (no unique decoy): 983 (0.1%) ← Could not generate valid decoy
```

This visibility helps detect potential issues:

- **High cycling rate** → Many palindromic sequences
- **High exclusion rate** → Small/redundant database
- **Zero exclusions** → Normal for typical proteome databases

### Why This Matters for FDR

Without collision detection:

```text
Score ranking:        With collisions:
1. PEPTIDEK (target)  ← True positive
2. EDITPEPK (decoy)   ← But this IS a real peptide!
3. ...                   Scores as well as target
                         FDR appears high even at top scores
```

With collision detection:

```text
Score ranking:        After excluding collisions:
1. PEPTIDEK (target)  ← True positive
2. ANOTHERR (target)  ← True positive
3. ...                   Decoys only score high when random match
                         FDR is accurate
```

## Best Practices

1. **Always use enzyme-aware reversal** - Maintains digestion specificity
2. **Recalculate fragment masses** - Essential for spectral matching
3. **Check collision statistics** - High exclusion rates may indicate issues
4. **Don't reuse decoys across experiments** - Generate fresh each time
5. **Check target/decoy ratio** - Should be ~1:1 in random data
