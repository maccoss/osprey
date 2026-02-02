# Decoy Generation

Osprey generates decoy peptides for FDR (False Discovery Rate) control using the target-decoy approach. The decoy generation follows the [pyXcorrDIA](https://github.com/maccoss/pyXcorrDIA) methodology.

## Overview

```
Target peptide: PEPTIDEK
                    ↓
         Enzyme-aware reversal
                    ↓
Decoy peptide:  EDITPEPK
                    ↓
      Fragment m/z recalculation
                    ↓
   Decoy spectrum (b↔y swapped)
```

## Why Decoys?

In proteomics, we need to estimate the false discovery rate (FDR) of our identifications. The target-decoy approach:

1. **Generates decoy sequences** that shouldn't exist in the sample
2. **Searches targets and decoys together** with the same scoring
3. **Uses decoy hits to estimate false positives**

```
FDR ≈ (decoy hits at threshold) / (target hits at threshold)
```

## Enzyme-Aware Sequence Reversal

### The Problem

Simple reversal breaks enzyme cleavage rules:

```
Trypsin cleaves after K/R

Target:  PEPTIDEK    ← ends with K (valid tryptic peptide)
Reverse: KEDIPTED    ← ends with D (NOT a valid tryptic peptide!)
```

### The Solution

Preserve the C-terminal residue (for C-terminal cleaving enzymes):

```
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

```
Original:  P  E  P  T  I  D  E  K
Position:  0  1  2  3  4  5  6  7

Reversed:  E  D  I  T  P  E  P  K
Old pos:   6  5  4  3  2  1  0  7

Position mapping: [6, 5, 4, 3, 2, 1, 0, 7]
```

## Fragment m/z Recalculation

### The Problem

When we reverse a sequence, the fragment ion m/z values change:

```
Target PEPTIDEK:          Decoy EDITPEPK:
  b1 = P                    b1 = E
  b2 = PE                   b2 = ED
  b3 = PEP                  b3 = EDI
  ...                       ...
```

If we don't recalculate, the decoy spectrum won't match any observed peaks!

### The Solution: Ion Type Swapping

When a sequence is reversed:
- **b-ions become y-ions** at complementary positions
- **y-ions become b-ions** at complementary positions

```
For sequence length N:
  b{i} → y{N-i}
  y{i} → b{N-i}
```

### Example

```
Target PEPTIDEK (length 8):
  b3 covers: P-E-P (positions 0-2)

  After reversal to EDITPEPK:
  The same amino acids (P-E-P) are now at positions 4-6
  That's a y4 ion (covers last 4 residues: P-E-P-K)

  b3 → y{8-3} = y5
```

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

```
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

Decoys are marked with:

1. **`is_decoy = true`** flag in the library entry
2. **`DECOY_` prefix** on modified sequence
3. **High bit set** on library ID: `id | 0x80000000`
4. **`DECOY_` prefix** on protein accessions

## Configuration

```yaml
decoy_method: Reverse    # Options: Reverse, Shuffle
decoys_in_library: false # Generate decoys or use existing ones
```

## Implementation

Key files:
- `crates/osprey-scoring/src/lib.rs` - DecoyGenerator, Enzyme enum
- `crates/osprey/src/pipeline.rs` - Decoy generation in pipeline
- `crates/osprey-core/src/config.rs` - DecoyMethod enum

## Target-Decoy Competition

After searching both targets and decoys:

```
Score threshold = 0.5

Targets above threshold: 1000
Decoys above threshold:  10

Estimated FDR = 10 / 1000 = 1%
```

The decoy hits represent false positives, allowing FDR estimation without ground truth.

## Best Practices

1. **Always use enzyme-aware reversal** - Maintains digestion specificity
2. **Recalculate fragment masses** - Essential for spectral matching
3. **Don't reuse decoys across experiments** - Generate fresh each time
4. **Check target/decoy ratio** - Should be ~1:1 in random data
