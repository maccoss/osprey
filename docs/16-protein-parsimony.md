# Protein Parsimony

Protein parsimony resolves the peptide-to-protein mapping from the spectral library into a minimal set of protein groups that can explain the observed peptides. It runs **always** (regardless of whether `--protein-fdr` is enabled) because the resulting peptide-to-group mapping is used by both protein-level FDR (when enabled) and downstream interpretation of results.

## Why Parsimony?

A peptide sequence can map to multiple proteins in the library (shared peptides across paralogs, isoforms, or homologous proteins). A single parsimony pass asks: *given the peptides we observed, what is the smallest set of proteins that can explain them?*

This is the classic **minimum set cover** problem, solved here with a combination of exact rules (identical-set merging, subset elimination) and a greedy heuristic for the remaining ambiguity (razor assignment).

## Algorithm

```
Input: library entries with protein_ids; optionally a detected-peptide set
Output: ProteinParsimonyResult {
    groups: Vec<ProteinGroup>,
    peptide_to_groups: HashMap<peptide, Vec<ProteinGroupId>>,
}

Steps:
  1. Build bipartite graph peptide <-> protein from target library entries.
     If a detected set is provided, exclude peptides that weren't detected.
  2. Group proteins with identical peptide sets (indistinguishable).
  3. Subset elimination: drop groups whose peptides are a strict subset of another.
  4. Classify remaining peptides as unique (1 group) or shared (2+ groups).
  5. Apply SharedPeptideMode (All | Razor | Unique).
```

### Step 1: Bipartite Graph

For each library entry with `protein_ids = [P1, P2, ...]`:
- Add edges peptide ↔ P1, peptide ↔ P2, ...

If a `detected_peptides` set is provided (typically peptides passing peptide-level FDR at experiment scope), only peptides in that set are added to the graph. This ensures protein groups reflect actual observations, not theoretical library content.

### Step 2: Identical-Set Merging

Proteins with identical peptide sets are indistinguishable based on the detected evidence. These are collapsed into a single **protein group** (with multiple `accessions`). For example:

```
P1 peptides: {A, B, C}
P2 peptides: {A, B, C}  ← identical to P1
P3 peptides: {A, B}

Group 1: accessions=[P1, P2], peptides={A, B, C}
Group 2: accessions=[P3],     peptides={A, B}
```

### Step 3: Subset Elimination

A protein group whose peptides are a strict subset of another group can be explained by that larger group — it provides no additional evidence. These subset groups are removed.

```
Group 2 (peptides {A, B}) is a strict subset of Group 1 (peptides {A, B, C})
→ Group 2 is eliminated
```

### Step 4: Unique vs Shared Classification

For each remaining peptide:
- **Unique**: belongs to exactly 1 group (proteotypic)
- **Shared**: belongs to 2+ groups

### Step 5: Shared Peptide Modes

Three modes control what happens to shared peptides:

| Mode | Behavior |
|------|----------|
| `All` (default) | Shared peptides contribute to **all** their groups. Each peptide inherits the best (lowest) protein q-value among its groups. Maximum sensitivity. |
| `Razor` | Shared peptides are assigned exclusively to **one** group each, chosen by iterative greedy set cover (see below). Matches MaxQuant's razor peptide logic. |
| `Unique` | Shared peptides are **excluded** from protein scoring and output entirely. Only unique peptides are used. Most conservative. |

## Razor: Iterative Greedy Set Cover

Razor mode resolves shared peptides to a single protein group each. Osprey uses a textbook iterative greedy algorithm:

```
while any shared peptides remain unassigned:
    find the group G with the MOST unique peptides that still has
        at least one unassigned shared peptide (tiebreak: lowest group ID)
    for each unassigned shared peptide of G (sorted alphabetically):
        assign to G (add to G.unique_peptides, remove from others' shared sets)
        mark as assigned
```

### Why Iterative?

Consider three groups:

```
Group 1: unique={A, B, C}, shared={X, Y}
Group 2: unique={D},        shared={X}
Group 3: unique={E},        shared={Y}
```

Naive single-pass razor might iterate shared peptides in arbitrary order and assign each to "the group with the most unique peptides right now." But the order of iteration matters — assigning X first bumps Group 1's unique count from 3 to 4, which could affect Y's decision.

The iterative greedy looks at the **global state** at every step. In this example:
- **Round 1**: Group 1 has 3 unique, Group 2 and Group 3 each have 1. Group 1 wins. It claims both X and Y in one batch. Result: Group 1 has {A, B, C, X, Y}, Groups 2 and 3 unchanged.

A single-pass approach could have produced a different result (e.g., X → Group 2 if X was processed before Group 1 had been strengthened). The iterative approach is path-independent.

### Determinism

All sources of non-determinism are eliminated:

1. **Collection of shared peptides** at the start is sorted alphabetically
2. **Winning group selection** uses `max_by_key((unique_count, Reverse(group_id)))` — ties broken by lowest group ID
3. **Claimed peptides per round** are sorted alphabetically before processing
4. **peptide_to_groups updates** are deterministic because they follow the sorted claim order

This means running the same parsimony on the same input always produces byte-identical output, regardless of HashMap iteration order or parallelism. See `test_shared_peptides_razor_deterministic` in `crates/osprey-fdr/src/protein.rs`.

## Examples

### Example 1: Simple isoforms

```
P1 peptides: {A, B, C, X}
P2 peptides: {D, X}

Unique: A, B, C → P1
        D       → P2
Shared: X       → P1, P2

All mode:    X contributes to both P1 and P2
Razor mode:  X → P1 (3 unique > 1 unique)
Unique mode: X excluded
```

### Example 2: Cascading razor

```
P1 peptides: {A, B, C, X, Y}
P2 peptides: {D, X, Z}
P3 peptides: {E, Y, Z}

Round 1: P1 has 3 unique (most) → claims X, Y
         P1 = {A, B, C, X, Y}, P2 = {D, Z}, P3 = {E, Z}
Round 2: P2 and P3 tie at 1 unique → P2 wins on lower group ID → claims Z
         P2 = {D, Z}, P3 = {E}
Done.
```

### Example 3: No shared peptides

```
P1 peptides: {A, B}
P2 peptides: {C, D}

All three modes produce identical output (no shared peptides to resolve).
```

## How Parsimony Feeds into Protein FDR

When `--protein-fdr` is set, the parsimony result is used as input to **two-pass picked-protein FDR** (Savitski 2015). The parsimony graph itself is built once per pass (same library, same shared-peptide mode, same result) but the peptide score pool differs between passes.

### First Pass (Pre-Compaction, Gating)

Runs before the compaction step in `pipeline.rs`. Produces `run_protein_qvalue` on each FdrEntry.

1. Build parsimony from peptides passing first-pass peptide FDR.
2. Call `collect_best_peptide_scores()` on the **full pre-compaction** `per_file_entries` — both targets and decoys, regardless of whether their base_id passed first-pass precursor FDR. This gives symmetric target + decoy pools for the picked-protein competition.
3. For each protein group:
   - Compute `target_score = max(SVM score over target peptides passing gate)`.
   - Compute `decoy_score = max(SVM score over DECOY_-prefixed peptides passing gate)`.
   - Gate: peptide's `run_peptide_qvalue <= config.run_fdr` (Savitski's 1% convention).
4. **Pairwise picking**: each group produces exactly one winner. Target wins if `target_score >= decoy_score`, else decoy wins.
5. Sort winners by score descending, compute cumulative FDR = `cum_decoys / max(1, cum_targets)`, backward sweep for monotonicity.
6. Propagate per-group q-values back to peptides (best/lowest across groups a peptide belongs to) and write into `run_protein_qvalue`.

First-pass protein q-values are used by:

- **Protein-aware compaction**: rescues borderline peptides whose protein passes first-pass FDR, even if their own peptide q-value is above the compaction gate.
- **Reconciliation consensus selection**: peptides from strong proteins can anchor consensus RT computation alongside peptides passing peptide-level FDR directly.

### Second Pass (Post-Reconciliation, Authoritative)

Runs after second-pass peptide FDR. Produces `experiment_protein_qvalue` on each FdrEntry.

1. Build parsimony from peptides passing second-pass peptide FDR.
2. Call `collect_best_peptide_scores()` on the compacted + reconciled + second-pass-scored `per_file_entries`. Scores now reflect reconciliation corrections.
3. Same picked-protein algorithm as first pass (single best peptide per protein, pairwise picking, cumulative FDR on winners).
4. Propagate into `experiment_protein_qvalue` (which feeds `--fdr-level protein` filtering and the protein CSV report).

### Shared Peptide Handling in Scoring

In `All` mode, shared peptides contribute to every protein group they belong to as a candidate for "best peptide". A protein's target-side score is `max(target SVM score)` over its unique peptides plus all shared peptides that touch it. A shared peptide with a high score can raise the best-peptide score of multiple groups simultaneously.

In `Razor` mode, shared peptides are first assigned to exactly one group via iterative greedy set cover (see the algorithm above), so they only contribute to their razor-assigned group's score.

In `Unique` mode, shared peptides are dropped from the parsimony graph entirely, so they never contribute to any protein's score.

### When `--protein-fdr` Is Not Set

Parsimony still runs. The peptide-to-group mapping is available for the blib writer and downstream features. No protein q-values are computed, so `run_protein_qvalue` and `experiment_protein_qvalue` remain at their default (1.0) on all FdrEntry stubs. Compaction falls back to using only the peptide-level q-value gate.

## Edge Cases

- **Peptides with empty `protein_ids`**: excluded from the parsimony graph. These peptides will not be in `peptide_to_groups` and cannot pass protein-level filtering.
- **Decoy entries**: excluded from the target parsimony graph. Decoys are handled separately by the picked-protein FDR via `DECOY_` prefix pairing.
- **Library entries not in `detected_peptides`**: when a detected set is provided, these entries are skipped. The resulting graph only contains peptides with at least one detection.

## Implementation

| File | Function | Purpose |
|------|----------|---------|
| `crates/osprey-fdr/src/protein.rs` | `build_protein_parsimony()` | Main entry point: builds bipartite graph, merges identical sets, eliminates subsets, applies shared peptide mode |
| `crates/osprey-fdr/src/protein.rs` | `compute_protein_fdr()` | Picked-protein q-value computation (only when `--protein-fdr` is set) |
| `crates/osprey-fdr/src/protein.rs` | `propagate_protein_qvalues()` | Writes protein q-values back into FdrEntry stubs |
| `crates/osprey-fdr/src/protein.rs` | `write_protein_report()` | CSV report writer (`*.proteins.csv`) |
| `crates/osprey/src/pipeline.rs` | protein FDR block | Orchestrates parsimony + optional FDR in the pipeline |

## Tests

| Test | Description |
|------|-------------|
| `test_shared_peptides_all_mode` | `All` mode preserves the full peptide-to-groups mapping |
| `test_shared_peptides_razor_mode` | Basic razor: shared peptide goes to the group with more unique peptides |
| `test_shared_peptides_razor_iterative_greedy` | Razor picks global best across a branching case |
| `test_shared_peptides_razor_cascading_assignment` | Two-round razor: cascading assignment after first-round claims |
| `test_shared_peptides_razor_deterministic` | 10 consecutive runs produce identical razor assignments |
| `test_shared_peptides_unique_mode` | `Unique` mode drops shared peptides from the mapping |
| `test_decoy_entries_excluded_from_parsimony` | Decoys are not in the target parsimony graph |

## References

- Protein parsimony general approach: Nesvizhskii AI, Aebersold R. "Interpretation of shotgun proteomic data: the protein inference problem" Mol Cell Proteomics. 2005;4(10):1419-1440.
- Picked-protein FDR: Savitski MM, et al. Mol Cell Proteomics. 2015;14(9):2394-2404.
- Picked-protein group FDR: The M, Käll L. J Proteome Res. 2016;15(4):1456-1461.
- Razor peptides as used in MaxQuant: Cox J, Mann M. Nat Biotechnol. 2008;26(12):1367-1372.
