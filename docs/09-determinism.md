# Determinism

Osprey is designed to produce **bit-identical results** across runs on the same input data, regardless of thread scheduling or platform. This is critical for reproducibility and debugging. This document describes the patterns used to maintain determinism and what to watch for when modifying the code.

## Sources of Non-Determinism

In a parallel Rust pipeline, four sources can break determinism:

1. **HashMap iteration order** — Rust's `HashMap` uses randomized hashing by default, so iteration order is non-deterministic across runs
2. **Rayon parallel iteration order** — `par_iter().flat_map()` and `par_iter().map().collect()` do not guarantee output order
3. **Floating-point accumulation order** — `a + (b + c)` may differ from `(a + b) + c` due to rounding, so parallel reductions with different thread scheduling can produce different sums
4. **Random number generation** — Any use of random sampling or shuffling must use a fixed seed

## Patterns Used in Osprey

### 1. Sort After HashMap Collection

Every time results are collected from a `HashMap`, they are immediately sorted into a deterministic order.

**Where this applies:**
- `deduplicate_pairs()`: Groups entries by `base_id` using HashMap, then sorts results by `entry_id` (`pipeline.rs`)
- `select_consensus_peaks()`: Groups entries by `modified_sequence` using HashMap, result sorted by `(entry_id, scan_number)` (`pipeline.rs`)
- `compete_calibration_pairs()`: Collects winners from HashMap, sorts by `(score, base_id)` (`calibration_ml.rs`)
- Best-per-precursor subsampling: Collects best target/decoy per `base_id` from HashMap, sorts `best_observations` by `(file_idx, local_idx)`, sorts peptide groups alphabetically before subsampling (`pipeline.rs`)

**Pattern:**
```rust
let mut groups: HashMap<K, Vec<V>> = HashMap::new();
// ... populate ...
let mut results: Vec<V> = groups.into_values().flatten().collect();
results.sort_by_key(|e| e.some_deterministic_key);
```

### 2. Sort After Parallel Collection

After any `par_iter()` that produces results, sort the output before further processing.

**Where this applies:**
- Main search: `all_entries` collected from parallel per-window processing, sorted by `entry_id` before FDR (`pipeline.rs`)
- Overlapping-window dedup: `deduped` sorted by `(entry_id, scan_number)` after HashMap-based dedup (`pipeline.rs`)
- Multi-charge consensus: Results re-sorted by `(entry_id, scan_number)` after merging kept and re-scored entries (`pipeline.rs`)

**Pattern:**
```rust
let results: Vec<T> = entries.par_iter().flat_map(|e| process(e)).collect();
let mut sorted = results;
sorted.sort_by(|a, b| a.key.cmp(&b.key));
```

### 3. Seeded Random Number Generation

All random operations use a fixed seed for reproducibility.

**Where this applies:**
- Percolator SVM cross-validation: seed=42 for fold assignment and weight initialization (`percolator.rs`)
- Library subsampling: Deterministic selection based on entry_id hash, not random sampling

### 4. Three-Level Tiebreaking with `total_cmp`

When comparing floating-point values, ties must be broken deterministically. Osprey uses multi-level sort keys:

**Pattern:**
```rust
// Sort by score descending, then entry_id ascending for deterministic tiebreaking
entries.sort_by(|a, b| {
    b.score.total_cmp(&a.score)       // Primary: score descending
        .then(a.entry_id.cmp(&b.entry_id))  // Tiebreak: entry_id ascending
});
```

**Why `total_cmp`?** Standard `partial_cmp` returns `None` for NaN comparisons and doesn't define a total order. `total_cmp` (stable since Rust 1.62) treats NaN as greater than all other values and provides a total order, making sorts deterministic even with NaN values.

### 5. Sort Window Groups

DIA isolation windows extracted from mzML are sorted by lower m/z bound before processing:

```rust
first_cycle_windows.sort_by(|a, b| a.0.total_cmp(&b.0));
```

This ensures windows are processed in the same order regardless of how they appear in the mzML file.

### 6. Deterministic Fold Assignment

Cross-validation fold assignment uses a deterministic hash of the peptide sequence:

```rust
// Fold assignment based on sequence hash, not random
let fold = hash(modified_sequence) % n_folds;
```

This ensures the same peptide always goes to the same fold across runs, making cross-validation reproducible.

### 7. Overlapping-Window Deduplication Tiebreaking

When the same precursor is detected in two overlapping DIA windows, the entry with the higher `coelution_sum` wins. Ties (same coelution_sum) are broken by `entry_id`:

```rust
// Iteration order is deterministic (entries sorted by apex_rt + entry_id),
// so the >= tiebreaker (outer loop wins ties) is also deterministic.
```

### 8. Parquet Round-Trip Determinism

Per-file score caches (`.scores.parquet`) are written once after scoring and reloaded on-demand for FDR, reconciliation, and blib output. The Parquet format preserves exact floating-point values (no rounding), so entries rehydrated from Parquet are bit-identical to the original in-memory entries. Column ordering and row ordering in the Parquet file match the deterministic entry order established by the sort-after-collect patterns above.

### 9. Iterative Greedy Set Cover for Razor Peptide Assignment

The `SharedPeptideMode::Razor` shared-peptide assignment (`crates/osprey-fdr/src/protein.rs`) resolves shared peptides to a single protein group each. A naive single-pass implementation over a `HashMap<peptide, groups>` would be non-deterministic because the iteration order leaks into the assignment order, and each assignment mutates the group's unique count that later iterations depend on.

Osprey uses an **iterative greedy set cover** that makes the decision globally at each step:

1. Collect shared peptides up-front and **sort alphabetically**
2. At each round, find the protein group with the most unique peptides that still has at least one unassigned shared peptide (tiebreak: lowest group ID)
3. That group claims **all** its unassigned shared peptides in one batch, processed in sorted order
4. Repeat until no shared peptides remain

This is path-independent: the same input always produces the same output, regardless of HashMap hashing or iteration order. Verified by `test_shared_peptides_razor_deterministic` which runs parsimony 10 times on the same input and asserts byte-identical assignments. See [16-protein-parsimony.md](16-protein-parsimony.md) for the full algorithm.

---

## Checklist for New Code

When adding new code to Osprey, verify these invariants:

### If you use a HashMap or HashSet:
- [ ] Results collected from the map are sorted before use
- [ ] The sort key is based on stable identifiers (`entry_id`, `modified_sequence`), not memory addresses or insertion order

### If you use `par_iter()` or `par_bridge()`:
- [ ] Collected results are sorted into deterministic order before downstream use
- [ ] No parallel floating-point reduction (use sequential sum, or sort-then-sum)

### If you compare floating-point values:
- [ ] Use `total_cmp` instead of `partial_cmp` for sort comparisons
- [ ] Add tiebreaker keys when scores can be equal (e.g., `.then(a.entry_id.cmp(&b.entry_id))`)
- [ ] Check for NaN: any division should guard against zero denominators

### If you sample or shuffle:
- [ ] Use a seeded RNG (seed=42 by convention)
- [ ] Document the seed location so it can be overridden if needed

### If you add cross-validation or fold splitting:
- [ ] Target-decoy pairs must stay together (same `base_id`)
- [ ] Same peptide, different charge states must stay together
- [ ] Subsampling operates on paired groups, not individual entries
- [ ] See [FDR Control](07-fdr-control.md) for the full grouping invariant

### If you write/read Parquet caches:
- [ ] Entries are written in a deterministic order (sorted before writing)
- [ ] Reloaded entries are used in the same order as they were written
- [ ] No floating-point rounding is introduced by the serialize/deserialize round-trip

---

## Known Determinism-Critical Paths

| Location | What | Why |
|----------|------|-----|
| `pipeline.rs:run_search()` | Sort `all_entries` by `entry_id` | Rayon parallel window processing |
| `pipeline.rs:run_search()` | Sort `deduped` by `(entry_id, scan_number)` | HashMap-based overlapping-window dedup |
| `pipeline.rs:run_search()` | Re-sort after multi-charge consensus | HashMap-based grouping + re-scoring merge |
| `pipeline.rs` | Best-per-precursor: sort `best_observations`, sort peptide groups | HashMap-based best target/decoy per base_id |
| `pipeline.rs:deduplicate_pairs()` | Sort by `entry_id` | HashMap-based target-decoy pairing |
| `calibration_ml.rs` | Sort winners by `(score, base_id)` | HashMap-based target-decoy competition |
| `calibration_ml.rs` | Fold assignment by sequence hash | Cross-validation reproducibility |
| `percolator.rs` | Seeded Xorshift64 (seed=42) | SVM weight initialization |
| `pipeline.rs` | Window groups sorted by lower m/z | Window processing order |
| `pipeline.rs` | FdrEntry stubs maintain entry order | Parquet write order = sorted entry order |
| `pipeline.rs` | Parquet round-trip preserves f64 values | No floating-point drift across reload |
| `pipeline.rs` | Sequential re-scoring (`iter_mut()`) | Avoids nondeterministic parallel float accumulation + OOM |
| `pipeline.rs:rt_mz_index` | Entries sorted by expected_rt within m/z bins | Binary search reproducibility |
| `protein.rs` | Iterative greedy set cover for Razor mode | Path-independent shared peptide assignment |

---

## Testing

The test `test_deduplicate_pairs_deterministic` in `pipeline.rs` verifies that `deduplicate_pairs()` produces deterministic output despite HashMap internals. This catches the most common source of non-determinism (HashMap iteration order affecting downstream row ordering for SVM training).

To verify full-pipeline determinism, run Osprey twice on the same input and compare outputs:
```bash
diff <(sqlite3 run1.blib "SELECT * FROM RetentionTimes ORDER BY 1,2") \
     <(sqlite3 run2.blib "SELECT * FROM RetentionTimes ORDER BY 1,2")
```

---

## Implementation

| File | Determinism-related code |
|------|--------------------------|
| `crates/osprey/src/pipeline.rs` | All sort-after-collect patterns, window group sorting, dedup tiebreaking |
| `crates/osprey-scoring/src/calibration_ml.rs` | Fold assignment, competition sorting |
| `crates/osprey-fdr/src/lib.rs` | Percolator seed, fold grouping |
| `crates/osprey-fdr/src/protein.rs` | Iterative greedy set cover, sorted peptide processing |
