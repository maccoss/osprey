# 18 - HPC Scoring Split (`--no-join` / `--join-at-pass`)

For large experiments (hundreds to thousands of mzML files), Osprey can split its pipeline into phases that scale very differently:

- **Per-file scoring** (Stages 1-4): per-file, embarrassingly parallel, CPU-bound. Each input mzML is scored independently and writes a `<stem>.scores.parquet` cache next to the input. No cross-file communication.
- **First join** (Stage 5): single-process, requires all per-file caches. Runs first-pass Percolator FDR, multi-charge consensus, cross-run consensus RT computation, per-file LOESS calibration refit, and reconciliation planning. Writes the **boundary file pair** per file: `<stem>.1st-pass.fdr_scores.bin` + `<stem>.reconciliation.json`.
- **Per-file rescore** (Stage 6): per-file again, embarrassingly parallel. Each file consumes its boundary file pair plus its own parquet, runs `run_search` with reconciled boundary overrides, executes the gap-fill two-pass, and rewrites a reconciled `<stem>.scores.parquet`.
- **Second join** (Stages 7-8): single-process. Runs second-pass FDR, optional protein FDR, and writes the final blib.

The CLI surface for orchestrating this is the **`--join-at-pass=<N>` family of flags**, with `--no-join` and `--join-only` as modifiers that pick which side of pass `N` to run.

## When to Use

- **Multi-node HPC clusters**: per-file scoring and rescore are embarrassingly parallel; the two join steps are not. Fanning each fan-out phase across N nodes gives near-linear speedup until the joins become the bottleneck.
- **Cross-implementation interop with OspreySharp**: the parquet caches and JSON boundary files are byte-identical between the Rust and C# implementations. Either binary can pick up the other's fan-out output.
- **Iterative reanalysis**: re-running FDR, reconciliation, or protein parameters without re-scoring is supported via parameter-hash-validated cache reuse (see "Parquet Footer Hash Validation" below).

For small experiments (< ~50 files on a single workstation) the default single-process pipeline is simpler and incurs no orchestration overhead.

## Pipeline Stage Mapping

| Stage | Phase | `--no-join` | Default | `--join-at-pass=1 --join-only` | `--join-at-pass=1 --no-join` | `--join-at-pass=1` |
|---|---|---|---|---|---|---|
| 1 | Library load + decoy generation | yes | yes | from parquet metadata | from parquet metadata | from parquet metadata |
| 2 | Per-file calibration discovery | yes | yes | from `.calibration.json` | from `.calibration.json` | from `.calibration.json` |
| 3 | Coelution search | yes | yes | from parquet | from parquet | from parquet |
| 4 | Feature extraction + per-file Parquet write | yes | yes | from parquet | from parquet | from parquet |
| 5 | First-pass Percolator FDR | no | yes | yes (writes boundary files) | from boundary files | yes |
| 5b | Multi-charge consensus | no | yes | yes (writes boundary files) | from boundary files | yes |
| 6 | Cross-run reconciliation per-file rescore | no | yes | no (exits after 5) | yes | yes |
| 6b | Second-pass Percolator FDR | no | yes | no | no (exits after 6) | yes |
| 6c | First-pass + second-pass protein FDR | no | yes | no | no | yes |
| 7 | blib output | no | yes | no | no | yes |

Reading the table: each column shows what one HPC mode does. The "yes / from X" cells show whether that stage actually runs or whether it's loaded from disk.

## Three HPC Modes

### 1. Per-file scoring — `--no-join`

Runs Stages 1-4 only and exits. Writes `<stem>.scores.parquet`, `<stem>.calibration.json`, and `<stem>.spectra.bin` next to each input mzML. No cross-file communication.

```bash
# One Slurm/PBS array job per file.
osprey --no-join -i file_001.mzML -l library.tsv --resolution hram
osprey --no-join -i file_002.mzML -l library.tsv --resolution hram
# ... fan out across N nodes ...
```

`--no-join` requires `--input <mzML...>` and rejects `--input-scores`. The `--output` flag is accepted for backward compatibility but logs a warning that no blib is written.

### 2. First join — `--join-at-pass=1 --join-only`

Runs Stage 5 only (first-pass FDR + multi-charge consensus + cross-run consensus + reconciliation planning) and exits. Writes the per-file boundary file pair (`<stem>.1st-pass.fdr_scores.bin` + `<stem>.reconciliation.json`).

```bash
osprey --join-at-pass=1 --join-only --input-scores /path/to/parquet/dir \
       -l library.tsv -o results.blib --resolution hram
```

The `--output` flag is required for argument validation but no blib is actually written (the run exits before Stage 7). Requires 2+ input parquets and `reconciliation.enabled=true`; the boundary file pair is only meaningful for multi-file fan-back-in with reconciliation.

### 3. Per-file rescore — `--join-at-pass=1 --no-join`

Runs Stage 6 only (per-file rescore + gap-fill + reconciled parquet write-back) and exits. Consumes the boundary file pair persisted by mode 2 plus the per-file parquet from mode 1, and rewrites each `<stem>.scores.parquet` with reconciled scores. Embarrassingly parallel: one Slurm/PBS array job per parquet.

```bash
osprey --join-at-pass=1 --no-join --input-scores file_001.scores.parquet \
       -l library.tsv -o results.blib --resolution hram
```

`--input` is rejected (mzML paths are derived from the parquet stems via the matching `<stem>.mzML` sibling). `--library` and `--output` are required. The boundary file siblings (`<stem>.1st-pass.fdr_scores.bin` + `<stem>.reconciliation.json`) must be present alongside each parquet.

The worker hydrates the in-memory state at the Stage 5 / Stage 6 seam, applies the same compaction predicate the in-process pipeline uses (`run_peptide_qvalue ≤ peptide_gate OR run_protein_qvalue ≤ protein_gate`), and drives the same `rescore_per_file_loop` that the in-process pipeline calls. **Bit-parity contract**: `Compare-Stage6-Worker.ps1` validates that Phase A in-process and Phase C worker-on-renamed-folder produce byte-identical reconciled parquets across all 40 columns.

### 4. Full join — `--join-at-pass=1`

Runs Stages 5 through 8 from the per-file `<stem>.scores.parquet` outputs (mode 1's output, or mode 3's reconciled outputs). Writes the final blib. This is the legacy single-host join mode, equivalent to the previous `--join-only` flag (which is now a modifier of `--join-at-pass`).

```bash
osprey --join-at-pass=1 --input-scores /path/to/parquet/dir \
       -l library.tsv -o results.blib --resolution hram
```

`--input-scores` accepts either a single directory (scanned **non-recursively** for `*.scores.parquet`) or an explicit list of paths.

### 5. Combined (single-machine, default)

Default invocation runs all stages in one process:

```bash
osprey -i *.mzML -l library.tsv -o results.blib --resolution hram
```

Output is identical to running mode 1 → mode 2 → mode 3 → final-join sequentially over the same files.

## End-to-End HPC Orchestration

A typical 240-file experiment fans out as:

```text
1. Per-file scoring (240 nodes, ~minutes each)
   for each mzML:
     osprey --no-join -i file.mzML -l library.tsv --parquet-compression snappy

2. First join (1 node, ~minutes)
   osprey --join-at-pass=1 --join-only --input-scores /parquet_dir \
          -l library.tsv -o results.blib

3. Per-file rescore (240 nodes, ~minutes each)
   for each parquet:
     osprey --join-at-pass=1 --no-join --input-scores file.scores.parquet \
            -l library.tsv -o results.blib

4. Final join (1 node, ~minutes)
   osprey --join-at-pass=1 --input-scores /parquet_dir \
          -l library.tsv -o results.blib
```

Steps 1 and 3 are fan-out; steps 2 and 4 are fan-in. Each fan-out is a Slurm/PBS array; each fan-in is a single head-node invocation. The `--parquet-compression snappy` flag in step 1 is needed only if the joins will run under OspreySharp (Parquet.Net 3.x supports Snappy only).

## Parquet Footer Hash Validation

The join steps are only safe if every input parquet was scored with **the same library, the same search parameters, and a compatible Osprey binary version**. Osprey enforces this via SHA-256 hashes stored in each parquet's footer metadata, validated **before any work is done** by `--join-at-pass=<N>`.

| Hash | Covers | Failure mode |
|---|---|---|
| `osprey.version` | Osprey binary major/minor version | Different binary versions may extract different feature sets — refuse to join |
| `osprey.search_hash` | Resolution mode, fragment/precursor tolerance, prefilter setting, decoy method, RT calibration parameters, `reconciliation.top_n_peaks` | Search parameter drift would mix incomparable scores in FDR — refuse to join |
| `osprey.library_hash` | Library file name + size + filesystem mtime | Library file changed between scoring and joining — refuse to join (the entry IDs in the parquet may not resolve correctly) |
| `osprey.reconciled` | Whether this parquet has already been through reconciliation | Used to short-circuit reconciliation on full-cache reruns |
| `osprey.reconciliation_hash` | Reconciliation parameters + sorted file set | Used together with `osprey.reconciled` for the skip-reconciliation fast path |

On hash mismatch, `--join-at-pass=<N>` **aborts with a clear, file-named error** before any FDR or reconciliation work begins. The error message names the offending parquet and the mismatched hash so the operator can re-score the affected file.

The same hashes are also written into the `<stem>.reconciliation.json` boundary file (`search_hash` and `library_hash` fields). The Stage 6 worker validates these match the per-file parquet before proceeding, so a boundary file accidentally paired with the wrong parquet is caught early.

For the full hash composition and the `CacheValidity` enum (`ValidReconciled` / `ValidFirstPass` / `Stale`), see [12 - Intermediate File Formats § Cache Invalidation via Parquet Metadata](12-intermediate-files.md#cache-invalidation-via-parquet-metadata).

## Cross-Implementation Compatibility

The parquet cache format and the `<stem>.reconciliation.json` boundary file are **shared between Rust Osprey and OspreySharp** (the C# port in ProteoWizard/pwiz). Either binary can pick up the other's output at any seam in the pipeline.

```bash
# Score with Snappy compression so OspreySharp can read the result.
osprey --no-join -i file.mzML -l library.tsv --parquet-compression snappy
```

Default compression is `zstd` (smaller files). OspreySharp uses Parquet.Net 3.x which supports Snappy only, so cross-impl interop requires `--parquet-compression snappy` on the **scoring** step. Reading auto-dispatches on per-column-chunk metadata, so this only governs writes. When Snappy is selected, the writer also disables `RLE_DICTIONARY` encoding (unsupported by Parquet.Net 3.x).

The reconciliation JSON uses **alphabetical key ordering at every nesting level** and routes every f64 through a shared `format_f64_roundtrip` helper so the file is byte-identical when written by either implementation. Numeric values use a single canonical fixed-point decimal form (no scientific notation) on both sides; the `serde_json` `float_roundtrip` feature flag ensures exact f64 round-trip on the Rust read path. See [docs/12 § Reconciliation Boundary File](12-intermediate-files.md#5-reconciliation-boundary-file-stemreconciliationjson) for the full schema.

Cross-impl byte-parity at each pipeline seam is validated via env-var-gated diagnostic dumps in [crates/osprey/src/diagnostics.rs](../crates/osprey/src/diagnostics.rs):

- **Stage 5 dumps**: `OSPREY_DUMP_STANDARDIZER`, `OSPREY_DUMP_SUBSAMPLE`, `OSPREY_DUMP_SVM_WEIGHTS`, `OSPREY_DUMP_PERCOLATOR`
- **Stage 6 planning dumps**: `OSPREY_DUMP_CALIBRATION`, `OSPREY_DUMP_PROTEIN_FDR`, `OSPREY_DUMP_CONSENSUS`, `OSPREY_DUMP_INV_PREDICT`, `OSPREY_DUMP_MULTICHARGE`, `OSPREY_DUMP_REFIT`, `OSPREY_DUMP_LOESS_FIT`, `OSPREY_DUMP_RECONCILIATION`
- **Stage 6 worker bisection dumps**: `OSPREY_DUMP_MP_INPUTS` (per-call inputs to `tukey_median_polish`), `OSPREY_DUMP_PREDICT_RT` (per-call inputs/outputs of `RTCalibration::predict` plus per-file calibration arrays)

All dumps use the `OSPREY_DUMP_<NAME>` / `OSPREY_<NAME>_ONLY` gate convention; the `_ONLY` companion exits after writing for fast bisection. Production runs see zero overhead with these unset.

## Persistent Files Per Input

Each input mzML produces persistent caches alongside it:

| File | Written by | Reused by |
|---|---|---|
| `<stem>.calibration.json` | Stage 2 (auto-calibration) | All later stages — loaded for inverse RT prediction in Stage 6 |
| `<stem>.scores.parquet` | Stage 4 (feature extraction); rewritten by Stage 6 with reconciled scores | Loaded selectively per phase: stubs, PIN features, CWT candidates, blib plan |
| `<stem>.spectra.bin` | Stage 1-2 (mzML decode) | Local-only fast-reload cache; not needed by joins |
| `<stem>.1st-pass.fdr_scores.bin` | Stage 5 (first-pass FDR + protein FDR) | Stage 6 worker (mandatory); skip-Percolator optimization on rerun |
| `<stem>.2nd-pass.fdr_scores.bin` | Stage 6b (second-pass FDR) | Skip-Percolator optimization on rerun |
| `<stem>.reconciliation.json` | Stage 5 (reconciliation planner) | Stage 6 worker (mandatory) |

`<stem>.spectra.bin` is a per-node fast-reload cache for mzML decoding. It does not need to travel between scoring and joining nodes; if you're using shared filesystem mounts it will be reused on subsequent scoring runs of the same file, but no join step depends on it.

The boundary file pair (`<stem>.1st-pass.fdr_scores.bin` + `<stem>.reconciliation.json`) is required for `--join-at-pass=1 --no-join` (Stage 6 worker mode). The persisted format includes `search_hash` and `library_hash` fields that the worker validates against the per-file parquet before proceeding.

For full cache file format documentation, see [docs/12 — Intermediate File Formats](12-intermediate-files.md).

## Failure Modes and Mutual Exclusion

The CLI rejects illegal flag combinations early, with specific error messages:

| Combination | Error |
|---|---|
| `--no-join --join-only` | "`--no-join` and `--join-only` are mutually exclusive." |
| `--input-scores` without `--join-at-pass` | "`--input-scores` requires `--join-at-pass=1` (with or without `--join-only` / `--no-join`)." |
| `--join-at-pass=1` without `--input-scores` | "`--join-at-pass=1` requires `--input-scores <path...>`." |
| `--join-at-pass=1 --no-join` without `--input-scores` | "`--join-at-pass=1 --no-join` (per-file rescore worker) requires `--input-scores <path...>`." |
| `--join-at-pass=1 --no-join` with `--input` | "`--join-at-pass=1 --no-join` cannot be combined with `--input`. Use `--input-scores`; mzML paths are derived from the parquet stems." |
| `--join-at-pass=1 --no-join` without `--library` or `--output` | "`--join-at-pass=1 --no-join` requires `--library` and `--output`." |
| `--no-join --input-scores` without `--join-at-pass=1` | "`--no-join` cannot be combined with `--input-scores`." |
| `--no-join` without `--input` | "`--no-join` requires `--input <mzML...>`." |
| `--input-scores <dir>` with no `*.scores.parquet` files | "No `*.scores.parquet` files found in `--input-scores` directory `<dir>`." |
| `--input-scores <path>` not found | "`--input-scores` path not found: `<path>`." |
| `--join-at-pass=1 --join-only` with single input parquet | "`--join-at-pass=1 --join-only` requires 2+ input parquets. The Stage 5 → Stage 6 boundary file pair is only meaningful for multi-file fan-back-in." |
| `--join-at-pass=1 --join-only` with `reconciliation.enabled=false` | "`--join-at-pass=1 --join-only` requires `reconciliation.enabled = true`." |
| Parquet hash mismatch under any join mode | Aborts before any work, naming the offending file and which hash mismatched. |

## Safe Concurrent Writes

All cache writers (parquet, calibration JSON, spectra binary, FDR score sidecars, reconciliation JSON) use the `copy_and_verify` pattern: the writer first writes to a local temp file, verifies the write completed, then moves to the final destination. This makes both fan-out steps safe on NAS / CIFS mounts where a process kill mid-write could otherwise leave a corrupted cache that downstream stages would have to reject and recompute. See [09 - Determinism § Safe NAS File Writes](09-determinism.md) for details.

## See Also

- [12 - Intermediate File Formats](12-intermediate-files.md) — full cache file documentation, including the v3 sidecar format and reconciliation JSON schema
- [07 - FDR Control](07-fdr-control.md) — what each join step does (first-pass, second-pass, protein FDR)
- [10 - Cross-Run Reconciliation](10-cross-run-reconciliation.md) — Stage 5 planner + Stage 6 rescore algorithm
- [crates/osprey/src/main.rs](../crates/osprey/src/main.rs) — `validate_hpc_args`, `normalize_hpc_args`, and `resolve_input_scores` logic
- [crates/osprey/src/rescore.rs](../crates/osprey/src/rescore.rs) — Stage 6 worker hydration and dispatch (`hydrate_for_rescore`, `run_rescore`)
- [crates/osprey/src/reconciliation_io.rs](../crates/osprey/src/reconciliation_io.rs) — boundary file reader/writer
