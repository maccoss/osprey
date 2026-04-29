# 18 - HPC Scoring Split (`--no-join` / `--join-only`)

For large experiments (hundreds to thousands of mzML files), Osprey can split its pipeline into two phases that scale very differently:

- **Scoring** (Stages 1-4): per-file, embarrassingly parallel, CPU-bound. Each input mzML is scored independently and writes a `<stem>.scores.parquet` cache next to the input. No cross-file communication.
- **Joining** (Stages 5+): single-process, requires all per-file caches. Runs Percolator FDR, cross-run reconciliation, second-pass FDR, optional protein FDR, and writes the final blib.

The `--no-join` and `--join-only` flags decouple these two phases so scoring can fan out across an HPC cluster while joining runs once on a head node.

## When to Use

- **Multi-node HPC clusters**: scoring time scales linearly with file count and dominates the total runtime on large experiments. Fanning scoring out across N nodes gives near-linear speedup until joining becomes the bottleneck.
- **Cross-implementation interop with OspreySharp**: the parquet caches written by Rust Osprey can be joined by either the Rust `osprey` binary or the C# `OspreySharp` port (see "Cross-Implementation Compatibility" below).
- **Iterative reanalysis**: if you need to re-run FDR or reconciliation with different parameters but keep the expensive scoring step, `--join-only` reuses existing caches without re-scoring (parquet footer hashes guard against parameter mismatch).

For small experiments (< ~50 files on a single workstation) the default single-process pipeline is simpler and incurs no orchestration overhead.

## Pipeline Stage Mapping

| Stage | Phase | `--no-join` | Default | `--join-only` |
|---|---|---|---|---|
| 1 | Library load + decoy generation | yes | yes | from parquet metadata |
| 2 | Per-file calibration discovery | yes | yes | from `.calibration.json` |
| 3 | Coelution search | yes | yes | from parquet |
| 4 | Feature extraction + per-file Parquet write | yes | yes | from parquet |
| 5 | First-pass Percolator FDR | no | yes | yes |
| 5b | Multi-charge consensus | no | yes | yes |
| 6 | Cross-run reconciliation + second-pass FDR | no | yes | yes |
| 6b | First-pass + second-pass protein FDR | no | yes | yes |
| 7 | blib output | no | yes | yes |

`--no-join` exits cleanly after Stage 4; `--join-only` enters at Stage 5.

## Usage

### Scoring step (per file, fan out across nodes)

```bash
# Score each file independently. Writes <stem>.scores.parquet,
# <stem>.calibration.json, and <stem>.spectra.bin next to the input mzML.
osprey --no-join -i file_001.mzML -l library.tsv --resolution hram
osprey --no-join -i file_002.mzML -l library.tsv --resolution hram
# ... fan out across N nodes ...
```

Each scoring invocation is fully independent: no shared state, no cross-file communication. A typical cluster orchestration is one Slurm/PBS array job per file.

`--no-join` requires `--input <mzML...>` and rejects `--input-scores`. The `--output` flag is accepted for backward compatibility but logs a warning that no blib is written.

### Joining step (single-process, on the head node)

```bash
# Once all per-file parquets are written, join once on a single node.
osprey --join-only --input-scores /path/to/parquet/dir \
       -l library.tsv -o results.blib --resolution hram
```

`--input-scores` accepts either:
- A single directory, scanned **non-recursively** for `*.scores.parquet`, or
- An explicit list of `.scores.parquet` paths.

`--join-only` requires `--library` and `--output`, rejects `--input` (the mzML inputs are no longer needed; they're discoverable from the parquet metadata if needed downstream). The parquet count after expansion determines whether reconciliation runs: a single parquet skips reconciliation; two or more enter the multi-file Stage 6 path.

### Combined (single-machine)

The default invocation (no `--no-join` or `--join-only`) runs both phases in one process:

```bash
osprey -i *.mzML -l library.tsv -o results.blib --resolution hram
```

This is identical in output to running `--no-join` for each file followed by `--join-only` over the resulting parquets.

## Parquet Footer Hash Validation

The join step is only safe if every input parquet was scored with **the same library, the same search parameters, and a compatible Osprey binary version**. Osprey enforces this via SHA-256 hashes stored in each parquet's footer metadata, validated **before any work is done** by `--join-only`.

| Hash | Covers | Failure mode |
|---|---|---|
| `osprey.version` | Osprey binary major/minor version | Different binary versions may extract different feature sets — refuse to join |
| `osprey.search_hash` | Resolution mode, fragment/precursor tolerance, prefilter setting, decoy method, RT calibration parameters, `reconciliation.top_n_peaks` | Search parameter drift would mix incomparable scores in FDR — refuse to join |
| `osprey.library_hash` | Library file path + size + filesystem mtime | Library file changed between scoring and joining — refuse to join (the entry IDs in the parquet may not resolve correctly) |
| `osprey.reconciled` | Whether this parquet has already been through reconciliation | Used to short-circuit reconciliation on full-cache reruns |
| `osprey.reconciliation_hash` | Reconciliation parameters + sorted file set | Used together with `osprey.reconciled` for the skip-reconciliation fast path |

On hash mismatch, `--join-only` **aborts with a clear, file-named error** before any scoring, FDR, or reconciliation work begins. The error message names the offending parquet and the mismatched hash so the operator can re-score the affected file.

For the full hash composition and the `CacheValidity` enum (`ValidReconciled` / `ValidFirstPass` / `Stale`), see [12 - Intermediate File Formats § Cache Invalidation via Parquet Metadata](12-intermediate-files.md#cache-invalidation-via-parquet-metadata).

## Cross-Implementation Compatibility

The parquet cache format is **shared between Rust Osprey and OspreySharp** (the C# port in ProteoWizard/pwiz). Either binary can join parquets written by the other, provided the compression codec is set to a format both implementations support.

```bash
# Score with Snappy compression so OspreySharp can read the result.
osprey --no-join -i file.mzML -l library.tsv --parquet-compression snappy
```

Default compression is `zstd` (smaller files). OspreySharp uses Parquet.Net 3.x which supports Snappy only, so cross-impl interop requires `--parquet-compression snappy` on the **scoring** step. Reading auto-dispatches on per-column-chunk metadata, so this only governs writes. When Snappy is selected, the writer also disables `RLE_DICTIONARY` encoding (unsupported by Parquet.Net 3.x).

Cross-impl byte-parity at the join step is validated via the env-var-gated diagnostic dumps in [crates/osprey/src/diagnostics.rs](../crates/osprey/src/diagnostics.rs). See `OSPREY_DUMP_PERCOLATOR`, `OSPREY_DUMP_CONSENSUS`, `OSPREY_DUMP_REFIT`, and the other Stage 5 / Stage 6 dumps. These are no-op when the gates are unset.

## Persistent Files Per Input

Each input mzML produces three persistent caches alongside it:

| File | Written by | Reused by |
|---|---|---|
| `<stem>.calibration.json` | `--no-join` (Stage 2) | `--join-only` (loaded for inverse RT prediction in Stage 6 reconciliation) |
| `<stem>.scores.parquet` | `--no-join` (Stage 4) | `--join-only` (loaded selectively per phase: stubs, PIN features, CWT candidates, blib plan) |
| `<stem>.spectra.bin` | `--no-join` (Stages 1-2) | Local-only — not needed by `--join-only` (reconciliation re-scoring loads spectra from mzML if needed) |

`<stem>.spectra.bin` is a per-node fast-reload cache for mzML decoding. It does not need to travel between scoring and joining nodes; if you're using shared filesystem mounts it will be reused on subsequent scoring runs of the same file, but the joining step does not depend on it.

The two FDR score sidecars (`<stem>.1st-pass.fdr_scores.bin`, `<stem>.2nd-pass.fdr_scores.bin`) are written by the joining step, not by `--no-join`. They enable skipping Percolator training entirely on subsequent `--join-only` reruns with matching parameters; see [12 - Intermediate File Formats § Skip-Percolator Optimization](12-intermediate-files.md#skip-percolator-optimization).

## Failure Modes and Mutual Exclusion

The CLI rejects illegal flag combinations early, with specific error messages:

| Combination | Error |
|---|---|
| `--no-join --join-only` | "`--no-join` and `--join-only` are mutually exclusive." |
| `--join-only --input` | "`--join-only` cannot be combined with `--input`. Use `--input-scores`." |
| `--join-only` without `--input-scores` | "`--join-only` requires `--input-scores <path...>`." |
| `--join-only` without `--library` or `--output` | "`--join-only` requires `--library` and `--output`." |
| `--no-join --input-scores` | "`--no-join` cannot be combined with `--input-scores`." |
| `--no-join` without `--input` | "`--no-join` requires `--input <mzML...>`." |
| `--input-scores <dir>` with no `*.scores.parquet` files | "No `*.scores.parquet` files found in `--input-scores` directory `<dir>`." |
| `--input-scores <path>` not found | "`--input-scores` path not found: `<path>`." |
| Parquet hash mismatch under `--join-only` | Aborts before any work, naming the offending file and which hash mismatched. |

## Safe Concurrent Writes

All cache writers (parquet, calibration JSON, spectra binary, FDR score sidecars) use the `copy_and_verify` pattern: the writer first writes to a local temp file, verifies the write completed, then moves to the final destination. This makes the scoring step safe on NAS / CIFS mounts where a process kill mid-write could otherwise leave a corrupted cache that `--join-only` would then have to reject and re-score. See [09 - Determinism § Safe NAS File Writes](09-determinism.md) for details.

## See Also

- [12 - Intermediate File Formats](12-intermediate-files.md) — full cache file documentation, metadata hashing, and `CacheValidity` enum
- [07 - FDR Control](07-fdr-control.md) — what the join step does in Stages 5 and 6
- [10 - Cross-Run Reconciliation](10-cross-run-reconciliation.md) — the multi-file Stage 6 algorithm
- [crates/osprey/src/main.rs](../crates/osprey/src/main.rs) — `validate_hpc_args` and `resolve_input_scores` logic
