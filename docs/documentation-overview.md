# Osprey Documentation Overview

## Docs Ordered by Pipeline Phase

The numbered docs roughly follow the Osprey pipeline (5 phases), with some cross-cutting concerns at the end:

### Phase 1: Initialization

| Doc | Summary |
|-----|---------|
| [01-decoy-generation.md](01-decoy-generation.md) | Enzyme-aware sequence reversal (preserves cleavage-site residues), fragment m/z recalculation for b/y ion swaps |

### Phase 2: Calibration Discovery

| Doc | Summary |
|-----|---------|
| [02-calibration.md](02-calibration.md) | Joint RT + MS1 + MS2 calibration from high-confidence matches; LOESS fitting, mass error stats, per-file caching |
| [03-spectral-scoring.md](03-spectral-scoring.md) | Calibration-phase scoring: XCorr, E-value, isotope cosine, fragment correlation, libcosine -- combined via LDA |
| [04-xcorr-scoring.md](04-xcorr-scoring.md) | Comet-style XCorr implementation details (preprocessing, bin sizes, flanking subtraction, BLAS sdot) |

### Phase 3: Coelution Search

| Doc | Summary |
|-----|---------|
| [05-peak-detection.md](05-peak-detection.md) | CWT consensus peak detection (Mexican Hat wavelet, median across transitions, zero-crossing boundaries, valley guard) |

### Phase 4: First-Pass FDR + Post-FDR Consensus

| Doc | Summary |
|-----|---------|
| [06-multi-charge-consensus.md](06-multi-charge-consensus.md) | Forces charge states of the same peptide to share peak boundaries; best SVM-scoring charge state leads |
| [07-fdr-control.md](07-fdr-control.md) | Two-level FDR (run + experiment), three methods (Percolator SVM, Mokapot, simple TDC), dual precursor+peptide q-values, `FdrLevel` filtering (precursor/peptide/protein/both) |
| [16-protein-parsimony.md](16-protein-parsimony.md) | Always-on parsimony: bipartite graph, identical-set merging, subset elimination, iterative greedy razor set cover; optional picked-protein FDR |

### Phase 5: Cross-Run Reconciliation (multi-file)

| Doc | Summary |
|-----|---------|
| [10-cross-run-reconciliation.md](10-cross-run-reconciliation.md) | Consensus RT computation, LOESS refit, per-entry reconciliation plan (Keep / UseCwtPeak / ForcedIntegration) |
| [13-boundary-overrides.md](13-boundary-overrides.md) | How `run_search()` accepts override boundaries for consensus and reconciliation re-scoring |
| [14-rt-alignment.md](14-rt-alignment.md) | LOWESS calibration, inverse prediction, weighted median consensus RT, peak imputation |

### Output

| Doc | Summary |
|-----|---------|
| [08-blib-output-schema.md](08-blib-output-schema.md) | BiblioSpec SQLite schema: standard tables + Osprey extensions (peak boundaries, run/experiment scores), nullable retentionTime for Skyline ID lines |

### Cross-Cutting / Development

| Doc | Summary |
|-----|---------|
| [09-determinism.md](09-determinism.md) | Patterns for bit-identical results across runs (thread ordering, float stability, fold assignment) |
| [11-testing.md](11-testing.md) | Testing guide: all 317 unit tests in `#[cfg(test)]` blocks, no separate integration tests |
| [12-intermediate-files.md](12-intermediate-files.md) | Calibration JSON, binary spectra cache, Parquet score caches -- all written alongside input mzML |
| [DIA-PeakDetectionPlan.md](DIA-PeakDetectionPlan.md) | Original technical spec for the CWT peak detection design (predates implementation) |

## How the Three Top-Level Docs Relate

### README.md (user-facing)

- Installation (Linux/Windows/macOS), quick start, CLI options, output format
- High-level algorithm overview (calibration strategy, coelution scoring, reconciliation, memory architecture)
- Links to [08-blib-output-schema.md](08-blib-output-schema.md), [10-cross-run-reconciliation.md](10-cross-run-reconciliation.md), and [docs/README.md](README.md)

### docs/README.md (algorithm reference)

- The most detailed doc: full 5-phase pipeline diagram, all 21 PIN features with descriptions, 26 removed features, Tukey median polish explanation, FDR control logic, memory architecture, crate-level architecture diagram, configuration reference
- Links to every numbered doc (01-14) for deep dives

### CLAUDE.md (developer context for Claude Code)

- Build commands, CI requirements (`cargo fmt`, `clippy -D warnings`, `cargo test`)
- Critical invariants (cross-validation fold grouping, target-decoy pairing)
- Key file paths and function names for navigation
- Points to specific docs:
  - [07-fdr-control.md](07-fdr-control.md) for two-level FDR details
  - [08-blib-output-schema.md](08-blib-output-schema.md) for blib format and Skyline integration
  - [12-intermediate-files.md](12-intermediate-files.md) for intermediate file formats
- Contains the "Recent Changes" changelog and "Current Status" not in the other docs
- Duplicates some content from docs/README.md (feature set, FDR logic, memory architecture) to keep it self-contained for Claude's context window

### Overlap / Redundancy

The FDR control logic, memory architecture, and feature set descriptions appear in all three files. CLAUDE.md is the most up-to-date (it tracks recent changes), docs/README.md is the most detailed, and README.md is the most concise. The numbered docs are the canonical deep-dive references -- CLAUDE.md and docs/README.md primarily link to them rather than duplicating their full content.
