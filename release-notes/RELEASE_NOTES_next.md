# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

<!-- none yet -->

## Bug Fixes

- **Carafe's per-peptide `_pepNNNNN` pseudo-protein accessions are stripped at library load**: Carafe emits one synthetic accession PER PEPTIDE (`sp|O95139_pep00019|NDUB6_HUMAN`), so every peptide became its own protein — 759,805 distinct accessions collapsing to 42,592 real ones on the 3-file Stellar entrapment library. The damage was not confined to the protein report: second-pass `protein-compact` q-value estimation builds its stratum from proteins with >= 2 *detected* peptides, and a protein owning exactly one peptide can never qualify, so the stratum collapsed, nothing recompeted, and every reported q-value was the carried-over first-pass value with nothing in the output saying so. The manifest's `proteins` column (see v26.6.0) already carried clean accessions, but the only path that applies it is the library-supplied-decoy path, which also hard-fails when the library holds no decoys — so a Carafe library searched with *generated* decoys had no route back to real accessions at all. Stripping at load closes that for every mode; on a manifest run the manifest still wins where it genuinely disagrees. Applied on both the fresh-parse path (before deduplication) and the `.libcache` path. Entrapment and decoy markers are preserved — only the `_pep` + digits token is removed, so `sp|Q9ULK4_p_target_pep00052|MED23_HUMAN_p_target` becomes `sp|Q9ULK4_p_target|MED23_HUMAN_p_target`. Mirrors the C# Osprey `CarafeProteinIdNormalizer` (ProteoWizard/pwiz#4573); the two implementations are validated in cross-impl parity at 1e-9 on `StellarLibraryDecoy`. **Note**: the strip assumes the default `--pep-suffix-format` (`_pep{:05d}`) of `scripts/build_entrapment_peptide_fasta.py`; a library built with a custom format is not recognized.

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
