# Testing Guide

Osprey uses Rust's built-in `#[test]` framework. All tests are unit tests located in `#[cfg(test)] mod tests` blocks within the source files. There are no integration tests in separate `tests/` directories.

## Running Tests

```bash
# Run all tests across all crates
cargo test --all

# Run tests for a specific crate
cargo test -p osprey-core
cargo test -p osprey-io
cargo test -p osprey-chromatography
cargo test -p osprey-scoring
cargo test -p osprey-fdr
cargo test -p osprey-ml
cargo test -p osprey

# Run a specific test by name
cargo test test_pearson_known_value

# Run tests matching a pattern
cargo test xcorr
cargo test competition

# Show output from passing tests (normally suppressed)
cargo test -- --nocapture

# Run tests in release mode (faster execution, slower compilation)
cargo test --release --all
```

### CI Requirements

Before committing, always run all three checks:

```bash
cargo fmt          # Format code — CI rejects formatting diffs
cargo clippy --all-targets --all-features -- -D warnings  # Lint — CI treats warnings as errors
cargo test --all   # Run all tests
```

## Test Overview

**Total: ~290 tests** across 7 crates.

| Crate | Tests | Focus |
|---|---|---|
| `osprey-core` | 28 | Types, config, isotope distributions |
| `osprey-io` | 26 | Library loading (DIA-NN, elib, blib), blib output |
| `osprey-chromatography` | 66 | Peak detection, CWT, RT/mass calibration |
| `osprey-scoring` | 70 | Spectral scoring, XCorr, coelution, decoy generation |
| `osprey-fdr` | 41 | FDR control, Percolator, Mokapot integration |
| `osprey-ml` | 30 | SVM, LDA, PEP estimation, q-value computation |
| `osprey` | 26 | Pipeline logic, candidate selection, consensus, reconciliation |

---

## Tests by Crate

### osprey-core (28 tests)

Core types, configuration, and isotope calculations.

#### Isotope Distribution — [crates/osprey-core/src/isotope.rs](../crates/osprey-core/src/isotope.rs)

| Test | Description |
|---|---|
| `test_amino_acid_compositions` | Elemental compositions for Ala, Cys, Met |
| `test_peptide_composition` | Peptide composition includes terminal H₂O |
| `test_binomial_coefficient` | Binomial coefficient for known values |
| `test_isotope_distribution_small_peptide` | Dominant M+0 peak for small peptides |
| `test_isotope_distribution_larger_peptide` | Significant M+1 for larger peptides |
| `test_isotope_cosine_perfect_match` | Identical distributions yield cosine = 1.0 |
| `test_isotope_cosine_orthogonal` | Orthogonal distributions yield near-zero |
| `test_peptide_isotope_cosine` | End-to-end peptide isotope scoring |
| `test_sulfur_effect` | Sulfur-containing peptides have higher M+2 ratios |

#### Types — [crates/osprey-core/src/types.rs](../crates/osprey-core/src/types.rs)

| Test | Description |
|---|---|
| `test_isolation_window_contains` | Half-open [lower, upper) convention |
| `test_neutral_loss_mass` | NeutralLoss H₂O and NH₃ masses |
| `test_ion_type_from_char` | IonType::from_char parses b/y/z ions |
| `test_ms1_spectrum_find_peak` | find_peak_ppm returns most intense within tolerance |
| `test_isotope_envelope_mz_calculation` | Isotope m/z spacing for charge 2 |
| `test_isotope_envelope_extraction` | Envelope extracts intensities and M+0 mass error |
| `test_isotope_envelope_missing_peaks` | Missing peaks reported as zero |

#### Configuration — [crates/osprey-core/src/config.rs](../crates/osprey-core/src/config.rs)

| Test | Description |
|---|---|
| `test_bin_config_unit_resolution` | Unit resolution BinConfig creates correct bins |
| `test_library_source_from_path` | Auto-detects format from file extension |
| `test_default_config` | Default OspreyConfig has expected values |
| `test_yaml_roundtrip` | YAML serialization preserves all fields |
| `test_config_merge` | CLI overrides correctly merge with config |
| `test_rt_calibration_config` | RTCalibrationConfig defaults and disabled() |
| `test_fragment_tolerance_ppm` | PPM tolerance varies correctly with m/z |
| `test_fragment_tolerance_da` | Da tolerance constant across m/z |

#### Error Handling — [crates/osprey-core/src/error.rs](../crates/osprey-core/src/error.rs)

| Test | Description |
|---|---|
| `test_error_display` | Error Display impl formats messages correctly |
| `test_error_helpers` | Helper constructors produce correct variants |

#### Traits — [crates/osprey-core/src/traits.rs](../crates/osprey-core/src/traits.rs)

| Test | Description |
|---|---|
| `test_mock_loader` | LibraryLoader::supports_format matches extension |

---

### osprey-io (26 tests)

File I/O for spectral libraries and mass spec data.

#### DIA-NN TSV — [crates/osprey-io/src/library/diann.rs](../crates/osprey-io/src/library/diann.rs)

| Test | Description |
|---|---|
| `test_strip_modifications` | Bracketed and parenthesized modifications stripped |
| `test_parse_modifications` | Mass shift modifications parsed correctly |
| `test_split_list` | Semicolon/comma-delimited lists split correctly |
| `test_parse_mod_mass` | Signed numeric and named modifications parsed |

#### BiblioSpec blib Input — [crates/osprey-io/src/library/blib.rs](../crates/osprey-io/src/library/blib.rs)

| Test | Description |
|---|---|
| `test_parse_blib_modifications_simple` | Unmodified peptide produces empty list |
| `test_parse_blib_modifications_carbamidomethyl` | Carbamidomethyl parsed correctly |
| `test_parse_blib_modifications_oxidation` | Oxidation modification parsed correctly |
| `test_identify_modification` | Mass delta identified as correct modification |

#### EncyclopeDIA elib — [crates/osprey-io/src/library/elib.rs](../crates/osprey-io/src/library/elib.rs)

| Test | Description |
|---|---|
| `test_parse_modified_sequence_simple` | Unmodified sequence parsed correctly |
| `test_parse_modified_sequence_with_mod` | Bracketed modification extracted |
| `test_parse_modified_sequence_nterm` | N-terminal modification at position 0 |
| `test_parse_modified_sequence_multiple` | Multiple modifications parsed |

#### blib Output — [crates/osprey-io/src/output/blib.rs](../crates/osprey-io/src/output/blib.rs)

| Test | Description |
|---|---|
| `test_create_blib` | blib created with spectrum, boundaries, metadata |
| `test_blob_compression` | compress_bytes produces valid zlib or raw fallback |
| `test_modifications_1based` | Modifications written with 1-based positions |
| `test_strip_flanking_chars` | Flanking amino acids/underscores/dots stripped |
| `test_convert_unimod_to_mass` | UniMod notation converted to mass shift |
| `test_fragment_roundtrip` | Fragment m/z and intensity survive write-read |
| `test_modseq_roundtrip` | UniMod converted to mass notation in database |
| `test_protein_mapping_roundtrip` | Protein mappings stored and retrievable |
| `test_retention_times_multirun` | Multi-run RT with correct Skyline flags |

#### Library Deduplication — [crates/osprey-io/src/library/mod.rs](../crates/osprey-io/src/library/mod.rs)

| Test | Description |
|---|---|
| `test_deduplicate_no_duplicates` | No changes when no duplicates present |
| `test_deduplicate_removes_duplicates` | Duplicates merged (best kept, RT averaged) |
| `test_deduplicate_different_charges_are_separate` | Different charges kept separate |
| `test_deduplicate_sequential_ids` | IDs renumbered sequentially |

#### mzML Parsing — [crates/osprey-io/src/mzml/parser.rs](../crates/osprey-io/src/mzml/parser.rs)

| Test | Description |
|---|---|
| `test_isolation_window` | Symmetric window has correct center and width |

---

### osprey-chromatography (66 tests)

Peak detection, CWT, calibration.

#### CWT Peak Detection — [crates/osprey-chromatography/src/cwt.rs](../crates/osprey-chromatography/src/cwt.rs)

| Test | Description |
|---|---|
| `test_kernel_zero_mean` | Mexican hat kernel is zero-mean |
| `test_kernel_symmetric` | Kernel is symmetric |
| `test_kernel_positive_center_negative_tails` | Center positive, tails negative |
| `test_kernel_size` | Kernel has correct size (2×width+1) |
| `test_convolve_same_length` | Convolution result matches signal length |
| `test_convolve_delta_function` | Delta convolved with kernel produces kernel |
| `test_convolve_gaussian_response` | CWT of Gaussian produces positive response |
| `test_estimate_scale_known_peak` | Estimated sigma near 5.0 for known peak |
| `test_estimate_scale_fallback` | Fallback to 4.0 for all-zero XICs |
| `test_consensus_single_gaussian_peak` | Single Gaussian peak detected correctly |
| `test_consensus_interference_rejection` | Interfering fragment rejected via median |
| `test_consensus_two_separated_peaks` | Two separated peaks detected |
| `test_consensus_degenerate_single_xic` | Empty for single XIC |
| `test_consensus_degenerate_short_xic` | Empty for too-short XICs |
| `test_consensus_all_zero` | Empty for all-zero XICs |
| `test_consensus_peak_bounds_valid` | Boundaries properly ordered |
| `test_consensus_noise_robustness` | Peak detected despite noise |

#### Peak Detection (classical) — [crates/osprey-chromatography/src/lib.rs](../crates/osprey-chromatography/src/lib.rs)

| Test | Description |
|---|---|
| `test_peak_detector` | Single Gaussian peak detected with correct apex |
| `test_find_best_peak` | Highest-coefficient peak selected within window |
| `test_fwhm_cap_symmetric_peak` | FWHM capping doesn't change symmetric peak |
| `test_fwhm_cap_tailing_peak` | FWHM capping tightens right boundary |
| `test_fwhm_cap_preserves_valley` | Valley boundary preserved between peaks |
| `test_fwhm_cap_fallback` | Graceful fallback when half-height not found |
| `test_asymmetric_half_widths` | Tailing peak shows larger right half-width |
| `test_detect_all_xic_peaks_finds_multiple` | Multiple candidates sorted by intensity |
| `test_detect_all_xic_peaks_single_peak` | Single peak for clean Gaussian |
| `test_detect_all_xic_peaks_noise_plus_real` | Tallest candidate over noise bump |

#### EMG Peak Fitting — [crates/osprey-chromatography/src/lib.rs](../crates/osprey-chromatography/src/lib.rs)

| Test | Description |
|---|---|
| `test_erfc_approx` | Approximation matches known values |
| `test_normal_cdf` | Normal CDF matches known values |
| `test_emg_pdf_basic` | EMG produces right-tailing peak |
| `test_emg_cdf_monotonic` | EMG CDF monotonically increasing |
| `test_emg_fit_gaussian` | EMG fit to Gaussian recovers parameters |
| `test_emg_fit_tailing_peak` | EMG fit to tailing recovers parameters |
| `test_emg_boundaries` | 95% boundaries capture specified area |
| `test_emg_fit_and_boundaries` | Full pipeline fits and verifies boundaries |

#### Mass Calibration — [crates/osprey-chromatography/src/calibration/mass.rs](../crates/osprey-chromatography/src/calibration/mass.rs)

| Test | Description |
|---|---|
| `test_ppm_error_calculation` | PPM error calculated correctly |
| `test_apply_calibration_negative_offset` | Negative offset shifts m/z upward |
| `test_apply_calibration_positive_offset` | Positive offset shifts m/z downward |
| `test_calculate_mz_calibration` | MS1/MS2 statistics computed correctly |
| `test_empty_qc_data` | Empty data produces uncalibrated parameters |
| `test_uncalibrated_passthrough` | Uncalibrated m/z returned unchanged |
| `test_within_calibrated_tolerance` | Tolerance accepts/rejects correctly |
| `test_histogram_generation` | Histogram bin structure and counts correct |
| `test_median_calculation` | Median for odd and even lengths |
| `test_unit_resolution_calibration` | Unit resolution produces correct labels |

#### RT Calibration — [crates/osprey-chromatography/src/calibration/rt.rs](../crates/osprey-chromatography/src/calibration/rt.rs)

| Test | Description |
|---|---|
| `test_rt_calibration_linear` | LOESS fits linear relationship with high R² |
| `test_rt_calibration_nonlinear` | LOESS captures sinusoidal nonlinearity |
| `test_stratified_sampler` | Samples distributed evenly across bins |
| `test_calibration_edge_cases` | Succeeds with minimum points, fails with fewer |
| `test_median` | Median computes correctly |
| `test_std_dev` | Sample standard deviation matches reference |
| `test_local_tolerance_varies_with_rt` | Tolerances vary with RT |
| `test_local_tolerance_minimum_floor` | Respects minimum floor |
| `test_local_tolerance_extrapolation` | Valid outside calibration range |
| `test_model_params_roundtrip_with_abs_residuals` | Export/import preserves residuals |
| `test_backwards_compatibility_no_abs_residuals` | Loads old format without residuals |
| `test_predict_with_duplicate_library_rts` | No NaN for duplicate points |
| `test_inverse_predict_roundtrip_linear` | Inverse predict roundtrips linear |
| `test_inverse_predict_roundtrip_nonlinear` | Inverse predict roundtrips nonlinear |
| `test_inverse_predict_duplicate_fitted` | Handles duplicate fitted values |

#### Calibration I/O — [crates/osprey-chromatography/src/calibration/io.rs](../crates/osprey-chromatography/src/calibration/io.rs)

| Test | Description |
|---|---|
| `test_calibration_filename` | Appends ".calibration.json" correctly |
| `test_calibration_filename_for_input` | Derives filename from mzML paths |
| `test_calibration_path_for_input` | Joins output directory correctly |

#### Calibration Config — [crates/osprey-chromatography/src/calibration/mod.rs](../crates/osprey-chromatography/src/calibration/mod.rs)

| Test | Description |
|---|---|
| `test_uncalibrated_params` | Uncalibrated reports all flags false |
| `test_effective_tolerance` | Returns base or adjusted tolerance |

---

### osprey-scoring (70 tests)

Spectral scoring, decoy generation, coelution features.

#### Spectral Scoring — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

| Test | Description |
|---|---|
| `test_spectral_scorer_lib_cosine` | LibCosine ~1.0 for perfect match |
| `test_spectral_scorer_partial_match` | Partial match reports correct coverage |
| `test_spectral_scorer_no_match` | Zero score for no overlapping peaks |
| `test_spectral_scorer_xcorr` | XCorr non-zero for matching spectra |

#### XCorr — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

| Test | Description |
|---|---|
| `test_xcorr_perfect_vs_partial_match` | Full match scores higher than partial |
| `test_xcorr_no_match_is_low` | Non-overlapping fragments produce near-zero XCorr |
| `test_xcorr_scaling_factor` | 0.005 scaling produces reasonable magnitude |
| `test_xcorr_empty_inputs` | Zero for empty spectrum/library |
| `test_xcorr_preprocessed_matches_direct` | Batch preprocessing matches one-shot path |

#### Pearson Correlation — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

Tests for `pearson_correlation_raw()`, the core function used for fragment co-elution scoring.

| Test | Description |
|---|---|
| `test_pearson_identical_vectors` | r = 1.0 for identical vectors |
| `test_pearson_perfect_negative` | r = -1.0 for anti-correlated |
| `test_pearson_uncorrelated` | r ≈ 0 for uncorrelated |
| `test_pearson_too_short` | Returns 0.0 for < 3 points |
| `test_pearson_constant_input` | Finite output for zero-variance input |
| `test_pearson_known_value` | Known reference (r ≈ 0.9820 for [1,2,3] vs [2,4,5]) |
| `test_pearson_linear_transform` | Invariant to positive linear scaling |

#### Coelution Sum (Pairwise Fragment Correlation) — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

Tests for the pairwise sum of Pearson correlations used to score candidate peaks. This is how peaks are ranked in the main search pipeline.

| Test | Description |
|---|---|
| `test_coelution_sum_perfect_coelution` | C(n,2) pairs × r=1.0 for identical XICs |
| `test_coelution_sum_scaled_fragments` | Scaled profiles still r=1.0 |
| `test_coelution_sum_one_interferer` | Interference reduces score |
| `test_coelution_sum_two_fragments` | Single pair equals pairwise r |
| `test_coelution_sum_noise_reduces_score` | Noise drops correlations below 1.0 |
| `test_coelution_sum_peak_selection` | Co-eluting peak selected over interfered peak |

#### Mean Pairwise Correlation (Pipeline Peak Selection) — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

Tests for the mean pairwise correlation metric used in the pipeline to choose the best peak among candidates.

| Test | Description |
|---|---|
| `test_mean_pairwise_correlation_peak_selection` | Best co-eluting peak selected |
| `test_mean_pairwise_correlation_two_fragments` | Mean equals single Pearson for 2 fragments |
| `test_mean_pairwise_correlation_degenerate` | Handles 0–1 fragment edge cases |

#### Decoy Generation — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

| Test | Description |
|---|---|
| `test_decoy_reverse_trypsin` | Trypsin-aware reversal preserves C-terminus |
| `test_decoy_reverse_lysn` | LysN-aware reversal preserves N-terminus |
| `test_decoy_preserves_precursor_mz` | Same amino acids preserve m/z |
| `test_decoy_protein_id_prefix` | Protein IDs prefixed "DECOY_" |
| `test_decoy_with_modification` | Modification positions remapped |
| `test_short_sequence` | Short sequences returned unchanged |
| `test_enzyme_detection` | Enzyme cleavage types detected correctly |
| `test_decoy_collision_detection` | Collision prevention works |
| `test_cycle_sequence` | Cyclic permutation shifts correctly |

#### Spectrum Aggregation — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

| Test | Description |
|---|---|
| `test_spectrum_aggregator` | Intensities summed across spectra |
| `test_spectrum_aggregator_weighted` | Coefficient weights applied correctly |

#### Fragment Matching — [crates/osprey-scoring/src/lib.rs](../crates/osprey-scoring/src/lib.rs)

| Test | Description |
|---|---|
| `test_get_top_n_fragment_indices` | Top N fragments selected by intensity |
| `test_has_match_within_tolerance_ppm` | PPM matching works |
| `test_has_match_within_tolerance_mz` | Da matching works |
| `test_has_match_empty_spectrum` | No match for empty spectrum |

#### Calibration ML (LDA) — [crates/osprey-scoring/src/calibration_ml.rs](../crates/osprey-scoring/src/calibration_ml.rs)

| Test | Description |
|---|---|
| `test_feature_matrix_extraction` | Matrix extracted with correct dimensions |
| `test_lda_training` | Targets score higher than decoys after LDA |
| `test_compete_calibration_pairs_target_wins` | Target wins when higher score |
| `test_compete_calibration_pairs_decoy_wins` | Decoy wins when higher score |
| `test_compete_calibration_pairs_tie_goes_to_decoy` | Tie conservative to decoy |
| `test_compete_calibration_pairs_singleton_auto_wins` | Unpaired entries auto-win |
| `test_compete_calibration_pairs_multiple_pairs` | Multiple pairs compete correctly |
| `test_compete_calibration_pairs_empty_input` | Empty input handled |
| `test_compete_calibration_pairs_subset` | Subset only considers specified |
| `test_compete_calibration_pairs_deterministic_with_ties` | Deterministic with ties |

#### Batch Scoring — [crates/osprey-scoring/src/batch.rs](../crates/osprey-scoring/src/batch.rs)

| Test | Description |
|---|---|
| `test_preprocess_library` | Library indexed and mapped correctly |
| `test_preprocess_spectra` | Spectra stored with retention times |
| `test_batch_scoring_perfect_match` | Score ~1.0 for identical peaks |
| `test_batch_scoring_no_match` | Score ~0.0 for no overlap |
| `test_find_best_matches` | Each library matched to best spectrum |
| `test_library_subset` | Subset retains only requested entries |
| `test_libcosine_perfect_match` | Score 1.0 for exact PPM match |
| `test_libcosine_no_match` | Score 0.0 for no PPM match |
| `test_libcosine_ppm_matching` | Matches within PPM, rejects outside |
| `test_libcosine_mass_errors` | Per-fragment errors accurate |
| `test_libcosine_da_tolerance` | Da tolerance matching works |
| `test_libcosine_multiple_fragments` | Multiple fragments scored correctly |
| `test_libcosine_closest_mz_selected` | Closest m/z selected over higher intensity |
| `test_count_top6_matched_at_apex` | Top-6 fragments counted correctly |
| `test_group_spectra_by_isolation_window_sorted` | Sorted isolation windows |

#### Preprocessing Pipeline — [crates/osprey-scoring/src/pipeline/](../crates/osprey-scoring/src/pipeline/)

| Test | File | Description |
|---|---|---|
| `test_preprocessing_worker` | preprocessor.rs | Metadata preserved, XCorr vector correct size |
| `test_empty_spectrum` | preprocessor.rs | Empty spectrum handled gracefully |
| `test_accumulator_grouping` | accumulator.rs | Spectra grouped by window |
| `test_drain_all` | accumulator.rs | All windows returned and emptied |
| `test_window_key` | accumulator.rs | Discretized bounds treated as equal |

---

### osprey-fdr (41 tests)

FDR control, target-decoy competition, Percolator, Mokapot.

#### FDR Controller — [crates/osprey-fdr/src/lib.rs](../crates/osprey-fdr/src/lib.rs)

| Test | Description |
|---|---|
| `test_qvalue_computation` | Q-values computed from target/decoy distributions |
| `test_filter_by_qvalue` | Items filtered at threshold |
| `test_count_at_thresholds` | Counts at 0.1%, 1%, 5%, 10% FDR |

#### Target-Decoy Competition — [crates/osprey-fdr/src/lib.rs](../crates/osprey-fdr/src/lib.rs)

| Test | Description |
|---|---|
| `test_competition_target_wins` | Target wins when higher score |
| `test_competition_decoy_wins` | Decoy wins when higher score |
| `test_competition_tie_goes_to_decoy` | Tied scores go to decoy (conservative) |
| `test_competition_fdr_calculation` | Cumulative FDR = decoy_wins/target_wins |
| `test_competition_multiple_pairs_all_targets_win` | All targets pass at 0% FDR |
| `test_competition_all_decoys_win` | No targets pass |
| `test_competition_keeps_best_score_per_peptide` | Best score per entry_id retained |
| `test_competition_empty_input` | Empty input handled |
| `test_competition_target_without_decoy` | Unpaired target wins by default |
| `test_competition_fdr_recovers_after_spike` | FDR recovery: max cumulative targets at valid FDR |

#### Percolator (SVM-based semi-supervised FDR) — [crates/osprey-fdr/src/percolator.rs](../crates/osprey-fdr/src/percolator.rs)

| Test | Description |
|---|---|
| `test_fold_assignment_peptide_grouping` | Same peptide targets share fold |
| `test_conservative_qvalues` | Conservative q-values with +1 correction |
| `test_compete_and_count` | Target-decoy competition counting |
| `test_percolator_basic` | Full Percolator run: targets score > decoys |
| `test_best_precursor_per_peptide` | Best score per peptide selected |
| `test_percolator_empty` | Empty input handled |
| `test_compete_from_indices_deterministic_with_ties` | Deterministic ordering with tied scores |

#### Multi-Level FDR — [crates/osprey-fdr/src/percolator.rs](../crates/osprey-fdr/src/percolator.rs)

Tests for run-level and experiment-level q-value computation, and the dual precursor + peptide FDR.

| Test | Description |
|---|---|
| `test_per_run_precursor_qvalues_independent_files` | Per-file q-values don't leak across files |
| `test_experiment_precursor_qvalues_cross_file` | Experiment-level picks best across files |
| `test_experiment_peptide_qvalues_aggregates_by_peptide` | Same q-value for all charge states |
| `test_effective_qvalue_is_max_of_precursor_and_peptide` | Dual FDR: `max(precursor_q, peptide_q)` |
| `test_per_run_peptide_qvalues_propagate_to_charges` | Peptide q-values propagate to charge states |

#### Cross-Validation Fold Integrity — [crates/osprey-fdr/src/percolator.rs](../crates/osprey-fdr/src/percolator.rs)

These tests ensure the [critical invariant](../CLAUDE.md) that target-decoy pairs, charge states, and peptide groups are never split across CV folds or subsamples.

| Test | Description |
|---|---|
| `test_subsample_keeps_target_decoy_pairs` | Target-decoy pairs never split during subsampling |
| `test_subsample_keeps_charge_states_together` | Multi-charge peptides stay together |
| `test_fold_assignment_multi_charge_with_decoys` | Multi-charge + decoy fold assignment |
| `test_subsample_no_reduction_when_under_limit` | No-op when under max_entries limit |

#### Mokapot PIN Format — [crates/osprey-fdr/src/mokapot.rs](../crates/osprey-fdr/src/mokapot.rs)

| Test | Description |
|---|---|
| `test_format_peptide` | Peptide wrapped in `-.PEPTIDE.-` format |
| `test_strip_flanking_chars_preserves_modification_dots` | Dots inside `[+57.02146]` preserved |
| `test_feature_header` | Expected feature names present, removed features absent |
| `test_charge_features` | One-hot encoding for charges 1–5 |
| `test_mokapot_runner_availability` | Availability check doesn't panic |
| `test_get_pin_feature_names_count` | Feature count matches header |
| `test_pin_feature_value_matches_format` | Values consistent with format |

#### Mokapot Result Parsing — [crates/osprey-fdr/src/mokapot.rs](../crates/osprey-fdr/src/mokapot.rs)

Tests `parse_results()` with mock TSV files covering the various header variants mokapot can produce.

| Test | Description |
|---|---|
| `test_parse_results_standard_headers` | Standard `PSMId/score/q-value/posterior_error_prob` |
| `test_parse_results_alternative_headers` | Alternative `specid/mokapot_score/mokapot_qvalue` |
| `test_parse_results_psmid_header` | `psmid` and `mokapot q-value` variants |
| `test_parse_results_empty_file` | Header-only file returns empty results |
| `test_parse_results_skips_malformed_lines` | Truncated lines skipped gracefully |

---

### osprey-ml (30 tests)

Machine learning primitives: SVM, LDA, PEP, q-values.

#### Linear Discriminant Analysis — [crates/osprey-ml/src/linear_discriminant.rs](../crates/osprey-ml/src/linear_discriminant.rs)

| Test | Description |
|---|---|
| `linear_discriminant` | Power method computes eigenvector correctly |

#### Posterior Error Probability — [crates/osprey-ml/src/pep.rs](../crates/osprey-ml/src/pep.rs)

| Test | Description |
|---|---|
| `test_isotonic_regression_already_decreasing` | Already decreasing unchanged |
| `test_isotonic_regression_single_violation` | Single violation fixed via PAVA |
| `test_isotonic_regression_all_increasing` | All increasing averaged |
| `test_isotonic_regression_empty_and_single` | Edge cases handled |
| `test_pep_well_separated` | High scores → low PEP, low scores → high PEP |
| `test_pep_monotonicity` | PEP monotonically non-increasing with score |
| `test_pep_empty_input` | Empty input returns 1.0 |
| `test_pep_all_targets` | All targets returns 1.0 |

#### Matrix Operations — [crates/osprey-ml/src/matrix.rs](../crates/osprey-ml/src/matrix.rs)

| Test | Description |
|---|---|
| `dotv` | Matrix-vector dot product |
| `tranpose` | Matrix transpose |
| `dot` | Matrix-matrix dot product |
| `slice` | Row slice extraction |

#### SVM — [crates/osprey-ml/src/svm.rs](../crates/osprey-ml/src/svm.rs)

| Test | Description |
|---|---|
| `test_linearly_separable` | Perfectly separated by hyperplane |
| `test_overlapping_classes` | Overlapping classes still discriminated |
| `test_weights_direction` | Discriminative feature weighted heavily |
| `test_empty_input` | Empty data handled |
| `test_deterministic_with_seed` | Same seed produces same model |
| `test_decision_function` | Decision function computes correctly |
| `test_feature_standardizer` | Mean/std computed correctly |
| `test_standardizer_zero_variance` | Zero variance handled |
| `test_grid_search_c` | Grid search selects reasonable C |
| `test_count_passing_targets` | Passing targets counted correctly |
| `test_count_passing_targets_clean` | Clean data all pass at FDR |
| `test_xorshift_deterministic` | RNG deterministic with seed |
| `test_shuffle_deterministic` | Fisher-Yates deterministic |
| `test_dual_cd_convergence` | Dual coordinate descent converges |
| `test_count_passing_targets_svm_deterministic` | Deterministic with tied scores |

#### Q-Value Computation — [crates/osprey-ml/src/qvalue.rs](../crates/osprey-ml/src/qvalue.rs)

| Test | Description |
|---|---|
| `test_qvalue_calculation` | Q-values computed from target/decoy mixture |
| `test_qvalue_mixed` | Monotonically non-increasing property |

---

### osprey (main crate, 26 tests)

Pipeline logic, candidate selection, multi-charge consensus, cross-run reconciliation.

#### Candidate Selection — [crates/osprey/src/pipeline.rs](../crates/osprey/src/pipeline.rs)

| Test | Description |
|---|---|
| `test_select_candidates` | Candidates within isolation window and RT selected |
| `test_build_mz_index` | Entries indexed in ±1 Da bins |
| `test_build_mz_index_empty` | Empty library handled |
| `test_select_candidates_max_limit` | Max candidates limit respected |

#### Fragment Overlap — [crates/osprey/src/pipeline.rs](../crates/osprey/src/pipeline.rs)

| Test | Description |
|---|---|
| `test_fragment_overlap_identical` | Identical fragments: 6/6 overlap |
| `test_fragment_overlap_disjoint` | Disjoint fragments: 0 overlap |
| `test_fragment_overlap_ppm_tolerance` | PPM tolerance used correctly |
| `test_top_n_fragments` | Top N fragments returned correctly |

#### DIA Isolation Scheme — [crates/osprey/src/pipeline.rs](../crates/osprey/src/pipeline.rs)

| Test | Description |
|---|---|
| `test_extract_isolation_scheme` | DIA window cycle detected |
| `test_extract_isolation_scheme_empty` | Empty spectra handled |

#### Multi-Charge Consensus — [crates/osprey/src/pipeline.rs](../crates/osprey/src/pipeline.rs)

Tests for the post-FDR consensus phase where charge states of the same peptide are reconciled to a single peak.

| Test | Description |
|---|---|
| `test_deduplicate_pairs_deterministic` | Deterministic order by entry_id |
| `test_consensus_single_charge_no_change` | Single charge not re-scored |
| `test_consensus_two_charges_same_peak` | Same peak: no re-scoring needed |
| `test_consensus_two_charges_different_peaks` | Different peaks: one re-scored |
| `test_consensus_uses_svm_score_not_coelution` | SVM score preferred |
| `test_consensus_skips_groups_where_none_pass_fdr` | Non-passing groups skipped |
| `test_consensus_three_charges_two_agree` | Multi-charge consensus |
| `test_consensus_decoys_separate_from_targets` | Separate target/decoy groups |
| `test_consensus_multiple_peptides_independent` | Independent per-peptide processing |
| `test_consensus_groups_by_file_name` | File-based grouping |

#### Cross-Run Reconciliation — [crates/osprey/src/reconciliation.rs](../crates/osprey/src/reconciliation.rs)

| Test | Description |
|---|---|
| `test_weighted_median_single` | Single value returns that value |
| `test_weighted_median_equal_weights` | Equal weights give regular median |
| `test_weighted_median_skewed_weights` | Heavy weight pulls median |
| `test_simple_median_odd` | Odd-length median correct |
| `test_simple_median_even` | Even-length median correct |
| `test_weighted_median_empty` | Empty returns 0.0 |

---

## Test Design Principles

### All tests are unit tests
Every test is co-located with the code it tests in `#[cfg(test)] mod tests` blocks. This keeps tests close to the implementation and avoids external test data dependencies.

### No external data required
Tests use synthetically constructed inputs (spectra, library entries, feature vectors). The `example_test_data/` directory contains real mass spec files for manual validation but is not used by `cargo test`.

### Critical invariants have dedicated tests
The [cross-validation grouping invariant](../CLAUDE.md#cross-validation-grouping) — target-decoy pairs, charge states, and peptide groups must never be split across folds — is tested by multiple dedicated tests in percolator.rs (`test_subsample_keeps_target_decoy_pairs`, `test_subsample_keeps_charge_states_together`, `test_fold_assignment_*`).

### Determinism is verified
Several tests run operations multiple times to verify deterministic output, especially for code involving HashMaps or tied scores (`test_compete_from_indices_deterministic_with_ties`, `test_deterministic_with_seed`).

### The mokapot tests use mock files
Mokapot result parsing tests write temporary TSV files with known content and verify the parser handles all known header variants. Mokapot itself (the Python tool) is not required to run the test suite.
