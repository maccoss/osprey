# Osprey Release Notes (Next Release)

Working draft for the next release. Append entries here as features and fixes land on the development branch. At release time this file is renamed to `RELEASE_NOTES_v{version}.md` and the workspace `Cargo.toml` version is updated to match.

## New Features

- **Widened RT penalty sigma from 3x to 5x calibration MAD.** The 3-sigma Gaussian RT penalty introduced in v26.3.0 was too aggressive for peptides whose true elution time deviates slightly from the LOESS prediction. A peptide with its correct peak 0.3 min from the predicted RT could lose to a shoulder at the predicted position because the 10% RT penalty (at 3-sigma) was enough to flip the selection when coelution scores were comparable. At 5-sigma, the penalty at 0.3 min is ~4% (essentially negligible), while the penalty at 1.0 min is still ~37% (substantial for genuine wrong-peak selections). This preserves the interferer-rejection benefit while avoiding over-penalization of peaks with slight RT deviations from the calibration.

- **Added intensity tiebreaker to CWT peak selection.** When two CWT candidates from the same chromatographic peak (main peak vs shoulder) have nearly identical coelution scores, the more intense peak now wins. The selection score is multiplied by `log(1 + apex_intensity)`, which keeps intensity as a secondary factor that breaks ties without dominating the coelution-based ranking. This fixes cases where a narrow shoulder near the calibration-predicted RT beat the main peak because mean pairwise correlation is similar for both (they share the same eluting species) and the RT penalty gave the shoulder a slight edge.

## Bug Fixes

<!-- none yet -->

## Performance

<!-- none yet -->

## Breaking Changes

<!-- none yet -->
