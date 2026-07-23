//! Learned per-platform peak-pick model.
//!
//! Ported from the C# `PickLdaModel` (`pwiz_tools/Osprey/Osprey.Scoring/PickLdaModel.cs`)
//! for cross-implementation parity. When active it REPLACES the product-form CWT
//! candidate rank score in `run_search` (`pipeline.rs`) with a standardized linear
//! combination of the same four raw terms the `OSPREY_PICK_DUMP_CANDIDATES` dump
//! captures:
//!
//! ```text
//! rank = w0*z(coelution) + w1*z(ln_intensity) + w2*z(rt_penalty) + w3*z(median_polish)
//! ```
//!
//! where `z(x_i) = (x_i - mean[i]) / scale[i]`. The argmax + IEEE-754 total-order
//! tie-break in the pick loop are unchanged.
//!
//! Selection precedence (matches C#), resolved once at first use:
//!   1. `OSPREY_PICK_LDA_MODEL` (path to a JSON model) → that model (test override);
//!   2. else `OSPREY_PICK_LDA` set and not `"0"` → the hardcoded resolution-keyed model
//!      (Stellar for unit, Astral for HRAM);
//!   3. else (DEFAULT) → the pure product pick (`None`; Rust parity + regression golden).

use std::sync::OnceLock;

use serde::Deserialize;

/// Number of pick features: coelution, ln_intensity, rt_penalty, median_polish.
const N: usize = 4;

/// The fixed feature order `score()` applies weights to. A loaded model MUST list these
/// names in this exact order (see `load_from_json`) so a re-ordered or older-schema file
/// cannot silently map weights to the wrong term. Matches the C# `ExpectedFeatures` and
/// pick_lda_train.py's FEATURES.
const EXPECTED_FEATURES: [&str; N] = ["coelution", "ln_intensity", "rt_penalty", "median_polish"];

/// A frozen linear pick model: per-feature weights over standardized (z-scored) terms.
#[derive(Clone, Debug)]
pub struct PickLdaModel {
    weights: [f64; N],
    means: [f64; N],
    scales: [f64; N],
}

impl PickLdaModel {
    /// The learned rank score for one candidate from its four raw terms. A zero scale
    /// standardizes to 0 for that feature (a constant feature contributes nothing)
    /// rather than dividing by zero. Feature order is fixed and MUST match the C#
    /// `Score` call and the dump columns: coelution, ln_intensity, rt_penalty,
    /// median_polish.
    pub fn score(
        &self,
        coelution: f64,
        ln_intensity: f64,
        rt_penalty: f64,
        median_polish: f64,
    ) -> f64 {
        self.weights[0] * self.z(coelution, 0)
            + self.weights[1] * self.z(ln_intensity, 1)
            + self.weights[2] * self.z(rt_penalty, 2)
            + self.weights[3] * self.z(median_polish, 3)
    }

    fn z(&self, x: f64, i: usize) -> f64 {
        if self.scales[i] != 0.0 {
            (x - self.means[i]) / self.scales[i]
        } else {
            0.0
        }
    }
}

// Resolution-keyed default peak-pick models, learned per platform via the paired
// target/decoy pick-LDA (see pick_lda_train.py / OSPREY_PICK_DUMP_CANDIDATES). Unit
// resolution -> Stellar-trained weights; HRAM -> Astral-trained weights. Values copied
// verbatim from the C# PickLdaModel StellarModel / AstralModel (pick-model-stellar.json
// / pick-model-astral.json), so the two implementations rank candidates identically.
const STELLAR_MODEL: PickLdaModel = PickLdaModel {
    weights: [
        0.9933168416485256,
        0.047052481253413006,
        0.027130393118192445,
        0.10184133676728513,
    ],
    means: [
        0.14931086687377143,
        9.15749607304815,
        0.9158212545538758,
        0.8904037620854307,
    ],
    scales: [
        0.2074610054197197,
        2.260347450504608,
        0.07333900267211861,
        0.08220791724094854,
    ],
};

const ASTRAL_MODEL: PickLdaModel = PickLdaModel {
    weights: [
        0.5348241578558818,
        0.0041302671426268105,
        0.3352868625222239,
        0.7755828652613985,
    ],
    means: [
        0.027393438120134818,
        6.585876043601798,
        0.939316453828307,
        0.6880222328717774,
    ],
    scales: [
        0.11714825645722571,
        3.9104476306002494,
        0.05338768554968953,
        0.1461117956225306,
    ],
};

/// The process-wide pick mode, resolved once from the environment.
enum PickMode {
    /// `OSPREY_PICK_LDA_MODEL` override — applies to every resolution.
    EnvModel(PickLdaModel),
    /// `OSPREY_PICK_LDA` — the hardcoded resolution-keyed model (opt-in).
    ResolutionModel,
    /// Default: the pure product-form pick (no model).
    Legacy,
}

static PICK_MODE: OnceLock<PickMode> = OnceLock::new();

fn pick_mode() -> &'static PickMode {
    PICK_MODE.get_or_init(|| {
        // (1) OSPREY_PICK_LDA_MODEL: a non-empty path to an existing JSON file wins.
        // Matches C# LoadFromEnv: empty / missing file → fall through (null); a present
        // but unparseable file is a hard error.
        if let Ok(path) = std::env::var("OSPREY_PICK_LDA_MODEL") {
            if !path.is_empty() && std::path::Path::new(&path).exists() {
                return PickMode::EnvModel(load_from_json(&path));
            }
        }
        // (2) OSPREY_PICK_LDA: set and not "0" → the resolution-keyed built-in model (opt-in).
        if let Ok(v) = std::env::var("OSPREY_PICK_LDA") {
            if !v.is_empty() && v != "0" {
                return PickMode::ResolutionModel;
            }
        }
        // (3) Default: the legacy pure product-form pick (Rust parity + regression golden).
        PickMode::Legacy
    })
}

/// The active pick model for the given resolution, or `None` when the legacy
/// product-form pick is in effect. `is_hram` mirrors the C# `Resolution.HasMs1Features`
/// (HRAM → Astral-trained, unit resolution → Stellar-trained).
pub fn active_model(is_hram: bool) -> Option<&'static PickLdaModel> {
    match pick_mode() {
        PickMode::EnvModel(m) => Some(m),
        PickMode::ResolutionModel => Some(if is_hram {
            &ASTRAL_MODEL
        } else {
            &STELLAR_MODEL
        }),
        PickMode::Legacy => None,
    }
}

/// JSON schema for `OSPREY_PICK_LDA_MODEL`, matching the C# `Dto`:
/// `{ "features": [...], "weights": [w..], "means": [m..], "scales": [s..] }`.
#[derive(Deserialize)]
struct Dto {
    #[serde(default)]
    features: Vec<String>,
    weights: Vec<f64>,
    means: Vec<f64>,
    scales: Vec<f64>,
}

fn load_from_json(path: &str) -> PickLdaModel {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("OSPREY_PICK_LDA_MODEL: failed to read '{path}': {e}"));
    let dto: Dto = serde_json::from_str(&text).unwrap_or_else(|e| {
        panic!("OSPREY_PICK_LDA_MODEL: failed to parse pick model JSON at '{path}': {e}")
    });
    if dto.weights.len() != N || dto.means.len() != N || dto.scales.len() != N {
        panic!(
            "OSPREY_PICK_LDA_MODEL: pick model JSON at '{path}' must define weights, means, \
             and scales as length-{N} arrays."
        );
    }
    // The weights are positional (score() applies weights[i] to the i-th raw term), so a file
    // whose features are in a different order -- or absent -- would silently score the wrong
    // term. Require the names and the exact expected order rather than trusting array position.
    if dto.features.len() != N
        || dto
            .features
            .iter()
            .zip(EXPECTED_FEATURES)
            .any(|(got, want)| got.as_str() != want)
    {
        panic!(
            "OSPREY_PICK_LDA_MODEL: pick model JSON at '{path}' must list features as \
             {EXPECTED_FEATURES:?} in that exact order; got {:?}.",
            dto.features
        );
    }
    let to_arr = |v: Vec<f64>| -> [f64; N] {
        let mut a = [0.0; N];
        a.copy_from_slice(&v);
        a
    };
    PickLdaModel {
        weights: to_arr(dto.weights),
        means: to_arr(dto.means),
        scales: to_arr(dto.scales),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_scale_feature_contributes_nothing() {
        let m = PickLdaModel {
            weights: [1.0, 1.0, 1.0, 1.0],
            means: [0.0, 0.0, 0.0, 0.0],
            scales: [1.0, 1.0, 1.0, 0.0],
        };
        // Fourth term has scale 0 → contributes 0 regardless of input.
        assert_eq!(m.score(1.0, 2.0, 3.0, 999.0), 1.0 + 2.0 + 3.0);
    }

    #[test]
    fn score_matches_standardized_linear_form() {
        // Reproduce the Stellar model's score on a hand-computed point.
        let coel = 0.5;
        let ln_i = 10.0;
        let rt = 0.9;
        let mp = 0.85;
        let z = |x: f64, mean: f64, scale: f64| (x - mean) / scale;
        let expected = 0.9933168416485256 * z(coel, 0.14931086687377143, 0.2074610054197197)
            + 0.047052481253413006 * z(ln_i, 9.15749607304815, 2.260347450504608)
            + 0.027130393118192445 * z(rt, 0.9158212545538758, 0.07333900267211861)
            + 0.10184133676728513 * z(mp, 0.8904037620854307, 0.08220791724094854);
        assert_eq!(STELLAR_MODEL.score(coel, ln_i, rt, mp), expected);
    }
}
