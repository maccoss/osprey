#!/usr/bin/env python3
"""
Inspect mokapot model weights to identify important features.

Usage:
    python inspect_mokapot_weights.py <mokapot.model.pkl>

This script loads a mokapot model file and displays the feature weights
from the trained linear model, ranked by absolute importance.
"""

import sys
import pickle
import argparse
from pathlib import Path


def load_mokapot_model(model_path: str):
    """Load a mokapot model from a pickle file."""
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    return model


def extract_feature_weights(model):
    """Extract feature names and weights from a mokapot model.

    Mokapot uses a linear SVM (Percolator-style) which stores weights
    in the model's estimator. The structure depends on mokapot version.
    """
    weights = {}

    # Handle list of models (one per cross-validation fold)
    if isinstance(model, list):
        print(f"Found {len(model)} cross-validation fold models")
        # Average weights across folds
        fold_weights = []
        for fold_model in model:
            fw = _get_weights_from_single_model(fold_model)
            if fw:
                fold_weights.append(fw)

        if fold_weights:
            # Get feature names from first fold
            feature_names = list(fold_weights[0].keys())
            for name in feature_names:
                values = [fw.get(name, 0) for fw in fold_weights]
                weights[name] = sum(values) / len(values)
    else:
        weights = _get_weights_from_single_model(model)

    return weights


def _get_weights_from_single_model(model):
    """Extract weights from a single mokapot model object."""
    weights = {}

    # Try different model structures based on mokapot version
    try:
        # Mokapot >= 0.10: Model wraps sklearn estimator
        if hasattr(model, 'estimator'):
            estimator = model.estimator
        elif hasattr(model, 'model'):
            estimator = model.model
        else:
            estimator = model

        # Get feature names
        feature_names = None
        if hasattr(model, 'features'):
            feature_names = model.features
        elif hasattr(model, 'feature_names'):
            feature_names = model.feature_names
        elif hasattr(model, '_feature_columns'):
            feature_names = model._feature_columns

        # Get coefficients from linear model
        if hasattr(estimator, 'coef_'):
            coef = estimator.coef_.flatten()
            if feature_names is None:
                feature_names = [f"feature_{i}" for i in range(len(coef))]
            for name, weight in zip(feature_names, coef):
                weights[name] = float(weight)

        # For SVC, coefficients might be in dual_coef_ or elsewhere
        elif hasattr(estimator, 'dual_coef_'):
            # This is more complex for kernel SVMs - weights aren't directly interpretable
            print("Warning: Model uses kernel SVM, weights are in dual form")

    except Exception as e:
        print(f"Warning: Could not extract weights: {e}")

    return weights


def print_weights_table(weights: dict, top_n: int = None):
    """Print feature weights as a formatted table, ranked by absolute importance."""
    if not weights:
        print("No weights found in model")
        return

    # Sort by absolute weight
    sorted_weights = sorted(weights.items(), key=lambda x: abs(x[1]), reverse=True)

    if top_n:
        sorted_weights = sorted_weights[:top_n]

    print("\n" + "=" * 70)
    print("FEATURE WEIGHTS (ranked by absolute importance)")
    print("=" * 70)
    print(f"{'Rank':<6} {'Feature':<40} {'Weight':>12} {'|Weight|':>10}")
    print("-" * 70)

    for rank, (name, weight) in enumerate(sorted_weights, 1):
        abs_weight = abs(weight)
        bar_len = min(int(abs_weight * 10), 20)  # Scale for visualization
        bar = "+" * bar_len if weight > 0 else "-" * bar_len
        print(f"{rank:<6} {name:<40} {weight:>12.4f} {abs_weight:>10.4f}  {bar}")

    print("=" * 70)

    # Summary statistics
    total_features = len(weights)
    nonzero_features = sum(1 for w in weights.values() if abs(w) > 1e-6)
    print(f"\nTotal features: {total_features}")
    print(f"Non-zero features: {nonzero_features}")

    # Identify potentially unimportant features (very small weights)
    threshold = 0.01
    small_weight_features = [name for name, w in weights.items() if abs(w) < threshold]
    if small_weight_features:
        print(f"\nFeatures with |weight| < {threshold} (potentially removable):")
        for name in small_weight_features:
            print(f"  - {name}: {weights[name]:.6f}")


def main():
    parser = argparse.ArgumentParser(
        description="Inspect mokapot model weights to identify important features"
    )
    parser.add_argument(
        "model_file",
        type=str,
        help="Path to mokapot model pickle file (mokapot.model.pkl)"
    )
    parser.add_argument(
        "--top",
        type=int,
        default=None,
        help="Show only top N features by importance"
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Optional: Save weights to TSV file"
    )

    args = parser.parse_args()

    model_path = Path(args.model_file)
    if not model_path.exists():
        print(f"Error: Model file not found: {model_path}")
        sys.exit(1)

    print(f"Loading model from: {model_path}")
    model = load_mokapot_model(model_path)

    # Debug: Show model structure
    print(f"Model type: {type(model)}")
    if hasattr(model, '__dict__'):
        print(f"Model attributes: {list(model.__dict__.keys())}")

    weights = extract_feature_weights(model)
    print_weights_table(weights, top_n=args.top)

    # Optionally save to TSV
    if args.output:
        output_path = Path(args.output)
        sorted_weights = sorted(weights.items(), key=lambda x: abs(x[1]), reverse=True)
        with open(output_path, 'w') as f:
            f.write("feature\tweight\tabs_weight\n")
            for name, weight in sorted_weights:
                f.write(f"{name}\t{weight:.6f}\t{abs(weight):.6f}\n")
        print(f"\nWeights saved to: {output_path}")


if __name__ == "__main__":
    main()
