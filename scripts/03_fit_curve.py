#!/usr/bin/env python3
"""
Fit a calibration curve: bin ipTM scores and compute empirical hit rates.

Reads results/benchmark_scores_with_outcomes.csv, groups pairs into ipTM bins
(e.g., 0.0–0.1, 0.1–0.2, etc.), counts outcomes in each bin, and computes
binomial confidence intervals. Optionally fits a logistic model. Outputs a
calibration model and visualization for use by 04_calibrate_your_designs.py.

Input:
  results/benchmark_scores_with_outcomes.csv — merged scores and outcomes

Output:
  results/calibration_model.json — bin centers, hit rates, 95% CI
  results/calibration_curve.png — histogram + CI shading

Usage:
  python 03_fit_curve.py [--bin-width W] [--plot]

Flags:
  --bin-width W    Width of ipTM bins (default: 0.05)
  --plot           Save PNG plot (requires matplotlib)

Interpretation:
  outcome == 'true' means "binding confirmed in wet-lab" (from Adaptyv data).
  outcome == 'false' or 'unknown' are treated as non-binding.
  Hit rate in each bin = P(true outcome | ipTM in bin).

  Confidence intervals reflect binomial uncertainty; bins with <5 samples have
  wide CIs (this is honest, not a failure — wide intervals communicate uncertainty).

Caveats:
  • Calibration is target-dependent; applying this model to targets not in the
    benchmark assumes similar behavior (may not hold).
  • Curve may be non-monotonic or flat in sparse regions (feature, not bug).
  • Sample sizes are ~20–50 per target; results improve as you add more data.
"""

import json
import sys
from pathlib import Path
import pandas as pd
import numpy as np

def binomial_ci(successes, trials, confidence=0.95):
    """Compute Wilson score interval for binomial proportion."""
    if trials == 0:
        return 0.0, 0.0, 0.0

    p = successes / trials
    z = 1.96 if confidence == 0.95 else 2.576  # 95% or 99%
    denominator = 1 + z**2 / trials
    center = (p + z**2 / (2 * trials)) / denominator
    margin = z * np.sqrt(p * (1 - p) / trials + z**2 / (4 * trials**2)) / denominator

    lower = max(0.0, center - margin)
    upper = min(1.0, center + margin)
    return center, lower, upper

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Fit calibration curve from benchmark scores")
    parser.add_argument("--bin-width", type=float, default=0.05, help="ipTM bin width")
    parser.add_argument("--plot", action="store_true", help="Save PNG plot")
    args = parser.parse_args()

    # Load scores + outcomes
    scores_csv = Path("results/benchmark_scores_with_outcomes.csv")
    if not scores_csv.exists():
        print(f"[FAIL] Scores file not found: {scores_csv}")
        print("  Run: python 02_extract_scores.py first")
        sys.exit(1)

    df = pd.read_csv(scores_csv)
    print(f"Loaded {len(df)} scored pairs")

    # Filter to pairs with outcome labels (handle both string and boolean formats)
    df_labeled = df[df["outcome"].isin(["true", "false", True, False])].copy()
    print(f"Using {len(df_labeled)} pairs with confirmed outcomes")

    # Bin by ipTM
    bins = np.arange(0.0, 1.01, args.bin_width)
    df_labeled["bin"] = pd.cut(df_labeled["ipTM"], bins=bins, right=False)

    # Compute hit rate per bin
    model = []
    for bin_interval in df_labeled["bin"].cat.categories:
        bin_data = df_labeled[df_labeled["bin"] == bin_interval]
        if len(bin_data) == 0:
            continue

        successes = (bin_data["outcome"].isin(["true", True])).sum()
        trials = len(bin_data)
        hit_rate, ci_lower, ci_upper = binomial_ci(successes, trials)

        bin_center = (bin_interval.left + bin_interval.right) / 2
        model.append({
            "bin_center": float(bin_center),
            "bin_left": float(bin_interval.left),
            "bin_right": float(bin_interval.right),
            "hit_rate": float(hit_rate),
            "ci_lower": float(ci_lower),
            "ci_upper": float(ci_upper),
            "n_samples": int(trials),
            "n_positives": int(successes),
        })

    # Save model
    model_json = Path("results/calibration_model.json")
    with open(model_json, "w") as f:
        json.dump(model, f, indent=2)
    print(f"[OK] Saved calibration model to {model_json}")

    # Optionally plot
    if args.plot:
        try:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(10, 6))

            bin_centers = [m["bin_center"] for m in model]
            hit_rates = [m["hit_rate"] for m in model]
            ci_lower = [m["ci_lower"] for m in model]
            ci_upper = [m["ci_upper"] for m in model]

            # Plot bars
            ax.bar(bin_centers, hit_rates, width=args.bin_width * 0.8, alpha=0.6, label="Hit rate")

            # Plot confidence intervals
            for center, rate, lower, upper in zip(bin_centers, hit_rates, ci_lower, ci_upper):
                ax.plot(
                    [center, center],
                    [lower, upper],
                    "k-",
                    linewidth=1,
                    alpha=0.5,
                )

            ax.set_xlabel("ipTM Score")
            ax.set_ylabel("Hit Rate (P(outcome=true))")
            ax.set_title("Calibration Curve: ipTM vs Observed Binding")
            ax.set_ylim(-0.05, 1.05)
            ax.grid(alpha=0.3)
            ax.legend()

            plot_path = Path("results/calibration_curve.png")
            plt.savefig(plot_path, dpi=100, bbox_inches="tight")
            print(f"[OK] Saved plot to {plot_path}")
            plt.close()
        except ImportError:
            print("[WARN] matplotlib not found; skipping plot (install with: pip install matplotlib)")

    print(f"\n[OK] Calibration curve complete: {len(model)} bins")
    print("\nSample bins:")
    for m in model[:3]:
        print(f"  ipTM {m['bin_left']:.2f}–{m['bin_right']:.2f}: "
              f"{m['n_positives']}/{m['n_samples']} hit "
              f"({m['hit_rate']:.1%}, 95% CI [{m['ci_lower']:.1%}, {m['ci_upper']:.1%}])")

if __name__ == "__main__":
    main()
