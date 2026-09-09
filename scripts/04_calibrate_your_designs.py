#!/usr/bin/env python3
"""
Apply the calibration curve to user-designed sequences.

Reads a CSV of your de novo designs (with ipTM and other metrics), matches each
to the appropriate bin in the calibration model, and returns an estimated hit rate
(along with 95% confidence interval). This gives you a probability-of-binding
estimate *relative to the benchmark*, not an absolute prediction.

Input:
  examples/example_user_designs.csv (or your own file) — columns:
    sequence, target, ipTM, pTM, mean_pLDDT, [optional: ipSAE, pDockQ2]
  results/calibration_model.json — the fitted curve from 03_fit_curve.py

Output:
  results/[your_file]_calibrated.csv — columns:
    [all input columns] + estimated_hit_rate, ci_lower, ci_upper, notes

Usage:
  python 04_calibrate_your_designs.py examples/example_user_designs.csv

Interpretation:
  estimated_hit_rate is the empirical hit rate from the benchmark bin your
  design falls into. A design with ipTM=0.8 gets the hit rate of the
  0.75–0.80 (or 0.80–0.85) bin from the benchmark curve.

  Confidence intervals show the range of plausible rates given the bin's sample size.
  Wider CIs = fewer benchmark samples in that region = less confident in the estimate.

  Notes column flags:
  • "target not in benchmark" — target differs from all benchmark targets; apply with caution
  • "ipTM out of range" — design's ipTM is above the highest binned score (extrapolation)
  • "no bin coverage" — too few samples in this bin; very wide CI

Important caveats:
  1. This is NOT a binding predictor; it's a *recalibration* of your pipeline's scores
     against empirical data.
  2. Assumes your designs have been scored with the same ColabFold settings (fast-screen
     or high-accuracy) as the benchmark.
  3. Target generalization: the model is fit on 4 targets (~215 binders); applying it
     to a novel target assumes similar scoring behavior (may not hold).
  4. Start with designs in high-confidence regions (e.g., ipTM > 0.2); extrapolation
     into low-ipTM regions is unreliable.
  5. Use this as a *prior*, not a point estimate. Update it with your own labeled data
     as you accumulate wet-lab results.
"""

import json
import sys
from pathlib import Path
import pandas as pd

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Calibrate your designs using the benchmark curve")
    parser.add_argument("input_csv", help="CSV of your designs (must have ipTM column)")
    args = parser.parse_args()

    input_path = Path(args.input_csv)
    if not input_path.exists():
        print(f"[FAIL] Input file not found: {input_path}")
        sys.exit(1)

    # Load user designs
    df_user = pd.read_csv(input_path)
    print(f"Loaded {len(df_user)} designs from {input_path}")

    # Load calibration model
    model_json = Path("results/calibration_model.json")
    if not model_json.exists():
        print(f"[FAIL] Calibration model not found: {model_json}")
        print("  Run: python 03_fit_curve.py first")
        sys.exit(1)

    with open(model_json) as f:
        model = json.load(f)

    # Check that model has bins
    if not model:
        print(f"[FAIL] Calibration model is empty (no bins fit from benchmark data)")
        print("       This may happen if the benchmark has no confirmed outcomes.")
        print("       Check: does benchmark/adaptyv_benchmark.csv have outcome='true' or 'false'?")
        sys.exit(1)

    # Apply calibration
    calibrated = []
    for idx, row in df_user.iterrows():
        iptm = row.get("ipTM")
        if pd.isna(iptm):
            print(f"[FAIL] Row {idx} missing ipTM score, skipping")
            continue

        # Find bin
        bin_match = None
        for bin_data in model:
            if bin_data["bin_left"] <= iptm < bin_data["bin_right"]:
                bin_match = bin_data
                break

        # Out of range?
        if bin_match is None:
            # Find closest bin
            if iptm > max(b["bin_right"] for b in model):
                bin_match = max(model, key=lambda b: b["bin_right"])
                note = "ipTM out of range (extrapolation)"
            elif iptm < min(b["bin_left"] for b in model):
                bin_match = min(model, key=lambda b: b["bin_left"])
                note = "ipTM out of range (extrapolation)"
            else:
                # Should not happen
                print(f"[FAIL] Row {idx}: no bin match for ipTM={iptm}, skipping")
                continue
        else:
            note = ""

        # Add sample-size warning
        if bin_match["n_samples"] < 5:
            note = note or "low bin coverage"

        result_row = dict(row)
        result_row["estimated_hit_rate"] = bin_match["hit_rate"]
        result_row["ci_lower"] = bin_match["ci_lower"]
        result_row["ci_upper"] = bin_match["ci_upper"]
        result_row["n_benchmark_samples_in_bin"] = bin_match["n_samples"]
        result_row["notes"] = note
        calibrated.append(result_row)

    # Save
    df_out = pd.DataFrame(calibrated)
    output_path = input_path.parent / f"{input_path.stem}_calibrated.csv"
    df_out.to_csv(output_path, index=False)

    print(f"[OK] Calibrated {len(df_out)} designs")
    print(f"[OK] Saved to {output_path}")
    print(f"\nSample (first 3 rows):")
    print(df_out[["sequence", "ipTM", "estimated_hit_rate", "ci_lower", "ci_upper", "notes"]].head(3).to_string())

if __name__ == "__main__":
    main()
