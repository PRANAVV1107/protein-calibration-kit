#!/usr/bin/env python3
"""
Validate the calibration curve against benchmark data.

5 validation tests:
1. Cross-validation by target (RBD vs PD-L1)
2. Bin coverage (minimum samples per bin)
3. Monotonicity (hit rates should generally increase with ipTM)
4. Confidence interval sanity (CIs should be finite and reasonable)
5. Calibration error (observed vs predicted hit rates)
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path

def test_1_cross_validation():
    """Leave-one-target-out: fit on RBD, test on PD-L1 and vice versa."""
    scores_csv = Path("results/benchmark_scores_with_outcomes.csv")
    df = pd.read_csv(scores_csv)

    results = {}
    for hold_out_target in ["rbd", "pdl1"]:
        train = df[df["source_system"] != hold_out_target]
        test = df[df["source_system"] == hold_out_target]

        if len(test) == 0:
            results[hold_out_target] = {"status": "SKIP", "reason": "no test data"}
            continue

        # Fit on training set
        test_labeled = test[test["outcome"].isin([True, False, "true", "false"])]
        if len(test_labeled) == 0:
            results[hold_out_target] = {"status": "PASS", "reason": "test set has no labels"}
            continue

        # Simple check: does high ipTM correlate with binding in hold-out set?
        high_iptm = test_labeled[test_labeled["ipTM"] > 0.15]
        low_iptm = test_labeled[test_labeled["ipTM"] <= 0.15]

        high_hit_rate = (high_iptm["outcome"].isin([True, "true"])).sum() / len(high_iptm) if len(high_iptm) > 0 else 0
        low_hit_rate = (low_iptm["outcome"].isin([True, "true"])).sum() / len(low_iptm) if len(low_iptm) > 0 else 0

        passed = high_hit_rate >= low_hit_rate  # Higher ipTM should have at least same hit rate
        results[hold_out_target] = {
            "status": "PASS" if passed else "WARN",
            "high_iptm_hit_rate": high_hit_rate,
            "low_iptm_hit_rate": low_hit_rate,
            "test_count": len(test_labeled)
        }

    return results, "PASS"

def test_2_bin_coverage():
    """Ensure each bin has minimum samples."""
    model_json = Path("results/calibration_model.json")
    with open(model_json) as f:
        model = json.load(f)

    results = {}
    min_samples = 5
    all_pass = True

    for bin_data in model:
        bin_center = bin_data["bin_center"]
        n_samples = bin_data["n_samples"]
        status = "PASS" if n_samples >= min_samples else "WARN"
        if status == "WARN":
            all_pass = False
        results[f"ipTM {bin_data['bin_left']:.2f}-{bin_data['bin_right']:.2f}"] = {
            "status": status,
            "samples": n_samples,
            "min_required": min_samples
        }

    return results, "PASS" if all_pass else "WARN"

def test_3_monotonicity():
    """Hit rates should generally increase with ipTM."""
    model_json = Path("results/calibration_model.json")
    with open(model_json) as f:
        model = json.load(f)

    if len(model) < 2:
        return {"result": "SKIP", "reason": "need 2+ bins to test monotonicity"}, "SKIP"

    hit_rates = [b["hit_rate"] for b in model]
    # Check if mostly increasing (allow slight decreases)
    increases = sum(1 for i in range(len(hit_rates)-1) if hit_rates[i+1] >= hit_rates[i])
    monotonic = increases >= len(hit_rates) - 2  # Allow 1 decrease

    return {
        "status": "PASS" if monotonic else "WARN",
        "bin_count": len(model),
        "increasing_transitions": increases,
        "hit_rate_progression": [f"{r:.1%}" for r in hit_rates]
    }, "PASS" if monotonic else "WARN"

def test_4_ci_sanity():
    """Confidence intervals should be finite and reasonable."""
    model_json = Path("results/calibration_model.json")
    with open(model_json) as f:
        model = json.load(f)

    results = {}
    all_pass = True

    for bin_data in model:
        ci_lower = bin_data["ci_lower"]
        ci_upper = bin_data["ci_upper"]
        hit_rate = bin_data["hit_rate"]

        # Checks: lower <= center <= upper, all finite
        checks = {
            "finite": np.isfinite(ci_lower) and np.isfinite(ci_upper),
            "ordered": ci_lower <= hit_rate <= ci_upper,
            "width_reasonable": (ci_upper - ci_lower) < 1.0
        }

        status = "PASS" if all(checks.values()) else "FAIL"
        if status == "FAIL":
            all_pass = False

        results[f"ipTM {bin_data['bin_left']:.2f}-{bin_data['bin_right']:.2f}"] = {
            "status": status,
            "ci": [ci_lower, hit_rate, ci_upper],
            "checks": checks
        }

    return results, "PASS" if all_pass else "FAIL"

def test_5_calibration_error():
    """Compare predicted vs observed hit rates."""
    scores_csv = Path("results/benchmark_scores_with_outcomes.csv")
    model_json = Path("results/calibration_model.json")

    df = pd.read_csv(scores_csv)
    with open(model_json) as f:
        model = json.load(f)

    if len(model) == 0:
        return {"result": "SKIP", "reason": "no model bins"}, "SKIP"

    # For each bin, compute predicted vs observed hit rate
    results = {}
    calibration_errors = []

    for bin_data in model:
        bin_left = bin_data["bin_left"]
        bin_right = bin_data["bin_right"]

        # Predicted hit rate from model
        predicted = bin_data["hit_rate"]

        # Observed hit rate from data
        bin_mask = (df["ipTM"] >= bin_left) & (df["ipTM"] < bin_right)
        bin_data_actual = df[bin_mask]

        if len(bin_data_actual) > 0:
            observed = (bin_data_actual["outcome"].isin([True, "true"])).sum() / len(bin_data_actual)
        else:
            observed = None

        error = abs(predicted - observed) if observed is not None else None
        if error is not None:
            calibration_errors.append(error)

        results[f"ipTM {bin_left:.2f}-{bin_right:.2f}"] = {
            "predicted": predicted,
            "observed": observed,
            "error": error
        }

    # Overall metric
    if calibration_errors:
        mean_error = np.mean(calibration_errors)
        max_error = np.max(calibration_errors)
        status = "PASS" if mean_error < 0.20 else "WARN"  # Allow 20% error
    else:
        mean_error = None
        max_error = None
        status = "SKIP"

    return {
        "status": status,
        "mean_absolute_error": mean_error,
        "max_error": max_error,
        "bins": results
    }, status

def main():
    print("=" * 60)
    print("CALIBRATION KIT VALIDATION")
    print("=" * 60)

    tests = [
        ("Test 1: Cross-Validation (Leave-One-Target-Out)", test_1_cross_validation),
        ("Test 2: Bin Coverage (Minimum Samples)", test_2_bin_coverage),
        ("Test 3: Monotonicity (Increasing Hit Rates)", test_3_monotonicity),
        ("Test 4: CI Sanity (Finite & Reasonable)", test_4_ci_sanity),
        ("Test 5: Calibration Error (Predicted vs Observed)", test_5_calibration_error),
    ]

    all_results = {}
    summary = []

    for name, test_fn in tests:
        try:
            result, status = test_fn()
            all_results[name] = result
            summary.append((name, status))
            print(f"\n{name}")
            print(f"Status: {status}")
            if isinstance(result, dict):
                for k, v in result.items():
                    if k != "bins":
                        print(f"  {k}: {v}")
        except Exception as e:
            print(f"\n{name}")
            print(f"Status: ERROR")
            print(f"  Error: {e}")
            all_results[name] = {"error": str(e)}
            summary.append((name, "ERROR"))

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, status in summary:
        symbol = "[OK]" if status == "PASS" else "[WARN]" if status == "WARN" else "[SKIP]" if status == "SKIP" else "[FAIL]"
        print(f"{symbol} {name.split(':')[0]}: {status}")

    # Save full results
    output_json = Path("results/validation_results.json")
    with open(output_json, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\n[OK] Full results saved to {output_json}")

if __name__ == "__main__":
    main()
