#!/usr/bin/env python3
"""
Assign synthetic outcomes to benchmark scores based on ipTM thresholds.
Simulates what real wet-lab validation might look like.

Logic: higher ipTM → higher probability of binding
- ipTM < 0.10: mostly false (5% hit rate)
- 0.10–0.15: low signal (15% hit rate)
- 0.15–0.25: moderate (40% hit rate)
- 0.25+: high (85% hit rate)
"""

import pandas as pd
import random
from pathlib import Path

def main():
    random.seed(42)

    # Load scores
    scores_csv = Path("results/benchmark_scores_with_outcomes.csv")
    df = pd.read_csv(scores_csv)
    print(f"Loaded {len(df)} scored pairs")

    # Assign mock outcomes based on ipTM
    outcomes = []
    for iptm in df["ipTM"]:
        if pd.isna(iptm):
            outcomes.append("unknown")
        elif iptm < 0.10:
            outcomes.append("true" if random.random() < 0.05 else "false")
        elif iptm < 0.15:
            outcomes.append("true" if random.random() < 0.15 else "false")
        elif iptm < 0.25:
            outcomes.append("true" if random.random() < 0.40 else "false")
        else:
            outcomes.append("true" if random.random() < 0.85 else "false")

    df["outcome"] = outcomes
    # Ensure outcomes stay as strings, not booleans
    df.to_csv(scores_csv, index=False, quoting=1)  # quoting=1 forces all strings to be quoted

    # Summary
    true_count = (df["outcome"] == "true").sum()
    false_count = (df["outcome"] == "false").sum()
    print(f"[OK] Created mock outcomes:")
    print(f"     {true_count} true (binders)")
    print(f"     {false_count} false (non-binders)")
    print(f"     Hit rate: {true_count/(true_count+false_count):.1%}")

if __name__ == "__main__":
    main()
