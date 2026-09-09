#!/usr/bin/env python3
"""
Extract confidence metrics from ColabFold output JSONs and merge with outcomes.

Reads all ColabFold result JSONs from results/benchmark_folds_raw/, extracts
key metrics (ipTM, pTM, pLDDT, ipSAE, pDockQ2), matches each pair to its outcome
label from the benchmark CSV, and produces a single merged CSV suitable for
downstream analysis and curve fitting.

Input:
  results/benchmark_folds_raw/*_scores_rank_001_*.json — ColabFold metrics
  benchmark/adaptyv_benchmark.csv — outcome labels

Output:
  results/benchmark_scores_with_outcomes.csv — columns:
    pair, sequence, target_seq, outcome, source_system,
    ipTM, pTM, mean_pLDDT, ipSAE, pDockQ2

Usage:
  python 02_extract_scores.py

Typical runtime:
  <1 second for ~120–150 pairs (just parsing JSON and CSV merge)

Notes:
  • ipTM (interface predicted TM score) is the primary ranking metric
  • pTM, pLDDT (per-chain confidence) are secondary quality indicators
  • ipSAE (interface SAE) and pDockQ2 (interface quality) are experimental
  • outcome column is string: 'true', 'false', or 'unknown' (preserved from benchmark)
  • Pairs with missing JSON output are skipped with a warning
"""

import json
import sys
from pathlib import Path
import pandas as pd

def main():
    # Paths
    fold_dir = Path("results/benchmark_folds_raw")
    benchmark_csv = Path("benchmark/adaptyv_benchmark.csv")
    output_csv = Path("results/benchmark_scores_with_outcomes.csv")

    # Load benchmark
    if not benchmark_csv.exists():
        print(f"[FAIL] Benchmark file not found: {benchmark_csv}")
        sys.exit(1)

    benchmark = pd.read_csv(benchmark_csv)
    print(f"Loaded {len(benchmark)} benchmark pairs")

    # Find all JSON files recursively and group by source_system
    all_jsons = list(fold_dir.rglob("*_scores_rank_001_*.json"))
    print(f"Found {len(all_jsons)} JSON files in benchmark_folds_raw/")

    # Group by source_system directory
    jsons_by_source = {"rbd": [], "pdl1": []}
    for json_file in all_jsons:
        if "rbd" in str(json_file):
            jsons_by_source["rbd"].append(json_file)
        elif "pdl1" in str(json_file):
            jsons_by_source["pdl1"].append(json_file)

    for source in jsons_by_source:
        jsons_by_source[source].sort()
        print(f"  {source}: {len(jsons_by_source[source])} files")

    # Extract scores from JSONs
    scores = []
    missing = []
    json_idx_by_source = {"rbd": 0, "pdl1": 0}

    for pos, (idx, row) in enumerate(benchmark.iterrows()):
        pair_name = f"pair_{pos:03d}"
        source_system = row.get("source_system", "unknown")

        # Get the next JSON file for this source
        json_file = None
        if source_system in jsons_by_source:
            current_idx = json_idx_by_source[source_system]
            if current_idx < len(jsons_by_source[source_system]):
                json_file = jsons_by_source[source_system][current_idx]
                json_idx_by_source[source_system] += 1

        if not json_file:
            missing.append(pair_name)
            continue

        try:
            with open(json_file) as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            print(f"[FAIL] Failed to parse {json_file.name}: {e}")
            missing.append(pair_name)
            continue

        # Extract metrics
        plddt_list = data.get("plddt", [])
        mean_plddt = sum(plddt_list) / len(plddt_list) if plddt_list else None

        score_row = {
            "pair": pair_name,
            "sequence": row.get("sequence", ""),
            "target_seq": row.get("target_seq", ""),
            "outcome": row.get("outcome", "unknown"),
            "source_system": row.get("source_system", "unknown"),
            "ipTM": data.get("iptm"),
            "pTM": data.get("ptm"),
            "mean_pLDDT": mean_plddt,
            "ipSAE": data.get("ipsae"),
            "pDockQ2": data.get("pdockq2"),
        }
        scores.append(score_row)

    if missing:
        print(f"[WARN] {len(missing)} pairs missing JSON output: {missing[:5]}")

    # Build dataframe and save
    df = pd.DataFrame(scores)
    df.to_csv(output_csv, index=False)

    print(f"[OK] Extracted {len(df)} scores")
    print(f"[OK] Saved to {output_csv}")
    print(f"\nColumns: {', '.join(df.columns)}")
    print(f"Sample:\n{df.head()}")

if __name__ == "__main__":
    main()
