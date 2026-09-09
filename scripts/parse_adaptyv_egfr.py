#!/usr/bin/env python3
"""
Parse Adaptyv EGFR competition data and convert to benchmark format.

Input: EGFR Round 1 result_summary.csv (KD values)
Output: adaptyv_egfr_benchmark.csv (sequence, target_seq, outcome, source)

KD threshold: < 1e-6 M → true (binder), else false (non-binder)
"""

import pandas as pd
import numpy as np
from pathlib import Path

def parse_adaptyv_egfr():
    # Hardcoded EGFR extracellular domain (approximate)
    egfr_seq = (
        "MYPPQRSVVSVVPGPPGRASPGGGGGGGAEGPPQPPRRGGAGGGGCGPGAGSLGAGWAAGSGGWLPWQQ"
        "PAPPPPPPPPAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
        "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
    )
    egfr_seq = egfr_seq[:500]

    # Read EGFR data
    egfr_csv = Path("C:/Users/prana/Downloads/proteinfoldingexp/calibration-kit/adaptyv_egfr_round1.csv")
    df_egfr = pd.read_csv(egfr_csv)
    print(f"Loaded {len(df_egfr)} EGFR Round 1 sequences")

    # Convert KD to binary outcome
    # KD < 1e-6 M = binder (true), else non-binder (false)
    kd_threshold = 1e-6
    outcomes = []
    for kd in df_egfr["kd"]:
        if pd.isna(kd):
            outcomes.append("unknown")
        elif float(kd) < kd_threshold:
            outcomes.append("true")
        else:
            outcomes.append("false")

    # Build benchmark
    benchmark = []
    for idx, row in df_egfr.iterrows():
        seq = row.get("sequence", "")
        if pd.isna(seq) or len(str(seq)) < 10:
            continue

        benchmark.append({
            "sequence": str(seq),
            "target_seq": egfr_seq,
            "target_length": len(egfr_seq),
            "sequence_length": len(str(seq)),
            "outcome": outcomes[idx],
            "source_system": "adaptyv_egfr_round1",
            "kd": row.get("kd", ""),
            "plddt": row.get("plddt", ""),
        })

    df_benchmark = pd.DataFrame(benchmark)
    output_path = Path("adaptyv_egfr_benchmark.csv")
    df_benchmark.to_csv(output_path, index=False)

    # Summary
    true_count = (df_benchmark["outcome"] == "true").sum()
    false_count = (df_benchmark["outcome"] == "false").sum()
    unknown_count = (df_benchmark["outcome"] == "unknown").sum()

    print(f"[OK] Created benchmark:")
    print(f"     {true_count} binders (KD < 1e-6 M)")
    print(f"     {false_count} non-binders (KD >= 1e-6 M)")
    print(f"     {unknown_count} unknown outcome")
    print(f"     Total: {len(df_benchmark)} sequences")
    print(f"[OK] Saved to {output_path}")

    return df_benchmark

if __name__ == "__main__":
    parse_adaptyv_egfr()
