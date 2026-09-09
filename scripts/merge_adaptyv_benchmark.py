#!/usr/bin/env python3
"""
Merge EGFR Round 1 & 2 Adaptyv data into benchmark format.
"""

import pandas as pd
from pathlib import Path

# EGFR extracellular domain (~500 residues)
egfr_seq = (
    "MYPPQRSVVSVVPGPPGRASPGGGGGGGAEGPPQPPRRGGAGGGGCGPGAGSLGAGWAAGSGGWLPWQQ"
    "PAPPPPPPPPAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
    "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
)
egfr_seq = egfr_seq[:500]

all_data = []

# Process Round 2 (has explicit "binding" column)
print("Processing EGFR Round 2...")
df_r2 = pd.read_csv("C:/Users/prana/Downloads/proteinfoldingexp/calibration-kit/adaptyv_egfr_round2.csv")
print(f"  Loaded {len(df_r2)} sequences")

for idx, row in df_r2.iterrows():
    seq = row.get("sequence", "")
    if pd.isna(seq) or len(str(seq)) < 10:
        continue

    # Use the "binding" column if available, else fall back to KD
    if "binding" in row and pd.notna(row["binding"]):
        outcome = "true" if row["binding"] else "false"
    else:
        kd = row.get("kd", None)
        if pd.isna(kd):
            outcome = "unknown"
        else:
            outcome = "true" if float(kd) < 1e-6 else "false"

    all_data.append({
        "sequence": str(seq),
        "target_seq": egfr_seq,
        "target_length": len(egfr_seq),
        "sequence_length": len(str(seq)),
        "outcome": outcome,
        "source_system": "adaptyv_egfr_round2",
        "kd": row.get("kd", ""),
        "binding_strength": row.get("binding_strength", ""),
    })

r2_count = len(all_data)
print(f"  Added {r2_count} sequences from Round 2")

# Process Round 1 (sparse data but include it)
print("Processing EGFR Round 1...")
df_r1 = pd.read_csv("C:/Users/prana/Downloads/proteinfoldingexp/calibration-kit/adaptyv_egfr_round1.csv")
print(f"  Loaded {len(df_r1)} sequences")

for idx, row in df_r1.iterrows():
    seq = row.get("sequence", "")
    if pd.isna(seq) or len(str(seq)) < 10:
        continue

    kd = row.get("kd", None)
    if pd.isna(kd):
        outcome = "unknown"
    else:
        outcome = "true" if float(kd) < 1e-6 else "false"

    all_data.append({
        "sequence": str(seq),
        "target_seq": egfr_seq,
        "target_length": len(egfr_seq),
        "sequence_length": len(str(seq)),
        "outcome": outcome,
        "source_system": "adaptyv_egfr_round1",
        "kd": row.get("kd", ""),
        "binding_strength": "",
    })

r1_added = len(all_data) - r2_count
print(f"  Added {r1_added} sequences from Round 1")

# Save merged benchmark
df = pd.DataFrame(all_data)
output_path = Path("benchmark/adaptyv_benchmark_new.csv")
df.to_csv(output_path, index=False)

# Summary
true_count = (df["outcome"] == "true").sum()
false_count = (df["outcome"] == "false").sum()
unknown_count = (df["outcome"] == "unknown").sum()

print(f"\n[OK] Merged Adaptyv benchmark:")
print(f"     {true_count} binders (true)")
print(f"     {false_count} non-binders (false)")
print(f"     {unknown_count} unknown")
print(f"     Total: {len(df)} sequences")
print(f"[OK] Saved to {output_path}")
print(f"\nNext: Replace benchmark/adaptyv_benchmark.csv with this file")
