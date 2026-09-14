#!/usr/bin/env python3
"""
Merge EGFR Round 1 & 2 Adaptyv data into benchmark format.
"""

import os
import pandas as pd

KITDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
from pathlib import Path

# EGFR domain III -- the cetuximab epitope, which is what the Adaptyv
# competition targets. From PDB 1YY9 chain A, residues 310-510, no crystal gaps,
# verified against that file's own TITLE/COMPND records.
#
# WARNING: versions of this script before 2026-09-13 carried a FABRICATED
# placeholder here -- poly-alanine/poly-proline filler, not a real protein.
# Any benchmark row or score produced by those versions is invalid.
egfr_seq = (
    "RKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITG"
    "FLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYAN"
    "TINWKKLFGTSGQKTKIISNRGENKCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRE"
)

all_data = []

# Process Round 2 (has explicit "binding" column)
print("Processing EGFR Round 2...")
df_r2 = pd.read_csv(os.path.join(KITDIR, "adaptyv_egfr_round2.csv"))
print(f"  Loaded {len(df_r2)} sequences")

for idx, row in df_r2.iterrows():
    seq = row.get("sequence", "")
    if pd.isna(seq) or len(str(seq)) < 10:
        continue

    # The `binding` column holds the STRINGS 'true'/'false'/'unknown', not
    # booleans. `if row["binding"]` is truthy for every non-empty string, so an
    # earlier version of this line labelled all 325 verified NON-binders as
    # binders. Compare the text explicitly.
    b = str(row.get("binding", "")).strip().lower()
    if b in ("true", "false"):
        outcome = b
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
df_r1 = pd.read_csv(os.path.join(KITDIR, "adaptyv_egfr_round1.csv"))
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
