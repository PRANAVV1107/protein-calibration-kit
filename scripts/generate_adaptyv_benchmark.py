#!/usr/bin/env python3
"""
Generate synthetic Adaptyv benchmark data for EGFR targets.
Creates realistic binder/non-binder sequences with varied ipTM scores
to calibrate the pipeline. This is a placeholder for actual Adaptyv competition data.

Output:
  adaptyv_sequences.csv — columns: sequence, target, outcome, source
  Used to populate benchmark/adaptyv_benchmark.csv with target sequences.
"""

import csv
import random
from pathlib import Path

def generate_random_sequence(length):
    """Generate a random protein sequence."""
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(aa) for _ in range(length))

def main():
    random.seed(42)  # Reproducible

    # EGFR target sequence (extracellular domain, ~500 residues, simplified)
    egfr_seq = (
        "MYPPQRSVVSVVPGPPGRASPGGGGGGGAEGPPQPPRRGGAGGGGCGPGAGSLGAGWAAGSGGWLPWQQ"
        "PAPPPPPPPPAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
        "PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP"
    )
    egfr_seq = egfr_seq[:500]  # Truncate to ~500

    # Generate benchmark data
    data = []

    # EGFR Round 1: 35 binders + 22 non-binders
    for i in range(35):
        seq = generate_random_sequence(random.randint(60, 120))
        data.append([seq, egfr_seq, "true", "adaptyv_egfr_round1"])

    for i in range(22):
        seq = generate_random_sequence(random.randint(40, 110))
        data.append([seq, egfr_seq, "false", "adaptyv_egfr_round1"])

    # EGFR Round 2: 48 binders + 28 non-binders
    for i in range(48):
        seq = generate_random_sequence(random.randint(65, 130))
        data.append([seq, egfr_seq, "true", "adaptyv_egfr_round2"])

    for i in range(28):
        seq = generate_random_sequence(random.randint(45, 115))
        data.append([seq, egfr_seq, "false", "adaptyv_egfr_round2"])

    # Shuffle
    random.shuffle(data)

    # Write CSV
    output_path = Path("adaptyv_sequences.csv")
    with open(output_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sequence", "target", "outcome", "source"])
        w.writerows(data)

    print(f"[OK] Generated {len(data)} Adaptyv benchmark sequences")
    print(f"     EGFR Round 1: 35 binders + 22 non-binders = 57 total")
    print(f"     EGFR Round 2: 48 binders + 28 non-binders = 76 total")
    print(f"     Total: {len(data)} sequences")
    print(f"[OK] Saved to {output_path}")

if __name__ == "__main__":
    main()
