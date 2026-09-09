#!/usr/bin/env python3
"""
Fold all benchmark sequence pairs and save ColabFold prediction outputs.

This script orchestrates structure prediction for a benchmark dataset of protein
complex pairs (target + designed binder). It reads sequence pairs from
benchmark/adaptyv_benchmark.csv, folds them using ColabFold (AF2-multimer),
and saves the output JSONs and PDB structures.

Input:
  benchmark/adaptyv_benchmark.csv — columns: sequence, target_seq, outcome, source_system

Output:
  results/benchmark_folds_raw/[pair_name]_scores_rank_001_*.json — ColabFold output
  results/benchmark_folds_raw/[pair_name]_unrelaxed_rank_001_*.pdb — structures

Usage:
  python 01_fold_benchmark.py [--resume] [--limit N]

Flags:
  --resume    Skip pairs that already have output files (checkpoint recovery)
  --limit N   Fold only first N pairs (for testing; default: all)

Typical runtime:
  • ~2–3 minutes per pair (single-sequence, 1 model, 3 recycles, RTX 3060)
  • 120–150 pairs ≈ 4–8 hours wall-clock (parallelization possible with batch submission)
  • Peak VRAM: ~3–4 GB per pair (with TF_FORCE_UNIFIED_MEMORY and XLA flags)

Environment:
  Requires ColabFold conda environment with JAX, jax-cuda12-pallas, and dependencies.
  Run: bash environment/setup.sh first to activate.

Checkpointing:
  After every 10 pairs, logs a recovery file. If interrupted, re-run with --resume
  to skip completed pairs and continue.

Caveats:
  • Single-sequence mode (no MSA lookup) is fast but lower-accuracy than real-MSA folds
  • Scores reported are from single-model ColabFold, suitable for ranking not final evaluation
  • Network calls to ColabFold's public MMseqs2 server (rate-limited, shared resource)
"""

import os
import json
import argparse
import subprocess
import sys
from pathlib import Path
import pandas as pd
from datetime import datetime

def main():
    parser = argparse.ArgumentParser(
        description="Fold benchmark sequence pairs using ColabFold"
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip pairs that already have output; checkpoint recovery"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Fold only first N pairs (for testing)"
    )
    args = parser.parse_args()

    # Paths
    benchmark_csv = Path("benchmark/adaptyv_benchmark.csv")
    output_dir = Path("results/benchmark_folds_raw")
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_file = output_dir / "checkpoint.txt"

    # Load benchmark
    if not benchmark_csv.exists():
        print(f"[FAIL] Benchmark file not found: {benchmark_csv}")
        sys.exit(1)

    df = pd.read_csv(benchmark_csv)
    print(f"Loaded {len(df)} benchmark pairs from {benchmark_csv}")

    # Limit if requested
    if args.limit:
        df = df.head(args.limit)
        print(f"Limited to first {len(df)} pairs")

    # Track completed
    completed = set()
    if args.resume and checkpoint_file.exists():
        with open(checkpoint_file) as f:
            completed = set(l.strip() for l in f if l.strip())
        print(f"Resume mode: skipping {len(completed)} already-completed pairs")

    # Fold each pair
    for pos, (idx, row) in enumerate(df.iterrows()):
        pair_name = f"pair_{pos:03d}"

        # Check if already done
        if args.resume:
            score_json = sorted(output_dir.glob(f"{pair_name}_scores_rank_001_*.json"))
            if score_json and pair_name in completed:
                print(f"[{pos+1}/{len(df)}] {pair_name} already done, skipping")
                continue

        # Build input FASTA
        target_seq = row.get("target_seq", "")
        binder_seq = row.get("sequence", "")
        fasta_content = f">{pair_name}\n{target_seq}:{binder_seq}\n"

        input_fasta = output_dir / f"{pair_name}.fasta"
        input_fasta.write_text(fasta_content)

        # Run ColabFold (low-VRAM, fast-screen settings)
        print(f"[{pos+1}/{len(df)}] Folding {pair_name}...")
        env = os.environ.copy()
        env.update({
            "TF_FORCE_UNIFIED_MEMORY": "1",
            "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
            "XLA_PYTHON_CLIENT_MEM_FRACTION": "4.0",
            "XLA_PYTHON_CLIENT_ALLOCATOR": "platform",
        })

        try:
            result = subprocess.run(
                [
                    "colabfold_batch",
                    "--msa-mode", "single_sequence",
                    "--model-type", "alphafold2_multimer_v3",
                    "--num-models", "1",
                    "--num-recycle", "3",
                    str(input_fasta),
                    str(output_dir),
                ],
                env=env,
                capture_output=True,
                text=True,
                timeout=600,
            )
            if result.returncode != 0:
                print(f"  [FAIL] Error: {result.stderr}")
                continue
            print(f"  [OK] Success")
        except subprocess.TimeoutExpired:
            print(f"  [FAIL] Timeout (>10min)")
            continue
        except FileNotFoundError:
            print(f"  [FAIL] colabfold_batch not found — run: bash environment/setup.sh")
            sys.exit(1)

        # Log completion every 10 pairs
        completed.add(pair_name)
        if (pos + 1) % 10 == 0:
            with open(checkpoint_file, "w") as f:
                f.write("\n".join(sorted(completed)))
            print(f"  Checkpoint: {len(completed)} pairs done")

    # Final checkpoint
    with open(checkpoint_file, "w") as f:
        f.write("\n".join(sorted(completed)))

    print(f"\n[OK] Complete: folded {len(completed)} pairs")
    print(f"Results in: {output_dir}")

if __name__ == "__main__":
    main()
