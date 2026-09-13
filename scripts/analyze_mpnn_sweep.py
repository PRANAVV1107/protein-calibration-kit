#!/usr/bin/env python3
"""Compare sequence composition across ProteinMPNN sampling configurations.

Reads the sweep output produced by sweep_mpnn_settings.sh and reports, per
configuration, how much low-complexity collapse remains -- using the same
composition penalty the RL reward applies.
"""
import glob
import os
import sys
from collections import Counter

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_composition import composition_penalty  # noqa: E402

SWEEP = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit/results/mpnn_sweep"

CONFIGS = [
    ("baseline_t0.1", "temp 0.1, no bias  (CURRENT)"),
    ("t0.2", "temp 0.2, no bias"),
    ("t0.3", "temp 0.3, no bias"),
    ("t0.1_biasA1", "temp 0.1, A bias -1.0"),
    ("t0.2_biasA1", "temp 0.2, A bias -1.0"),
    ("t0.2_biasA2", "temp 0.2, A bias -2.0"),
    ("t0.2_biasA1G", "temp 0.2, A -1.0 / G -0.5"),
]


def read_designed_seqs(cfg_dir):
    """ProteinMPNN .fa layout: line0 header, line1 original, line2 header, line3 designed."""
    seqs = []
    for fa in sorted(glob.glob(os.path.join(cfg_dir, "seqs", "*.fa"))):
        with open(fa) as f:
            lines = f.readlines()
        if len(lines) >= 4:
            seqs.append(lines[3].strip())
    return seqs


def main():
    rows = []
    for cfg, label in CONFIGS:
        cfg_dir = os.path.join(SWEEP, cfg)
        seqs = read_designed_seqs(cfg_dir)
        if not seqs:
            print(f"[WARN] no sequences found for {cfg}")
            continue

        pens, runs, fracs, ala = [], [], [], []
        for s in seqs:
            p, d = composition_penalty(s, detail=True)
            pens.append(p)
            runs.append(d["max_run"])
            fracs.append(d["max_aa_fraction"])
            ala.append(s.count("A") / len(s))

        rows.append({
            "config": label,
            "n": len(seqs),
            "mean_penalty": sum(pens) / len(pens),
            "max_penalty": max(pens),
            "n_penalised>0.1": sum(1 for p in pens if p > 0.1),
            "mean_Ala%": 100 * sum(ala) / len(ala),
            "max_Ala%": 100 * max(ala),
            "mean_maxrun": sum(runs) / len(runs),
            "max_maxrun": max(runs),
        })

    df = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print("\nProteinMPNN sampling sweep — 50 IL-7Ra backbones\n")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\nLower mean_penalty / max_Ala% / max_maxrun is better,")
    print("but very high temperature can also degrade design quality — check ipTM separately.")


if __name__ == "__main__":
    main()
