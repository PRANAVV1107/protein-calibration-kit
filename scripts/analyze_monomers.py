#!/usr/bin/env python3
"""Does a designed binder fold on its own, and does that predict the complex score?

A binder that only folds when bound to its target is likely disordered in
solution -- poor expression, aggregation, or simply not the molecule the complex
prediction implied. Complex predictions cannot reveal this, because the target
scaffolds the binder in every one of them.

Each monomer here was folded alone (single sequence, no target). Every one is
paired with the complex ipTM of the same sequence, so the two can be correlated.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
sys.path.insert(0, os.path.join(KIT, "scripts"))
from seq_composition import composition_penalty  # noqa: E402


def main():
    os.chdir(KIT)
    man = pd.read_csv("results/queue/monomers_manifest.csv").set_index("monomer_name")

    rows = []
    for f in glob.glob("results/queue/monomers/*_scores_rank_001_*.json"):
        name = os.path.basename(f).split("_scores_rank_001_")[0]
        if name not in man.index:
            continue
        d = json.load(open(f))
        p = np.array(d.get("plddt", [0.0]))
        r = man.loc[name]
        rows.append({
            "name": name,
            "target": r["target"],
            "source": r["source"],
            "length": r["length"],
            "mono_pLDDT": p.mean(),
            "mono_pLDDT_min": p.min(),
            "mono_pTM": d.get("ptm", np.nan),
            "complex_ipTM": r["complex_ipTM"],
            "penalty": composition_penalty(r["sequence"]),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("[FAIL] no monomer folds matched the manifest")

    pd.set_option("display.width", 200)
    print("=" * 74)
    print("MONOMER FOLDABILITY -- %d binders folded alone" % len(df))
    print("=" * 74)

    print("\nDoes the binder fold without its target?")
    print("  mono pLDDT : mean=%.1f  median=%.1f  range %.0f-%.0f"
          % (df.mono_pLDDT.mean(), df.mono_pLDDT.median(), df.mono_pLDDT.min(), df.mono_pLDDT.max()))
    for lo, hi, lab in [(90, 101, "very high (>90)"), (70, 90, "confident (70-90)"),
                        (50, 70, "low (50-70)"), (0, 50, "DISORDERED (<50)")]:
        n = ((df.mono_pLDDT >= lo) & (df.mono_pLDDT < hi)).sum()
        print("    %-18s %3d  (%4.1f%%)" % (lab, n, 100 * n / len(df)))

    print("\nBy target:")
    for t, g in df.groupby("target"):
        print("  %-7s n=%3d  mono pLDDT %.1f   complex ipTM %.3f   disordered %.0f%%"
              % (t, len(g), g.mono_pLDDT.mean(), g.complex_ipTM.mean(),
                 100 * (g.mono_pLDDT < 50).mean()))

    print("\n" + "-" * 74)
    print("Does monomer foldability PREDICT the complex score?")
    print("-" * 74)
    for t, g in list(df.groupby("target")) + [("POOLED", df)]:
        if len(g) < 5:
            continue
        rho, p = spearmanr(g.mono_pLDDT, g.complex_ipTM)
        print("  %-7s n=%3d   rho=%+.3f  p=%.4f%s"
              % (t, len(g), rho, p, "   <- significant" if p < 0.05 else ""))

    # practical framing: would a monomer-pLDDT filter cost you good designs?
    print("\n" + "-" * 74)
    print("If you filtered on monomer pLDDT, what would it cost?")
    print("-" * 74)
    good = df.complex_ipTM >= 0.7
    print("  designs with complex ipTM >= 0.7 : %d" % good.sum())
    for thr in (50, 60, 70, 80):
        kept = df.mono_pLDDT >= thr
        lost = (good & ~kept).sum()
        print("    keep mono_pLDDT >= %2d : retains %3d/%3d designs, "
              "loses %d of the %d good ones" % (thr, kept.sum(), len(df), lost, good.sum()))

    df.to_csv("results/monomer_analysis.csv", index=False)
    print("\n[OK] results/monomer_analysis.csv")


if __name__ == "__main__":
    main()
