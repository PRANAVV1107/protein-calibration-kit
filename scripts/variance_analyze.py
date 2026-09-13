#!/usr/bin/env python3
"""Analyse the sequence-draw variance study.

Reports, per target and pooled:
  * within-backbone spread of high-accuracy ipTM across independent draws
  * the gain from best-of-N over a single draw (the quantity that matters if
    you are deciding how many sequences to sample per backbone)
  * how much of the total ipTM variance sits within backbones vs between them
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd


def collect(study_dir):
    rows = []
    for target in ("il7ra", "rbd", "pdl1"):
        man_path = os.path.join(study_dir, target, "manifest.csv")
        fold_dir = os.path.join(study_dir, target, "folds")
        if not (os.path.exists(man_path) and os.path.isdir(fold_dir)):
            continue
        man = pd.read_csv(man_path)
        keyed = {f"{r.backbone}__{r['sample']}": r for _, r in man.iterrows()}

        for f in glob.glob(os.path.join(fold_dir, "*_scores_rank_001_*.json")):
            name = os.path.basename(f).split("_scores_rank_001_")[0]
            if name not in keyed:
                continue
            d = json.load(open(f))
            plddt = d.get("plddt", [0])
            r = keyed[name]
            rows.append({
                "target": target,
                "backbone": r.backbone,
                "sample": r["sample"],
                "ipTM": d["iptm"],
                "pTM": d.get("ptm"),
                "mean_pLDDT": sum(plddt) / len(plddt),
                "penalty": r.penalty,
                "length": r.length,
            })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-dir", required=True)
    args = ap.parse_args()

    df = collect(args.study_dir)
    if df.empty:
        raise SystemExit("[FAIL] no fold results matched any manifest")

    out_csv = os.path.join(args.study_dir, "variance_results.csv")
    df.to_csv(out_csv, index=False)
    pd.set_option("display.width", 200)

    print(f"\nCollected {len(df)} folds across {df['target'].nunique()} targets, "
          f"{df.groupby(['target','backbone']).ngroups} backbones\n")

    # ---- per-backbone spread ----
    g = df.groupby(["target", "backbone"])["ipTM"]
    per_bb = g.agg(n="size", best="max", worst="min", mean="mean", std="std").reset_index()
    per_bb["spread"] = per_bb["best"] - per_bb["worst"]
    per_bb["best_minus_mean"] = per_bb["best"] - per_bb["mean"]

    print("=" * 78)
    print("WITHIN-BACKBONE ipTM SPREAD ACROSS DRAWS")
    print("=" * 78)
    for t, grp in per_bb.groupby("target"):
        print(f"\n{t}  ({len(grp)} backbones, {int(grp['n'].mean())} draws each)")
        print(f"  spread (best-worst): mean={grp['spread'].mean():.3f}  "
              f"median={grp['spread'].median():.3f}  max={grp['spread'].max():.3f}")
        print(f"  within-backbone std: mean={grp['std'].mean():.3f}")
        print(f"  best-of-N vs mean draw: +{grp['best_minus_mean'].mean():.3f}")

    print(f"\nPOOLED ({len(per_bb)} backbones)")
    print(f"  spread (best-worst): mean={per_bb['spread'].mean():.3f}  "
          f"median={per_bb['spread'].median():.3f}  max={per_bb['spread'].max():.3f}")
    print(f"  within-backbone std: mean={per_bb['std'].mean():.3f}")
    print(f"  best-of-N vs mean draw: +{per_bb['best_minus_mean'].mean():.3f}")

    # ---- variance decomposition ----
    print("\n" + "=" * 78)
    print("VARIANCE DECOMPOSITION (is the backbone or the sequence in control?)")
    print("=" * 78)
    for t, grp in df.groupby("target"):
        bb_means = grp.groupby("backbone")["ipTM"].mean()
        between = bb_means.var(ddof=1)
        within = grp.groupby("backbone")["ipTM"].var(ddof=1).mean()
        tot = between + within
        if tot > 0:
            print(f"  {t:6s}: between-backbone {between/tot:5.1%}   "
                  f"within-backbone (sequence) {within/tot:5.1%}")

    # ---- expected value of sampling N draws ----
    print("\n" + "=" * 78)
    print("EXPECTED ipTM FROM BEST-OF-N (bootstrap over observed draws)")
    print("=" * 78)
    rng = np.random.default_rng(0)
    max_n = int(per_bb["n"].min())
    print(f"{'N draws':>8s}  " + "  ".join(f"{t:>8s}" for t in sorted(df['target'].unique())) + "    pooled")
    for n in range(1, max_n + 1):
        line = f"{n:>8d}  "
        pooled = []
        for t in sorted(df["target"].unique()):
            vals = []
            for _, grp in df[df["target"] == t].groupby("backbone"):
                arr = grp["ipTM"].to_numpy()
                draws = rng.choice(arr, size=(400, n), replace=True).max(axis=1)
                vals.append(draws.mean())
            line += f"{np.mean(vals):8.3f}  "
            pooled.extend(vals)
        line += f"  {np.mean(pooled):.3f}"
        print(line)

    per_bb.to_csv(os.path.join(args.study_dir, "per_backbone_summary.csv"), index=False)
    print(f"\n[OK] {out_csv}")
    print(f"[OK] {os.path.join(args.study_dir, 'per_backbone_summary.csv')}")


if __name__ == "__main__":
    main()
