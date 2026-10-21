#!/usr/bin/env python3
"""Is 0.86 a ceiling of the pipeline, or a ceiling of binder_47?

Why this exists
---------------
Three searches from the root `binder_47` topped out at 0.820-0.860, and a
bootstrap over 2,104 pooled folds put expected best-of-N at 0.849 by N=1000,
rising ~0.01 per doubling after. That was reported as "the pipeline's reachable
distribution caps near 0.86".

It is not. Those 2,104 folds are overwhelmingly binder_47 and its descendants,
so the number describes ONE root. The variance study already hinted at the
problem: per-backbone best-of-6 on RBD ran from 0.460 to 0.870 across ten
backbones, a spread wider than anything search recovered inside one of them.

This measures the between-root spread directly, and compares it against what
search bought inside a single root. If between-root variance dominates, then
the cheapest available improvement is to draw more roots rather than to search
harder from one -- and every search result in this project is conditional on a
root that was never chosen on evidence.

Instrument note -- this bit matters
-----------------------------------
An earlier version of this script read `results/rl/noise_poc`, which folds with
a DONATED target alignment at 1 model / 3 recycles. That is the configuration
§3.2 measured at AUC 0.545, the weakest of the three, and on the ten backbones
both datasets cover it disagrees with high accuracy badly:

    backbone     hybrid   high-acc
    binder_42     0.330     0.870
    binder_5      0.410     0.860
    binder_47     0.160     0.520

Ranking roots on that instrument ranks them on noise. This version reads the
variance study instead -- full MSA, 3 models, 8 recycles, 6 draws per backbone.
Fewer draws, but measured with the instrument the 0.870 bar was measured on.

No new folds. Everything here is already on disk.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"


def mn(x):
    if isinstance(x, dict):
        v = [float(i) for i in x.values()]
        return min(v) if v else np.nan
    return float(x) if x is not None else np.nan


def load_manifest_dir(folds, manifest, target):
    """noise_poc / variance_study layout: manifest.csv maps tag -> backbone."""
    man = pd.read_csv(manifest).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(folds, "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        r = man.loc[tag]
        if target and r["target"] != target:
            continue
        d = json.load(open(f))
        rows.append({"backbone": r["backbone"], "ipTM": d.get("iptm"),
                     "ipSAE_min": mn(d.get("ipsae"))})
    return pd.DataFrame(rows)


def best_of_k(vals, k, rng, reps=2000):
    """Expected max of k draws, with replacement, from this backbone's pool."""
    v = np.asarray(vals, float)
    v = v[~np.isnan(v)]
    if len(v) == 0:
        return np.nan
    return float(np.mean([rng.choice(v, k, replace=True).max() for _ in range(reps)]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", default="ipTM", choices=["ipTM", "ipSAE_min"])
    ap.add_argument("--target", default="rbd")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.chdir(KIT)
    rng = np.random.default_rng(a.seed)

    v = pd.read_csv("results/variance_study/variance_results.csv")
    df = v[v.target == a.target][["backbone", "ipTM"]].copy()
    if "ipSAE_min" not in df.columns:
        df["ipSAE_min"] = np.nan
    if df.empty:
        raise SystemExit("[FAIL] no high-accuracy folds for target %s" % a.target)

    print("=" * 74)
    print("BETWEEN-ROOT SPREAD vs WITHIN-ROOT SEARCH GAIN  (%s, %s)"
          % (a.target.upper(), a.metric))
    print("=" * 74)
    per = df.groupby("backbone")[a.metric]
    n_per = per.size()
    print("  %d backbones, %d-%d draws each, %d folds total\n"
          % (df.backbone.nunique(), n_per.min(), n_per.max(), len(df)))

    K = int(n_per.min())
    print("  Per-root, at K=%d draws (the number every root actually has):\n" % K)
    print("  %-12s %7s %7s %7s %8s" % ("root", "median", "mean", "max", "E[best-of-%d]" % K))
    print("  " + "-" * 56)
    rows = []
    for b, g in df.groupby("backbone"):
        v = g[a.metric].dropna().to_numpy()
        e = best_of_k(v, K, rng)
        rows.append({"backbone": b, "median": np.median(v), "mean": v.mean(),
                     "max": v.max(), "ebest": e})
    rows = sorted(rows, key=lambda r: -r["max"])
    for r in rows:
        print("  %-12s %7.3f %7.3f %7.3f %8.3f"
              % (r["backbone"], r["median"], r["mean"], r["max"], r["ebest"]))

    mx = np.array([r["max"] for r in rows])
    eb = np.array([r["ebest"] for r in rows])
    print("\n  " + "-" * 56)
    print("  between-root spread of best-of-%d:  %.3f to %.3f  (range %.3f, sd %.3f)"
          % (K, eb.min(), eb.max(), eb.max() - eb.min(), eb.std(ddof=1)))

    print("\n" + "=" * 74)
    print("THE COMPARISON THAT MATTERS")
    print("=" * 74)
    print("""
  Search, from ONE root (binder_47), 1012 folds:
      run 1  0.860     run 2  0.820     run 3  0.830
      blind sampling on that root, 60 folds: 0.870

  Drawing a DIFFERENT root costs %d folds to evaluate at K=%d.
""" % (K, K))
    # what does picking the best of j random roots buy?
    print("  %-26s %8s" % ("strategy", "E[best]"))
    print("  " + "-" * 38)
    for j in (1, 2, 3, 5, len(eb)):
        vals = [np.max(rng.choice(eb, j, replace=True)) for _ in range(4000)]
        print("  %-26s %8.3f"
              % ("best of %d random roots" % j if j > 1 else "one random root",
                 float(np.mean(vals))))
    print("  %-26s %8.3f" % ("best root observed", eb.max()))

    print("""
  Read this as: if the spread across roots at a FIXED small budget exceeds what
  search recovered inside one root at 1012 folds, then root selection is the
  cheaper lever and every search conclusion here is conditional on an arbitrary
  starting point.
""")
    out = pd.DataFrame(rows)
    out.to_csv("results/root_sweep_%s.csv" % a.target, index=False)
    print("[OK] results/root_sweep_%s.csv" % a.target)


if __name__ == "__main__":
    main()
