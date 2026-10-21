#!/usr/bin/env python3
"""Given a fixed fold budget, how many roots should it be spread across?

The question
------------
Three searches spent ~1000 folds each descending from ONE root. The root sweep
showed between-root spread (E[best-of-6] from 0.371 to 0.827) is larger than
anything search recovered inside a root. That suggests breadth beats depth --
but "suggests" is not a measurement, because the two were never compared at a
matched budget.

This simulates the allocation directly. For a budget of B folds, split as
R roots x K draws each, what is E[best]?

Instrument note
---------------
Uses `results/rl/noise_poc` -- 10 RBD backbones, 16 draws each, folded with a
donated alignment at 1 model / 3 recycles. That is deliberately the SAME
configuration the MCTS runs used, so the search numbers (run 1 0.860, run 2
0.820, run 3 0.830) can be placed on the same axis. Everything here is
therefore on the low-reading instrument; see reports/EXPERIMENTS_QUEUED.md 0a.

The high-accuracy variance study is used as a cross-check on the SHAPE of the
curve only. It has 6 draws per backbone, so bootstrapping a deep single-root
arm from it caps at that root's observed max and is biased against depth.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"


def load_noise_poc(target="rbd"):
    man = pd.read_csv(os.path.join(KIT, "results/rl/noise_poc/manifest.csv")).set_index("tag")
    pools = {}
    for f in glob.glob(os.path.join(KIT, "results/rl/noise_poc/folds",
                                    "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index or man.loc[tag, "target"] != target:
            continue
        v = json.load(open(f)).get("iptm")
        if v is not None:
            pools.setdefault(man.loc[tag, "backbone"], []).append(float(v))
    return {k: np.asarray(v) for k, v in pools.items()}


def load_variance(target="rbd"):
    v = pd.read_csv(os.path.join(KIT, "results/variance_study/variance_results.csv"))
    v = v[v.target == target]
    return {b: g.ipTM.dropna().to_numpy() for b, g in v.groupby("backbone")}


def e_best(pools, R, K, rng, reps=4000):
    """E[max] when spending K draws on each of R randomly chosen roots."""
    names = list(pools)
    if R > len(names):
        return np.nan
    out = np.empty(reps)
    for i in range(reps):
        chosen = rng.choice(len(names), R, replace=False)
        best = -1.0
        for c in chosen:
            p = pools[names[c]]
            best = max(best, rng.choice(p, K, replace=True).max())
        out[i] = best
    return float(out.mean())


def table(pools, label, budgets, rng):
    names = list(pools)
    print("\n" + "=" * 70)
    print(label)
    print("=" * 70)
    print("  %d roots, %d-%d draws each\n"
          % (len(names), min(len(p) for p in pools.values()),
             max(len(p) for p in pools.values())))
    header = "  %-8s" % "budget" + "".join("%10s" % ("R=%d" % r) for r in (1, 2, 5, 10))
    print(header)
    print("  " + "-" * (len(header) - 2))
    rows = []
    for B in budgets:
        line = "  %-8d" % B
        best_r, best_v = None, -1
        for R in (1, 2, 5, 10):
            if R > len(names) or B % R:
                line += "%10s" % "-"
                continue
            K = B // R
            v = e_best(pools, R, K, rng)
            line += "%10.3f" % v
            if v > best_v:
                best_v, best_r = v, R
            rows.append({"budget": B, "roots": R, "draws": K, "e_best": v})
        print(line + "   <- best: R=%d" % best_r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="rbd")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    hyb = load_noise_poc(a.target)
    if not hyb:
        raise SystemExit("[FAIL] no noise_poc folds for %s" % a.target)

    print("=" * 70)
    print("MULTI-ROOT vs SINGLE-ROOT AT MATCHED BUDGET  (%s)" % a.target.upper())
    print("=" * 70)
    print("""
  E[best ipTM] when a budget of B folds is split R roots x (B/R) draws.
  Roots drawn at random each replicate, so this is what you would get
  WITHOUT knowing in advance which root is good.
""")
    df = table(hyb, "HYBRID instrument -- same configuration the searches used",
               [10, 20, 40, 80, 160], rng)

    print("""
  Search, on this same instrument and the same target:
      run 1   0.860   1012 folds, 1 root (binder_47)
      run 2   0.820    955 folds, 1 root
      run 3   0.830   1012 folds, 1 root
""")
    b160 = df[(df.budget == 160)]
    if len(b160):
        r10 = b160[b160.roots == 10].e_best.values
        r1 = b160[b160.roots == 1].e_best.values
        if len(r10) and len(r1):
            print("  At 160 folds: 10 roots gives %.3f, 1 root gives %.3f (%+.3f)"
                  % (r10[0], r1[0], r10[0] - r1[0]))
            print("  The searches spent 6.3x that budget on one root to reach 0.82-0.86.")

    try:
        hi = load_variance(a.target)
        table(hi, "HIGH-ACCURACY cross-check -- shape only, depth arm is capped",
              [10, 20, 40, 60], rng)
        print("""
  Caveat: 6 draws per root here, so the R=1 column bootstraps from 6 values and
  cannot exceed that root's observed max. It is biased AGAINST depth. Read the
  direction, not the magnitude.""")
    except Exception as e:
        print("\n[skip high-accuracy cross-check: %s]" % e)

    df.to_csv(os.path.join(KIT, "results/root_allocation_%s.csv" % a.target), index=False)
    print("\n[OK] results/root_allocation_%s.csv" % a.target)


if __name__ == "__main__":
    main()
