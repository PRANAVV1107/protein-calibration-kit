#!/usr/bin/env python3
"""How many sequence draws does a backbone need before its score means anything?

The question and why it is the prerequisite
--------------------------------------------
Guided search over backbones failed here: 1012 folds reached ipTM 0.860 against
0.870 for blind sampling at 60 folds. The diagnosis from that run was that node
scores were noise -- parent score did not predict child quality, expanding the
best node failed repeatedly, and depth carried no information.

This measures that noise directly. Scoring a backbone means drawing K sequences
with ProteinMPNN, folding them, and taking the max. If two disjoint sets of K
draws from the SAME backbone disagree by more than the differences the search
was trying to resolve, then the search was selecting on measurement error.

From the variance study (6 draws per backbone) the disagreement at K=3 is
already known and is target-dependent:

    RBD      mean 0.192   median 0.160   p90 0.481
    IL-7Ra   mean 0.027   median 0.010
    PD-L1    mean 0.024   median 0.010

RBD -- the target the search ran on -- is roughly 7x noisier than the others.
What is NOT known is how fast this falls with K, and that is what decides
whether the approach is fixable.

Two outcomes, both worth having
-------------------------------
  noise falls as ~1/sqrt(K)   sampling error; a larger K fixes it and the cost
                              is calculable
  noise plateaus              irreducible backbone-level spread; no compute
                              budget fixes it, and tree search over backbones is
                              not viable for this pipeline

What this does NOT show
-----------------------
It does not test whether search beats random sampling. That needs a separate
run with a matched baseline. This establishes only whether node evaluation can
be made reliable enough for such a test to be meaningful.
"""
import argparse
import glob
import itertools
import json
import os
import re

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"


def mn(x):
    if isinstance(x, dict):
        v = [float(i) for i in x.values()]
        return min(v) if v else np.nan
    return float(x) if x is not None else np.nan


def disagreement(vals, k, rng, reps=400):
    """|max(A) - max(B)| for disjoint random subsets of size k."""
    v = np.asarray(vals, float)
    v = v[~np.isnan(v)]
    if len(v) < 2 * k:
        return np.nan, np.nan, np.nan
    out = []
    idx = np.arange(len(v))
    for _ in range(reps):
        p = rng.permutation(idx)
        out.append(abs(v[p[:k]].max() - v[p[k:2 * k]].max()))
    o = np.array(out)
    return o.mean(), float(np.median(o)), float(np.percentile(o, 90))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="noise_poc output dir")
    ap.add_argument("--metric", default="ipTM", choices=["ipTM", "ipSAE_min"])
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.chdir(KIT)
    rng = np.random.default_rng(args.seed)

    man = pd.read_csv(os.path.join(args.dir, "manifest.csv")).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(args.dir, "folds", "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        d = json.load(open(f))
        r = man.loc[tag]
        rows.append({"target": r["target"], "backbone": r["backbone"],
                     "ipTM": d.get("iptm"), "ipSAE_min": mn(d.get("ipsae"))})
    if not rows:
        raise SystemExit("[FAIL] no folds found in %s/folds" % args.dir)
    df = pd.DataFrame(rows)

    print("=" * 72)
    print("NODE-SCORE NOISE vs NUMBER OF SEQUENCE DRAWS")
    print("=" * 72)
    print("  metric %s" % args.metric)
    for t, g in df.groupby("target"):
        per = g.groupby("backbone").size()
        print("  %-7s %d backbones, %d-%d draws each, %d folds"
              % (t, g.backbone.nunique(), per.min(), per.max(), len(g)))
    print()
    print("  |max(K draws) - max(K other draws)| on the SAME backbone.")
    print("  Lower is better. Compare against the score differences a search")
    print("  would try to resolve -- typically 0.05-0.15.\n")

    print("%-8s %5s %8s %8s %8s   %s" % ("target", "K", "mean", "median", "p90", "vs K=2"))
    print("-" * 72)
    res = []
    for t, g in df.groupby("target"):
        base = None
        for k in (2, 4, 8):
            vals = [disagreement(h[args.metric].to_numpy(), k, rng)
                    for _, h in g.groupby("backbone")]
            vals = [v for v in vals if not np.isnan(v[0])]
            if not vals:
                continue
            m = float(np.mean([v[0] for v in vals]))
            md = float(np.mean([v[1] for v in vals]))
            p9 = float(np.mean([v[2] for v in vals]))
            if base is None:
                base = m
            print("%-8s %5d %8.3f %8.3f %8.3f   %s"
                  % (t, k, m, md, p9,
                     "--" if k == 2 else "%.2fx" % (m / base if base else float("nan"))))
            res.append({"target": t, "K": k, "mean": m, "median": md, "p90": p9})

    print("\n" + "=" * 72)
    print("READ THIS AS")
    print("=" * 72)
    print("  ~1/sqrt(K) decay  -> sampling error, a larger K fixes it")
    print("                       (K=8 should be about 0.50x the K=2 noise)")
    print("  plateau           -> irreducible backbone spread; no budget fixes it,")
    print("                       and tree search over backbones is not viable here")
    print("\n  This does NOT show whether search beats random sampling.")
    print("  It shows whether node evaluation can be made reliable enough")
    print("  for that question to be worth asking.")
    if res:
        pd.DataFrame(res).to_csv("results/noise_curve.csv", index=False)
        print("\n[OK] results/noise_curve.csv")


if __name__ == "__main__":
    main()
