#!/usr/bin/env python3
"""Did the search actually lose? Re-scored on the bar's own instrument.

Reads the high-accuracy re-folds and answers three things, in order:

1. The instrument mapping. Using the RANDOM arm only -- designs not selected for
   their hybrid score -- how do the two instruments relate? A systematic offset
   here is the thing that would have invalidated the original comparison.
2. The primary, redone. Best high-accuracy ipTM among search designs against
   0.870, the blind-sampling bar, now measured the same way.
3. Whether hybrid ranking survives. Spearman between the two instruments on the
   same designs. If it is low, the search was optimising a metric only loosely
   related to the one it was judged on, which is a finding in itself.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

BAR = 0.870   # blind sampling, 60 folds, 10 RBD backbones, HIGH ACCURACY
PRIOR = {"run1": 0.860, "run2": 0.820, "run3": 0.830}   # as scored on hybrid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/rescore")
    a = ap.parse_args()

    man = pd.read_csv(os.path.join(a.dir, "manifest.csv")).set_index("id")
    rows = []
    for f in glob.glob(os.path.join(a.dir, "folds", "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        d = json.load(open(f))
        r = man.loc[tag]
        rows.append({"id": tag, "run": r["run"], "arm": r["arm"],
                     "hybrid": r["hybrid_ipTM"], "highacc": d.get("iptm"),
                     "pLDDT": np.mean(d.get("plddt", [np.nan]))})
    if not rows:
        raise SystemExit("[FAIL] no re-folds found in %s/folds" % a.dir)
    df = pd.DataFrame(rows).dropna(subset=["highacc"])

    print("=" * 74)
    print("SEARCH DESIGNS RE-SCORED ON THE BAR'S OWN INSTRUMENT")
    print("=" * 74)
    print("  %d designs re-folded at full MSA / 3 models / 8 recycles" % len(df))
    print("  (originally scored with a donated alignment at 1 model / 3 recycles)\n")

    print("-" * 74)
    print("1. INSTRUMENT MAPPING  (random arm only -- not selected on either score)")
    print("-" * 74)
    rnd = df[df.arm == "random"]
    if len(rnd) >= 5:
        d_mean = (rnd.highacc - rnd.hybrid).mean()
        rho, p = stats.spearmanr(rnd.hybrid, rnd.highacc)
        print("  n = %d" % len(rnd))
        print("  hybrid   mean %.3f" % rnd.hybrid.mean())
        print("  high-acc mean %.3f" % rnd.highacc.mean())
        print("  offset        %+.3f   <- how much the instruments differ" % d_mean)
        print("  Spearman rho  %+.3f (p = %.3f)  <- whether ranking survives" % (rho, p))
    else:
        print("  too few random-arm designs (%d) to estimate" % len(rnd))
        d_mean = np.nan

    print("\n" + "-" * 74)
    print("2. THE PRIMARY, REDONE")
    print("-" * 74)
    print("  %-6s %10s %12s %12s   %s" % ("run", "n", "hybrid best", "high-acc best", "vs 0.870"))
    for run, g in df.groupby("run"):
        hb = g.highacc.max()
        print("  %-6s %10d %12.3f %12.3f   %s"
              % (run, len(g), g.hybrid.max(), hb,
                 "CLEARS" if hb > BAR else "below by %.3f" % (BAR - hb)))
    overall = df.highacc.max()
    print("\n  best across all runs, high accuracy: %.3f" % overall)
    print("  blind-sampling bar (same instrument):  %.3f" % BAR)
    print("  -> %s" % ("SEARCH CLEARS THE BAR" if overall > BAR
                       else "search still below by %.3f" % (BAR - overall)))

    print("\n" + "-" * 74)
    print("3. DID THE SEARCH OPTIMISE THE RIGHT THING?")
    print("-" * 74)
    rho, p = stats.spearmanr(df.hybrid, df.highacc)
    print("  Spearman(hybrid, high-acc) over all %d designs: %+.3f (p = %.4f)"
          % (len(df), rho, p))
    top = df[df.arm == "top"]
    if len(top) and len(rnd):
        print("  top arm   (chosen for high hybrid): high-acc mean %.3f" % top.highacc.mean())
        print("  random arm                        : high-acc mean %.3f" % rnd.highacc.mean())
        t, tp = stats.mannwhitneyu(top.highacc, rnd.highacc, alternative="greater")
        print("  Mann-Whitney one-sided p = %.4f" % tp)
        print("\n  If the top arm is not above the random arm, the hybrid score the")
        print("  search maximised does not identify designs that are good on the")
        print("  instrument that matters -- which would make the search's objective,")
        print("  not its algorithm, the thing that failed.")

    df.sort_values("highacc", ascending=False).to_csv(
        os.path.join(a.dir, "rescore_results.csv"), index=False)
    print("\n[OK] %s/rescore_results.csv" % a.dir)


if __name__ == "__main__":
    main()
