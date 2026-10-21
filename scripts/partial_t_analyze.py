#!/usr/bin/env python3
"""Correlation length of the score landscape under partial re-diffusion.

One number per partial_T: the Spearman correlation between a parent's value and
its best child's, pooled over parents. Compare against run 3's rho = 0.168
(p = 0.384, n = 42) measured at partial_T = 10.
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd
from scipy import stats


def run_values(d):
    """Parent own-value and each child's own-value, parsed from a run's log."""
    log = os.path.join(d, "search.log")
    if not os.path.isfile(log):
        return None, []
    root, kids = None, []
    for line in open(log, errors="ignore"):
        m = re.search(r"root score ([0-9.]+)", line)
        if m:
            root = float(m.group(1))
        m = re.search(r"child \d+ -> value ([0-9.]+)", line)
        if m:
            kids.append(float(m.group(1)))
    return root, kids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/rl/partial_t")
    a = ap.parse_args()

    rows = []
    for d in sorted(glob.glob(os.path.join(a.dir, "t*_p*"))):
        m = re.search(r"t(\d+)_p(\d+)$", os.path.basename(d))
        if not m:
            continue
        T, p = int(m.group(1)), int(m.group(2))
        root, kids = run_values(d)
        if root is None or not kids:
            continue
        rows.append({"T": T, "parent_id": p, "parent": root,
                     "best_child": max(kids), "mean_child": float(np.mean(kids)),
                     "n_kids": len(kids)})
    if not rows:
        raise SystemExit("[FAIL] no completed sweeps in %s" % a.dir)
    df = pd.DataFrame(rows)

    print("=" * 72)
    print("SCORE-LANDSCAPE LOCALITY vs partial_T")
    print("=" * 72)
    print("  Does a parent's value predict its children's?")
    print("  Reference: run 3 at partial_T=10 gave rho = 0.168, p = 0.384, n = 42\n")
    print("  %-6s %7s %10s %10s %9s %9s" %
          ("T", "parents", "rho", "p", "mean gap", "kid>parent"))
    print("  " + "-" * 56)
    for T, g in df.groupby("T"):
        if len(g) < 3:
            print("  %-6d %7d   (too few parents)" % (T, len(g)))
            continue
        rho, p = stats.spearmanr(g.parent, g.best_child)
        gap = (g.best_child - g.parent).mean()
        win = 100.0 * (g.best_child > g.parent).mean()
        print("  %-6d %7d %10.3f %10.3f %9.3f %8.0f%%" % (T, len(g), rho, p, gap, win))

    # ------------------------------------------------------------------
    # Does the gain from expanding depend on how good the parent already is?
    #
    # Run 3 measured this inside one tree: top-quartile parents averaged 0.342
    # and their best children 0.279 -- they LOST -- while average parents gained
    # about +0.05. If that holds across a range of starting qualities then a
    # search seeded from good roots gains nothing, and multi-root search
    # (Experiment 5) is answered without running it.
    #
    # This is the cross-parent version: parents span whatever quality the fixed
    # index rule drew, rather than the narrow range inside a single tree.
    # ------------------------------------------------------------------
    print("")
    print("=" * 72)
    print("DOES GAIN DEPEND ON PARENT QUALITY?")
    print("=" * 72)
    print("  Reference, run 3 within one tree (n=42):")
    print("    top-quartile parents 0.342 -> best child 0.279   (lost 0.063)")
    print("    average parents                                  (gained ~0.05)")
    print("")
    print("  %-6s %7s %10s %10s %11s" % ("T", "parents", "rho", "p", "gain slope"))
    print("  " + "-" * 50)
    for T, g in df.groupby("T"):
        if len(g) < 4:
            print("  %-6d %7d   (too few parents)" % (T, len(g)))
            continue
        gain = g.best_child - g.parent
        rho, p = stats.spearmanr(g.parent, gain)
        slope = np.polyfit(g.parent, gain, 1)[0] if g.parent.nunique() > 1 else float("nan")
        print("  %-6d %7d %10.3f %10.3f %11.3f" % (T, len(g), rho, p, slope))
    if len(df) >= 6:
        gain_all = df.best_child - df.parent
        rho, p = stats.spearmanr(df.parent, gain_all)
        print("")
        print("  pooled over all T: rho = %+.3f (p = %.4f, n = %d)" % (rho, p, len(df)))
        hi = df.parent >= df.parent.quantile(0.5)
        print("  better-than-median parents: mean gain %+.3f" % gain_all[hi].mean())
        print("  worse-than-median parents : mean gain %+.3f" % gain_all[~hi].mean())
        print("")
        print("  Strongly negative rho means the search only helps when the starting")
        print("  point is bad, so seeding it from good roots buys nothing and")
        print("  Experiment 5 would be measuring regression to the mean. Flat or")
        print("  positive rho means the gain is real across the quality range and")
        print("  multi-root search is worth building.")
    print("")
    print("\n" + "=" * 72)
    print("READ THIS AS")
    print("=" * 72)
    print("""  rho climbing as T falls   -> the landscape has a correlation length and
                               partial_T = 10 was simply past it. Re-run the
                               search at the T where rho is highest.
  rho flat and null at all T -> no locality at any step size tested. Tree search
                               over these backbones cannot work, and the reason
                               is measured rather than inferred.

  Caveat that belongs in any writeup: with %d parents per T this detects only
  large correlations. A null here bounds the effect, it does not prove zero.""" % df.groupby("T").size().min())
    df.to_csv("results/partial_t_pairs.csv", index=False)
    print("\n[OK] results/partial_t_pairs.csv")


if __name__ == "__main__":
    main()
