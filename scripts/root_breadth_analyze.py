#!/usr/bin/env python3
"""Does E[best] climb faster with more ROOTS than with more DRAWS per root?

Four questions, in order of what they settle:

1. **Instrument calibration.** Ten of the fifty backbones are also in the
   variance study at 3 models / 8 recycles. Comparing them tells us whether
   this configuration reads on the same scale, before anything is pooled.
2. **The prediction.** E[best] at fixed total budget, spent as R roots x K
   draws. If the tail argument is right, spending on roots beats spending on
   draws, and the R axis should not flatten by R=50 the way it did by R=10.
3. **The distribution.** Between-root spread at n=50 rather than n=10.
4. **Free root screening.** Whether any cheap structural feature predicts root
   quality -- hopeless at n=10, properly powered at n=50.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"


def load(d):
    man = pd.read_csv(os.path.join(d, "manifest.csv")).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(d, "folds", "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        r = man.loc[tag]
        rows.append({"backbone": r["backbone"], "draw": r["draw"],
                     "overlap": bool(r["overlap"]), "length": r["length"],
                     "ipTM": json.load(open(f)).get("iptm")})
    return pd.DataFrame(rows).dropna(subset=["ipTM"])


def e_best(pools, R, K, rng, reps=3000):
    names = list(pools)
    if R > len(names):
        return np.nan
    out = np.empty(reps)
    for i in range(reps):
        chosen = rng.choice(len(names), R, replace=False)
        out[i] = max(rng.choice(pools[names[c]], K, replace=True).max() for c in chosen)
    return float(out.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="results/root_breadth")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    os.chdir(KIT)
    rng = np.random.default_rng(a.seed)

    df = load(a.dir)
    if df.empty:
        raise SystemExit("[FAIL] no folds in %s/folds" % a.dir)
    pools = {b: g.ipTM.to_numpy() for b, g in df.groupby("backbone")}
    per = df.groupby("backbone").size()

    print("=" * 74)
    print("ROOT-BREADTH SWEEP")
    print("=" * 74)
    print("  %d backbones, %d-%d draws each, %d folds"
          % (len(pools), per.min(), per.max(), len(df)))
    print("  full MSA, 1 model, 3 recycles\n")

    # ---- 1. instrument calibration on the overlap set ----
    print("-" * 74)
    print("1. INSTRUMENT CALIBRATION (the 10 backbones also in the variance study)")
    print("-" * 74)
    v = pd.read_csv("results/variance_study/variance_results.csv")
    v = v[v.target == "rbd"].groupby("backbone").ipTM.agg(["max", "mean"])
    ov = df[df.overlap].groupby("backbone").ipTM.agg(["max", "mean"])
    j = ov.join(v, lsuffix="_new", rsuffix="_var", how="inner")
    if len(j) >= 4:
        print("  %-12s %10s %10s   %10s %10s"
              % ("backbone", "max new", "max var", "mean new", "mean var"))
        for b, r in j.sort_values("max_var", ascending=False).iterrows():
            print("  %-12s %10.3f %10.3f   %10.3f %10.3f"
                  % (b, r["max_new"], r["max_var"], r["mean_new"], r["mean_var"]))
        rho, p = stats.spearmanr(j["max_new"], j["max_var"])
        print("\n  n=%d  offset on max %+.3f   Spearman %+.3f (p=%.3f)"
              % (len(j), (j["max_new"] - j["max_var"]).mean(), rho, p))
        print("  A near-zero offset and a high rho mean the two configurations")
        print("  are interchangeable and results may be pooled. Otherwise not.")
    else:
        print("  only %d overlap backbones folded so far" % len(j))

    # ---- 2. the prediction: roots versus draws ----
    print("\n" + "-" * 74)
    print("2. E[best] AT FIXED BUDGET: spend on ROOTS or on DRAWS?")
    print("-" * 74)
    Kmax = int(per.min())
    budgets = [20, 40, 80, 160, 200]
    combos = [(R, K) for K in (1, 2, 4) for R in (5, 10, 25, 50)]
    print("  %-8s" % "budget" + "".join("%12s" % ("R=%d,K=%d" % c) for c in combos[:6]))
    print("  " + "-" * 78)
    rows = []
    for B in budgets:
        line = "  %-8d" % B
        for R, K in combos[:6]:
            if R * K != B or K > Kmax or R > len(pools):
                line += "%12s" % "-"
                continue
            val = e_best(pools, R, K, rng)
            line += "%12.3f" % val
            rows.append({"budget": B, "roots": R, "draws": K, "e_best": val})
        print(line)

    print("\n  R axis at fixed K (the prediction):")
    for K in (1, 2, 4):
        if K > Kmax:
            continue
        vals = []
        for R in (5, 10, 25, 50):
            if R > len(pools):
                continue
            vals.append("R=%-3d %.3f" % (R, e_best(pools, R, K, rng)))
        print("    K=%d   %s" % (K, "   ".join(vals)))
    print("\n  K axis at fixed R (the comparison):")
    for R in (5, 10, 25, 50):
        if R > len(pools):
            continue
        vals = []
        for K in range(1, Kmax + 1):
            vals.append("K=%-2d %.3f" % (K, e_best(pools, R, K, rng)))
        print("    R=%-3d %s" % (R, "   ".join(vals)))
    print("""
  If E[best] rises faster along R than along K at equal budget, breadth beats
  depth and the tail argument holds. If the R axis flattens by R=50, the
  recommendation has a ceiling and the paper must say where it is.""")

    # ---- 3. the distribution ----
    print("\n" + "-" * 74)
    print("3. BETWEEN-ROOT DISTRIBUTION (n=%d)" % len(pools))
    print("-" * 74)
    bb = df.groupby("backbone").ipTM.agg(["max", "mean", "std"]).sort_values("max", ascending=False)
    print("  best-of-%d per root: min %.3f  median %.3f  max %.3f  sd %.3f"
          % (Kmax, bb["max"].min(), bb["max"].median(), bb["max"].max(), bb["max"].std()))
    print("  fraction of roots reaching >0.7: %.0f%%   >0.8: %.0f%%"
          % (100 * (bb["max"] > 0.7).mean(), 100 * (bb["max"] > 0.8).mean()))
    print("\n  top 8 roots: %s" % ", ".join("%s(%.2f)" % (b, r["max"]) for b, r in bb.head(8).iterrows()))
    print("  bottom 5   : %s" % ", ".join("%s(%.2f)" % (b, r["max"]) for b, r in bb.tail(5).iterrows()))

    # ---- 4. can a good root be spotted for free? ----
    print("\n" + "-" * 74)
    print("4. DOES ANY CHEAP FEATURE PREDICT ROOT QUALITY? (n=%d, now powered)" % len(bb))
    print("-" * 74)
    fp = "results/root_features_all.csv"
    if os.path.exists(fp):
        feat = pd.read_csv(fp)
        feat = feat[feat.batch == 0].set_index("backbone")
        m = bb.join(feat, how="inner")
        if len(m) >= 20:
            print("  %-14s %9s %9s   %s" % ("feature", "rho", "p", ""))
            print("  " + "-" * 46)
            for c in ("length", "rog", "contacts", "hotspot_c", "helix_frac", "contact_dens"):
                if c not in m:
                    continue
                rho, p = stats.spearmanr(m[c], m["max"])
                flag = "**" if p < 0.05 / 6 else ("*" if p < 0.05 else "")
                print("  %-14s %9.3f %9.3f   %s" % (c, rho, p, flag))
            print("\n  * p<0.05 uncorrected, ** survives Bonferroni over 6 features.")
            print("  At n=%d, |rho| >= %.2f reaches p<0.05 uncorrected."
                  % (len(m), 1.96 / np.sqrt(len(m) - 1)))
        else:
            print("  only %d backbones joined to features" % len(m))
    else:
        print("  run scripts/root_anatomy.py first to generate %s" % fp)

    bb.to_csv("results/root_breadth_per_root.csv")
    if rows:
        pd.DataFrame(rows).to_csv("results/root_breadth_allocation.csv", index=False)
    print("\n[OK] results/root_breadth_per_root.csv")


if __name__ == "__main__":
    main()
