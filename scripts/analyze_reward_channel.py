#!/usr/bin/env python3
"""Does fast-screen ipTM rank SEQUENCES within a fixed backbone?

The RL refinement loop uses fast-screen ipTM as its reward. The fast screen was
shown to carry no signal ACROSS designs (top-10 vs bottom-10, p=0.970) -- but
that varied both backbone and sequence. The RL loop varies only sequence on one
backbone, which is a different question.

This correlates fast-screen against high-accuracy ipTM for the same sequences,
computed WITHIN each backbone (so between-backbone differences cannot inflate
it) and then pooled across backbones.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def load(study_dir, fast_dir, target):
    man = pd.read_csv(os.path.join(study_dir, target, "manifest.csv"))
    key = {f"{r['backbone']}__{r['sample']}": r for _, r in man.iterrows()}

    def scores(d, strip_prefix=""):
        out = {}
        for f in glob.glob(os.path.join(d, "*_scores_rank_001_*.json")):
            n = os.path.basename(f).split("_scores_rank_001_")[0]
            if strip_prefix and n.startswith(strip_prefix):
                n = n[len(strip_prefix):]
            out[n] = json.load(open(f))["iptm"]
        return out

    hi = scores(os.path.join(study_dir, target, "folds"))
    fa = scores(fast_dir)

    rows = []
    for n in set(hi) & set(fa):
        if n in key:
            rows.append({"name": n, "backbone": key[n]["backbone"],
                         "fast": fa[n], "highacc": hi[n]})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-dir", required=True)
    ap.add_argument("--fast-dir", required=True)
    ap.add_argument("--target", required=True)
    args = ap.parse_args()

    df = load(args.study_dir, args.fast_dir, args.target)
    if df.empty:
        raise SystemExit("[FAIL] no matched sequences")

    print(f"\n{'='*74}\nREWARD-CHANNEL TEST -- {args.target}\n{'='*74}")
    print(f"Matched {len(df)} sequences across {df.backbone.nunique()} backbones\n")

    print(f"fast-screen  : mean={df.fast.mean():.3f} sd={df.fast.std():.3f} "
          f"range {df.fast.min():.2f}-{df.fast.max():.2f}")
    print(f"high-accuracy: mean={df.highacc.mean():.3f} sd={df.highacc.std():.3f} "
          f"range {df.highacc.min():.2f}-{df.highacc.max():.2f}")

    # --- the question the RL loop depends on: within-backbone ranking ---
    print(f"\n{'-'*74}\nWITHIN-BACKBONE rank correlation (this is what the RL reward needs)\n{'-'*74}")
    rhos = []
    for bb, g in df.groupby("backbone"):
        if g.fast.nunique() < 2:
            print(f"  {bb:26s} n={len(g)}  fast screen gave IDENTICAL scores -- no ranking possible")
            continue
        rho, p = spearmanr(g.fast, g.highacc)
        rhos.append(rho)
        print(f"  {bb:26s} n={len(g)}  rho={rho:+.3f}  p={p:.3f}   "
              f"fast {g.fast.min():.2f}-{g.fast.max():.2f} -> high {g.highacc.min():.2f}-{g.highacc.max():.2f}")

    if rhos:
        print(f"\n  mean within-backbone rho = {np.mean(rhos):+.3f}  (n={len(rhos)} backbones)")
        print(f"  backbones with rho > 0.5 : {sum(1 for r in rhos if r > 0.5)}/{len(rhos)}")
    else:
        print("\n  No backbone had variable fast-screen scores -- the channel carries no ranking at all.")

    # --- pooled, for reference ---
    rho, p = spearmanr(df.fast, df.highacc)
    print(f"\nPOOLED across all sequences: rho={rho:+.3f}  p={p:.3f}")
    print("  (pooled mixes between-backbone differences in; the within-backbone number above is the relevant one)")

    # --- practical framing: would picking the fast-screen best work? ---
    print(f"\n{'-'*74}\nPRACTICAL: pick the best sequence per backbone BY FAST SCREEN\n{'-'*74}")
    gains = []
    for bb, g in df.groupby("backbone"):
        chosen = g.loc[g.fast.idxmax(), "highacc"]
        best = g.highacc.max()
        avg = g.highacc.mean()
        gains.append({"backbone": bb, "fast_pick": chosen, "oracle_best": best,
                      "random_draw": avg, "lost_vs_oracle": best - chosen})
    gd = pd.DataFrame(gains)
    print(f"  picking by fast screen : mean high-acc {gd.fast_pick.mean():.3f}")
    print(f"  picking at random      : mean high-acc {gd.random_draw.mean():.3f}")
    print(f"  oracle (true best)     : mean high-acc {gd.oracle_best.mean():.3f}")
    edge = gd.fast_pick.mean() - gd.random_draw.mean()
    print(f"\n  fast screen's edge over random selection: {edge:+.3f} ipTM")
    if abs(edge) < 0.02:
        print("  => no better than choosing at random. As an RL reward this is noise.")


if __name__ == "__main__":
    main()
