#!/usr/bin/env python3
"""Paired comparison: does the tuned model sample better sequences?

Paired BY BACKBONE. Backbone quality varies more than the effect being measured
-- this project put between-backbone variance at 22-41% of total -- so an
unpaired comparison of 60 sequences against 60 would mostly measure which
backbones happened to land in which arm. Each backbone contributes one paired
difference (tuned mean minus base mean) and the test is on those differences.
"""
import argparse, glob, json, os
import numpy as np, pandas as pd

def mn(x):
    if isinstance(x, dict):
        v=[float(i) for i in x.values()]; return min(v) if v else np.nan
    return float(x) if x is not None else np.nan

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    args = ap.parse_args()
    man = pd.read_csv(os.path.join(args.dir, "manifest.csv")).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(args.dir, "folds", "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index: continue
        d = json.load(open(f)); r = man.loc[tag]
        p = np.array(d.get("plddt",[0.0]))
        rows.append({"arm": r["arm"], "backbone": r["backbone"],
                     "ipTM": d.get("iptm"), "ipSAE_min": mn(d.get("ipsae")),
                     "pDockQ2_min": mn(d.get("pdockq2")), "pLDDT": float(p.mean())})
    if not rows: raise SystemExit("[FAIL] no folds found")
    df = pd.DataFrame(rows)
    print("="*74); print("EXPERT ITERATION ROUND 1 -- tuned vs untuned SAMPLING"); print("="*74)
    print("  folds %d   base %d   tuned %d   backbones %d\n"
          % (len(df), (df.arm=="base").sum(), (df.arm=="tuned").sum(), df.backbone.nunique()))
    for m in ("pDockQ2_min","ipSAE_min","ipTM","pLDDT"):
        piv = df.pivot_table(index="backbone", columns="arm", values=m, aggfunc="mean")
        if not {"base","tuned"}.issubset(piv.columns): continue
        piv = piv.dropna()
        d = (piv["tuned"] - piv["base"]).to_numpy(float)
        n = len(d)
        if n < 3: continue
        mean = d.mean(); se = d.std(ddof=1)/np.sqrt(n)
        from math import erfc, sqrt
        t = mean/se if se>0 else 0.0
        p = erfc(abs(t)/sqrt(2))
        wins = int((d>0).sum())
        print("%-14s base %.3f  tuned %.3f  paired diff %+.3f  se %.3f  p %.3f  (tuned better on %d/%d backbones)"
              % (m, piv["base"].mean(), piv["tuned"].mean(), mean, se, p, wins, n))
    print("\n  Paired by backbone; n = number of backbones, not sequences.")
    print("  A positive difference means the tuned model samples better sequences")
    print("  on structures it was never trained on.")

if __name__ == "__main__":
    main()
