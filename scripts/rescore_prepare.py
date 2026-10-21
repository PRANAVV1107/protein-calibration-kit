#!/usr/bin/env python3
"""Collect search designs for re-scoring at high accuracy.

The problem this exists to test
-------------------------------
All three searches folded with a DONATED target alignment at 1 model / 3
recycles ("hybrid"). The pre-registered bar they were judged against -- 0.870,
60 blind folds across 10 RBD backbones -- was measured in the variance study at
FULL MSA, 3 models, 8 recycles. Those are different instruments, and on the ten
backbones both have scored they disagree badly:

    backbone     hybrid   high-acc
    binder_42     0.330     0.870
    binder_5      0.410     0.860
    binder_47     0.160     0.520
    binder_13     0.790     0.810

High accuracy reads systematically higher (+0.18 mean on these ten) and reorders
the backbones outright -- binder_42 is third-worst on hybrid and best on
high-accuracy. So "search reached 0.860 against a bar of 0.870" compares a
number from the low-reading instrument against one from the high-reading one.

This does NOT mean the search actually won. The search *optimised* the hybrid
metric, and the two instruments disagree, so designs selected for high hybrid
score need not be high at full MSA. That is exactly why it has to be measured.

What this selects
-----------------
Per run: the top N designs by hybrid ipTM, plus N drawn at random across the
whole score range. The random arm is what lets the instrument mapping be
estimated rather than assumed -- with only the winners we could not tell a real
gain from a systematic offset.
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
RBD = ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
       "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
       "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG")
RUNS = {"run1": "results/rl/mcts_run1",
        "run2": "results/rl/mcts_run2_mean",
        "run3": "results/rl/mcts_run3_uct"}


def collect(run_dir):
    """Every folded design in a run: sequence + its hybrid ipTM."""
    rows = []
    for ev in sorted(glob.glob(os.path.join(run_dir, "eval", "*"))):
        fa = os.path.join(ev, "binders.fasta")
        if not os.path.isfile(fa):
            continue
        seqs, name = {}, None
        for line in open(fa):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif name:
                seqs[name] = line
                name = None
        for f in glob.glob(os.path.join(ev, "folds", "*_scores_rank_001_*.json")):
            tag = os.path.basename(f).split("_scores_rank_001_")[0]
            if tag not in seqs:
                continue
            d = json.load(open(f))
            rows.append({"node_dir": os.path.basename(ev), "tag": tag,
                         "seq": seqs[tag], "hybrid_ipTM": d.get("iptm")})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/rescore")
    ap.add_argument("--top", type=int, default=15, help="top designs per run")
    ap.add_argument("--rand", type=int, default=5, help="random designs per run")
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    os.chdir(KIT)
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    picks = []
    for run, d in RUNS.items():
        df = collect(d)
        if df.empty:
            print("[warn] %s: no designs found" % run)
            continue
        df = df.dropna(subset=["hybrid_ipTM"]).drop_duplicates("seq")
        df["run"] = run
        top = df.nlargest(a.top, "hybrid_ipTM").copy()
        top["arm"] = "top"
        rest = df.drop(top.index)
        k = min(a.rand, len(rest))
        rnd = rest.iloc[rng.choice(len(rest), k, replace=False)].copy()
        rnd["arm"] = "random"
        picks.append(pd.concat([top, rnd]))
        print("  %-5s %4d designs on disk, hybrid ipTM %.3f-%.3f -> %d top + %d random"
              % (run, len(df), df.hybrid_ipTM.min(), df.hybrid_ipTM.max(), len(top), k))

    if not picks:
        raise SystemExit("[FAIL] nothing to re-score")
    sel = pd.concat(picks, ignore_index=True)
    sel = sel.drop_duplicates("seq").reset_index(drop=True)
    sel["id"] = ["%s_%s_%03d" % (r.run, r.arm, i) for i, r in sel.iterrows()]

    fa = os.path.join(a.out, "to_fold.fasta")
    with open(fa, "w") as fh:
        for _, r in sel.iterrows():
            fh.write(">%s\n%s:%s\n" % (r["id"], r["seq"], RBD))
    sel.to_csv(os.path.join(a.out, "manifest.csv"), index=False)

    print("\n[OK] %d unique designs -> %s" % (len(sel), fa))
    print("     %s" % os.path.join(a.out, "manifest.csv"))
    print("\n     at ~240 s/design high-accuracy: about %.1f hours"
          % (len(sel) * 240 / 3600.0))


if __name__ == "__main__":
    main()
