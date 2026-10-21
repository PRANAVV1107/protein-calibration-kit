#!/usr/bin/env python3
"""Cross-target validation of the screening result on an independent dataset.

Why
---
Section 3.1 rests on 164 designs from one competition against one target family
(Adaptyv, EGFR). A single AUC on a single target family is an anecdote. This
re-tests it on a published, independent set.

The dataset
-----------
Overath, Rygaard et al., "Predicting Experimental Success in De Novo Binder
Design: A Meta-Analysis of 3,766 Experimentally Characterised Binders"
(bioRxiv 2025, doi 10.1101/2025.08.14.670059). Data ships in the paper's repo:

    github.com/DigBioLab/de_novo_binder_scoring
      analysis/data/prepared_training_dataset.csv            3,760 designs, 15 targets
      analysis/data/all_affinity_data_with_predicted_metrics.csv   637 designs, 98 K_D

**No folding required.** They re-predicted every complex with AF2 (initial
guess), ColabFold, AF3 and Boltz-1 and publish per-design confidence metrics,
so the comparison is pure analysis against our own numbers.

Two things make this a stronger check than an arbitrary second dataset:

  1. **Our own data is inside it.** `source` includes "Adaptyv binder comp
     round 1" and "round 2", so the Adaptyv subset is a direct replication
     rather than an independent comparison.
  2. **All three of our targets are in it** -- SARS-CoV-2 RBD, IL-7Ra and
     PD-L1 each appear as their own study.

What is compared
----------------
Our figure is ipTM AUC = 0.627 on 164 designs balanced 54/110. Pooled AUC on a
heterogeneous set with an 11.6% binder rate is not comparable to that, so this
reports three things: the pooled number, the per-study spread, and a
base-rate-matched number built by subsampling each study to our 1:2 ratio.
"""
import argparse
import os
import urllib.request

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
RAW = "https://raw.githubusercontent.com/DigBioLab/de_novo_binder_scoring/HEAD/analysis/data"
FILES = {"training": "prepared_training_dataset.csv",
         "affinity": "all_affinity_data_with_predicted_metrics.csv"}
OURS = {"ipTM AUC, full MSA": 0.627, "ipSAE_min AUC, full MSA": 0.618,
        "pDockQ2_min AUC, full MSA": 0.623, "ipSAE_min AUC, cheap screen": 0.466}


def fetch(which, data_dir):
    os.makedirs(data_dir, exist_ok=True)
    p = os.path.join(data_dir, FILES[which])
    if not os.path.exists(p):
        print("  downloading %s ..." % FILES[which])
        urllib.request.urlretrieve("%s/%s" % (RAW, FILES[which]), p)
    return p


def auc(score, label):
    """Rank-sum AUC, the same identity used elsewhere in this project."""
    s = pd.DataFrame({"s": score, "y": label}).dropna()
    if s.y.nunique() < 2:
        return np.nan, len(s), int(s.y.sum())
    r = s.s.rank()
    n1 = int(s.y.sum())
    n0 = len(s) - n1
    return (r[s.y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0), len(s), n1


def balanced_auc(score, label, ratio=2.0, reps=500, seed=0):
    """AUC after subsampling non-binders to `ratio` per binder.

    Our 164-design set is 54 binders to 110 non-binders. Pooled AUC on an
    11.6%-binder set is not the same quantity, so this matches the design.
    """
    rng = np.random.default_rng(seed)
    s = pd.DataFrame({"s": score, "y": label}).dropna()
    pos = s[s.y == 1]
    neg = s[s.y == 0]
    if len(pos) < 3 or len(neg) < 3:
        return np.nan
    k = min(len(neg), int(round(ratio * len(pos))))
    out = []
    for _ in range(reps):
        sub = pd.concat([pos, neg.iloc[rng.choice(len(neg), k, replace=False)]])
        a, _, _ = auc(sub.s, sub.y)
        if not np.isnan(a):
            out.append(a)
    return float(np.mean(out)) if out else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="results/external")
    ap.add_argument("--metric", default="colab_iptm_avg")
    a = ap.parse_args()
    os.chdir(KIT)

    print("=" * 78)
    print("EXTERNAL VALIDATION -- Overath/Rygaard et al. 2025, 3,760 designs")
    print("=" * 78)
    d = pd.read_csv(fetch("training", a.data_dir), low_memory=False)
    print("  %d designs, %d targets, %d studies, %.1f%% binders\n"
          % (len(d), d.target_id.nunique(), d.source.nunique(), 100 * d.binder.mean()))

    # ---- 1. pooled, and why it is not our number ----
    print("-" * 78)
    print("1. POOLED AUC, AND WHY IT IS NOT COMPARABLE TO OURS")
    print("-" * 78)
    pooled, n, nb = auc(d[a.metric], d.binder)
    print("  %-40s %.3f   n=%4d  binders=%3d" % ("pooled over 15 targets", pooled, n, nb))
    print("  %-40s %.3f   n= 164  binders= 54" % ("ours (Adaptyv EGFR, balanced 54/110)", OURS["ipTM AUC, full MSA"]))
    print("""
  Pooled AUC over heterogeneous studies is inflated by between-study
  differences in base rate: this set runs 11.6% binders overall but 1/127 to
  42/91 per study. Our set was deliberately balanced, which lowers AUC. The
  base-rate-matched figure below is the comparable one.""")

    # ---- 2. the overlap: our own data is in theirs ----
    print("\n" + "-" * 78)
    print("2. THE OVERLAP -- our dataset appears in theirs")
    print("-" * 78)
    for src in sorted(s for s in d.source.unique() if "Adaptyv" in str(s)):
        g = d[d.source == src]
        v, n, nb = auc(g[a.metric], g.binder)
        print("  %-34s AUC %.3f   n=%4d  binders=%3d" % (src, v, n, nb))
    ad = d[d.source.str.contains("Adaptyv", na=False)]
    v, n, nb = auc(ad[a.metric], ad.binder)
    bal = balanced_auc(ad[a.metric], ad.binder)
    print("  %-34s AUC %.3f   n=%4d  binders=%3d" % ("both rounds pooled", v, n, nb))
    print("  %-34s AUC %.3f   (subsampled to 1:2, as ours)" % ("both rounds, base-rate matched", bal))

    # ---- 3. per-study spread: the real finding ----
    print("\n" + "-" * 78)
    print("3. PER-STUDY SPREAD -- the figure that matters")
    print("-" * 78)
    rows = []
    for src, g in d.groupby("source"):
        v, n, nb = auc(g[a.metric], g.binder)
        if n >= 20 and not np.isnan(v):
            rows.append({"study": src, "auc": v, "n": n, "binders": nb,
                         "balanced": balanced_auc(g[a.metric], g.binder)})
    r = pd.DataFrame(rows).sort_values("auc", ascending=False)
    print("  %-34s %6s %6s %6s %8s" % ("study", "AUC", "bal", "n", "binders"))
    for _, x in r.iterrows():
        print("  %-34s %6.3f %6s %6d %8d"
              % (x["study"], x["auc"],
                 "%.3f" % x["balanced"] if not np.isnan(x["balanced"]) else "  -  ",
                 x["n"], x["binders"]))
    print("\n  median %.3f   range %.3f - %.3f   across %d studies"
          % (r.auc.median(), r.auc.min(), r.auc.max(), len(r)))
    bal = r.balanced.dropna()
    if len(bal):
        print("  base-rate matched: median %.3f   range %.3f - %.3f"
              % (bal.median(), bal.min(), bal.max()))
    print("  studies below chance (AUC < 0.5): %d of %d" % ((r.auc < 0.5).sum(), len(r)))

    # ---- 4. our three targets, in their data ----
    print("\n" + "-" * 78)
    print("4. OUR THREE TARGETS, MEASURED BY THEM")
    print("-" * 78)
    for key in ("RBD", "IL7Ra", "Pdl1", "EGFR"):
        hits = [s for s in d.source.unique() if key.lower() in str(s).lower()]
        for s_ in sorted(hits):
            g = d[d.source == s_]
            v, n, nb = auc(g[a.metric], g.binder)
            print("  %-34s AUC %.3f   n=%4d  binders=%3d" % (s_, v, n, nb))

    # ---- 5. metric comparison where columns allow ----
    print("\n" + "-" * 78)
    print("5. WHICH METRIC RANKS BEST, ON THEIR DATA")
    print("-" * 78)
    cands = [c for c in ("colab_iptm_avg", "colab_ipae", "af2_pae_interaction",
                         "colab_dockQ", "af3_iptm_avg", "boltz1_iptm_avg")
             if c in d.columns]
    print("  %-26s %6s   %s" % ("metric", "AUC", "note"))
    for c in cands:
        sc = -d[c] if "pae" in c or "ipae" in c else d[c]   # lower pAE is better
        v, n, _ = auc(sc, d.binder)
        print("  %-26s %6.3f   n=%d%s" % (c, v, n, "  (sign flipped)" if sc is not d[c] else ""))
    print("""
  Their headline is that ipSAE outperforms ipTM and ipAE. Our section 3.1 used
  ipSAE_min (0.618) and pDockQ2_min (0.623) alongside ipTM (0.627) and found
  them interchangeable on our set. Those are consistent: a metric difference
  they resolve on 3,760 designs across 15 targets is not resolvable on 164 from
  one target family.""")

    r.to_csv("results/external_validation_per_study.csv", index=False)
    print("\n[OK] results/external_validation_per_study.csv")


if __name__ == "__main__":
    main()
