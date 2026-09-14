#!/usr/bin/env python3
"""The controlled comparison: cheap screen vs high accuracy on the SAME designs.

Why this script exists rather than reading two AUCs off two reports:

The cheap arm covers 365 designs, the high-accuracy arm 165. Comparing 0.518
against whatever the 165 gives is not a comparison of accuracy settings -- it
is partly a comparison of two different design sets, which differ in class
balance (310:55 against 110:55). So both arms are restricted here to the 165
designs that have folds in both, and the difference is tested PAIRED, by
resampling designs rather than scores, which is what makes the CI on the
difference meaningful.

Three baselines are reported alongside, and each exists to kill a specific
alternative explanation:

  length        Binders in this set average 126.7 aa against 100.5 aa for
                non-binders. Any metric that tracks chain length inherits an
                AUC above 0.5 that has nothing to do with modelling the
                interface. If a metric cannot beat length, it has not earned
                its compute.
  adaptyv_iptm  The competition's own reported ipTM -- an external pipeline on
                the same designs. This is the bar.
  random        0.5, stated explicitly so the reader is not asked to eyeball
                whether an interval clears it.

Labels are normalised by explicit text comparison. pandas reads the strings
'true'/'false' back from CSV as booleans, and a truthiness test on that has
already silently inverted every label in this project once.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"


def mn(x):
    if isinstance(x, dict):
        v = [float(i) for i in x.values()]
        return min(v) if v else np.nan
    return float(x) if x is not None else np.nan


def auc(pos, neg):
    """ROC AUC via the rank-sum identity; ties get average ranks."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    m = ~np.isnan(pos)
    pos = pos[m]
    neg = neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank(method="average").to_numpy()
    r_pos = ranks[:len(pos)].sum()
    return (r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def auc_df(df, col):
    return auc(df.loc[df.binding == "true", col], df.loc[df.binding == "false", col])


def load(fold_dir, manifest, prefix):
    """Scores from one arm, keyed by tag, columns prefixed."""
    man = pd.read_csv(manifest).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(fold_dir, "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        d = json.load(open(f))
        p = np.array(d.get("plddt", [0.0]))
        rows.append({
            "tag": tag,
            prefix + "ipTM": d.get("iptm"),
            prefix + "pTM": d.get("ptm"),
            prefix + "ipSAE_min": mn(d.get("ipsae")),
            prefix + "pDockQ2_min": mn(d.get("pdockq2")),
            prefix + "pLDDT": p.mean(),
        })
    return pd.DataFrame(rows).set_index("tag")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cheap-dir", default="results/queue/adaptyv_val_folds")
    ap.add_argument("--cheap-man", default="results/queue/adaptyv_val_manifest.csv")
    ap.add_argument("--high-dir", default="results/queue/adaptyv_highacc_folds")
    ap.add_argument("--high-man", default="results/queue/adaptyv_highacc_manifest.csv")
    ap.add_argument("--out", default="results/arm_comparison.csv")
    ap.add_argument("--boot", type=int, default=4000)
    ap.add_argument("--min-folds", type=int, default=40,
                    help="refuse to report below this many high-accuracy folds")
    args = ap.parse_args()

    os.chdir(KIT)
    hman = pd.read_csv(args.high_man).set_index("tag")
    cheap = load(args.cheap_dir, args.cheap_man, "cheap_")
    high = load(args.high_dir, args.high_man, "high_")

    meta = hman[["binding", "length", "composition_penalty",
                 "adaptyv_iptm", "adaptyv_plddt"]].copy()
    meta["binding"] = meta["binding"].map(
        lambda v: "true" if str(v).strip().lower() in ("true", "1")
        else ("false" if str(v).strip().lower() in ("false", "0") else "unknown"))

    df = meta.join(cheap, how="inner").join(high, how="inner")
    df = df[df.binding.isin(["true", "false"])]

    n_pos = (df.binding == "true").sum()
    n_neg = (df.binding == "false").sum()
    print("=" * 78)
    print("CONTROLLED COMPARISON -- identical designs, cheap screen vs high accuracy")
    print("=" * 78)
    print("  paired designs: %d   binders %d   non-binders %d"
          % (len(df), n_pos, n_neg))
    print("  high-accuracy arm is %d / 165 complete" % len(high))
    if len(df) < args.min_folds or n_pos < 10 or n_neg < 10:
        print("\n  TOO EARLY -- need >=%d paired and >=10 per class. Rerun later."
              % args.min_folds)
        return
    if len(high) < 165:
        print("  NOTE: arm incomplete, so this is an interim read on whichever")
        print("        designs finished first, not a random subsample.")

    df.to_csv(args.out)

    # ---- baselines and metrics on the same paired set ----
    cand = [("length (confound baseline)", "length"),
            ("adaptyv_iptm (external)", "adaptyv_iptm"),
            ("cheap  ipTM", "cheap_ipTM"),
            ("cheap  ipSAE_min", "cheap_ipSAE_min"),
            ("cheap  pDockQ2_min", "cheap_pDockQ2_min"),
            ("HIGH   ipTM", "high_ipTM"),
            ("HIGH   ipSAE_min", "high_ipSAE_min"),
            ("HIGH   pDockQ2_min", "high_pDockQ2_min"),
            ("HIGH   pLDDT", "high_pLDDT")]
    cand = [(lab, c) for lab, c in cand if c in df.columns and df[c].notna().any()]

    rng = np.random.default_rng(0)
    idx = np.arange(len(df))
    # one set of resampled design indices, reused for every metric -- this is
    # what makes the differences paired rather than independent
    boots = [rng.choice(idx, len(idx), True) for _ in range(args.boot)]
    bdist = {}
    for lab, c in cand:
        vals = []
        for b in boots:
            s = df.iloc[b]
            if (s.binding == "true").sum() < 2 or (s.binding == "false").sum() < 2:
                vals.append(np.nan)
            else:
                vals.append(auc_df(s, c))
        bdist[c] = np.array(vals, float)

    print("\n%-28s %6s   %-16s" % ("metric", "AUC", "95% CI"))
    print("-" * 78)
    for lab, c in cand:
        a = auc_df(df, c)
        d = bdist[c][~np.isnan(bdist[c])]
        lo, hi = np.percentile(d, [2.5, 97.5])
        flag = "" if lo > 0.5 else "   (CI includes 0.5)"
        print("%-28s %6.3f   [%.3f, %.3f]%s" % (lab, a, lo, hi, flag))

    # ---- the paired question: does high accuracy beat the cheap screen? ----
    print("\n" + "=" * 78)
    print("PAIRED DIFFERENCES (same designs, same resamples)")
    print("=" * 78)
    pairs = [("high_ipTM", "cheap_ipTM", "high accuracy - cheap screen (ipTM)"),
             ("high_ipSAE_min", "cheap_ipSAE_min", "high accuracy - cheap screen (ipSAE)"),
             ("high_ipTM", "adaptyv_iptm", "high accuracy - Adaptyv's own (ipTM)"),
             ("high_ipTM", "length", "high accuracy ipTM - length alone")]
    for a_c, b_c, lab in pairs:
        if a_c not in bdist or b_c not in bdist:
            continue
        d = bdist[a_c] - bdist[b_c]
        d = d[~np.isnan(d)]
        obs = auc_df(df, a_c) - auc_df(df, b_c)
        lo, hi = np.percentile(d, [2.5, 97.5])
        # two-sided bootstrap p for "no difference"
        p = 2 * min((d <= 0).mean(), (d >= 0).mean())
        verdict = "significant" if lo > 0 or hi < 0 else "not distinguishable"
        print("  %-38s %+0.3f  [%+0.3f, %+0.3f]  p=%.3f  %s"
              % (lab, obs, lo, hi, min(p, 1.0), verdict))

    # ---- the decision gate the user set ----
    gate = 0.564  # within 10% of Adaptyv's 0.627
    best_high = max((auc_df(df, c), lab) for lab, c in cand if c.startswith("high_"))
    print("\n" + "=" * 78)
    print("DECISION GATE: best high-accuracy metric >= %.3f justifies the RL layer" % gate)
    print("=" * 78)
    a, lab = best_high
    print("  best high-accuracy metric: %s  AUC %.3f" % (lab.strip(), a))
    print("  gate %s (%.3f vs %.3f)" % ("CLEARED" if a >= gate else "NOT cleared", a, gate))
    if a >= gate:
        c = dict((l, cc) for l, cc in cand)[lab]
        d = bdist[c][~np.isnan(bdist[c])]
        print("  fraction of bootstraps above the gate: %.1f%%" % (100 * (d >= gate).mean()))
    print("\n[OK] %s" % args.out)


if __name__ == "__main__":
    main()
