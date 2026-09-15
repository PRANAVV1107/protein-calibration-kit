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


def auc_stratified(df, col, by="length", k=3):
    """AUC computed within length strata, then pooled.

    The reason this is the headline number and not plain AUC: binders in this
    set are longer than non-binders (127 aa against 100), so any metric that
    tracks length scores above 0.5 without modelling anything. Comparing only
    within strata of similar length removes that channel -- a metric that
    still separates binders inside a narrow length band is separating them on
    something other than size.

    Strata are weighted by n_pos * n_neg, the number of comparable pairs each
    contributes, which is the Mantel-Haenszel weighting for this statistic.
    """
    d = df[df[col].notna() & df[by].notna()].copy()
    if len(d) < 3 * k:
        return np.nan
    try:
        d["_s"] = pd.qcut(d[by], k, labels=False, duplicates="drop")
    except ValueError:
        return np.nan
    num = den = 0.0
    for _, g in d.groupby("_s"):
        p = g.loc[g.binding == "true", col]
        n = g.loc[g.binding == "false", col]
        if len(p) < 1 or len(n) < 1:
            continue
        a = auc(p, n)
        if np.isnan(a):
            continue
        w = len(p) * len(n)
        num += w * a
        den += w
    return num / den if den else np.nan


def residualize(df, col, by="length"):
    """The metric with its linear dependence on length removed.

    A second, independent check on the same confound. Stratifying and
    residualizing fail in different ways -- stratification is robust to a
    nonlinear length effect but throws away resolution, regression keeps
    resolution but assumes linearity -- so agreement between them is
    informative and disagreement is a flag to investigate.
    """
    d = df[[col, by]].dropna()
    if len(d) < 5 or d[by].std() == 0:
        return pd.Series(np.nan, index=df.index)
    b, a = np.polyfit(d[by].to_numpy(float), d[col].to_numpy(float), 1)
    return df[col] - (a + b * df[by])


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
    if len(high) < len(hman):
        # ColabFold sorts its input queries by length and folds the shortest
        # first, so a partial arm is not a random subsample -- it is the short
        # tail. At 70/165 the folded set spanned 27-58 aa against a manifest
        # spanning 27-250. Every AUC below is then computed inside a narrow
        # length band and does not transfer to the full set. This prints the
        # actual coverage so no interim file can be mistaken for a result.
        cov_lo, cov_hi = df.length.min(), df.length.max()
        all_lo, all_hi = hman.length.min(), hman.length.max()
        frac = (df.length.max() - df.length.min()) / max(all_hi - all_lo, 1)
        print("\n  *** INTERIM -- NOT A RESULT ***")
        print("  ColabFold folds shortest-first, so this is the short tail,")
        print("  not a random subsample.")
        print("  length covered here : %d-%d aa (sd %.0f)"
              % (cov_lo, cov_hi, df.length.std()))
        print("  length in full set  : %d-%d aa" % (all_lo, all_hi))
        print("  covering %.0f%% of the length range; %d designs still unfolded."
              % (100 * frac, len(hman) - len(high)))
        if frac < 0.5:
            print("  Length-adjusted columns are unreliable at this coverage:")
            print("  there is too little length variation left to adjust for.")

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

    # ---- the confound check that decides whether any of the above is real ----
    print("\n" + "=" * 78)
    print("LENGTH-ADJUSTED (binders %.0f aa vs non-binders %.0f aa in this set)"
          % (df.loc[df.binding == "true", "length"].mean(),
             df.loc[df.binding == "false", "length"].mean()))
    print("=" * 78)
    print("  raw    = plain AUC, inherits whatever length predicts on its own")
    print("  strat  = pooled within length tertiles")
    print("  resid  = AUC after regressing length out of the metric")
    print("  A metric is only real if strat and resid both stay above 0.5.\n")
    print("%-28s %7s %7s %7s" % ("metric", "raw", "strat", "resid"))
    print("-" * 78)
    adj = {}
    for lab, c in cand:
        if c == "length":
            continue
        raw = auc_df(df, c)
        st = auc_stratified(df, c)
        tmp = df.copy()
        tmp["_r"] = residualize(df, c)
        rs = auc_df(tmp, "_r")
        adj[c] = (raw, st, rs)
        drop = ""
        if not np.isnan(st) and not np.isnan(rs):
            if max(st, rs) < 0.55:
                drop = "   <- collapses; was length"
            elif min(st, rs) > 0.60:
                drop = "   <- survives"
        print("%-28s %7.3f %7.3f %7.3f%s" % (lab, raw, st, rs, drop))

    # ---- the decision gate ----
    #
    # The bar is "within 10% of Adaptyv", and which Adaptyv number that means
    # matters more than it sounds. Their AUC on the full 365-design set is
    # 0.627, giving a bar of 0.564. On the subset folded here they score
    # higher, so the bar derived from THAT is stricter. Only the second is
    # like-for-like: same designs, same labels, same class balance, which is
    # the entire reason this script restricts both arms to a common set.
    # Judging our subset number against their full-set number re-introduces
    # exactly the confound the script exists to remove.
    #
    # So the bar is computed from the measured same-design baseline rather
    # than hardcoded, and both sides are taken length-controlled -- correcting
    # ours while crediting theirs uncorrected would flatter us.
    MARGIN = 0.90
    gate = 0.564  # legacy: 0.9 * 0.627, the full-365 figure, kept for continuity
    print("\n" + "=" * 78)
    print("DECISION GATE")
    print("=" * 78)
    highs = [(lab, c) for lab, c in cand if c.startswith("high_")]
    best_raw = max(highs, key=lambda lc: auc_df(df, lc[1]))
    # judged on the weaker of the two adjusted figures: a metric has to clear
    # the gate on the confound-controlled number, not the flattering one
    def adjmin(c):
        _, st, rs = adj[c]
        v = [x for x in (st, rs) if not np.isnan(x)]
        return min(v) if v else np.nan
    best_adj = max(highs, key=lambda lc: (adjmin(lc[1]) if not np.isnan(adjmin(lc[1])) else -1))
    a_raw = auc_df(df, best_raw[1])
    a_adj = adjmin(best_adj[1])
    print("  best on raw AUC       : %-18s %.3f" % (best_raw[0].strip(), a_raw))
    print("  best length-controlled: %-18s %.3f  (weaker of strat/resid)"
          % (best_adj[0].strip(), a_adj))
    base_adj = adjmin("adaptyv_iptm") if "adaptyv_iptm" in adj else np.nan
    if not np.isnan(base_adj):
        bar = MARGIN * base_adj
        print("  Adaptyv, SAME designs, adj: %-18s %.3f  -> bar %.3f"
              % ("adaptyv_iptm", base_adj, bar))
        ok = (not np.isnan(a_adj)) and a_adj >= bar
        print("\n  LIKE-FOR-LIKE GATE %s  (%.3f vs %.3f, %+.3f)"
              % ("CLEARED" if ok else "NOT CLEARED", a_adj, bar, a_adj - bar))
    else:
        bar, ok = np.nan, False
        print("\n  Adaptyv baseline missing on this subset -- gate cannot be judged.")
    legacy_ok = (not np.isnan(a_adj)) and a_adj >= gate
    print("  legacy gate vs the full-365 figure (%.3f): %s"
          % (gate, "cleared" if legacy_ok else "not cleared"))
    if legacy_ok and not ok and not np.isnan(bar):
        print("\n  These disagree, and the like-for-like gate is the one that counts.")
        print("  The legacy gate judges our subset number against their full-set")
        print("  number, which is not the same comparison.")
    if a_raw >= (bar if not np.isnan(bar) else gate) and not ok:
        print("  Raw AUC would have cleared it; the gap is the length confound.")
        print("  Do not report the raw figure as the result.")
    d = bdist[best_adj[1]][~np.isnan(bdist[best_adj[1]])]
    print("  (raw-AUC bootstraps above the like-for-like bar: %.1f%%)"
          % (100 * (d >= (bar if not np.isnan(bar) else gate)).mean()))
    print("\n[OK] %s" % args.out)


if __name__ == "__main__":
    main()
