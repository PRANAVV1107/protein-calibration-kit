#!/usr/bin/env python3
"""Do our metrics predict MEASURED binding?

The Adaptyv EGFR round-2 set is the only ground truth in this project: designs
that were synthesised and assayed, with 55 binders and 310 non-binders. Every
other number here is AlphaFold's opinion checked against AlphaFold's opinion.

This asks, for each metric, whether it separates real binders from real
non-binders -- reported as ROC AUC, which is threshold-free and equals the
probability that a randomly chosen binder outscores a randomly chosen
non-binder. 0.5 is coin-flipping.

AUC is computed from ranks (the Mann-Whitney identity) rather than via sklearn,
which is not installed in this environment.

Also fits a REAL calibration curve. The one shipped in this repository is fit on
labels synthesised from ipTM thresholds and the README disclaims it. These
labels are measured, so the resulting curve means something.
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
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    allv = np.concatenate([pos, neg])
    order = allv.argsort()
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks within ties
    s = pd.Series(allv)
    ranks = s.rank(method="average").to_numpy()
    r_pos = ranks[: len(pos)].sum()
    return (r_pos - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def boot_ci(pos, neg, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if len(pos) < 2 or len(neg) < 2:
        return (np.nan, np.nan)
    vals = [auc(rng.choice(pos, len(pos), True), rng.choice(neg, len(neg), True))
            for _ in range(n)]
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold-dir", default="results/queue/adaptyv_val_folds")
    ap.add_argument("--manifest", default="results/queue/adaptyv_val_manifest.csv")
    ap.add_argument("--out", default="results/adaptyv_validation.csv")
    args = ap.parse_args()

    os.chdir(KIT)
    man = pd.read_csv(args.manifest).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(args.fold_dir, "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        d = json.load(open(f))
        r = man.loc[tag]
        p = np.array(d.get("plddt", [0.0]))
        rows.append({
            "tag": tag, "binding": r["binding"], "kd": r.get("kd"),
            "length": r["length"], "composition_penalty": r["composition_penalty"],
            "adaptyv_iptm": r.get("adaptyv_iptm"),
            "ipTM": d.get("iptm"), "pTM": d.get("ptm"),
            "ipSAE_min": mn(d.get("ipsae")), "pDockQ2_min": mn(d.get("pdockq2")),
            "pLDDT": p.mean(),
        })

    if not rows:
        raise SystemExit("[FAIL] no folds matched the manifest")

    df = pd.DataFrame(rows)
    df["shaped_pDockQ2"] = df["pDockQ2_min"] - df["composition_penalty"]
    df["shaped_ipSAE"] = df["ipSAE_min"] - df["composition_penalty"]

    # pandas reads the strings 'true'/'false' back from CSV as booleans, so a
    # string comparison silently matches nothing. This is the same failure that
    # inverted every label in the Adaptyv merge; normalise explicitly.
    df["binding"] = df["binding"].map(
        lambda v: "true" if str(v).strip().lower() in ("true", "1")
        else ("false" if str(v).strip().lower() in ("false", "0") else "unknown"))
    df.to_csv(args.out, index=False)

    pos = df[df.binding == "true"]
    neg = df[df.binding == "false"]
    print("=" * 74)
    print("ADAPTYV EGFR VALIDATION -- measured outcomes")
    print("=" * 74)
    print("  folded: %d of %d labelled   binders %d   non-binders %d"
          % (len(df), len(man), len(pos), len(neg)))
    if len(pos) < 5 or len(neg) < 5:
        print("  too few of one class so far; rerun when more folds land")
        return

    print("\n%-20s %7s %7s %8s %s" % ("metric", "AUC", "", "95% CI", "binder vs non-binder mean"))
    print("-" * 74)
    metrics = ["ipTM", "ipSAE_min", "pDockQ2_min", "pTM", "pLDDT",
               "shaped_pDockQ2", "shaped_ipSAE", "adaptyv_iptm"]
    res = []
    for m in metrics:
        if m not in df.columns or df[m].isna().all():
            continue
        a = auc(pos[m], neg[m])
        lo, hi = boot_ci(pos[m], neg[m])
        res.append((m, a, lo, hi))
        star = "  <-- best" if a == max(r[1] for r in res) else ""
        print("%-20s %7.3f        [%.3f, %.3f]   %.3f vs %.3f%s"
              % (m, a, lo, hi, pos[m].mean(), neg[m].mean(), star))

    best = max(res, key=lambda r: r[1])
    print("\n  best separator: %s (AUC %.3f)" % (best[0], best[1]))
    print("  note: adaptyv_iptm is the competition's OWN reported ipTM, so it is")
    print("        a baseline for what re-scoring under our settings adds.")

    # ---- real calibration curve, on measured labels ----
    m = best[0]
    print("\n" + "=" * 74)
    print("CALIBRATION on measured outcomes -- %s" % m)
    print("=" * 74)
    sub = df[df[m].notna()].copy()
    try:
        sub["bin"] = pd.qcut(sub[m], 5, duplicates="drop")
    except ValueError:
        sub["bin"] = pd.cut(sub[m], 5)
    print("%-26s %6s %8s %s" % ("bin", "n", "binders", "observed hit rate"))
    for b, g in sub.groupby("bin", observed=True):
        k = (g.binding == "true").sum()
        print("%-26s %6d %8d   %5.1f%%" % (str(b), len(g), k, 100 * k / len(g)))
    print("\n  overall base rate: %.1f%%" % (100 * len(pos) / len(df)))
    print("\n[OK] %s" % args.out)


if __name__ == "__main__":
    main()
