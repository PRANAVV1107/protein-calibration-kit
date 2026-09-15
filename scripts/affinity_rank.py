#!/usr/bin/env python3
"""Do any of these metrics rank AFFINITY, or only presence of binding?

Every analysis in this project so far has treated the Adaptyv set as a binary
problem: binder against non-binder, scored by AUC. But all 55 binders carry a
measured KD spanning 1.2 nM to 10 uM -- almost four decades -- and that
information has gone unused.

The distinction matters for the reward function. A reward that maximises the
probability of binding at all is not the same as one that maximises affinity.
If no metric ranks KD among designs that already bind, then a policy trained
to maximise those metrics will drift toward many weak binders rather than few
strong ones, and the metric cannot tell it apart. That is a Goodhart failure
that would only show up on the bench.

Method. Restricted to the 55 measured binders, so this asks a strictly harder
question than the AUC work: given that a design binds, does the metric know
how well. Correlation is Spearman on log10(KD) -- rank-based because the
relationship need not be linear, and log because KD spans decades and an
untransformed correlation would be dominated by the weakest few.

Sign convention: lower KD is tighter binding, so a metric that predicts
affinity correlates NEGATIVELY with log10(KD). Reported as rho_vs_logKD with
that sign preserved, and interpreted in the printed verdict, rather than
silently flipped.

n = 55 is small. A Spearman rho needs roughly |rho| > 0.27 to clear p = 0.05
here, so anything below that is not distinguishable from noise and is labelled
as such rather than reported as a weak trend.
"""
import argparse
import os

import numpy as np
import pandas as pd

RAW = ("/mnt/c/Users/prana/Downloads/proteinfoldingexp/experiments/"
       "calibration_2026/raw/egfr2_result_summary.csv")


def spearman(x, y):
    """Spearman rho and a two-sided p from the t approximation."""
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    n = len(x)
    if n < 5:
        return np.nan, np.nan, n
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    rho = np.corrcoef(rx, ry)[0, 1]
    if abs(rho) >= 1.0:
        return rho, 0.0, n
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    # two-sided p via the survival function of |t|, normal approx for n>30
    from math import erfc, sqrt
    p = erfc(abs(t) / sqrt(2))
    return float(rho), float(p), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=RAW)
    ap.add_argument("--out", default="results/affinity_rank.csv")
    args = ap.parse_args()

    d = pd.read_csv(args.raw)
    d["binding_n"] = d["binding"].map(lambda v: str(v).strip().lower())
    d["kd_n"] = pd.to_numeric(d["kd"], errors="coerce")

    b = d[(d.binding_n == "true") & d.kd_n.notna()].copy()
    b["log_kd"] = np.log10(b.kd_n)
    b["length"] = b.sequence.astype(str).str.len()

    print("=" * 74)
    print("AFFINITY RANKING among measured binders")
    print("=" * 74)
    print("  binders with KD : %d" % len(b))
    print("  KD range        : %.3g to %.3g M (%.2f decades)"
          % (b.kd_n.min(), b.kd_n.max(), b.log_kd.max() - b.log_kd.min()))
    print("  median KD       : %.3g M" % b.kd_n.median())
    print()
    print("  Lower KD = tighter. A metric that predicts affinity should show a")
    print("  NEGATIVE rho against log10(KD) if higher metric = better.")
    print()

    cands = [("iptm", "adaptyv ipTM", -1),
             ("plddt", "adaptyv pLDDT", -1),
             ("pae_interaction", "adaptyv PAE_int", +1),
             ("esm_pll", "ESM pseudo-LL", -1),
             ("length", "binder length", -1)]

    print("%-20s %9s %9s %6s  %s" % ("metric", "rho", "p", "n", "verdict"))
    print("-" * 74)
    rows = []
    for col, lab, better in cands:
        if col not in b.columns:
            continue
        x = pd.to_numeric(b[col], errors="coerce").to_numpy(float)
        rho, p, n = spearman(x, b.log_kd.to_numpy(float))
        if np.isnan(rho):
            continue
        # expected sign for "predicts tighter binding":
        #   better=-1 (higher is better) -> expect rho < 0
        #   better=+1 (lower is better, e.g. PAE) -> expect rho > 0
        expected = "neg" if better < 0 else "pos"
        got = "neg" if rho < 0 else "pos"
        if p >= 0.05:
            verdict = "not distinguishable from noise"
        elif got == expected:
            verdict = "predicts affinity"
        else:
            verdict = "WRONG DIRECTION"
        rows.append({"metric": lab, "rho_vs_logKD": rho, "p": p, "n": n,
                     "verdict": verdict})
        print("%-20s %9.3f %9.3f %6d  %s" % (lab, rho, p, n, verdict))

    print()
    print("  (n=%d needs |rho| > ~0.27 for p<0.05; smaller values are noise)"
          % len(b))

    # binding_strength is a coarser label but an independent view of the same
    # thing -- if a metric ranks affinity it should separate these buckets too
    if "binding_strength" in b.columns:
        print("\n" + "=" * 74)
        print("KD by reported binding_strength (sanity check on the labels)")
        print("=" * 74)
        for s, g in b.groupby(b.binding_strength.astype(str).str.lower()):
            print("  %-8s n=%2d   median KD %.3g M   ipTM %.3f"
                  % (s, len(g), g.kd_n.median(),
                     pd.to_numeric(g.get("iptm"), errors="coerce").mean()))

    if rows:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        pd.DataFrame(rows).to_csv(args.out, index=False)
        print("\n[OK] %s" % args.out)


if __name__ == "__main__":
    main()
