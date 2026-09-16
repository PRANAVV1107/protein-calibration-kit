#!/usr/bin/env python3
"""Do AlphaFold's confidence metrics predict physical interface stability?

Why this is worth asking separately from everything else
--------------------------------------------------------
ipTM, ipSAE, pDockQ2 and pLDDT are one network's opinion re-expressed four
ways; they cannot disagree with each other informatively. Molecular dynamics
asks a physically different question -- propagate the predicted complex under a
force field and see whether the interface survives -- so it is the only
orthogonal evidence available here without a wet lab.

What it is NOT: a binding assay. A complex can be stable in simulation and not
bind, or bind and drift in 2 ns. Stability correlation is a secondary axis, and
the measured Adaptyv outcomes remain the stronger ground truth. Do not let a
strong correlation here be read as validation against binding.

The panel
---------
37 complexes drawn across il7ra, PD-L1 and RBD, deliberately spanning the full
ipTM range (0.13-0.92) rather than selecting good ones -- picking only
high-scoring complexes would truncate the predictor and drive any correlation
toward zero however good the metric was. 2 ns each, implicit solvent GBn2,
Amber14, 300 K.

Stability measures, and why more than one
-----------------------------------------
  contacts_retained   final interface contacts / initial. Near 1 or above means
                      the interface held; decay toward 0 means it came apart.
  binder_rmsd_final   binder displacement after superposing on the TARGET, so
                      it reports drift relative to the partner rather than
                      overall tumbling.
  com_drift           change in centre-of-mass separation. Climbing means
                      dissociation.
They fail differently: an interface can lose contacts while the binder stays
put (local unfolding) or keep contacts while sliding (rigid-body drift).
Agreement between them is informative; disagreement is a flag.

Correlations are Spearman. Reported pooled AND within-target, because with
three targets of very different difficulty a pooled correlation can be driven
entirely by between-target differences -- the same trap that made
ProteinMPNN's likelihood look like a usable value function (pooled rho -0.52,
within-backbone -0.06).

At n=37 a Spearman rho needs |rho| > ~0.32 for p<0.05. Within a single target
(n about 12) it needs |rho| > ~0.58, so within-target figures are reported for
direction only and should not be called significant.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"


def spearman(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    n = len(x)
    if n < 5:
        return np.nan, np.nan, n
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return np.nan, np.nan, n
    rho = float(np.corrcoef(rx, ry)[0, 1])
    if abs(rho) >= 1:
        return rho, 0.0, n
    from math import erfc, sqrt
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    return rho, float(erfc(abs(t) / sqrt(2))), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md-dir", default="results/md_gpu")
    ap.add_argument("--panel", default="results/md_panel.csv")
    ap.add_argument("--out", default="results/md_correlation.csv")
    args = ap.parse_args()
    os.chdir(KIT)

    rows = []
    for f in glob.glob(os.path.join(args.md_dir, "*_md.json")):
        d = json.load(open(f))
        if "error" in d:
            continue
        rows.append({
            "name": d["name"],
            "contacts_retained": d.get("contacts_retained"),
            "binder_rmsd": d.get("binder_rmsd_final_nm"),
            "com_drift": (d.get("com_final_nm", np.nan) - d.get("com_start_nm", np.nan)),
            "dissociated": bool(d.get("dissociated", False)),
            "ns": d.get("ns"),
        })
    md = pd.DataFrame(rows)
    if md.empty:
        raise SystemExit("[FAIL] no MD results found")

    panel = pd.read_csv(args.panel)
    panel["name"] = panel.pdb.apply(
        lambda p: os.path.basename(p).replace(".pdb", ""))
    df = panel.merge(md, on="name", how="inner")
    if df.empty:
        raise SystemExit("[FAIL] MD results did not join to the panel")

    print("=" * 78)
    print("PREDICTED CONFIDENCE vs SIMULATED STABILITY")
    print("=" * 78)
    print("  complexes %d over %d targets, %.0f ns each"
          % (len(df), df.target.nunique(), df.ns.median()))
    print("  dissociated in simulation: %d" % int(df.dissociated.sum()))
    print("  ipTM range %.2f-%.2f (panel spans the full study range by design)"
          % (df.ipTM.min(), df.ipTM.max()))
    print("\n  contacts_retained: higher = interface held")
    print("  binder_rmsd, com_drift: LOWER = more stable, so a good predictor")
    print("  correlates NEGATIVELY with these two.\n")

    metrics = [c for c in ("ipTM", "pTM", "mean_pLDDT", "penalty", "length")
               if c in df.columns]
    stabs = ["contacts_retained", "binder_rmsd", "com_drift"]
    out = []
    for s in stabs:
        print("-" * 78)
        print("%s" % s)
        for m in metrics:
            rho, p, n = spearman(df[m].to_numpy(float), df[s].to_numpy(float))
            if np.isnan(rho):
                continue
            sig = "significant" if p < 0.05 else "noise"
            print("  %-12s rho %+.3f  p %.3f  n %2d   %s" % (m, rho, p, n, sig))
            out.append({"stability": s, "metric": m, "rho": rho, "p": p, "n": n,
                        "scope": "pooled"})
        for t, g in df.groupby("target"):
            r_i, p_i, n_i = spearman(g["ipTM"].to_numpy(float),
                                     g[s].to_numpy(float))
            if not np.isnan(r_i):
                print("    within %-6s ipTM rho %+.3f  n %d  (direction only, "
                      "underpowered)" % (t, r_i, n_i))
                out.append({"stability": s, "metric": "ipTM", "rho": r_i,
                            "p": p_i, "n": n_i, "scope": "within_%s" % t})

    print("\n" + "=" * 78)
    print("RANKING: which metric best predicts stability (pooled |rho|)")
    print("=" * 78)
    o = pd.DataFrame(out)
    pooled = o[o.scope == "pooled"].copy()
    pooled["absrho"] = pooled.rho.abs()
    rank = (pooled.groupby("metric").absrho.mean()
            .sort_values(ascending=False))
    for m, v in rank.items():
        print("  %-12s mean |rho| %.3f across the three stability measures" % (m, v))
    print("\n  n=%d needs |rho| > ~0.32 for p<0.05. Anything below that is not"
          % len(df))
    print("  distinguishable from noise however it ranks.")
    print("\n  Stability is NOT binding. The Adaptyv measured outcomes remain the")
    print("  stronger ground truth; this is an orthogonal secondary axis.")

    o.to_csv(args.out, index=False)
    print("\n[OK] %s" % args.out)


if __name__ == "__main__":
    main()
