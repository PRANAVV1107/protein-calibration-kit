#!/usr/bin/env python3
"""Choose the MD panel for correlating predicted scores against stability.

Two selection rules, both of which exist to stop the correlation being
meaningless before it is measured.

Spread, not quality. The question is whether a metric PREDICTS stability, so
the panel must span the metric's range. Picking the top-scoring complexes --
the obvious instinct, and what every other selection in this project does --
would truncate the predictor and drive the correlation toward zero no matter
how good the metric is. Complexes are therefore drawn evenly across score
quantiles, deliberately including ones expected to fall apart.

Three targets, not one. The previously armed panel was fifteen il7ra
refinements, which would have produced a single-target result that generalises
to nothing. The variance study folded 60 complexes per target at high accuracy
across il7ra, PD-L1 and RBD, so the panel is drawn from all three and target is
recorded for a per-target breakdown.

Sample size. At n=15 a Spearman rho needs |rho| > 0.51 for p < 0.05, so a
ranking of five metrics against 15 points would be reading noise. At n=40 the
threshold falls to ~0.31, which is attainable. That is why the panel is 40 at
2 ns rather than 15 at 5 ns for the same wall time -- for a correlation, spread
across designs buys more than trajectory length per design.
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
VS = os.path.join(KIT, "results/variance_study")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=40, help="total complexes")
    ap.add_argument("--metric", default="ipTM",
                    help="metric whose range the panel should span")
    ap.add_argument("--out", default="results/md_panel.csv")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.chdir(KIT)
    res = pd.read_csv(os.path.join(VS, "variance_results.csv"))
    rng = np.random.default_rng(args.seed)

    # map (target, backbone, sample) -> fold PDB
    rows = []
    for t in res.target.unique():
        for f in glob.glob(os.path.join(VS, t, "folds", "*_unrelaxed_rank_001_*.pdb")):
            stem = os.path.basename(f).split("_unrelaxed_")[0]
            rows.append({"target": t, "stem": stem, "pdb": f})
    pdbs = pd.DataFrame(rows)
    if pdbs.empty:
        raise SystemExit("[FAIL] no variance-study fold PDBs found")

    res["stem"] = res.backbone.astype(str) + "__" + res["sample"].astype(str)
    df = res.merge(pdbs, on=["target", "stem"], how="inner")
    if df.empty:
        raise SystemExit("[FAIL] scores and PDBs did not match on (target, stem)")

    print("=" * 74)
    print("MD PANEL SELECTION")
    print("=" * 74)
    print("  scored complexes with structures: %d across %d targets"
          % (len(df), df.target.nunique()))

    per = max(1, args.n // df.target.nunique())
    picks = []
    for t, g in df.groupby("target"):
        g = g[g[args.metric].notna()].copy()
        if len(g) <= per:
            picks.append(g)
            continue
        # even coverage of the metric's range: one draw per quantile bin
        g["_q"] = pd.qcut(g[args.metric], per, labels=False, duplicates="drop")
        sel = (g.groupby("_q", group_keys=False)
                 .apply(lambda x: x.sample(1, random_state=int(rng.integers(1e6)))))
        picks.append(sel)
    out = pd.concat(picks, ignore_index=True).drop(columns=["_q"], errors="ignore")

    print()
    print("%-8s %6s %8s %8s %8s %8s" % ("target", "n", "min", "median", "max", "range"))
    print("-" * 74)
    for t, g in out.groupby("target"):
        m = g[args.metric]
        print("%-8s %6d %8.3f %8.3f %8.3f %8.3f"
              % (t, len(g), m.min(), m.median(), m.max(), m.max() - m.min()))
    m = out[args.metric]
    print("-" * 74)
    print("%-8s %6d %8.3f %8.3f %8.3f %8.3f"
          % ("TOTAL", len(out), m.min(), m.median(), m.max(), m.max() - m.min()))
    print("\n  panel spans %.0f%% of the full %s range in the study"
          % (100 * (m.max() - m.min()) / (df[args.metric].max() - df[args.metric].min()),
             args.metric))

    out.to_csv(args.out, index=False)
    with open("results/md_panel_paths.txt", "w") as f:
        for p in out.pdb:
            f.write(p + "\n")
    print("\n[OK] %d complexes -> %s" % (len(out), args.out))
    print("     paths -> results/md_panel_paths.txt")


if __name__ == "__main__":
    main()
