#!/usr/bin/env python3
"""What IS a root, and does anything cheap about it predict its quality?

Context
-------
Root choice turned out to matter more than search: between-root spread of
E[best-of-6] ran 0.371 to 0.827 on the high-accuracy instrument, while three
searches from one root gained little. Evaluating a root costs K folds. If some
property of the BACKBONE FILE -- computable in milliseconds, before any folding
-- predicts root quality, then roots can be pre-screened for free and the whole
budget spent on the promising ones.

This measures the backbones the project actually has and tests that.

Features, all from the PDB alone
--------------------------------
  length        binder chain residue count. RFdiffusion sampled 70-100 here, so
                this varies by design rather than by accident.
  rog           radius of gyration of the binder CA trace; compactness.
  contacts      binder CA within 10 A of any target CA. Interface size proxy.
  hotspot_c     contacts against the six hotspots the generation was steered to
                (E449, E493, E496, E498, E501, E505). Whether the backbone
                actually engaged what it was asked to engage.
  helix_frac    fraction of residues whose CA(i)-CA(i+4) distance is 5-7 A,
                the helical signature. Crude, but it separates helical bundles
                from loops without a DSSP dependency.
  contact_dens  contacts / length. Interface per residue spent.

Power warning, stated up front
------------------------------
Only ten RBD backbones have high-accuracy quality measured (variance study,
6 draws each). Six features against ten points will produce a strong-looking
correlation by chance alone -- with n=10, |rho| must exceed about 0.63 for
p<0.05 uncorrected, and testing six features makes one such hit likely even if
nothing is real. Treat every number here as hypothesis-generating. The honest
use is to pick ONE feature to test prospectively on new backbones, not to
select the best-looking one and report its p-value.
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
HOTSPOTS = {449, 493, 496, 498, 501, 505}


def read_ca(path):
    """CA coordinates per chain, plus target residue numbers."""
    chains = {}
    resnum = {}
    for line in open(path):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        ch = line[21]
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        chains.setdefault(ch, []).append(xyz)
        resnum.setdefault(ch, []).append(int(line[22:26]))
    return {k: np.asarray(v) for k, v in chains.items()}, resnum


def features(path):
    ch, resnum = read_ca(path)
    if len(ch) < 2:
        return None
    # binder is the designed chain; RFdiffusion writes it as A, target as B
    bkey = "A" if "A" in ch else sorted(ch)[0]
    tkey = [k for k in ch if k != bkey][0]
    b, t = ch[bkey], ch[tkey]
    if len(b) < 10 or len(t) < 10:
        return None

    d = np.linalg.norm(b[:, None, :] - t[None, :, :], axis=-1)
    contacts = int((d < 10.0).sum())

    tnum = np.asarray(resnum[tkey])
    hs = np.isin(tnum, sorted(HOTSPOTS))
    hotspot_c = int((d[:, hs] < 10.0).sum()) if hs.any() else 0

    cen = b - b.mean(axis=0)
    rog = float(np.sqrt((cen ** 2).sum(axis=1).mean()))

    if len(b) > 4:
        d14 = np.linalg.norm(b[4:] - b[:-4], axis=1)
        helix = float(((d14 > 5.0) & (d14 < 7.0)).mean())
    else:
        helix = np.nan

    return {"length": len(b), "rog": rog, "contacts": contacts,
            "hotspot_c": hotspot_c, "helix_frac": helix,
            "contact_dens": contacts / len(b)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batches", default="0,1,2,3")
    a = ap.parse_args()

    rows = []
    for bt in a.batches.split(","):
        d = os.path.join(ROOT, "real_test_results/rbd_200_batch",
                         "batch_%s" % bt, "rfdiffusion")
        for p in sorted(glob.glob(os.path.join(d, "binder_*.pdb"))):
            f = features(p)
            if f:
                f["batch"] = bt
                f["backbone"] = os.path.basename(p)[:-4]
                rows.append(f)
    if not rows:
        raise SystemExit("[FAIL] no backbones parsed")
    df = pd.DataFrame(rows)

    print("=" * 74)
    print("THE ROOT SYSTEM: what RFdiffusion actually produced")
    print("=" * 74)
    print("  %d backbones across %d batches of 50, target RBD (6M0J chain E)\n"
          % (len(df), df.batch.nunique()))
    print("  %-14s %8s %8s %8s %8s" % ("feature", "min", "median", "max", "sd"))
    print("  " + "-" * 50)
    for c in ("length", "rog", "contacts", "hotspot_c", "helix_frac", "contact_dens"):
        print("  %-14s %8.2f %8.2f %8.2f %8.2f"
              % (c, df[c].min(), df[c].median(), df[c].max(), df[c].std()))

    print("\n  Binder length was SAMPLED, not fixed -- the generation contig")
    print("  allowed a range, so roots differ in size before anything else.")
    print("  length histogram:")
    h = df.length.value_counts().sort_index()
    for L in range(int(df.length.min()), int(df.length.max()) + 1, 5):
        n = int(h[(h.index >= L) & (h.index < L + 5)].sum())
        print("    %3d-%3d  %s (%d)" % (L, L + 4, "#" * n, n))

    # ---- join to measured quality where it exists ----
    v = pd.read_csv(os.path.join(KIT, "results/variance_study/variance_results.csv"))
    v = v[v.target == "rbd"]
    q = v.groupby("backbone").ipTM.agg(["max", "mean"]).rename(
        columns={"max": "best_of_6", "mean": "mean_of_6"})
    m = df[df.batch == "0"].set_index("backbone").join(q, how="inner")

    print("\n" + "=" * 74)
    print("DOES ANY CHEAP FEATURE PREDICT ROOT QUALITY?")
    print("=" * 74)
    print("  %d backbones with high-accuracy quality (6 draws each)\n" % len(m))
    if len(m) < 5:
        print("  too few to correlate")
        return
    print("  %-14s %9s %9s   %s" % ("feature", "rho", "p", "verdict"))
    print("  " + "-" * 56)
    crit = 0.63   # |rho| for p<0.05 at n=10, uncorrected
    for c in ("length", "rog", "contacts", "hotspot_c", "helix_frac", "contact_dens"):
        rho, p = stats.spearmanr(m[c], m.best_of_6)
        flag = "nothing"
        if abs(rho) >= crit:
            flag = "LOOKS STRONG -- but see power warning"
        elif abs(rho) >= 0.4:
            flag = "weak hint"
        print("  %-14s %9.3f %9.3f   %s" % (c, rho, p, flag))

    print("""
  With n=%d and six features tested, one apparently-significant hit is expected
  by chance. Nothing here is a result. The correct next step is to pick the
  single most mechanistically plausible feature and test it PROSPECTIVELY on
  backbones that played no part in choosing it.""" % len(m))

    m.sort_values("best_of_6", ascending=False).to_csv(
        os.path.join(KIT, "results/root_anatomy.csv"))
    df.to_csv(os.path.join(KIT, "results/root_features_all.csv"), index=False)
    print("\n[OK] results/root_anatomy.csv  (the %d with quality)" % len(m))
    print("[OK] results/root_features_all.csv  (all %d backbones)" % len(df))


if __name__ == "__main__":
    main()
