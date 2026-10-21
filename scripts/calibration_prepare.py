#!/usr/bin/env python3
"""Proper instrument calibration: the SAME sequences under both configurations.

Why the existing overlap set failed
-----------------------------------
The root-breadth sweep included 10 backbones that the variance study had also
folded, intended as a calibration. It cannot serve as one: the two datasets
differ in configuration (1 model/3 recycles versus 3 models/8 recycles), in
sequence set (only 4 of 40 position-matched draws are identical, because
variance_prepare selects a seeded-random subset and folds three targets in one
MPNN call, so the RNG stream differs), and in K (4 versus 6). Three variables,
one measurement.

The observed offset is -0.447 on max and -0.227 on means, with rank correlation
+0.148. Part is explained -- reporting `rank_001` of `--num-models 3` is a
max-over-3 and inflates by a measured +0.132 -- but the rest cannot be
attributed.

What this does
--------------
Re-folds the variance study's OWN 60 RBD sequences at 1 model / 3 recycles.
Same designs, same K, same everything except the configuration. That is what a
calibration means, and it is the only way to license pooling the breadth sweep
with the rest of the project or comparing it to the 0.870 bar.
"""
import argparse
import os

import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
RBD = ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
       "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
       "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/calibration")
    a = ap.parse_args()
    os.chdir(KIT)
    os.makedirs(a.out, exist_ok=True)

    vs = pd.read_csv("results/variance_study/rbd/manifest.csv")
    vs = vs.dropna(subset=["sequence"]).copy()
    vs["tag"] = ["cal__%s__%s" % (r.backbone, str(getattr(r, "sample")))
                 for r in vs.itertuples()]

    with open(os.path.join(a.out, "to_fold.fasta"), "w") as fh:
        for _, r in vs.iterrows():
            fh.write(">%s\n%s:%s\n" % (r["tag"], r["sequence"], RBD))
    vs.to_csv(os.path.join(a.out, "manifest.csv"), index=False)

    print("=" * 70)
    print("INSTRUMENT CALIBRATION -- fold set")
    print("=" * 70)
    print("  %d sequences over %d backbones" % (len(vs), vs.backbone.nunique()))
    print("  identical designs to results/variance_study/rbd, which was folded")
    print("  at 3 models / 8 recycles. This folds them at 1 model / 3 recycles.")
    print("  Only the configuration differs.")
    print("\n[OK] %s/to_fold.fasta" % a.out)


if __name__ == "__main__":
    main()
