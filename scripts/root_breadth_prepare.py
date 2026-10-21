#!/usr/bin/env python3
"""Fold set for the root-breadth sweep: 50 RBD backbones, 4 sequence draws each.

The prediction being tested
---------------------------
ALLOCATION_RESULT.md found that no policy beats uniform and that an oracle
spending every fold on the best-mean root loses to spreading folds across ten.
The tail-statistic argument in the paper explains why: the objective is a
maximum, and maxima come from sampling many distributions rather than the best
one harder.

That argument makes a sharp prediction nobody has tested: **E[best] should keep
climbing with the NUMBER of roots, faster than it climbs with draws per root.**
With ten roots it cannot be tested -- uniform flattens at 0.869 simply because
there are only ten distributions to draw from. Fifty roots can test it.

Why this configuration
----------------------
Full MSA at 1 model / 3 recycles: the configuration section 3.2 measured at AUC
0.627 against 0.614 for 3 models / 8 recycles, at a quarter the cost. It is what
the paper recommends, so the sweep should use it.

**Ten of the fifty backbones are already measured in the variance study** at
3 models / 8 recycles (binder_0, 5, 9, 10, 13, 14, 15, 34, 42, 47). They are
kept in deliberately: they calibrate this instrument against that one directly,
which is exactly the check whose absence invalidated three search runs. Never
run a new configuration without an overlap set.
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
# The campaign's own ProteinMPNN output holds only ONE draw per backbone, which
# is not enough to estimate a per-root distribution. root_breadth.sh re-samples
# 4 draws per backbone into results/root_breadth/mpnn/seqs first.
SEQS_DEFAULT = os.path.join(KIT, "results/root_breadth/mpnn/seqs")
RBD = ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
       "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
       "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG")
OVERLAP = {"binder_0", "binder_5", "binder_9", "binder_10", "binder_13",
           "binder_14", "binder_15", "binder_34", "binder_42", "binder_47"}


def read_fa(path):
    """ProteinMPNN .fa: first record is the native/input, the rest are draws."""
    seqs, cur = [], None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cur is not None:
                seqs.append(cur)
            cur = ""
        elif cur is not None:
            cur += line
    if cur:
        seqs.append(cur)
    return seqs[1:] if len(seqs) > 1 else seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/root_breadth")
    ap.add_argument("--draws", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--seqs-dir", default=None,
                    help="ProteinMPNN seqs/ directory (default: the re-sampled one)")
    a = ap.parse_args()
    os.chdir(KIT)
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    seqs_dir = a.seqs_dir or SEQS_DEFAULT
    if not glob.glob(os.path.join(seqs_dir, "*.fa")):
        raise SystemExit("[FAIL] no .fa in %s -- run the MPNN step first" % seqs_dir)
    rows = []
    short = []
    for fa in sorted(glob.glob(os.path.join(seqs_dir, "*.fa")),
                     key=lambda p: int(re.search(r"binder_(\d+)", p).group(1))):
        bb = os.path.basename(fa)[:-3]
        draws = [re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s.upper()) for s in read_fa(fa)]
        draws = [s for s in draws if s]
        if not draws:
            continue
        if len(draws) < a.draws:
            short.append((bb, len(draws)))
        # take the first `draws` rather than sampling: MPNN writes them in
        # generation order, and picking by score would bias the pool upward.
        for i, s in enumerate(draws[:a.draws]):
            rows.append({"tag": "%s__d%d" % (bb, i + 1), "backbone": bb,
                         "draw": i + 1, "sequence": s, "length": len(s),
                         "overlap": bb in OVERLAP})

    if not rows:
        raise SystemExit("[FAIL] no sequences found in %s" % seqs_dir)
    df = pd.DataFrame(rows)

    fa_out = os.path.join(a.out, "to_fold.fasta")
    with open(fa_out, "w") as fh:
        for _, r in df.iterrows():
            fh.write(">%s\n%s:%s\n" % (r["tag"], r["sequence"], RBD))
    df.to_csv(os.path.join(a.out, "manifest.csv"), index=False)

    print("=" * 66)
    print("ROOT-BREADTH SWEEP -- fold set")
    print("=" * 66)
    print("  %d backbones, %d draws each -> %d folds"
          % (df.backbone.nunique(), a.draws, len(df)))
    print("  overlap with variance study: %d backbones (instrument calibration)"
          % df[df.overlap].backbone.nunique())
    print("  binder length %d-%d aa (median %d)"
          % (df.length.min(), df.length.max(), int(df.length.median())))
    if short:
        print("  [warn] %d backbones had fewer than %d draws: %s"
              % (len(short), a.draws, ", ".join("%s(%d)" % s for s in short[:6])))
    print("\n[OK] %s" % fa_out)
    print("[OK] %s" % os.path.join(a.out, "manifest.csv"))


if __name__ == "__main__":
    main()
