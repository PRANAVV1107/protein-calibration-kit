#!/usr/bin/env python3
"""Prospective test of the radius-of-gyration predictor. REGISTERS BEFORE FOLDING.

The claim under test
--------------------
ROOT_BREADTH_RESULT.md found that radius of gyration predicts per-root quality
at rho = -0.482 (p = 0.0004, n = 50, surviving Bonferroni over six features) --
compact binders score better. That is the project's only cheap, positive,
forward-looking result.

It was also **selected post hoc from six features on the data it was measured
on**, which is exactly the situation where this project has been wrong before:
at n = 10 the leading feature was `contacts` (+0.511, "the one to test
prospectively") and at n = 50 it was dead (+0.191).

So rog is currently a hypothesis, not a finding. This tests it the only way that
counts: rank backbones that played no part in selecting it, write the prediction
down, then fold.

Design
------
Batches 1-3 hold 150 RFdiffusion backbones that have never been folded. None
contributed to the rog correlation, which came entirely from batch 0. Rank all
150 by rog, take the 10 most compact and the 10 least compact, fold 4
ProteinMPNN draws each.

  80 folds, full MSA at 1 model / 3 recycles -- the same configuration the
  n=50 correlation was measured on, so no instrument question arises.

The registration this writes contains the two groups, the predicted direction,
the effect size to beat, and the decision rule. It is written before a single
fold exists and must not be edited afterwards.
"""
import argparse
import glob
import json
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
# The campaign's ProteinMPNN output holds ONE sequence per backbone in every
# batch. The first run of this experiment silently took whatever existed and
# produced 20 folds at K=1 instead of the registered 80 at K=4 -- the same
# trap that root_breadth.sh already works around. overnight.sh now re-samples
# into results/rog_prospective/mpnn/seqs before this script runs, and this
# script REFUSES to proceed on short input rather than quietly degrading.
RESAMPLED = os.path.join(KIT, "results/rog_prospective/mpnn/seqs")
RBD = ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
       "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
       "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG")


def read_fa(path):
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
    ap.add_argument("--out", default="results/rog_prospective")
    ap.add_argument("--n-per-group", type=int, default=10)
    ap.add_argument("--draws", type=int, default=4)
    a = ap.parse_args()
    os.chdir(KIT)
    os.makedirs(a.out, exist_ok=True)

    reg_path = os.path.join(a.out, "REGISTRATION.md")
    if os.path.exists(reg_path):
        print("[ok] registration already exists, not overwriting: %s" % reg_path)

    feat = pd.read_csv("results/root_features_all.csv")
    held = feat[feat.batch.isin([1, 2, 3])].copy()
    if len(held) < 2 * a.n_per_group:
        raise SystemExit("[FAIL] only %d held-out backbones" % len(held))

    held = held.sort_values("rog")
    compact = held.head(a.n_per_group)          # low rog -> predicted BETTER
    extended = held.tail(a.n_per_group)         # high rog -> predicted WORSE

    # ---- sequences ----
    rows = []
    for grp, sub in (("compact", compact), ("extended", extended)):
        for _, r in sub.iterrows():
            # Batch-keyed: binder_0 exists in batches 1, 2 and 3 and the
            # selection can contain more than one of them.
            fa = os.path.join(RESAMPLED, "b%d_%s.fa" % (r["batch"], r["backbone"]))
            if not os.path.exists(fa):   # fall back to the campaign output
                fa = os.path.join(ROOT, "real_test_results/rbd_200_batch",
                                  "batch_%d" % r["batch"], "proteinmpnn/output/seqs",
                                  "%s.fa" % r["backbone"])
            if not os.path.exists(fa):
                print("  [warn] no sequences for batch_%d/%s" % (r["batch"], r["backbone"]))
                continue
            draws = [re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s.upper()) for s in read_fa(fa)]
            draws = [s for s in draws if s]
            for i, s in enumerate(draws[:a.draws]):
                rows.append({"tag": "b%d_%s__d%d" % (r["batch"], r["backbone"], i + 1),
                             "batch": int(r["batch"]), "backbone": r["backbone"],
                             "group": grp, "rog": float(r["rog"]),
                             "draw": i + 1, "sequence": s, "length": len(s)})
    if not rows:
        raise SystemExit("[FAIL] no sequences assembled")
    df = pd.DataFrame(rows)

    # The registered protocol is N backbones x --draws. Anything less is a
    # different experiment, and running it silently is how the first attempt
    # produced 16 backbones at K=1 against a registration specifying 20 at K=4.
    per = df.groupby(["batch", "backbone"]).size()
    want_bb, want_k = 2 * a.n_per_group, a.draws
    if len(per) < want_bb or int(per.min()) < want_k:
        msg = [
            "[FAIL] protocol not satisfiable.",
            "  have: %d backbones with %d-%d draws each" % (len(per), int(per.min()), int(per.max())),
            "  registration requires: %d backbones x %d draws" % (want_bb, want_k),
            "  re-sample first:  bash scripts/overnight.sh  (it runs ProteinMPNN)",
            "  short backbones:  %s" % ", ".join("b%d/%s" % k for k in per[per < want_k].index[:8]),
        ]
        raise SystemExit(chr(10).join(msg))

    with open(os.path.join(a.out, "to_fold.fasta"), "w") as fh:
        for _, r in df.iterrows():
            fh.write(">%s\n%s:%s\n" % (r["tag"], r["sequence"], RBD))
    df.to_csv(os.path.join(a.out, "manifest.csv"), index=False)

    # ---- the registration, written before any fold exists ----
    if not os.path.exists(reg_path):
        with open(reg_path, "w", newline="\n") as fh:
            fh.write("""# Prospective test of radius of gyration — REGISTERED %s

Written **before any fold in this experiment existed**. Do not edit below.

## The claim

Radius of gyration of the binder chain predicts per-root design quality:
**more compact backbones score higher.** Measured on batch 0, n = 50,
rho = -0.482, p = 0.0004, surviving Bonferroni over six features
(ROOT_BREADTH_RESULT.md section 4).

## Why it needs this test

`rog` was chosen by inspecting six features on the same data that measured it.
The project has already been wrong this way once: at n = 10 the leading feature
was `contacts` at +0.511, named as "the one worth testing prospectively"; at
n = 50 it was +0.191, nothing. A post-hoc feature is a hypothesis.

## Design

150 backbones in batches 1-3, none of which contributed to the correlation.
Ranked by rog; the %d most compact and %d most extended selected.
%d ProteinMPNN draws each, full MSA at 1 model / 3 recycles — the configuration
the original correlation was measured on.

    compact group   rog %.2f - %.2f   (predicted BETTER)
    extended group  rog %.2f - %.2f   (predicted WORSE)
    separation      %.2f angstroms

    compact:  %s
    extended: %s

## Predictions, fixed now

1. **Primary.** Mean best-of-%d ipTM is higher in the compact group than the
   extended group. One-sided Mann-Whitney, alpha = 0.05.
2. **Effect size to beat.** The batch-0 correlation implies a gap of roughly
   **0.15-0.25 ipTM** between these groups. A significant result with a gap
   below 0.10 confirms direction but means the feature is too weak to screen on.
3. **Secondary.** Spearman(rog, best-of-%d) across all %d backbones should be
   negative and of similar magnitude to -0.482.

## Outcomes, decided now

| result | conclusion |
|---|---|
| compact > extended, p < 0.05, gap >= 0.10 | rog is a validated free pre-screen; keep recommendation #2 |
| compact > extended, p < 0.05, gap < 0.10 | direction real, effect too small to act on; demote to a note |
| not significant | rog joins `contacts` as a post-hoc feature that did not replicate; **remove recommendation #2 from the paper** |
| extended > compact | the batch-0 correlation was an artifact; report it as such |

**Declared expectation:** I expect this to replicate but weaker than -0.482,
because the original estimate was selected as the largest of six and is
therefore biased upward by selection. A gap of 0.10-0.15 seems most likely.

## Limits, acknowledged in advance

- One target (RBD), one generator configuration.
- 20 backbones, 4 draws each. Powered only for a large effect.
- Extreme-group design inflates the apparent effect relative to a full-range
  correlation; the secondary test across all 20 is the fairer estimate.
""" % (datetime.now().strftime("%Y-%m-%d %H:%M"),
       a.n_per_group, a.n_per_group, a.draws,
       compact.rog.min(), compact.rog.max(),
       extended.rog.min(), extended.rog.max(),
       extended.rog.min() - compact.rog.max(),
       ", ".join("b%d/%s" % (r.batch, r.backbone) for r in compact.itertuples()),
       ", ".join("b%d/%s" % (r.batch, r.backbone) for r in extended.itertuples()),
       a.draws, a.draws, 2 * a.n_per_group))

    print("=" * 70)
    print("ROG PROSPECTIVE TEST -- registered and prepared")
    print("=" * 70)
    print("  held-out pool: %d backbones (batches 1-3)" % len(held))
    print("  compact  rog %.2f - %.2f  (predicted better)" % (compact.rog.min(), compact.rog.max()))
    print("  extended rog %.2f - %.2f  (predicted worse)" % (extended.rog.min(), extended.rog.max()))
    print("  separation %.2f A" % (extended.rog.min() - compact.rog.max()))
    print("  %d backbones x %d draws = %d folds" % (df.backbone.nunique(), a.draws, len(df)))
    print("\n[OK] %s" % reg_path)
    print("[OK] %s/to_fold.fasta" % a.out)


if __name__ == "__main__":
    main()
