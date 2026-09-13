#!/usr/bin/env python3
"""Build a null panel: every designed binder against targets it was NOT designed for.

Why this is needed. The calibration curve in this repository is currently fit on
labels synthesised from ipTM thresholds -- circular, and useless for deciding
whether a score is good. What it actually needs is an empirical null: the
distribution each metric takes on pairs that should NOT bind.

Without that, "ipSAE_min 0.717 is high" is unanchored. High relative to what?
The field's 0.61 cutoff comes from other people's datasets, with different chain
lengths, targets and design methods -- and ipTM in particular is sensitive to
chain length, which is the whole reason ipSAE exists.

What a null panel gives you:
  * a percentile for any candidate against ~N non-binding pairs from YOUR pipeline
  * negatives labelled BY CONSTRUCTION rather than derived from the metric being
    calibrated, which removes the circularity
  * a comparison of which metric (ipTM / ipSAE_min / pDockQ2_min) separates
    designed-for from not-designed-for pairs most sharply

Construction. Three targets here have binders designed against them and MSAs
already on disk. Every cross-pairing -- an IL-7Ra binder against RBD, say -- is a
negative by design intent: unrelated proteins do not bind by default.

Honest limitation: these are PREDICTED non-binders. No assay has been run. The
prior is strong but not certainty, and any genuine cross-reactivity would sit in
this panel mislabelled -- which is itself worth discovering.

Cost. Folded through the hybrid-MSA route (~50 s each) rather than high accuracy
(~4.5 min), because characterising a null needs volume, not precision.
"""
import argparse
import glob
import os
import random
import sys

import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
sys.path.insert(0, os.path.join(KIT, "scripts"))
from make_hybrid_msa import build_hybrid, extract_target_msa, read_a3m, split_a3m_row  # noqa: E402

# a source complex a3m per target -- donates that target's alignment
SOURCE_A3M = {
    "il7ra": "results/variance_study/il7ra/folds/*.a3m",
    "rbd": "results/variance_study/rbd/folds/*.a3m",
    "pdl1": "results/variance_study/pdl1/folds/*.a3m",
}

# where each target's designed binders live (multi-draw ProteinMPNN output)
BINDER_SEQS = {
    "il7ra": "results/il7ra_seqrefine/mpnn/seqs/*.fa",
    "rbd": "results/variance_study/rbd/mpnn/seqs/*.fa",
    "pdl1": "results/variance_study/pdl1/mpnn/seqs/*.fa",
}


def load_binders(pattern, per_backbone=2):
    out = []
    for fa in sorted(glob.glob(pattern)):
        bb = os.path.splitext(os.path.basename(fa))[0]
        lines = [l.strip() for l in open(fa) if l.strip()]
        draws = [(lines[i], lines[i + 1]) for i in range(2, len(lines) - 1, 2)]
        for k, (hdr, seq) in enumerate(draws[:per_backbone]):
            out.append((f"{bb}__d{k}", seq))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-backbone", type=int, default=2,
                    help="draws per backbone to include")
    ap.add_argument("--max-per-pairing", type=int, default=70,
                    help="cap on binders per (binder-target, decoy-target) pairing")
    ap.add_argument("--out-dir", default="results/queue/null_panel")
    ap.add_argument("--seed", type=int, default=11)
    args = ap.parse_args()

    os.chdir(KIT)
    rng = random.Random(args.seed)

    # target MSAs, extracted once each
    msa = {}
    for t, pat in SOURCE_A3M.items():
        files = sorted(glob.glob(pat))
        if not files:
            print(f"  [skip] {t}: no source a3m"); continue
        t_len, _, rows = extract_target_msa(files[0])
        _, entries = read_a3m(files[0])
        q_first, _ = split_a3m_row(entries[0][1], t_len)
        tseq = "".join(c for c in q_first if not c.islower()).replace("-", "")
        msa[t] = (t_len, rows, tseq)
        print(f"  {t:6s} target len {t_len}, MSA {len(rows)} rows")

    binders = {t: load_binders(p, args.per_backbone) for t, p in BINDER_SEQS.items()}
    for t, b in binders.items():
        print(f"  {t:6s} binders available: {len(b)}")

    os.makedirs(args.out_dir, exist_ok=True)
    manifest, n = [], 0
    for src in binders:
        for dec in msa:
            if dec == src:
                continue                       # cross-pairings only
            pool = binders[src]
            if len(pool) > args.max_per_pairing:
                pool = rng.sample(pool, args.max_per_pairing)
            t_len, rows, tseq = msa[dec]
            for name, seq in pool:
                tag = f"{src}__{name}--NULLvs--{dec}"
                with open(os.path.join(args.out_dir, f"{tag}.a3m"), "w") as f:
                    f.write(build_hybrid(rows, t_len, tseq, seq, tag))
                manifest.append({"tag": tag, "binder_target": src, "decoy_target": dec,
                                 "binder": name, "sequence": seq, "label": "non_binder_by_design"})
                n += 1

    pd.DataFrame(manifest).to_csv(os.path.join(args.out_dir, "..", "null_panel_manifest.csv"),
                                  index=False)
    print(f"\n[OK] {n} cross-pairings -> {args.out_dir}")
    print(f"     ~{n * 50 / 3600:.1f} h at ~50 s/fold via the hybrid route")
    print("     label: non_binder_by_design (PREDICTED, not assayed)")


if __name__ == "__main__":
    main()
