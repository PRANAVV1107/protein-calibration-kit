#!/usr/bin/env python3
"""Build the v2 refold: the complete design space, scored by the corrected method.

Why refold at all
-----------------
The IL-7Ra results currently rest on 12 of 50 backbones, and those 12 were
selected by the fast screen -- which was measured to select WORSE than random
(-0.162 ipTM against random choice, because without an MSA the target folds at
pLDDT 27 rather than 84). The other 38 backbones were never evaluated. On this
target the fast screen's bottom decile contained the joint-best design in the
whole campaign, at 0.870.

So the current "best candidate" is the best of an arbitrarily chosen quarter of
the design space. This rebuilds the whole thing.

What makes it affordable
------------------------
The hybrid route -- precomputed target alignment, binder as a single sequence --
runs at roughly 50 s per fold against 4.5 min for high accuracy, and makes no
MSA-server calls because the alignment is computed once per target and reused.
At that price the entire IL-7Ra design space is about one overnight run.

Every choice below has a measurement behind it
----------------------------------------------
  draws per backbone    8, not 1     59-78% of interface-score variance is
                                     within-backbone, i.e. attributable to which
                                     sequence was drawn
  alanine bias -1.0     applied      raising temperature alone left mean alanine
                                     at ~20% with a 16-residue run; the bias
                                     brought it to 8.8%, natural abundance
  composition filter    before fold  design costs ~1.5 s, folding ~50 s; filter
                                     first
  fast screen           removed      no selective signal (p=0.970 across designs;
                                     rho=-0.14 within a backbone)
  monomer filter        removed      0 of 183 designs disordered, and monomer
                                     pLDDT does not predict complex score
  ranking metric        deferred     set from the Adaptyv validation, which tests
                                     each metric against 365 MEASURED binding
                                     outcomes rather than against itself
"""
import argparse
import glob
import os
import sys

import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
sys.path.insert(0, os.path.join(KIT, "scripts"))
from make_hybrid_msa import build_hybrid, extract_target_msa, read_a3m, split_a3m_row  # noqa: E402
from seq_composition import composition_penalty  # noqa: E402

TARGETS = {
    "il7ra": {
        "seqs": "results/il7ra_seqrefine/mpnn/seqs/*.fa",
        "a3m": "results/variance_study/il7ra/folds/*.a3m",
    },
    "pdl1": {
        "seqs": "results/variance_study/pdl1/mpnn/seqs/*.fa",
        "a3m": "results/variance_study/pdl1/folds/*.a3m",
    },
    "rbd": {
        "seqs": "results/variance_study/rbd/mpnn/seqs/*.fa",
        "a3m": "results/variance_study/rbd/folds/*.a3m",
    },
}


def draws(pattern):
    out = []
    for fa in sorted(glob.glob(pattern)):
        bb = os.path.splitext(os.path.basename(fa))[0]
        lines = [l.strip() for l in open(fa) if l.strip()]
        for k in range(2, len(lines) - 1, 2):
            hdr, seq = lines[k], lines[k + 1]
            s = "s0"
            for fld in hdr.lstrip(">").split(","):
                if "sample=" in fld:
                    s = "s" + fld.split("=")[1].strip()
            out.append((bb, s, seq))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="il7ra",
                    help="comma-separated: il7ra, pdl1, rbd")
    ap.add_argument("--max-penalty", type=float, default=0.10,
                    help="composition ceiling; free, applied before any folding")
    ap.add_argument("--out-root", default="results/queue/v2refold")
    args = ap.parse_args()

    os.chdir(KIT)
    total = 0
    summary = []
    for t in [x.strip() for x in args.targets.split(",")]:
        if t not in TARGETS:
            print("  [skip] unknown target %r" % t); continue
        cfg = TARGETS[t]
        a3ms = sorted(glob.glob(cfg["a3m"]))
        if not a3ms:
            print("  [skip] %s: no source alignment" % t); continue
        t_len, _, rows = extract_target_msa(a3ms[0])
        _, entries = read_a3m(a3ms[0])
        q_first, _ = split_a3m_row(entries[0][1], t_len)
        tseq = "".join(c for c in q_first if not c.islower()).replace("-", "")

        all_draws = draws(cfg["seqs"])
        if not all_draws:
            print("  [skip] %s: no sequences" % t); continue

        kept, dropped = [], 0
        for bb, s, seq in all_draws:
            p = composition_penalty(seq)
            if p <= args.max_penalty:
                kept.append((bb, s, seq, p))
            else:
                dropped += 1

        out_dir = os.path.join(args.out_root, t)
        os.makedirs(out_dir, exist_ok=True)
        man = []
        for bb, s, seq, p in kept:
            tag = "%s__%s__%s" % (t, bb, s)
            with open(os.path.join(out_dir, tag + ".a3m"), "w") as f:
                f.write(build_hybrid(rows, t_len, tseq, seq, tag))
            man.append({"tag": tag, "target": t, "backbone": bb, "sample": s,
                        "sequence": seq, "length": len(seq), "composition_penalty": p})
        pd.DataFrame(man).to_csv(os.path.join(args.out_root, "%s_manifest.csv" % t), index=False)

        n_bb = len({b for b, _, _, _ in kept})
        print("  %-6s %4d draws over %2d backbones -> %4d kept, %3d dropped on composition"
              % (t, len(all_draws), n_bb, len(kept), dropped))
        print("         target %d aa, MSA %d rows (computed once, reused)" % (t_len, len(rows)))
        total += len(kept)
        summary.append((t, len(kept), n_bb))

    print("\n[OK] %d folds queued -> %s" % (total, args.out_root))
    print("     ~%.1f h at ~50 s/fold via the hybrid route, no MSA-server calls" % (total * 50 / 3600))
    print("\n     For comparison, the current IL-7Ra results rest on 12 of 50")
    print("     backbones, chosen by a screen that selected worse than random.")


if __name__ == "__main__":
    main()
