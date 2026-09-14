#!/usr/bin/env python3
"""Build the candidate x family-decoy panel, reusing each decoy's alignment.

Stage 1 of queue_v6 folds one binder against each decoy at high accuracy purely
to obtain that decoy's MSA. This script harvests those alignments and builds a
hybrid a3m for every remaining pairing, so the rest of the panel costs ~50 s per
fold and no further MSA-server calls.
"""
import glob
import os
import sys

import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
sys.path.insert(0, os.path.join(KIT, "scripts"))
from make_hybrid_msa import build_hybrid, extract_target_msa, read_a3m, split_a3m_row  # noqa: E402

DECOYS = ["gammaC", "IL4Ra", "IL21R", "IgGFc", "albumin"]


def main():
    os.chdir(KIT)
    out_dir = "results/queue/family_panel"
    os.makedirs(out_dir, exist_ok=True)

    # candidates: the IL-7Ra finalists we have on-target scores for
    man = pd.read_csv("results/il7ra_seqrefine/manifest.csv")
    cands = [(f"{r['backbone']}__{r['sample']}", r["sequence"]) for _, r in man.iterrows()]
    print("  candidates: %d" % len(cands))

    rows, n = [], 0
    for dec in DECOYS:
        hits = sorted(glob.glob("results/queue/family_seed/*%s*.a3m" % dec))
        if not hits:
            print("  [skip] %s: no alignment produced in stage 1" % dec)
            continue
        t_len, _, msa_rows = extract_target_msa(hits[0])
        _, entries = read_a3m(hits[0])
        q_first, _ = split_a3m_row(entries[0][1], t_len)
        tseq = "".join(c for c in q_first if not c.islower()).replace("-", "")
        print("  %-8s target %3d res, MSA %d rows" % (dec, t_len, len(msa_rows)))

        for name, seq in cands:
            tag = "%s--FAMvs--%s" % (name, dec)
            with open(os.path.join(out_dir, tag + ".a3m"), "w") as f:
                f.write(build_hybrid(msa_rows, t_len, tseq, seq, tag))
            rows.append({"tag": tag, "candidate": name, "decoy": dec,
                         "sequence": seq, "label": "non_binder_by_design"})
            n += 1

    pd.DataFrame(rows).to_csv("results/queue/family_panel_manifest.csv", index=False)
    print("[OK] %d pairings -> %s  (~%.1f h at ~50 s/fold)" % (n, out_dir, n * 50 / 3600))


if __name__ == "__main__":
    main()
