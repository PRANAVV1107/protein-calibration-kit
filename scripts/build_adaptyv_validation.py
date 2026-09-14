#!/usr/bin/env python3
"""Build the Adaptyv EGFR validation set -- the only measured outcomes here.

Every other result in this repository validates AlphaFold's opinion against
AlphaFold's opinion. This set has experimental binding calls: 380 designs, 55
binders and 325 non-binders, with K_D between 1.2 nM and 10 uM for the binders.

It answers the question the rest of the project cannot: does the selection
criterion (on-target minus worst off-target minus composition penalty) predict
MEASURED binding, or only itself?

Target: EGFR domain III, PDB 1YY9 chain A residues 310-510 (201 aa, no crystal
gaps). Correct because 1YY9 is EGFR bound to cetuximab, the cetuximab epitope
spans EGFR 349-473 by contact analysis, and Adaptyv's own reference positive is
Cetuximab_scFv. Structure verified against its own TITLE/COMPND records.

Folded through the hybrid-MSA route: the EGFR alignment is computed once and
reused across all 380 designs, so this costs ~50 s per design and no per-design
MSA-server calls.

Note on heterogeneity: these designs come from many different methods (the
competition's `design_models` column), not one pipeline. That makes them a
harder and more honest test than re-scoring our own RFdiffusion output, because
it tests generalisation rather than fit to one generator's quirks.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed-a3m", default="",
                    help="a3m containing the EGFR target alignment (from the seed fold)")
    ap.add_argument("--out-dir", default="results/queue/adaptyv_val")
    ap.add_argument("--min-len", type=int, default=25,
                    help="skip designs shorter than this; a few are 13 aa fragments")
    args = ap.parse_args()

    os.chdir(KIT)
    df = pd.read_csv("adaptyv_egfr_round2.csv")
    df["binding"] = df["binding"].astype(str).str.strip().str.lower()
    lab = df[df["binding"].isin(["true", "false"])].copy()
    lab = lab[lab["sequence"].notna()]
    lab["length"] = lab["sequence"].str.len()
    dropped = (lab["length"] < args.min_len).sum()
    lab = lab[lab["length"] >= args.min_len]

    print("Adaptyv EGFR round 2")
    print("  labelled designs      : %d  (%d binders / %d non-binders)"
          % (len(lab), (lab.binding == "true").sum(), (lab.binding == "false").sum()))
    print("  dropped as too short  : %d  (< %d aa)" % (dropped, args.min_len))
    print("  length range          : %d-%d aa" % (lab.length.min(), lab.length.max()))

    if not args.seed_a3m:
        # stage 1 hasn't run yet: emit the seed FASTA and stop
        target = open("benchmark/decoys/EGFRd3_seq.txt").read().strip()
        os.makedirs("results/queue", exist_ok=True)
        with open("results/queue/adaptyv_seed.fasta", "w") as f:
            f.write(">seed_EGFRd3\n%s:%s\n" % (target, lab.iloc[0]["sequence"]))
        print("\n[stage 1] wrote results/queue/adaptyv_seed.fasta")
        print("          one high-accuracy fold to generate the EGFR alignment")
        return

    t_len, _, rows = extract_target_msa(args.seed_a3m)
    _, entries = read_a3m(args.seed_a3m)
    q_first, _ = split_a3m_row(entries[0][1], t_len)
    tseq = "".join(c for c in q_first if not c.islower()).replace("-", "")
    print("\n  EGFR alignment: %d rows, target %d aa (reused for every design)"
          % (len(rows), t_len))

    os.makedirs(args.out_dir, exist_ok=True)
    man = []
    for i, r in enumerate(lab.itertuples()):
        tag = "adaptyv_%03d" % i
        with open(os.path.join(args.out_dir, tag + ".a3m"), "w") as f:
            f.write(build_hybrid(rows, t_len, tseq, r.sequence, tag))
        man.append({"tag": tag, "name": r.name, "sequence": r.sequence,
                    "length": len(r.sequence), "binding": r.binding,
                    "kd": getattr(r, "kd", None),
                    "binding_strength": getattr(r, "binding_strength", None),
                    "expression": getattr(r, "expression", None),
                    "adaptyv_iptm": getattr(r, "iptm", None),
                    "adaptyv_plddt": getattr(r, "plddt", None),
                    "composition_penalty": composition_penalty(r.sequence)})
    pd.DataFrame(man).to_csv("results/queue/adaptyv_val_manifest.csv", index=False)
    print("[OK] %d hybrid inputs -> %s  (~%.1f h at ~50 s/fold)"
          % (len(man), args.out_dir, len(man) * 50 / 3600))


if __name__ == "__main__":
    main()
