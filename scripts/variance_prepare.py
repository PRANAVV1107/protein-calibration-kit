#!/usr/bin/env python3
"""Build the fold set for the sequence-draw variance study.

Selects a seeded-random subset of backbones per target (no cherry-picking by
score -- that would bias the variance estimate) and writes one FASTA per target
pairing every draw against that target's sequence.
"""
import argparse
import glob
import os
import random
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_composition import composition_penalty  # noqa: E402

TARGET_SEQS = {
    "il7ra": ("DYSFSCYSQLEVNGSQHSLTCAFEDPDVNTTNLEFEICGALVEVKCLNFRKLQEIYFIETKKFLLIGKSNICVKVGEKSL"
              "TCKKIDLTTIVKPEAPFDLSVVYREGANDFVVTFNTSHLQKKYVKVLMHDVAYRQEKDENKWTHVNLSSTKLTLLQRKLQ"
              "PAAMYEIKVRSIPDHYFKGFWSEWSPSYYFRTP"),
    "rbd": ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
            "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
            "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG"),
    "pdl1": ("AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNA"
             "ALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA"),
}


def read_draws(seqs_dir):
    """Parse multi-sample ProteinMPNN .fa files."""
    rows = []
    for fa in sorted(glob.glob(os.path.join(seqs_dir, "*.fa"))):
        backbone = os.path.splitext(os.path.basename(fa))[0]
        lines = [l.strip() for l in open(fa) if l.strip()]
        for i in range(2, len(lines) - 1, 2):
            header, seq = lines[i], lines[i + 1]
            sample = "s0"
            for fld in header.lstrip(">").split(","):
                if "sample=" in fld:
                    sample = "s" + fld.split("=")[1].strip()
            rows.append({"backbone": backbone, "sample": sample, "sequence": seq})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--study-dir", required=True)
    ap.add_argument("--n-backbones", type=int, default=10)
    ap.add_argument("--draws", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1234)
    args = ap.parse_args()

    total = 0
    for name, target_seq in TARGET_SEQS.items():
        seqs_dir = os.path.join(args.study_dir, name, "mpnn", "seqs")
        if not os.path.isdir(seqs_dir):
            print(f"[WARN] {name}: no seqs dir, skipping")
            continue

        df = read_draws(seqs_dir)
        if df.empty:
            print(f"[WARN] {name}: no sequences parsed, skipping")
            continue

        # keep only backbones that produced the full number of draws, so every
        # backbone contributes equally to the variance estimate
        counts = df.groupby("backbone").size()
        complete = counts[counts >= args.draws].index.tolist()

        rng = random.Random(args.seed)
        chosen = sorted(rng.sample(complete, min(args.n_backbones, len(complete))))

        sub = df[df["backbone"].isin(chosen)].copy()
        # if a backbone somehow has extra draws, take the first N deterministically
        sub = sub.sort_values(["backbone", "sample"]).groupby("backbone", as_index=False).head(args.draws)
        sub["penalty"] = sub["sequence"].apply(composition_penalty)
        sub["length"] = sub["sequence"].str.len()
        sub["target"] = name

        out_fa = os.path.join(args.study_dir, name, "to_fold.fasta")
        with open(out_fa, "w") as f:
            for _, r in sub.iterrows():
                f.write(f">{r['backbone']}__{r['sample']}\n{target_seq}:{r['sequence']}\n")
        sub.to_csv(os.path.join(args.study_dir, name, "manifest.csv"), index=False)

        total += len(sub)
        print(f"  {name:6s}: {len(chosen)} backbones x {args.draws} draws = {len(sub)} folds  "
              f"(binder len {sub['length'].min()}-{sub['length'].max()}, "
              f"mean composition penalty {sub['penalty'].mean():.3f})")

    print(f"\n[OK] total folds queued: {total}   estimated {total*4.5/60:.1f} h at 4.5 min/fold")


if __name__ == "__main__":
    main()
