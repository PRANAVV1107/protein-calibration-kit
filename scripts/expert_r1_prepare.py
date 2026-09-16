#!/usr/bin/env python3
"""Turn sampled sequences into hybrid-MSA fold inputs for the round-1 comparison.

Only the 10 RBD backbones that appear in the variance study are used: those are
the ones held out of DPO training, so the comparison is on structures neither
model was tuned on.

The target's precomputed alignment is reused and the binder is given as a single
sequence. That is the configuration this project measured at 63 s/design -- and
critically, the same one used for every other fold in the comparison, so any
difference between arms is the weights and not the scoring setup.
"""
import argparse, glob, json, os, re
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
RBD_TARGET = ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
              "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
              "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG")

def read_fa(p):
    """ProteinMPNN writes the native first, then sampled designs."""
    seqs, cur = [], None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            if cur is not None: seqs.append(cur)
            cur = ""
        elif cur is not None:
            cur += line
    if cur: seqs.append(cur)
    return seqs[1:] if len(seqs) > 1 else seqs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    keep = set(pd.read_csv(os.path.join(KIT, "results/variance_study/rbd/manifest.csv"))
               .backbone.astype(str).unique())
    indir = os.path.join(args.out, "fold_in"); os.makedirs(indir, exist_ok=True)
    rows, n = [], 0
    for arm in ("base", "tuned"):
        for fa in sorted(glob.glob(os.path.join(args.out, "mpnn_%s/seqs/*.fa" % arm))):
            bb = os.path.basename(fa)[:-3]
            if bb not in keep: continue
            for i, s in enumerate(read_fa(fa)):
                s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s.upper())
                if not s: continue
                tag = "%s__%s__s%d" % (arm, bb, i + 1)
                with open(os.path.join(indir, tag + ".a3m"), "w") as f:
                    f.write("#%d,%d\t1,1\n" % (len(RBD_TARGET), len(s)))
                    f.write(">101\n%s%s\n" % (RBD_TARGET, s))
                    f.write(">target\n%s%s\n" % (RBD_TARGET, "-" * len(s)))
                rows.append({"tag": tag, "arm": arm, "backbone": bb,
                             "draw": i + 1, "sequence": s, "length": len(s)})
                n += 1
    if not rows: raise SystemExit("[FAIL] no sequences matched the held-out backbones")
    pd.DataFrame(rows).to_csv(os.path.join(args.out, "manifest.csv"), index=False)
    d = pd.DataFrame(rows)
    print("[OK] %d inputs  base=%d tuned=%d  backbones=%d"
          % (n, (d.arm=="base").sum(), (d.arm=="tuned").sum(), d.backbone.nunique()))

if __name__ == "__main__":
    main()
