#!/usr/bin/env python3
"""Fold inputs for the noise PoC: the 10 variance-study backbones per target.

Those backbones are used because each already has 6 high-accuracy draws
recorded, so the new 16 low-cost draws can be compared against an independent
measurement of the same structure if that is ever wanted.
"""
import argparse, glob, os, re
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
SRC_A3M = {
    "rbd":   KIT + "/results/queue/v2refold/rbd/rbd__binder_0__s1.a3m",
    "il7ra": KIT + "/results/queue/v2refold/il7ra/il7ra__il7ra_binder_0__s1.a3m",
}

def read_fa(p):
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
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    indir = os.path.join(a.out, "fold_in"); os.makedirs(indir, exist_ok=True)
    rows = []
    for t in ("rbd", "il7ra"):
        mp = os.path.join(KIT, "results/variance_study/%s/manifest.csv" % t)
        if not os.path.exists(mp): continue
        keep = set(pd.read_csv(mp).backbone.astype(str).unique())
        tmp = os.path.join(a.out, "tmp_%s.fasta" % t)
        n = 0
        with open(tmp, "w") as fh:
            for fa in sorted(glob.glob(os.path.join(a.out, "mpnn_%s/seqs/*.fa" % t))):
                bb = os.path.basename(fa)[:-3]
                if bb not in keep: continue
                for i, s in enumerate(read_fa(fa)):
                    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s.upper())
                    if not s: continue
                    tag = "%s__%s__d%02d" % (t, bb, i + 1)
                    fh.write(">" + tag + "\n" + s + "\n")
                    rows.append({"tag": tag, "target": t, "backbone": bb,
                                 "draw": i + 1, "sequence": s, "length": len(s)})
                    n += 1
        if n:
            os.system("%s/bin/python %s/scripts/make_hybrid_msa.py --source-a3m %s "
                      "--binders %s --out-dir %s > /dev/null 2>&1"
                      % ("/root/miniconda3/envs/colabfold", KIT, SRC_A3M[t], tmp, indir))
        print("[%s] %d sequences over %d backbones" % (t, n, len(keep)))
    if not rows: raise SystemExit("[FAIL] nothing prepared")
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "manifest.csv"), index=False)
    print("[OK] %d inputs" % len(rows))

if __name__ == "__main__":
    main()
