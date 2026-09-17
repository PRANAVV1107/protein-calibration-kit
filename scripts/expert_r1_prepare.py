#!/usr/bin/env python3
"""Fold inputs for the round-1 comparison, via make_hybrid_msa.py.

The first version hand-wrote a two-line a3m and every fold failed with
"MSA 0 must contain at least one sequence". ColabFold splits a paired a3m into
per-chain MSAs and rejects a chain with no alignment, and the query header must
be paired (">101<tab>102"). The same defect was fixed in mcts_design.py and not
propagated here, which cost this job a full run.

make_hybrid_msa.py already handles the paired header and a3m lowercase
insertion columns, so the target alignment is donated from an existing complex
a3m rather than reconstructed.

Only the 10 RBD backbones present in the variance study are used: those were
held out of DPO training, so the comparison is on structures neither model saw.
"""
import argparse, glob, os, re, subprocess
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
CF = "/root/miniconda3/envs/colabfold"
SRC_A3M = KIT + "/results/queue/v2refold/rbd/rbd__binder_0__s1.a3m"


def read_fa(p):
    """ProteinMPNN writes the native first, then the sampled designs."""
    seqs, cur = [], None
    for line in open(p):
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
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    keep = set(pd.read_csv(os.path.join(KIT, "results/variance_study/rbd/manifest.csv"))
               .backbone.astype(str).unique())
    indir = os.path.join(args.out, "fold_in")
    os.makedirs(indir, exist_ok=True)
    fasta = os.path.join(args.out, "binders.fasta")

    rows = []
    with open(fasta, "w") as fh:
        for arm in ("base", "tuned"):
            for fa in sorted(glob.glob(os.path.join(args.out, "mpnn_%s/seqs/*.fa" % arm))):
                bb = os.path.basename(fa)[:-3]
                if bb not in keep:
                    continue
                for i, s in enumerate(read_fa(fa)):
                    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s.upper())
                    if not s:
                        continue
                    tag = "%s__%s__s%d" % (arm, bb, i + 1)
                    fh.write(">" + tag + "\n" + s + "\n")
                    rows.append({"tag": tag, "arm": arm, "backbone": bb,
                                 "draw": i + 1, "sequence": s, "length": len(s)})
    if not rows:
        raise SystemExit("[FAIL] no sequences matched the held-out backbones")

    if not os.path.exists(SRC_A3M):
        raise SystemExit("[FAIL] source a3m missing: %s" % SRC_A3M)
    subprocess.run([CF + "/bin/python", KIT + "/scripts/make_hybrid_msa.py",
                    "--source-a3m", SRC_A3M, "--binders", fasta,
                    "--out-dir", indir], check=True,
                   env={**os.environ, "PYTHONPATH": ""})

    pd.DataFrame(rows).to_csv(os.path.join(args.out, "manifest.csv"), index=False)
    made = len(glob.glob(os.path.join(indir, "*.a3m")))
    d = pd.DataFrame(rows)
    print("[OK] %d sequences (base=%d tuned=%d, %d backbones) -> %d hybrid a3m"
          % (len(d), (d.arm == "base").sum(), (d.arm == "tuned").sum(),
             d.backbone.nunique(), made))


if __name__ == "__main__":
    main()
