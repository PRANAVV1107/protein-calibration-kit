#!/usr/bin/env python3
"""Build monomer-foldability inputs: each designed binder folded ALONE.

Why this matters: a de novo binder that only folds when bound to its target is
likely disordered in solution -- it will express poorly, aggregate, or simply
not be the thing the complex prediction implied. Complex ipTM cannot see this,
because the target scaffolds the binder in every complex prediction.

Why it is cheap: a designed sequence has no natural homologs, so an MSA search
returns essentially nothing useful. single_sequence mode is the CORRECT setting
here rather than a shortcut -- which also means zero MSA-server load. And an
80-120 residue monomer is far smaller than a 250-315 residue complex.

Outputs one FASTA per source set, plus a manifest mapping each monomer back to
its complex-fold score so the two can be correlated.
"""
import glob
import json
import os

import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
QUEUE = os.path.join(KIT, "results", "queue")


def variance_sequences():
    """All variance-study sequences that already have a complex ipTM."""
    rows = []
    for t in ("il7ra", "rbd", "pdl1"):
        man_p = os.path.join(KIT, "results", "variance_study", t, "manifest.csv")
        fold_d = os.path.join(KIT, "results", "variance_study", t, "folds")
        if not os.path.exists(man_p):
            continue
        man = pd.read_csv(man_p)
        scored = {}
        for f in glob.glob(os.path.join(fold_d, "*_scores_rank_001_*.json")):
            n = os.path.basename(f).split("_scores_rank_001_")[0]
            scored[n] = json.load(open(f))["iptm"]
        for _, r in man.iterrows():
            name = f"{r['backbone']}__{r['sample']}"
            if name in scored:
                rows.append({"monomer_name": f"{t}__{name}", "target": t,
                             "sequence": r["sequence"], "complex_ipTM": scored[name],
                             "length": len(r["sequence"]), "source": "variance_study"})
    return rows


def candidate_sequences():
    """The IL-7Ra finalists from the sequence-refinement run."""
    man_p = os.path.join(KIT, "results", "il7ra_seqrefine", "manifest.csv")
    ranked_p = os.path.join(KIT, "results", "il7ra_seqrefine", "ranked.csv")
    if not (os.path.exists(man_p) and os.path.exists(ranked_p)):
        return []
    man = pd.read_csv(man_p)
    ranked = pd.read_csv(ranked_p).set_index("name")
    rows = []
    for _, r in man.iterrows():
        name = f"{r['backbone']}__{r['sample']}"
        if name in ranked.index:
            rows.append({"monomer_name": f"cand__{name}", "target": "il7ra",
                         "sequence": r["sequence"],
                         "complex_ipTM": ranked.loc[name, "ipTM"],
                         "length": len(r["sequence"]), "source": "candidate"})
    return rows


def main():
    os.makedirs(QUEUE, exist_ok=True)
    rows = variance_sequences() + candidate_sequences()
    if not rows:
        raise SystemExit("[FAIL] no sequences found")

    df = pd.DataFrame(rows).drop_duplicates(subset=["sequence", "target"])
    out_fa = os.path.join(QUEUE, "monomers.fasta")
    with open(out_fa, "w") as f:
        for _, r in df.iterrows():
            f.write(f">{r['monomer_name']}\n{r['sequence']}\n")
    df.to_csv(os.path.join(QUEUE, "monomers_manifest.csv"), index=False)

    print(f"[OK] {len(df)} monomers -> {out_fa}")
    for t, g in df.groupby("target"):
        print(f"     {t:6s}: {len(g):3d}  (length {g['length'].min()}-{g['length'].max()})")
    # monomers are small and single-sequence; ~1 min each rather than ~2.5
    print(f"     estimated ~{len(df)*1.0/60:.1f} h at ~1 min/monomer, no MSA server")


if __name__ == "__main__":
    main()
