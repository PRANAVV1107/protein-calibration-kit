#!/usr/bin/env python3
"""Pick finalists from the v2 hybrid screen for high-accuracy rescoring.

The hybrid screen is cheap but imperfect: on RBD it recovered about 47% of the
gap between random selection and an oracle. So it is used to narrow, not to
decide -- the top of its ranking gets high-accuracy folds, and those are what a
candidate is judged on.

Ranking metric: pDockQ2_min by default, because on the within-backbone
sequence-selection task it beat both ipTM (+0.166 against +0.141 over random
choice) and ipSAE (+0.077). If the Adaptyv validation -- which tests each metric
against 365 MEASURED binding outcomes rather than against itself -- prefers a
different metric, pass --metric to match it. That file is the only ground truth
available here, so it should win any disagreement.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"


def mn(x):
    if isinstance(x, dict):
        v = [float(i) for i in x.values()]
        return min(v) if v else np.nan
    return float(x) if x is not None else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", default="pDockQ2_min",
                    choices=["pDockQ2_min", "ipSAE_min", "ipTM"])
    ap.add_argument("--per-target", type=int, default=20,
                    help="finalists per target for high-accuracy rescoring")
    ap.add_argument("--per-backbone", type=int, default=1,
                    help="cap per backbone, so one lucky backbone cannot fill the list")
    args = ap.parse_args()

    os.chdir(KIT)
    rows = []
    for t in ("il7ra", "pdl1", "rbd"):
        man_p = "results/queue/v2refold/%s_manifest.csv" % t
        fold_d = "results/queue/v2refold_folds_%s" % t
        if not (os.path.exists(man_p) and os.path.isdir(fold_d)):
            continue
        man = pd.read_csv(man_p).set_index("tag")
        for f in glob.glob(os.path.join(fold_d, "*_scores_rank_001_*.json")):
            tag = os.path.basename(f).split("_scores_rank_001_")[0]
            if tag not in man.index:
                continue
            d = json.load(open(f))
            r = man.loc[tag]
            rows.append({
                "tag": tag, "target": t, "backbone": r["backbone"],
                "sequence": r["sequence"], "length": r["length"],
                "composition_penalty": r["composition_penalty"],
                "ipTM": d.get("iptm"), "pTM": d.get("ptm"),
                "ipSAE_min": mn(d.get("ipsae")), "pDockQ2_min": mn(d.get("pdockq2")),
            })

    if not rows:
        raise SystemExit("[FAIL] no v2 fold results found")

    df = pd.DataFrame(rows)
    df["score"] = df[args.metric] - df["composition_penalty"]
    df.to_csv("results/queue/v2_screen_all.csv", index=False)

    print("v2 screen: %d folds across %d targets" % (len(df), df.target.nunique()))
    for t, g in df.groupby("target"):
        print("  %-6s n=%4d  %s median %.3f  max %.3f  backbones %d"
              % (t, len(g), args.metric, g[args.metric].median(),
                 g[args.metric].max(), g.backbone.nunique()))

    picks = (df.sort_values("score", ascending=False)
               .groupby(["target", "backbone"], as_index=False).head(args.per_backbone)
               .sort_values("score", ascending=False)
               .groupby("target", as_index=False).head(args.per_target))

    targets = {
        "il7ra": open("benchmark/il7ra_target_seq.txt").read().strip(),
        "pdl1": ("AFTVTVPKDLYVVEYGSNMTIECKFPVEKQLDLAALIVYWEMEDKNIIQFVHGEEDLKVQHSSYRQRARLLKDQLSLGNA"
                 "ALQITDVKLQDAGVYRCMISYGGADYKRITVKVNA"),
        "rbd": ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
                "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
                "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG"),
    }
    with open("results/queue/v2_finalists.fasta", "w") as f:
        for _, r in picks.iterrows():
            f.write(">%s\n%s:%s\n" % (r["tag"], targets[r["target"]], r["sequence"]))
    picks.to_csv("results/queue/v2_finalists.csv", index=False)

    print("\n[OK] %d finalists (<=%d per target, <=%d per backbone) ranked on %s - composition"
          % (len(picks), args.per_target, args.per_backbone, args.metric))
    print("     ~%.1f h of high-accuracy folding" % (len(picks) * 4.5 / 60))


if __name__ == "__main__":
    main()
