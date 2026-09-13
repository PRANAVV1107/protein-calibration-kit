#!/usr/bin/env python3
"""Build the FASTA inputs for the overnight queue.

Two independent questions, both cheap relative to what they decide:

  A. reward-channel test -- does fast-screen ipTM rank SEQUENCES within a
     fixed backbone? The RL loop uses fast-screen ipTM as its reward, and the
     fast screen was shown to carry no signal ACROSS designs (p=0.970). Those
     are different questions. This re-scores, in fast-screen mode, the exact
     sequences the variance study already scored at high accuracy, so the two
     can be correlated per backbone. Uses single_sequence mode -- no MSA
     server calls, so it cannot be throttled.

  B. fast-screen null replication on RBD -- the null was measured on IL-7Ra
     only (n=20). This folds RBD's fast-screen top-10 and bottom-10 ORIGINAL
     campaign sequences at high accuracy, giving a second, independent target.
"""
import argparse
import glob
import json
import os

import pandas as pd

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")

RBD_TARGET = ("TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
              "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
              "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG")


def build_reward_channel(targets):
    """A: re-score variance-study sequences in fast-screen mode."""
    total = 0
    for t in targets:
        man_path = os.path.join(KIT, "results", "variance_study", t, "manifest.csv")
        fold_dir = os.path.join(KIT, "results", "variance_study", t, "folds")
        if not os.path.exists(man_path):
            print(f"  [skip] {t}: no manifest")
            continue
        man = pd.read_csv(man_path)

        # only sequences that actually have a high-accuracy score to compare against
        have = set()
        for f in glob.glob(os.path.join(fold_dir, "*_scores_rank_001_*.json")):
            have.add(os.path.basename(f).split("_scores_rank_001_")[0])

        target_seq = man["target_seq"].iloc[0] if "target_seq" in man.columns else None
        if target_seq is None:
            # manifest stores only binder seqs; recover target from the fold fasta
            fa = os.path.join(KIT, "results", "variance_study", t, "to_fold.fasta")
            target_seq = open(fa).readlines()[1].split(":")[0].strip()

        out = os.path.join(KIT, "results", "queue", f"rewardchannel_{t}.fasta")
        n = 0
        with open(out, "w") as f:
            for _, r in man.iterrows():
                name = f"{r['backbone']}__{r['sample']}"
                if name not in have:
                    continue
                f.write(f">{name}\n{target_seq}:{r['sequence']}\n")
                n += 1
        total += n
        print(f"  {t:6s}: {n} sequences -> {os.path.basename(out)}")
    print(f"  [A] total {total} fast-screen folds (~{total*2.5/60:.1f} h, no MSA server)")


def build_rbd_null():
    """B: RBD fast-screen top-10 vs bottom-10, folded at high accuracy."""
    d = json.load(open(os.path.join(ROOT, "real_test_results/rbd_200_batch/ranking_all_200.json")))
    top = sorted(d, key=lambda x: -x["iptm"])[:10]
    bot = sorted(d, key=lambda x: x["iptm"])[:10]

    out = os.path.join(KIT, "results", "queue", "rbd_null.fasta")
    meta = []
    with open(out, "w") as f:
        for grp, items in (("top", top), ("bot", bot)):
            for x in items:
                p = os.path.join(ROOT, "real_test_results/rbd_200_batch",
                                 x["batch"], "proteinmpnn/output/seqs", f"{x['design']}.fa")
                seq = open(p).readlines()[3].strip()
                # design names repeat across batches -- qualify with batch
                name = f"{grp}__{x['batch']}__{x['design']}"
                f.write(f">{name}\n{RBD_TARGET}:{seq}\n")
                meta.append({"name": name, "group": grp, "batch": x["batch"],
                             "design": x["design"], "fast_ipTM": x["iptm"],
                             "fast_rank": x["rank"], "length": len(seq)})
    pd.DataFrame(meta).to_csv(os.path.join(KIT, "results", "queue", "rbd_null_manifest.csv"), index=False)
    print(f"  [B] 20 high-accuracy folds -> {os.path.basename(out)} (~1.5 h)")
    print(f"      top-10 fast ipTM {top[-1]['iptm']}-{top[0]['iptm']}, "
          f"bottom-10 {bot[0]['iptm']}-{bot[-1]['iptm']}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="rbd,il7ra")
    args = ap.parse_args()

    os.makedirs(os.path.join(KIT, "results", "queue"), exist_ok=True)
    print("A. reward-channel test (fast-screen re-score of variance-study sequences)")
    build_reward_channel([t.strip() for t in args.targets.split(",")])
    print("\nB. fast-screen null replication on RBD")
    build_rbd_null()


if __name__ == "__main__":
    main()
