#!/usr/bin/env python3
"""Print the state of every fold set in this project.

A file rather than an inline heredoc: nested quoting through
wsl.exe -> bash -c -> python mangles f-strings and shell variables.
"""
import glob
import os
import subprocess
import time

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"

JOBS = [
    ("variance il7ra", "results/variance_study/il7ra/folds", 60),
    ("variance rbd", "results/variance_study/rbd/folds", 60),
    ("variance pdl1", "results/variance_study/pdl1/folds", 60),
    ("seqrefine candidates", "results/il7ra_seqrefine/folds", 15),
    ("decoys (5 cand)", "results/il7ra_decoys", 10),
    ("rbd_null", "results/queue/rbd_null", 20),
    ("rewardchannel rbd", "results/queue/rewardchannel_rbd", 60),
    ("rewardchannel il7ra", "results/queue/rewardchannel_il7ra", 60),
    ("hybrid rbd", "results/queue/hybrid_rbd_folds", 60),
    ("template rbd", "results/queue/template_rbd_folds", 60),
    ("monomers", "results/queue/monomers", 183),
    ("il7ra rescreen", "results/queue/il7ra_rescreen_folds", 114),
    ("decoys full (15 cand)", "results/queue/il7ra_decoys_full_folds", 30),
]

INPUTS = [
    "results/queue/il7ra_rescreen.fasta",
    "results/queue/il7ra_decoys_full.fasta",
    "results/queue/rewardchannel_rbd.fasta",
    "results/queue/templates_rbd/rbdt.pdb",
]


def main():
    os.chdir(KIT)
    total = done = 0
    print("=== fold sets ===")
    for name, d, exp in JOBS:
        n = len(glob.glob(os.path.join(d, "*_scores_rank_001_*.json")))
        total += exp
        done += min(n, exp)
        filled = int(16 * n / exp) if exp else 0
        bar = "#" * min(filled, 16) + "." * max(0, 16 - filled)
        flag = "" if n >= exp else ("  <- pending" if n == 0 else "  <- partial")
        print("  %-22s [%s] %3d/%-4d%s" % (name, bar, n, exp, flag))
    print("  %-22s %d/%d folds complete" % ("TOTAL", done, total))

    print("\n=== prepared inputs ===")
    for p in INPUTS:
        if not os.path.exists(p):
            print("  %-32s MISSING" % os.path.basename(p))
            continue
        entries = sum(1 for l in open(p) if l.startswith(">")) if p.endswith("fasta") else "-"
        print("  %-32s %6.0f KB  entries=%s" % (os.path.basename(p), os.path.getsize(p) / 1024, entries))

    print("\n=== running ===")
    env = dict(os.environ, LINES="50", COLUMNS="300")
    ps = subprocess.run(["ps", "-eo", "args"], capture_output=True, text=True, env=env).stdout
    folds = [l for l in ps.splitlines() if "envs/colabfold/bin/colabfold_batch" in l]
    queues = [l for l in ps.splitlines() if "queue_v" in l and "grep" not in l]
    print("  colabfold jobs : %d" % len(folds))
    print("  queue scripts  : %d" % len(queues))
    for q in queues:
        print("    " + q.strip()[:80])
    print("  time           : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))


if __name__ == "__main__":
    main()
