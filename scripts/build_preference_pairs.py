#!/usr/bin/env python3
"""Build within-backbone preference pairs for offline preference tuning.

The idea this rests on
----------------------
Online RL over this pipeline is unaffordable: every rollout needs a reward,
the only predictive reward we have is a high-accuracy ColabFold run at 4-8
minutes, and the cheap screen that would be fast enough scores BELOW chance
against measured binding (ipSAE 0.466, pDockQ2 0.461 on 164 Adaptyv designs).
Training a policy against an anti-predictive reward would select against
binders.

Offline preference tuning sidesteps that entirely. The scores already exist.
What is needed is pairs, and the variance study already generated exactly the
right shape: several independent ProteinMPNN draws per backbone, each folded
at high accuracy. Within one backbone the structure is held constant, so a
score difference between two sequences isolates the SEQUENCE choice -- which
is the same observation as this project's headline result, that within-backbone
variance dominates between-backbone variance (59-78% of total). The finding
and the training signal are one thing seen twice.

What counts as a pair
---------------------
Two sequences on the SAME backbone whose scores differ by more than --margin.
Pairs closer than that are dropped rather than labelled: our metrics reach only
AUC ~0.62 against measured binding, so a small score gap carries no reliable
preference and would teach noise. The margin is the one knob that trades pair
count against label quality.

Never pair across backbones. A cross-backbone comparison conflates sequence
quality with backbone quality, and backbone is the factor this design is trying
to hold constant.

Metric choice
-------------
Default pDockQ2_min, with composition penalty subtracted. On the Adaptyv set
the high-accuracy interface metrics were the ones that survived length
adjustment (ipTM 0.605, ipSAE 0.565, pDockQ2 0.581) while pLDDT did not
(0.551, CI containing 0.5). Pass --metric to use a different one; the pair
file records which was used so a model cannot be evaluated against a metric it
was not trained on without that being visible.
"""
import argparse
import glob
import itertools
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


def load_target(t):
    """Scored designs for one target, joined to their manifest."""
    man_p = os.path.join(KIT, "results/queue/v2refold/%s_manifest.csv" % t)
    fold_d = os.path.join(KIT, "results/queue/v2refold_folds_%s" % t)
    if not (os.path.exists(man_p) and os.path.isdir(fold_d)):
        return pd.DataFrame()
    man = pd.read_csv(man_p).set_index("tag")
    rows = []
    for f in glob.glob(os.path.join(fold_d, "*_scores_rank_001_*.json")):
        tag = os.path.basename(f).split("_scores_rank_001_")[0]
        if tag not in man.index:
            continue
        d = json.load(open(f))
        r = man.loc[tag]
        p = np.array(d.get("plddt", [0.0]))
        rows.append({
            "tag": tag, "target": t, "backbone": r["backbone"],
            "sequence": r["sequence"], "length": r["length"],
            "composition_penalty": float(r["composition_penalty"]),
            "ipTM": d.get("iptm"), "pTM": d.get("ptm"),
            "ipSAE_min": mn(d.get("ipsae")), "pDockQ2_min": mn(d.get("pdockq2")),
            "pLDDT": float(p.mean()),
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metric", default="pDockQ2_min",
                    choices=["pDockQ2_min", "ipSAE_min", "ipTM", "pLDDT"])
    ap.add_argument("--margin", type=float, default=0.05,
                    help="minimum score gap for a pair to be kept")
    ap.add_argument("--no-composition", action="store_true",
                    help="rank on the raw metric without the composition penalty")
    ap.add_argument("--targets", default="il7ra,pdl1,rbd")
    ap.add_argument("--holdout-target", default=None,
                    help="exclude this target entirely, for cross-target generalisation")
    ap.add_argument("--out", default="results/rl/preference_pairs.jsonl")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.chdir(KIT)
    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    if args.holdout_target:
        targets = [t for t in targets if t != args.holdout_target]

    frames = [load_target(t) for t in targets]
    frames = [f for f in frames if len(f)]
    if not frames:
        raise SystemExit("[FAIL] no v2 refold folds found yet -- run the queue first")
    df = pd.concat(frames, ignore_index=True)
    df = df[df[args.metric].notna()]
    df["score"] = df[args.metric] if args.no_composition else \
        df[args.metric] - df["composition_penalty"]

    print("=" * 76)
    print("PREFERENCE PAIRS -- within-backbone only")
    print("=" * 76)
    print("  metric  : %s%s" % (args.metric,
                                "" if args.no_composition else " - composition_penalty"))
    print("  margin  : %.3f" % args.margin)
    print("  targets : %s%s" % (", ".join(targets),
                                "   (held out: %s)" % args.holdout_target
                                if args.holdout_target else ""))
    print()

    pairs, skipped = [], 0
    for (t, bb), g in df.groupby(["target", "backbone"]):
        if len(g) < 2:
            continue
        for a, b in itertools.combinations(g.itertuples(), 2):
            gap = abs(a.score - b.score)
            if gap < args.margin:
                skipped += 1
                continue
            win, lose = (a, b) if a.score > b.score else (b, a)
            pairs.append({
                "target": t, "backbone": str(bb),
                "chosen": win.sequence, "rejected": lose.sequence,
                "chosen_tag": win.tag, "rejected_tag": lose.tag,
                "chosen_score": float(win.score), "rejected_score": float(lose.score),
                "gap": float(gap), "metric": args.metric,
                "composition_subtracted": not args.no_composition,
            })

    if not pairs:
        raise SystemExit("[FAIL] no pairs cleared the margin; lower --margin")

    pd.set_option("display.width", 100)
    pp = pd.DataFrame(pairs)
    print("%-8s %8s %10s %10s %10s" % ("target", "pairs", "backbones", "med gap", "max gap"))
    print("-" * 76)
    for t, g in pp.groupby("target"):
        print("%-8s %8d %10d %10.3f %10.3f"
              % (t, len(g), g.backbone.nunique(), g.gap.median(), g.gap.max()))
    print("-" * 76)
    print("%-8s %8d %10d %10.3f %10.3f"
          % ("TOTAL", len(pp), pp.backbone.nunique(), pp.gap.median(), pp.gap.max()))
    print("\n  dropped below margin: %d  (kept %.0f%% of candidate pairs)"
          % (skipped, 100 * len(pp) / (len(pp) + skipped)))

    # A sequence appearing as 'chosen' in one pair and 'rejected' in another is
    # normal and correct -- preference is relative within a backbone, not an
    # absolute label. Reported because it looks alarming if unexplained.
    both = set(pp.chosen_tag) & set(pp.rejected_tag)
    print("  sequences appearing on both sides: %d (expected; preference is relative)"
          % len(both))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        for r in pairs:
            f.write(json.dumps(r) + "\n")
    print("\n[OK] %d pairs -> %s" % (len(pairs), args.out))


if __name__ == "__main__":
    main()
