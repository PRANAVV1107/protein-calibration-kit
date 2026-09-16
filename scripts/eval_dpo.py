#!/usr/bin/env python3
"""Preference accuracy of a DPO-tuned ProteinMPNN on pairs it never saw.

Why a separate script. dpo_finetune.py reports validation accuracy on a random
split of its own training pairs -- same targets, often the same backbones. That
measures whether the model memorised the ordering it was shown, which it will,
and says nothing about generalisation. This evaluates on a target excluded from
training entirely.

What this still does NOT show. Preference accuracy asks whether the tuned model
assigns higher likelihood to the sequence that scored better. It does not show
that sequences SAMPLED from the tuned model are better, which is the claim that
matters and which needs fresh generation followed by high-accuracy folding. A
model can rank existing candidates well while sampling no better than before.
Treat a good number here as permission to run that experiment, not as its
result.

Chance is 0.5. With n pairs the standard error is about 0.5/sqrt(n), so on ~98
pairs anything inside roughly 0.40-0.60 is not distinguishable from chance.
"""
import argparse
import copy
import json
import os
import sys

import numpy as np
import torch

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
MPNN_DIR = os.path.join(ROOT, "external_tools/ProteinMPNN")
sys.path.insert(0, MPNN_DIR)
sys.path.insert(0, os.path.join(KIT, "scripts"))

from dpo_finetune import (load_model, load_backbones, seq_logprob,  # noqa: E402
                          SOURCES)
from protein_mpnn_utils import tied_featurize  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--tuned", required=True, help="DPO checkpoint .pt")
    ap.add_argument("--base",
                    default=os.path.join(MPNN_DIR, "vanilla_model_weights/v_48_020.pt"))
    ap.add_argument("--orders", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.chdir(KIT)
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    pairs = [json.loads(l) for l in open(args.pairs)]
    targets = sorted({p["target"] for p in pairs})
    bb = load_backbones(targets)
    pairs = [p for p in pairs if (p["target"], p["backbone"]) in bb]
    if not pairs:
        raise SystemExit("[FAIL] no pairs resolvable to a backbone")

    tuned = load_model(args.tuned, device)
    base = load_model(args.base, device)
    tuned.eval()
    base.eval()

    cache = {}
    for key, (protein, chains) in bb.items():
        cd = {protein["name"]: (chains[0], chains[1])}
        o = tied_featurize([copy.deepcopy(protein)], device, cd,
                           None, None, None, None, None, ca_only=False)
        feats = (o[0], o[1], o[2], o[4], o[10], o[12], o[5])
        cache[key] = (feats, o[2] * o[4] * o[10])

    def acc(model):
        """Fraction of pairs where the model prefers the better sequence."""
        hits = []
        with torch.no_grad():
            for p in pairs:
                feats, dmask = cache[(p["target"], p["backbone"])]
                rnd = [torch.randn(feats[3].shape, device=device)
                       for _ in range(args.orders)]
                lw = seq_logprob(model, feats, p["chosen"], rnd, dmask)
                ll = seq_logprob(model, feats, p["rejected"], rnd, dmask)
                hits.append(int(float(lw - ll) > 0))
        return np.array(hits)

    h_t, h_b = acc(tuned), acc(base)
    n = len(pairs)
    se = 0.5 / np.sqrt(n)

    print("=" * 72)
    print("HELD-OUT PREFERENCE ACCURACY")
    print("=" * 72)
    print("  pairs %d   targets: %s" % (n, ", ".join(targets)))
    print("  chance 0.500, se ~%.3f, so ~%.3f-%.3f is indistinguishable\n"
          % (se, 0.5 - 1.96 * se, 0.5 + 1.96 * se))
    print("  untuned ProteinMPNN : %.3f" % h_b.mean())
    print("  DPO-tuned           : %.3f" % h_t.mean())
    print("  difference          : %+.3f" % (h_t.mean() - h_b.mean()))

    # paired: the same pairs under both models, so McNemar is the right test
    b = int(((h_t == 1) & (h_b == 0)).sum())
    c = int(((h_t == 0) & (h_b == 1)).sum())
    print("\n  tuned right / base wrong : %d" % b)
    print("  tuned wrong / base right : %d" % c)
    if b + c > 0:
        from math import erfc, sqrt
        z = abs(b - c) / sqrt(b + c)
        print("  McNemar z %.2f  p %.3f" % (z, erfc(z / sqrt(2))))
    else:
        print("  models disagree on nothing")

    print("\n  Ranking existing candidates is not the same as SAMPLING better")
    print("  ones. That needs fresh generation and high-accuracy folding.")


if __name__ == "__main__":
    main()
