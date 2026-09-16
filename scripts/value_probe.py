#!/usr/bin/env python3
"""Is a cheap value function learnable? The load-bearing test for AlphaZero-style search.

Why this is the experiment that decides the architecture
--------------------------------------------------------
A naive MCTS over this pipeline costs ~5 min per leaf, because every leaf needs
a fold. That is unaffordable on one consumer GPU. But AlphaZero does not roll
out to terminal states -- it evaluates leaves with a VALUE FUNCTION, and that is
what makes the search affordable. So the whole proposal reduces to one
question: can something cheap predict the expensive score?

This tests the cheapest possible candidate first, before building anything
learned: ProteinMPNN's OWN log-likelihood of the sequence given the backbone.
It costs ~0.1 s, it requires no training, and the model is already loaded as
the policy. If it correlates with the high-accuracy fold score, the value
function is free. If it does not, a learned value net is needed and the 180
examples available here are thin for that -- which is itself worth knowing
before promising the approach in a proposal.

Design of the test
------------------
Correlation is Spearman, rank-based, because the search only needs the ORDERING
of candidates, not calibrated values.

Reported three ways, and the differences matter:

  pooled          all 180 designs. Flattered by between-target differences: a
                  metric that only knows "RBD is hard" scores well here while
                  being useless for choosing among candidates.
  within-target   pooled across backbones inside one target.
  within-backbone THE relevant number. MCTS chooses between sequences on the
                  same backbone, so this is the task the value function
                  actually faces. It is also the hardest, and this project has
                  already measured that within-backbone variance dominates --
                  which is exactly why the easier numbers must not be quoted
                  in its place.

Held-out target reporting is included because a value function that works only
on targets it was tuned on is not a value function.

Sign convention: ProteinMPNN reports negative log-likelihood, so LOWER is a
sequence the model likes more. A useful value function should therefore
correlate NEGATIVELY with a fold score where higher is better. The sign is
preserved in the output rather than silently flipped.
"""
import argparse
import copy
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
import torch

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
MPNN_DIR = os.path.join(ROOT, "external_tools/ProteinMPNN")
sys.path.insert(0, MPNN_DIR)

from protein_mpnn_utils import ProteinMPNN, tied_featurize, _scores  # noqa: E402

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
AA_IDX = {a: i for i, a in enumerate(ALPHABET)}
VS = os.path.join(KIT, "results/variance_study")

SOURCES = {
    "il7ra": os.path.join(KIT, "results/il7ra_real_pilot50/proteinmpnn"),
    "pdl1":  os.path.join(ROOT, "real_test_results/pdl1_pilot50/proteinmpnn"),
    "rbd":   os.path.join(ROOT, "real_test_results/rbd_200_batch/batch_0/proteinmpnn"),
}


def mn(x):
    if isinstance(x, dict):
        v = [float(i) for i in x.values()]
        return min(v) if v else np.nan
    return float(x) if x is not None else np.nan


def spearman(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    n = len(x)
    if n < 5:
        return np.nan, np.nan, n
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    if rx.std() == 0 or ry.std() == 0:
        return np.nan, np.nan, n
    rho = float(np.corrcoef(rx, ry)[0, 1])
    if abs(rho) >= 1:
        return rho, 0.0, n
    from math import erfc, sqrt
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    return rho, float(erfc(abs(t) / sqrt(2))), n


def load_model(weights, device):
    ckpt = torch.load(weights, map_location=device)
    m = ProteinMPNN(num_letters=21, node_features=128, edge_features=128,
                    hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                    augment_eps=0.0, k_neighbors=48)   # noise off: see dpo_finetune.py
    m.load_state_dict(ckpt["model_state_dict"])
    m.to(device).eval()
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights",
                    default=os.path.join(MPNN_DIR, "vanilla_model_weights/v_48_020.pt"))
    ap.add_argument("--orders", type=int, default=5,
                    help="decoding orders averaged per likelihood")
    ap.add_argument("--out", default="results/rl/value_probe.csv")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.chdir(KIT)
    torch.manual_seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = load_model(args.weights, device)

    rows = []
    for t, src in SOURCES.items():
        man_p = os.path.join(VS, t, "manifest.csv")
        if not os.path.exists(man_p):
            continue
        man = pd.read_csv(man_p)
        # high-accuracy fold scores
        scores = {}
        for f in glob.glob(os.path.join(VS, t, "folds", "*_scores_rank_001_*.json")):
            stem = os.path.basename(f).split("_unrelaxed_")[0].split("_scores_rank_001_")[0]
            d = json.load(open(f))
            p = np.array(d.get("plddt", [0.0]))
            scores[stem] = {"ipTM": d.get("iptm"), "ipSAE_min": mn(d.get("ipsae")),
                            "pDockQ2_min": mn(d.get("pdockq2")), "pLDDT": float(p.mean())}
        # backbones
        with open(os.path.join(src, "assigned.jsonl")) as fh:
            assigned = json.loads(fh.read().strip().splitlines()[0])
        proteins = {}
        for line in open(os.path.join(src, "parsed.jsonl")):
            pr = json.loads(line)
            if pr["name"] in assigned:
                proteins[pr["name"]] = pr

        cache = {}
        for _, r in man.iterrows():
            bb, smp = str(r["backbone"]), str(r["sample"])
            stem = "%s__%s" % (bb, smp)
            if stem not in scores or bb not in proteins:
                continue
            if bb not in cache:
                cd = {bb: (assigned[bb][0], assigned[bb][1])}
                o = tied_featurize([copy.deepcopy(proteins[bb])], device, cd,
                                   None, None, None, None, None, ca_only=False)
                cache[bb] = (o[0], o[1], o[2], o[4], o[10], o[12], o[5])
            X, S_t, mask, chain_M, chain_M_pos, residue_idx, chain_enc = cache[bb]
            dmask = mask * chain_M * chain_M_pos
            seq = str(r["sequence"])
            idx = torch.tensor([AA_IDX.get(c, 20) for c in seq], device=device,
                               dtype=torch.long)
            S = S_t.clone()
            n = min(len(idx), S.shape[1])
            S[0, :n] = idx[:n]
            with torch.no_grad():
                vals = []
                for _ in range(args.orders):
                    rnd = torch.randn(chain_M.shape, device=device)
                    lp = model(X, S, mask, chain_M * chain_M_pos, residue_idx,
                               chain_enc, rnd)
                    vals.append(float(_scores(S, lp, dmask)))
            rows.append({"target": t, "backbone": bb, "sample": smp,
                         "mpnn_nll": float(np.mean(vals)),
                         "mpnn_nll_sd": float(np.std(vals)),
                         "penalty": float(r["penalty"]), "length": int(r["length"]),
                         **scores[stem]})

    if not rows:
        raise SystemExit("[FAIL] no variance-study designs could be matched")
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)

    print("=" * 78)
    print("VALUE-FUNCTION PROBE -- can ProteinMPNN's own likelihood predict the fold?")
    print("=" * 78)
    print("  designs %d over %d targets, %d backbones"
          % (len(df), df.target.nunique(), df.backbone.nunique()))
    print("  mpnn_nll: LOWER = model likes it more, so a useful value function")
    print("  correlates NEGATIVELY with ipTM / ipSAE / pDockQ2.\n")

    for metric in ("ipSAE_min", "pDockQ2_min", "ipTM"):
        print("-" * 78)
        print("%s" % metric)
        rho, p, n = spearman(df.mpnn_nll.to_numpy(float), df[metric].to_numpy(float))
        print("  pooled            rho %+.3f  p %.3f  n %d   (flattered by target effects)"
              % (rho, p, n))
        for t, g in df.groupby("target"):
            rho, p, n = spearman(g.mpnn_nll.to_numpy(float), g[metric].to_numpy(float))
            print("  within %-8s   rho %+.3f  p %.3f  n %d" % (t, rho, p, n))
        # within-backbone: the task MCTS actually faces
        rs, ws = [], []
        for (_t, _b), g in df.groupby(["target", "backbone"]):
            if len(g) < 3:
                continue
            rho, _p, n = spearman(g.mpnn_nll.to_numpy(float), g[metric].to_numpy(float))
            if not np.isnan(rho):
                rs.append(rho)
                ws.append(n)
        if rs:
            mean_rho = float(np.average(rs, weights=ws))
            print("  WITHIN-BACKBONE   rho %+.3f  (mean over %d backbones)  <-- the real task"
                  % (mean_rho, len(rs)))

    print("\n" + "=" * 78)
    print("READ THIS AS")
    print("=" * 78)
    print("  The within-backbone number is the one that matters: MCTS chooses")
    print("  between sequences on a fixed backbone. |rho| below ~0.2 there means")
    print("  ProteinMPNN's likelihood is NOT a usable free value function and a")
    print("  learned value net is required -- on 180 examples, which is thin.")
    print("\n[OK] %s" % args.out)


if __name__ == "__main__":
    main()
