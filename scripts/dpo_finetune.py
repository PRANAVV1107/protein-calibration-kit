#!/usr/bin/env python3
"""Direct Preference Optimization fine-tuning of ProteinMPNN.

The argument for this approach
------------------------------
Online RL over this pipeline is unaffordable, and that is a measured result
rather than an assumption. Each rollout needs a reward; the only reward that
predicts measured binding is a high-accuracy ColabFold run at 4-8 minutes; and
the cheap screen fast enough to substitute scores BELOW chance against the 164
labelled Adaptyv designs (ipSAE 0.466, pDockQ2 0.461). A policy trained against
an anti-predictive reward would learn to select against binders.

DPO removes the cost rather than routing around it: the rewards already exist,
as scores on designs already folded. What it needs is preference PAIRS, and the
variance study produced exactly that shape -- several ProteinMPNN draws per
backbone, each folded at high accuracy. Holding the backbone fixed makes the
comparison isolate the sequence choice, which is the same fact as this
project's headline result that within-backbone variance dominates
between-backbone variance.

The loss
--------
    L = -log sigmoid( beta * [ (logp_t(w) - logp_r(w)) - (logp_t(l) - logp_r(l)) ] )

for chosen w and rejected l, where logp_t is the model being trained and logp_r
a frozen copy of the starting weights. The reference term is what keeps this
from collapsing: without it the objective is maximised by assigning all mass to
whatever sequences appear as 'chosen', destroying the model's general
distribution. beta sets how far the policy may move from the reference; low
beta permits more drift and more reward hacking.

Two details specific to ProteinMPNN
-----------------------------------
1. Decoding order is random. ProteinMPNN factorises p(sequence|backbone) over a
   random permutation, so a single evaluation is a noisy estimate of the
   sequence log-likelihood. This averages over --orders permutations, and uses
   THE SAME permutations for the trained and reference models within a pair.
   Sharing them matters more than the count: the loss depends on a difference
   of log-probabilities, so common random numbers cancel most of the variance
   that would otherwise swamp the signal.

2. Only the designed chain is scored. mask * chain_M * chain_M_pos restricts
   the likelihood to the binder; including the target chain would average the
   signal against a constant the model cannot change.

Evaluation must be prospective
------------------------------
The accuracy printed here is on held-out PAIRS, which measures only whether the
model reproduces the preference ordering it was trained on. That is not
evidence the designs are better. The real test is to generate fresh sequences
from the tuned model on unseen backbones, fold them at high accuracy, and
compare against the untuned baseline on those same backbones. --holdout-target
supports the harder version of that, training on two targets and leaving one
untouched.
"""
import argparse
import copy
import json
import os
import random
import sys

import numpy as np
import torch

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
MPNN_DIR = os.path.join(ROOT, "external_tools/ProteinMPNN")
sys.path.insert(0, MPNN_DIR)

from protein_mpnn_utils import ProteinMPNN, tied_featurize, _scores  # noqa: E402

ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
AA_IDX = {a: i for i, a in enumerate(ALPHABET)}

# Where each target's backbones and chain assignments live, and which chain
# ProteinMPNN designed. Taken from variance_study.sh, which generated them.
SOURCES = {
    "il7ra": (os.path.join(KIT, "results/il7ra_real_pilot50/proteinmpnn"), "A"),
    "pdl1":  (os.path.join(ROOT, "real_test_results/pdl1_pilot50/proteinmpnn"), "B"),
    "rbd":   (os.path.join(ROOT, "real_test_results/rbd_200_batch/batch_0/proteinmpnn"), "A"),
}


def load_model(weights, device, hidden=128, layers=3, k=48):
    """Load ProteinMPNN with backbone noise DISABLED.

    augment_eps is training-time coordinate noise: the forward pass perturbs X
    by that much every call. The checkpoint carries 0.2, and inheriting it
    breaks DPO in a way that is easy to miss -- the policy and the frozen
    reference then see different structures on every evaluation, so their
    log-probabilities differ even when the weights are identical. The dry-run
    margin was -6.19 where it must be exactly 0.

    protein_mpnn_run.py sets this to 0.0 for inference for the same reason.
    Zero here makes a likelihood depend only on the weights and the decoding
    order, which is what lets shared random permutations cancel between the two
    models and leave only the signal DPO is meant to train on.
    """
    ckpt = torch.load(weights, map_location=device)
    m = ProteinMPNN(num_letters=21, node_features=hidden, edge_features=hidden,
                    hidden_dim=hidden, num_encoder_layers=layers,
                    num_decoder_layers=layers, augment_eps=0.0, k_neighbors=k)
    m.load_state_dict(ckpt["model_state_dict"])
    return m.to(device)


def load_backbones(targets):
    """name -> (protein dict, chain_id_dict entry, target)."""
    out = {}
    for t in targets:
        d, _chain = SOURCES[t]
        pj, aj = os.path.join(d, "parsed.jsonl"), os.path.join(d, "assigned.jsonl")
        if not (os.path.exists(pj) and os.path.exists(aj)):
            print("  [warn] %s: parsed/assigned jsonl missing, skipping" % t)
            continue
        assigned = json.load(open(aj)) if aj.endswith("jsonl") else {}
        # assigned.jsonl here is a single JSON object, not line-delimited
        with open(aj) as fh:
            assigned = json.loads(fh.read().strip().splitlines()[0])
        for line in open(pj):
            p = json.loads(line)
            if p["name"] in assigned:
                out[(t, p["name"])] = (p, assigned[p["name"]])
    return out


def seq_logprob(model, feats, seq, randns, designed_mask):
    """Total log p(seq | backbone), averaged over the given decoding orders."""
    X, S_template, mask, chain_M, chain_M_pos, residue_idx, chain_enc = feats
    idx = torch.tensor([AA_IDX.get(c, 20) for c in seq],
                       device=X.device, dtype=torch.long)
    S = S_template.clone()
    n = min(len(idx), S.shape[1])
    S[0, :n] = idx[:n]
    tot = 0.0
    for r in randns:
        log_probs = model(X, S, mask, chain_M * chain_M_pos, residue_idx, chain_enc, r)
        # _scores returns MEAN negative log-likelihood per scored residue, so
        # the total log-likelihood is -(mean NLL) * (number of residues).
        nll = _scores(S, log_probs, designed_mask)
        tot = tot + -(nll * designed_mask.sum(dim=-1))
    return tot / len(randns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default="results/rl/preference_pairs.jsonl")
    ap.add_argument("--weights",
                    default=os.path.join(MPNN_DIR, "vanilla_model_weights/v_48_020.pt"))
    ap.add_argument("--out-dir", default="results/rl/dpo")
    ap.add_argument("--beta", type=float, default=0.1,
                    help="lower = more drift from reference allowed")
    ap.add_argument("--lr", type=float, default=1e-5)
    ap.add_argument("--epochs", type=int, default=2)
    ap.add_argument("--orders", type=int, default=3,
                    help="decoding orders averaged per likelihood estimate")
    ap.add_argument("--val-frac", type=float, default=0.15)
    ap.add_argument("--holdout-target", default=None,
                    help="train without this target entirely")
    ap.add_argument("--max-pairs", type=int, default=0, help="0 = all")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true",
                    help="load everything, score one pair, then stop")
    args = ap.parse_args()

    os.chdir(KIT)
    torch.manual_seed(args.seed)
    random.seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    pairs = [json.loads(l) for l in open(args.pairs)]
    if args.holdout_target:
        pairs = [p for p in pairs if p["target"] != args.holdout_target]
    random.shuffle(pairs)
    if args.max_pairs:
        pairs = pairs[:args.max_pairs]
    targets = sorted({p["target"] for p in pairs})

    print("=" * 76)
    print("DPO FINE-TUNING of ProteinMPNN")
    print("=" * 76)
    print("  device   : %s" % device)
    print("  pairs    : %d   targets: %s%s"
          % (len(pairs), ", ".join(targets),
             "   (held out: %s)" % args.holdout_target if args.holdout_target else ""))
    print("  beta %.3g   lr %.3g   epochs %d   orders %d"
          % (args.beta, args.lr, args.epochs, args.orders))

    bb = load_backbones(targets)
    print("  backbones loaded: %d" % len(bb))
    pairs = [p for p in pairs if (p["target"], p["backbone"]) in bb]
    print("  pairs with a resolvable backbone: %d" % len(pairs))
    if not pairs:
        raise SystemExit("[FAIL] no pairs could be matched to a backbone structure")

    policy = load_model(args.weights, device)
    ref = load_model(args.weights, device)
    ref.eval()
    for q in ref.parameters():
        q.requires_grad_(False)
    # eval() on the POLICY too, which disables dropout without touching
    # gradients -- eval only changes dropout/batchnorm behaviour, autograd
    # still flows. ProteinMPNN carries dropout 0.1 on every attention layer, so
    # leaving the policy in train() while the reference sits in eval() makes
    # their log-probabilities differ by pure noise: the dry-run margin was 4.26
    # where identical weights must give exactly 0. DPO's objective is already
    # a small difference of large log-likelihoods, and dropout noise of that
    # size would dominate the gradient.
    policy.eval()
    opt = torch.optim.AdamW(policy.parameters(), lr=args.lr)

    # featurise each backbone once
    cache = {}
    for key, (protein, chains) in bb.items():
        designed, fixed = chains[0], chains[1]
        cd = {protein["name"]: (designed, fixed)}
        out = tied_featurize([copy.deepcopy(protein)], device, cd,
                             None, None, None, None, None, ca_only=False)
        X, S, mask, _l, chain_M, chain_enc = out[0], out[1], out[2], out[3], out[4], out[5]
        chain_M_pos, residue_idx = out[10], out[12]
        feats = (X, S, mask, chain_M, chain_M_pos, residue_idx, chain_enc)
        cache[key] = (feats, mask * chain_M * chain_M_pos)

    n_val = max(1, int(len(pairs) * args.val_frac))
    val, train = pairs[:n_val], pairs[n_val:]
    print("  train %d   val %d\n" % (len(train), len(val)))

    def pair_terms(p, model_a, model_b):
        feats, dmask = cache[(p["target"], p["backbone"])]
        randns = [torch.randn(feats[3].shape, device=device) for _ in range(args.orders)]
        la_w = seq_logprob(model_a, feats, p["chosen"], randns, dmask)
        la_l = seq_logprob(model_a, feats, p["rejected"], randns, dmask)
        with torch.no_grad():
            lb_w = seq_logprob(model_b, feats, p["chosen"], randns, dmask)
            lb_l = seq_logprob(model_b, feats, p["rejected"], randns, dmask)
        return (la_w - lb_w) - (la_l - lb_l)

    if args.dry_run:
        p = train[0]
        with torch.no_grad():
            d = pair_terms(p, policy, ref)
        print("[dry-run] %s/%s  margin=%.4f (should be ~0 before training)"
              % (p["target"], p["backbone"], float(d)))
        print("[dry-run] plumbing OK")
        return

    os.makedirs(args.out_dir, exist_ok=True)
    hist = []
    for ep in range(1, args.epochs + 1):
        random.shuffle(train)
        run_loss, run_acc, seen = 0.0, 0, 0
        for i, p in enumerate(train, 1):
            delta = pair_terms(p, policy, ref)
            loss = -torch.nn.functional.logsigmoid(args.beta * delta).mean()
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(policy.parameters(), 1.0)
            opt.step()
            run_loss += float(loss)
            run_acc += int(float(delta) > 0)
            seen += 1
            if i % 50 == 0:
                print("  ep%d %4d/%d  loss %.4f  train-pref-acc %.3f"
                      % (ep, i, len(train), run_loss / seen, run_acc / seen), flush=True)

        with torch.no_grad():
            v_acc = sum(int(float(pair_terms(p, policy, ref)) > 0) for p in val) / len(val)
        print("  [epoch %d] train loss %.4f  train acc %.3f  VAL pref-acc %.3f"
              % (ep, run_loss / max(seen, 1), run_acc / max(seen, 1), v_acc))
        hist.append({"epoch": ep, "train_loss": run_loss / max(seen, 1),
                     "train_acc": run_acc / max(seen, 1), "val_pref_acc": v_acc})
        torch.save({"model_state_dict": policy.state_dict(),
                    "num_edges": 48, "noise_level": 0.2},
                   os.path.join(args.out_dir, "dpo_epoch%d.pt" % ep))

    with open(os.path.join(args.out_dir, "history.json"), "w") as f:
        json.dump({"args": vars(args), "history": hist}, f, indent=2)

    print("\n  Val pref-acc measures only whether the model reproduces the")
    print("  training preference ordering. It is NOT evidence the designs are")
    print("  better. Generate fresh sequences, fold at high accuracy, and")
    print("  compare against the untuned baseline on the same backbones.")
    print("\n[OK] %s" % args.out_dir)


if __name__ == "__main__":
    main()
