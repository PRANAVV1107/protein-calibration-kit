#!/usr/bin/env python3
"""Would allocating folds across roots beat spreading them uniformly?

The question and why it is not "should we build multi-root MCTS"
----------------------------------------------------------------
The partial_T sweep turned up two regimes:

    across independent backbones   rho +1.000  (n=6 clustered, p=0.0014)
    within one lineage             rho +0.185  (n=62, p=0.146)

A forest MCTS makes selection decisions in both places, but only the first has
measurable signal. Strip the tree away and what is left is a bandit: evaluate
several roots, then spend the remaining folds on the ones that look good.

This replays that bandit against the alternatives on folds already on disk, so
the 14-hour forest run only happens if a policy actually beats uniform.

Policies
--------
  single     all folds on one randomly chosen root -- what runs 1-3 did
  uniform    B/R folds on each of R roots -- the breadth baseline
  halving    successive halving: evaluate all, drop the worst half, repeat
  greedy     seed each root with `warm` draws, then all remaining folds on the
             current leader
  ucb        seed each root with `warm` draws, then allocate one at a time by
             UCB1 on the running mean
  oracle     all folds on the root with the best true mean -- an upper bound
             that no policy can reach, included to show the headroom

What this cannot show
---------------------
Every draw is resampled from that backbone's observed pool, so no policy can
exceed a root's observed maximum. Deep single-root allocation is therefore
penalised relative to reality. The comparison between MULTI-root policies is
fair, since they all draw from the same pools; the `single` arm should be read
as a floor, not an estimate.
"""
import argparse
import glob
import json
import math
import os

import numpy as np
import pandas as pd

KIT = os.environ.get("CALIB_KIT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load_pools(source, target):
    """{backbone: array of ipTM draws}."""
    if source == "noise_poc":
        man = pd.read_csv(os.path.join(KIT, "results/rl/noise_poc/manifest.csv")).set_index("tag")
        pools = {}
        for f in glob.glob(os.path.join(KIT, "results/rl/noise_poc/folds",
                                        "*_scores_rank_001_*.json")):
            tag = os.path.basename(f).split("_scores_rank_001_")[0]
            if tag not in man.index or man.loc[tag, "target"] != target:
                continue
            v = json.load(open(f)).get("iptm")
            if v is not None:
                pools.setdefault(man.loc[tag, "backbone"], []).append(float(v))
        return {k: np.asarray(v) for k, v in pools.items()}
    v = pd.read_csv(os.path.join(KIT, "results/variance_study/variance_results.csv"))
    v = v[v.target == target]
    return {b: g.ipTM.dropna().to_numpy() for b, g in v.groupby("backbone")}


def draw(pool, rng, n=1):
    return rng.choice(pool, n, replace=True)


def run_policy(policy, pools, B, rng, warm=2):
    """Return the best score found after spending B folds."""
    names = list(pools)
    R = len(names)
    seen = {k: [] for k in names}
    best = -1.0

    def spend(k, n=1):
        nonlocal best
        d = draw(pools[k], rng, n)
        seen[k].extend(d.tolist())
        best = max(best, float(d.max()))

    if policy == "single":
        spend(names[rng.integers(R)], B)

    elif policy == "oracle":
        k = max(names, key=lambda n_: pools[n_].mean())
        spend(k, B)

    elif policy == "uniform":
        per = max(1, B // R)
        for k in names:
            spend(k, per)

    elif policy == "halving":
        alive, left = list(names), B
        while left > 0 and alive:
            per = max(1, left // (len(alive) * max(1, int(math.log2(len(alive))) + 1)))
            for k in alive:
                if left <= 0:
                    break
                n = min(per, left)
                spend(k, n)
                left -= n
            if len(alive) == 1:
                if left > 0:
                    spend(alive[0], left)
                    left = 0
                break
            alive.sort(key=lambda k_: -np.mean(seen[k_]))
            alive = alive[:max(1, len(alive) // 2)]

    elif policy in ("greedy", "ucb"):
        left = B
        for k in names:
            n = min(warm, left)
            if n <= 0:
                break
            spend(k, n)
            left -= n
        if policy == "greedy":
            if left > 0:
                k = max(names, key=lambda k_: np.mean(seen[k_]) if seen[k_] else -1)
                spend(k, left)
        else:
            t = sum(len(v) for v in seen.values())
            while left > 0:
                k = max(names, key=lambda k_: np.mean(seen[k_])
                        + 0.3 * math.sqrt(math.log(t + 1) / max(1, len(seen[k_]))))
                spend(k, 1)
                left -= 1
                t += 1
    else:
        raise ValueError(policy)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default="noise_poc", choices=["noise_poc", "variance"])
    ap.add_argument("--target", default="rbd")
    ap.add_argument("--reps", type=int, default=3000)
    ap.add_argument("--warm", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    pools = load_pools(a.source, a.target)
    if not pools:
        raise SystemExit("[FAIL] no folds for %s/%s" % (a.source, a.target))
    sizes = [len(v) for v in pools.values()]

    print("=" * 76)
    print("ALLOCATING FOLDS ACROSS ROOTS: does a bandit beat uniform?")
    print("=" * 76)
    print("  source %s, target %s" % (a.source, a.target))
    print("  %d roots, %d-%d draws each, %d folds total\n"
          % (len(pools), min(sizes), max(sizes), sum(sizes)))
    print("  true per-root means: %s\n"
          % "  ".join("%.2f" % pools[k].mean() for k in sorted(pools, key=lambda k_: -pools[k_].mean())))

    policies = ["single", "uniform", "halving", "greedy", "ucb", "oracle"]
    budgets = [20, 40, 80, 160]
    print("  E[best ipTM found], %d replicates, warm-start %d draws/root\n" % (a.reps, a.warm))
    hdr = "  %-8s" % "budget" + "".join("%9s" % p for p in policies)
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    rows = []
    for B in budgets:
        line = "  %-8d" % B
        vals = {}
        for p in policies:
            v = float(np.mean([run_policy(p, pools, B, rng, a.warm) for _ in range(a.reps)]))
            vals[p] = v
            line += "%9.3f" % v
            rows.append({"budget": B, "policy": p, "e_best": v})
        gain = vals["ucb"] - vals["uniform"]
        best_bandit = max(vals["halving"], vals["greedy"], vals["ucb"])
        line += "   best bandit %+0.3f vs uniform" % (best_bandit - vals["uniform"])
        print(line)

    print("""
  READ THIS AS
  ------------
  bandit >> uniform   allocation is worth building; the across-backbone signal
                      is exploitable and a forest MCTS is worth its 14 hours
  bandit ~= uniform   the signal exists but a fixed budget cannot exploit it --
                      spread folds evenly and skip the machinery entirely
  bandit << uniform   early draws mislead; committing to a leader is worse than
                      not choosing

  The `single` column is a floor, not an estimate: resampling caps it at one
  root's observed maximum. `oracle` is the unreachable ceiling.""")

    out = os.path.join(KIT, "results/bandit_sim_%s_%s.csv" % (a.source, a.target))
    pd.DataFrame(rows).to_csv(out, index=False)
    print("\n[OK] %s" % out)


if __name__ == "__main__":
    main()
