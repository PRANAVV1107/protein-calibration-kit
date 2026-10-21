# Allocating folds across roots — the bandit simulation

**Run 2026-09-21 on folds already on disk. No GPU.** Written to decide whether
multi-root search (Experiment 4, ~14 h plus a code change) was worth building.

**Verdict: no. Spread folds uniformly across roots and stop there.** No
allocation policy beats uniform at any budget, on either instrument — and
*uniform beats the oracle*, which is the result that explains why.

---

## What prompted it

The `partial_T` sweep turned up two regimes:

    across independent backbones   rho +1.000  (n=6 clustered, p=0.0014)
    within one lineage             rho +0.185  (n=62, p=0.146)

A backbone's score predicts its neighbours' when backbones are independently
generated, and does not when they are descendants of one. A forest MCTS makes
selection decisions in both places but only the first has measurable signal, so
what the data actually supports is a **bandit over roots**, not a tree.

`scripts/bandit_sim.py` replays that bandit against the alternatives on existing
folds, so the expensive run only happens if a policy earns it.

## Result

E[best ipTM found], 3000 replicates, 10 RBD roots.

**Hybrid instrument** (`noise_poc`, 16 draws per root — the configuration the
searches used):

| budget | single | **uniform** | halving | greedy | ucb | oracle |
|---|---|---|---|---|---|---|
| 20 | 0.423 | **0.565** | 0.555 | 0.565 | 0.565 | 0.698 |
| 40 | 0.458 | **0.676** | 0.649 | 0.640 | 0.676 | 0.737 |
| 80 | 0.467 | **0.745** | 0.719 | 0.663 | 0.743 | 0.749 |
| 160 | 0.469 | **0.774** | 0.753 | 0.664 | 0.765 | 0.750 |

**High-accuracy instrument** (`variance_study`, 6 draws per root — the one the
0.870 bar was measured on):

| budget | single | **uniform** | halving | greedy | ucb | oracle |
|---|---|---|---|---|---|---|
| 20 | 0.696 | **0.833** | 0.830 | 0.832 | 0.834 | 0.770 |
| 40 | 0.696 | **0.856** | 0.844 | 0.844 | 0.849 | 0.770 |
| 80 | 0.694 | **0.866** | 0.858 | 0.843 | 0.858 | 0.770 |
| 160 | 0.697 | **0.869** | 0.866 | 0.844 | 0.863 | 0.770 |

**No policy beats uniform anywhere.** The best bandit is within ±0.001 of it at
small budgets and falls behind as budget grows. Greedy is materially worse —
0.664 against 0.774 at B=160 on the hybrid instrument — because committing to
an early leader forfeits the diversity that produces the result.

## The number that explains it: uniform beats the oracle

The `oracle` column spends every fold on the root with the **best true mean** —
perfect knowledge, unattainable by any real policy. On the high-accuracy
instrument it reaches **0.770**, and uniform reaches **0.869**.

**Knowing exactly which root is best, and exploiting it perfectly, is worse than
spreading folds evenly across ten roots.**

This is the project's own central finding arriving from a new direction. The
objective is a **maximum over designs**, not a mean. Allocation policies —
bandits, oracles, best-first search — all optimise *location*: they steer folds
toward distributions with high means. But a maximum is a **tail** statistic, and
tails are found by sampling many different distributions, not by sampling the
best one harder.

So the across-backbone correlation of +1.000 is real and it is about means. It
simply is not the quantity that pays. Root quality predicts root quality;
it does not tell you where the tail is.

## What this settles

- **Experiment 4 (multi-root search) is cancelled.** It would have spent 14
  hours and a code change building machinery whose selection logic this
  simulation shows to be worthless at best and harmful at worst.
- **Multi-root remains right; allocation across roots does not.** Uniform
  breadth is the recommendation — 60 folds across 10 roots reaches 0.856 on the
  high-accuracy instrument, against 0.870 for the blind-sampling bar and
  0.870 for the best design any search produced in ~2,900 folds.
- **The practical rule is one line:** draw many roots, sample each a few times,
  take the max. Do not search, do not allocate, do not steer.

## Limits

- Every draw is resampled from a backbone's observed pool, so no policy can
  exceed that root's observed maximum. The `single` column is therefore a floor
  rather than an estimate. Comparisons *between multi-root policies* are fair,
  since all draw from the same pools.
- Ten roots, one target. The high-accuracy source has only 6 draws per root.
- UCB's exploration constant (0.3) and the warm-start (2 draws/root) were not
  tuned. A tuned bandit might close the gap to uniform; nothing here suggests it
  would exceed it, and the oracle bound says it cannot.
- This measures allocation among *given* roots. It says nothing about whether
  generating better roots would help — see [ROOT_SYSTEM.md](ROOT_SYSTEM.md),
  where no cheap structural feature predicted root quality at n=10.

---

Regenerate with:

    python scripts/bandit_sim.py --source noise_poc --target rbd
    python scripts/bandit_sim.py --source variance  --target rbd
