# Result — guided search over denoising trajectories

Run completed 2026-09-17. Evaluated against `MCTS_PREREGISTRATION.md`, written
before any result existed. Nothing below deviates from the criteria fixed there.

## Outcome

| pre-registered outcome | result |
|---|---|
| **primary** — best ipTM > 0.87 | **0.860 — NOT cleared** |
| improved on root best-of-6 (0.52) | yes, 0.860 |
| **secondary** — children beat parent above chance | **26.2%, CI [20.8, 31.6] — below chance** |

**1012 folds, 253 backbones, 21 hours, one RTX 3060.**

The headline comparison, stated as the pre-registration requires:

    guided search   0.860  in 1012 folds
    blind sampling  0.870  in   60 folds

Guided search did not beat blind sampling, at roughly **17x the compute**. The
0.010 difference is far inside the resolution at which these scores mean
anything; the honest reading is that the two are indistinguishable, and the
compute ratio is the finding.

## Why it failed — the mechanism

Three independent observations from this run converge on one explanation.

**Parent score does not predict child quality.** Every improvement came from
expanding a *weak* node. The lineage that produced the best design reads
0.150 -> 0.140 -> 0.740 -> 0.840 -> 0.860: two of its five steps were from
parents scoring below 0.16. Meanwhile expanding the best-known node produced
nothing above 0.27 on three separate occasions (nodes 85, 99, 105).

**Depth carries no information.** With the full sample the profile is flat from
depth 2 onward — means 0.306 to 0.416 across depths 2-6, within noise of each
other. Two earlier readings of this run claimed a peak at depth 3 and then a
decline past it; both were artifacts of n=4 and n=8 and are withdrawn.

**Children score below their parents at 26.2%.** This metric is confounded --
UCT deliberately expands high-scoring nodes, so regression to the mean predicts
below 50% on its own. But combined with the two observations above, the simplest
account is that `partial_T` neighbours are drawn from a distribution largely
independent of where the re-noising started.

If a backbone's score does not predict its neighbours' scores, then the value
estimates UCT propagates are noise, and tree search degenerates into stratified
random sampling with bookkeeping. That is consistent with what the run produced:
six designs at 0.81+ from 253 backbones (2.4%), against blind sampling's one at
0.87 from 60 (1.7%) — the same band at a similar rate, reached far more
expensively.

## Convergence with the rest of the project

This is the third independent result pointing the same way:

- `value_probe.py` — ProteinMPNN's own likelihood cannot rank sequences within a
  fixed backbone (within-backbone rho about -0.06), so there is no cheap value
  function
- `dpo_finetune.py` / `eval_dpo.py` — preference tuning reaches 0.844 in
  distribution and 0.490 on a held-out target, i.e. chance
- this run — no cheap signal predicts which neighbour will fold well

Nothing available cheaply predicts design quality. That is why search, learned
preference, and heuristic value all fail here for the same underlying reason.

## What this does not establish

- **n = 1.** One root, one target, one seed, one hyperparameter setting. This
  cannot separate "search does not work" from "this configuration did not work".
- **Baselines are historical**, folded at 3 models / 8 recycles against this
  run's 1 model / 3 recycles. That asymmetry was checked directly and is small:
  on 85 designs folded both ways, Spearman 0.971 and identical maxima, so it
  does not explain the outcome.
- **No wet-lab validation.** ipTM is a prediction.
- A different `partial_T`, a deeper tree, a learned value function, or a
  different target could all change this. Proteina-Complexa (arXiv:2603.27950)
  reports MCTS beating best-of-N on hard targets with far larger budgets.

## Permitted claim

> In a single search run on one target, guided best-first search over denoising
> trajectories reached ipTM 0.860 within 1012 folds, against 0.870 for blind
> sampling at 60 folds. Search did not outperform random sampling at ~17x the
> compute. One root, one target, no wet-lab validation.
