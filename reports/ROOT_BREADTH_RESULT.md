# Root-breadth sweep — does E[best] climb faster with roots than with draws?

**Complete 2026-09-22. 50 RBD backbones × 4 ProteinMPNN draws = 200 folds**,
full MSA at 1 model / 3 recycles. No errors, no rate-limiting.

**Verdict: breadth beats depth at every budget tested, and the ceiling is set by
the best root in the pool rather than by how many roots you draw.**

---

## Why this was run

[ALLOCATION_RESULT.md](ALLOCATION_RESULT.md) found that no allocation policy
beats uniform, and that spending every fold on the root with the best *true
mean* (0.770) loses to spreading folds across ten roots (0.869). The
tail-statistic argument explains it: the objective is a maximum, and maxima come
from sampling many distributions rather than the best one harder.

That argument predicts something the ten-root data could not test — **E[best]
should keep climbing with the number of roots**. Uniform flattened at 0.869 there
only because there were ten distributions available. Fifty can distinguish "the
recommendation has a ceiling" from "we ran out of roots".

## 1. The prediction, at matched budget

E[best ipTM] for a budget of B folds spent as R roots × K draws, 3,000
replicates, roots drawn at random:

| budget | spend on **roots** | spend on **draws** | breadth gain |
|---|---|---|---|
| 20 | R=10, K=2 → **0.521** | R=5, K=4 → 0.446 | **+0.075** |
| 50 | R=50, K=1 → **0.763** | R=25, K=2 → 0.719 | **+0.044** |
| 100 | R=50, K=2 → **0.835** | R=25, K=4 → 0.786 | **+0.049** |

**Breadth wins at every budget.** The full surface:

| | K=1 | K=2 | K=3 | K=4 |
|---|---|---|---|---|
| R=5 | 0.337 | 0.398 | 0.429 | 0.453 |
| R=10 | 0.440 | 0.517 | 0.562 | 0.594 |
| R=25 | 0.629 | 0.720 | 0.764 | 0.783 |
| R=50 | **0.765** | **0.836** | **0.852** | **0.856** |

Reading down a column moves far more than reading across a row. Going from 5
roots to 50 at K=1 gains **+0.428**; going from 1 draw to 4 at R=5 gains
**+0.116**. Roots are worth roughly four times as much as draws, fold for fold.

## 2. Where the ceiling actually is

The R axis does **not** flatten at low K — at K=1 it is still climbing steeply
(0.337 → 0.440 → 0.629 → 0.765, the last step worth +0.136). At K=4 it does
saturate, at **0.856**.

But that saturation is not a property of root count. **The best root in the pool
has a maximum of 0.860**, and E[best] cannot exceed the best value any root can
produce. The curve is converging on its own sample maximum.

**So the ceiling is set by the best root available, not by how many roots you
sample.** Adding roots helps until you have found the good one; after that only
a better *generator* helps. This is the cleanest statement of the tail argument
the project has produced, and it is the practical form of the recommendation:
draw roots until the best-so-far stops improving, then stop — more draws on
what you have will not rescue a pool with no good root in it.

## 3. The distribution, n = 50

    best-of-4 per root:  min 0.120   median 0.210   max 0.860   sd 0.180
    roots above 0.7: 6%      roots above 0.8: 6%

    top:    binder_23 (0.86), binder_22 (0.86), binder_9 (0.83), binder_28 (0.65)
    bottom: binder_29 (0.14), binder_24 (0.13), binder_36 (0.13), binder_39 (0.12)

**Three roots in fifty carry essentially all the reachable quality.** The
median root reaches 0.210 and the best reaches 0.860 — a factor of four. That
concentration is why breadth wins and why steering toward high-*mean* roots
fails: the payoff lives in a handful of outliers, and the mean does not locate
them.

Note `binder_47` — the root all three MCTS runs descended from — sits in the
**bottom five at 0.13**.

## 4. A cheap feature that finally survives correction

At n = 50, six structural features computed from the PDB alone against per-root
best-of-4:

| feature | ρ | p | |
|---|---|---|---|
| **radius of gyration** | **−0.482** | **0.0004** | ** |
| helix fraction | −0.254 | 0.075 | |
| hotspot contacts | +0.203 | 0.158 | |
| contacts | +0.191 | 0.185 | |
| binder length | +0.148 | 0.304 | |
| contacts per residue | +0.030 | 0.837 | |

`**` survives Bonferroni over six features (threshold p < 0.0083).

**Compact binders score better, and radius of gyration costs milliseconds.**
This is the first cheap predictor of root quality the project has found.

**How it got there matters as much as the result.** Tracking the two leading
candidates as n grew:

| n | contacts | rog |
|---|---|---|
| 10 | **+0.511** "the one to test prospectively" | −0.286 |
| 22 | +0.057 | −0.497 |
| 27 | +0.095 | −0.413 |
| 29 | +0.178 | −0.442 |
| **50** | **+0.191** (dead) | **−0.482** (survives) |

At n = 10 I named `contacts` as the feature worth testing and `rog` as
"nothing". Exactly the reverse held. **The feature that survived is the one that
was stable across sample sizes, not the one that looked best when the sample was
small** — which is the same lesson as every other small-n result in this project.

**It still needs prospective validation.** `rog` was selected by looking at these
data, on one target, on 50 backbones from one generation batch. The correct next
test is to predict quality on a fresh batch before folding it.

---

## What this changes

- **§3.4 of the preprint gains its positive half.** It currently says steering
  fails. It can now say what to do instead, with a curve: spend on roots, not
  draws, worth about 4× fold-for-fold, until the best-so-far stops improving.
- **Recommendation #1 does not need qualifying, but it needs a stopping rule.**
  "Spread folds across many roots" is confirmed; "and stop when the best root
  stops improving, because the ceiling is the best root you drew, not the number
  you drew" is the missing half.
- **A new, free screening step is available** — rank candidate backbones by
  radius of gyration before folding any of them. Pending prospective validation.

## Limits

- **No instrument calibration.** The ten overlap backbones cannot separate
  configuration from sequence set from K (see
  [EXP1_EXP3_LOG.md](EXP1_EXP3_LOG.md), 2026-09-21 18:00). Offset reads −0.447
  with ρ = 0.148, but three variables are confounded. **Nothing here may be
  pooled with the variance study or compared to the 0.870 bar.** All comparisons
  above are internal to this sweep and unaffected. The fix — re-folding the
  variance study's own 60 sequences at this configuration — is ~3.2 h and has
  not been run.
- One target (RBD), one generation batch, 50 backbones, K=4.
- E[best] is estimated by resampling each root's four observed draws, so no
  policy can exceed a root's observed maximum. This compresses the high-K end
  and is why the K=4 column saturates so cleanly; the R-versus-K *comparison* is
  fair because both arms draw from the same pools.
- `rog` was selected post hoc from six features on these data.

---

Regenerate with:

    python scripts/root_breadth_analyze.py --dir results/root_breadth

Per-root table: `results/root_breadth_per_root.csv`.
