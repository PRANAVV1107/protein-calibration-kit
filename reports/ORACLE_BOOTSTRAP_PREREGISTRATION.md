# Pre-registration — how often does uniform beat the oracle?

**Written before the analysis was run.** Committed to fix the decision rule
before any bootstrap output existed.

## The claim under test

ALLOCATION_RESULT.md and PREPRINT_DRAFT.md §3.4 lead on:

> An oracle spending every fold on the root with the best true mean reaches
> 0.770 where uniform across ten roots reaches 0.869.

That is **one draw of ten RBD roots on one instrument**. Two facts motivate
re-testing it:

1. The same simulation on the other two targets gives oracle − uniform of
   **+0.001 (PD-L1)** and **+0.008 → +0.001 (IL-7Rα)**. The effect is RBD-only
   in the existing data.
2. In the 50-root breadth pool for RBD, the best-mean root (`binder_22`,
   mean 0.597) **also holds the pool maximum (0.860)** — so an oracle confined
   to it would tie, not lose.

The oracle is capped at its single root's observed maximum. It therefore loses
to uniform **iff the best-mean root does not contain the pool maximum** — an
argmax discordance. The headline may be one favourable draw.

## Hypothesis

H1: Across random 10-root subsets of the 50-root RBD pool, the best-mean root
fails to contain the subset maximum in a **minority** of draws, and the
headline gap (−0.099 at B=160) sits in the **upper tail** of the gap
distribution rather than near its centre.

Declared expectation before running: discordance in **30–50%** of draws, with a
median gap materially smaller than 0.099.

## Method, fixed in advance

- Source: `results/root_breadth_per_root.csv` (50 RBD backbones, 4 draws each).
- Draw 10 roots without replacement, 10,000 replicates, seed 0.
- Per replicate record: (a) whether argmax(mean) == argmax(max); (b) oracle
  value = max of the best-mean root; (c) uniform value = subset max as the
  ceiling uniform approaches at high budget; (d) gap = oracle − uniform.
- Report: discordance rate, median and IQR of gap, and the percentile of
  −0.099 in the gap distribution.
- Ties in per-root mean or max are broken by lowest backbone index, declared
  now because `binder_22` and `binder_23` tie at max 0.860.

## Decision rule, fixed in advance

| outcome | what happens to the paper |
|---|---|
| discordance ≥ 50% **and** median gap ≤ −0.05 | the headline stands as a general effect; report the distribution alongside it |
| discordance 20–50% | headline is **re-stated as conditional**: "when the best-mean root does not contain the best design, which happens in X% of 10-root draws" |
| discordance < 20% **or** −0.099 beyond the 90th percentile | the headline is **withdrawn as a general claim** and reported as one favourable draw; §3.4 leads on the roots-beat-draws result instead |

## Limits, acknowledged in advance

- 4 draws per root is thin; a root's observed max underestimates its true max,
  and unevenly across roots.
- Resampling from one target's 50 backbones. Says nothing about IL-7Rα or
  PD-L1, where the existing data already shows no effect.
- Subsets are drawn from one pool, so replicates are not independent.
- This tests the **ceiling** the oracle and uniform approach, not the full
  budget curve. A budget-resolved version requires re-running `bandit_sim.py`
  per subset.
