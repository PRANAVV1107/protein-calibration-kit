# How often does uniform beat the oracle? A bootstrap, pre-registered

**Run 2026-10-29, 10,000 replicates, no GPU.** Judged against
[ORACLE_BOOTSTRAP_PREREGISTRATION.md](ORACLE_BOOTSTRAP_PREREGISTRATION.md),
committed before this analysis was run.

**Verdict: withdrawn as a general claim.** The headline in
[ALLOCATION_RESULT.md](ALLOCATION_RESULT.md) and PREPRINT_DRAFT.md §3.4 — "an
oracle spending every fold on the best-mean root reaches 0.770 where uniform
across ten roots reaches 0.869" — is one favourable draw of ten roots, not a
property of the oracle-vs-uniform comparison in general. In the median draw
from the same 50-root pool, **the oracle does not lose to uniform at all.**

---

## Why this needed testing

The oracle result was measured once, on one draw of ten RBD roots. Two facts
already on record made it worth re-testing before it carried the paper's
headline:

1. The same comparison on the other two targets gives oracle − uniform of
   +0.001 (PD-L1) and +0.001–0.008 (IL-7Rα) — no effect.
2. In the full 50-root RBD pool, the best-mean root (`binder_22`, mean 0.5975)
   **ties the pool maximum** (0.86, shared with `binder_23`). An oracle that
   draws the full pool wins outright. The reported loss is therefore a
   property of which *ten* of the fifty roots were drawn, not of the method.

The oracle is mechanically capped at its own root's observed maximum. It can
only lose to uniform when the best-mean root in a given draw fails to also
contain that draw's best design — an argmax discordance between `mean` and
`max`. The pre-registration fixed the test of how often that happens before
running it.

## Method

10,000 replicates, 10-of-50 roots without replacement, seed 0
(`numpy.random.default_rng`), ties broken by lowest backbone index. Per
replicate: whether argmax(mean) = argmax(max); oracle value = max of the
best-mean root; uniform value = the draw's true max (uniform's ceiling at high
budget); gap = oracle − uniform. Full method and decision thresholds fixed in
advance in the pre-registration; not altered after seeing output.

## Result

    discordance rate              19.63%        (10,000 replicates)
    gap, median                   +0.0000
    gap, IQR                      [+0.0000, +0.0000]
    gap, mean                     -0.0230
    gap, among discordant reps    median -0.080, IQR [-0.220, -0.030]
    headline (-0.099) percentile  92.75th percentile of severity

**In 80.4% of 10-root draws, the oracle and the true-best root are the same
root — gap is exactly zero, and the oracle does not lose at all.** The median
and both quartiles of the full gap distribution are zero. Loss only occurs in
the minority where the best-mean root is not the pool's best root, and even
then the typical loss (median −0.080 among discordant draws) is smaller than
the reported −0.099.

The reported gap of −0.099 sits at the **92.75th percentile of severity** —
only 7.25% of 10,000 resampled draws produce a gap this bad or worse. It is a
tail outcome, not a typical one.

## Against the decision rule, fixed in advance

| condition | threshold | observed | met? |
|---|---|---|---|
| discordance ≥ 50% and median gap ≤ −0.05 | → headline stands | 19.63%, median 0.000 | no |
| discordance 20–50% | → conditional restatement | 19.63% | no (just under) |
| discordance < 20% **or** −0.099 beyond 90th percentile | → **withdrawn** | 19.63% < 20%, **and** 92.75th percentile | **yes, on both independent conditions** |

Both legs of the withdrawal condition hold on their own; this is not a
borderline call decided by one close number.

**The declared expectation was also wrong, in the same direction as the
finding.** Before running, I expected discordance of 30–50%. The measured
19.63% is below even that — the oracle result is more RBD-root-specific, not
less, than the pre-registration anticipated.

## What this changes

- **PREPRINT_DRAFT.md §3.4** no longer leads on the oracle-vs-uniform number.
  The section's lead becomes the roots-beat-draws result (R-vs-K, measured
  across the full 50-root pool and not contradicted by this test), with the
  oracle finding reported as a single-draw illustration and the bootstrap
  statistics given alongside it.
- **The abstract's closing claim** ("spending every fold on the root with the
  best true mean reaches 0.770 ... loses") is rewritten to state what is
  actually general: no allocation *policy* (bandit, greedy, UCB) beats uniform
  at any tested budget on any of three targets. That claim is untouched by this
  result — it was never about the oracle's single-root ceiling, and nothing
  here tests it further (see Limits).
- **README.md's headline hook** is the same oracle number and needs the same
  downgrade: either replaced with the roots-beat-draws number, or kept but
  explicitly labelled as one draw with the bootstrap context attached.
- **ALLOCATION_RESULT.md** keeps its "no policy beats uniform" finding, which
  this result does not touch, but its oracle framing needs the same caveat
  inline rather than left as the section's explaining mechanism.

## What this does not change

- **"No allocation policy beats uniform" stands.** Single, halving, greedy and
  UCB were all ≤ uniform at every tested budget on the original ten-root draw
  for all three targets. This bootstrap only tests the oracle's ceiling against
  uniform's ceiling; it does not re-run the policy simulation across resampled
  subsets. That remains a real gap — see Limits.
- **Roots-beat-draws stands, unexamined by this test.** The R-vs-K result (roots
  worth ~4× single draws fold-for-fold) comes from the full 50-root pool
  directly, not from a 10-root subsample, and this bootstrap does not bear on
  it either way.

## Limits, as declared in advance

- 4 draws per root is thin; a root's observed max underestimates its true
  max, unevenly across roots, so both the discordance rate and the gap
  distribution are measured on a noisy proxy for the true per-root ceiling.
- Resampling from one target's 50 backbones. Says nothing new about IL-7Rα or
  PD-L1, where the oracle effect was already absent in the original data.
- The 10,000 replicates are draws from one fixed pool of 50 roots, so they are
  not independent samples of "RBD root pools in general" — they are a
  bootstrap over one realised set of 50 RFdiffusion outputs.
- This tests the **ceiling** oracle and uniform approach as budget grows, not
  the full budget curve or the bandit policies at finite budget. A
  budget-resolved re-run of `bandit_sim.py` per resampled subset would be
  needed to test whether "no policy beats uniform" is similarly draw-dependent,
  and that has not been done.

---

Regenerate with:

    python scripts/oracle_bootstrap.py

Per-replicate audit trail: `results/oracle_bootstrap_replicates.csv` (10,000
rows: discordant, oracle_val, uniform_val, gap).
