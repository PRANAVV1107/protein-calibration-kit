# Experiment 3 — does the search space have locality at any step size?

**Complete 2026-09-21. 24 cells, 480 folds, no errors.** 6 parent backbones ×
4 values of `partial_T`, 4 children each, 4 ProteinMPNN draws per node.

**Verdict on the question asked: no. `rho` does not rise as the diffusion step
shrinks, so there is no step size at which guided search would have worked.**

**Verdict on the question it accidentally answered: the predictive signal lives
*between* independently generated backbones, not *within* a lineage.** That is
the more interesting result and it reframes the whole search programme.

---

## Why this was run

Run 3 measured parent→best-child Spearman **rho = 0.168, p = 0.384 (n = 42)** at
`partial_T = 10`, with every record-setting design descending from a
bottom-decile parent. If expanding a node returns an uncorrelated draw, the tree
has no edges worth following and best-first search *is* random search.

But `partial_T = 10` was never chosen on evidence. This was the last hypothesis
under which the three search runs might have been salvageable: that locality
exists at a shorter step and the hyperparameter was simply past it.

## 1. Locality versus step size

| T | parents | rho | p |
|---|---|---|---|
| 2 | 6 | 0.257 | 0.623 |
| 4 | 6 | 0.600 | 0.208 |
| 10 | 6 | **1.000** | 0.000 |
| 20 | 6 | 0.771 | 0.072 |
| *10 (run 3, within one lineage)* | *42* | *0.168* | *0.384* |

**There is no trend with T, and what trend there is runs the wrong way** —
`rho` is *lowest* at the shortest step and highest at T = 10. A correlation
length would produce the opposite ordering.

At n = 6, Spearman needs |rho| ≥ 0.83 for p < 0.05, so only T = 10 clears it,
and that value is the across-backbone effect described below rather than
anything specific to T = 10.

**The "fixable hyperparameter" hypothesis is closed.** `partial_T = 10` was not
a bad choice; no choice available would have helped.

### The sweep falsified its own power, from the inside

At an interim check T = 10 read **rho = 1.000 on four parents** while run 3 had
read **0.168 on forty-two at the same T in the same pipeline**. A perfect rank
correlation on four points occurs by chance about 4% of the time, so it was an
ordinary small-sample accident — but seeing the same T land at both extremes
calibrates the whole sweep better than any caveat could.

**A properly powered version needs about 30 parents per T** to detect rho = 0.5
at 80% power: roughly 2,400 folds and 40 hours, five times this budget. That is
the number to quote rather than a vague "underpowered".

## 2. The two regimes — what the experiment actually found

Investigating the T = 10 discrepancy turned up a distinction the design had not
anticipated:

| | rho | Pearson | n |
|---|---|---|---|
| **across independent backbones** (clustered, one point per backbone) | **+0.771** (p = 0.072) | **+0.989** (p = 0.0002) | 6 |
| **within one lineage** (run 3) | +0.185 (p = 0.146) | — | 62 |

    binder_24   parent 0.137 -> best_child 0.183
    binder_10   parent 0.142 -> best_child 0.182
    binder_38   parent 0.173 -> best_child 0.245
    binder_3    parent 0.185 -> best_child 0.296
    binder_31   parent 0.196 -> best_child 0.243
    binder_17   parent 0.367 -> best_child 0.655

Near-monotone with a single adjacent swap, and near-perfectly linear.

**This is not restriction of range.** Parent spread is 0.098 in the sweep and
0.087 in run 3 — comparable. The regimes differ in *what varies*, not how much.

**Why.** Within a lineage every node is a re-noised variant of one original
backbone, so all of them share whatever makes that backbone good or bad. The
only thing left to vary is measurement noise, and noise correlates with nothing.
Across independent backbones the variation *is* backbone quality, and that
quality survives re-diffusion.

**A tree search descending within one lineage operates precisely in the regime
where there is no signal to exploit.**

## 3. Does gain depend on parent quality?

Pooled **rho = +0.364, p = 0.080, n = 24**; better-than-median parents gain
+0.144 against worse-than-median +0.057. Per-T: +0.257, −0.314, +0.886, +0.257.

Positive rather than negative — the opposite of run 3's within-tree result,
where top-quartile parents *lost* 0.063. But the n = 24 is pseudo-replicated
(six backbones at four T values; effective n = 6), the pooled p exceeds 0.05,
and a sign flip between adjacent T values is what noise looks like. Read
alongside section 2, it says **good backbones produce good children** — a
statement about backbone quality, not about search.

## What this settles, and what it does not

**Settled.** No step size rescues guided search. The three search runs were not
victims of a bad hyperparameter.

**Settled elsewhere.** The two-regimes finding argued for a bandit over roots.
That was simulated on existing folds in
[ALLOCATION_RESULT.md](ALLOCATION_RESULT.md): **no allocation policy beats
uniform.** (A specific oracle-vs-uniform number cited there at the time was
later pre-registered and re-tested across 10,000 resampled root draws —
[ORACLE_BOOTSTRAP_RESULT.md](ORACLE_BOOTSTRAP_RESULT.md) — and held in only
19.6% of draws, so it is withdrawn as a general claim; the no-policy-beats-
uniform finding itself is untouched.) The across-backbone signal is real, is
about *means*, and still cannot reliably be exploited by an allocation policy
— because the objective is a maximum, and maxima are tail statistics.
**Experiment 4 (multi-root search) is cancelled.**

**Not settled.** Whether locality exists below the detection threshold of six
parents. This bounds the effect; it does not prove it is zero.

## Limits

- 6 parents per T, drawn by a fixed index rule from one batch of 50.
- One target (RBD), hybrid instrument throughout — see
  [RESCORE_RESULT.md](RESCORE_RESULT.md) for what that instrument is worth.
- `best_child` is a max over 4 children of a mean-of-4 while `parent` is a
  single mean-of-4, so the mean gap and the kid>parent rate are biased upward at
  every T by construction. **Only the correlations are reported**; a constant
  upward shift leaves rank correlation alone.
- The same backbone re-evaluated with fresh draws moved by up to **0.335**
  (binder_17: 0.170, 0.392, 0.505). Parent values carry large noise, which
  attenuates correlation rather than creating it.

---

Regenerate with:

    python scripts/partial_t_analyze.py --dir results/rl/partial_t

Raw pairs: `results/partial_t_pairs.csv`.
