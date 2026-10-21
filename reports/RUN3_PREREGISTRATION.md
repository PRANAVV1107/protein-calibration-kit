# Run 3 — pre-registration

Committed **2026-09-19, before the run started**, in the same form as
[MCTS_PREREGISTRATION.md](MCTS_PREREGISTRATION.md).

---

## Why not PD-L1, as run 2's writeup recommended

The recommendation was to move to a low-skew target because node evaluation
measures ~7× cleaner there. Checking the distributions before spending the GPU
showed that recommendation is confounded:

| target | all folds >0.8 | blind best, 60 folds | headroom to 1.0 | skew |
|---|---|---|---|---|
| PD-L1 | **73.3%** | 0.920 | 0.080 | −1.88 |
| IL-7Rα | 58.3% | 0.890 | 0.110 | −1.74 |
| **SARS-CoV-2 RBD** | **6.7%** | 0.870 | 0.130 | +0.51 |

**PD-L1 and IL-7Rα are left-skewed because they are saturated, not because they
are well-behaved.** Three quarters of PD-L1 draws already clear 0.8 and sixty
blind folds reach 0.92. A search there would post a high number and demonstrate
nothing, because there is nothing left to find.

This is the tension the project has been circling, stated plainly:

> **Search is only worth running on targets where blind sampling struggles, and
> those are exactly the targets where node evaluation is noisy. The targets
> where node evaluation is reliable are the targets that do not need search.**

RBD is the only one of the three with real headroom — 6.7% of blind draws clear
0.8 — so it is the only one where a search win would mean anything. It is also
the target both previous runs used, which makes run 3 a controlled A/B.

## The run

Identical to runs 1 and 2 in every respect — same root (`binder_47`), target,
`--partial-t 10`, `--children 4`, `--seqs 4`, `--max-depth 6`, `--metric ipTM`,
`--node-value mean` — changing only the corrected selection rule:

    --budget 1000        # now genuinely 1000 folds; the old unit was backbones

| | run 1 | run 2 | run 3 |
|---|---|---|---|
| node value | `max` | `mean` | `mean` |
| UCT selection | broken | broken | **fixed** |
| folds | 1012 | 955 | 1000 |

Two fixes are in the corrected rule: descend the tree so candidates are siblings
sharing a reference count, and normalise node values to their observed range so
the exploration bonus (~1.2 at `c_uct` = 0.6) no longer exceeds the entire
spread of the values it is added to. Both are argued in
[RUN2_RESULT.md](RUN2_RESULT.md).

## What counts as what — fixed in advance

**Primary.** Best single ipTM found, against **0.870** — blind sampling across
10 RBD backbones at 60 folds. The same bar runs 1 and 2 were judged against.

- **≥ 0.870** — search clears blind sampling. The selection bug was the binding
  constraint, and both earlier runs understated what this method can do.
- **0.820 – 0.869** — indistinguishable from runs 1 and 2. The selection bug was
  real but not what was holding the search back. This is the outcome the scaling
  curve predicts and the one I expect.
- **< 0.820** — worse than a search with a broken selection rule. Would mean the
  correction is wrong or the root is exhausted; report as such.

**Secondary, decided now so it cannot be chosen afterwards.** These say whether
the fix did what it was argued to do, independent of the primary outcome:

1. **Children beating their own parent.** 26% (run 1) → 40% (run 2) → ? A
   further rise means selection is now concentrating folds on better regions.
2. **Depth reached.** Both earlier runs were argued to be breadth-first in
   effect. A correctly weighted UCT should go deeper on the same budget. Run 2
   reached depth 6 with 48 nodes there; more mass at depth 5–6 is the signal.
3. ~~**Node revisits.** Every node in runs 1 and 2 sat at `n = 1` forever. Under
   the corrected rule the tree should show a visit distribution, not a
   constant.~~ **RETRACTED 2026-09-20 00:40, before the outcome was known.**
   The premise is false and the test cannot discriminate: `backprop` increments
   every ancestor, so internal nodes always have `n > 1` and leaves always have
   `n = 1` until expanded. All three runs show 23–25% of nodes above 1 — a
   function of tree shape, not of the selection rule. **Replaced by subtree
   concentration:** the share of nodes falling in each of the root's four
   branches, which was 50/42/7/2 in run 1 and 51/43/4/2 in run 2. A less skewed
   split is what the corrected rule should produce, compared at matched node
   counts. Reasoning in [RUN3_LOG.md](RUN3_LOG.md).

**Declared before the fact:** I expect the middle outcome. The pooled scaling
curve over 2,104 RBD folds saturates near 0.86 and run 2 reached 0.780 by its
32nd fold, which says the reachable distribution from this root — not the
selection rule — is what bounds the result. Run 3 tests whether that reading is
right. If it clears 0.870 I was wrong about which constraint binds, and that is
the more interesting outcome.
