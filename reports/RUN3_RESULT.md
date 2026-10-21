# Run 3 — mean-valued nodes with a corrected UCT selection rule

**Complete 2026-09-20, 1012 folds, 253 nodes** — the same counts as run 1,
making the three-way comparison as clean as it could be. Judged against
[RUN3_PREREGISTRATION.md](RUN3_PREREGISTRATION.md), committed before launch.

**Verdict: the selection bug was real and is fixed; fixing it changed how the
search explores and did not change what it finds.**

> **Read [RESCORE_RESULT.md](RESCORE_RESULT.md) alongside this.** Every score
> below is on the hybrid instrument, and the 0.870 bar is not. Re-scored at full
> MSA, run 3's best design is 0.820 and the best across all three runs is 0.870
> — a tie with the bar, not a loss. The primary outcome recorded here was
> assessed on an invalid comparison.

---

## Result

| | run 1 (`max`, broken UCT) | run 2 (`mean`, broken UCT) | run 3 (`mean`, fixed UCT) |
|---|---|---|---|
| folds | 1012 | 955 | **1012** |
| nodes | 253 | 237 | **253** |
| best single sequence (hybrid) | 0.860 | 0.820 | **0.830** |
| best re-scored at full MSA | **0.870** | 0.820 | 0.820 |
| best node value | — | 0.427 | 0.463 |
| children beating their parent | 26% | 40% | 45.6% [39.5, 51.8] |
| max depth reached | 6 | 6 | 6 |

Pre-registered primary: best single sequence against **0.870**. On the hybrid
instrument run 3 reached 0.830 — the **middle outcome**, 0.820–0.869, which is
what I declared in advance that I expected.

**The result was determined early.** 0.820 by fold 52; 0.830 at fold 820. The
remaining 960 folds bought **+0.010**.

    node   6   0.670
    node  13   0.820     <- fold 52
    node 201   0.830     <- fold 820

---

## The two fixes, and what each did

### The UCT selection rule was wrong in two independent ways

Both present in runs 1 and 2, both corrected for run 3.

**Wrong comparison set.** `uct()` took a global argmax over every expandable
node, scoring each against *its own parent's* visit count. Candidates for
expansion are always unexpanded leaves, so `n.n` was permanently 1 and the
exploration term reduced to `sqrt(log(parent_n + 1))` — a function of the parent
alone. Fixed by descending the tree so candidates at each step are siblings
sharing one reference count.

**Unscaled exploration constant.** UCT assumes rewards in [0,1]; ipTM node values
span ~0.10–0.45. At `c_uct` = 0.6 the bonus ran 0.50 to 1.41 against a value
spread of 0.31, so selection was decided by *having a heavily-visited parent*.
Fixed by normalising node values to their observed range.

**Effect: measurable, and in the right direction.** Top-2 subtree share at
matched node counts:

| nodes | run 1 | run 2 | run 3 |
|---|---|---|---|
| 27 | 92% | 62% | 62% |
| 64 | 97% | 78% | 65% |
| 136 | 84% | 90% | 69% |
| 180 | 88% | 92% | **74%** |

Both broken runs converged to 84–94%; run 3 climbed at roughly half the rate and
stayed 14–18 points below. The corrected rule also expanded its own best nodes —
node 75 (0.463, the highest value in any run) was chosen on the third
opportunity, where run 2 passed over its 4th- and 5th-ranked nodes more than
thirty times each.

### The `--budget` unit was wrong

It counted backbones while claiming to count folds, so run 2's `--budget 1000`
meant 4000 real folds. Fixed; run 3's 1000 is genuinely 1000.

---

## The finding that outlasts the fixes

**A node's score carries no information about what expanding it will produce.**
Three independent views, all from run 3 where selection demonstrably worked:

1. **Parent→best-child Spearman ρ = 0.168, p = 0.384** (n = 42). No correlation,
   and fixing selection did not change it — run 1 gave 0.058, run 2 gave 0.162.
2. **Top-quartile parents produce children worse than themselves**: 0.342 →
   0.279. Average parents gain about +0.05. Expanding the best node is the worst
   use of a fold, which inverts the premise of best-first search.
3. **Every record-setting design descended from a bottom-decile parent**:

   | record | best seq | parent value |
   |---|---|---|
   | node 6 | 0.670 | 0.150 |
   | node 13 | 0.820 | 0.133 |
   | node 201 | 0.830 | 0.135 |

   Population mean node value is 0.19. Meanwhile the two highest-valued nodes in
   the tree, 75 (0.463) and 98 (0.440), were both expanded and produced nothing.

In runs 1 and 2 this could be blamed on the broken selection rule steering the
search into bad regions. Run 3 rules that out: **selection chose the best node,
and expansion still produced nothing.** It is a property of the search space —
value does not survive one round of partial re-diffusion at `partial_T = 10`.

`scripts/partial_t_sweep.sh` tests whether any step size preserves it.

---

## Pre-registered secondaries: both failed, for opposite reasons

Recorded because the failures are informative and because the test that *did*
work was chosen after seeing data.

**Secondary 1 (node revisit counts) — retracted mid-run, before the outcome was
known.** The premise was false: `backprop` increments every ancestor, so
internal nodes always have `n > 1` and leaves always have `n = 1` until
expanded. All three runs show 23–25% of nodes above 1 — a function of tree
shape, not of the selection rule. It could not discriminate.

**Secondary 2 (depth distribution) — predicted backwards.** I expected a correct
UCT to push mass to depth 5–6. At matched 136 nodes both broken runs were at
depth 6 while run 3 had not passed depth 4. A correct UCT goes *breadth before
depth*; the broken rule's depth-6 reach was the rich-get-richer loop tunnelling,
not healthy search.

**The replacement (subtree concentration) was chosen after seeing data** and
carries correspondingly less weight than the two registered tests that failed.
It is also n = 1 per arm: runs 1 and 2 differ by 30 points at 27 nodes under the
*identical* rule, so between-seed variance is comparable to the effect.

---

## Claims walked back during this run

All three erred the same way — a trend read from early points and stated more
firmly than the sample supported. Logged in full in [RUN3_LOG.md](RUN3_LOG.md).

| claim | correction |
|---|---|
| "both runs were breadth-first searches wearing UCT's clothing" | measured the opposite: 92–94% of nodes in two of four root branches |
| "only 39% of node-value variance is real", with point predictions | the one pre-committed prediction missed (0.282 predicted, 0.380 observed); the model beats naive by 16%, not decisively |
| "run 3's concentration has gone flat at 65–70%" | it was climbing, just half as fast — 74% by 180 nodes |

**The descriptive measurements in this project have held up; the diagnostic
claims have not.** That belongs in the paper's limitations.

---

## Limits

- n = 1: one root (binder_47), one target, one seed, one hyperparameter setting.
- binder_47 ranks **9th of 10** on measured root quality — see
  [ROOT_SYSTEM.md](ROOT_SYSTEM.md). Every search result is conditional on an
  arbitrary and poor starting point.
- Scores here are hybrid-instrument; see [RESCORE_RESULT.md](RESCORE_RESULT.md).
- No design has been tested at a bench.
