# Run 2 — does valuing nodes by `mean` instead of `max` fix the search?

> **Superseded in part, 2026-09-21.** Every score in this report is on the
> hybrid instrument (donated alignment, 1 model / 3 recycles); the 0.870 bar is
> not. Re-scored at full MSA, run 2's best design is 0.820 and the best across
> all three runs is 0.870 — **a tie with the bar, not a loss**. See
> [RESCORE_RESULT.md](RESCORE_RESULT.md). The within-run comparisons below
> (run 1 vs run 2, children-beating-parent, depth profiles) are unaffected,
> since both runs used the same instrument.


**Status: interim, read at matched budget.** Run 2 is still running. This report
is written at the point where its fold count matches run 1's, because that is
the only point at which the two are comparable. See *Budget accounting* below
for why the run did not stop there on its own.

**Verdict: the diagnosis holds, the fix does not rescue the search.** Valuing
nodes by the mean makes node scores behave like estimates — depth now carries
information and children beat parents far more often — but the best design found
is no better. Higher floors, same ceiling.

---

## Setup

Identical to the pre-registered run except for one flag:

    scripts/mcts_design.py --root-pdb .../binder_47.pdb --target rbd \
      --children 4 --seqs 4 --partial-t 10 --max-depth 6 --metric ipTM \
      --node-value mean

Run 1 used `--node-value max`. Same root, same target, same partial-diffusion
depth, same 4 sequences per node.

---

## Result

| | run 1 (`max`) | run 2 (`mean`) | bar |
|---|---|---|---|
| real folds | 1012 | 955 | 60 |
| nodes evaluated | 253 | 237 | — |
| **best single sequence** | **0.860** | **0.820** | **0.870** |
| best node value | — | 0.427 | — |
| children beating their own parent | 26% | **39.8%** [33.6, 46.1] | 50% = chance |

Neither run cleared the pre-registered bar of **0.870**, which is what blind
best-of-N sampling reached in 60 folds. Run 2 spent roughly 16× the compute of
that baseline and came in below it.

**Own score by depth** (run 1's profile was flat):

| depth | n | mean | max |
|---|---|---|---|
| 0 | 1 | 0.117 | 0.117 |
| 1 | 4 | 0.144 | 0.170 |
| 2 | 16 | 0.200 | 0.358 |
| 3 | 36 | 0.186 | 0.325 |
| 4 | 64 | 0.200 | 0.355 |
| 5 | 68 | 0.192 | 0.398 |
| 6 | 48 | **0.220** | **0.427** |

Mean own score rises monotonically from depth 1 to depth 6 apart from one dip at
depth 3, and the best node value climbs 0.170 → 0.427. Under `max` valuation
this profile was flat. **The search is now moving somewhere**; it is the
selection statistic that changed, and nothing else did.

**Where the best design actually came from:**

    node   1   0.120
    node   2   0.260
    node   8   0.780     <- 32 folds
    node  20   0.790
    node 192   0.820     <- 768 folds

The search found 95% of its final quality in its first 32 folds, then bought
**+0.040 for the next ~900**. Run 1 has the same shape. This is the
diminishing-returns curve measured separately over 2,104 pooled RBD folds, seen
from inside a single run.

---

## What this does and does not establish

**Confirmed.** The tail-statistic diagnosis in [NOISE_RESULT.md](NOISE_RESULT.md)
was correct about node *values*. Switching to `mean` produced exactly the two
changes it predicts: parent score began to predict child quality (26% → 40% of
children beat their parent) and a depth signal appeared where there had been
none. That is a real, measured improvement in the search's internal machinery.

**Not established.** The paper's implied claim — that the mean/max correction is
a "free fix" for the search — does not survive. Estimable node values did not
produce a better design. §3.3 of the preprint must say this: the correction
fixes value estimation and value estimation was not the binding constraint.

**Still confounded.** Two limitations apply to both runs and neither was fixed
before run 2 launched:

1. **The UCT selection was wrong, in two independent ways.** Both are now fixed
   in `scripts/mcts_design.py`; neither fix has been run.

   *Wrong comparison set.* `uct()` took a global argmax over every expandable
   node in the tree, scoring each against **its own parent's** visit count.
   Because expansion evaluates a node exactly once, every candidate sits at
   `n = 1` forever, so the exploration term collapsed to
   `sqrt(log(parent_n + 1))` — a function of the parent alone. A node hanging
   off a heavily-visited parent outranked a better node hanging off a lightly
   visited one, which inverts the intent. Fixed by descending the tree, so the
   candidates at each step are siblings and share a reference count.

   *Unscaled exploration constant.* This one survives the first fix and is
   arguably the larger problem. UCT's constant assumes rewards in [0, 1]; ipTM
   node values here span ~0.10–0.45. Candidates for expansion are always
   unexpanded leaves, so `n.n` was always exactly 1 and the bonus reduced to
   `0.6·sqrt(log(parent_n + 1))` — 0.50 at `parent_n` = 1, rising to 1.41 at
   250. That is a span of **0.91 against a value spread of 0.31**, so what
   decided selection was *having an already heavily-visited parent*, not being
   good.

   **This is a rich-get-richer loop, and it is measurable.** Both runs funnelled
   **92–94% of their nodes into two of the root's four branches** (run 1:
   50/42/7/2; run 2: 51/43/4/2) and drove to maximum depth, while node value
   only ever broke ties between siblings. It explains node 192 — passed over for
   hanging off a starved branch, not for being deep — and it explains run 2's
   depth signal arriving with no gain in the best design: the search drilled to
   the bottom of two branches and exhausted them while the reachable set stayed
   narrow. Fixed by normalising node values to their observed range, which keeps
   `c_uct` meaning what it is documented to mean across targets and metrics.

   *(An earlier version of this section called both runs "breadth-first". That
   was wrong in the opposite direction and is corrected above; the subtree
   shares are what settled it.)*

   A regression test covering both — reproducing run 2's tree, where node 192 at
   0.427 was passed over roughly ten times — now selects it. Whatever the two
   runs measure, it is not a correctly implemented UCT search, and this belongs
   ahead of the mean/max result in any writeup.
2. **n = 1.** One root, one target, one seed, one hyperparameter setting.

---

## Budget accounting — a unit bug

`--budget` is documented as "total folds to spend" and is compared against
`self.folds_used`, but `folds_used` is incremented once per *node*, while each
node costs `--seqs` = 4 folds. The budget is therefore in node units, 4× coarser
than its name.

Run 1 was launched with a budget that happened to land at 1012 real folds. Run 2
was launched with `--budget 1000`, which means **4000 real folds** — roughly 40
more GPU-hours, and an unmatched comparison against run 1. It is at 948 folds
now, which is why this report is written now.

Given the measured returns to compute (expected best-of-N gains ~+0.01 per
doubling past N≈1000 on this target), the remaining 3000 folds are worth
approximately +0.02 in expectation. That is not worth 40 hours, and it destroys
budget matching. **Recommendation: stop run 2 here and treat 948 folds as the
result.** Fix the `--budget` unit before any run 3.

---

## If there is a run 3

In priority order, and only one variable at a time:

1. ~~Fix the UCT selection term~~ — **done**, both halves, with a test.
2. ~~Fix the `--budget` unit~~ — **done**. It now counts folds, as documented.
   Budgets written for the old behaviour must be multiplied by `--seqs`.
3. Run on a **low-skew** target (PD-L1, skew −1.98; IL-7Rα, −1.83) rather than
   RBD (+3.38). The noise measurement says node evaluation is ~7× cleaner there.
   This is the only remaining variable worth changing, and it should change
   alone.

A fourth run on RBD at a larger budget is the one thing not worth doing. The
scaling curve and both runs agree that this root's reachable distribution tops
out near 0.86 and more sampling does not breach it.
