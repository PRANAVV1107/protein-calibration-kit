# Run 3 — running log

Appended as the run proceeds, newest entry at the bottom. Observations only:
anything that turns out to matter graduates to
[RUN3_RESULT.md](RUN3_RESULT.md) at the end, judged against
[RUN3_PREREGISTRATION.md](RUN3_PREREGISTRATION.md).

**What run 3 is.** Mean-valued nodes plus a corrected UCT selection rule, on
RBD, budget-matched to runs 1 and 2. One variable changed from run 2.

| | run 1 | run 2 | run 3 |
|---|---|---|---|
| node value | `max` | `mean` | `mean` |
| UCT selection | broken | broken | **fixed** |
| folds | 1012 | 955 | 1000 (target) |
| best single sequence | 0.860 | 0.820 | — |

Bar: **0.870**, blind sampling at 60 folds. Declared expectation before launch:
the middle outcome, 0.820–0.869.

Watching for, in priority order:

1. **Node revisit counts.** Every node in runs 1 and 2 sat at `n = 1` forever.
   A visit distribution appearing is the most direct evidence the fix bit, and
   it is visible long before the primary outcome is.
2. **Depth distribution** against run 2's 48 nodes at depth 6.
3. **Children beating their own parent**: 26% → 40% → ?
4. Best single sequence against 0.870.

---

## 2026-09-19 23:49 — launched

Driver up, pid 351, root `binder_47`, contig `[E333-526/0 71-71]` (binder 71 aa,
target E 194 aa) — identical to both previous runs. 4 root sequences submitted.

An earlier launch at 10:38 was paused at the user's request two minutes in,
before its first iteration; no `state.json` and 0 folds on disk, so this is a
clean start rather than a resume.

At ~14 h for 1000 folds (run 2's observed rate) this finishes around 14:00 on
2026-09-20. Hourly watcher armed. A background check is waiting on the first
`iter: expanding` line — the first call into the rewritten `uct()`, and the code
path most likely to fail outright.

## 2026-09-19 23:54 — first expansion clean; budget fix confirmed

    [23:53:57] root score 0.130 (folds 4)
    [23:53:57] iter: expanding node 0 (depth 0, mean 0.130, n 1) [folds 4/1000]

The rewritten `uct()` returned a node without throwing — the code path most
likely to fail outright is past. Nothing yet about whether it selects *better*;
at one evaluated node there is nothing to choose between.

**The budget fix is visible in the same line.** The root cost 4 folds and the
counter reads `4/1000`. Under the old accounting this would have read `1/1000`,
and the run would have continued to 4000 real folds. Run 3 is genuinely
budget-matched to run 1 (1012) and run 2 (955) for the first time.

Root scored **0.130**, against run 2's 0.117 on the same backbone with the same
`--seqs 4`. Both are draws of 4 sequences from a stochastic sampler, and the
K=4 disagreement on RBD measured 0.167 in NOISE_RESULT.md, so a gap of 0.013 is
well inside noise and is not evidence of anything. Noting it only so it is not
read as a difference later.

## 2026-09-20 00:40 — advancing; secondary 1 retracted, and "breadth-first" was wrong

**Progress.** 52/1000 folds, 13 nodes, driver alive. Best node value 0.283
(node 6, depth 2); best single sequence **0.670**. Depth 2 reached. Too early
for the primary — run 2 was at 0.780 by fold 32, and a gap this size at 5% of
budget is noise.

Children beating their own parent: 8 of 12 = 66.7%, CI [40.0, 93.3] — not
distinguishable from chance at n=12. Ignore until a few hundred folds.

---

### Secondary 1 cannot discriminate. Retracted.

The pre-registration said "every node in runs 1 and 2 sat at `n = 1` forever"
and made a visit distribution appearing the headline evidence the fix bit. That
was **wrong**, and checking run 3's histogram against the finished runs shows
why:

| run | nodes | leaves with n>1 | internal with n>1 |
|---|---|---|---|
| run 1 (broken UCT) | 253 | **0 of 190** | 63 of 63 |
| run 2 (broken UCT) | 237 | **0 of 178** | 59 of 59 |
| run 3 (fixed UCT) | 13 | **0 of 10** | 3 of 3 |

`backprop` increments a node and all its ancestors, so an internal node has
`n > 1` automatically and a leaf has `n = 1` until it is expanded. The visit
count is a deterministic function of subtree size. It is identical in all three
runs (23–25% above 1) and says nothing whatever about the selection rule. I
should have derived this before pre-registering it.

One thing it does do is *strengthen* the bug diagnosis. Candidates for expansion
are always unexpanded leaves, so `max(1, n.n)` was always exactly 1, and the old
exploration term really did reduce to `sqrt(log(parent_n + 1))` with no
dependence on the candidate at all.

### And the mechanism was the opposite of what I wrote

RUN2_RESULT.md and the preprint both say runs 1 and 2 were "breadth-first
searches wearing UCT's clothing". Measuring it directly:

| run | node share across the root's 4 subtrees | max depth |
|---|---|---|
| run 1 | **50% / 42% / 7% / 2%** | 6 |
| run 2 | **51% / 43% / 4% / 2%** | 6 |
| run 3 (13 nodes so far) | 42% / 42% / 8% / 8% | 2 |

Both runs put **92–94% of their nodes into two of four root branches** and drove
to maximum depth. That is not breadth-first; it is close to its opposite.

The old rule's bonus, `0.6·sqrt(log(parent_n + 1))`, runs from 0.50 at
`parent_n = 1` to 1.41 at `parent_n = 250` — a span of 0.91 against a value
spread of about 0.31. So the bonus dominated, and what it rewarded was *having
an already heavily-visited parent*. That is a rich-get-richer loop: whichever
branch took an early lead collected nearly every subsequent fold, and node value
only ever broke ties between siblings sharing a parent.

This accounts for node 192 better than the old story did — it was not passed
over for being deep, but for hanging off a starved branch. And it accounts for
run 2's depth signal arriving with no gain in the best design: the search drilled
to the bottom of two branches and exhausted them while the reachable set stayed
narrow.

Corrected in RUN2_RESULT.md and PREPRINT_DRAFT.md. **Subtree concentration
replaces secondary 1** as the test of whether the fix bit — it discriminates
cleanly between the two rules, where visit counts do not. Run 3 at 13 nodes is
far too early to read, and the comparison is only fair at similar node counts.

## 2026-09-20 01:35 — 0.820 at fold 52; and the replacement secondary is weak too

**Progress.** 116/1000 folds, 29 nodes, advancing. Best node value 0.333
(node 14, depth 2). Depth 3 reached.

**Best single sequence 0.820 at node 13 — 52 folds in.** That equals run 2's
*final* best, which took 955 folds and 192 nodes to reach.

| | run 2 (broken UCT) | run 3 (fixed UCT) |
|---|---|---|
| first ≥0.78 | node 8 | node 6 |
| **reached 0.820** | **node 192** | **node 13** |

Worth writing down, and worth not over-reading. It is a single draw, the K=4
disagreement on RBD measured 0.167 in NOISE_RESULT.md, and run 2 got most of the
way there by node 8 as well. The honest statement is that run 3 is ahead of run
2's pace at the same budget and that this is one sample. It becomes interesting
only if the run goes on to clear 0.870, which nothing here yet suggests.

Children beating their own parent: 15 of 28 = 53.6%, CI [35.1, 72.0]. Still
uninformative.

### The replacement secondary is also underpowered. Say so now.

Comparing subtree concentration at a **matched node count of 27**, rather than
run 3's 27 against the finished runs' ~250:

| run | shares | top-2 combined |
|---|---|---|
| run 1 (broken) | 74/19/4/4 | **93%** |
| run 2 (broken) | 33/30/19/19 | **63%** |
| run 3 (fixed) | 44/19/19/19 | **63%** |
| *run 1, final (252)* | *50/42/7/2* | *91%* |
| *run 2, final (236)* | *51/43/4/2* | *94%* |

Two problems, both visible only because the comparison was made at matched
counts:

1. **Run 3 is indistinguishable from run 2 here** — both 63%. The metric has no
   discriminating power at this stage.
2. **Run 1 and run 2 differ enormously from each other at 27 nodes (93% vs
   63%) despite running the identical broken rule.** Run 1 concentrated early;
   run 2 concentrated late, ending at 94%. So between-seed variance in this
   metric is at least as large as the effect it is meant to detect.

With one run per arm, subtree concentration cannot separate the rules. It may
still be readable at ~200 nodes, where both broken runs converged to 91–94% and
a fixed rule ought not to — but that is a single comparison against a quantity
whose spread we cannot estimate, and it will not carry much weight.

**Both pre-registered mechanism tests have now failed on power, for different
reasons: secondary 1 was structurally uninformative, and its replacement is
swamped by between-run variance.** Run 3 can still answer the primary question —
does a correctly implemented search beat 0.870 — but the claim "the fix changed
how the search behaves" is not going to be demonstrable from n=1 per arm. That
limitation belongs in the writeup rather than being discovered by a reader.

## 2026-09-20 02:39 — flat since fold 52; best-node test favours the fix; node 192 claim verified

**Progress.** 196/1000 folds, 49 nodes, advancing, driver alive.

**Best single sequence still 0.820, unchanged since node 13.** 144 folds and 35
nodes have bought nothing. Best node value also flat at 0.333 (node 14). The
same shape as runs 1 and 2 — nearly all of the result arrives early, and the
rest of the budget is spent confirming it.

Children beating their own parent has drifted down to 19 of 48 = **39.6%**,
CI [25.7, 53.4] — now sitting on run 2's final 40%, and still not separable from
chance.

Depth is spreading rather than descending: 16 at depth 2, 24 at depth 3, only 4
at depth 4. Run 2 at this stage was also shallow, so nothing to read yet.

### A mechanism test that does discriminate

Rather than subtree shares, ask directly what failed in run 2: **does the search
expand its own best-known nodes?**

| | #1 node expanded? | top-10 expanded |
|---|---|---|
| run 2 (broken), at 48 expansions | yes (node 184) | 8 of 10 |
| run 2 (broken), final | **no** (node 192, 0.427) | **5 of 10** |
| run 3 (fixed), at 13 expansions | **yes** (node 14, 0.333) | 4 of 10 |

Run 3 has expanded its #1, #3 and #4 nodes. It is the right question to ask, but
13 expansions is too few to compare against run 2's 60, and run 2 looked fine at
48 as well. Recheck at ~60 expansions, which is where run 2's behaviour fell
apart.

### The node-192 claim holds, and 112/113 are better evidence

RUN2_RESULT.md and the preprint both lean on "node 192, the highest-valued node
in the tree, was never selected across roughly ten subsequent opportunities". A
node born near the end of a run has few chances by construction, so this was
worth checking rather than repeating. Counting expansions after each node's
birth in run 2's log:

| node | value | expansions after it was created | ever chosen |
|---|---|---|---|
| 192 | 0.427 | 12 | never |
| 218 | 0.398 | 5 | never |
| **113** | **0.375** | **31** | **never** |
| **112** | **0.370** | **32** | **never** |

The claim survives as stated — 12 chances, never taken — and nodes 112 and 113
are the stronger cases: ranked 4th and 5th in the whole tree, each passed over
**more than thirty times**. Worth substituting into the writeup, since neither
can be waved away as an end-of-run artifact.

## 2026-09-20 03:33 — 260/1000, still flat at 0.820; concentration trajectories diverging

**Progress.** 260/1000 folds, 65 nodes, advancing. Best node value 0.333 (node
14), best single sequence **0.820** — both unchanged since node 13. That is now
**208 folds with no gain**, four times the budget that produced the result.

Children beating their own parent back up to 29 of 64 = 45.3%, CI [33.1, 57.5].
It has wandered 67% → 54% → 40% → 45% as n grew, which is what a statistic with
no signal looks like. Depth 4 now has 8 nodes with mean 0.218, the highest of
any depth, but on 8 nodes that is nothing.

**Best-node test, matched at 16 expansions: no discrimination.** Both runs
expanded 5 of their top 10. Run 3's #2, #3 and #4 were created in the two most
recent expansions and have had no opportunity yet. Run 2's divergence only
appeared near its end, so this needs ~60 expansions to mean anything; run 3 is
at 16.

**Subtree concentration is the one thing separating them, and it is the
trajectory rather than the value that looks informative:**

| nodes | run 1 (broken) | run 2 (broken) | run 3 (fixed) |
|---|---|---|---|
| 27 | 93% | 63% | 63% |
| 64 | 97% | **78%** | **65%** |
| final (~250) | 91% | 94% | — |

Run 2 climbed 63 → 78 → 94 as its rich-get-richer loop took hold. Run 3 has gone
63 → 65. A flat trajectory where run 2 climbed would be better evidence than any
single matched-count comparison, because it does not depend on the two runs
having started from the same place.

Still not conclusive: run 1 and run 2 differ by 19 points at 64 nodes under the
*identical* rule, so run 3's 13-point gap below run 2 sits inside the spread
between two runs of the same arm. Worth tracking to ~200 nodes, where both
broken runs had converged above 90%. If run 3 is still near 65% there, that is
the strongest statement this design can support — and it is still n=1 per arm.

## 2026-09-20 04:27 — best node value in any run (0.463), and it has not been chosen

**Progress.** 324/1000 folds, 81 nodes, advancing.

**Node 75 scored 0.463 — the highest node value in any of the three runs**, past
run 2's final best of 0.427. Best single sequence is still **0.820**, unchanged
since node 13 (272 folds ago), now tied by node 58.

Children beating their own parent: 42 of 80 = **52.5%**, CI [41.6, 63.4] — above
50% for the first time, CI still spanning chance.

### Watch item: the run-2 pathology may be recurring

| rank | node | value | best seq | expanded? |
|---|---|---|---|---|
| **#1** | **75** | **0.463** | 0.790 | **not yet — 2 chances** |
| #2 | 14 | 0.333 | 0.630 | expanded |
| #3 | 42 | 0.325 | 0.800 | never |
| #4 | 63 | 0.323 | 0.530 | never |

Since node 75 was created the search has expanded nodes **19 and 9** — node 9's
value is 0.150, less than a third of node 75's. This is precisely what the fix
was meant to prevent, and it is the first sign of it under the corrected rule.

Two chances is far too few to call, and node 75 is three minutes old. But a
plausible mechanism is already visible in the depth histogram: **48 nodes at
depth 3, only 12 at depth 4.** Within a sibling set the exploration bonus is
constant and value decides, which is the part the fix repaired — but *along the
descent* an unexpanded mediocre child (n = 1, full bonus) still outranks an
expanded good one (n > 1, reduced bonus). So the descent keeps stopping at the
first level it finds unexpanded nodes, filling level 3 broadly instead of
pushing into level 4. Even after normalising values to [0, 1], `c_uct` = 0.6
gives a bonus of 0.80 at `parent_n` = 5 and 1.05 at 20 — still comparable to or
larger than the entire value range.

If node 75 is still unexpanded in ~10 more expansions, **the constant is still
too large after normalisation** and that is a third finding about this
implementation, not a repeat of the second.

### A clean live illustration of mean vs max

The two statistics have fully decoupled in this run:

| | best node value | best single sequence |
|---|---|---|
| node 75 | **0.463** (rank 1) | 0.790 (rank 4) |
| node 13 | 0.302 (rank 12) | **0.820** (rank 1) |

The node best worth *searching from* and the node carrying the best design to
*synthesise* are different nodes, and neither ranking predicts the other. This is
NOISE_RESULT.md's argument appearing directly in a live run rather than in a
resampling experiment — worth a sentence in §3.3.

## 2026-09-20 05:34 — node 75 WAS expanded, and produced nothing. The cleanest result yet.

**Progress.** 404/1000 folds, 101 nodes, advancing. Best single sequence still
**0.820** — unchanged for 352 folds. Best node value still 0.463 (node 75).
Children beating their own parent back to 46.0%, CI [36.2, 55.8].

### The watch item from 04:27 resolves in the fix's favour

Node 75 was chosen on the **third** expansion after it appeared — order since its
birth: 19, 9, **75**, 27, 33, 54, 55. The corrected rule does select its
best-known node. That is the behaviour run 2 failed at, where nodes ranked 4th
and 5th were passed over more than thirty times each.

The "constant may still be too large after normalisation" hypothesis from the
last entry is **not supported**. Withdrawn.

### And expanding it produced nothing

Node 75's four children: **0.140, 0.230, 0.180, 0.138**, against the parent's
0.463. The best node in the tree yielded four children all well below the tree
average, and none came close to the parent.

**This is the finding.** Run 1 logged the same anomaly — "the best-known node
yielded nothing on three separate expansions" — and RUN2_RESULT.md attributed
it, along with everything else, to the broken selection rule. That attribution
is now refuted from inside a run where selection worked correctly:

> Selection chose the best node. Expansion still produced nothing.

So "expanding the best node is not rewarded" is **not an artifact of the broken
UCT**. It is a property of this search space: a node's value does not survive
one round of partial re-diffusion at `partial_T = 10`. Whatever made node 75
score 0.463 is not inherited by its structural neighbours.

Corroborating, from the same window: the second-best node in the tree, **node 98
at 0.440, is a child of node 54, whose own value is 0.285** — below the top six.
Good children keep arriving from ordinary parents. Run 1 recorded exactly this
("the largest gains arrived from low-scoring parents") and it survives the fix.

**What this means for the paper.** It separates two claims that RUN2_RESULT.md
had merged. The selection bug was real and is fixed. But the reason guided
search does not beat blind sampling here is upstream of selection: **the value
signal does not propagate across the edges the search moves along.** A search
cannot exploit structure that is not there, however correctly it selects. This
belongs in §3.3 as the primary explanation, with the selection bug demoted to an
implementation note — the reverse of the current draft.

Provisional until the run finishes: one parent, four children. Node 75 is a
single case, and the run-1 precedent is three more. Worth counting properly at
the end — for every expanded node, the gap between its own value and its best
child's — which is a whole-run statistic rather than an anecdote.

## 2026-09-20 06:28 — the whole-run number, and a reliability estimate that explains it

**Progress.** 468/1000 folds, 117 nodes, advancing. Best single sequence
**0.820** (flat for 416 folds), best node value 0.463. Children beating their
own parent 44.8%, CI [35.8, 53.9]. Depth 4 now has 40 nodes.

### Parent value does not predict child value — in any run

Promised last entry: the whole-run statistic rather than the node-75 anecdote.
For every expanded node, its own value against its best child's:

| run | pairs | Spearman ρ | p | mean gap | best child > parent |
|---|---|---|---|---|---|
| run 1 (max, broken) | 63 | 0.058 | 0.649 | +0.037 | 54% |
| run 2 (mean, broken) | 59 | 0.162 | 0.219 | +0.044 | 71% |
| **run 3 (mean, fixed)** | 29 | **0.168** | **0.384** | +0.050 | 72% |

**No correlation anywhere, and fixing selection did not change it.** The search
moves along edges that carry no information about value.

And the top of the distribution is actively misleading — top-decile parents
produce children *worse than themselves* in all three runs, while the average
parent improves:

| run | n | parent mean | best child |
|---|---|---|---|
| run 1 | 9 | 0.816 | **0.530** |
| run 2 | 6 | 0.365 | **0.311** |
| run 3 | 4 | 0.352 | **0.253** |

Ordinary parents gain (+0.04 to +0.05 on average); the best parents lose. That
is regression to the mean, and it means **expanding the best node is the worst
use of a fold** — which inverts the premise of best-first search.

### Why: less than half the node-value variance is real

Run 3's 116 node values have mean 0.183 and **sd 0.067**. The measured K=4
mean-value disagreement on RBD is 0.059 (NOISE_RESULT.md); for a difference of
two independent means, `E|A−B| = σ_e·1.128`, so the per-node measurement error
is **σ_e = 0.052**.

    observed sd        0.067
    measurement error  0.052
    true-signal sd     0.042
    RELIABILITY        0.39

**Only 39% of the variance in node values is real.** The search is selecting on
differences that are 61% noise, so a node selected for a high value is mostly
selected for a lucky draw, and its children re-measure without that luck.

The arithmetic predicts node 75 quantitatively. Regressing its value toward the
mean by the reliability, `0.183 + 0.39×(0.463 − 0.183) = 0.291`; observed best
child **0.230**, on four draws. That is the right size.

### This is the paper's §3.3, and it is a different argument than the draft makes

The draft says max-of-K is a tail statistic and mean fixes it. True, and the
mean did help: K=4 noise falls from 0.167 (max) to 0.059 (mean). **But it did
not help enough.** At reliability 0.39 the node values are still mostly noise,
and every downstream symptom follows from that one number — no parent-child
correlation, best parents regressing, depth carrying nothing, and search losing
to blind sampling at 16× the cost.

The three runs now separate the causes cleanly:

| | fixed? | did it help? |
|---|---|---|
| wrong statistic (max → mean) | yes, run 2 | noise 0.167 → 0.059; not enough |
| broken selection rule | yes, run 3 | best node now gets expanded; no gain |
| **node values 61% noise** | **no** | **the binding constraint** |

The remedy for the third is more folds per node, and it is priced: halving σ_e
needs K = 16, quadrupling the cost of every node. At 1000 folds that buys 62
nodes instead of 250 — and the reliability would still only reach about 0.7.
Worth stating explicitly in the discussion, because it converts "search did not
work" into "here is what a search would cost to work, and it is not affordable
on this hardware".

## 2026-09-20 07:36 — the fix demonstrably changed behaviour; secondary 2 falsified

**Progress.** 548/1000 folds, 137 nodes, advancing. Best single sequence
**0.820** (flat 496 folds), best node value 0.463. Children beating their own
parent drifting down again: 42.6%, CI [34.3, 51.0].

### Concentration: the first secondary that discriminates

Top-2 subtree share against node count, all three runs:

| nodes | 27 | 64 | 100 | 136 | 180 | 236 |
|---|---|---|---|---|---|---|
| run 1 (broken) | 92% | 97% | 91% | 84% | 88% | 91% |
| run 2 (broken) | 62% | 78% | 86% | 90% | 92% | 94% |
| **run 3 (fixed)** | 62% | 65% | **70%** | **69%** | — | — |

The variance problem I flagged at 01:35 has resolved itself with more data.
Runs 1 and 2 began 30 points apart (92 vs 62) and **converged to 84–90% by 136
nodes**; run 3 has gone flat at 65–70%. At matched count it now sits 15–21
points below both broken runs, and the separation is in the *trajectory*, which
does not depend on the runs starting alike.

**The corrected selection rule measurably changed how the search explores.** It
has not changed what the search finds. Those two facts together are a cleaner
result than either alone.

### Secondary 2 was predicted backwards

The pre-registration said "a correctly weighted UCT should push more mass to
depth 5–6". At matched 136 nodes:

| run | depth histogram | max depth |
|---|---|---|
| run 1 (broken) | 1/4/16/32/32/31/20 | **6** |
| run 2 (broken) | 1/4/16/27/36/36/16 | **6** |
| **run 3 (fixed)** | 1/4/16/**60/55** | **4** |

Run 3 is **shallower**, not deeper — it has not reached depth 5 at all, while
both broken runs were at depth 6 by this point. My prediction was wrong, and in
hindsight it was wrong for a reason I had already written down at 04:27: along
the descent an unexpanded child (n = 1, full exploration bonus) outranks an
expanded one (n > 1, reduced bonus), so a correct UCT stops at the shallowest
level still holding unexpanded nodes. Breadth before depth is what UCT is
supposed to do. The broken rule's depth-6 reach was the rich-get-richer loop
tunnelling, not healthy search.

So both pre-registered mechanism secondaries failed, and for opposite reasons:
secondary 1 could not discriminate at all, and secondary 2 predicted the wrong
direction. The unregistered replacement — subtree concentration — is the only
one that worked. **That ordering has to be reported honestly: the test that
succeeded was chosen after seeing data, and carries correspondingly less
weight than the two registered ones that failed.**

The primary is unaffected — it was fixed in advance and is still 0.870.

## 2026-09-20 08:30 — a prediction, written before the data exists

**Progress.** 612/1000 folds, 153 nodes, advancing. Best single sequence
**0.820** (flat 560 folds), best node value 0.463. Children beating their own
parent 42.8%, CI [34.9, 50.6]. Depth still capped at 4 (72 nodes there now).

Selection continues to do its job: having expanded node 75 (#1), the search is
now expanding **node 98, value 0.440 — the #2 node in the tree**. Node 149 was
also just created at 0.435, third-highest in the run.

### Testable prediction, logged now

The reliability model from 06:28 says node values are 38% signal. Population as
of this moment: n = 152, mean 0.183, sd 0.067, reliability 0.38.

**Node 98 is being expanded as this is written. Its four children do not yet
exist.** The model predicts the best of them at:

    0.183 + 0.38 x (0.440 - 0.183)  =  0.282

against **0.440** if node values were real. The gap between those two numbers is
the whole claim, and it will be settled by the next check with no opportunity to
adjust either figure.

Precedent: node 75 (value 0.463) predicted **0.291**, observed **0.230**.

A hit near 0.28 would make this the second quantitative confirmation. A best
child near 0.44 would falsify the reliability account and mean node 75 was a
fluke — in which case the search is working and the shortfall is elsewhere.
Either way it is decided by a number that already exists on the GPU and not by
anything chosen afterwards.

*(Caveat, stated now rather than after: n = 2 predictions. The proper version is
the whole-run regression of best-child on parent value across all expansions,
which is what RUN3_RESULT.md should carry. This is a spot check that happens to
be available and pre-committable.)*

## 2026-09-20 09:38 — the prediction missed; regression is real but weaker than claimed

**Progress.** 692/1000 folds, 173 nodes, advancing. Best single sequence
**0.820** (flat 640 folds), best node value 0.463. **Depth 5 reached** — 16
nodes, mean 0.255, the highest of any depth. Children beating their own parent
41.9%, CI [34.5, 49.2] — now formally **below** chance.

### The node-98 prediction was wrong

Logged at 08:30, before the children existed:

    predicted best child        0.282
    naive (node value is real)  0.440
    OBSERVED                    0.380    <- children 0.380 / 0.307 / 0.260 / 0.150

Model error 0.098; naive error 0.060. **The naive model won this one.** I said
at 08:30 that a result near 0.44 would falsify the reliability account and one
near 0.28 would confirm it. 0.380 is neither, and it is closer to naive.

Score on pre-committed predictions: **1 hit, 1 miss.** Node 75 predicted 0.291,
observed 0.230 (model error 0.061, naive 0.233 — model won decisively). Node 98
the other way.

### What the whole-run numbers actually support

The proper version, all 42 expansions with scored children:

    reliability (re-estimated)     0.47   <- was 0.38 at 06:28
    mean |predicted - observed|    0.081
    mean |parent    - observed|    0.097  (naive)
    top-quartile parents: 0.342 -> predicted 0.261, observed 0.279

Two things to correct from the 06:28 entry:

1. **The reliability estimate is not stable.** It moved 0.38 → 0.47 as 60 more
   nodes arrived, because it is derived from the sample sd of node values, which
   is still growing. Any number quoted from it carries that wobble, and 06:28
   quoted 0.39 to two figures as though it were measured.
2. **The model beats naive by 0.081 against 0.097 — about 16%.** That is a real
   edge but a modest one, not the decisive account 06:28 implied. Observed
   top-quartile children (0.279) land *above* the model's prediction (0.261) and
   well below the parent (0.342).

**What survives.** Regression to the mean is real and reproducible: top-quartile
parents average 0.342 and their best children 0.279, a drop of 0.063, and the
parent-child correlation is still null. So "expanding the best node is not
rewarded" stands on direct measurement.

**What does not.** The claim that 61% of node-value variance is noise, and that
this quantitatively explains the shortfall, is over-stated. The reliability
figure is unstable and its point predictions are only slightly better than
assuming node values are real. RUN3_RESULT.md should report the regression
effect and the null correlation, and treat the variance decomposition as a
plausible mechanism rather than a measured one.

This is the second time in this run that a mechanism I found compelling has not
survived its own test — the first was "both runs were breadth-first". Worth
noting as a pattern: the diagnostic claims in this project have a much worse
track record than the descriptive measurements.

## 2026-09-20 10:33 — 756/1000; "flat" was wrong, it is climbing slowly

**Progress.** 756/1000 folds, 189 nodes, advancing. Best single sequence
**0.820** (flat 704 folds), best node value 0.463. Children beating their own
parent 42.6%, CI [35.5, 49.6], still below chance. ~244 folds left, roughly
3 hours.

### Correction to the 07:36 entry

I wrote that run 3 "has gone flat at 65–70%" on subtree concentration. Updated:

| nodes | 27 | 64 | 100 | 136 | 180 |
|---|---|---|---|---|---|
| run 1 (broken) | 92% | 97% | 91% | 84% | 88% |
| run 2 (broken) | 62% | 78% | 86% | 90% | 92% |
| run 3 (fixed) | 62% | 65% | 70% | 69% | **74%** |

It is not flat — it is **climbing, just more slowly**: 62 → 65 → 70 → 69 → 74.
The correct statement is that run 3 concentrates at roughly half the rate of the
broken runs and sits 14–18 points below both at matched node counts, not that it
avoids concentrating. With about 60 nodes of budget left it will not be known
whether it converges to their 90%+ or plateaus lower.

This is the third claim in this run that needed walking back after more data
(after "breadth-first" and the reliability point-prediction), and all three erred
the same way: I read a trend from the first few points and stated it more firmly
than the sample supported. **For RUN3_RESULT.md the rule should be that no
trajectory claim gets made below ~150 nodes.**

Depth 5 also illustrates it — its mean was 0.255 on 16 nodes at 09:38 and is
0.227 on 24 nodes now. Early means on small n run high here, every time.

Nothing else changed.

## 2026-09-20 11:27 — new best 0.830, and every record came from a bad parent

**Progress.** 820/1000 folds, 205 nodes, advancing. **New best single sequence
0.830** (node 201, depth 4) — the first improvement since node 13 at fold 52.
Best node value still 0.463. Children beating their own parent 42.6%, below
chance.

    0.820 at fold 52   ->   0.830 at fold 820
    +0.010 for 768 folds

Run 3 now sits above run 2 (0.820), below run 1 (0.860), and below the bar
(0.870) — inside the pre-registered middle band.

### Every record-setter in this run descends from a low-value parent

| record | best seq | parent | parent's value |
|---|---|---|---|
| node 6 | 0.670 | node 3 | **0.150** |
| node 13 | 0.820 | node 1 | **0.133** |
| **node 201** | **0.830** | node 22 | **0.135** |

Population mean node value is 0.19. **All three parents sit at or below 0.15 —
the bottom of the distribution.** Meanwhile the two best-valued nodes in the
tree, 75 (0.463) and 98 (0.440), were both expanded and produced nothing that
set a record.

Node 201's own value is 0.427 — near the top of the tree — and it came from a
parent at 0.135. The value signal does not propagate across a single
partial-diffusion step in either direction.

This is run 1's "the largest gains arrived from low-scoring parents" reproduced
under a **correct** selection rule, which is what makes it worth something.
In run 1 it could be dismissed as an artifact of the broken UCT steering the
search into bad regions. Here selection demonstrably worked — it expanded its
top two nodes — and the pattern held anyway.

Combined with the null parent-child correlation (ρ = 0.168, p = 0.384) and
top-quartile parents regressing (0.342 → 0.279), three independent views now say
the same thing: **a node's score carries no information about what expanding it
will produce.** That is the cleanest statement this project can make about why
guided search does not beat blind sampling here, and unlike the reliability
decomposition it rests on direct measurement rather than a modelled variance
split.

~180 folds left, roughly 2 hours.

## 2026-09-20 12:35 — 900/1000, nothing changed

Best single sequence **0.830**, best node value 0.463, children beating their
own parent 43.3% [36.8, 49.8]. Depth 5 has 28 nodes; depth 6 never reached.
~100 folds left, finishing around 14:00.

## 2026-09-20 13:29 — 964/1000, final stretch

Best single sequence **0.830**, best node value 0.463, unchanged. Children
beating their own parent drifted up to 45.8% [39.5, 52.1] — back inside chance
rather than below it, which is what that statistic has done all run.

Two or three expansions left; the loop runs until `folds_used >= 1000`, so it
will stop near 1012 — coincidentally the exact fold count of run 1. Depth 6 will
not be reached.
