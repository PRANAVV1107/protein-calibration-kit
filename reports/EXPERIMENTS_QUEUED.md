# Queued experiments — prepared 2026-09-20, none run

Everything below is written and ready. Nothing has been launched. The GPU is idle.

**Two findings emerged while preparing these, and they reorder the queue.** Both
are in section 0 because they change how the existing results should be read.

---

## 0. What preparation turned up

### 0a. The pre-registered bar and the search were measured on different instruments

All three searches folded with a **donated target alignment at 1 model / 3
recycles** ("hybrid"). The bar they were judged against — 0.870, from 60 blind
folds across 10 RBD backbones — was measured in the variance study at **full
MSA, 3 models, 8 recycles**. `reports/MCTS_PREREGISTRATION.md` line 51 records
the bar's provenance as "variance study, high accuracy"; `scripts/mcts_design.py`
line 310 folds at `--num-models 1 --num-recycle 3`.

On the ten backbones both datasets cover, the instruments disagree badly:

| backbone | hybrid | high-acc | diff |
|---|---|---|---|
| binder_42 | 0.330 | **0.870** | +0.54 |
| binder_5 | 0.410 | **0.860** | +0.45 |
| binder_47 | 0.160 | 0.520 | +0.36 |
| binder_10 | 0.180 | 0.540 | +0.36 |
| binder_14 | 0.150 | 0.460 | +0.31 |
| binder_34 | 0.480 | 0.710 | +0.23 |
| binder_13 | 0.790 | 0.810 | +0.02 |
| binder_9 | 0.750 | 0.770 | +0.02 |
| binder_0 | 0.750 | 0.730 | -0.02 |
| binder_15 | 0.730 | 0.710 | -0.02 |

High accuracy reads **+0.18 higher on average** and **reorders the set** —
binder_42 is third-worst on hybrid and best on high-accuracy. So "search reached
0.860 against a bar of 0.870" compared a number from the low-reading instrument
against one from the high-reading one.

**This does not mean the search won.** The search *optimised* the hybrid metric,
and the two instruments disagree, so its winners need not score well at full
MSA. It means the comparison was never valid in either direction, and the
primary result of all three runs is currently unknown rather than null.
Experiment 1 settles it.

### 0b. The searches started from the second-worst root available

Re-running the root comparison on the high-accuracy instrument (no new folds —
`results/variance_study`, 10 RBD backbones x 6 draws):

| root | best-of-6 | E[best-of-6] |
|---|---|---|
| binder_5 | 0.860 | **0.827** |
| binder_42 | 0.870 | 0.814 |
| binder_9 | 0.770 | 0.764 |
| ... | | |
| **binder_47** | **0.520** | **0.485**  <- every search ran from here |
| binder_14 | 0.460 | 0.371 |

Between-root spread of E[best-of-6]: **0.371 to 0.827**, sd 0.154.

| strategy | folds | E[best] |
|---|---|---|
| one random root | 6 | 0.644 |
| best of 5 random roots | 30 | 0.789 |
| **best of 10 random roots** | **60** | **0.816** |
| search from binder_47 | **1012** | 0.830 *(hybrid instrument)* |

**Sixty folds spread across ten roots reaches roughly what 1012 folds of search
from one root reached.** Root selection is ~17x cheaper than search and was
never done on evidence — binder_47 was picked arbitrarily and is 9th of 10.

Run: `python scripts/root_sweep.py --target rbd` (instant, no GPU, already run;
output in `results/root_sweep_rbd.csv`).

### 0c. Multi-root beats single-root at every budget, and the gap widens

Simulated directly from `results/rl/noise_poc` (10 RBD backbones x 16 draws,
folded in the **same hybrid configuration the searches used**, so the numbers sit
on one axis). E[best ipTM] for a budget of B folds split R roots x B/R draws,
roots drawn at random each replicate:

| budget | R=1 | R=2 | R=5 | R=10 |
|---|---|---|---|---|
| 10 | 0.365 | 0.400 | 0.433 | **0.435** |
| 20 | 0.417 | 0.493 | 0.545 | **0.563** |
| 40 | 0.456 | 0.572 | 0.650 | **0.676** |
| 80 | 0.485 | 0.615 | 0.720 | **0.747** |
| 160 | 0.469 | 0.632 | 0.750 | **0.774** |

Breadth wins at every budget and the margin grows: **+0.070 at B=10, +0.305 at
B=160**. More importantly, **single-root saturates** (0.365 -> 0.485 -> 0.469,
flat past B~80) while multi-root is still climbing at the edge of the data.
The high-accuracy cross-check agrees on direction.

Run: `python scripts/root_allocation.py --target rbd` (no GPU).

### 0d. The search was compared against a multi-root baseline while being single-root

This is the reframe that matters. On the hybrid instrument:

| | folds | backbones | best ipTM |
|---|---|---|---|
| blind, binder_47 alone | 16 | 1 | **0.160** |
| **search from binder_47** | **1012** | **253** | **0.830** |
| blind, 10 *fresh* backbones | 160 | 10 | **0.790** |

**Against a matched single-root baseline the search wins overwhelmingly** —
0.830 against 0.160, an improvement of +0.67 over the root it started from. The
method is not inert. It escapes a bad root, which is exactly what it is for.

What it loses to is *root diversity*. The 0.870 bar was 60 blind folds across
**ten** backbones — a multi-root baseline — so all three runs were judged
single-root against multi-root. That is not the comparison the pre-registration
intended to make, and it is not a fair one.

And the reason breadth wins is measurable: **253 backbones produced by partial
re-diffusion from one root are worth about the same as 10 freshly generated
ones** (0.830 against 0.790). Partial diffusion at `partial_T = 10` yields far
less diversity per backbone than fresh generation does. Which is the same
conclusion the locality measurements reached from the other direction, and it
makes Experiment 3 (`partial_T` sweep) more interesting rather than less: if
locality exists at low T, then high T is the only regime where re-diffusion
buys diversity, and there it destroys the value signal. That tension may be
fundamental to the approach.


---

## Experiment 1 — re-score the search designs on the bar's instrument

**Why first:** 0a makes the headline result of three runs and ~2,900 folds
uninterpretable. Nothing else is worth running until this is settled.

    scripts/rescore_prepare.py     select designs  (no GPU)
    scripts/rescore_winners.sh     fold at full MSA / 3 models / 8 recycles
    scripts/rescore_analyze.py     compare against 0.870

**Design.** Per run: 15 top designs by hybrid ipTM, plus 5 drawn at random
across the whole range. The random arm is essential — with only the winners, a
real gain cannot be told apart from the systematic +0.18 offset.

**Cost:** 60 designs x ~240 s = **about 4 hours**.

**Outcomes.** If the best re-scored design clears 0.870, the search did beat
blind sampling and three reports need rewriting. If it does not, the null
survives on a valid comparison and is much harder to attack. If the top arm is
no better than the random arm, then the hybrid objective the search maximised
does not identify good designs at all — which would make the search's
*objective* the thing that failed, not its algorithm.

    bash scripts/rescore_winners.sh

## Experiment 2 — root sweep, extended

**Status:** the analysis in 0b already answers the headline question from
existing folds. A GPU extension is only worth it to raise K from 6 to 16 on the
same 10 roots (160 folds) or to add roots beyond the variance study's ten.
**Recommend deferring** until Experiment 1 resolves which instrument future work
should use.

## Experiment 3 — partial_T sweep: does the search space have locality?

**Why it matters:** run 3 measured parent-to-best-child rho = 0.168 (p = 0.384)
at `partial_T = 10`, with every record-setting design descending from a
bottom-decile parent. But `partial_T = 10` was never chosen on evidence.

    scripts/partial_t_sweep.sh      T in {2, 4, 10, 20}
    scripts/partial_t_analyze.py

**Design.** 6 parents per T, chosen by a fixed index rule rather than by score —
selecting good parents would restrict range and hide a real correlation.
4 children x 4 sequences per node.

**Cost:** 4 T-values x 6 parents x 5 nodes x 4 seqs = **480 folds, about 8
hours** at the search's own configuration.

**Outcomes.** rho climbing as T falls means the landscape has a correlation
length and every negative search result here is conditional on one badly chosen
hyperparameter. rho null at every T means tree search over RFdiffusion backbones
is structurally dead for this pipeline, with a measured mechanism — the stronger
paper. Stated limitation: 6 parents per T detects only large correlations.

    bash scripts/partial_t_sweep.sh

## Experiment 4 — second validation dataset (no GPU)

**Why:** the headline AUROC of 0.627 rests on 164 designs from one competition
against one target family (EGFR). A second independent dataset turns it from an
anecdote into a claim, and this is data wrangling rather than compute.

**Candidates found 2026-09-20, none yet verified as downloadable:**

1. **ProtDBench** (arXiv:2605.04118, May 2026) — the best fit on paper. States
   it includes designs "validated without prior AlphaFold-based filtering",
   which is exactly what is needed: AF-filtered sets are selected on the very
   metric under test, so any AUROC computed on them is biased toward the filter.
   Also benchmarks AF2-IG, ColabFold, Protenix, Boltz-1/2 and Chai-1, so our
   numbers would land beside published ones on the same designs.
2. **"Predicting Experimental Success in De Novo Binder Design: A Meta-Analysis
   of 3,766 Experimentally Characterised Binders"** (bioRxiv, Aug 2025) —
   3,766 binders across multiple targets, roughly 23x the current set.
3. **Cao et al. 2022** (*Nature* 605:551, doi 10.1038/s41586-022-04654-9) —
   large minibinder campaign with SSM data in the supplement; older and likely
   more work to parse into a binder/non-binder table.

**Next step is a person reading the data availability sections**, not a script.
Once a table exists, `scripts/compare_arms.py` and the AUC machinery in
`scripts/analyze_adaptyv.py` already do the rest — they take a
design/sequence/outcome frame and are not EGFR-specific.

**Watch for:** whether outcomes are binary or K_D, whether designs were
AF-filtered before synthesis (biases AUROC upward), and target-family overlap
with EGFR (would not be independent).

## Experiment 5 — multi-root search (NEW, from 0c/0d)

**Why:** 0c shows breadth beats depth at every budget tested; 0d shows the
search is good at escaping a bad root but was judged against root diversity it
never had. The obvious protocol is both at once: seed the tree with R fresh
roots rather than one.

**Code change needed.** `scripts/mcts_design.py --root-pdb` takes a single PDB
and `run()` builds one root node. It would need a forest: R roots evaluated up
front, then UCT selection across all of them. `uct()` already descends from
`self.root`, so the change is a virtual super-root whose children are the R real
roots — roughly 20 lines, and the existing descent logic then works unmodified.

**Design.** 10 roots drawn without reference to score, 1000 folds total, against
two baselines at the same 1000 folds: blind over the same 10 roots, and
single-root search (runs 1-3, already done). That is the three-way comparison
the project has never had.

**Cost:** ~1000 folds, about 14 hours, plus the code change.

**Run this only after Experiment 1**, which decides which instrument the
comparison should use.

---

## Priority order — corrected 2026-09-21

An earlier version of this file ranked multi-root search third, on the appeal of
a positive result. That was wrong. Run 3 already measured the confound it would
walk into: **gain runs inversely to starting quality** — top-quartile parents
averaged 0.342 and their best children 0.279, a loss, while average parents
gained about +0.05. The search's headline +0.67 came from starting at
binder_47, the second-worst root available. Seeded from good roots there may be
nothing to gain, so Experiment 5 risks spending 14 hours plus a code change
measuring regression to the mean.

| # | experiment | GPU | why here |
|---|---|---|---|
| **1** | **rescore at high accuracy** | **4 h** | three reports and the preprint state a conclusion that may be an instrument artifact; nothing else is interpretable until this lands |
| **2** | **second validation dataset** | **none** | makes the project's strongest and most novel claim cross-target instead of single-competition; costs no compute, so it runs alongside everything |
| **3** | **partial_T sweep** | **8 h** | publishable either way, and it *gates* 5 — if no T has both locality and diversity, multi-root search is multi-root blind with tree overhead |
| ~~4~~ | ~~multi-root search~~ | ~~14 h + code~~ | **CANCELLED 2026-09-21** — [ALLOCATION_RESULT.md](ALLOCATION_RESULT.md): no allocation policy beats uniform, and uniform beats the oracle |
| 5 | root sweep extension | — | largely answered in 0b/0c |

**Added to Experiment 3 at no cost:** `partial_t_analyze.py` now also tests
whether child gain correlates with parent quality, across the 6 parents per T
rather than within one tree. A strongly negative correlation answers Experiment
4 without running it.

---

## Status

    2026-09-21  Experiment 1  COMPLETE -> RESCORE_RESULT.md
    2026-09-21  Experiment 3  COMPLETE -> PARTIAL_T_RESULT.md
    2026-09-21  Experiment 4  CANCELLED -> ALLOCATION_RESULT.md
    2026-09-21  Experiment 5  answered from existing folds -> ROOT_SYSTEM.md

    Remaining: Experiment 2 (second validation dataset, no GPU).
    All GPU work is finished. The card is idle.
