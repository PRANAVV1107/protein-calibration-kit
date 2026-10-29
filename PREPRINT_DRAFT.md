# Preprint draft — working document

**Format target:** short methods paper, ~4 pages + figures (MLSB / Bioinformatics
Applications Note length). Everything below is real measured data from this
repository.

**Convention notes are in blockquotes like this** — they explain why a section is
written the way it is. Delete them before submitting.

> **Restructured 2026-09-21.** The 2026-09-18 draft was three loosely-related
> claims about cheap shortcuts. Two experiments run on 2026-09-20/21 — re-scoring
> the search's designs on the baseline's own instrument, and simulating
> allocation policies over roots — turned those into **one thesis with four
> independent measurements**. The escalating structure below (screen → search →
> allocate → mechanism) is the change; most individual numbers survive.
>
> An earlier version of §3.3 reported that the search *lost* to blind sampling.
> That comparison was invalid: the search and the baseline were scored with
> different AlphaFold configurations. Re-measured properly it is a tie. The
> correction is in the text, not hidden.
>
> **Updated 2026-10-29.** §3.4 led on an oracle-vs-uniform result (0.770
> versus 0.869) measured on one draw of ten RBD roots. Pre-registered a
> bootstrap of 10,000 such draws from the full 50-root pool before running it;
> the oracle matched uniform in 80.4% of draws and the reported gap sat at the
> 92.75th percentile of severity. Withdrawn as a general claim per the
> decision rule fixed in advance — see
> reports/ORACLE_BOOTSTRAP_PREREGISTRATION.md and
> reports/ORACLE_BOOTSTRAP_RESULT.md. §3.4 and the abstract now lead on
> roots-beat-draws, which this test does not touch, and report the oracle
> number as a single-draw illustration with the bootstrap attached.

---

## Title

Working options, strongest first:

1. **Cleverness does not pay: test-time search and adaptive allocation both lose
   to uniform sampling in AlphaFold-filtered binder design**
2. The objective is a maximum: why selection, search and allocation all fail in
   de novo binder design
3. Use your cheap screen as a filter, never as an estimator

> Option 1 states the result and is the one an ML audience forwards. Option 3 is
> the practitioner's version and makes a good subtitle. The previous draft's
> "where cheap shortcuts fail" is now too weak — we can say something sharper
> than "these things fail".

---

## Abstract

> ~150–200 words. Numbers are what make people read on.

De novo binder design pipelines generate backbones with RFdiffusion, assign
sequences with ProteinMPNN, and score complexes with AlphaFold2. A natural
instinct is to spend that scoring budget cleverly — screen cheaply, search over
promising regions, allocate toward promising starting points. We measure three
such strategies against uniform sampling on one consumer GPU and find that none
of them pays.

A single-sequence screen selects **worse than random** (−0.162 ipTM), because
the *target* fails to fold without an alignment (pLDDT 27.2 versus 84.3). A
pre-registered test-time search over denoising trajectories, run three times
across ~2,900 folds and two corrected implementation bugs, **ties** 60 folds of
blind sampling (0.870 versus 0.870). Simulated allocation over roots is no
better: **no bandit policy beats uniform at any tested budget**. An oracle
given the best-mean root is capped at that root's own ceiling and lost to
uniform on one draw of ten roots (0.770 versus 0.869); a pre-registered
bootstrap of 10,000 such draws from the full root pool found that loss in only
19.6% of draws, with a median gap of zero, so we report it as an illustration of
the mechanism rather than a general result.

One mechanism explains all three. The objective is a **maximum** over designs, a
tail statistic, while every clever strategy optimises a **mean**. Selecting on a
cheap score is worth +0.194 as a filter and −0.26 as an estimate of the same
quantity.

---

## 1. Introduction

> Short. Three paragraphs.

**¶1 — the pipeline.** RFdiffusion generates backbones; ProteinMPNN assigns
sequences; AlphaFold2-multimer scores the complex and top-scoring designs go
forward. Cite Watson 2023, Dauparas 2022, Evans 2021, Bennett 2023.

**¶2 — the gap.** AlphaFold scoring dominates the compute budget, so every
practical pipeline economises: screening without an MSA, fewer models and
recycles, test-time search in place of exhaustive sampling, adaptive allocation
toward promising regions. Each is adopted on the premise that a cheap score
identifies where to spend next. Whether it does is largely unmeasured.

**¶3 — contribution.** We measure four such strategies against a uniform
baseline at matched budget, with pre-registered success criteria, against 164
designs with experimental binding outcomes.

---

## 2. Methods

- **Targets.** SARS-CoV-2 Spike RBD (6M0J chain E, 333–526, 194 aa); PD-L1
  (5O45 chain A, 18–132, 115 aa); IL-7Rα (3DI3 chain B, 17–209, 193 aa).
- **Ground truth.** 365 designs from the Adaptyv EGFR competition, synthesised
  and assayed by their submitters; 55 binders, 310 non-binders. A balanced
  subset of 164 (54/110) was folded under every configuration compared. One
  positive control (Cetuximab, an approved anti-EGFR therapeutic) was excluded.
- **Scoring configurations.** ColabFold 1.6.2, AlphaFold2-multimer v3.
  *single-sequence*; *hybrid* (precomputed target alignment, single-sequence
  binder, 1 model / 3 recycles); *full MSA* at 1 model / 3 recycles; *full MSA*
  at 3 models / 8 recycles. **These are different instruments and we treat them
  as such** — §3.3 measures the offset between them.
- **Statistics.** ROC AUC via the rank-sum identity. Paired bootstrap resampling
  *designs*, not scores. Because binders in this set are longer than
  non-binders (125 aa against 101), every AUC is also reported stratified within
  length tertiles and after regressing length out.
- **Hardware.** One RTX 3060, 6 GB, under WSL2. ~6,000 AlphaFold predictions.

---

## 3. Results

> The arc is deliberate: each subsection tries to be cleverer than the last, and
> each fails by a larger margin. §3.5 explains why all of them had to.

### 3.1 A cheap screen is a good filter and a bad estimator

**[FIGURE 1]** target pLDDT against selection quality (four conditions), beside
the cheap-versus-accurate regression line.

Scoring the same 60 RBD complexes under four conditions differing only in how
the target's structure is supplied:

| condition | target pLDDT | selection vs random |
|---|---|---|
| no MSA (cheap screen) | 27.2 | **−0.162** |
| target structure as template | 45.8 | −0.013 |
| precomputed target MSA | 79.2 | +0.141 |
| full MSA (reference) | 84.3 | — |

pLDDT below 50 indicates a disordered prediction. **The cheap screen was scoring
interfaces against an unfolded target**, and selecting on it is worse than
choosing at random.

Against measured binding outcomes on the 164-design set, cheap-screen interface
metrics fall *below* chance while the same metrics at full accuracy do not:

| metric | cheap screen | full accuracy |
|---|---|---|
| ipSAE_min | **0.466** | 0.618 |
| pDockQ2_min | **0.461** | 0.623 |

Paired on identical designs, full accuracy beats the cheap screen by **+0.155
on ipSAE (p = 0.001)**.

**But a screen that is useless as an estimator can still be useful as a filter,
and conflating the two is the field's error as much as ours.** We separated them
by re-folding 60 designs at full accuracy: 45 chosen for high cheap-screen score
and 15 drawn at random across the range.

| | value |
|---|---|
| as a **filter** — selected versus unselected designs | **+0.194 ipTM**, p = 0.0027 |
| as an **estimator** — over-optimism above cheap score 0.75 | **−0.260** |

The error regresses on the score itself:

    Δ = −0.694 × (cheap score) + 0.282,   crossover at 0.406,   p = 0.0001

Designs the screen rates below 0.41 **improve** when measured properly; those
above it decline. Two designs the screen ranked near the bottom (0.120, 0.140)
scored **0.750** at full accuracy.

**Replicated on an independent dataset, and the single number is the wrong
summary.** Our AUC rests on one competition against one target family, so we
re-tested it against 3,760 experimentally characterised designs across 15
targets from an independent meta-analysis (Overath, Rygaard et al. 2025), which
publishes per-design ColabFold, AF2, AF3 and Boltz-1 confidence metrics.

Our own competition data appears in theirs, giving a direct replication rather
than a comparison: on the Adaptyv subset their pipeline gives ipTM AUC **0.659**
(round 2, n = 175) where we measure **0.627**. Base-rate matched to our 1:2
design, their 21 studies have median AUC **0.676** with a range of **0.489 to
1.000**, and **two fall below chance** — including *Cao et al EGFR* at 0.489, the
same target family as the competition.

**ipTM's usefulness is target-dependent to a degree no single AUC conveys.** Our
0.627 is unremarkable inside that distribution, which is the point: reporting one
number from one target family — ours or anyone's — overstates how well this is
known. Their larger set does resolve a metric ordering we could not (ipSAE above
ipTM and ipAE); on 164 designs from one family those three were
indistinguishable, which is consistent rather than contradictory.

**Use the screen to choose what to measure. Never report its numbers.**

> The four-condition dose–response turns the first half into a mechanism. The
> filter/estimator split is the half that is genuinely new, and it is the
> practitioner's takeaway.

### 3.2 Models and recycles are nearly free to cut

**[FIGURE 2]** AUC against compute, three configurations.

Same 164 designs, same alignment, varying only model count and recycles:

> **Cost column excludes MSA acquisition.** Measured end-to-end against the
> public mmseqs2 server, a fresh alignment costs roughly 130 s per design:
> 1 model / 3 recycles runs ~190 s and 3 models / 8 recycles ~600 s. The 4×
> ratio holds; the absolute figures below do not, and must be relabelled or
> re-measured before submission.

| configuration | cost/design | ipTM AUC | length-controlled |
|---|---|---|---|
| hybrid a3m, 1 model, 3 recycle | ~120 s | 0.545 | 0.521 |
| **full MSA, 1 model, 3 recycle** | **~63 s** | **0.627** | **0.620** |
| full MSA, 3 models, 8 recycles | ~240 s | 0.614 | 0.605 |

The cheap-model configuration matches full accuracy at a quarter of the compute.
Combined with §3.1, **the alignment is doing the work; sampling depth is not.**

**But "matches" means one specific thing, and we checked.** Re-folding 60
identical de novo RBD designs under both configurations, they are **not
interchangeable**:

    1 model / 3 recycles     mean 0.158   sd 0.104   max 0.580
    3 models / 8 recycles    mean 0.396   sd 0.252   max 0.870
    Spearman -0.053 (p = 0.69)    fit: cheap = 0.066 x expensive + 0.132

Thirteen designs score above 0.7 under the expensive configuration; **one of
them scores even above 0.3 under the cheap one.** The best backbone in the set
(0.870) is scored 0.160. Roughly half the offset is arithmetic — ColabFold's
`rank_001` under `--num-models 3` is a max over three models, worth a measured
+0.132 — but a constant offset would preserve ranking, and what remains is rank
disagreement, which cannot be subtracted away.

Two weak signals can each beat chance against binding outcomes without agreeing
with each other, and at AUC ≈ 0.62 that is what happens. **The licensed claim is
narrow: for separating characterised binders from non-binders, the cheap
configuration is as good and four times cheaper. It is not licensed for ranking
de novo designs against one another.**

> The paper's one unambiguously positive, immediately actionable result. One
> sentence: compute the target's alignment once per campaign and fold at 1 model
> / 3 recycles.

### 3.3 Test-time search ties blind sampling at 17× the compute

**[FIGURE 3]** best-so-far against folds, three search runs versus the blind
baseline.

We pre-registered a guided best-first search over RFdiffusion denoising
trajectories (`diffuser.partial_T`), with nodes evaluated by folding ProteinMPNN
sequences. Success criteria and the comparison baseline were fixed in advance and
committed before any result existed. Three runs, ~2,900 folds:

| | run 1 | run 2 | run 3 |
|---|---|---|---|
| node value statistic | `max` | `mean` | `mean` |
| UCT selection | broken | broken | **corrected** |
| folds | 1012 | 955 | 1012 |
| children beating their own parent | 26% | 40% | 46% |
| **best design, baseline instrument** | **0.870** | 0.820 | 0.820 |

Blind sampling across ten backbones at **60 folds** reaches **0.870**. The search
matches it at seventeen times the budget and does not exceed it.

**Establishing even that required rebuilding the comparison.** The searches
scored designs with a donated target alignment at 1 model / 3 recycles; the
baseline was measured at full MSA, 3 models, 8 recycles. On ten backbones both
had scored, the instruments differ by a clean **+0.139 offset on unselected
designs** and disagree in rank. Re-folding the search's designs on the
baseline's instrument is what produces the table above.

Three explanations were tested and eliminated:

- **Wrong statistic.** Scoring nodes by `max`-of-K rather than `mean` (§3.5).
  Corrected in run 2; node values became estimable and the best design did not
  improve.
- **Broken selection.** The UCT exploration term reduced to a function of the
  parent's visit count alone, producing a rich-get-richer loop that funnelled
  92–94% of nodes into two of the root's four branches. Corrected in run 3;
  subtree concentration fell to 74%, the search began expanding its own
  best-valued nodes, and the best design did not improve.
- **Wrong step size.** `partial_T = 10` was never chosen on evidence. Sweeping
  T ∈ {2, 4, 10, 20}, parent→child correlation shows **no trend with T**
  (ρ = 0.257, 0.600, 1.000, 0.771 at n = 6 each) — and is *lowest* at the
  shortest step, the opposite of what a correlation length predicts.

**What the sweep did find is a distinction between two regimes:**

| | ρ | n |
|---|---|---|
| across **independently generated** backbones | **+0.771** (Pearson +0.989, p = 0.0002) | 6 |
| within **one lineage** of re-diffused variants | +0.185 (p = 0.146) | 62 |

Parent spreads are comparable (0.098, 0.087), so this is not restriction of
range. Within a lineage every node shares the original backbone's quality and
only noise varies; across backbones the variation is real. **A tree search
descending within one lineage operates precisely where there is no signal.**

> Report the pre-registration, the two corrected bugs and the tie plainly. A
> negative result that survived two implementation fixes and an instrument
> correction is worth far more than one that did not, and reviewers can tell.

### 3.4 Roots beat draws, and allocation policies do not improve on uniform

**[FIGURE 4]** E[best] against budget, roots (R) versus draws (K), with the
single-root oracle's bootstrapped severity distribution inset.

**The robust result: spend a fixed budget on more roots, not more draws per
root.** Extending the root pool to 50 roots × 4 draws (200 folds) gives E[best]
for a budget spent as R roots × K draws:

| | K=1 | K=2 | K=3 | K=4 |
|---|---|---|---|---|
| R=5 | 0.337 | 0.398 | 0.429 | 0.453 |
| R=10 | 0.440 | 0.517 | 0.562 | 0.594 |
| R=25 | 0.629 | 0.720 | 0.764 | 0.783 |
| R=50 | **0.765** | **0.836** | **0.852** | **0.856** |

At matched budget, roots beat draws every time — +0.075 at B=20, +0.044 at
B=50, +0.049 at B=100 — and **roots are worth roughly 4× draws fold-for-fold**
(5→50 roots gains +0.428 at K=1; 1→4 draws gains +0.116 at R=5).

**The ceiling is the best root in the pool, not the number of roots.** The R
axis is still climbing steeply at K=1, and saturates at K=4 only because the
best of these 50 backbones tops out at 0.860. Three roots in fifty carry
essentially all reachable quality (median 0.210, max 0.860); the root all three
search runs descended from ranks in the bottom five at 0.13.

**Allocation policies do not improve on uniform.** Section 3.3 found signal
*between* independently generated backbones, which motivates spending a budget
adaptively — evaluate several roots, then commit more folds to the ones
scoring well. Simulated on 3,000 replicates over ten RBD roots at full
accuracy, no policy beats uniform at any tested budget:

| budget | single root | **uniform** | halving | greedy | UCB |
|---|---|---|---|---|---|
| 20 | 0.696 | **0.833** | 0.830 | 0.832 | 0.834 |
| 80 | 0.694 | **0.866** | 0.858 | 0.843 | 0.858 |
| 160 | 0.697 | **0.869** | 0.866 | 0.844 | 0.863 |

Greedy is materially worse at every budget. This held on all three targets in
the original simulation, and we have not re-tested it against resampled root
subsets (see limitations).

**An oracle given perfect knowledge of the best-mean root is capped at that
root's own ceiling, and on one draw of ten RBD roots that ceiling was 0.770
against uniform's 0.869 — a result we pre-registered to re-test before
reporting it as general, and withdrew.** Bootstrapping 10,000 draws of ten
roots from the full 50-root pool: the oracle matches uniform exactly in 80.4%
of draws (median gap 0.000, both quartiles 0.000), and the reported −0.099 gap
sits at the 92.75th percentile of severity — a tail outcome from one
favourable draw, not a property of the comparison. The mechanism is still
real — a maximum cannot be found by perfectly locating the best mean when the
best-mean root does not contain the best design — but how often that
mechanism bites on a given ten-root draw is itself root-selection-dependent,
and here it bit less than a fifth of the time.

**A cheap predictor of root quality, tested prospectively and found too weak to
use.** Of six structural features computed from the backbone PDB alone, radius
of gyration correlated with per-root quality at ρ = −0.482 (p = 0.0004, n = 50,
surviving Bonferroni) — compact binders score better. Because it was selected as
the largest of six on the data that measured it, we pre-registered a prospective
test on 150 held-out backbones: the 10 most compact against the 10 most
extended, 4 draws each.

**It replicated and fell short of the threshold we had fixed in advance.**
Compact 0.310 against extended 0.227, gap +0.083, p = 0.019; ρ = −0.374 across
all 20. The registered rule required a gap ≥ 0.10 to call it a usable screen, so
we report it as a real but unusable effect rather than a recommendation.

This is the second feature we selected post hoc and then tested. The first,
interface contact count, led at n = 10 (+0.511) and was nothing at n = 50
(+0.191). **Two for two, the selected effect shrank on honest re-testing**, by
the amount selection bias predicts. The oracle bootstrap above is a third
instance of the same discipline applied to a headline result rather than a
feature: pre-register before re-testing, and withdraw what does not survive it.

> Lead on roots-beat-draws, which is measured directly from the 50-root pool
> and not contradicted by anything here. The oracle number was the more
> dramatic claim and the one we led the abstract with; it is also the one that
> turned out to be a single favourable draw. Report the bootstrap plainly
> rather than quietly dropping the number — a reviewer who finds the original
> draw in the repo and reruns it will get 0.770 vs 0.869 exactly as stated, and
> should also find the test that shows it is not typical.

### 3.5 Why: the objective is a maximum and every strategy optimises a mean

**[FIGURE 5]** node-score disagreement against K, two targets, max versus mean.

Scoring a backbone by the maximum of K sequence draws, the disagreement between
two disjoint sets of K draws from the *same* backbone:

| target | K=2 | K=4 | K=8 |
|---|---|---|---|
| IL-7Rα, `max` | 0.141 | 0.068 | **0.034** |
| SARS-CoV-2 RBD, `max` | 0.120 | 0.167 | **0.214** |
| SARS-CoV-2 RBD, `mean` | 0.071 | 0.059 | **0.044** |

RBD's score distribution has skew **+3.38** (median 0.14; 91% of draws below
0.3, 3% above 0.7); IL-7Rα's is −0.66. The maximum is a **tail statistic**, so
on a right-skewed target its value depends on whether a subset caught an
excursion — and adding draws makes two subsets *more* likely to disagree. The
mean decays normally on both.

**This single fact accounts for §3.1, §3.3 and §3.4 together.** What a designer
wants is the best single design: a maximum. What every clever strategy estimates
and steers by — a screen score, a node value, a root's expected quality — is a
location statistic. Optimising location concentrates sampling on the
highest-*mean* region, but tails are found by sampling **many different
distributions**, not the best one harder. Hence:

- selecting on a cheap score inflates it by 0.26 and loses what it gains (§3.1);
- best-first search concentrates on high-mean branches and ties blind sampling
  (§3.3);
- an oracle that perfectly identifies the highest-mean root loses to uniform
  (§3.4).

Consistent with this, returns to compute are logarithmic and saturate: pooling
2,104 RBD folds, expected best-of-N rises from 0.727 at N=60 to 0.849 at N=1000
and gains ~+0.01 per doubling after; 0.33% of draws exceed 0.80 and none exceed
0.90.

### 3.6 Supporting: sequence choice dominates backbone choice

Within-backbone (sequence) share of total interface-score variance, 180
high-accuracy predictions over 30 backbones, 6 draws each:

| target | ipTM | ipSAE_min | pDockQ2_min |
|---|---|---|---|
| IL-7Rα | 59.0% | 43.4% | 53.4% |
| PD-L1 | 68.9% | 55.1% | 66.4% |
| SARS-CoV-2 RBD | 78.3% | 79.8% | 78.9% |

Best-of-six rather than one raises mean ipTM by 0.155 pooled and 0.246 on RBD.

### 3.7 Supporting: confidence scores do not rank affinity

All 55 binders carry a measured K_D spanning 1.2 nM to 10 µM. Spearman ρ against
log10(K_D): ipTM +0.150, pLDDT +0.139, PAE_interaction −0.164, ESM
pseudo-likelihood +0.226, length −0.155. At n = 55, significance requires
|ρ| > 0.27; none reaches it.

> Frame as confirmation, not discovery, and cite.

### 3.8 Supporting: orthogonal validation by molecular dynamics

37 complexes across three targets, 2 ns implicit-solvent MD, panel spanning the
full ipTM range by construction. ipTM predicts binder RMSD relative to the
target (ρ = −0.486, p = 0.001). Contacts retained and centre-of-mass drift show
nothing — 0 of 37 complexes dissociated, so those measures have no variance at
this timescale and are uninformative rather than null.

> The only result here not measured by AlphaFold against AlphaFold.

---

## 4. Discussion

**Recommendations, in order of what they save.**

1. **Spread folds across many starting points, uniformly, and stop when the
   best-so-far stops improving.** Do not search, do not allocate, do not steer —
   an oracle loses to uniform. Roots are worth ~4× draws fold-for-fold, and the
   ceiling is the best root you drew, not the number you drew: once the pool
   contains a good root, further draws on it do not help and a pool with no good
   root cannot be rescued by depth (§3.4).
2. **Never screen without the target's alignment.** Compute it once per
   campaign; it is the single largest quality lever measured here (§3.1).
3. **Fold at 1 model / 3 recycles to triage, not to rank.** The extra models and
   recycles buy nothing for separating binders from non-binders, and four times
   the compute. They buy a great deal for ranking de novo designs against each
   other, where the two configurations do not agree at all (§3.2).
4. **Use cheap scores as filters and never as estimates.** Worth +0.194 for the
   first purpose and −0.26 for the second (§3.1).
5. **Check your target's score skew before designing any selection scheme.** A
   pilot of 16 draws on 10 backbones takes 40 minutes and tells you whether
   max-of-K is estimable at all (§3.5).
6. Draw at least three sequences per backbone (§3.6).
7. Do not treat any of these metrics as an affinity estimate (§3.7).

**Limitations — first in the paragraph, not last.**

- Three targets; the search result is three runs sharing one root, one target,
  one seed. That root ranks **9th of 10** on measured quality, chosen before any
  of this was measured.
- No wet-lab validation of our own designs. Every "accurate" measurement here is
  a *different AlphaFold configuration*, not reality. The MD in §3.8 is the only
  orthogonal check and it is 2 ns.
- 54 binders bounds every null. "No detectable relationship" is not "no
  relationship".
- The `partial_T` sweep uses 6 parents per T and detects only large
  correlations; ~30 per T would be needed for ρ = 0.5 at 80% power.
- Allocation results are simulated by resampling observed pools, so no policy
  can exceed a root's observed maximum. Comparisons between multi-root policies
  are fair; the single-root arm is a floor.
- One consumer GPU constrained which configurations could be tested, not only
  how many.
- **Our diagnostic claims have a worse track record than our measurements.**
  Three mechanisms we found compelling during this work did not survive their
  own tests: that both broken search runs were breadth-first (they were the
  opposite), that a variance decomposition explained the shortfall (its one
  pre-committed prediction missed), and that search gain concentrates in poor
  starting points (the sweep found the reverse). The descriptive numbers held.
  We report this because it is the reason to trust the tables and not the
  stories.

**Future work.** Every negative result here is bounded by one consumer GPU.
Proteina-Complexa (arXiv:2603.27950) reports MCTS beating best-of-N on hard
targets at far larger budgets, so this does not generalise upward. The
tail-statistic argument does predict the shape of any such win: it should come
from sampling breadth, not from steering. The decisive experiment is a validated
value function, which requires measured binding data rather than more compute.

---

## 5. Data and code availability

github.com/PRANAVV1107/protein-calibration-kit — all scripts, fold outputs,
analysis code, the two search pre-registrations, a running log of the
experiments as they were performed, and a log of errors found and corrected
during the work.

---

## Notes on writing this

> Delete before submitting.

- **One thesis, four measurements, escalating.** Screen → search → allocate →
  why. Resist adding findings that dilute it; §3.6–3.8 are already at the limit.
- **The oracle result is the hook.** Lead the abstract and the talk with it.
- **The error log is an asset.** A random-sequence stub that passed three
  validation stages, a label inversion, a positive control scored as a design, a
  sampling bias from fold ordering, a process guard that matched itself five
  times, and an invalid instrument comparison that stood for three search runs.
  Reporting them is what makes the rest credible.
- **Cite what you confirm.** For the affinity result and for ipSAE, say
  "consistent with" and cite.
- **Figures carry the paper.** Five good ones; Figure 4 (the oracle) is the one
  people will remember.
