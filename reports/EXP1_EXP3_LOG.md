# Experiments 1 and 3 — running log and interpretation

Appended as they proceed, newest entry at the bottom. Experiment 1 is running;
Experiment 3 starts automatically when it finishes, via
`scripts/chain_exp1_exp3.sh`.

---

## What these two are testing

**Experiment 1 — re-score search designs at high accuracy.** All three MCTS runs
folded with a donated target alignment at 1 model / 3 recycles ("hybrid"). The
pre-registered 0.870 bar was measured at full MSA / 3 models / 8 recycles. Those
are different instruments and they disagree badly, so the headline result of
~2,900 folds is currently *unknown* rather than null. Re-folds 15 top + 5 random
designs per run on the bar's own instrument.

The random arm is the load-bearing part. With winners only, a real gain cannot
be separated from a flat instrument offset.

**Experiment 3 — `partial_T` sweep.** Whether the search space has any locality:
does a parent's score predict its children's, at T in {2, 4, 10, 20}? Run 3
measured ρ = 0.168 (p = 0.384) at T = 10, but T = 10 was never chosen on
evidence. Also now tests whether gain depends on parent quality, which would
answer the multi-root question without running it.

**Three outcomes for Experiment 1, fixed before the data:**

| outcome | meaning |
|---|---|
| best re-scored design > 0.870 | search beat blind sampling; three reports need rewriting |
| below 0.870, top arm > random arm | the null survives on a valid comparison and is much harder to attack |
| **top arm ≈ random arm** | the hybrid objective does not identify good designs at all — the *objective* failed, not the algorithm |

---

## 2026-09-21 00:03 — Experiment 1 launched

60 designs queued (`results/rescore/to_fold.fasta`), 15 top + 5 random from each
of runs 1–3, deduplicated by sequence.

Two bugs caught at launch, both worth recording:

**The process guard self-matched.** The first launch aborted silently because
the guard was

    grep -cE '...|envs/colabfold/bin/colabfold_batch|[m]d_stability\.py'

`ps` prints the grep's own argv, so the pattern matched itself and the guard
always reported the GPU busy. **Fourth occurrence of this exact bug in this
project**, and written on a line where `[m]d_stability` was already correctly
bracketed. Fixed by bracketing every alternative. Then repeated it a fifth time
minutes later with `pkill -f 'chain_exp1_exp3.sh'`, which matched the shell
command containing that string and killed the wrong process.

**The runtime estimate was 2.5x low.** Quoted 4 hours from 240 s/design. The
measured rate is **~10 min/design** — the remote mmseqs2 MSA fetch plus 3 models
× 8 recycles is far heavier than the search's 1 model / 3 recycles. 60 designs
is about **10 hours**. The chain's `MAX_WAIT` was set to exactly 10 h and would
have timed out just short of the finish; raised to 24 h.

## 2026-09-21 00:09 — first design re-scored, and it collapses

    run1_top_000    hybrid ipTM 0.860  ->  high-accuracy ipTM 0.227
                    ipSAE 0.0133   pDockQ2 0.0143   pLDDT 78.4

This is **run 1's single best design**, the number the whole project reported as
its search result. On the instrument the bar was measured on, it drops **0.63**.

One design; hold it loosely. But note the direction. The instrument offset
measured on *blind* backbones was **+0.18** — high accuracy reads *higher*. This
design moved 0.63 the other way, and that is consistent rather than
contradictory: `run1_top_000` was selected for maximum hybrid score across 1012
folds, which is precisely where selection-on-a-noisy-metric bites hardest.

**This is the winner's curse.** The top design by a noisy measure regresses hard
when re-measured, and the more designs you select from, the worse it gets. The
project has now met this mechanism four separate ways — max-of-K noise growing
with K, top-quartile parents producing worse children, the reliability estimate,
and now this.

If it holds across the top arm, it is the third outcome in the table above: the
hybrid score the search spent 2,900 folds maximising does not identify designs
that survive a better measurement. That would make every search conclusion a
statement about the wrong quantity — and it would also mean the searches were
never really being compared to blind sampling at all.

**What settles it is the random arm.** If random-arm designs show the same
collapse, it is a flat instrument offset and the comparison merely needs
rescaling. If they do *not* collapse while the top arm does, the drop is
selection-driven and the objective is the problem.

Status at 00:13: 2 of 60 folded, no errors, no mmseqs2 rate-limiting, chain
armed and waiting.

## 2026-09-21 00:45 — 13/60. The collapse is bimodal, not an offset

Advancing, ~10 min/design, no mmseqs2 rate-limiting, no errors. Chain armed.
All 13 so far are run 1's top arm — ColabFold folds in FASTA order, so the
random arm (positions 15-19) has not been reached. **Nothing here is yet
interpretable as the primary result.**

    id                hybrid  high-acc     delta
    run1_top_002       0.840     0.870    +0.030   <- best so far, equals the bar
    run1_top_011       0.760     0.810    +0.050
    run1_top_006       0.800     0.790    -0.010
    run1_top_005       0.810     0.780    -0.030
    run1_top_003       0.820     0.760    -0.060
    run1_top_007       0.790     0.700    -0.090
    -----------------------------------------------
    run1_top_009       0.790     0.330    -0.460
    run1_top_010       0.760     0.300    -0.460
    run1_top_012       0.740     0.250    -0.490
    run1_top_000       0.860     0.230    -0.630   <- run 1's reported best
    run1_top_001       0.850     0.230    -0.620
    run1_top_004       0.810     0.160    -0.650
    run1_top_008       0.790     0.140    -0.650

**The designs split cleanly into two groups.** Six hold up (delta -0.09 to
+0.05); seven collapse by 0.46 to 0.65. There is nothing in between. This is not
the smooth +0.18 instrument offset measured on blind backbones — it is bimodal,
and averaging it (mean delta -0.313) hides the structure entirely.

The dispersions say the same thing:

    hybrid    range 0.740 - 0.860    sd 0.036
    high-acc  range 0.140 - 0.870    sd 0.292

**The hybrid metric compressed a 0.73-wide range of true quality into a
0.12-wide band.** Within its top band it cannot tell a 0.87 design from a 0.14
one. That is the direct, design-level version of what section 3.2 measured
indirectly — hybrid a3m scored AUC 0.545 against full MSA's 0.627.

Spearman(hybrid, high-acc) = -0.042, p = 0.89 on these 13. **Do not read that as
a result yet**: all 13 are top-arm with hybrid sd 0.036, and a null correlation
on a range that narrow is exactly what restriction of range produces. The
bimodality and the magnitude of the collapse are the findings here; the
correlation is not.

**One genuinely good number.** `run1_top_002` scores **0.870** at high accuracy
— equal to the blind-sampling bar. So the search did produce at least one design
that holds up. Whether that beats blind sampling depends on the random arm and
on the other two runs.

### What the random arm now has to settle

If random-arm designs show the same bimodal collapse, it is a property of the
instrument pair and the comparison needs rescaling, not rethinking. If they do
**not** collapse while the top arm does, then selection on hybrid score actively
concentrates designs that fail at high accuracy — the objective would be worse
than uninformative, it would be adversarial. Those two readings differ enough
that no conclusion should be drawn before positions 15-19 land.

## 2026-09-21 01:45 — 29/60. The random arm settles it: the collapse is selection, not instrument

Advancing, no rate-limiting. The five random-arm designs landed and they answer
the question the arm was built to answer.

    arm        n   mean delta   collapsed (delta < -0.3)
    random     5      +0.172     0 of 5
    top       24      -0.310    13 of 24

**Random designs go UP by +0.172. Selected designs go DOWN by -0.310.** A swing
of 0.48 between arms that differ only in whether the hybrid score was used to
pick them.

The random arm's +0.172 offset matches the **+0.18** measured independently
across ten blind backbones (EXPERIMENTS_QUEUED.md 0a) — two unrelated estimates
of the same instrument gap agreeing to within 0.01. That is a real consistency
check, and it means the instrument mapping is a clean positive offset. Every
bit of the top arm's collapse is therefore attributable to selection.

**This is the winner's curse, measured.** Not inferred from skew, not modelled
from a reliability estimate — measured directly, by re-scoring selected and
unselected designs from the same runs on the same instrument. It is worth
**0.48 ipTM** here.

### The number that makes the point on its own

    run1_random_019    hybrid 0.140  ->  high-accuracy 0.750

A design the hybrid screen ranked near the *bottom* scores 0.750 when measured
properly. This is the same phenomenon section 3.1 reported from the fast screen
("the single lowest-ranked design of fifty scored 0.870 at high accuracy"), now
reproduced inside the search's own output.

### The primary, redone

    best across all runs, high accuracy   0.870   (run1_top_002)
    blind-sampling bar, same instrument   0.870
    -> dead tie. Not cleared.

Run 2's best re-scores to 0.820. So on a valid like-for-like comparison the
search **matches** blind sampling rather than losing to it — an improvement on
the invalid comparison, but still not a win, and at ~17x the fold budget.

### Is the hybrid score informative at all?

    top arm    high-acc mean 0.478
    random arm high-acc mean 0.320
    Mann-Whitney one-sided p = 0.0784

Directionally yes — selected designs really are better on average. But
**p = 0.078 at n = 24 vs 5 is not significance**, and the top arm's advantage
(+0.158) is smaller than the amount it over-promises (-0.310). The honest
reading: the hybrid score carries weak real signal, and is wildly
over-optimistic at the top of its range.

Spearman(hybrid, high-acc) across all 29 is +0.172 (p = 0.373) — still null, and
still partly restriction of range.

### Caveats, before this gets quoted

- **The random arm is n = 5 and entirely from run 1.** Runs 2 and 3's random
  designs are at FASTA positions 35-39 and 55-59 and have not been folded. The
  +0.172 estimate rests on five designs.
- One of those five (`run1_random_019`, +0.610) dominates the random arm's mean.
  Without it the offset is +0.063. The median delta is +0.040. **The mean is
  carried by one design** — report the median alongside it.
- 31 designs still to fold, including both remaining random arms.

## 2026-09-21 02:45 — 46/60. The hybrid score IS informative. Two earlier reads corrected.

Advancing, no rate-limiting, all three runs now represented.

**Both significance tests flipped as n grew, and they flipped toward the
screen being useful:**

    n=29 (01:45)   top vs random  p = 0.078   Spearman(hybrid,high-acc) p = 0.373
    n=46 (now)     top vs random  p = 0.0034  Spearman(hybrid,high-acc) p = 0.021

    top arm    n=36  high-acc mean 0.541  median 0.595
    random arm n=10  high-acc mean 0.306  median 0.220

The top arm beats the random arm by **+0.235, p = 0.0034**. The hybrid objective
the search maximised **does** identify designs that hold up at full MSA. The
third outcome in the table above — "the objective failed" — is **ruled out**.

I overstated this twice on smaller samples and both need correcting:

- At 13 designs I wrote that the hybrid metric "cannot tell a 0.87 design from a
  0.14 one". That was a restricted-range artifact read as a finding. It can.
- At 5 random designs the offset read +0.172; at 10 it is **+0.116**, and the
  median is **+0.040**. The single +0.610 design I flagged as carrying the mean
  was indeed carrying it.

### What is actually going on: plain regression to the mean

    Spearman(hybrid, delta) = -0.423, p = 0.0034
    delta = -0.638 x hybrid + 0.253
    crossover at hybrid = 0.397

| hybrid band | n | mean delta | mean high-acc |
|---|---|---|---|
| 0.00-0.20 | 8 | **+0.189** | 0.331 |
| 0.20-0.50 | 2 | -0.175 | 0.205 |
| 0.50-0.75 | 8 | -0.219 | 0.506 |
| 0.75-1.00 | 28 | **-0.247** | 0.550 |

Designs scoring below **0.397** on hybrid go up when re-measured; above it they
go down, monotonically. There is no separate "collapse" mechanism and no
adversarial behaviour — it is one regression line. The bimodality I described at
13 designs was the top band of that line sampled sparsely.

### The practical number

Filtering at hybrid >= 0.75 gives mean high-accuracy **0.550** against the
unfiltered random arm's **0.306**. The cheap screen is worth roughly **+0.24
ipTM as a filter** while being **~0.25 too optimistic as an estimate**. Those
are two different claims about the same number and the project has been
conflating them.

### The primary, unchanged

    run 1   hybrid best 0.860  ->  high-acc best 0.870
    run 2   hybrid best 0.820  ->  high-acc best 0.820
    run 3   hybrid best 0.830  ->  high-acc best 0.820   (6 of 20 folded)

    best overall 0.870   bar 0.870   -> dead tie, not cleared

14 designs left, including run 3's random arm.

## 2026-09-21 04:52 — Experiment 1 closed, Experiment 3 running

**Experiment 1: 60/60, complete.** Written up in
[RESCORE_RESULT.md](RESCORE_RESULT.md). Final: best high-accuracy 0.870 against
a bar of 0.870 (tie); top arm beats random arm by +0.194 (p = 0.0027);
delta = -0.694 x hybrid + 0.282, crossover 0.406 (p = 0.0001). No rate-limiting
across 10 hours. The "errors" reported at earlier checks were my own grep
pattern `429` matching `ipTM=0.429` — there were none.

**Experiment 3 started automatically at 03:47** via the chain. 5 of 24 cells,
88 of 480 folds, no errors. The script sweeps T in the outer loop, so all five
completed cells are T = 2. Running ~12 min/cell, so completion around 08:30 —
faster than the 8 h estimate, since each cell is only 20 folds.

### A construction bias in my own analysis, caught before it was reported

The first interim read showed, at T = 2, n = 4:

    rho 0.400 (p = 0.600)   mean gap +0.134   kid>parent 100%

The 100% and the +0.134 are **artifacts, not signal**. `best_child` is the
**max over 4 children of a mean-of-4**; `parent` is a **single mean-of-4**. A
max of four draws beats one draw of the same thing most of the time regardless
of whether the landscape has any locality at all, so both statistics are biased
upward at every T and cannot discriminate between T values on their own.

This is the project's own max-versus-mean finding reappearing — this time in the
analysis code rather than in the search. Only the **correlations** are unbiased,
because a constant upward bias shifts every point equally and leaves rank
correlation alone. The watcher has been instructed not to report the gap as a
finding.

rho = 0.400 at n = 4 is meaningless either way; it needs all six parents per T.

## 2026-09-21 05:53 — 9/24 cells. T=2 complete, and it shows no locality gain

172/480 folds, no errors, advancing at ~12 min/cell. T = 2 is finished with all
six parents; T = 4 is three cells in. Completion around 09:00.

**T = 2, all six parents: rho = 0.257, p = 0.623.**

Against run 3's T = 10 reference of rho = 0.168, p = 0.384 (n = 42). Both null.
At n = 6, rho would have to exceed about 0.83 to reach p < 0.05, so 0.257 is not
distinguishable from zero *or* from 0.168. **There is no evidence so far that
shortening the diffusion step recovers locality** — which is the "structurally
dead" outcome if it holds at T = 4 and 20.

Reporting only the correlations, per the construction bias noted at 04:52. The
mean gap (+0.117) and kid>parent (100%) are artifacts of comparing a max-of-4
children against a single parent value and say nothing about T.

**Gain-vs-parent-quality, first data: rho = +0.257 at T = 2, pooled +0.238
(p = 0.570, n = 8).**

    better-than-median parents  mean gain +0.100
    worse-than-median parents   mean gain +0.087

This runs *opposite* to run 3's within-tree finding, where top-quartile parents
lost 0.063 while average parents gained. Here good and bad parents gain
equally. Two caveats before anyone reads that as a contradiction: n = 8, and
these parents were drawn by a fixed index rule across the batch rather than
being the top of one tree, so the quality range differs. If it holds at full n,
it would mean the regression-to-the-mean effect seen inside a single tree does
not generalise across independently generated roots — and multi-root search
would be worth building after all.

Too early to say. Needs T = 4, 10 and 20.

## 2026-09-21 06:56 — 13/24. T=4 is the highest so far, and probably shouldn't be believed

249/480 folds, no errors, advancing. T = 2 and T = 4 complete at six parents
each; T = 10 has started. Completion around 09:40.

    T        parents     rho        p       (reference: T=10 run 3, rho 0.168, n=42)
    2              6   0.257    0.623
    4              6   0.600    0.208
    10            42   0.168    0.384   <- run 3

**T = 4 gives the highest correlation seen at any step size.** It is also not
significant — at n = 6, rho must exceed about 0.83 for p < 0.05 — and there are
two reasons to hold it loosely beyond that:

1. **Four T values are being tested.** The maximum of four independent null
   correlations at n = 6 is expected to land near 0.6 by chance alone. Seeing
   one value stand out is the null hypothesis's most likely appearance here, not
   evidence against it.
2. **It is not monotone.** The locality story predicts rho falling smoothly as T
   rises. Observed so far: 0.257 at T = 2, 0.600 at T = 4, 0.168 at T = 10. A
   peak in the middle with a dip below it fits noise better than it fits a
   correlation length.

If T = 20 comes in high as well, that is decisive evidence of no structure. If
T = 20 is low and T = 4 holds, a follow-up at T = 4 with more parents would be
worth running — but on its own, at n = 6, this is not a finding.

**Gain-vs-parent-quality is going nowhere**, as expected at this n:

    T = 2    rho +0.257
    T = 4    rho -0.314
    pooled   rho +0.077  (p = 0.812, n = 12)

Opposite signs at adjacent T values, both null. The multi-root question is not
answered yet in either direction.

Correlations only, per the construction bias. The gaps (+0.117, +0.094) and the
100% kid>parent rates at both T values are artifacts of comparing a max-of-4
against a single value.

## 2026-09-21 07:57 — 17/24. The experiment just falsified its own power, from the inside

323/480 folds, no errors, advancing. T = 10 is four parents in and has produced
the number that settles how much any of this can be trusted.

    T        parents     rho        p
    2              6   0.257    0.623
    4              6   0.600    0.208
    10             4   1.000    0.000    <- this sweep
    10            42   0.168    0.384    <- run 3, same T, same pipeline

**The same partial_T gives rho = 1.000 here and rho = 0.168 in run 3.** A
perfect rank correlation on four points occurs by chance with probability
1/4! ~ 4%, so this is an ordinary small-sample accident, not a discovery. But
it is the cleanest possible internal calibration of the sweep: **at n = 4-6, the
estimates this experiment produces are noise**, and the 0.600 at T = 4 that
looked interesting an hour ago should be read the same way.

The design anticipated this — `partial_t_analyze.py` says in its own output that
six parents per T "detects only large correlations" — but seeing the same T
land at 0.168 and 1.000 makes the point far better than the caveat did. With
six parents, the sampling distribution of rho under no correlation is wide
enough to cover essentially the whole range.

**What this means for the experiment's conclusion.** The sweep can still
support one statement: *no T shows a correlation large enough to be detected
with six parents.* It cannot rank the T values against each other, and it cannot
distinguish "no locality" from "locality too weak for n = 6". Run 3's T = 10
estimate at n = 42 remains the only adequately powered number in the project,
and it is 0.168, p = 0.384.

A properly powered version needs roughly **30 parents per T** to detect rho =
0.5 at 80% power — five times this budget, about 2,400 folds and 40 hours. That
is a real number to put in the writeup rather than a vague "underpowered".

**Gain-vs-parent-quality** now pools to rho = +0.412 (p = 0.113, n = 16), with
better-than-median parents gaining +0.146 against worse-than-median +0.058 —
still trending opposite to run 3's within-tree result. Same power problem
applies; per-T values are +0.257, -0.314, +0.800 at T = 2, 4, 10. Sign flips
between adjacent T values are what noise looks like.

Correlations only; the gaps and 100% kid>parent rates remain construction
artifacts.

## 2026-09-21 08:57 — 20/24. The sweep answered a better question than it asked

386/480 folds, no errors. T = 2, 4 and 10 complete at six parents; T = 20 is one
cell in. The per-T correlations remain uninterpretable at n = 6 — that stands
from the 07:57 entry. But investigating *why* T = 10 read 1.000 here and 0.168
in run 3 turned up something the experiment was not designed to find.

### Two regimes, and only one of them has signal

    across independent backbones, clustered   n=6    rho +1.000   Pearson +0.997
    across independent backbones, pooled      n=19   rho +0.732   (pseudo-replicated)
    within one lineage (run 3)                n=62   rho +0.185   p = 0.146

**A backbone's score predicts its neighbours' scores when the backbones are
independently generated. It does not when they are descendants in one lineage.**

The clustered figure is one point per backbone, averaging its 3-4 T-values, so
the pseudo-replication in the pooled n = 19 is removed:

    binder_24   parent 0.135 -> best_child 0.169
    binder_10   parent 0.141 -> best_child 0.186
    binder_38   parent 0.164 -> best_child 0.219
    binder_31   parent 0.177 -> best_child 0.244
    binder_3    parent 0.185 -> best_child 0.296
    binder_17   parent 0.356 -> best_child 0.605

Perfectly monotone across all six. Under the null that ordering has probability
1/6! = **0.0014**.

**This is not restriction of range.** Parent spread is 0.098 in the sweep and
0.087 in run 3 — comparable. The two regimes differ in *what varies*, not in how
much.

### Why this happens, and what it means for search

Within a lineage every node is a re-noised variant of the same original
backbone, so they share whatever makes that backbone good or bad. The only thing
left to vary is measurement noise, and noise does not correlate with anything —
hence rho = 0.185.

Across independent backbones the variation *is* backbone quality, and that
quality survives re-diffusion. Hence rho = 1.000.

**A tree search descending within one lineage operates precisely in the regime
where there is no signal to exploit.** A search spanning independent roots
operates in the regime where there is. That is a mechanistic argument for
multi-root search, and it is the opposite of the conclusion I gave when
Experiment 4 was deprioritised — I argued then that gain would concentrate in
bad roots and good roots would contribute nothing. The gain-vs-quality data
does not support that either: pooled rho is +0.391, positive, with
better-than-median parents gaining *more* (+0.121 against +0.058).

### Before this gets quoted

- **n = 6 backbones.** Significant by permutation, but six is six, and these
  were drawn by a fixed index rule from one batch.
- Run 3's within-lineage figure is one tree, one root.
- The same backbone re-evaluated with fresh MPNN draws moves by up to **0.335**
  (binder_17: 0.170, 0.392, 0.505 across three runs). Parent values are very
  noisy, which *attenuates* correlation — it cannot manufacture one. A
  correlation this strong despite that noise means the between-backbone signal
  is large.
- The sweep's own question — whether locality depends on T — remains
  unanswerable at this n, and T = 20 is still running.

## 2026-09-21 09:56 — 23/24. All four T values in; the sweep's own question is answered

449/480 folds, no errors. One cell left (t20_p5), running now.

    T        parents     rho        p
    2              6   0.257    0.623
    4              6   0.600    0.208
    10             6   1.000    0.000
    20             4   0.400    0.600
    10            42   0.168    0.384   <- run 3, within one lineage

**rho does not climb as T falls.** The sequence is 0.257, 0.600, 1.000, 0.400 —
non-monotone, with the *highest* value at T = 10, which is the opposite of what
a correlation length predicts. Three of the four are null at n = 6; the T = 10
value is the across-backbone effect already logged at 08:57, not a property of
T = 10 specifically.

**Conclusion for the sweep as designed: no step size tested shows locality that
n = 6 can detect, and there is no trend with T.** `partial_T = 10` was not a
bad choice; no choice was available that would have helped. This closes the
"fixable hyperparameter" hypothesis — the one remaining explanation under which
the three search runs might have been salvageable.

**Gain-vs-parent-quality, final: pooled rho = +0.466 (p = 0.029, n = 22).**
Better-than-median parents gain +0.145 against worse-than-median +0.057 — so
gain rises with starting quality, the opposite of run 3's within-tree result
where top-quartile parents lost 0.063.

Two things temper it. The n = 22 is pseudo-replicated (six backbones at three
to four T values each; effective n = 6). And per-T values are +0.257, -0.314,
+0.886, +0.400 — a sign flip between adjacent T values is what noise looks
like. Taken together with the across-backbone correlation, the reading is that
**good backbones produce good children**, which is a statement about backbone
quality rather than about search.

**It no longer decides anything.** The multi-root question was settled by
[ALLOCATION_RESULT.md](ALLOCATION_RESULT.md): no allocation policy beats
uniform, and uniform beats the oracle. Whether gain rises or falls with parent
quality does not matter if steering toward good parents loses to not steering.

Correlations only; gaps and kid>parent rates remain construction artifacts.

## 2026-09-21 14:55 — root-breadth sweep: 19/200, and a cost figure in the paper looks wrong

Advancing, no mmseqs2 rate-limiting, no errors. 50 backbones x 4 draws, full MSA
at 1 model / 3 recycles, 10-backbone overlap set included for calibration.

**Measured rate: 189 s/design, so ~10.5 h — three times the 3.5 h I estimated.**
That estimate came from section 3.2's "~63 s" for this exact configuration.

This is the second time this session a runtime estimate taken from section 3.2's
cost column has been ~3x low. Experiment 1 ran full MSA at 3 models / 8 recycles
and measured **~600 s/design** against the paper's **~240 s**.

| configuration | section 3.2 | measured here, incl. MSA fetch |
|---|---|---|
| full MSA, 1 model, 3 recycles | ~63 s | **189 s** |
| full MSA, 3 models, 8 recycles | ~240 s | **~600 s** |

**Section 3.2's per-design costs appear to exclude MSA acquisition.** The
benchmark that produced them almost certainly reused precomputed alignments,
while both runs this session fetch a fresh MSA per design from the remote
mmseqs2 server. That also explains the oddity noted earlier — section 3.2 lists
the *hybrid* configuration (donated alignment, no fetch) at 120 s, **slower**
than full MSA at 63 s, which is impossible if the full-MSA figure included a
fetch.

**What this does and does not change.** The *ratio* survives: 189 vs 600 is
3.2x, against the paper's 240/63 = 3.8x, so "a quarter of the compute" is
roughly right and the recommendation stands. The **absolute figures do not** —
anyone budgeting a campaign from that column will be out by 3x.

Before submission, section 3.2's cost column must either be re-measured
end-to-end or labelled explicitly as excluding MSA acquisition, with a note that
a fresh alignment costs roughly 130 s per design on the public server. Added to
the pre-submission checklist rather than silently corrected, since the original
measurement may have been right for its own context.

Analysis needs 40 folds; next check should have it, including the instrument
calibration on the overlap set.

## 2026-09-21 17:00 — 56/200. The overlap calibration fails, and I built it wrong

Advancing steadily (~186 s/design, ~7.7 h left), no rate-limiting, no errors.
Six of the ten overlap backbones have folded, and the calibration they produce
is alarming:

    backbone    max new    max var     delta
    binder_5      0.190      0.860     -0.67
    binder_13     0.140      0.810     -0.67
    binder_0      0.140      0.730     -0.59
    binder_10     0.150      0.540     -0.39
    binder_14     0.160      0.460     -0.30
    binder_9      0.830      0.770     +0.06

    n=6   offset on max -0.427   Spearman +0.145 (p = 0.784)

Offset of nearly half a unit and **no rank agreement at all**. Taken at face
value this says the sweep cannot be compared to the existing ten-root numbers.

**But the comparison is confounded three ways, and two of them are my fault.**

1. **Different K.** The new pools have 2-4 draws so far; the variance study has
   6. `max` is a tail statistic — this project's own central finding — so
   max-of-2 is systematically below max-of-6 on a right-skewed target by
   a large margin. This will shrink as the remaining draws land, but it alone
   could account for much of the offset.
2. **Different sequences.** I re-sampled 4 fresh ProteinMPNN draws (seed 37)
   rather than re-folding the variance study's own sequences. The overlap
   backbones therefore share only the *backbone*, not the designs. So this
   measures "does this backbone score similarly under a fresh sample *and* a
   different configuration" — it cannot isolate the instrument.
3. **Different configuration**, which is the only thing it was meant to measure.

**A calibration that cannot separate its own variables is not a calibration.**
I described this overlap set as the check whose absence invalidated three search
runs; as built, it does not do that job. The honest statement is that the sweep
currently has **no instrument calibration**, not that the instruments disagree.

**The fix is cheap and specific.** Re-fold the variance study's own 60 RBD
sequences at 1 model / 3 recycles. Same designs, same K, only the configuration
differs — which is what a calibration means. 60 folds at ~190 s is about
3.2 hours, and it can run after the sweep rather than competing with it.

**What is still usable meanwhile.** Everything internal to the sweep: the
R-versus-K comparison, the between-root distribution, and the feature
correlations all live on one instrument and are unaffected. Only statements
that pool this sweep with variance-study numbers, or compare it to the 0.870
bar, are blocked.

Early R-versus-K at 15 backbones, equal budget of 10 folds: R=10,K=1 gives
0.474 against R=5,K=2 at 0.416. Breadth ahead by +0.058, same direction as
before and still inside the regime the ten-root result already covered.

## 2026-09-21 18:00 — 75/200. Calibration failure diagnosed: ColabFold's default reporting is a max

Advancing (~193 s/design, ~7 h left), no rate-limiting, no errors. Seven overlap
backbones folded; offset stable at **-0.447 on max, Spearman +0.148 (p = 0.75)**.

Last entry blamed the offset on different K (max-of-2 versus max-of-6). **That
was wrong** — the offset on *means*, which is K-insensitive, is **-0.227**, and
the max offset grew slightly as draws accumulated rather than shrinking. Three
checks now decompose it properly.

### 1. ColabFold's rank_001 is a maximum over models, worth +0.132

`variance_study` ran `--num-models 3` and reports `rank_001`, the best-scoring
of the three. `root_breadth` runs `--num-models 1`. Computing both from the
variance study's own per-model files (60 designs, rank_001/002/003 all on disk):

    rank_001 (best of 3 models)   0.396
    mean across the 3 models      0.264
    inflation from taking the max +0.132

**This accounts for 58% of the -0.227 mean offset**, and it is this project's
central finding appearing inside ColabFold's default output. `--num-models N`
with `rank_001` reported is a **max-of-N**, a tail statistic, and it inflates
scores accordingly. Anyone comparing runs at different `--num-models` is
comparing different estimators, not the same one at different precision.

### 2. The sequences also differ

Both runs used identical ProteinMPNN settings — temp 0.2, seed 37, alanine bias
-1.0 — so I expected the first four draws to match. They do not: **4 of 40
position-matched sequences are identical, and only 1 of 10 backbones has its
breadth draws as a subset of its variance draws.** `variance_prepare.py` selects
a seeded-random backbone subset and folds three targets in one MPNN call, so the
RNG state at each backbone differs. Same seed, different stream, different
sequences.

### 3. What remains

    total mean offset            -0.227
    explained by max-of-3-models +0.132   (58%)
    remainder                    ~-0.095   recycles (3 vs 8) + different sequences

**The calibration as built cannot be repaired by waiting.** Two variables are
confounded with the one being measured. The fix remains: re-fold the variance
study's own 60 RBD sequences at 1 model / 3 recycles — same designs, same K,
only the configuration varying. ~3.2 h, runnable after the sweep.

### This changes something in the paper

Section 3.2 reports that 1 model / 3 recycles "matches" 3 models / 8 recycles at
AUC 0.627 against 0.614. **That is a statement about ranking within each
configuration, not about scale.** The two differ by roughly 0.23 on means,
mostly because one reports a max over three models. Both claims are true and
they are easy to confuse — the paper should say so explicitly, since a reader
who cuts to 1 model expecting the same numbers will see everything drop by
about a quarter of the scale.

Internal comparisons within the sweep are unaffected; only cross-dataset
statements are blocked.

---

## Interpretation guide for whoever reads this next

Three numbers decide what Experiment 1 means, in this order:

1. **Random-arm offset.** `high-acc mean − hybrid mean` on the 15 random
   designs. This is the instrument mapping, uncontaminated by selection. The
   blind-backbone estimate was +0.18; if the random arm reproduces that, the
   instruments differ by a constant and the top arm's behaviour is selection.
2. **Top arm versus random arm** (Mann–Whitney, one-sided). If the top arm is
   not higher, the hybrid objective carries no information about high-accuracy
   quality. That is the most consequential possible outcome and it invalidates
   more than the search.
3. **Best high-accuracy ipTM against 0.870.** The primary, now like-for-like.

Do not read the top arm's raw drop as the headline. It is the *difference*
between the arms that carries the meaning, which is why the random arm exists.
