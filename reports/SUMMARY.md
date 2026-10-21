# Summary

Index and one-line verdict for every report in this folder. Written 2026-09-18.

**What this project is.** A de novo binder design pipeline (RFdiffusion →
ProteinMPNN → AlphaFold2) run end to end on one consumer GPU, with every cheap
shortcut in it measured against 365 designs that other groups synthesised and
assayed. **No design of ours has been tested at a bench.**

---

## The reports

| file | what it covers | verdict |
|---|---|---|
| [RESULTS.md](RESULTS.md) | all nine findings, limits, and the error log | the main document |
| [NOISE_RESULT.md](NOISE_RESULT.md) | why the search failed: max-of-K is a tail statistic | **strongest result** |
| [MCTS_PREREGISTRATION.md](MCTS_PREREGISTRATION.md) | success criteria fixed before the search ran | committed before data existed |
| [MCTS_RESULT.md](MCTS_RESULT.md) | the pre-registered search outcome | search lost, cause diagnosed |
| [RUN2_RESULT.md](RUN2_RESULT.md) | does the mean-value fix work? | diagnosis holds, fix does not rescue |

| [ROOT_SYSTEM.md](ROOT_SYSTEM.md) | what a root is, how 200 were made, what varies | **root choice beats search** |
| [PARTIAL_T_RESULT.md](PARTIAL_T_RESULT.md) | is there locality at any diffusion step size? | no trend with T; two regimes found |
| [ALLOCATION_RESULT.md](ALLOCATION_RESULT.md) | should folds be allocated across roots? | **no — uniform beats the oracle** |
| [ROOT_BREADTH_RESULT.md](ROOT_BREADTH_RESULT.md) | roots or draws? 50 backbones x 4 | **roots, by ~4x** |
| [ROG_PROSPECTIVE_RESULT.md](ROG_PROSPECTIVE_RESULT.md) | does the rog predictor replicate? | real, too weak to use |
| [CALIBRATION_RESULT.md](CALIBRATION_RESULT.md) | do the two AF configs agree? | **no — uncorrelated rankings** |
| [EXPERIMENTS_QUEUED.md](EXPERIMENTS_QUEUED.md) | the next four experiments, prioritised | two findings that reorder everything |
| [RESCORE_RESULT.md](RESCORE_RESULT.md) | search designs re-scored on the bar's instrument | **search ties, does not lose** |
| [RUN3_RESULT.md](RUN3_RESULT.md) | corrected UCT: what it changed and did not | selection fixed, outcome unchanged |
| [EXP1_EXP3_LOG.md](EXP1_EXP3_LOG.md) | rescore + partial_T, running log | both complete |

Supporting numeric output lives in `../results/*_RESULT.txt`.

---

## The search result changed on 2026-09-21

The three MCTS runs were scored with a donated alignment at 1 model / 3
recycles; the 0.870 bar was measured at full MSA / 3 models / 8 recycles. That
comparison was invalid in both directions. Sixty designs were re-folded on the
bar's own instrument ([RESCORE_RESULT.md](RESCORE_RESULT.md)):

    best across all three runs, high accuracy   0.870
    blind sampling, 60 folds, same instrument   0.870
    -> dead tie. Not a loss.

**"The search lost to blind sampling" is withdrawn. It tied, at ~17x the fold
budget.** Two further results came out of the same experiment:

- **The cheap screen is a good filter and a bad estimator.** Selected designs
  beat unselected ones by +0.194 (p = 0.0027) — so the objective carries real
  signal — while over-promising by 0.26 at the top of its range. The project had
  been treating a demonstration of the second as evidence against the first.
- **It is one regression line, not a collapse.** delta = -0.694 x hybrid +
  0.282, crossover at 0.406 (p = 0.0001). Designs below 0.406 improve when
  re-measured; above it they decline. Two designs the screen ranked near the
  bottom scored 0.750.

**Third result, 2026-09-21 — and it retires a whole line of work.** The
`partial_T` sweep showed a backbone's score predicts its neighbours' when
backbones are independently generated (rho +1.000, n=6 clustered, p=0.0014) but
not within one lineage (rho +0.185, n=62). That argued for a bandit over roots.
Simulated on existing folds ([ALLOCATION_RESULT.md](ALLOCATION_RESULT.md)), **no
allocation policy beats uniform at any budget on either instrument** — and
spending every fold on the root with the best true mean reaches 0.770 where
uniform across ten roots reaches 0.869.

**Knowing which root is best and exploiting it perfectly is worse than not
choosing at all.** The objective is a maximum over designs; allocation optimises
a mean. Tails are found by sampling many distributions, not the best one harder.
The practical rule is one line: **draw many roots, sample each a few times, take
the max — do not search, do not allocate, do not steer.**


**Fourth result, 2026-09-23 — and it narrows an earlier recommendation.** The
two AlphaFold configurations the project treats as cheap and expensive were
compared on 60 **identical** designs. They do not agree:

    1 model / 3 recycles     mean 0.158   sd 0.104   max 0.580
    3 models / 8 recycles    mean 0.396   sd 0.252   max 0.870
    Spearman -0.053 (p = 0.69)

Thirteen designs clear 0.7 under the expensive configuration; **one clears even
0.3 under the cheap one.** §3.2's "the cheap configuration matches full
accuracy" is true for separating characterised binders from non-binders (AUC
0.627 vs 0.614) and **false for ranking de novo designs** — two weak signals can
each beat chance without agreeing. Half the gap is arithmetic (ColabFold's
`rank_001` under `--num-models 3` is a max-of-3, worth +0.132); the rest is rank
disagreement, which cannot be corrected away.

Separately, the radius-of-gyration predictor was **pre-registered and tested
prospectively** on 150 held-out backbones. It replicated (gap +0.083, p = 0.019)
but missed the 0.10 threshold fixed in advance, so it is reported as a real but
unusable effect rather than a screening recommendation. Both features this
project selected post hoc have now shrunk on honest re-testing.


Separately, [ROOT_SYSTEM.md](ROOT_SYSTEM.md) shows all three searches ran from
**binder_47, 9th of 10** on measured quality, and that spreading the same folds
across 10 roots beats deepening one at every budget tested.

---

## The three claims that lead

Confirmed as not-obviously-published by a literature check on 2026-09-17.

**1. Single-sequence screening selects worse than random.** −0.162 ipTM against
random choice. Traced to the *target* failing to fold without an alignment
(pLDDT 27.2 versus 84.3), so the interface is scored against a disordered
chain. Four-condition dose-response. Against measured binding on 164 designs,
cheap-screen interface metrics fall below chance (ipSAE 0.466) while the same
metrics at full accuracy reach 0.618 — paired gain **+0.155, p = 0.001**.

**2. Models and recycles are nearly free to cut.** Once the target alignment is
reused, 1 model / 3 recycles gives ipTM AUC **0.627** (0.620 length-controlled)
at ~63 s/design, against 0.614 (0.605) for 3 models / 8 recycles at ~240 s.
**A quarter of the compute, no loss.** The alignment does the work.

**3. `max`-of-K is a tail statistic and it breaks value estimation.** On a
right-skewed target (RBD, skew +3.38 under cheap scoring), the error in a
backbone's score *grows* with the number of draws — 0.120 at K=2 to 0.214 at
K=8 — while `mean`-of-K falls 0.071 → 0.044. `max` is correct for choosing a
design to synthesise, `mean` for valuing a region to search. Conflating them
explains every anomaly in the failed search.

## Supporting findings

| | |
|---|---|
| sequence choice dominates backbone choice | 59–78% of variance; best-of-6 gains +0.155 pooled, +0.246 on RBD |
| our scoring ≈ the competition's own pipeline | ΔAUC −0.002, p = 0.979 — underpowered at 54 binders, not equivalence |
| no metric ranks affinity | all \|ρ\| < 0.23 across 3.9 decades of K_D, n = 55 — confirms published commentary |
| ipTM predicts physical interface rigidity | ρ = −0.486, p = 0.001, 37 complexes, 2 ns MD — the only orthogonal check |
| skew tracks target difficulty | PD-L1 −1.98, IL-7Rα −1.83, RBD +0.54; n = 3, illustrative only |
| 7 of 15 candidates cross-react with homologs | decoy panel; abundance tier clean |
| preference tuning does not transfer | DPO 0.844 in-distribution, 0.490 held-out; sampling +0.018, p = 0.372 |
| monomer foldability | clean null, ρ = +0.081, p = 0.27 |

---

## Where the search stands

**Search run 2** — identical to the pre-registered run except nodes are valued
by `mean` rather than `max`. Read at matched budget; full writeup in
[RUN2_RESULT.md](RUN2_RESULT.md).

    blind sampling, 60 folds       0.870   <- the bar
    run 1, max-valued, 1012 folds  0.860
    run 2, mean-valued, 955 folds  0.820

The middle outcome of the three that were pre-registered. The correction does
what NOISE_RESULT.md predicted *to the node values* — children beat their own
parent 40% of the time against run 1's 26%, and a depth signal appeared where
run 1's profile was flat — but the best design found did not improve. Higher
floors, same ceiling. **The preprint's "the fix is free" line in §3.3 has to
become "the fix makes node values estimable, and value estimation was not the
binding constraint."**

Two caveats bound both runs, **both now fixed in `scripts/mcts_design.py` and
neither yet re-run**: the UCT selection compared candidates against their own
parents rather than siblings *and* used an exploration constant ~3× the entire
spread of the values it was added to — so both runs were effectively
breadth-first, which is the simplest account of a depth signal arriving with no
gain in the best design; and `--budget` counted nodes while claiming to count
folds, which is why run 2 was headed for 4000 folds against run 1's 1012.

**Returns to compute, measured.** Pooling 2,104 RBD folds, expected best-of-N
rises 0.727 (N=60) → 0.849 (N=1000) and gains ~+0.01 per doubling thereafter;
0.33% of draws exceed 0.80 and none exceed 0.90. Run 2 shows the same curve from
the inside: 0.780 by its 32nd fold, +0.040 over the next ~900. More compute is
not the constraint — the absence of a signal that predicts which design is good
is, and this project measured that absence three independent ways (value probe
within-backbone ρ −0.042; DPO held-out 0.490; search 0.820–0.860 against blind
0.870).

## Standing limits on everything here

- No wet-lab validation of our own designs
- Three targets; the search result is n = 1 (one root, one seed, one config)
- 54 binders bounds every null — "no detectable relationship" is not "no
  relationship"
- MD is 2 ns: long enough for rigidity, too short for dissociation
- One consumer NVIDIA GPU constrained which configurations could be tested

## Errors found and corrected

Kept because they shaped the results and because the corrections are part of
the method. Full table in [RESULTS.md](RESULTS.md).

A stub emitting `random.choice()` sequences passed three validation stages, one
reaching ipTM 0.615. A PDB was chosen by name similarity (1ILR is not IL-7Rα).
A pandas truthiness bug inverted all 325 non-binder labels. Cetuximab — an
approved therapeutic — was scored as a de novo design. A length confound that
looked decisive at n=52 (AUC 0.812) vanished at full N (0.569). ColabFold folds
shortest-first, so every interim AUC was computed on a biased short tail. And a
speculative JAX memory change stalled every fold and was misdiagnosed three
times before a four-minute A/B settled it.
