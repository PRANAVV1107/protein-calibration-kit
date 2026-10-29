# Where cleverness fails in AlphaFold-filtered binder design

Measuring whether it pays to be clever about *where* to spend AlphaFold
predictions in a de novo binder design pipeline — screening cheaply, searching
over promising regions, allocating toward promising starting points.

**It does not. Every strategy tested loses to, or ties, spending the same
compute uniformly.** One mechanism explains all of them: the objective is a
**maximum** over designs, a tail statistic, while every clever strategy
optimises a **mean**.

All of it on one RTX 3060 (6 GB) under WSL2. ~6,700 AlphaFold2 predictions.

---

## The one result

> **Spend a fixed fold budget on more starting points, not more draws per
> starting point — roots are worth about 4× draws, fold for fold. The ceiling
> is set by the best root in your pool, not by how many you draw.**

Measured directly across 50 backbones × 4 draws (200 folds): at every matched
budget, splitting folds across more roots beats drawing more sequences from
fewer roots, and the gain from adding roots only stops once the pool's best
design has already been found.

A sharper, less stable version of the same idea: an oracle given the single
best-mean starting point is capped at that point's own ceiling, and on one
draw of ten roots it lost to spreading folds uniformly (0.770 against 0.869
expected best ipTM). We pre-registered a re-test before leaning on that number
— bootstrapping 10,000 such draws from the full 50-root pool — and it holds in
only 19.6% of draws (median gap: zero). The mechanism is real; that particular
number was one favourable draw. Both the claim and its own re-test are in the
repo. [ORACLE_BOOTSTRAP_RESULT.md](reports/ORACLE_BOOTSTRAP_RESULT.md)

## The four measurements

| strategy | result | report |
|---|---|---|
| **screen cheaply** | selects **worse than random** (−0.162 ipTM); the *target* fails to fold without an alignment (pLDDT 27.2 vs 84.3) | [RESULTS.md](reports/RESULTS.md) |
| **search over trajectories** | **ties** 60 folds of blind sampling (0.870 vs 0.870) across 3 runs, ~2,900 folds and two corrected implementation bugs | [RUN3_RESULT.md](reports/RUN3_RESULT.md), [RESCORE_RESULT.md](reports/RESCORE_RESULT.md) |
| **allocate adaptively** | no bandit policy beats uniform on 3 targets, at any tested budget | [ALLOCATION_RESULT.md](reports/ALLOCATION_RESULT.md) |
| **why** | max-of-K error *grows* with K on a right-skewed target (0.120 → 0.214); the mean decays normally | [NOISE_RESULT.md](reports/NOISE_RESULT.md) |

## What does pay

- **Spread folds across many starting points.** Roots are worth ~4× draws
  fold-for-fold, and the ceiling is the best root you drew, not how many you
  drew — [ROOT_BREADTH_RESULT.md](reports/ROOT_BREADTH_RESULT.md)
- **Compute the target's alignment once per campaign**, then fold at 1 model /
  3 recycles — same AUC for triage at a quarter the compute (§3.2)
- **Use cheap scores as filters, never as estimates.** Worth **+0.194** for the
  first and **−0.26** for the second, on the same number —
  [RESCORE_RESULT.md](reports/RESCORE_RESULT.md)
- **Check your two scorers agree before pooling them.** Ours do not: identical
  designs, Spearman −0.053 — [CALIBRATION_RESULT.md](reports/CALIBRATION_RESULT.md)

## Validation

- **164 designs with measured binding outcomes** from the Adaptyv EGFR
  competition, synthesised and assayed by their submitters.
- **Replicated on 3,760 independent designs across 15 targets** (Overath,
  Rygaard et al. 2025). Our ipTM AUC of 0.627 against their 0.659 on the same
  competition data; their 21 studies span **0.489 to 1.000**, two below chance.
  `scripts/external_validation.py`
- **2 ns implicit-solvent MD** on 37 complexes — the only measurement here not
  made by AlphaFold about AlphaFold (ρ = −0.486, p = 0.001).
- **No design of ours has been tested at a bench.**

## Method

Three pre-registrations, all committed before any data existed, and all three
changed what got reported:

- [MCTS_PREREGISTRATION.md](reports/MCTS_PREREGISTRATION.md) — the search was
  judged against a bar fixed in advance
- [RUN3_PREREGISTRATION.md](reports/RUN3_PREREGISTRATION.md) — two of its three
  mechanism tests failed, one of them retracted mid-run before the outcome was
  known
- [results/rog_prospective/REGISTRATION.md](results/rog_prospective/REGISTRATION.md)
  — a post-hoc feature registered and then tested; it replicated and **missed
  the threshold set in advance**, so the recommendation built on it was removed
- [ORACLE_BOOTSTRAP_PREREGISTRATION.md](reports/ORACLE_BOOTSTRAP_PREREGISTRATION.md)
  — the project's own headline result, re-tested against a rule fixed before
  the bootstrap ran, and **withdrawn as a general claim** when it only held in
  19.6% of resampled draws

**[reports/EXP1_EXP3_LOG.md](reports/EXP1_EXP3_LOG.md) is the running log**,
written as the work happened rather than afterwards. It contains an invalid
instrument comparison that stood for three search runs before being caught, four
mechanisms that failed their own tests, a void experiment archived with an
explanation of why, and a process-guard bug hit five times. It is kept because
the corrections are part of the method.

A pattern worth stating plainly: **the descriptive measurements here held up and
the diagnostic explanations mostly did not.** Trust the tables over the stories.

## Layout

    PREPRINT_DRAFT.md     the 4-page writeup
    reports/SUMMARY.md    index and one-line verdict for every report
    reports/              one document per finding
    scripts/              every analysis, runnable
    results/              manifests, score JSONs, analysis CSVs, logs

Fold payloads (PDBs, PAE matrices, MSAs) are gitignored and regenerable from the
scripts plus the tracked manifests.

## Reproducing

RFdiffusion, ProteinMPNN and ColabFold 1.6.2 in separate conda environments.
Analyses need only `numpy`, `pandas`, `scipy`:

    python scripts/root_allocation.py --target rbd                # roots vs draws
    python scripts/bandit_sim.py --source variance --target rbd   # allocation vs uniform vs oracle
    python scripts/oracle_bootstrap.py                            # the oracle re-test
    python scripts/external_validation.py                         # cross-target check
    python scripts/rescore_analyze.py --dir results/rescore       # filter vs estimator

None of those need a GPU.
