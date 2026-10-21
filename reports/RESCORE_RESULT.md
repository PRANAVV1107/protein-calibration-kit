# Experiment 1 — search designs re-scored on the bar's own instrument

**Complete 2026-09-21, 60/60 designs.** Written against the three questions
fixed in `scripts/rescore_analyze.py` before the run.

**Verdict: the search ties blind sampling, it does not lose to it — and the
cheap screen is a good filter while being a bad estimator.** Both conclusions
required this experiment; neither was available from the original comparison,
which was invalid in both directions.

---

## Why this had to be run

All three MCTS runs folded with a **donated target alignment at 1 model / 3
recycles** ("hybrid"). The pre-registered bar of 0.870 was measured in the
variance study at **full MSA / 3 models / 8 recycles**. Different instruments,
never reconciled. On ten backbones both had scored, high accuracy read +0.18
higher *and reordered them*.

So "search reached 0.860 against a bar of 0.870" compared a number from the
low-reading instrument to one from the high-reading one. This re-folds 45 top
designs (15 per run, by hybrid ipTM) and 15 random designs (5 per run, drawn
across the whole score range) on the bar's instrument.

**The random arm is what makes the experiment work.** With winners only, a real
effect cannot be separated from a flat instrument offset.

---

## 1. The instrument mapping

Random arm only — designs not selected on either score:

    n = 15
    hybrid   mean 0.190
    high-acc mean 0.329
    offset        +0.139   (median +0.040)
    collapsed (delta < -0.3):  0 of 15

**+0.139 against the +0.18 measured independently on ten blind backbones.** Two
unrelated estimates of the same gap, agreeing within 0.04. The instrument
mapping is a clean positive offset.

The mean is skewed by two designs that gained +0.61 and +0.63; the **median
+0.040** is the more honest summary, and both are reported wherever this is
quoted.

## 2. The primary, redone

| run | hybrid best | high-accuracy best | vs 0.870 |
|---|---|---|---|
| run 1 (max, broken UCT) | 0.860 | **0.870** | **tie** |
| run 2 (mean, broken UCT) | 0.820 | 0.820 | −0.050 |
| run 3 (mean, fixed UCT) | 0.830 | 0.820 | −0.050 |

    best across all runs   0.870
    blind-sampling bar     0.870
    -> dead tie. Bar not cleared.

**The pre-registered conclusion changes from "search loses" to "search ties".**
It still does not win, and it still costs ~17x the folds of the 60-fold blind
baseline. But the three reports that state the search *lost* are overstating a
comparison that was never valid.

## 3. Did the search optimise the right thing? Yes.

    top arm    n=45  high-acc mean 0.523  median 0.550
    random arm n=15  high-acc mean 0.329  median 0.230
    Mann-Whitney one-sided  p = 0.0027

    Spearman(hybrid, high-acc) over all 60:  +0.320  (p = 0.013)

**The hybrid objective does identify designs that hold up at full MSA**, by
+0.194 on the mean. The worst available outcome — that the search spent 2,900
folds maximising a quantity unrelated to design quality — is ruled out.

---

## The mechanism: one regression line, not a collapse

    Spearman(hybrid, delta) = -0.478,  p = 0.00011
    delta = -0.694 x hybrid + 0.282
    crossover at hybrid = 0.406

| hybrid band | n | mean delta | median delta | mean high-acc |
|---|---|---|---|---|
| 0.00–0.20 | 11 | **+0.206** | +0.100 | 0.345 |
| 0.20–0.50 | 4 | −0.048 | −0.110 | 0.282 |
| 0.50–0.75 | 10 | −0.231 | −0.220 | 0.497 |
| 0.75–1.00 | 35 | **−0.260** | −0.230 | 0.531 |

Designs scoring below **0.406** on hybrid go *up* when re-measured; above it
they go *down*, monotonically. There is no separate collapse mechanism and
nothing adversarial — it is textbook regression to the mean, and it is highly
significant (p = 0.0001).

**Two designs make the point on their own:**

    run3_random_055   hybrid 0.120  ->  high-accuracy 0.750
    run1_random_019   hybrid 0.140  ->  high-accuracy 0.750

Both were ranked near the bottom by the cheap screen. Both are among the best
designs in the set. This reproduces, inside the search's own output, what
section 3.1 found in the fast screen — the single lowest-ranked design of fifty
scored 0.870 at high accuracy.

## The distinction the project has been conflating

| | value |
|---|---|
| cheap screen as a **filter** — top arm vs random arm | **+0.194**, p = 0.0027 |
| cheap screen as an **estimator** — over-optimism at hybrid > 0.75 | **−0.260** |

Filtering at hybrid ≥ 0.75 raises mean high-accuracy quality from 0.475
(unfiltered) to 0.531, and the fraction above 0.7 from 30% to 40%.

**The screen is worth using and its numbers are not worth believing.** Those are
two different claims about the same quantity, and this project has repeatedly
treated a demonstration of the second as evidence against the first.

---

## What this changes

- **[MCTS_RESULT.md](MCTS_RESULT.md), [RUN2_RESULT.md](RUN2_RESULT.md),
  [RUN3_RESULT.md](RUN3_RESULT.md)** — "the search lost" becomes "the search
  tied, on a comparison that had to be rebuilt to be meaningful".
- **[PREPRINT_DRAFT.md](../PREPRINT_DRAFT.md) §3.3** — the search narrative
  needs the instrument mismatch stated, and the tie reported instead of a loss.
- **§3.1 and §3.2 are strengthened, not weakened.** This is a third independent
  demonstration that cheap screening scores are miscalibrated at the top, now
  with a fitted regression line and a crossover point.

## Limits

- 60 designs, 20 per run; the random arm is 15.
- One target (RBD), one root (binder_47) — every design here descends from it.
- ipTM only. ipSAE and pDockQ2 were recorded but not analysed here.
- "High accuracy" is itself AlphaFold. This measures agreement between two
  AlphaFold configurations, not agreement with reality. No design in this set
  has been tested at a bench.
- The regression line is fitted on designs spanning hybrid 0.12–0.86; the
  crossover at 0.406 should not be extrapolated outside that.

---

Regenerate with:

    python scripts/rescore_analyze.py --dir results/rescore

Raw per-design table: `results/rescore/rescore_results.csv`.
