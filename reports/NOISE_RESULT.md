# Why the search failed: max-of-K is the wrong statistic

320 folds — 10 backbones x 16 ProteinMPNN draws x 2 targets, folded at the
search's own settings. Run to answer one question: how many draws does a
backbone need before its score means anything?

## The measurement

Disagreement between two disjoint sets of K draws from the **same** backbone,
scored by `max`:

| target | K=2 | K=4 | K=8 | vs K=2 |
|---|---|---|---|---|
| IL-7Rα | 0.141 | 0.068 | **0.034** | 0.24x |
| SARS-CoV-2 RBD | 0.120 | 0.167 | **0.214** | **1.78x** |

IL-7Rα behaves as expected — noise falls faster than 1/sqrt(K), so more draws
fix it. **RBD gets worse.** Neither a 1/sqrt(K) decay nor a plateau; noise
nearly doubles from K=2 to K=8.

## Why

The two targets have very different score distributions:

| | skew | median | frac < 0.3 | frac > 0.7 |
|---|---|---|---|---|
| IL-7Rα | −0.66 | 0.74 | 21% | 52% |
| **RBD** | **+3.38** | **0.14** | **91%** | **3%** |

RBD is a spike at 0.14 with a rare tail reaching 0.75+. `max` is a **tail
statistic**: its value depends on whether a sample happened to catch an
excursion. Adding draws raises the chance that one subset catches the tail
while the other does not, so disagreement *grows* with K. On IL-7Rα, where
half of all draws already exceed 0.7, `max` behaves like a location statistic
and converges normally.

## The fix is the statistic, not the budget

Same test, scoring by `mean` instead:

| target | K=2 | K=4 | K=8 |
|---|---|---|---|
| IL-7Rα | 0.170 | 0.125 | **0.089** |
| RBD | 0.071 | 0.059 | **0.044** |

Monotone decay on both. Node values become estimable at a cost that scales
predictably.

## What this explains

The pre-registered search scored nodes by **max of 4 draws** on RBD — skew
+3.38, the worst possible case for that statistic. Every node value was
effectively a lottery ticket on whether its four draws caught a tail event.

That accounts for every anomaly in that run, none of which needed a separate
explanation:

- parent score did not predict child quality
- expanding the best-known node produced nothing above 0.27, three times
- the largest gains came from expanding *weak* nodes
- the depth profile was flat
- children beat parents only 26% of the time

UCT was not misled by a weak signal. It was reading tail noise and propagating
it as value.

## The distinction that was missed

`max` and `mean` answer different questions, and conflating them was the bug:

- **Selecting a design to synthesise** — you want the best one, so `max` is
  correct
- **Valuing a region to search** — you want to know whether the region is
  good, which is `mean`

The search used the design-selection statistic as a search value. On a
left-skewed target that happens to work; on a right-skewed one it cannot.

## What to do differently

1. **Score nodes by mean, select finalists by max.** Free — same folds, two
   statistics.
2. **Check the target's score distribution before searching.** Skew is
   computable from a pilot of ~16 draws on a handful of backbones, which is
   40 minutes. A target with skew above ~2 needs a different value statistic.
3. **Expect targets to differ.** RBD and IL-7Rα required opposite treatments
   from the same pipeline. A method validated on one says little about the
   other.

## Limits

- Two targets, 10 backbones each, 16 draws. The skew estimate per backbone
  rests on 16 points.
- `mean` was evaluated only on this resampling test, not inside a live search.
  That it is estimable does not prove a search using it would succeed.
- ipTM only. Whether ipSAE or pDockQ2 are as heavy-tailed is untested.
- 1 model / 3 recycles with a donated target MSA throughout.
