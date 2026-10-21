# Instrument calibration — do the two AlphaFold configurations agree?

**Complete 2026-09-23. 60 designs, identical sequences, folded under both
configurations.** The only difference is `--num-models` and `--num-recycle`.

**Verdict: they do not agree, and the disagreement is not a scale offset. The
cheap configuration compresses its score range and fails to identify the designs
the expensive one rates highest.**

---

## Why this was needed

The root-breadth sweep included ten backbones the variance study had also
folded, intended as a calibration. It could not serve as one: the two datasets
differed in configuration, in sequence set, and in K simultaneously. Three
variables, one measurement.

This re-folds the variance study's **own 60 RBD sequences** at 1 model / 3
recycles. Same designs, same K, only the configuration varies.

## Result

    1 model / 3 recycles     mean 0.158   sd 0.104   range 0.100 - 0.580
    3 models / 8 recycles    mean 0.396   sd 0.252   range 0.130 - 0.870

    offset      -0.238   (median -0.115)
    Spearman    -0.053   (p = 0.69)
    Pearson     +0.160   (p = 0.22)
    fit         cheap = 0.066 x expensive + 0.132

**No rank agreement on identical designs.** The slope of 0.066 is nearly flat:
the cheap configuration's score is almost independent of the expensive one's.

Per backbone, best-of-6 under each:

| backbone | cheap | expensive |
|---|---|---|
| binder_42 | **0.160** | **0.870** |
| binder_5 | **0.210** | **0.860** |
| binder_13 | 0.570 | 0.810 |
| binder_9 | 0.580 | 0.770 |
| binder_0 | 0.440 | 0.730 |
| binder_15 | **0.140** | **0.710** |
| binder_34 | **0.160** | **0.710** |
| binder_10 | 0.130 | 0.540 |
| binder_47 | 0.160 | 0.520 |
| binder_14 | 0.420 | 0.460 |

Per-backbone Spearman +0.326 (p = 0.358, n = 10) — no better.

**The sharpest statement:**

    13 of 60 designs score above 0.7 under the expensive configuration.
    Exactly 1 of those 13 scores even above 0.3 under the cheap one.
    The cheap configuration's maximum over all 60 designs is 0.580.

The cheap configuration cannot reach the scores the expensive one assigns to its
best designs, and does not rank them highly either. `binder_42`, the best
backbone in the set at 0.870, is scored 0.160 — indistinguishable from the worst.

## Where the difference comes from

Part is arithmetic and was measured separately: ColabFold's `rank_001` under
`--num-models 3` is the **best of three models**, a max-of-3, and taking it
inflates scores by a measured **+0.132** relative to the mean across models.
That accounts for roughly half the offset. The rest is recycling depth — these
complexes move substantially between recycle 3 and recycle 8.

But the offset is the less important half. **A constant offset would preserve
ranking and could simply be subtracted.** What is measured here is rank
disagreement, which cannot be corrected away.

## What this changes

**§3.2 of the preprint needs revising.** It reports that 1 model / 3 recycles
"matches" 3 models / 8 recycles at ipTM AUC 0.627 against 0.614, and recommends
cutting to the cheap configuration. Both numbers are correct and the conclusion
drawn from them is too strong.

Two weak signals can each beat chance against binding data without agreeing with
each other. AUC ~0.62 is weak, and that is exactly what has happened: the
configurations rank binders above non-binders about equally well *on the Adaptyv
set*, while ranking **de novo RBD designs** in mutually uncorrelated orders.

The correct claim is narrower: *for separating experimentally-characterised
binders from non-binders, the cheap configuration is as good and four times
cheaper.* It is **not** licensed to say the two configurations are
interchangeable, and it should not be used to rank de novo designs against each
other.

**It also conditions the root-breadth sweep**, which used the cheap
configuration throughout. Its finding that breadth beats depth is internal to
one instrument and survives — the comparison is between allocation strategies
scored the same way. But its *per-root rankings* should not be treated as
identifying which backbones are genuinely good, and
[ROOT_BREADTH_RESULT.md](ROOT_BREADTH_RESULT.md) already carries that caveat.

## Limits

- 60 designs, 10 backbones, one target.
- The cheap configuration's range on this set is narrow (sd 0.104), and
  restriction of range attenuates rank correlation. **That compression is itself
  the finding** — a configuration that collapses its output range cannot
  discriminate, whatever the cause — but the Spearman of −0.053 should be read
  as "no usable signal here", not as a precise estimate of zero.
- Both configurations are AlphaFold. Neither is ground truth. The expensive one
  is treated as the reference because the 0.870 bar and the variance study were
  measured on it, not because it is known to be more accurate.

---

Regenerate from `results/calibration/` and `results/variance_study/rbd/`.
