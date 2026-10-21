# Prospective test of radius of gyration

**Complete 2026-09-23. 20 backbones × 4 draws = 80 folds, exactly as
registered.** Judged against
[REGISTRATION.md](../results/rog_prospective/REGISTRATION.md), written before
any fold existed and not edited since.

**Verdict: the middle outcome. The effect is real and replicated, and it is too
small to screen on.**

---

## The claim under test

Radius of gyration predicts per-root design quality — compact backbones score
higher. Measured on batch 0 at ρ = −0.482, p = 0.0004, n = 50, surviving
Bonferroni over six features ([ROOT_BREADTH_RESULT.md](ROOT_BREADTH_RESULT.md)).

It was selected by inspecting six features on the data that measured it, which
is the situation this project has already been wrong in once: at n = 10 the
leading feature was `contacts` (+0.511, named as the one to test
prospectively) and at n = 50 it was +0.191, nothing.

## Result

150 backbones in batches 1–3, none of which contributed to the correlation.
Ranked by rog; the 10 most compact against the 10 most extended, groups
separated by 8.62 Å with no overlap.

    compact    n=10   mean best-of-4  0.310   median 0.270   max 0.600
    extended   n=10   mean best-of-4  0.227   median 0.180   max 0.500

    PRIMARY    gap +0.083    Mann-Whitney one-sided  p = 0.0187
    SECONDARY  Spearman(rog, best-of-4) -0.374  (p = 0.104, n = 20)

## Against the registered decision rule

| registered outcome | fires? |
|---|---|
| p < 0.05 **and** gap ≥ 0.10 → validated pre-screen | no |
| **p < 0.05 and gap < 0.10 → direction real, too small to act on; demote to a note** | **yes** |
| not significant → did not replicate; remove recommendation | no |
| extended > compact → batch-0 correlation was an artifact | no |

**`rog` replicates.** It is a real property of the design space, not a
selection artifact. But at a gap of 0.083 between the two extremes of a
150-backbone distribution, it cannot carry a screening decision — and the
extreme-group design flatters it, since a full-range correlation would be
smaller still.

The attenuation from −0.482 to −0.374 is what selection bias predicts: the
original estimate was the largest of six and is inflated. The declared
expectation in the registration was a gap of 0.10–0.15; the observed 0.083 is
just below that, so the prediction was slightly optimistic in the direction the
bias predicts.

## A void first attempt, recorded

The first run of this experiment folded **16 backbones at 1–2 draws each, 20
folds**, against a registration specifying 20 × 4 = 80. Every batch's
ProteinMPNN output holds one sequence per backbone, and the prepare script took
what existed instead of failing. Archived at
`results/rog_prospective/VOID_k1_run/` with a note.

That run read **gap +0.039, p = 0.064**, and was briefly reported as "rog failed
to replicate". It was retracted once the fold count was checked. At the
registered K = 4 the same comparison reads +0.083, p = 0.019.

**Best-of-1 is a much noisier statistic than best-of-4, and the noise
attenuated a real effect into a null.** This is the project's own tail-statistic
finding operating on its own experiment, and it is the second time a premature
read was reported and withdrawn. The prepare script now refuses to run when the
registered protocol is not satisfiable.

## What this changes

- **Recommendation #2 in the discussion — "rank candidate backbones by radius of
  gyration before folding" — comes out**, replaced by a note in §3.4 that rog
  predicts root quality weakly (ρ = −0.374 prospectively) and is not strong
  enough to screen on.
- The project now has **two post-hoc features tested prospectively**: `contacts`
  did not replicate at all, `rog` replicated at roughly 60% of its selected
  effect size and below the actionable threshold. Together they are a cleaner
  methodological point than either alone.

## Limits

- **Extreme-group design.** Sampling the tails of the rog distribution
  (10.9–11.9 against 20.5–30.4, batch-0 median 15.3) makes both groups atypical
  backbones and inflates the apparent contrast relative to a full-range
  correlation. A correlation on 40 random held-out backbones would have been the
  better design at identical cost. That was my error in constructing the test
  and it is why the secondary (ρ = −0.374 across all 20) is the fairer number.
- 20 backbones, one target, one generator configuration.
- Scored at 1 model / 3 recycles, the configuration
  [CALIBRATION_RESULT.md](CALIBRATION_RESULT.md) shows ranks de novo designs
  differently from 3 models / 8 recycles. The test is internally consistent —
  the batch-0 correlation was measured the same way — but whether rog predicts
  quality under the expensive configuration is untested.
