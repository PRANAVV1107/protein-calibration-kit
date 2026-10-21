# Prospective test of radius of gyration — REGISTERED 2026-09-22 12:09

Written **before any fold in this experiment existed**. Do not edit below.

## The claim

Radius of gyration of the binder chain predicts per-root design quality:
**more compact backbones score higher.** Measured on batch 0, n = 50,
rho = -0.482, p = 0.0004, surviving Bonferroni over six features
(ROOT_BREADTH_RESULT.md section 4).

## Why it needs this test

`rog` was chosen by inspecting six features on the same data that measured it.
The project has already been wrong this way once: at n = 10 the leading feature
was `contacts` at +0.511, named as "the one worth testing prospectively"; at
n = 50 it was +0.191, nothing. A post-hoc feature is a hypothesis.

## Design

150 backbones in batches 1-3, none of which contributed to the correlation.
Ranked by rog; the 10 most compact and 10 most extended selected.
4 ProteinMPNN draws each, full MSA at 1 model / 3 recycles — the configuration
the original correlation was measured on.

    compact group   rog 10.88 - 11.91   (predicted BETTER)
    extended group  rog 20.53 - 30.39   (predicted WORSE)
    separation      8.62 angstroms

    compact:  b3/binder_32, b1/binder_32, b2/binder_44, b3/binder_23, b2/binder_0, b1/binder_43, b2/binder_6, b3/binder_48, b2/binder_35, b3/binder_0
    extended: b2/binder_15, b3/binder_18, b3/binder_45, b2/binder_41, b2/binder_39, b1/binder_45, b3/binder_21, b2/binder_18, b2/binder_42, b2/binder_31

## Predictions, fixed now

1. **Primary.** Mean best-of-4 ipTM is higher in the compact group than the
   extended group. One-sided Mann-Whitney, alpha = 0.05.
2. **Effect size to beat.** The batch-0 correlation implies a gap of roughly
   **0.15-0.25 ipTM** between these groups. A significant result with a gap
   below 0.10 confirms direction but means the feature is too weak to screen on.
3. **Secondary.** Spearman(rog, best-of-4) across all 20 backbones should be
   negative and of similar magnitude to -0.482.

## Outcomes, decided now

| result | conclusion |
|---|---|
| compact > extended, p < 0.05, gap >= 0.10 | rog is a validated free pre-screen; keep recommendation #2 |
| compact > extended, p < 0.05, gap < 0.10 | direction real, effect too small to act on; demote to a note |
| not significant | rog joins `contacts` as a post-hoc feature that did not replicate; **remove recommendation #2 from the paper** |
| extended > compact | the batch-0 correlation was an artifact; report it as such |

**Declared expectation:** I expect this to replicate but weaker than -0.482,
because the original estimate was selected as the largest of six and is
therefore biased upward by selection. A gap of 0.10-0.15 seems most likely.

## Limits, acknowledged in advance

- One target (RBD), one generator configuration.
- 20 backbones, 4 draws each. Powered only for a large effect.
- Extreme-group design inflates the apparent effect relative to a full-range
  correlation; the secondary test across all 20 is the fairer estimate.
