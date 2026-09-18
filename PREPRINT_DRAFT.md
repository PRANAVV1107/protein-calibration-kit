# Preprint draft — working document

**Format target:** short methods paper, ~4 pages + figures (MLSB / Bioinformatics
Applications Note length). Everything below is real measured data from this
repository unless marked `[PENDING]`.

**Convention notes are in blockquotes like this** — they explain why a section is
written the way it is. Delete them before submitting.

> **Restructured 2026-09-18** after a literature check. The previous draft led
> with a variance decomposition whose premise — "roughly eight sequences per
> backbone is standard" — could not be confirmed; published sources describe
> 1–5, and some pipelines use one. That finding survives as supporting evidence
> but can no longer carry the paper. The three claims that survived the check
> lead instead, and they share a theme: every cheap shortcut in this pipeline
> fails for a specific, measurable reason, and two of the three have free fixes.

---

## Title

Working options, strongest first:

1. **Where cheap shortcuts fail in de novo binder design, and how to fix two of them**
2. Three failure modes of low-cost screening in AlphaFold-filtered binder design
3. Tail statistics break value estimation in test-time search over binder backbones

> Option 1 promises a practitioner something actionable. Option 3 is narrower
> but more striking to an ML audience — worth splitting into a second short
> paper if the first grows past four pages.

---

## Abstract

> ~150–200 words. Numbers are what make people read on.

De novo binder design pipelines generate backbones with RFdiffusion, assign
sequences with ProteinMPNN, and filter with AlphaFold2. Each stage admits a
cheaper approximation, and we measure where three of them fail. First,
single-sequence AlphaFold screening — a common triage step — selects *worse than
random* (−0.162 ipTM); we trace this to the **target** failing to fold without an
alignment (pLDDT 27.2 versus 84.3), so the interface is scored against a
disordered chain. Second, on 164 experimentally characterised designs, reusing a
precomputed target alignment at 1 model and 3 recycles matches full-accuracy
scoring (ipTM AUC 0.627 against 0.614) at **one quarter the compute** — models
and recycles are nearly free to cut once the alignment is present. Third, in a
pre-registered test-time search over denoising trajectories, scoring nodes by the
*maximum* of K sequence draws causes estimation error to **increase** with K on a
right-skewed target (skew +3.38); the mean decays normally. That single error
explains why the search lost to blind sampling at 17× the compute. All
computation ran on one consumer GPU.

---

## 1. Introduction

> Short. Three paragraphs.

**¶1 — the pipeline.** RFdiffusion generates backbones; ProteinMPNN assigns
sequences; AlphaFold2-multimer scores the complex and top-scoring designs go
forward. Cite Watson 2023, Dauparas 2022, Evans 2021, Bennett 2023.

**¶2 — the gap.** Each stage has a cheaper variant adopted for throughput:
screening without an MSA, fewer models and recycles, and test-time search in
place of exhaustive sampling. Whether these approximations preserve the ranking
they are used to produce is largely unmeasured.

**¶3 — contribution.** We measure three, against 164 designs with experimental
binding outcomes and a pre-registered search protocol.

---

## 2. Methods

- **Targets.** SARS-CoV-2 Spike RBD (6M0J chain E, 333–526, 194 aa); PD-L1
  (5O45 chain A, 18–132, 115 aa); IL-7Rα (3DI3 chain B, 17–209, 193 aa).
- **Ground truth.** 365 designs from the Adaptyv EGFR competition, synthesised
  and assayed by their submitters; 55 binders, 310 non-binders. A balanced
  subset of 164 (54/110) was folded under every configuration compared. One
  positive control (Cetuximab, an approved anti-EGFR therapeutic) was excluded.
- **Scoring.** ColabFold 1.6.2, AlphaFold2-multimer v3. Configurations compared:
  single-sequence; hybrid (precomputed target alignment, single-sequence
  binder); full MSA at 1 model / 3 recycles; full MSA at 3 models / 8 recycles.
- **Statistics.** ROC AUC via the rank-sum identity. Paired bootstrap resampling
  *designs*, not scores. Because binders in this set are longer than
  non-binders (125 aa against 101), every AUC is also reported stratified
  within length tertiles and after regressing length out.
- **Hardware.** One RTX 3060, 6 GB, under WSL2.

---

## 3. Results

### 3.1 Single-sequence screening selects worse than random

**[FIGURE 1]** target pLDDT against selection quality, four conditions.

Scoring the same 60 RBD complexes under four conditions differing only in how
the target's structure is supplied:

| condition | target pLDDT | selection vs random |
|---|---|---|
| no MSA (cheap screen) | 27.2 | **−0.162** |
| target structure as template | 45.8 | −0.013 |
| precomputed target MSA | 79.2 | +0.141 |
| full MSA (reference) | 84.3 | — |

pLDDT below 50 indicates a disordered prediction. The cheap screen was scoring
interfaces against an unfolded target.

Independently, folding the fast screen's top and bottom deciles at high accuracy
gave means of 0.670 and 0.671 (Mann–Whitney p = 0.970, n = 20). The single
lowest-ranked design of fifty scored 0.870 at high accuracy — the joint best in
that campaign.

Against measured outcomes on the 164-design set, cheap-screen interface metrics
fall *below* chance while the same metrics at full accuracy do not:

| metric | cheap screen | full accuracy |
|---|---|---|
| ipSAE_min | **0.466** | 0.618 |
| pDockQ2_min | **0.461** | 0.623 |

Paired on identical designs, full accuracy beats the cheap screen by **+0.155
on ipSAE (p = 0.001)**.

> The dose–response across four conditions is what turns this from an
> observation into a mechanism. Lead with it.

### 3.2 Models and recycles are nearly free to cut

**[FIGURE 2]** AUC against compute, three configurations.

Same 164 designs, same alignment, varying only model count and recycles:

| configuration | cost/design | ipTM AUC | length-controlled |
|---|---|---|---|
| hybrid a3m, 1 model, 3 recycle | ~120 s | 0.545 | 0.521 |
| **full MSA, 1 model, 3 recycle** | **~63 s** | **0.627** | **0.620** |
| full MSA, 3 models, 8 recycles | ~240 s | 0.614 | 0.605 |

The cheap-model configuration matches full accuracy at a quarter of the compute.
Combined with §3.1, the alignment is doing the work; sampling depth is not.

> The practical recommendation is one sentence: compute the target's alignment
> once per campaign and fold at 1 model / 3 recycles. Readers want a number to
> act on.

### 3.3 Tail statistics break value estimation in test-time search

**[FIGURE 3]** node-score disagreement against K, two targets, max vs mean.

We pre-registered a guided best-first search over RFdiffusion denoising
trajectories (`diffuser.partial_T`), with nodes evaluated by folding ProteinMPNN
sequences. Success criteria and the comparison baseline were fixed in advance
and committed before any result existed.

**The search lost.** 1012 folds over 253 backbones reached ipTM 0.860, against
0.870 for blind sampling at 60 folds — indistinguishable, at 17× the compute.

The cause is a statistics error. Scoring a backbone by the maximum of K sequence
draws, the disagreement between two disjoint sets of K draws from the *same*
backbone:

| target | K=2 | K=4 | K=8 |
|---|---|---|---|
| IL-7Rα, `max` | 0.141 | 0.068 | **0.034** |
| SARS-CoV-2 RBD, `max` | 0.120 | 0.167 | **0.214** |
| SARS-CoV-2 RBD, `mean` | 0.071 | 0.059 | **0.044** |

RBD's score distribution has skew **+3.38** (median 0.14; 91% of draws below
0.3, 3% above 0.7); IL-7Rα's is −0.66. The maximum is a **tail statistic**, so
on a right-skewed target its value depends on whether a subset happened to catch
an excursion — and adding draws makes two subsets *more* likely to disagree. The
mean decays normally on both.

The search scored nodes by max-of-4 on the right-skewed target. That single fact
accounts for every anomaly observed: parent score failing to predict child
quality, the best-known node yielding nothing on three separate expansions, the
largest gains arriving from low-scoring parents, a flat score-versus-depth
profile, and children beating parents only 26% of the time.

**The fix is free.** `max` is the correct statistic for selecting a design to
synthesise; `mean` is the correct statistic for valuing a region to search. The
same folds serve both.

> Report the pre-registration and the loss plainly. A negative result with a
> diagnosed cause and a costless remedy is worth more than a marginal win, and
> reviewers can tell the difference.

### 3.4 Supporting: sequence choice dominates backbone choice

Within-backbone (sequence) share of total interface-score variance, 180
high-accuracy predictions over 30 backbones, 6 draws each:

| target | ipTM | ipSAE_min | pDockQ2_min |
|---|---|---|---|
| IL-7Rα | 59.0% | 43.4% | 53.4% |
| PD-L1 | 68.9% | 55.1% | 66.4% |
| SARS-CoV-2 RBD | 78.3% | 79.8% | 78.9% |

Best-of-six rather than one raises mean ipTM by 0.155 pooled and 0.246 on RBD.

> Previously the paper's headline. Demoted because the framing rested on "eight
> sequences per backbone is standard", which could not be confirmed — published
> descriptions give 1–5. State the variance result on its own terms and make no
> claim about ablating a convention.

### 3.5 Supporting: confidence scores do not rank affinity

All 55 binders carry a measured K_D spanning 1.2 nM to 10 µM. Spearman ρ against
log10(K_D): ipTM +0.150, pLDDT +0.139, PAE_interaction −0.164, ESM
pseudo-likelihood +0.226, length −0.155. At n = 55, significance requires
|ρ| > 0.27; none reaches it.

> Frame as confirmation, not discovery. That confidence scores work as negative
> filters and not as rankings is already discussed in the field. The
> contribution is measuring it against real K_D across 3.9 decades.

### 3.6 Supporting: orthogonal validation by molecular dynamics

37 complexes across three targets, 2 ns implicit-solvent MD, panel spanning the
full ipTM range by construction. ipTM predicts binder RMSD relative to the
target (ρ = −0.486, p = 0.001). Contacts retained and centre-of-mass drift show
nothing — 0 of 37 complexes dissociated, so those measures have no variance at
this timescale and are uninformative rather than null.

> The only result here not measured by AlphaFold against AlphaFold.

---

## 4. Discussion

**Recommendations.**
1. Never screen without the target's alignment. Compute it once per campaign.
2. Fold at 1 model / 3 recycles; the extra models and recycles buy nothing.
3. In any search or selection loop, value regions by the **mean** of sampled
   scores and select final designs by the **max**. Check the target's score
   skew first — a pilot of 16 draws on 10 backbones takes 40 minutes.
4. Draw at least three sequences per backbone.
5. Do not treat any of these metrics as an affinity estimate.

**Limitations — first in the paragraph, not last.**
- Three targets. The search result is n = 1: one root, one target, one seed,
  one hyperparameter setting.
- No wet-lab validation of our own designs; the only measured outcomes are
  other groups' designs re-scored.
- 54 binders bounds every null here. "No detectable relationship" is not "no
  relationship".
- MD is 2 ns — sufficient for rigidity, too short for dissociation.
- One consumer GPU constrained which configurations could be tested, not only
  how many.

**Future work.** The mean-versus-max correction is untested inside a live
search; that it makes node values estimable does not prove a search using it
succeeds. Proteina-Complexa (arXiv:2603.27950) reports MCTS beating best-of-N on
hard targets at far larger budgets, and re-running with corrected node values on
a low-skew target is the obvious next experiment.

---

## 5. Data and code availability

github.com/PRANAVV1107/protein-calibration-kit — all scripts, fold outputs,
analysis code, the search pre-registration, and a log of errors found and
corrected during the work.

---

## Notes on writing this

> Delete before submitting.

- **One theme, three claims.** Cheap shortcuts fail for measurable reasons; two
  have free fixes. Resist adding findings that dilute it.
- **The error log is an asset, not an embarrassment.** A random-sequence stub
  that passed three validation stages, a label inversion, a positive control
  scored as a design, a sampling bias from fold ordering. Reporting them is
  what makes the rest credible.
- **Cite what you confirm.** For the affinity result and for ipSAE, say
  "consistent with" and cite.
- **Figures carry the paper.** Three good ones beat eight mediocre ones.
