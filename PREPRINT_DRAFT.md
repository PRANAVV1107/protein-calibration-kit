# Preprint draft — working document

**Format target:** short methods paper, ~4 pages + figures (MLSB / Bioinformatics
Applications Note length). Everything below is real measured data from this
repository unless marked `[PENDING]`.

**Convention notes are in blockquotes like this** — they explain why a section is
written the way it is. Delete them before submitting.

---

## Title

Working options, strongest first:

1. **Sequence choice dominates backbone choice in de novo binder design**
2. Where does the variance come from? Sequence versus backbone in AlphaFold-filtered binder design
3. How many sequences per backbone? A variance decomposition for de novo binder design

> Option 1 states the finding. Titles that state a finding get read more than
> titles that state a topic. Avoid a colon-plus-explainer; it reads as filler.

---

## Abstract

> ~150–200 words. Structure: context (1 sentence), gap (1), what you did (1–2),
> what you found (2–3 with numbers), what it means (1). Numbers in the abstract
> are what make people read on.

De novo binder design pipelines generate backbones with RFdiffusion, assign
sequences with ProteinMPNN, and filter with AlphaFold2. Standard practice samples
roughly eight sequences per backbone, a convention adopted without a published
ablation. We decomposed AlphaFold2 interface-score variance into between-backbone
and within-backbone (sequence) components across three targets — SARS-CoV-2 RBD,
PD-L1 and IL-7Rα — using 180 high-accuracy predictions over 30 backbones with six
independent sequence draws each. Within-backbone variance exceeded
between-backbone variance on every target (59.0%, 68.9% and 78.3% of total
variance by ipTM), meaning the sequence draw, not the backbone, dominates the
score a design receives. Drawing the best of six sequences rather than one raised
mean ipTM by 0.155 pooled and by 0.246 on the hardest target, with no additional
backbone generation. We further show that single-sequence AlphaFold screening — a
common cheap triage step — selects worse than random choice (−0.162 ipTM), and
trace the cause: without a multiple sequence alignment the *target* fails to fold
(pLDDT 27.2 versus 84.3), so the interface is scored against a disordered chain.

---

## 1. Introduction

> Short. Three paragraphs: what the pipeline is, what's unexamined, what we did.
> Resist the urge to review the whole field.

**¶1 — the pipeline.** RFdiffusion generates backbones; ProteinMPNN assigns
sequences; AlphaFold2-multimer scores the resulting complex and the top-scoring
designs go forward. Cite Watson 2023, Dauparas 2022, Evans 2021, Bennett 2023.

**¶2 — the gap.** Two parameters in this pipeline are set by convention rather
than measurement: how many sequences to draw per backbone (commonly eight), and
whether a cheap low-cost AlphaFold pass is a useful triage step. Neither appears
to have a published ablation. `[verify this claim once more before submission —
it is the novelty claim and the thing a reviewer will check first]`

**¶3 — contribution.** We measure both.

---

## 2. Methods

> Keep to what a reader needs to reproduce it. Full parameters go in supplement.

- **Targets.** SARS-CoV-2 Spike RBD (PDB 6M0J chain E, residues 333–526, 194 aa);
  PD-L1 (5O45 chain A, 18–132, 115 aa); IL-7Rα (3DI3 chain B, 17–209, 193 aa).
- **Design.** RFdiffusion with hotspot conditioning; ProteinMPNN at sampling
  temperature 0.2 with an alanine logit bias of −1.0 (justified in §3.4).
- **Scoring.** ColabFold 1.6.2, AlphaFold2-multimer v3. High accuracy = full MSA
  search, 3 models, 8 recycles. Cheap screen = single-sequence mode, 1 model,
  3 recycles.
- **Variance study.** 10 backbones per target selected at random (seeded, not by
  score, to avoid biasing the variance estimate), 6 independent ProteinMPNN draws
  each, all folded at high accuracy. 180 predictions.
- **Statistics.** Variance decomposed as between-backbone (variance of
  per-backbone means) against within-backbone (mean of per-backbone variances).
  Bootstrap confidence intervals; Spearman for rank correlation; Mann–Whitney for
  group comparisons.

---

## 3. Results

### 3.1 Sequence choice dominates backbone choice

**[FIGURE 1]** stacked bars, within- vs between-backbone variance per target,
with the three metrics side by side.

Within-backbone (sequence) share of total interface-score variance:

| target | ipTM | ipSAE_min | pDockQ2_min |
|---|---|---|---|
| IL-7Rα | 59.0% | 43.4% | 53.4% |
| PD-L1 | 68.9% | 55.1% | 66.4% |
| SARS-CoV-2 RBD | 78.3% | 79.8% | 78.9% |

> Reporting three metrics rather than one is deliberate: it shows the result is
> not an artifact of ipTM's known chain-length sensitivity. Note honestly that
> the effect narrows under ipSAE on the two easier targets — a reviewer will find
> that, so find it first.

Within-backbone spread across six draws averaged 0.228 (IL-7Rα), 0.256 (PD-L1)
and 0.540 (RBD), with a maximum of 0.700 — on geometry held fixed.

### 3.2 The cost of drawing one sequence

**[FIGURE 2]** best-of-N curves, three targets, N = 1…6.

| N draws | IL-7Rα | PD-L1 | RBD | pooled |
|---|---|---|---|---|
| 1 | 0.766 | 0.805 | 0.400 | 0.657 |
| 3 | 0.820 | 0.870 | 0.565 | 0.751 |
| 6 | 0.835 | 0.883 | 0.646 | 0.788 |

Most of the gain arrives by N = 3. On the hardest target a single draw averages
0.400 while the same backbones reach 0.646 at six draws.

> State the practical recommendation explicitly. Readers want a number to act on.

### 3.3 Single-sequence screening selects worse than random

**[FIGURE 3]** target pLDDT against selection quality, four conditions.

Scoring the same 60 RBD complexes under four conditions differing only in how the
target's structure is supplied:

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

> The dose–response across four conditions is what turns this from an observation
> into a mechanism. Lead with it.

### 3.4 Composition collapse is fixed by residue bias, not temperature

| ProteinMPNN config | mean Ala | max Ala | longest run |
|---|---|---|---|
| temp 0.1 (default) | 20.6% | 59.3% | 16 |
| temp 0.2 | 20.4% | 60.2% | 16 |
| temp 0.2 + Ala bias −1.0 | **8.8%** | 15.8% | 5 |

> Frame as confirmation, not discovery — a negative alanine bias is already the
> community remedy. The contribution is the measurement that temperature alone
> does nothing.

### 3.5 Validation against measured binding `[PARTIAL]`

Re-scoring 365 experimentally characterised designs from the Adaptyv EGFR
competition (55 binders, 310 non-binders):

| metric | AUC | 95% CI |
|---|---|---|
| ipTM (our cheap screen) | 0.518 | [0.439, 0.600] |
| ipSAE_min | 0.458 | [0.377, 0.544] |
| ipTM (competition's own pipeline) | **0.627** | **[0.544, 0.707]** |

`[PENDING: high-accuracy arm, same 165 designs — the controlled comparison]`

---

## 4. Discussion

**What to recommend.** Draw at least three sequences per backbone. Do not use
single-sequence AlphaFold as a triage step; if a cheap screen is needed, reuse
the target's alignment, which costs one MSA search per campaign.

**Limitations — state these plainly and first in the paragraph.**
- Three targets, ten backbones each. The variance direction is consistent but the
  magnitude rests largely on one hard target.
- The effect narrows under ipSAE on the two easier targets (43.4%, 55.1%).
- Everything is prediction. No binding assay was performed on our own designs.
- Six draws per backbone; the best-of-N curve is bootstrapped from those six and
  is biased low at N = 6.

> Reviewers trust papers that find their own weaknesses. Every limitation above
> is one someone would otherwise raise.

---

## 5. Data and code availability

github.com/PRANAVV1107/protein-calibration-kit — all scripts, fold outputs and
analysis code.

---

## Notes on writing this

> Delete this section before submitting; it is scaffolding.

- **One claim per paper.** The headline is the variance decomposition. Everything
  else is support. Resist adding findings that dilute it.
- **Numbers in the abstract.** Most readers read only that.
- **Report what contradicts you.** The ipSAE narrowing, the single hard target,
  the p = 0.222 RBD non-replication. Finding these yourself is what makes the
  rest credible.
- **Cite what you confirm.** For composition collapse and target-MSA screening,
  say "consistent with" and cite — claiming novelty on known practice is the
  fastest route to rejection.
- **Passive voice is fine in methods, active elsewhere.** "We measured" beats
  "it was measured."
- **Figures carry the paper.** Three good ones beat eight mediocre ones.
