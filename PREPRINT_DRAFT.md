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
Against 164 experimentally characterised designs, cheap-screen interface metrics
score below chance (ipSAE AUC 0.466) while the same metrics at full accuracy
reach 0.618, a paired gain of +0.153 (p = 0.001). Finally, among the 55 designs
that do bind, no metric we examined ranks affinity across four decades of K_D
(all |ρ| < 0.23, n = 55) — so these scores support triage but not optimisation
for binding strength.

> The abstract now carries three findings, which is one more than ideal. If it
> runs long, the affinity null is the one to compress, not cut: it is the
> result a reader is least likely to have assumed.

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

### 3.5 Validation against measured binding

**[FIGURE 4]** paired AUC, cheap against high accuracy, same designs.

The Adaptyv EGFR competition released 365 designs that were synthesised and
assayed. We re-scored a balanced subset of 164 (54 binders, 110 non-binders)
under both settings, so the comparison isolates the accuracy setting rather
than the design set.

| metric | AUC | 95% CI |
|---|---|---|
| chain length (confound baseline) | 0.569 | [0.466, 0.667] |
| cheap ipTM | 0.545 | [0.451, 0.638] |
| cheap ipSAE_min | 0.466 | [0.371, 0.564] |
| cheap pDockQ2_min | 0.461 | [0.367, 0.558] |
| high-accuracy ipTM | 0.614 | [0.528, 0.697] |
| high-accuracy ipSAE_min | 0.618 | [0.532, 0.701] |
| high-accuracy pDockQ2_min | 0.623 | [0.536, 0.706] |
| competition's own reported ipTM | 0.630 | [0.534, 0.720] |

Paired on identical designs and identical bootstrap resamples:

| comparison | ΔAUC | 95% CI | p |
|---|---|---|---|
| high accuracy − cheap (ipSAE) | **+0.153** | [+0.057, +0.249] | **0.001** |
| high accuracy − cheap (ipTM) | +0.069 | [−0.009, +0.150] | 0.080 |
| high accuracy − competition's own | −0.016 | [−0.121, +0.090] | 0.791 |

Two results. First, the cheap screen's interface metrics sit **below chance**
(0.466, 0.461) while the same metrics at full accuracy reach 0.618 and 0.623;
the paired gain is significant. Second, our high-accuracy scoring is
statistically indistinguishable from the competition's own production pipeline
on these designs.

> The second claim needs its caveat stated in the same breath: "not
> distinguishable" at 54 binders reflects limited power and is not evidence of
> equivalence. Write it that way or a reviewer will write it for you.

Because binders in this set are longer than non-binders (125 aa against 101),
every metric is also reported stratified within length tertiles and after
regressing length out. High-accuracy ipTM holds (0.615 stratified, 0.605
residualised) and the cheap interface metrics fall further (0.477, 0.466).

> One positive control, Cetuximab_scFv, was excluded. It is an approved
> anti-EGFR therapeutic spiked into the competition as a benchmark, not a de
> novo design, and this study evaluates metrics for de novo designs. Say so
> explicitly; silent exclusions are how papers lose trust.

### 3.6 No metric ranks affinity among binders

**[FIGURE 5]** metric against log K_D, 55 binders.

All 55 binders carry a measured K_D spanning 1.2 nM to 10 µM, nearly four
decades. Among those binders, Spearman ρ against log10(K_D):

| metric | ρ | p |
|---|---|---|
| competition ipTM | +0.150 | 0.28 |
| competition pLDDT | +0.139 | 0.32 |
| PAE_interaction | −0.164 | 0.24 |
| ESM pseudo-likelihood | +0.226 | 0.10 |
| binder length | −0.155 | 0.25 |

At n=55, significance requires |ρ| > 0.27. None of these reach it. Every metric
examined separates binders from non-binders to some degree and none ranks how
tightly a binder binds.

> This is the most consequential negative result here and it belongs in the
> abstract. The practical implication is concrete: an objective built from
> these metrics can optimise hit rate but not affinity, so a search or
> fine-tuning procedure driven by them is free to converge on many weak
> binders while the objective registers success.

---

## 4. Discussion

**What to recommend.** Draw at least three sequences per backbone. Do not use
single-sequence AlphaFold as a triage step; if a cheap screen is needed, reuse
the target's alignment, which costs one MSA search per campaign. Do not treat
any of these metrics as an affinity estimate.

**Limitations — state these plainly and first in the paragraph.**
- Three targets, ten backbones each. The variance direction is consistent but the
  magnitude rests largely on one hard target.
- The effect narrows under ipSAE on the two easier targets (43.4%, 55.1%).
- Everything is prediction. No binding assay was performed on our own designs;
  the only measured outcomes here are other groups' designs, re-scored.
- Six draws per backbone; the best-of-N curve is bootstrapped from those six and
  is biased low at N = 6.
- The validation set is 54 binders. Null results in §3.6 are bounded by that
  power, and "no detectable relationship" is not "no relationship".
- All computation ran on one consumer GPU (RTX 3060, 6 GB). That constrains
  which configurations could be tested, not only how many.

**Future work — search at generation time.** The variance decomposition says
sequence choice dominates, and best-of-N exploits that only greedily. Test-time
search over the generation trajectory is the principled generalisation, and
recent work reports it succeeds specifically on hard targets where single-pass
generation fails (Proteína-Complexa, arXiv:2603.27950). Our RBD target is that
regime: single-draw mean ipTM 0.400 against 0.766 and 0.805 for the other two,
with the largest best-of-N gain (+0.246). An earlier partial-diffusion
experiment here is a crude instance of the same idea — re-noising a known-good
backbone gave a 40% hit rate from 15 attempts against near-zero from 200 blind
designs.

> Two things make this a credible proposal rather than a wish. §3.5 measures
> why the obvious cheap reward fails, which is the failure mode such a search
> would hit first. §3.6 bounds what any such objective can optimise. State the
> compute honestly: at a measured ~5 min per node evaluation, a 500-node search
> is roughly 42 GPU-hours on this hardware, and whether that is reducible is
> itself an open question worth answering.

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
