# Results

Everything measured in this project, with the numbers that support it and the
limits on each claim. All computation on one RTX 3060 (6 GB) under WSL2.

The only ground truth here is the Adaptyv EGFR competition set — 365 designs
that were synthesised and assayed by other groups, re-scored by us. **No design
of ours has been tested at a bench.** Every other number is a prediction
checked against another prediction, and is labelled as such.

---

## 1. Sequence choice dominates backbone choice

Within-backbone (sequence) share of total interface-score variance, 180
high-accuracy predictions over 30 backbones, 6 sequence draws each:

| target | ipTM | ipSAE_min | pDockQ2_min |
|---|---|---|---|
| IL-7Rα | 59.0% | 43.4% | 53.4% |
| PD-L1 | 68.9% | 55.1% | 66.4% |
| SARS-CoV-2 RBD | 78.3% | 79.8% | 78.9% |

Drawing the best of six sequences rather than one raises mean ipTM by **0.155
pooled** and **0.246 on RBD**, with no additional backbone generation.

> Narrows under ipSAE on the two easier targets (43.4%, 55.1%). Three targets,
> ten backbones each; the magnitude rests largely on one hard target.

## 2. Cheap screening is anti-predictive, and the mechanism is known

On 164 designs with measured binding outcomes, folded under both settings:

| metric | cheap screen | high accuracy |
|---|---|---|
| ipSAE_min | **0.466** | 0.618 |
| pDockQ2_min | **0.461** | 0.623 |
| ipTM | 0.545 | 0.614 |

Paired on identical designs: high accuracy beats cheap by **+0.153 on ipSAE
(p = 0.001)**. The cheap screen's interface metrics sit *below chance*.

The cause, from a four-condition dose-response: without an alignment the
**target** fails to fold (pLDDT 27.2 against 84.3), so the interface is scored
against a disordered chain. Reusing a precomputed target alignment recovers
about 47% of the gap.

## 3. Our scoring is indistinguishable from the competition's own pipeline

Same 164 designs, paired bootstrap:

    high accuracy - competition's own ipTM   -0.016  [-0.121, +0.090]  p = 0.791

> "Not distinguishable" at 54 binders reflects limited power and is **not**
> evidence of equivalence. Say it that way or a reviewer will.

## 4. No metric ranks affinity among designs that bind

All 55 binders carry a measured K_D spanning 1.2 nM to 10 µM — 3.9 decades.
Spearman ρ against log10(K_D):

| metric | ρ | p |
|---|---|---|
| ipTM | +0.150 | 0.28 |
| pLDDT | +0.139 | 0.32 |
| PAE_interaction | −0.164 | 0.24 |
| ESM pseudo-likelihood | +0.226 | 0.10 |
| binder length | −0.155 | 0.25 |

At n = 55, significance needs |ρ| > 0.27. **Nothing reaches it.**

These metrics separate binders from non-binders and carry no information about
how tightly a binder binds. An objective built from them can optimise hit rate
but not affinity — a search or fine-tuning procedure driven by them is free to
converge on many weak binders while the objective registers success.

## 5. ipTM predicts physical interface stability

The one *orthogonal* check — every other result above is AlphaFold measured
against AlphaFold. 37 complexes, three targets, 2 ns molecular dynamics
(OpenMM, Amber14, GBn2 implicit solvent), panel spanning ipTM 0.13–0.92 by
design so the predictor is not truncated:

| stability measure | ipTM ρ | p |
|---|---|---|
| **binder RMSD** | **−0.486** | **0.001** |
| contacts retained | −0.112 | 0.504 |
| centre-of-mass drift | −0.113 | 0.502 |

Higher confidence means less binder movement relative to the target. Strongest
within RBD (ρ = −0.637).

> 0 of 37 complexes dissociated, so the other two measures have almost no
> variance — they are **uninformative at 2 ns**, not null. That follows from
> cutting 50 ns to 2 ns to buy n = 37 over n = 15.

## 6. Nothing cheap predicts design quality

Three independent attempts, three failures, one underlying cause.

**The generator cannot rank its own draws.** ProteinMPNN's sequence
log-likelihood against the fold score it predicts:

    pooled           rho -0.522   (looks excellent)
    within-backbone  rho -0.042   (the task that matters)

The pooled figure is entirely between-target and between-backbone structure.
Reporting it alone would have justified building on a value function that does
nothing.

**Preference tuning does not transfer.** DPO on 215 within-backbone pairs from
two targets, tested on a third held out entirely:

    in-distribution validation   0.844
    held-out target, untuned     0.449
    held-out target, tuned       0.490   (+0.041, McNemar p = 0.248)

**Guided search does not beat random sampling.** Pre-registered, 1012 folds,
253 backbones, 21 hours:

    guided search   0.860  in 1012 folds
    blind sampling  0.870  in   60 folds

Every improvement came from expanding a *weak* node; expanding the best-known
node produced nothing above 0.27 on three occasions; the depth profile is flat
from depth 2 onward. A backbone's score does not predict its neighbours', so
UCT propagates noise.

**Taken together:** search, learned preference, and heuristic value all fail
for the same reason. Nothing available cheaply predicts which design will fold
well, which is why the expensive evaluation cannot be avoided.

## 7. Specificity is not free

A tiered decoy panel (homologs / abundant serum proteins / unrelated folds)
found **7 of 15** IL-7Rα candidates cross-reacting with cytokine-receptor
relatives. The best on-target scorer hits IL-21R at 0.480 while passing the
unrelated-fold decoys. The abundance tier (IgG Fc, albumin) was clean.

## 8. Monomer foldability is a clean null

0 of 183 designs predicted disordered; ρ = +0.081, p = 0.27 against interface
score. The filter was removed.

---

## What this project has not done

- No wet-lab validation of any design of ours
- Three targets, one hard (RBD) carrying much of the variance result
- The search result is n = 1: one root, one target, one seed, one configuration
- MD is 2 ns — long enough for rigidity, too short for dissociation
- The affinity and transfer nulls are bounded by 54–55 binders; "no detectable
  relationship" is not "no relationship"

## Errors found and corrected during this work

Recorded because they shaped the results and because the corrections are part
of the method.

| error | consequence |
|---|---|
| stub emitting `random.choice()` sequences as designs | 100 fake designs passed three validation stages, one reaching ipTM 0.615 |
| PDB chosen by name similarity (1ILR ≠ IL-7Rα) | wrong target; caught by reading the TITLE record |
| pandas truthiness on the string `'false'` | inverted all 325 non-binder labels |
| Cetuximab scored as a de novo design | approved therapeutic in the validation set; now excluded |
| length confound at n = 52 (AUC 0.812) | vanished at full N (0.569) — sampling artifact |
| ColabFold folds shortest-first | every interim AUC computed on a biased short tail |
| speculative JAX `MEM_FRACTION` change | stalled every fold; misdiagnosed three times before an A/B settled it |

The quarantined bad inputs are kept in `QUARANTINE/` with their provenance
rather than deleted.
