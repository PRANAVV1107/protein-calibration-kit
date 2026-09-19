# Summary

Index and one-line verdict for every report in this folder. Written 2026-09-18.

**What this project is.** A de novo binder design pipeline (RFdiffusion →
ProteinMPNN → AlphaFold2) run end to end on one consumer GPU, with every cheap
shortcut in it measured against 365 designs that other groups synthesised and
assayed. **No design of ours has been tested at a bench.**

---

## The reports

| file | what it covers | verdict |
|---|---|---|
| [RESULTS.md](RESULTS.md) | all nine findings, limits, and the error log | the main document |
| [NOISE_RESULT.md](NOISE_RESULT.md) | why the search failed: max-of-K is a tail statistic | **strongest result** |
| [MCTS_PREREGISTRATION.md](MCTS_PREREGISTRATION.md) | success criteria fixed before the search ran | committed before data existed |
| [MCTS_RESULT.md](MCTS_RESULT.md) | the pre-registered search outcome | search lost, cause diagnosed |
| RUN2_RESULT.md | `[PENDING]` — does the mean-value fix work? | running, ~21 h |

Supporting numeric output lives in `../results/*_RESULT.txt`.

---

## The three claims that lead

Confirmed as not-obviously-published by a literature check on 2026-09-17.

**1. Single-sequence screening selects worse than random.** −0.162 ipTM against
random choice. Traced to the *target* failing to fold without an alignment
(pLDDT 27.2 versus 84.3), so the interface is scored against a disordered
chain. Four-condition dose-response. Against measured binding on 164 designs,
cheap-screen interface metrics fall below chance (ipSAE 0.466) while the same
metrics at full accuracy reach 0.618 — paired gain **+0.155, p = 0.001**.

**2. Models and recycles are nearly free to cut.** Once the target alignment is
reused, 1 model / 3 recycles gives ipTM AUC **0.627** (0.620 length-controlled)
at ~63 s/design, against 0.614 (0.605) for 3 models / 8 recycles at ~240 s.
**A quarter of the compute, no loss.** The alignment does the work.

**3. `max`-of-K is a tail statistic and it breaks value estimation.** On a
right-skewed target (RBD, skew +3.38 under cheap scoring), the error in a
backbone's score *grows* with the number of draws — 0.120 at K=2 to 0.214 at
K=8 — while `mean`-of-K falls 0.071 → 0.044. `max` is correct for choosing a
design to synthesise, `mean` for valuing a region to search. Conflating them
explains every anomaly in the failed search.

## Supporting findings

| | |
|---|---|
| sequence choice dominates backbone choice | 59–78% of variance; best-of-6 gains +0.155 pooled, +0.246 on RBD |
| our scoring ≈ the competition's own pipeline | ΔAUC −0.002, p = 0.979 — underpowered at 54 binders, not equivalence |
| no metric ranks affinity | all \|ρ\| < 0.23 across 3.9 decades of K_D, n = 55 — confirms published commentary |
| ipTM predicts physical interface rigidity | ρ = −0.486, p = 0.001, 37 complexes, 2 ns MD — the only orthogonal check |
| skew tracks target difficulty | PD-L1 −1.98, IL-7Rα −1.83, RBD +0.54; n = 3, illustrative only |
| 7 of 15 candidates cross-react with homologs | decoy panel; abundance tier clean |
| preference tuning does not transfer | DPO 0.844 in-distribution, 0.490 held-out; sampling +0.018, p = 0.372 |
| monomer foldability | clean null, ρ = +0.081, p = 0.27 |

---

## What is running now

**Search run 2** — identical to the pre-registered run except nodes are valued
by `mean` rather than `max`. Tests whether the correction implied by
NOISE_RESULT.md actually helps, which the paper currently asserts without
demonstrating.

    blind sampling, 60 folds      0.870   <- the bar
    run 1, max-valued, 1012 folds 0.860
    run 2, mean-valued            [PENDING]

Three outcomes, all reportable: clearing 0.870 means diagnosis and fix both
hold; landing near 0.860 means the statistic was not the whole story; doing
worse means the noise measurement does not transfer to live search and the
"free fix" claim comes out of the paper.

---

## Standing limits on everything here

- No wet-lab validation of our own designs
- Three targets; the search result is n = 1 (one root, one seed, one config)
- 54 binders bounds every null — "no detectable relationship" is not "no
  relationship"
- MD is 2 ns: long enough for rigidity, too short for dissociation
- One consumer NVIDIA GPU constrained which configurations could be tested

## Errors found and corrected

Kept because they shaped the results and because the corrections are part of
the method. Full table in [RESULTS.md](RESULTS.md).

A stub emitting `random.choice()` sequences passed three validation stages, one
reaching ipTM 0.615. A PDB was chosen by name similarity (1ILR is not IL-7Rα).
A pandas truthiness bug inverted all 325 non-binder labels. Cetuximab — an
approved therapeutic — was scored as a de novo design. A length confound that
looked decisive at n=52 (AUC 0.812) vanished at full N (0.569). ColabFold folds
shortest-first, so every interim AUC was computed on a biased short tail. And a
speculative JAX memory change stalled every fold and was misdiagnosed three
times before a four-minute A/B settled it.
