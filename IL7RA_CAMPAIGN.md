# IL-7Rα Binder Campaign — 50-Design Pilot + Sequence-Sampling Correction

**Date:** September 10–12, 2026
**Target:** PDB 3DI3 chain B — human IL-7 receptor α ectodomain, residues 17–209 (193 aa)
**Pipeline:** RFdiffusion → ProteinMPNN → ColabFold (fast screen) → high-accuracy rescore → negative controls → multi-draw sequence refinement
**Status:** Complete. Top candidate `il7ra_binder_6__s2` at **ipTM 0.88 / pTM 0.89 / pLDDT 94.8**, clean composition. Exceeds the RBD campaign's best confirmed result (0.799).

> **Two false starts preceded this run and are documented in §1 — they are the
> most transferable lesson here, more so than the final score.**

---

## 1. False starts (read this first)

### 1.1 Wrong target protein

The campaign initially ran against **PDB 1ILR**, believed to be IL-7Rα. It is not — `1ILR` is
**interleukin-1 receptor antagonist (IL-1Ra)**, an unrelated protein. The error came from an
unverified PDB ID and survived structure download, design, folding, and scoring without
tripping anything, because every downstream stage happily accepts any valid PDB.

Caught by reading the file header:

```
TITLE     CRYSTAL STRUCTURE OF THE INTERLEUKIN-1 RECEPTOR ANTAGONIST
```

The correct target is **3DI3** (McElroy, Dohm & Walsh, *Structure* 2009), chain B, whose header
reads `...GLYCOSYLATED HUMAN INTERLEUKIN-7 RECEPTOR ALPHA ECTODOMAIN`.

**Rule adopted:** grep `HEADER`/`TITLE`/`COMPND` of any newly downloaded PDB and confirm the
molecule name before any compute is spent on it.

### 1.2 Designs that were never designed

The first 100 "IL-7Rα binders" were **random amino-acid strings**. The generating script
contained a silent fallback:

```python
print("[SKIP] Skipping actual RFdiffusion run (not available in this environment)")
...
seq = "".join(random.choice(aa) for _ in range(length))   # <- this is what ran
```

100 random peptides were folded, scored, and one (`pair_026`, ipTM 0.42 fast-screen → 0.61
high-accuracy → 0.615 after 5-model ensemble + Amber relax) advanced through three validation
stages before the stub was noticed. Nothing in the scoring flagged it.

**Rule adopted:** a design stage that cannot run must fail loudly, never substitute placeholder
data. Verify backbone `.pdb` files exist and that ProteinMPNN headers report the expected
`fixed_chains`/`designed_chains` before trusting any downstream score.

### 1.3 Wrong benchmark file folded

`01_fold_benchmark.py` hard-coded `benchmark/adaptyv_benchmark.csv`. Launched for IL-7Rα, it
spent ~2 h re-folding 158 pairs of the *previous* RBD/PD-L1 campaign. Fixed by adding
`--benchmark` and `--output-dir` flags.

---

## 2. Target setup (verified)

| Property | Value |
|---|---|
| PDB | 3DI3 (IL-7 : IL-7Rα complex, glycosylated form) |
| Target chain | B (chain A is the IL-7 ligand, stripped) |
| Residue range | 17–209, contiguous, **no crystal gaps** |
| Length | 193 aa |
| Hotspots | **S31, K77, K138, Y192** — all four confirmed present with the expected residue identity in chain B |
| Binder length | 50–120 aa (contig `B17-209/0 50-120`) |

Hotspots are the IL-7-contacting residues from the 3DI3 paper's H-bond table. PDB numbering
matched the paper's construct numbering directly (no offset), verified residue-by-residue.

Target sequence in the benchmark CSV was checked byte-for-byte against the sequence extracted
from the PDB used by RFdiffusion — exact match.

---

## 3. Fast screen (50 designs)

| Metric | ipTM | pTM |
|---|---|---|
| Mean | 0.091 | 0.309 |
| Std dev | 0.015 | 0.071 |
| Min–Max | 0.06–0.12 | 0.16–0.42 |

Flat and featureless, matching the RBD 50-design pilot (mean 0.090, max 0.123) and PD-L1 pilot
(mean 0.103, max 0.170). No standout design.

### 3.1 The fast screen has no ranking signal

Folding the fast-screen top 10 at high accuracy gave:

| Design | Fast ipTM | High-acc ipTM |
|---|---|---|
| pair_042 | 0.11 | 0.86 |
| pair_043 | 0.11 | 0.85 |
| pair_046 | 0.12 | 0.84 |
| pair_002 | 0.12 | 0.83 |
| pair_028 | 0.12 | 0.82 |
| pair_021 | 0.11 | 0.75 |
| pair_024 | 0.11 | 0.71 |
| pair_013 | 0.12 | 0.49 |
| pair_041 | 0.12 | 0.30 |
| pair_008 | 0.10 | 0.25 |

**Spearman ρ = 0.060 (p = 0.87, n = 10).** Designs sharing fast ipTM 0.12 landed between 0.30
and 0.84 at high accuracy.

**Scope limit:** this is measured within a narrow band (fast 0.10–0.12), and range restriction
deflates correlation. It shows the fast screen cannot rank *within* its top band. It does **not**
establish whether the screen usefully separates that band from the bottom 40 — that remains
**untested**. The experiment that would settle it: high-accuracy fold the 10 *lowest* fast-screen
designs and compare distributions (~45 min of GPU).

---

## 4. The 40%-alanine problem

The top high-accuracy hit, `pair_042` (ipTM 0.858, pTM 0.871, pLDDT 93.0), was **40% alanine**
with a 7-residue poly-Ala run — against ~8% alanine in natural proteins. A separate design in
the same batch reached **59% alanine with a 16-residue run**.

Cause: ProteinMPNN at `--sampling_temp 0.1` (near-greedy) collapsing into low-complexity
sequence. Across all 50 designs, mean Shannon entropy was 3.01 bits (random 20-aa ≈ 4.32) and
mean max-single-residue fraction 27.7%.

### 4.1 Is the score real? — negative controls

| Test | ipTM |
|---|---|
| `pair_042` on-target (IL-7Rα) | **0.858** |
| Solo re-fold, independent of batch | 0.858 (bit-identical — determinism confirmed) |
| **Scrambled** — identical residue composition, order destroyed | **0.240** |
| **Off-target** — same binder vs PD-L1 | **0.320** |
| **Generic Ala-rich helix** (`AAAAEAAAKA…`) vs IL-7Rα | **0.120** |

Per-chain breakdown of the on-target fold: binder pLDDT 91.4 (not a disordered blob),
intra-binder PAE 3.0 Å, interface PAE 5.8 Å, ipSAE 0.71, pDockQ2 0.61.

**Conclusion: the score is not a composition artifact.** A generic alanine-rich helix earns
0.120, and scrambling `pair_042` while preserving its exact residue counts costs 0.62 ipTM. The
sequence *order* carries the signal, and the binder is target-specific (2.7× on- vs off-target).
`pair_042` was a real designed interface expressed in a degenerate sequence — a synthesizability
liability, not a fake hit.

---

## 5. Fixing composition

### 5.1 Temperature alone does not work

ProteinMPNN sweep over all 50 backbones:

| Config | mean Ala% | max Ala% | max run | mean penalty |
|---|---|---|---|---|
| temp 0.1, no bias *(original)* | 20.6% | 59.3% | 16 | 0.083 |
| temp 0.2, no bias | 20.4% | 60.2% | 16 | 0.072 |
| temp 0.3, no bias | 20.8% | 44.4% | 10 | 0.057 |
| temp 0.1, **A bias −1.0** | 9.1% | 15.9% | 6 | 0.058 |
| **temp 0.2, A bias −1.0** | **8.8%** | **15.8%** | **5** | **0.051** |
| temp 0.2, A bias −2.0 | 4.2% | 8.4% | 5 | 0.051 |
| temp 0.2, A −1.0 / G −0.5 | 9.3% | 15.8% | 6 | 0.056 |

Raising temperature 0.1 → 0.2 left mean alanine at ~20% and a 16-residue run intact. The
**alanine logit bias** is what works, landing mean Ala at 8.8% ≈ natural abundance. A −2.0 bias
over-corrects to 4.2%, below natural — alanine is legitimately useful in helices, so that risks
real design quality.

### 5.2 Designs are fragile to sequence identity

Re-drawing sequences on **unchanged** backbones:

| backbone | orig (t0.1) | t0.2 **no bias** | t0.2 **bias −1.0** | mutations orig→nobias |
|---|---|---|---|---|
| il7ra_binder_47 (pair_042) | 0.86 | 0.82 | 0.81 | 18/120 |
| il7ra_binder_48 (pair_043) | 0.85 | **0.58** | 0.29 | 11/95 |
| il7ra_binder_6 (pair_046) | 0.84 | 0.84 | 0.75 | 14/83 |
| il7ra_binder_10 (pair_002) | 0.83 | 0.74 | 0.71 | 9/111 |
| il7ra_binder_34 (pair_028) | 0.82 | 0.81 | 0.79 | 18/94 |
| **mean** | **0.840** | **0.758** | **0.670** | |

The no-bias column is the control that matters: with **no bias at all**, just a different
ProteinMPNN draw, `pair_043` fell 0.85 → 0.58. So ~half the apparent "cost of the bias"
(0.840 → 0.670) was ordinary **sampling variance**, not the bias.

Changing 9–18 residues (~10–15% of positions) on an unchanged backbone swings ipTM by up to
0.27. **The backbone is not the whole story — the specific sequence is load-bearing.**

### 5.3 Multi-draw selection removes the cost entirely

8 draws per backbone (temp 0.2, A bias −1.0) → composition filter → high-accuracy fold the 3
cleanest → select on shaped reward:

| backbone | A) orig, 1 draw | B) 1 biased draw | **C) best of 8** | orig Ala% | C penalty |
|---|---|---|---|---|---|
| il7ra_binder_6 | 0.84 | 0.75 | **0.88** ⬆ | 15.7% | 0.000 |
| il7ra_binder_34 | 0.82 | 0.79 | **0.89** ⬆ | 13.8% | 0.019 |
| il7ra_binder_47 | 0.86 | 0.81 | **0.86** = | **40.0%** | **0.000** |
| il7ra_binder_48 | 0.85 | **0.29** | **0.84** ⬆ | 16.8% | 0.001 |
| il7ra_binder_10 | 0.83 | 0.71 | **0.83** = | 15.3% | 0.000 |
| **mean** | **0.840** | **0.670** | **0.860** | | |

Clean composition at **no cost in score** — mean ipTM slightly exceeds the original. The
`il7ra_binder_48` backbone recovered 0.29 → 0.84, confirming that collapse was sampling
variance. `il7ra_binder_47` holds 0.86 while dropping from 40% alanine to 16.7% max-residue.

Within-backbone spread across just 3 kept draws: `binder_48` 0.62/0.81/0.84, `binder_34`
0.64/0.86/0.89, `binder_10` 0.70/0.80/0.83 — a ~0.2–0.25 swing. Single-draw design was a lottery.

---

## 5b. Metric correction: ipTM → ipSAE

After the sections below were written, a literature check surfaced Dunbrack 2025
("Rēs ipSAE loquuntur: What's wrong with AlphaFold's ipTM score and how to fix
it"). ipTM is computed across whole chains, so unequal chain lengths shift it
spuriously — **this is the published fix for the size-asymmetry artifact recorded
independently in this project's own notes**, where the true ACE2/RBD complex
scored lower than an off-target pair. ipSAE restricts the calculation to residue
pairs below a pAE threshold. Overath et al. 2025 (meta-analysis, 3,766
experimentally characterised binders) found AF3/Boltz-1 with **ipSAE_min > 0.61**
beat AF2-initial-guess on ipAE by ~1.4× average precision.

ipSAE was already present in every ColabFold JSON here and had simply been
ignored. Re-analysis required no new folding.

### Three things changed

**a. The variance result narrows.** Within-backbone (sequence) share of variance:

| target | ipTM | **ipSAE_min** | pDockQ2_min |
|---|---|---|---|
| IL-7Rα | 59.0% | **43.4%** | 53.4% |
| PD-L1 | 68.9% | **55.1%** | 66.4% |
| RBD | 78.3% | **79.8%** | 78.9% |

RBD holds at ~79% under every metric. But under ipSAE, IL-7Rα falls **below half**
(43.4%), meaning backbone choice slightly dominates on that target. The defensible
claim is therefore: *sequence dominates decisively on a hard target, and is
roughly co-equal with backbone on easy targets.* Not the universal claim implied
by the ipTM-only numbers.

**b. pDockQ2 selects better than either ipTM or ipSAE.** Mean high-accuracy ipTM
of the sequence each cheap metric picks (random 0.396, oracle 0.698):

| selecting by | fast screen | hybrid MSA |
|---|---|---|
| ipTM | 0.234 (−0.162) | 0.537 (+0.141) |
| ipSAE_min | 0.292 (−0.104) | 0.473 (+0.077) |
| **pDockQ2_min** | **0.341 (−0.055)** | **0.562 (+0.166)** |

Note this is a *different task* from the literature's recommendation: Overath
ranks ipSAE best for predicting experimental binding success, whereas this
measures within-backbone sequence selection. No contradiction.

**c. The specificity ranking reorders, and off-target signal mostly evaporates.**

| candidate | off ipTM | off ipSAE_min | ipSAE margin |
|---|---|---|---|
| **binder_47__s1** | 0.320 | **0.052** | **0.665** |
| binder_34__s6 | 0.300 | 0.064 | 0.663 |
| binder_34__s5 | 0.560 | 0.192 | 0.478 |
| binder_6__s2 | 0.610 | 0.218 | 0.457 |
| binder_47__s2 | 0.690 | 0.338 | 0.318 |

Off-target ipTM of 0.32 looked like meaningful residual binding; the same pair at
ipSAE_min 0.052 is essentially no interface at all. The apparent promiscuity was
substantially a chain-length artifact — precisely what ipSAE exists to remove.

### Epitope analysis (no GPU required)

Contacts computed from existing predicted structures, against IL-7's own footprint
in 3DI3:

- IL-7's natural epitope on IL-7Rα is **19 residues** (31, 33, 57–59, 77–83, 102,
  104, 138–139, 191–193).
- Designed binders overlap it by **15.9 residues on average**; `binder_47__s3`
  reaches 19/19.
- So the designs occupy the real ligand site and would plausibly **compete with
  IL-7** — the intended therapeutic mechanism.
- But epitope overlap does **not** predict specificity (ρ = −0.36, p = 0.55, n=5);
  `binder_47__s2` covers 17/19 and has the worst margin. Promiscuity appears to
  come from the binder's own surface properties, not from where it binds.

---

## 6. Final result

> **Superseded by §5b.** Under ipSAE_min the top candidate is
> **`il7ra_binder_47__s1`** (ipTM 0.86, pTM 0.89, pLDDT 92.5, **ipSAE_min 0.717**,
> pDockQ2_min 0.487, 120 aa), which also engages all four requested hotspots and
> covers 18/19 of the natural IL-7 epitope. Its off-target ipSAE_min is 0.052.
>
> ```
> MTPACLKLAELVKSKSPEEAGELAGKAIACLASSFINPSNLEVCDPELAEQVAQLTTVEEKIECMKLLADMAAEIGKEAAGKVYIGAAVGASYLEDEGGDEKLIKGLKEIAEKAYEAYKK
> ```
>
> The ipTM-selected candidate below remains accurate as recorded, but ipTM was the
> wrong metric to rank on.

**`il7ra_binder_6__s2`**

| Metric | Value |
|---|---|
| ipTM | **0.88** |
| pTM | 0.89 |
| mean pLDDT | 94.8 |
| Composition penalty | 0.000 |
| Max single-residue fraction | 22.9% |
| Max homopolymer run | 2 |
| Length | 83 aa |

```
MEKKEEIKKLLEETEKRMEKVCKEAGEKGNKELIDKCLEARQEVLGCFENASYYIKKGDLEKAMEEAKKAQEIVDKLEEIVKK
```

For comparison, the RBD campaign's best confirmed candidate (`neighbor_4`) reached ipTM 0.799
after partial-diffusion refinement.

### Not yet done
- Solo verification of `il7ra_binder_6__s2` independent of its batch
- Negative controls (scramble, off-target) on this specific sequence
- 5-model ensemble + Amber relaxation (CPU relax — `--use-gpu-relax` fails on this setup with
  "Minimization failed after 100 attempts", an OpenMM/CUDA gap under WSL2)
- Partial-diffusion backbone refinement around `il7ra_binder_6`

**These are required before calling this a validated hit.** §4.1 exists precisely because a
0.858 number survived three stages of validation before controls were run.

---

## 7. Pipeline changes adopted

`scripts/run_campaign_pipeline.sh` now implements, by default:

1. **ProteinMPNN `--num_seq_per_target 8`** instead of 1 — justified by §5.2
2. **`--bias_AA_jsonl {"A": -1.0}`, `--sampling_temp 0.2`** — justified by §5.1
3. **Composition filter before folding** — design is ~1.5 s, high-accuracy folding ~4.5 min
4. **Ranking by shaped reward** = ipTM − composition penalty

`external_tools/rl/seq_composition.py` implements the penalty (homopolymer runs beyond 4,
single-residue fraction beyond 25%, Shannon entropy below 3.2 bits; capped at 0.50). Calibration:
IL-7, PD-L1, and RBD native sequences all score **0.000**; `pair_042` scores 0.307.

It is wired into the RL refinement reward (`rl_refine_binder_rbd.py`), so the policy cannot gain
reward by collapsing into poly-alanine.

---

## 8. Cost

| Stage | Wall-clock |
|---|---|
| RFdiffusion, 50 backbones | ~3 h 10 m |
| ProteinMPNN, 50 × 1 seq | ~4 m |
| ColabFold fast screen, 50 designs | ~2 h |
| High-accuracy rescore, top 10 | ~1 h 20 m |
| Negative controls, 3 folds | ~15 m |
| ProteinMPNN sweep, 7 configs × 50 | ~6 m |
| Multi-draw refinement, 15 high-acc folds | ~70 m (4.5 min/fold) |

Wasted: ~2 h re-folding the wrong benchmark (§1.3), plus the entire random-sequence run (§1.2).

---

## 8b. Where this sits in the literature

Checked after the fact, which was the wrong order — several of these would have
changed how the campaign was run.

**Probably novel.** No published work decomposes AF2 interface-score variance into
within-backbone (sequence) versus between-backbone components. ProtDBench samples
8 sequences per backbone "to reduce variance" without reporting the split; the
RFdiffusion protocol fixes 8 with no ablation justifying it. The measurement here
— across three targets and three metrics — appears unpublished.

**Already standard, reinvented here.** Bennett et al. 2023 (`dl_binder_design`)
already runs AF2 with an MSA for the target and single sequence for the binder,
precisely because the target needs coevolutionary signal. The "hybrid MSA" built
in this project is that practice, rediscovered. What is *not* documented anywhere
found is the quantitative failure it avoids (target pLDDT 27.2 vs 84.3, binder
unaffected at 70.1 vs 72.3) — useful as a practitioner note, not a discovery.

**Superseded tooling.** AF2-initial-guess remains the reference baseline, but
Overath et al. 2025 report AF3 and Boltz-1 outperforming it, with ipSAE_min the
better filter. BindCraft (Pacesa et al., *Nature* 2025) reports 10–100% success
screening ~10 designs, co-folding the target each iteration.

**Contested, currently.** Specificity screening is not standardised: decoy panels
are not routine, and no published campaign was found using ProteinMPNN's
`--pos_neg_chain_betas` multi-state flags. Two 2026 preprints (RedNet; Odin-Multi)
address binder selectivity directly.

**Known folklore.** ProteinMPNN's low-complexity collapse at low sampling
temperature is discussed in RFdiffusion issue #102, with a negative alanine bias
as the standard remedy — the same fix arrived at here empirically.

## 9. Transferable lessons

1. **Verify the target's identity from the PDB header** before spending compute. A wrong PDB ID
   propagates silently through every stage.
2. **A skipped pipeline stage must fail loudly.** Placeholder data that looks plausible will
   survive multiple validation stages.
3. **Draw multiple sequences per backbone.** ipTM varies by ~0.2–0.25 across draws on an
   unchanged backbone; one draw is a lottery ticket.
4. **Check sequence composition, not just score.** AF2 scores a 40%-alanine design at 0.858 with
   pLDDT 93 — the score is real but the sequence is not synthesis-grade.
5. **Run negative controls before believing a hit**: scramble (same composition, destroyed order),
   off-target, and a generic-fold control. All three were needed here to distinguish "real
   interface, bad sequence" from "composition artifact".
6. **Isolate one variable at a time.** The no-bias control (§5.2) is what revealed that half the
   apparent cost of the alanine bias was sampling variance.
7. **Check which metric you are optimising before optimising it.** ipTM was used throughout, while
   ipSAE and pDockQ2 sat unread in the same JSON files. Switching metric reordered the specificity
   ranking, changed the top candidate, and cut the headline variance figure on one target from 59%
   to 43%. Zero new compute.
8. **Read the literature first, not last.** The target-MSA fix has been standard practice since
   2023, and the ipTM size-asymmetry artifact this project rediscovered by experiment has a
   published correction (ipSAE). Both were findable in an hour.
9. **Decompose a failing stage before discarding it.** The fast screen changed three variables at
   once (no MSA, 1 model, 3 recycles). Only one mattered: holding the other two fixed while
   restoring the target MSA moved selection from worse-than-random (−0.162) to better-than-random
   (+0.166 with pDockQ2).
