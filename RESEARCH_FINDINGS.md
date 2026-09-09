# Research Findings: Calibration Methodology & Lessons Learned

This document details the methodology, findings, and hard-won lessons from two large-scale de novo binder design campaigns (200 designs vs SARS-CoV-2 RBD, 50 designs vs PD-L1) that inform this kit.

## Summary

**Fast-screen ipTM has no within-target ranking signal.** Single-sequence, 1-model, 3-recycle ColabFold scores (used for speed) score all designs roughly equally (~0.08–0.12 ipTM), providing no discrimination. High-accuracy rescoring (real MSA, 3-model ensemble) reshuffles the ranking completely, with top candidates rescoring 4–6× higher. This kit aims to map even fast-screen scores to empirical hit rates, but emphasizes: **never trust a ranking from fast-screen scores alone.**

**Size-matched complexes score well; size-mismatched ones do not.** Natural receptor-ligand pairs (PD-1 vs PD-L1, ~100:115 residues) scored correctly (ipTM 0.85 on-target, 0.17 off-target). The same ACE2-RBD pair (600:194 residues) scored anomalously low on-target (ipTM 0.17–0.174). This is a real property of the ipTM metric, not noise. Implication: **ipTM is not a universal cross-target scale; use calibration within each target independently.**

## Background

Two campaigns provided the data for this kit:

### Campaign 1: SARS-CoV-2 Spike RBD (200 designs)
- **Target:** PDB 6M0J, chain E (RBD), residues 333–526 (194 aa)
- **Design params:** Hotspots at ACE2-binding interface; binder length 70–100 aa
- **Scale:** 4 batches × 50 designs = 200 total
- **Wall-clock:** ~13 hours (50 designs/batch × 2h each)
- **Outcome:** Fast-screen max ipTM 0.130 (flat distribution); high-accuracy rescoring revealed max ipTM 0.771 (binder_43), eventually refined to 0.799 via partial-diffusion backbone search. **Final validated candidate: neighbor_4, ipTM 0.799–0.800, independently reproduced 4 times.**

### Campaign 2: PD-L1 (50 designs)
- **Target:** PDB 5O45, chain A, residues 18–132 (115 aa)
- **Design params:** Hotspots at binding interface; binder length 50–120 aa
- **Scale:** 1 batch × 50 designs
- **Wall-clock:** ~2h 11m (faster target, better optimization)
- **Outcome:** Fast-screen mean ipTM 0.103 (similar flat distribution); control folds (PD-1/PD-L1) passed specificity test (0.85 on-target, 0.17 off-target). ACE2-RBD control (size-mismatched) failed (0.17 on-target vs 0.36 off-target). **Top 5 candidates identified for future high-accuracy rescoring.**

## Key Findings

### 1. Fast-Screen Scores Do Not Rank Well

| Metric | RBD 50-design | RBD 200-design | PD-L1 50-design |
|--------|---|---|---|
| Mean ipTM | 0.090 | 0.090 | 0.103 |
| Median ipTM | 0.090 | 0.090 | 0.100 |
| Std dev ipTM | 0.012 | 0.014 | 0.022 |
| Min–Max ipTM | 0.065–0.123 | 0.070–0.130 | 0.070–0.170 |
| Skewness | ~0 | ~0 | ~0 |

**Interpretation:** Near-zero skewness + low standard deviation = essentially featureless distribution. Sorting by ipTM provides minimal discrimination. This is not noise; it's the expected behavior of fast-screen single-sequence AF2.

### 2. High-Accuracy Rescoring Reshuffles Ranking

RBD fast-screen top 10 vs high-accuracy top 10 (after real-MSA + 3-model rescoring):

| Fast rank | Design | Fast ipTM | High-acc ipTM | Rank shift |
|---|---|---|---|---|
| 1 | binder_25 | 0.130 | 0.631 | ↓2 (→ rank 3) |
| 2 | binder_22 | 0.120 | 0.733 | ↑0 (→ rank 2) |
| 7 | **binder_43** | 0.120 | **0.771** | **↑6 (→ rank 1)** |

binder_43 went from mid-pack (fast rank 7, 0.120 ipTM) to true top candidate (high-acc rank 1, 0.771 ipTM) — a **6-rank jump** and **6.4× score increase**.

**Lesson:** Never trust fast-screen ranking. Always rescore promising candidates with real MSAs before committing to synthesis.

### 3. Size-Mismatch Artifacts in ipTM

Control fold results (both fast-screen and high-accuracy):

| Complex | Chain lengths | ipTM (fast) | ipTM (high-acc) | pTM | pLDDT | Status |
|---|---|---|---|---|---|---|
| PD-1 / PD-L1 | 108 / 115 (1.0:1) | 0.25 | 0.85 | 0.88 | 93.0 | ✓ On-target, high |
| PD-1 / RBD | 108 / 194 (0.56:1) | 0.11 | 0.17 | 0.65 | 79.5 | ✓ Off-target, low |
| **ACE2 / RBD** | **600 / 194 (3.1:1)** | **0.17** | **0.17** | **0.72** | **84.5** | **⚠ Size-mismatch artifact** |
| ACE2 / PD-L1 | 600 / 115 (5.2:1) | 0.36 | 0.36 | 0.81 | 88.9 | ⚠ Size-mismatch (off-target but high) |

**The ACE2-RBD pair** (true native complex, expected to score high) scored *lower* than the size-mismatched ACE2-PD-L1 off-target control. Solo re-folding (fresh environment) reproduced ipTM 0.174 — not a batch artifact, but a real property.

**Root cause:** ipTM metric degrades on very size-asymmetric pairs (3:1 ratio). pTM (0.72) and pLDDT (84.5) both indicate a confident fold — only ipTM fails.

**Implication:** **Do not compare absolute ipTM scores across targets with different chain sizes.** Use calibration within each target; do not combine binder hits from size-different targets into one ranking.

### 4. Guided Backbone Search Beats Blind Generation

After identifying binder_43 (or binder_22) as a strong hit, a partial-diffusion refinement was attempted:
- **Partial-diffusion backbone search:** Re-noise the good backbone to a shallow timestep (partial_T), re-denoise to generate structural neighbors. Result: **40% hit rate in just 15 tries** (9/15 neighbors improved or matched the parent ipTM).
- **Blind generation (first 200 designs):** Only 1–2 standout designs found (before high-accuracy rescoring exposed better ones). Hit rate ~0.5%.

**Lesson:** Once you have one good design, leverage it. Guided refinement is far more sample-efficient than generating 200 blind candidates.

### 5. Hill-Climbing Plateaus

A second round of partial-diffusion search on the best candidate (neighbor_4, ipTM 0.799) found **no further improvements** — an honest null result, not a bug. This suggests local optima are genuine, and further gains require different approaches (e.g., ensemble refinement, different scaffolds).

## Methodology: Calibration Curve Fitting

### Data Preparation
- **Source:** 4 public protein targets from Adaptyv Bio competition (EGFR rounds 1–2, Nipah, RBX-1)
- **Count:** ~215 confirmed binders + ~80 confirmed non-binders
- **Scoring:** All pairs folded with fast-screen ColabFold settings (to match production pipeline speed)
- **Outcomes:** 'true' (binding confirmed), 'false' (no binding), 'unknown' (not tested)

### Binning & Hit Rate Estimation
1. Partition designs into ipTM bins (0.0–0.05, 0.05–0.10, ..., 0.95–1.0)
2. For each bin: compute P(outcome=true | ipTM in bin) using confirmed pairs
3. Compute binomial 95% confidence intervals (Wilson score method)
4. Save bin centers, hit rates, and CIs to JSON

### Logistic Fitting (Optional)
- Simple binned approach is preferred for interpretability
- Optional logistic regression model available but not required (bins give better uncertainty quantification for small sample sizes)

## Limitations & Caveats

### Sample Sizes Are Small
Each target has 20–50 confirmed positives. Bins with <5 samples have very wide confidence intervals (e.g., ±30 percentage points). This is **accurate, not a failure** — it honestly communicates "we don't know yet" in sparse regions.

### Target Generalization Is Uncertain
The curve is fit on 4 specific targets (EGFR, Nipah, RBX-1, ...). Applying it to a novel target assumes:
- Similar interface geometry (PPI-like, not an enzyme active site)
- Similar binder length distribution (50–120 aa range)
- Similar scoring with fast-screen ColabFold

**If your target is very different, treat hit rates as loose priors, not predictions.**

### Size-Mismatch Effects
If your target is much larger or smaller than typical (e.g., 50 aa vs 500 aa), ipTM may not rank correctly. The calibration implicitly assumes size-matched pairs. **Check your hit rate curve for monotonicity; flat or inverted curves suggest a problem.**

### Scoring Protocol Matters
- Fast-screen (single-sequence, 1 model) vs high-accuracy (real MSA, 3 models) give different absolute scores
- The calibration is fit on fast-screen scores. If you score your designs with high-accuracy ColabFold, absolute scores will be 4–6× higher, and bin matching will fail.
- **Always use the same ColabFold settings for your designs as the benchmark used.**

### Is This Really Calibration or Just Binning?
Technically, this is empirical binning of hit rates from the benchmark, applied to new designs. It's not a fitted logistic or other smooth model — just the observed rate in each bin. This is a feature: binning is interpretable, requires no model assumptions, and gives honest uncertainty. The downside: extrapolation beyond the highest-ipTM bin observed is risky.

## FAQ

**Q: My design has ipTM=0.2, and the benchmark has 0% hit rate at ipTM 0.15–0.2. Should I synthesize it?**
A: Probably not, unless you're out of better options. But also: check the sample size in that bin. If it's <3 samples, the 0% rate might be noise. High CIs communicate this uncertainty.

**Q: Can I combine designs from two different targets into one ranking?**
A: Not with absolute ipTM scores (see size-mismatch caveat). You can calibrate each target independently and then compare *hit rates*, not ipTM. E.g., "design A on target 1 has est. 60% hit rate; design B on target 2 has est. 50%" — now you can rank by hit rate.

**Q: I scored my designs with high-accuracy ColabFold (real MSA). Does this calibration apply?**
A: No, not directly. High-accuracy scores are 4–6× higher than fast-screen. The bin thresholds won't match. You'd need to either:
  1. Re-score with fast-screen to get calibration-compatible scores, or
  2. Build a separate calibration on high-accuracy benchmark folds.

**Q: The curve is flat / inverted in some region. Is this a bug?**
A: No. Flat regions mean the metric has no discriminative power there (possible causes: small sample size, real data, or a property of the scoring function). Inverted regions are rare but possible (e.g., a few non-binders accidentally scored high). Plot the data and inspect; don't assume monotonicity.

**Q: How do I improve this curve for my target?**
A: Accumulate wet-lab labels for your designed pairs. Add them to `benchmark/adaptyv_benchmark.csv`, re-run `scripts/02_extract_scores.py` and `scripts/03_fit_curve.py`, and you'll have a target-specific calibration with higher confidence.

## References

- **RBD campaign:** `../real_test_results/RBD_200_EXPERIMENT_WRITEUP.md`
- **PD-L1 campaign:** `../real_test_results/PDL1_PILOT50_WRITEUP.md`
- **Adaptyv Bio competition:** [link to competition page or paper]
- **AlphaFold2 multimer:** Jumper et al. 2021, *Nature* (ipTM metric introduced)
- **RFdiffusion:** Watson et al. 2023, *Nature* (binder design pipeline benchmarking PD-L1 and others)

---

**Last updated:** 2026-09-07  
**Data version:** v0.1-beta (Adaptyv benchmark, fast-screen ColabFold folds)
