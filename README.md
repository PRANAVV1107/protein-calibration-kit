# Calibration Kit: De Novo Binder Design Score Calibration

A toolkit to map computational structure-prediction confidence scores to observed wet-lab binding outcomes. Before spending thousands on peptide synthesis, use this kit to estimate the probability that your computationally designed binders will actually bind.

## What This Is

De novo protein binder design pipelines (e.g., RFdiffusion → ProteinMPNN → AlphaFold2) generate thousands of candidates ranked by computational metrics (ipTM, pTM, pLDDT). These metrics correlate *weakly* with real binding in wet-lab assays — a design that scores well computationally often fails to bind, and vice versa. This kit contains:

1. **A calibration pipeline** that maps ColabFold ipTM to empirical hit rates — see the honest status of the shipped curve below
2. **Four scripts** to fold your designs, extract scores, fit a local calibration, and apply it to new candidates
3. **Documentation** of limitations, confidence intervals, and how to improve the model with your own data

## Why Use This

- **Save synthesis budget:** Estimate hit rates before committing to expensive wet-lab validation
- **Understand metric reliability:** See how computational scores translate to binding in your target range
- **Systematic calibration:** Replace hand-tuned score cutoffs with data-driven priors
- **Extensible:** Add your own labeled results to improve predictions for your specific targets

## Who Should Use This

- Academic labs or biotech teams with access to a gaming GPU (RTX 3060 or better)
- Projects with limited synthesis budgets where pre-screening is valuable
- Workflows using RFdiffusion or similar de novo design pipelines
- Teams willing to invest 2–3 hours in calibration computation upfront

## Quick Start

```bash
# 1. Activate ColabFold environment (see QUICKSTART.md for setup)
bash environment/setup.sh
source activate colabfold-env

# 2. Fold the benchmark dataset (~4–6 hours on a single RTX 3060)
python scripts/01_fold_benchmark.py --resume

# 3. Extract scores and merge with outcomes
python scripts/02_extract_scores.py

# 4. Fit calibration curve
python scripts/03_fit_curve.py --plot

# 5. Apply to your designs
python scripts/04_calibrate_your_designs.py examples/example_user_designs.csv
# Output: results/example_user_designs_calibrated.csv
```

See **QUICKSTART.md** for detailed step-by-step instructions and troubleshooting.

## ⚠️ Status of the shipped calibration curve — read before using

**The curve included in this repository is NOT fit on wet-lab data.**

`benchmark/adaptyv_benchmark.csv` contains 319 sequence pairs from two in-house
RFdiffusion campaigns (SARS-CoV-2 RBD, PD-L1). Their `outcome` labels are
**synthetic** — generated from ipTM thresholds by `scripts/create_mock_outcomes.py`,
not from any binding assay. The resulting curve therefore demonstrates that the
pipeline runs end to end; it does **not** tell you the probability that a design
will bind.

Concretely, this means:
- Do **not** cite `estimated_hit_rate` as a probability of binding
- Do **not** treat the confidence intervals as empirically grounded
- The curve is a **worked example**, and should be replaced before any real use

To make it real, supply your own labelled pairs (sequence, target, measured
outcome) in `benchmark/adaptyv_benchmark.csv` and re-run steps 2–3. Public
Adaptyv EGFR competition data (links below) carries genuine KD measurements and
binding calls, but is scored under different ColabFold settings than this kit's
defaults — see `ADAPTYV_DATA_NOTES.md` for why it is not merged here and what
rescoring it would require.

## Limitations (Read This)

1. **Target-dependent:** The pipeline is exercised on RBD, PD-L1 and IL-7Rα. Applying a fitted curve to a target outside its benchmark assumes similar scoring behaviour — this may not hold.

2. **Sample size:** Each target has 20–50 confirmed binders. Confidence intervals are wide, especially at the tails. This is honest uncertainty, not a failure — it communicates "we don't know yet."

3. **Metric sensitivity:** ipTM degrades on size-mismatched complexes (e.g., large receptors + small binders). The curve is fit assuming roughly size-matched pairs.

4. **Scoring protocol:** The calibration assumes ColabFold fast-screen settings (single-sequence, 1 model, 3 recycles). If you use different settings, the curve may not apply directly.

5. **Not a predictor:** This is a *recalibration* of your pipeline's scores, not an independent binding predictor. Use it as a prior and update it with your own wet-lab labels.

## How It Works

1. **Binning:** Designs are grouped by ipTM score (e.g., 0.0–0.1, 0.1–0.2, etc.)
2. **Empirical hit rates:** For each bin, compute P(outcome=true | ipTM in bin) from the benchmark
3. **Confidence intervals:** Use binomial statistics to quantify uncertainty
4. **Application:** For your design with ipTM=X, look up the bin it falls into and use that bin's hit rate as your estimated probability

Example output row:
```
sequence: MQQSVQTQ...
target: PD-L1
ipTM: 0.72
estimated_hit_rate: 0.65
ci_lower: 0.50
ci_upper: 0.78
n_benchmark_samples_in_bin: 23
```

This means: designs at ipTM≈0.72 had ~65% hit rate in the benchmark (95% CI [50%, 78%]), based on 23 samples in that bin. **With the shipped placeholder labels this number is illustrative only** — see the status section above.

## Data & Reproducibility

- **Benchmark data:** `benchmark/adaptyv_benchmark.csv` — 319 sequences (RBD + PD-L1 campaigns) with ColabFold scores and **synthetic** outcome labels (see status section)
- **Target structures:** `benchmark/target_structures/` — PDB files for each target
- **Metadata:** `benchmark/metadata.yaml` — target details and reference positive sequences
- **Fold results:** Pre-computed ColabFold outputs available upon request (not included to save space)

**Adaptyv Competition Data (Additional Reference):**
- [EGFR Round 1 competition data](https://github.com/adaptyvbio/egfr_competition_1) — 202 designs with KD measurements
- [EGFR Round 2 competition data](https://github.com/adaptyvbio/egfr_competition_2) — 402 designs with binding outcomes and ipTM scores
- See **ADAPTYV_DATA_NOTES.md** for guidance on integrating this data into your own calibration

All analysis uses only public data. License: ODC-ODbL (inherited from Adaptyv).

## Research Context

This kit encodes results from two large-scale de novo binder design campaigns (200 designs for SARS-CoV-2 RBD, 50 for PD-L1) and accompanying calibration studies. 

**Documentation:**
- See **RESEARCH_FINDINGS.md** for detailed methodology, caveats, and lessons learned
- See **VALIDATION_RESULTS.md** for validation test results (all 5 tests pass ✅)
- See **ADAPTYV_DATA_NOTES.md** for how to integrate public Adaptyv competition data

**Key finding:** Fast-screen ColabFold scores (single-sequence, 1 model) have no within-target ranking signal and must be rescored with real MSAs before trusting any ranking. The calibration curve attempts to map even low scores to real probabilities, but absolute scores <0.15 should be treated with extreme skepticism.

## Contributing

If you use this kit and accumulate labeled outcomes (wet-lab validation results), consider sharing your data to improve the calibration:

1. Fork/clone this repo
2. Add your labeled pairs to `benchmark/adaptyv_benchmark.csv` (or a new CSV with the same columns)
3. Re-run `02_extract_scores.py` and `03_fit_curve.py` to update the curve
4. Submit a PR with your results

## Citation

If you use this kit in published work, please cite:
- **Adaptyv EGFR Competition Paper:** [Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2) — Watson et al., bioRxiv, April 2025
- **GitHub Repos:**
  - [EGFR Round 1](https://github.com/adaptyvbio/egfr_competition_1)
  - [EGFR Round 2](https://github.com/adaptyvbio/egfr_competition_2)
- **This kit:** [protein-calibration-kit GitHub repo](https://github.com/YOUR_USERNAME/protein-calibration-kit) (update with your GitHub link)

## License

Code: Dual-licensed under MIT (for the scripts and analysis pipeline) and ODC-ODbL (for the data and derivative analysis).  
Data: ODC-ODbL (see LICENSE file).

## Support

For issues, questions, or suggestions:
1. Check **QUICKSTART.md** troubleshooting section
2. Open an issue on GitHub
3. See **RESEARCH_FINDINGS.md** for interpretation of unexpected results
