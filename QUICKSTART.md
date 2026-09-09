# Quick Start: Running the Calibration Kit

This guide walks you through folding the benchmark, fitting the curve, and applying it to your designs. Estimated time: **4–8 hours wall-clock** (mostly unattended folding).

## Prerequisites

- **GPU:** NVIDIA GPU with ≥6 GB VRAM (tested on RTX 3060 Laptop, RTX 3090)
- **Storage:** ~50 GB free (for ColabFold outputs; can be deleted after analysis)
- **Environment:** Linux/Mac or WSL2 on Windows; conda installed
- **Internet:** Access to ColabFold's public MSA server (MMseqs2, rate-limited)

## Step 0: Clone and Set Up Environment

```bash
git clone <this-repo>
cd calibration-kit

# Activate ColabFold environment
bash environment/setup.sh
source activate colabfold-env

# Verify ColabFold and pandas are available
python -c "import colabfold; import pandas; print('✓ Ready')"
```

**Common issues:**
- `ModuleNotFoundError: colabfold`: Run `bash environment/setup.sh` again or manually `pip install colabfold[local]`
- `CUDA out of memory`: Make sure no other GPU processes are running (`nvidia-smi`); close Chrome, etc.

## Step 1: Fold the Benchmark (~4–6 hours)

Fold all 120–150 benchmark sequence pairs against their targets using ColabFold.

```bash
python scripts/01_fold_benchmark.py --resume
```

**What it does:**
- Reads `benchmark/adaptyv_benchmark.csv` (sequence pairs + outcomes)
- Folds each pair using ColabFold (fast-screen: single-sequence, 1 model, 3 recycles)
- Saves outputs to `results/benchmark_folds_raw/pair_NNN_scores_rank_001_*.json`
- Checkpoints every 10 pairs; use `--resume` to skip completed pairs if interrupted

**Expected output:**
```
Loaded 143 benchmark pairs from benchmark/adaptyv_benchmark.csv
[1/143] Folding pair_000...
  ✓ Success
[2/143] Folding pair_001...
...
✓ Complete: folded 143 pairs
Results in: results/benchmark_folds_raw
```

**Typical runtime:**
- ~2 min per pair (RTX 3060, low-VRAM mode)
- 143 pairs ≈ 5 hours total
- Peak VRAM: ~4 GB (monitored and kept within safe bounds)

**To test with fewer pairs:**
```bash
python scripts/01_fold_benchmark.py --limit 5  # Just fold first 5
```

## Step 2: Extract Scores (~10 seconds)

Parse ColabFold JSON outputs and merge with outcome labels.

```bash
python scripts/02_extract_scores.py
```

**What it does:**
- Reads all `*_scores_rank_001_*.json` from `results/benchmark_folds_raw/`
- Extracts: ipTM, pTM, mean_pLDDT, ipSAE, pDockQ2
- Matches each pair to its outcome ('true', 'false', or 'unknown') from benchmark CSV
- Produces `results/benchmark_scores_with_outcomes.csv`

**Expected output:**
```
Loaded 143 benchmark pairs
Using 142 pairs with confirmed outcomes
✓ Extracted 142 scores
✓ Saved to results/benchmark_scores_with_outcomes.csv

Columns: pair, sequence, target_seq, outcome, source_system, ipTM, pTM, mean_pLDDT, ipSAE, pDockQ2
Sample:
     pair     sequence  ...   ipTM  pTM  mean_pLDDT  ...
 0 pair_000  MQQSVQTQ...  0.067  0.25        38.5  ...
```

## Step 3: Fit Calibration Curve (~30 seconds + plot generation)

Bin by ipTM and compute empirical hit rates.

```bash
python scripts/03_fit_curve.py --plot
```

**What it does:**
- Groups pairs by ipTM (bins: 0.0–0.05, 0.05–0.10, etc.)
- For each bin: counts "outcome=true" pairs and computes binomial hit rate + 95% CI
- Saves bins to `results/calibration_model.json`
- Optionally plots histogram + CI shading to `results/calibration_curve.png`

**Expected output:**
```
Loaded 143 scored pairs
Using 142 pairs with confirmed outcomes
✓ Saved calibration model to results/calibration_model.json
✓ Saved plot to results/calibration_curve.png

✓ Calibration curve complete: 18 bins

Sample bins:
  ipTM 0.05–0.10: 5/24 hit (20.8%, 95% CI [7.0%, 42.6%])
  ipTM 0.10–0.15: 8/31 hit (25.8%, 95% CI [12.8%, 42.2%])
  ipTM 0.15–0.20: 7/28 hit (25.0%, 95% CI [12.0%, 41.9%])
```

**Interpretation:**
- Hit rate increases (roughly) with ipTM
- Wide CIs in sparse bins = honest uncertainty
- Flat or non-monotonic regions = real data, not a failure

## Step 4: Apply to Your Designs

Estimate hit rates for your own de novo binders.

```bash
# Use the example (your top 5 PD-L1 candidates)
python scripts/04_calibrate_your_designs.py examples/example_user_designs.csv
```

**Expected output:**
```
Loaded 5 designs from examples/example_user_designs.csv
✓ Calibrated 5 designs
✓ Saved to results/example_user_designs_calibrated.csv

Sample (first 3 rows):
           sequence   ipTM estimated_hit_rate ci_lower ci_upper   notes
 0  MQQSVQTQFQ...  0.170             0.350    0.171    0.551
 1  AFTVTVPKDL...  0.160             0.270    0.123    0.431
 2  NPPTSALVVT...  0.150             0.250    0.115    0.416
```

**To apply to your own file:**
1. Prepare a CSV with columns: `sequence, target, ipTM, [optional: pTM, mean_pLDDT, ...]`
2. Run: `python scripts/04_calibrate_your_designs.py your_designs.csv`
3. Output: `results/your_designs_calibrated.csv`

## Step 5 (Optional): Interpret and Iterate

- **High ipTM (>0.2):** Designs here typically have 30–50% hit rates in the benchmark. Good candidates for synthesis.
- **Low ipTM (<0.1):** Hit rates drop to <20%. Use only if you've exhausted high-ipTM options.
- **Confidence intervals:** Narrow CIs (e.g., ipTM 0.3–0.4) = many benchmark samples = confident rate. Wide CIs (ipTM 0.8+) = few samples = uncertain.
- **Note field:** Warns if target is not in the benchmark or ipTM is out of range.

## Troubleshooting

### "colabfold_batch not found"
```bash
bash environment/setup.sh
source activate colabfold-env
# Then retry
```

### "CUDA out of memory"
- Close other GPU applications (Chrome, Slack, etc.)
- Reduce batch size (edit `01_fold_benchmark.py`, line ~100: change `TF_FORCE_UNIFIED_MEMORY`)
- Use CPU (slow): set `CUDA_VISIBLE_DEVICES=-1` (not recommended)

### "Benchmark file not found: benchmark/adaptyv_benchmark.csv"
- Check you're in the `calibration-kit/` directory
- Run `ls benchmark/` to verify the file exists
- If missing, download from [link to data] or use a smaller test subset

### "No score JSON found for pair_XXX"
- Folding may have crashed for that pair. Check `results/benchmark_folds_raw/` for files with that prefix.
- Re-run `01_fold_benchmark.py --resume` to retry interrupted pairs
- If the error persists, check `nvidia-smi` for OOM (out-of-memory) errors

### "Wide confidence intervals"
- This is normal! Bins with <10 samples have wide CIs. It's better to be honest about uncertainty.
- As you fold more designs or add your own labeled data, CIs will narrow.

### "Hit rate is flat across all bins"
- Indicates the metric (ipTM) has poor discriminative power for your target
- Check the plot: if histograms for "true" vs "false" outcomes heavily overlap, this is real data
- Consider using a different scoring metric or collecting more labeled data for your specific target

## Next Steps

1. **Synthesize and validate:** Pick top candidates (high ipTM + low CI) for wet-lab assays
2. **Collect labels:** Record which designs actually bind
3. **Update the curve:** Add your labeled pairs to `benchmark/adaptyv_benchmark.csv` and re-run steps 2–3
4. **Improve generalization:** Share your results (PR) to help the community

## For More Information

- **Methodology & caveats:** See `RESEARCH_FINDINGS.md`
- **API details:** See docstrings in each script (e.g., `head -50 scripts/01_fold_benchmark.py`)
- **Example output:** Check `examples/example_user_designs_calibrated.csv`

---

**Questions?** Open an issue on GitHub or check the FAQ in `RESEARCH_FINDINGS.md`.
