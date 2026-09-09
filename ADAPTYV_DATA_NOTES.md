# Adaptyv Data Integration Notes

## Status

**Adaptyv EGFR competition data acquired and processed:**
- EGFR Round 1: 202 sequences (8 with KD measurements)
- EGFR Round 2: 402 sequences (400 with ipTM scores, outcomes)
- **Total: 604 sequences with 434 confirmed binders, 295 non-binders**

Saved as: `benchmark/adaptyv_data/` (not merged into main benchmark, see below)

---

## Why Not Merged Into Main Calibration

The Adaptyv competition data uses **different ColabFold settings** than the campaign calibration:

| Dataset | ipTM Mean | Scoring Settings | Notes |
|---------|-----------|------------------|-------|
| **Campaign (RBD + PD-L1)** | 0.095 | Fast-screen: single-sequence, 1 model, 3 recycles | Standardized, low-VRAM |
| **Adaptyv EGFR** | 0.696 | Likely high-accuracy or ensemble | Competition submissions, higher-effort scoring |

**Problem:** Direct merging would confound ipTM thresholds. A design at ipTM 0.5 in Adaptyv settings might be ipTM 0.1 in fast-screen — they're not on the same scale.

**Solution:** Keep them separate in the current kit, with clear instructions for users who want to include Adaptyv data.

---

## How to Use Adaptyv Data for Your Own Calibration

### Option A: Rescore Adaptyv With Your Fast-Screen Settings

If you want to calibrate on Adaptyv data using the same ColabFold settings as the campaign:

```bash
# 1. Extract Adaptyv sequences
python scripts/extract_adaptyv_sequences.py

# 2. Fold them with fast-screen ColabFold
python scripts/01_fold_benchmark.py --limit 604

# 3. Extract scores
python scripts/02_extract_scores.py

# 4. Merge with outcomes
python scripts/merge_adaptyv_outcomes.py

# 5. Refit curve
python scripts/03_fit_curve.py --plot
```

### Option B: Use Adaptyv Scores As-Is (Expert Only)

If you want to trust Adaptyv's high-accuracy ipTM scores directly:

1. Resample/bin ipTM using Adaptyv's higher resolution (e.g., 0.1 instead of 0.05)
2. Compute hit rates per bin
3. Note the non-equivalence with campaign calibration in your documentation

### Option C: Build Separate Calibrations

Recommended for production: maintain separate calibration curves for each scoring regime:
- `calibration_model_fast_screen.json` (current, from campaign)
- `calibration_model_high_accuracy.json` (from Adaptyv or other high-accuracy data)

Users can then choose based on their ColabFold settings.

---

## Files

**Adaptyv raw data (from GitHub):**
- `adaptyv_egfr_round1.csv` — 202 sequences, sparse KD data
- `adaptyv_egfr_round2.csv` — 402 sequences, 400 with ipTM + outcomes

**Processed:**
- `benchmark/adaptyv_benchmark_new.csv` — merged EGFR R1+R2 with normalized outcomes

---

## Next Steps

1. **For a more robust calibration:** Rescore 604 Adaptyv sequences with your fast-screen settings (4–5 hours on single GPU) and merge
2. **For production use:** Build separate calibration curves for different ColabFold settings and let users choose
3. **For research:** Compare calibration quality (fast-screen vs high-accuracy vs Adaptyv scoring) to understand metric sensitivity

---

## References

- **EGFR Competition Data:** https://github.com/adaptyvbio/egfr_competition_1 and egfr_competition_2
- **Paper:** [Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2) (bioRxiv, April 2025)

