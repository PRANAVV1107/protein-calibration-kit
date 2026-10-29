#!/usr/bin/env python3
"""Pre-registered test: how often does uniform beat the oracle?

Method is fixed in reports/ORACLE_BOOTSTRAP_PREREGISTRATION.md, committed
before this script was run. Do not alter the method after seeing output.

RNG: numpy.random.default_rng(0) (the PCG64 generator), declared explicitly
since "seed 0" alone is ambiguous across RNG families and this must be
exactly reproducible by someone re-running it.
"""
import numpy as np
import pandas as pd

df = pd.read_csv("results/root_breadth_per_root.csv")
assert len(df) == 50, f"expected 50 backbones, got {len(df)}"

# Tie-break key: lowest backbone index, declared in advance.
df["idx"] = df["backbone"].str.replace("binder_", "", regex=False).astype(int)
df = df.reset_index(drop=True)

means = df["mean"].to_numpy()
maxs = df["max"].to_numpy()
idxs = df["idx"].to_numpy()

N_ROOTS = len(df)
SUBSET = 10
N_REPS = 10_000
SEED = 0

rng = np.random.default_rng(SEED)

def argbest_tiebreak(values, sub_idx, sub_backbone_idx):
    """argmax over `values[sub_idx]`, ties broken by lowest backbone index."""
    v = values[sub_idx]
    b = sub_backbone_idx
    best_v = v.max()
    candidates = sub_idx[v == best_v]
    cand_b = b[v == best_v]
    return candidates[np.argmin(cand_b)]

discordant = np.zeros(N_REPS, dtype=bool)
gaps = np.zeros(N_REPS, dtype=float)
oracle_vals = np.zeros(N_REPS, dtype=float)
uniform_vals = np.zeros(N_REPS, dtype=float)

for i in range(N_REPS):
    sub = rng.choice(N_ROOTS, size=SUBSET, replace=False)
    sub_b = idxs[sub]

    best_mean_root = argbest_tiebreak(means, sub, sub_b)
    best_max_root = argbest_tiebreak(maxs, sub, sub_b)

    discordant[i] = best_mean_root != best_max_root
    oracle_val = maxs[best_mean_root]          # oracle capped at its root's own max
    uniform_val = maxs[sub].max()              # subset's true max (uniform's ceiling)

    oracle_vals[i] = oracle_val
    uniform_vals[i] = uniform_val
    gaps[i] = oracle_val - uniform_val

discordance_rate = discordant.mean()
median_gap = np.median(gaps)
q1, q3 = np.percentile(gaps, [25, 75])

headline = -0.099
pct_le_headline = (gaps <= headline).mean()          # fraction at least as bad as headline
pct_rank = (gaps <= headline).mean() * 100           # percentile position of -0.099 (lower tail)

print("=" * 70)
print("ORACLE BOOTSTRAP -- pre-registered, 10,000 reps, seed=0, 10-of-50 draws")
print("=" * 70)
print(f"N replicates:              {N_REPS}")
print(f"Discordance rate:          {discordance_rate:.4f}  ({discordance_rate*100:.2f}%)")
print(f"Gap (oracle - uniform):")
print(f"  median:                  {median_gap:+.4f}")
print(f"  IQR:                     [{q1:+.4f}, {q3:+.4f}]")
print(f"  mean:                    {gaps.mean():+.4f}")
print(f"  min / max:               {gaps.min():+.4f} / {gaps.max():+.4f}")
print(f"  fraction of reps with gap <= -0.099 (as bad as or worse than headline): {pct_le_headline:.4f} ({pct_le_headline*100:.2f}%)")
print()
# Percentile of the headline gap, expressed in SEVERITY terms (severity = -gap,
# so severity=0 is best/no-loss and larger severity is worse). This avoids the
# sign confusion of "percentile of a negative number": we ask what fraction of
# replicates are AT LEAST AS GOOD as the headline, and the complement is the
# headline's percentile rank in the severity distribution.
severity = -gaps
headline_severity = -headline
frac_better_or_equal = (severity <= headline_severity).mean()
severity_percentile = frac_better_or_equal * 100
print(f"Headline severity (0.099) percentile within the simulated severity distribution: {severity_percentile:.2f}th percentile")
print(f"(i.e. only {100-severity_percentile:.2f}% of replicates are as bad as or worse than the headline -0.099)")
print()

# Among discordant reps only, how bad is the gap? (conditional severity)
if discordant.sum() > 0:
    cond_gaps = gaps[discordant]
    print(f"Among discordant reps only (n={discordant.sum()}):")
    print(f"  median gap: {np.median(cond_gaps):+.4f}   IQR [{np.percentile(cond_gaps,25):+.4f}, {np.percentile(cond_gaps,75):+.4f}]")
print()

# Decision rule application
print("-" * 70)
print("DECISION RULE (fixed in advance)")
print("-" * 70)
if discordance_rate >= 0.50 and median_gap <= -0.05:
    verdict = "STANDS -- headline generalizes; report distribution alongside it"
elif 0.20 <= discordance_rate < 0.50:
    verdict = "CONDITIONAL -- restate as: 'when best-mean root lacks the best design (X% of draws)'"
elif discordance_rate < 0.20 or severity_percentile >= 90:
    verdict = "WITHDRAWN as general claim -- one favourable draw; lead section 3.4 on roots-beat-draws instead"
else:
    verdict = "FALLS BETWEEN DECLARED BUCKETS -- report exact numbers, use judgment, do not force a bucket"
print(f"discordance_rate = {discordance_rate:.4f}  ({'>=0.50' if discordance_rate>=0.5 else ('0.20-0.50' if discordance_rate>=0.2 else '<0.20')})")
print(f"median_gap       = {median_gap:+.4f}  ({'<=-0.05' if median_gap<=-0.05 else '>-0.05'})")
print(f"headline pctile  = {severity_percentile:.2f}  ({'>=90th' if severity_percentile>=90 else '<90th'})")
print()
print("VERDICT:", verdict)

# save per-replicate data for audit
out = pd.DataFrame({
    "discordant": discordant,
    "oracle_val": oracle_vals,
    "uniform_val": uniform_vals,
    "gap": gaps,
})
out.to_csv("results/oracle_bootstrap_replicates.csv", index=False)
print()
print("[saved] results/oracle_bootstrap_replicates.csv (10,000 rows, full audit trail)")
