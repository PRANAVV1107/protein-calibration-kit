# Pre-registration — guided search over denoising trajectories

**Written before the run produced any result.** Committed to git so the
timestamp is checkable. Nothing below may be revised after seeing outcomes; if
the analysis changes, the change is recorded as an amendment with its own
timestamp and reason.

---

## Why pre-register this

Search will almost certainly return *some* design scoring well. That is not
evidence the search worked, because blind sampling already does it: 60 random
folds across 10 RBD backbones reached ipTM 0.87 in the variance study. A result
of 0.85 would look like success in isolation while being **worse than doing
nothing** at comparable budget.

The failure mode this guards against is choosing the comparison after seeing
the number. Fixing the bar in advance is the only thing that makes a single
run interpretable.

---

## What is being run

Guided best-first search over RFdiffusion's denoising trajectory. A node is a
backbone; expansion re-noises it with `diffuser.partial_T=10` and re-denoises
into 4 structural neighbours; each node is evaluated by 4 ProteinMPNN sequences
folded with ColabFold (`alphafold2_multimer_v3`, 1 model, 3 recycles, target
alignment reused). Selection is UCT with c = 0.6, max depth 6.

- **Target:** SARS-CoV-2 RBD — the hard target here (single-draw mean ipTM
  0.400, against 0.766 and 0.805 for IL-7Rα and PD-L1)
- **Root:** `binder_15`, chosen before the run for having a *middling* measured
  baseline, so improvement has headroom to be visible
- **Budget:** 1000 folds
- **Metric:** ipTM, chosen to match the pre-existing baselines below
- **Hardware:** one RTX 3060, 6 GB

**This is not AlphaZero.** There is no value network — `value_probe.py`
measured the obvious free candidate (ProteinMPNN's own sequence likelihood) at
within-backbone ρ ≈ −0.06, meaning it cannot rank sequences on a fixed backbone.
Every leaf therefore costs a real fold and the tree stays shallow. It should be
described as guided best-first search, not MCTS-with-learned-value.

## Baselines, both measured before this run

| baseline | ipTM | source |
|---|---|---|
| root `binder_15`, best of 6 sequences | **0.71** | variance study, high accuracy |
| 60 folds blind across 10 RBD backbones | **0.87** | variance study, high accuracy |

The second is the one that matters: comparable fold budget, blind against
guided. The first only asks whether search improved on where it started.

> Neither baseline is perfectly matched. The variance-study folds used 3 models
> and 8 recycles; this search uses 1 model and 3 recycles, which scores
> somewhat lower on the same structure. That biases the comparison **against**
> the search, so clearing 0.87 would be a conservative result — and failing to
> clear it is correspondingly weaker evidence of failure. This asymmetry is
> stated now rather than discovered later.

---

## Primary outcome, fixed in advance

| best ipTM found within 1000 folds | conclusion |
|---|---|
| **> 0.87** | search beats blind sampling at matched budget — a real result |
| 0.71 – 0.87 | improves on its root, does not beat blind sampling — honest null |
| **≤ 0.71** | search did not help |

## Secondary outcome

**Does score improve with depth?** Mean node score by tree depth, and the
fraction of children that beat their own parent. If children systematically
exceed parents, the mechanism works even when the headline number does not
clear 0.87. Reported regardless of the primary outcome.

Chance is 50% for "child beats parent"; with ~250 children the standard error
is about 3%, so anything inside roughly 44–56% is not distinguishable from
noise.

---

## What this run cannot establish, whatever it returns

- **n = 1.** One root, one target, one seed, one hyperparameter setting. A
  single run cannot separate "search works" from "this root was lucky".
- **No wet-lab validation.** ipTM is a prediction. The only measured binding
  in this project is other groups' designs, re-scored.
- **Baselines are historical, not concurrent.** They were folded under a
  different setting (see the asymmetry note above).
- **No statistical test on the primary outcome.** It is a single number against
  a fixed threshold, not a comparison of distributions.

Claims must therefore stay at the level of *"in one run on one target"*.
Reporting a high absolute ipTM without the blind-sampling comparison beside it
would be misleading, because blind sampling already reaches 0.87.

## Language committed to in advance

Permitted, if the bar is cleared:

> In a single search run on one target, guided search over denoising
> trajectories reached ipTM X within 1000 folds, against 0.87 for blind
> sampling at comparable budget. One root, one target, no wet-lab validation.

Not permitted, whatever the outcome:

- any absolute score quoted without the 0.87 comparison alongside it
- "outperforms", "state of the art", or any claim of generality from n = 1
- calling it AlphaZero or MCTS-with-value-network
- reporting the secondary depth result as the headline if the primary fails

---

## Amendments

### 2026-09-16 — root changed, for a hardware reason, before any result

`binder_15` cannot be folded on this machine. Its binder chain is 96 aa, giving
a 290-residue complex, and a 6 GB card under WSL2 stalls on it: the GPU
allocation is carved out of host RAM and appears as `shared`, so a 7.6 GB host
is left with nothing and thrashes rather than failing. Two smoke attempts ran
14+ minutes each and produced zero folds.

The size ceiling was then measured rather than guessed. Everything this machine
has successfully folded is **282 residues or smaller**:

| evidence | complex size |
|---|---|
| reward benchmark, 85 successes | 221–282 |
| variance study RBD (75 aa binder) | 269 |
| `binder_15` (96 aa binder) | **290 — fails** |

This also explains the reward benchmark's OOM at exactly 85/165: ColabFold
folds shortest-first, and 85 is where it ran out of designs under the ceiling.
That was read at the time as a batch-size problem; it was a size limit.

**New root: `binder_47`** — 71 aa binder, 265-residue complex, comfortably
inside the ceiling. Chosen for being middling with headroom, the same criterion
originally applied to `binder_15`.

Baselines update accordingly:

| baseline | ipTM | note |
|---|---|---|
| root `binder_47`, best of 6 | **0.52** | was 0.71 for `binder_15` |
| 60 folds blind across 10 RBD backbones | **0.87** | **unchanged** |

The primary bar is unchanged at 0.87, which is the comparison that matters and
the one a lower root baseline cannot flatter. The root's own best-of-6 falling
from 0.71 to 0.52 does make "improves on its root" easier to clear, so that
secondary band is correspondingly weaker evidence than originally written, and
should be reported as such.

No results existed for either root when this amendment was made.
