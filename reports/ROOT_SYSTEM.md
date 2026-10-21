# The root system

What a "root" actually is in this pipeline, how the 200 of them were made, how
much they differ, and whether any of that can be predicted without folding.

Written 2026-09-21. All numbers from `scripts/root_anatomy.py` and
`scripts/root_allocation.py`, both of which run on existing files with no GPU.

---

## 1. Where roots come from

A root is one RFdiffusion output: a backbone-only PDB with two chains — the
designed binder (chain A, no sequence yet) and the target held fixed (chain B,
RBD from 6M0J chain E, residues 333–526).

They were generated 2026-09-04 in **four batches of 50 = 200 backbones**
(`real_test_results/rbd_200_batch/batch_{0,1,2,3}/rfdiffusion/`), from the
`Complex_base_ckpt.pt` checkpoint, steered toward six hotspots on the RBD:
**E449, E493, E496, E498, E501, E505** — the ACE2-facing patch.

Nothing downstream re-generates roots. ProteinMPNN puts sequences on them;
AlphaFold scores the pair; the MCTS re-noises one of them. The 200 are the
entire supply.

## 2. What varies between roots

All 200 backbones, measured from the PDB files alone:

| feature | min | median | max | sd |
|---|---|---|---|---|
| binder length (aa) | 70 | 85 | 100 | 9.1 |
| radius of gyration (Å) | 10.9 | 15.3 | 30.4 | 2.9 |
| binder–target contacts (CA < 10 Å) | 83 | 160 | 254 | 28.9 |
| contacts to the six hotspots | 17 | 39 | 61 | 8.2 |
| helical fraction | 0.80 | 0.92 | 1.00 | 0.04 |
| contacts per residue | 1.04 | 1.85 | 2.83 | 0.36 |

**Length was sampled, not fixed.** The generation contig allowed a range, and
the result is near-uniform across 70–100 aa:

    70- 74  ############################ (28)
    75- 79  ################################# (33)
    80- 84  ################################# (33)
    85- 89  ############################ (28)
    90- 94  ################################# (33)
    95- 99  ################################### (35)
    100-104  ########## (10)

This matters for one specific reason: the MCTS reads its contig **from the root**
(`[E333-526/0 71-71]` for binder_47), so every descendant in a search tree is
locked to the root's length. A search from a 71-aa root can never produce a
90-aa binder. Length is a root-level choice that search cannot revisit.

**Everything is a helical bundle.** Helical fraction runs 0.80–1.00 with sd
0.04 — there is essentially no fold-type diversity in the supply. Whatever the
between-root spread in quality comes from, it is not a choice between helices
and sheets, because sheets were never generated.

## 3. How much root choice is worth

Measured on the high-accuracy instrument (variance study, 10 RBD backbones,
6 draws each — the same instrument the 0.870 bar came from):

| root | best-of-6 | E[best-of-6] |
|---|---|---|
| binder_5 | 0.860 | **0.827** |
| binder_42 | 0.870 | 0.814 |
| binder_9 | 0.770 | 0.764 |
| binder_13 | 0.810 | 0.694 |
| binder_0 | 0.730 | 0.694 |
| binder_34 | 0.710 | 0.660 |
| binder_15 | 0.710 | 0.658 |
| binder_10 | 0.540 | 0.469 |
| **binder_47** | **0.520** | **0.485** |
| binder_14 | 0.460 | 0.371 |

Spread of E[best-of-6]: **0.371 to 0.827**, sd 0.154. Choosing well among these
ten is worth about **0.46 ipTM** — larger than anything three searches recovered
inside one root.

**All three MCTS runs started from binder_47, which ranks 9th of 10.** It was
chosen arbitrarily, before any of this was measured.

### Allocation: breadth beats depth at every budget

E[best ipTM] for a budget of B folds split R roots × B/R draws, roots drawn at
random (so: what you get *without* knowing which root is good). Hybrid
instrument, matching the searches:

| budget | R=1 | R=2 | R=5 | R=10 |
|---|---|---|---|---|
| 10 | 0.365 | 0.400 | 0.433 | **0.435** |
| 20 | 0.417 | 0.493 | 0.545 | **0.563** |
| 40 | 0.456 | 0.572 | 0.650 | **0.676** |
| 80 | 0.485 | 0.615 | 0.720 | **0.747** |
| 160 | 0.469 | 0.632 | 0.750 | **0.774** |

Breadth wins everywhere and the margin grows — **+0.070 at B=10, +0.305 at
B=160** — because **single-root saturates** (flat past B≈80) while multi-root is
still climbing at the edge of the data.

## 4. Can a good root be spotted without folding it?

This is the question that would make root selection free. Six cheap structural
features against measured quality, n = 10:

| feature | Spearman ρ | p |
|---|---|---|
| contacts | 0.511 | 0.132 |
| contacts per residue | 0.480 | 0.160 |
| hotspot contacts | −0.355 | 0.314 |
| radius of gyration | −0.286 | 0.424 |
| helical fraction | 0.170 | 0.638 |
| **binder length** | **0.018** | **0.960** |

**Nothing reaches significance, and this is badly underpowered.** At n = 10,
|ρ| must exceed about 0.63 for p < 0.05 uncorrected, and testing six features
makes one apparent hit likely even if nothing is real. Treat all of it as
hypothesis-generating.

Two things are still worth saying:

- **Interface size is the only direction with any hint** (ρ ≈ 0.5 for both
  contact measures, neither significant). Mechanistically plausible — a bigger
  interface has more to score — and it is the one feature worth testing
  prospectively on backbones that played no part in choosing it.
- **Length predicts nothing** (ρ = 0.018), despite being the most obvious knob
  and the one thing a search irrevocably inherits from its root.
- **Hotspot engagement trends the wrong way** (ρ = −0.355). Not significant,
  but there is certainly no evidence that hitting the residues the generation
  was steered toward produces better binders.

So root quality is currently **only knowable by folding**. At 6 draws that is
6 folds per root, which the allocation table says is a good trade anyway — 10
roots × 6 draws = 60 folds reaches E[best] ≈ 0.816.

## 5. What this means for multi-root search

Running it is cheap to justify and the design is already written
(Experiment 4 in [EXPERIMENTS_QUEUED.md](EXPERIMENTS_QUEUED.md)). Three things
to build in, all of which come out of the numbers above:

1. **Draw roots without reference to score.** Selecting good roots first and
   then searching from them confounds the two effects; the comparison has to be
   against blind sampling over the *same* roots.
2. **Expect the gain to be concentrated in the bad roots.** Run 3 measured
   top-quartile parents at 0.342 producing best children at 0.279 — a loss —
   while average parents gained ~0.05. The search's headline +0.67 came from
   starting at the second-worst root available. Seeded from binder_5 there may
   be nothing to escape.
3. **Length is fixed per root, so a multi-root forest is also a length sweep.**
   That is a free extra axis, and since length predicts nothing on its own it
   is worth carrying rather than controlling away.

## 6. Limits

- Quality is measured on **10 of 200** backbones. The other 190 have never been
  folded at high accuracy.
- One target (RBD). PD-L1 and IL-7Rα roots exist but were not analysed here.
- The feature set is crude by design — no DSSP, no SASA, no Rosetta terms —
  because the point was features computable in milliseconds. A richer set might
  find what these missed.
- `helix_frac` is a CA(i)–CA(i+4) distance heuristic, not real secondary
  structure assignment.

---

Regenerate with:

    python scripts/root_anatomy.py        # features + the n=10 correlation
    python scripts/root_allocation.py     # the breadth-vs-depth table
    python scripts/root_sweep.py --target rbd
