# Glossary

Terms used across this project, from molecular basics up to the metrics and
statistics the pipeline reports. Where a term has a concrete value in our own
data, that value is given — abstract definitions are harder to hold onto than
real numbers.

---

## 1. Molecular basics

**Protein** — a chain of amino acids that folds into a specific 3D shape. The
shape determines the function. Essentially all of biology's machinery is protein.

**Amino acid** — the 20 building blocks proteins are made from. Each has a
common backbone part and a distinct **side chain** that gives it character:
charged (K, R, D, E), hydrophobic (L, I, V, F), small (A, G), reactive (C), etc.

**Residue** — one amino acid *once it is inside a chain*. "Residue 77" means the
77th amino acid along. Used interchangeably with "position."

**One-letter code** — each amino acid has a letter: A=alanine, C=cysteine,
K=lysine, and so on. A sequence like `MEKKEEIKK...` is a protein written out.

**Backbone (chemistry)** — the repeating N–Cα–C=O atoms running along the chain,
identical in every amino acid. The side chains hang off it. When RFdiffusion
"designs a backbone," it is placing this chain in 3D without deciding the side
chains yet.

**N-terminus / C-terminus** — the two ends of a protein chain. Sequences are
written N→C, left to right.

**Chain** — one continuous protein molecule. A complex has several, labelled
A, B, E… in structure files.

**Fold** — the overall 3D shape a sequence adopts. Many different sequences can
share a fold.

**Secondary structure** — local repeating shapes. **α-helix** (a coil, very
common in designed binders) and **β-sheet** (extended strands side by side).
Alanine is the strongest helix-former, which is why a 40%-alanine design still
folded confidently in silico.

**Domain** — a self-contained folding unit within a larger protein. IL-7Rα's
**ectodomain** is the part outside the cell, which is the part we target.

**Hydrophobic / hydrophilic** — water-avoiding vs water-liking. Hydrophobic
residues normally bury themselves in a protein's core; when they sit exposed on
the surface, the protein tends to aggregate — a real liability for a designed
binder.

**Disulfide bond** — a covalent link between two **cysteine** (C) residues.
Useful when intended, a problem when not: unpaired cysteines cause scrambling
and aggregation. Standard practice is to exclude C from designed binders
(`--omit_AAs C` in ProteinMPNN).

**Glycosylation** — sugars attached to a protein, often at an **N-X-S/T sequon**
(asparagine, any residue, then serine or threonine). Affects expression.

---

## 2. Structures and files

**PDB** — the Protein Data Bank, the public archive of experimentally determined
structures; also the file format. A **PDB ID** is a 4-character code like `3DI3`.

> A PDB ID identifies a specific structure, not a protein. `1ILR` is
> interleukin-1 receptor *antagonist*; `3DI3` is the interleukin-**7** receptor α
> complex. Confusing the two cost this project a full campaign — always check the
> `TITLE` line.

**Crystal structure** — a structure determined experimentally (usually X-ray
crystallography). Ground truth, as opposed to a prediction.

**Chain ID** — the letter labelling each chain in a structure file. Conventions
vary per campaign here: RBD used target=E / binder=A; PD-L1 used target=A /
binder=B; IL-7Rα used target=B / binder=A.

**Residue numbering** — structures keep the *biological* numbering, which rarely
starts at 1. IL-7Rα chain B runs 17–209.

**Ångström (Å)** — 0.1 nanometres. A carbon atom is ~1.5 Å wide; a small protein
~30 Å across.

**Complex** — two or more chains bound together. AlphaFold-multimer predicts these.

**Interface** — the surface where two chains touch. **Epitope** is the patch on
the *target* side specifically.

**Buried surface area** — how much surface gets hidden when two proteins bind.
Bigger generally means a stronger, more specific interaction.

---

## 3. Binding

**Ligand / receptor** — the small partner and the big one. IL-7 is the ligand;
IL-7Rα is the receptor.

**Binder** — any molecule designed or found to stick to a target.

**Affinity** — how tightly two molecules bind, usually as **K_D** (dissociation
constant) in molar units. **Lower K_D = tighter binding.** 1 nM (10⁻⁹ M) is
strong; 1 µM (10⁻⁶ M) is weak. The Adaptyv EGFR data uses K_D as its outcome.

> **ipTM is not affinity.** Every number in this project is a *prediction
> confidence*, not a measured interaction. No binding assay has been run here.

**Specificity** — binding the intended target and *not* others. Measured here as
the margin between on-target and worst off-target ipTM. Our candidates ranged
from 0.15 to 0.59 — a real spread.

**Off-target / decoy** — an unrelated protein used to test specificity. We used a
panel (PD-L1 and RBD), because candidates failed on *different* decoys — one
bound PD-L1 at 0.61 but RBD at 0.15, another the reverse.

---

## 4. The proteins in this project

**IL-7Rα** (interleukin-7 receptor α) — immune receptor controlling T-cell
development. Target of the current campaign. PDB **3DI3** chain B, residues
17–209.

**IL-7** — its natural ligand. Chain A of 3DI3; the residues it contacts became
our hotspots.

**SARS-CoV-2 Spike RBD** — the receptor-binding domain of the COVID virus's spike
protein, the part that grabs human **ACE2** to enter cells. PDB **6M0J** chain E,
residues 333–526. The "hard" target in our variance study.

**ACE2** — the human receptor the RBD binds. Used as a control.

**PD-L1 / PD-1** — immune checkpoint pair; blocking their interaction is the
mechanism of several cancer drugs. PD-L1 is PDB **5O45** chain A, residues 18–132.

**EGFR** — epidermal growth factor receptor, a major cancer target. The Adaptyv
competition data covers it.

**IL-1Ra** — interleukin-1 receptor antagonist, PDB `1ILR`. Unrelated to anything
here. Included only because we mistakenly designed against it for a full round.

---

## 5. Prediction

**AlphaFold2 (AF2)** — DeepMind's structure predictor. **AlphaFold-multimer** is
the variant for complexes. `alphafold2_multimer_v3` is the version used here.

**ColabFold** — a fast, accessible wrapper around AlphaFold2. `colabfold_batch`
is its command-line tool.

**MSA** (Multiple Sequence Alignment) — AlphaFold's main information source: a
stack of evolutionarily related sequences from across species. Positions that
mutate *together* are probably touching in 3D — this **co-evolution** signal is
most of what AlphaFold reads.

**a3m** — the file format an MSA is stored in.

**single-sequence mode** — running with no MSA. Fast, no database search.

> This is the project's central mechanical finding. A *designed* binder has no
> evolutionary relatives, so an MSA adds nothing to it — binder pLDDT was 70.1
> without an MSA and 72.3 with one. But the **target** is a natural protein with
> hundreds of relatives, and stripping its MSA dropped it from pLDDT 84.3 to
> **27.2** — random coil. The fast screen was scoring interfaces against an
> unfolded target, which is why it carried no signal.

**Template** — supplying an actual 3D structure for AlphaFold to work from, as an
alternative to an MSA.

**Initial guess** — the field-standard version of that idea for binder design
(Bennett et al. 2023): seed AlphaFold with the designed complex's coordinates.
Released in the `dl_binder_design` repository.

**Recycle** — AlphaFold feeds its own output back in to refine. 3 recycles =
quick, 8 = thorough.

**Model 1–5** — five independently-trained AlphaFold networks. Running several
gives a reliability check: our top candidate scored 0.88 / 0.86 / 0.85 / 0.85 /
**0.31** — four agree, one dissents.

**Amber relaxation** — a physics cleanup pass fixing clashing atoms and strained
bonds. Runs *after* scoring, so it doesn't change the metrics. GPU relaxation
fails on this setup ("Minimization failed after 100 attempts", an OpenMM/CUDA
issue under WSL2); CPU relaxation works in ~17 s.

---

## 6. Design tools

**De novo design** — building a protein from scratch rather than discovering one
in nature.

**RFdiffusion** — a diffusion model that generates protein **backbones**. Like an
image diffusion model, it starts from noise and denoises into a structure —
except the output is 3D coordinates conditioned on your target's surface.

**Diffusion timestep / `partial_T`** — how far back toward noise you push an
existing structure before re-denoising. Small values produce *structural
neighbours* of a known-good backbone, which is what **partial diffusion**
refinement does.

**ProteinMPNN** — **inverse folding**: given a backbone shape, predict which
amino acid sequence would fold into it. The step that was collapsing into
poly-alanine.

**Sampling temperature** — how randomly ProteinMPNN picks residues. 0.1 is nearly
greedy (always the most likely residue), 0.3 is more varied. Raising it did *not*
fix the alanine collapse; a logit bias did.

**Logit bias** — pushing a specific amino acid's probability up or down directly.
`{"A": -1.0}` made alanine ~2.7× less likely, cutting mean alanine content from
20.6% to 8.8% (natural abundance is ~8%).

**Contig** — RFdiffusion's specification of what to build.
`B17-209/0 50-120` = "hold chain B residues 17–209 fixed, chain break, then
design a new 50–120 residue chain."

**Hotspot residues** — target residues you steer the designed binder toward.
Ours: S31, K77, K138, Y192 on IL-7Rα.

**Negative design / multi-state design** — optimising a sequence to fold well on
one structure *and badly* on another. ProteinMPNN supports this via
`--pos_neg_chain_betas` (positive weight for the target, negative for the decoy).
The principled fix for poor specificity, versus merely filtering for it after.

---

## 7. Metrics

**pLDDT** (0–100, per residue) — local confidence. >90 very high, 70–90
confident, 50–70 low, **<50 effectively disordered**. Report it **per chain**:
averaging across a complex lets a well-folded target mask a bad binder (our
target is 70% of the residues).

**pTM** (0–1) — confidence in the overall fold of the whole complex.

**ipTM** (0–1) — confidence in the **interface** specifically: how the chains sit
relative to each other. The primary metric used throughout. >0.8 is conventionally
a confident interface.

**PAE** (Ångströms, lower better) — Predicted Aligned Error. A matrix: for each
residue pair, the expected positional error. The **cross-chain block** is the
interface PAE — our validated candidate's was 5.8 Å.

**ipSAE / pDockQ / pDockQ2** — newer interface-quality scores, present in the same
output. **pDockQ2 discriminated our controls far better than ipTM**: 0.61
on-target vs 0.013 scrambled (47×), where ipTM managed only 3.3×.

**Shannon entropy** (bits) — sequence diversity. A random 20-amino-acid sequence
scores ~4.32; natural proteins ~4.0–4.1; our designs averaged 3.01, and the
degenerate 40%-alanine one 3.15.

---

## 8. Statistics

**Spearman ρ** (−1 to +1) — rank correlation. Do two measurements order things
the same way? Our fast-screen vs high-accuracy within-backbone ρ was **−0.14**:
no usable ordering.

**p-value** — probability of seeing a result this extreme if nothing real were
happening. <0.05 is the conventional threshold. Our fast-screen null returned
**p = 0.970**.

**Mann–Whitney U** — tests whether two groups differ without assuming a normal
distribution. Used for top-10 vs bottom-10.

**Bootstrap** — resampling data many times to estimate behaviour under
repetition. Used for the best-of-N curve. Resampling *with replacement* from only
6 draws biases the N=6 endpoint **low**.

**Variance decomposition** — splitting total variation by source.
"78% within-backbone" means 78% of RBD's ipTM variation came from *which sequence
was drawn* and only 22% from *which backbone was used*.

**Range restriction** — correlations computed inside a narrow band look weaker
than they are. Worth naming because our first fast-screen result was measured
only within the top 10 — the follow-up with top *and* bottom deciles removed the
concern.

---

## 9. Reinforcement learning

**Policy** — the thing being learned. Here, a table of per-position amino acid
preferences.

**Reward** — the number being maximised. Originally raw ipTM; now ipTM minus a
composition penalty.

**REINFORCE** — a basic policy-gradient algorithm: sample, score, push the policy
toward what scored well.

**Reward shaping** — adding terms to steer behaviour, e.g. subtracting a penalty
for poly-alanine so the policy can't win by degenerating.

**Reward channel** — the measurement the reward is computed from. If the channel
can't rank, the optimiser is climbing noise. Ours was fast-screen ipTM, and
selecting by it performed **worse than random** (−0.105 ipTM per backbone).

**Best-of-N** — draw N candidates, keep the best. Simple, parallel, no
hyperparameters — and given how expensive and noisy our evaluations are, likely a
better use of budget than policy gradient.

---

## 10. Project-specific shorthand

**Fast screen** — 1 model, 3 recycles, no MSA. Cheap, and shown here to carry no
selective signal.

**High accuracy** — 3 models, 8 recycles, full MSA search. ~4.5 min/fold, used as
ground truth.

**Hybrid MSA** — our variant: precomputed target alignment, binder single-sequence.
Fast-screen settings otherwise, and no MSA-server calls.

**Composition penalty** — a score subtracted for low-complexity sequence
(homopolymer runs, one residue dominating, low entropy). 0.000 on natural
proteins, 0.307 on the 40%-alanine design.

**Shaped reward** — ipTM minus the composition penalty.

**Specificity margin** — on-target ipTM minus the *worst* off-target ipTM. Max,
not mean: binding one wrong thing is the failure mode.

**Monomer foldability** — does the binder hold its shape alone, without the
target scaffolding it? If not, it will likely be disordered in solution.

**Solo verification** — re-folding a surprising result in isolation, to rule out
batch artifacts. A project rule after a batched score of 0.771 collapsed to 0.164
when re-run alone.
