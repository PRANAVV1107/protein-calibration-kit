#!/bin/bash
# Family-decoy specificity panel for the IL-7Ra candidates.
#
# The decoys already used (PD-L1, RBD) were chosen by what happened to be on
# disk. They are unrelated folds, so they test the EASY case. The informative
# decoys are structural homologs: IL-7Ra is a type-I cytokine receptor, and a
# binder that also hits the common gamma chain or IL-4Ra is a real specificity
# failure with clinical consequence. One that hits a viral RBD is a curiosity.
#
# Panel (all verified by reading each PDB's own TITLE/COMPND records):
#   gammaC   191 res  2B5I chain C  common gamma chain, IL-7Ra's own partner
#   IL4Ra    188 res  1IAR chain B  same family fold
#   IL21R    209 res  3TGX chain C  same family fold
#   IgGFc    207 res  1FC1 chain A  abundance / nonspecific-binding indicator
#   albumin  578 res  1AO6 chain A  most abundant serum protein
#
# Note on albumin: 578 vs 193 residues is a 3:1 ratio, the regime where ipTM is
# known to misbehave (this project hit that independently with ACE2/RBD). Score
# albumin pairings on ipSAE, not ipTM.
#
# Stage 1 pays for one MSA per decoy (5 high-accuracy folds). Stage 2 reuses
# those alignments through the hybrid route for every remaining pairing, at
# ~50 s each and no further server calls.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
D=$KIT/benchmark/decoys
LOG=$KIT/results/queue_v6.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }

{
say "=== QUEUE V6 (family decoys) -- waiting for earlier queues ==="
while pgrep -f 'queue_v4.sh' > /dev/null || pgrep -f 'queue_v5_null.sh' > /dev/null; do sleep 120; done
while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
say "GPU free."

# ---- stage 1: one high-accuracy fold per decoy, to generate its MSA ----
mkdir -p "$Q/family_seed"
say "START family_seed (5 folds, one per decoy, generates the alignments)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3 \
  "$Q/family_seed.fasta" "$Q/family_seed" > "$Q/family_seed.log" 2>&1 \
  && say "DONE  family_seed ($(ls "$Q/family_seed"/*.a3m 2>/dev/null | wc -l) alignments)" \
  || say "FAIL  family_seed (see log)"

# ---- stage 2: build hybrid inputs for every candidate x decoy pairing ----
say "Building hybrid inputs from the new alignments..."
cd "$KIT"
python scripts/build_family_panel.py >> "$LOG" 2>&1 || say "FAIL  build_family_panel"

# ---- stage 3: fold the panel ----
mkdir -p "$Q/family_panel_folds"
say "START family_panel"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$Q/family_panel" "$Q/family_panel_folds" > "$Q/family_panel.log" 2>&1 \
  && say "DONE  family_panel ($(ls "$Q/family_panel_folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l))" \
  || say "FAIL  family_panel (partial results kept)"

say "=== QUEUE V6 COMPLETE ==="
} >> "$LOG" 2>&1
