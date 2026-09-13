#!/bin/bash
# Four-way screen comparison, then the remaining backlog.
#
# The same 60 RBD sequences are scored under conditions that differ in exactly
# one thing: how the TARGET's structure is supplied. The fast screen failed
# because without an MSA the target does not fold (pLDDT 27.2 vs 84.3), so the
# interface was being scored against a random coil.
#
#   fast screen   1 model, 3 recycles, nothing          MEASURED: pLDDT 27, rho -0.14
#   hybrid MSA    1 model, 3 recycles, target alignment  <- job 1
#   template      1 model, 3 recycles, target structure  <- job 2  (field-standard idea)
#   high accuracy 3 models, 8 recycles, full MSA         MEASURED: pLDDT 84, ground truth
#
# Jobs 1 and 2 make no MSA-server calls: the alignment is already on disk and
# the template is a local PDB.
#
# New file rather than an edit of the running queue -- bash reads a script by
# byte offset as it executes, so editing one mid-run can execute garbage.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v3.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

{
say "=== QUEUE V3 START ==="
if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
  say "waiting for in-flight fold (rbd_null) to finish..."
  while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
  say "GPU free."
fi

# ---- 1. hybrid MSA: target alignment, binder single-sequence ----
mkdir -p "$Q/hybrid_rbd_folds"
say "START hybrid_rbd (60, target MSA from disk)"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$Q/hybrid_rbd" "$Q/hybrid_rbd_folds" > "$Q/hybrid_rbd.log" 2>&1 \
  && say "DONE  hybrid_rbd ($(count "$Q/hybrid_rbd_folds"))" || say "FAIL  hybrid_rbd"

# ---- 2. template: target structure as a custom template ----
# --templates is required to enable --custom-template-path. Sequence input, so
# --msa-mode single_sequence keeps the binder MSA-free and the run server-free.
mkdir -p "$Q/template_rbd_folds" "$Q/template_cache"
say "START template_rbd (60, target PDB as template)"
colabfold_batch --msa-mode single_sequence --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3 \
  --templates --custom-template-path "$Q/templates_rbd" \
  --custom-template-cache-path "$Q/template_cache" \
  "$Q/rewardchannel_rbd.fasta" "$Q/template_rbd_folds" > "$Q/template_rbd.log" 2>&1 \
  && say "DONE  template_rbd ($(count "$Q/template_rbd_folds"))" || say "FAIL  template_rbd (see log)"

# ---- 3. reward channel on the easy target, contrast with RBD ----
mkdir -p "$Q/rewardchannel_il7ra"
say "START rewardchannel_il7ra (60)"
colabfold_batch --msa-mode single_sequence --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3 \
  "$Q/rewardchannel_il7ra.fasta" "$Q/rewardchannel_il7ra" > "$Q/rewardchannel_il7ra.log" 2>&1 \
  && say "DONE  rewardchannel_il7ra ($(count "$Q/rewardchannel_il7ra"))" || say "FAIL  rewardchannel_il7ra"

# ---- 4. finish the variance study to 10 backbones on all three targets ----
say "START pdl1_variance (12 remaining, uses MSA server)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 3 --num-recycle 8 \
  "$KIT/results/variance_study/pdl1/to_fold.fasta" \
  "$KIT/results/variance_study/pdl1/folds" > "$Q/pdl1_variance.log" 2>&1 \
  && say "DONE  pdl1_variance ($(count "$KIT/results/variance_study/pdl1/folds"))" || say "FAIL  pdl1_variance"

# ---- 5. monomer foldability ----
mkdir -p "$Q/monomers"
say "START monomers (183)"
colabfold_batch --msa-mode single_sequence --model-type alphafold2_ptm \
  --num-models 1 --num-recycle 3 \
  "$Q/monomers.fasta" "$Q/monomers" > "$Q/monomers.log" 2>&1 \
  && say "DONE  monomers ($(count "$Q/monomers"))" || say "FAIL  monomers"

say "=== QUEUE V3 COMPLETE ==="
} >> "$LOG" 2>&1
