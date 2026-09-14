#!/bin/bash
# Reordered queue: the hybrid-MSA test runs FIRST, because it is the most
# decisive open question (does a precomputed target MSA resurrect the cheap
# screen?) and costs only ~45 min.
#
# Written as a NEW file rather than editing the running queue: bash reads a
# script by byte offset while executing it, so editing one mid-run can make it
# jump into the middle of a line.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v2.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }

wait_gpu(){
  if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
    say "waiting for in-flight fold (rbd_null) to finish..."
    while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
    say "GPU free."
  fi
}

fast(){ # name, input, outdir  -- single_sequence, no MSA server
  say "START $1"
  colabfold_batch --msa-mode single_sequence --model-type alphafold2_multimer_v3 \
    --num-models 1 --num-recycle 3 "$2" "$3" > "$Q/$1.log" 2>&1 \
    && say "DONE  $1 ($(ls "$3"/*_scores_rank_001_*.json 2>/dev/null | wc -l))" \
    || say "FAIL  $1 (partial results kept)"
}

acc(){ # name, input, outdir -- full MSA, uses server
  say "START $1"
  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
    --num-models 3 --num-recycle 8 "$2" "$3" > "$Q/$1.log" 2>&1 \
    && say "DONE  $1 ($(ls "$3"/*_scores_rank_001_*.json 2>/dev/null | wc -l))" \
    || say "FAIL  $1 (partial results kept)"
}

{
say "=== QUEUE V2 START ==="
wait_gpu

# 1 -- HYBRID TEST (the decisive one). a3m inputs already carry the target MSA,
#      so no --msa-mode flag and no server calls. Same 1 model / 3 recycles as
#      the fast screen, so the ONLY changed variable is the target alignment.
mkdir -p "$Q/hybrid_rbd_folds"
say "START hybrid_rbd  (60 folds, fast-screen settings + precomputed target MSA)"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$Q/hybrid_rbd" "$Q/hybrid_rbd_folds" > "$Q/hybrid_rbd.log" 2>&1 \
  && say "DONE  hybrid_rbd ($(ls "$Q/hybrid_rbd_folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l))" \
  || say "FAIL  hybrid_rbd (partial results kept)"

# 2 -- reward channel on the easy target, for contrast with RBD
mkdir -p "$Q/rewardchannel_il7ra"
fast rewardchannel_il7ra "$Q/rewardchannel_il7ra.fasta" "$Q/rewardchannel_il7ra"

# 3 -- finish the variance study to 10 backbones on all three targets
acc  pdl1_variance "$KIT/results/variance_study/pdl1/to_fold.fasta" \
                   "$KIT/results/variance_study/pdl1/folds"

# 4 -- monomer foldability (183 small monomers, no MSA server)
mkdir -p "$Q/monomers"
say "START monomers (183 monomer folds)"
colabfold_batch --msa-mode single_sequence --model-type alphafold2_ptm \
  --num-models 1 --num-recycle 3 "$Q/monomers.fasta" "$Q/monomers" \
  > "$Q/monomers.log" 2>&1 \
  && say "DONE  monomers ($(ls "$Q/monomers"/*_scores_rank_001_*.json 2>/dev/null | wc -l))" \
  || say "FAIL  monomers (partial results kept)"

say "=== QUEUE V2 COMPLETE ==="
} >> "$LOG" 2>&1
