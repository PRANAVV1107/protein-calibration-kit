#!/bin/bash
# Overnight GPU queue. Runs jobs strictly in sequence -- concurrent ColabFold
# on a 6GB card has OOM'd before.
#
# Ordered by what each job decides, most decisive first:
#   1. reward-channel test, RBD   fast screen, no MSA server  ~2.5 h
#   2. RBD fast-screen null       high accuracy, needs MSA    ~1.5 h
#   3. PD-L1 variance completion  high accuracy, needs MSA    ~1.0 h
#   4. reward-channel test, IL-7Ra fast screen, no MSA server ~2.5 h
#
# Jobs 1 and 4 use single_sequence mode and make no MSA server calls, so they
# are immune to the throttling that stalled the variance study. They are placed
# around the MSA-dependent jobs deliberately.
#
# A failing job does not stop the queue: partial results are still analysable,
# and a throttled MSA job should not block the two that need no server.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/overnight_queue.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

# Match the interpreter path, not the bare binary name: a watcher running
# `pgrep -f colabfold_batch` carries that string in its own command line, so a
# bare-name match makes watchers match each other and deadlock.
REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'

say(){ echo "[$(date '+%F %T')] $*"; }

wait_for_gpu(){
  if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
    say "waiting for in-flight fold to finish..."
    while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
    say "GPU free."
  fi
}

fast(){   # name, fasta, outdir
  say "START $1  (fast screen, no MSA server)"
  colabfold_batch --msa-mode single_sequence --model-type alphafold2_multimer_v3 \
    --num-models 1 --num-recycle 3 "$2" "$3" > "$Q/$1.log" 2>&1 \
    && say "DONE  $1  ($(ls "$3"/*_scores_rank_001_*.json 2>/dev/null | wc -l) folds)" \
    || say "FAIL  $1 (continuing; partial results in $3)"
}

acc(){    # name, fasta, outdir
  say "START $1  (high accuracy, uses MSA server)"
  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
    --num-models 3 --num-recycle 8 "$2" "$3" > "$Q/$1.log" 2>&1 \
    && say "DONE  $1  ($(ls "$3"/*_scores_rank_001_*.json 2>/dev/null | wc -l) folds)" \
    || say "FAIL  $1 (continuing; partial results in $3)"
}

{
say "=== OVERNIGHT QUEUE START ==="
wait_for_gpu

mkdir -p "$Q/rewardchannel_rbd" "$Q/rbd_null" "$Q/rewardchannel_il7ra"

# 1 -- decides whether the RL reward signal is valid at all
fast rewardchannel_rbd   "$Q/rewardchannel_rbd.fasta"   "$Q/rewardchannel_rbd"

# 2 -- second independent target for the fast-screen null
acc  rbd_null            "$Q/rbd_null.fasta"            "$Q/rbd_null"

# 3 -- finishes the variance study to 10 backbones on all three targets
acc  pdl1_variance       "$KIT/results/variance_study/pdl1/to_fold.fasta" \
                         "$KIT/results/variance_study/pdl1/folds"

# 4 -- reward channel on the easy target, for contrast with RBD
fast rewardchannel_il7ra "$Q/rewardchannel_il7ra.fasta" "$Q/rewardchannel_il7ra"

# 5 -- monomer foldability: does each binder fold WITHOUT its target?
#      Designed sequences have no homologs, so single_sequence is the correct
#      setting here, not a shortcut. Small monomers, no MSA server.
mkdir -p "$Q/monomers"
say "START monomers  (183 monomer folds, fast screen, no MSA server)"
colabfold_batch --msa-mode single_sequence --model-type alphafold2_ptm \
  --num-models 1 --num-recycle 3 "$Q/monomers.fasta" "$Q/monomers" \
  > "$Q/monomers.log" 2>&1 \
  && say "DONE  monomers  ($(ls "$Q/monomers"/*_scores_rank_001_*.json 2>/dev/null | wc -l) folds)" \
  || say "FAIL  monomers (continuing; partial results in $Q/monomers)"

say "=== OVERNIGHT QUEUE COMPLETE ==="
} >> "$LOG" 2>&1
