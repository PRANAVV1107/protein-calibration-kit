#!/bin/bash
# Hybrid-MSA test: fast-screen SETTINGS with a real target MSA.
#
# The fast screen failed because the target does not fold without an MSA
# (target pLDDT 27.2 vs 84.3). This isolates that single variable: identical
# model count and recycles as the fast screen (1 model, 3 recycles), the only
# difference being that the target carries its precomputed alignment.
#
#   fast screen    : 1 model, 3 recycles, NO msa           -> target pLDDT 27
#   hybrid  (this) : 1 model, 3 recycles, target MSA only  -> ?
#   high accuracy  : 3 models, 8 recycles, full MSA        -> target pLDDT 84
#
# The target MSA is precomputed once and reused, so this makes ZERO MSA-server
# calls despite using a real alignment.
#
# Runs as a separate process from overnight_queue.sh and waits for it, so the
# running queue is never edited mid-execution (doing that corrupts bash's
# byte offset and can execute garbage).
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/hybrid_test.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

# Full interpreter path, not the bare binary name: a watcher running
# `pgrep -f colabfold_batch` carries that string in its own command line.
REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'

say(){ echo "[$(date '+%F %T')] $*"; }

{
say "=== HYBRID-MSA TEST: waiting for the overnight queue ==="
while pgrep -f 'overnight_queue.sh' > /dev/null; do sleep 60; done
say "overnight queue finished."
while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
say "GPU free."

mkdir -p "$Q/hybrid_rbd_folds"
say "START hybrid_rbd  (60 folds, 1 model / 3 recycles, precomputed target MSA)"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$Q/hybrid_rbd" "$Q/hybrid_rbd_folds" > "$Q/hybrid_rbd.log" 2>&1 \
  && say "DONE  hybrid_rbd ($(ls "$Q/hybrid_rbd_folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l) folds)" \
  || say "FAIL  hybrid_rbd (partial results in $Q/hybrid_rbd_folds)"

say "=== HYBRID-MSA TEST COMPLETE ==="
} >> "$LOG" 2>&1
