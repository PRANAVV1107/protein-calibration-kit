#!/bin/bash
# Null panel: 420 designed binders folded against targets they were NOT designed
# for. Every pairing is a negative by construction -- unrelated proteins do not
# bind by default.
#
# Purpose: give the calibration curve real negatives. It is currently fit on
# labels synthesised from ipTM thresholds, which is circular. A null panel
# provides an empirical distribution for each metric on non-binding pairs, so a
# candidate's score can be expressed as a percentile against this pipeline's own
# data rather than against a threshold borrowed from another dataset.
#
# Folded through the hybrid-MSA route (~50 s each). Characterising a null needs
# volume, not precision, and the hybrid keeps the target folded (pLDDT 79 vs 27
# for the naive fast screen) while making no MSA-server calls.
#
# Separate file, waits for queue_v4: bash reads a running script by byte offset,
# so editing one mid-execution can make it execute garbage.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v5.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }

{
say "=== QUEUE V5 (null panel) -- waiting for queue_v4 ==="
while pgrep -f 'queue_v4.sh' > /dev/null; do sleep 120; done
say "queue_v4 finished."
while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
say "GPU free."

mkdir -p "$Q/null_panel_folds"
say "START null_panel (420 cross-pairings, hybrid MSA route)"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$Q/null_panel" "$Q/null_panel_folds" > "$Q/null_panel.log" 2>&1 \
  && say "DONE  null_panel ($(ls "$Q/null_panel_folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l))" \
  || say "FAIL  null_panel ($(ls "$Q/null_panel_folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l) done; see log)"

say "=== QUEUE V5 COMPLETE ==="
} >> "$LOG" 2>&1
