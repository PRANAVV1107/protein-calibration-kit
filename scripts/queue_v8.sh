#!/bin/bash
# Priority order, with the OOM fixed.
#
# The Adaptyv validation was killed at 166/365 by the kernel OOM killer:
#
#   Out of memory: Killed process (colabfold_batch) anon-rss:7034832kB
#   Mem: 7.6Gi total
#
# Cause: one colabfold_batch invocation over 365 inputs accumulates memory
# across queries, and XLA_PYTHON_CLIENT_MEM_FRACTION=4.0 lets JAX oversubscribe
# GPU memory into system RAM. The setting that rescues a 6GB GPU is what
# exhausted a 7.6GB host.
#
# Fix: process in chunks, one colabfold_batch call per chunk, so memory is
# released between them. Completed folds are skipped on re-entry, so this also
# resumes cleanly from the 166 already done.
#
# Order, highest impact first:
#   1  adaptyv validation  199 remaining  the ONLY measured outcomes here
#   2  il7ra rescreen       61 remaining  resumes from 53
#   3  null panel          420            constructed negatives, demoted since
#                                         325 measured ones now exist
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v8.log
CHUNK=${CHUNK:-40}

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
# lowered from 4.0: less oversubscription of GPU memory into system RAM, which
# is the resource that actually ran out.
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

# fold a directory of a3m inputs in chunks, skipping anything already scored
chunked(){                      # name, in_dir, out_dir
  local name="$1" src="$2" dst="$3"
  mkdir -p "$dst" "$Q/.chunk_$name"
  local todo=()
  for f in "$src"/*.a3m; do
    local b; b=$(basename "$f" .a3m)
    if ! ls "$dst/${b}_scores_rank_001_"*.json >/dev/null 2>&1; then todo+=("$f"); fi
  done
  say "START $name: ${#todo[@]} remaining of $(ls "$src"/*.a3m | wc -l), chunks of $CHUNK"
  local i=0 n=0
  while [ $i -lt ${#todo[@]} ]; do
    rm -rf "$Q/.chunk_$name"; mkdir -p "$Q/.chunk_$name"
    local j=0
    while [ $j -lt "$CHUNK" ] && [ $i -lt ${#todo[@]} ]; do
      cp "${todo[$i]}" "$Q/.chunk_$name/"; i=$((i+1)); j=$((j+1))
    done
    n=$((n+1))
    colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
      "$Q/.chunk_$name" "$dst" >> "$Q/$name.log" 2>&1 \
      || say "  chunk $n exited nonzero (continuing)"
    say "  chunk $n done -- $(count "$dst") total"
  done
  rm -rf "$Q/.chunk_$name"
  say "DONE  $name ($(count "$dst"))"
}

{
say "=== QUEUE V8 START (OOM fixed: chunked, MEM_FRACTION 4.0 -> 2.0) ==="
if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
  say "waiting for in-flight fold..."
  while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 20; done
fi
say "GPU free."

# ---- 1. the only measured outcomes in this project ----
chunked adaptyv_val "$Q/adaptyv_val" "$Q/adaptyv_val_folds"

# ---- 2. the 38 backbones the broken fast screen discarded ----
say "START il7ra_rescreen (resumes from $(count "$Q/il7ra_rescreen_folds"))"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 3 --num-recycle 8 \
  "$Q/il7ra_rescreen.fasta" "$Q/il7ra_rescreen_folds" > "$Q/il7ra_rescreen.log" 2>&1 \
  && say "DONE  il7ra_rescreen ($(count "$Q/il7ra_rescreen_folds"))" \
  || say "FAIL  il7ra_rescreen ($(count "$Q/il7ra_rescreen_folds") done)"

# ---- 3. constructed negatives ----
chunked null_panel "$Q/null_panel" "$Q/null_panel_folds"

say "=== QUEUE V8 COMPLETE ==="
} >> "$LOG" 2>&1
