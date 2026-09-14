#!/bin/bash
# Final ordering. The v2 refold changes the structure, so this supersedes the
# tail of queue_v8 rather than appending to it.
#
# Why the rescreen is dropped: it folds 38 IL-7Ra backbones x 3 draws at high
# accuracy. The v2 refold covers ALL backbones x ALL draws through the hybrid
# route at a ninth of the cost per fold. Screening everything cheaply and then
# spending high accuracy only on the winners is strictly better than spending
# high accuracy on a subset chosen in advance -- which is the same mistake the
# fast screen made, just with a better metric.
#
#   1  adaptyv validation   (running in queue_v8)  decides WHICH metric to rank on
#   2  v2 refold      881    the complete design space, corrected method
#   3  v2 finalists   ~60    high accuracy on the top of that screen
#   4  GPU MD         15     50 ns, orthogonal to any structure predictor
#   5  null panel     420    constructed negatives, demoted since 325 measured
#                            ones now exist
#
# Steps 2 and 5 make no MSA-server calls; every alignment is precomputed.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v9.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

CHUNK=${CHUNK:-40}
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

# Fold a directory of a3m inputs in chunks, skipping anything already scored.
# Chunking is not cosmetic: one call over hundreds of inputs accumulates memory
# across queries and was OOM-killed at 7.0GB on this 7.6GB host.
chunked(){                      # name, in_dir, out_dir
  local name="$1" src="$2" dst="$3"
  mkdir -p "$dst"
  local todo=()
  for f in "$src"/*.a3m; do
    local b; b=$(basename "$f" .a3m)
    ls "$dst/${b}_scores_rank_001_"*.json >/dev/null 2>&1 || todo+=("$f")
  done
  say "START $name: ${#todo[@]} remaining, chunks of $CHUNK"
  local i=0 n=0
  while [ $i -lt ${#todo[@]} ]; do
    rm -rf "$Q/.c_$name"; mkdir -p "$Q/.c_$name"
    local j=0
    while [ $j -lt "$CHUNK" ] && [ $i -lt ${#todo[@]} ]; do
      cp "${todo[$i]}" "$Q/.c_$name/"; i=$((i+1)); j=$((j+1))
    done
    n=$((n+1))
    colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
      "$Q/.c_$name" "$dst" >> "$Q/$name.log" 2>&1 || say "  chunk $n nonzero (continuing)"
    say "  chunk $n done -- $(count "$dst") total"
  done
  rm -rf "$Q/.c_$name"
  say "DONE  $name ($(count "$dst"))"
}

{
say "=== QUEUE V9 START ==="
while pgrep -f 'queue_v8.sh' > /dev/null; do sleep 60; done
say "queue_v8 finished."
# interpreter PLUS arguments: a bare path fragment matches any process that
# merely mentions it, which has deadlocked this project twice.
while pgrep -f 'colabfold_batch --model' > /dev/null; do sleep 30; done
say "GPU free."

# ---- 2. the complete design space ----
for t in il7ra pdl1 rbd; do
  [ -d "$Q/v2refold/$t" ] || continue
  chunked "v2_$t" "$Q/v2refold/$t" "$Q/v2refold_folds_$t"
done

# ---- 3. high accuracy on the winners of that screen ----
cd "$KIT"
say "selecting v2 finalists"
python scripts/select_v2_finalists.py >> "$LOG" 2>&1 || say "FAIL select_v2_finalists"
if [ -f "$Q/v2_finalists.fasta" ]; then
  mkdir -p "$Q/v2_finalists_folds"
  say "START v2_finalists (high accuracy)"
  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
    --num-models 3 --num-recycle 8 \
    "$Q/v2_finalists.fasta" "$Q/v2_finalists_folds" > "$Q/v2_finalists.log" 2>&1 \
    && say "DONE  v2_finalists ($(count "$Q/v2_finalists_folds"))" || say "FAIL  v2_finalists"
fi

# ---- 5. constructed negatives, lowest priority ----
chunked null_panel "$Q/null_panel" "$Q/null_panel_folds"

say "=== QUEUE V9 COMPLETE ==="
} >> "$LOG" 2>&1
