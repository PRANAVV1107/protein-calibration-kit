#!/bin/bash
# Finish the reward benchmark after MD, in chunks this time.
#
# The first attempt ran one colabfold_batch call over all 165 inputs and was
# OOM-killed at 85 ("Killed", line 57). That is the third time this project has
# hit the same failure: a single call holds every query's features resident, and
# 7.6 GB of host RAM will not take 165 of them. queue_v11 used chunks of 40 over
# far more inputs and never OOM'd. So: chunks of 30 here, smaller than 40
# because MD may still be holding memory when this starts.
#
# Completed folds persist via .done.txt and are skipped, so this picks up the
# remaining ~80 rather than redoing 165.
#
# Waits for MD first. MD is paper-critical and more than half done; the
# benchmark only needs its AUC, which nothing downstream is blocked on tonight.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/resume_bench.log
IN=$Q/reward_bench_in
OUT=$Q/reward_bench_m1r3
CHUNK=${CHUNK:-30}

export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$OUT"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }
md_running(){ ps -eo args --width 400 2>/dev/null | grep -c '[m]d_stability\.py'; }

{
say "=== RESUME BENCH START (at $(count)/165) ==="
say "waiting for MD to finish (running: $(md_running))"
while [ "$(md_running)" -gt 0 ]; do sleep 120; done
say "MD done. $(ls $KIT/results/md_gpu/*_md.json 2>/dev/null | wc -l)/37 systems."

todo=()
for f in "$IN"/*.a3m; do
  b=$(basename "$f" .a3m)
  ls "$OUT/${b}_scores_rank_001_"*.json >/dev/null 2>&1 || todo+=("$f")
done
say "${#todo[@]} remaining, chunks of $CHUNK"
[ ${#todo[@]} -eq 0 ] && { say "already complete"; exit 0; }

i=0; n=0
while [ $i -lt ${#todo[@]} ]; do
  rm -rf "$Q/.rb_chunk"; mkdir -p "$Q/.rb_chunk"
  j=0
  while [ $j -lt "$CHUNK" ] && [ $i -lt ${#todo[@]} ]; do
    cp "${todo[$i]}" "$Q/.rb_chunk/"; i=$((i+1)); j=$((j+1))
  done
  n=$((n+1))
  colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
    "$Q/.rb_chunk" "$OUT" >> "$Q/reward_bench.log" 2>&1 || say "  chunk $n nonzero (continuing)"
  say "  chunk $n done -- $(count)/165"
done
rm -rf "$Q/.rb_chunk"
say "=== RESUME BENCH COMPLETE ($(count)/165) ==="
} >> "$LOG" 2>&1
