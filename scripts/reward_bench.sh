#!/bin/bash
# Find the cheapest reward that still predicts measured binding.
#
# This is the MCTS critical path. Node cost is dominated by the reward, so the
# whole question of whether MCTS fits on a 6GB card is "how cheap can the
# reward get before it stops predicting". Two configurations are already
# measured on the 164 labelled Adaptyv designs:
#
#   hybrid a3m, 1 model, 3 recycle   ~2 min   ipTM 0.545, ipSAE 0.466  (chance)
#   full MSA,   3 models, 8 recycle  ~4 min   ipTM 0.614, ipSAE 0.618  (works)
#
# The untested middle is full MSA with the cheap model settings. If that keeps
# the predictive power at a third of the cost, MCTS goes from ~42 h to ~15 h
# for 500 nodes. If it collapses to the hybrid's level, then the MSA is not
# what mattered and the recycles are, which is equally worth knowing.
#
# MSAs are reused from the completed high-accuracy run, so no MMseqs2 search is
# repeated -- that was the dominant cost there and it is already paid.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/reward_bench.log
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform
cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }

{
say "=== REWARD BENCH START ==="

# 1. stop the chance-level screen and reclaim the GPU
for p in queue_v11.sh; do
  pkill -f "bash scripts/$p" && say "stopped $p" || say "$p not running"
done
sleep 3
pkill -f 'envs/colabfold/bin/colabfold_batch' && say "stopped in-flight fold" || true
sleep 5

# 2. cheap-model / full-MSA arm, reusing the a3m files already computed
IN=$Q/reward_bench_in; OUT=$Q/reward_bench_m1r3
mkdir -p "$IN" "$OUT"
n=0
for a in "$Q"/adaptyv_highacc_folds/*.a3m; do
  b=$(basename "$a"); cp "$a" "$IN/$b" 2>/dev/null && n=$((n+1))
done
say "reusing $n precomputed MSAs (no MMseqs2 re-search)"

t0=$(date +%s)
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$IN" "$OUT" >> "$Q/reward_bench.log" 2>&1 || say "colabfold returned nonzero"
t1=$(date +%s)
done_n=$(ls "$OUT"/*_scores_rank_001_*.json 2>/dev/null | wc -l)
say "m1r3 arm: $done_n folds in $(( (t1-t0)/60 )) min  ($(( (t1-t0) / (done_n>0?done_n:1) )) s/design)"

say "=== REWARD BENCH DONE -- analyse with compare_arms.py --high-dir $OUT ==="
} >> "$LOG" 2>&1
