#!/bin/bash
# The two jobs displaced by the search. Serial, chunked, resumable.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/chain_remaining.log
cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }
{
say "=== CHAIN: reward benchmark -> expert iteration r1 ==="
say "--- 1. reward benchmark (80 remaining of 165) ---"
bash scripts/resume_bench.sh
say "benchmark now $(ls results/queue/reward_bench_m1r3/*_scores_rank_001_*.json 2>/dev/null | wc -l)/165"
say "--- 2. expert iteration round 1 (120 folds) ---"
bash scripts/expert_iter_round1.sh
say "=== CHAIN COMPLETE ==="
} >> "$LOG" 2>&1
