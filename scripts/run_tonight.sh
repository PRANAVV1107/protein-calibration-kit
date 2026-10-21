#!/bin/bash
# Reward benchmark, then MD. Chained so tonight produces both.
#
# The benchmark must finish before MD starts: they would otherwise contend for
# 6GB of VRAM, and concurrent GPU work has OOM'd in this project before.
# queue_md_gpu3 has its own guard against ColabFold, but running them in
# sequence here makes the ordering explicit rather than relying on it.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/run_tonight.log
cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }
{
say "=== RUN TONIGHT START ==="
say "--- stage 1: reward benchmark ---"
bash scripts/reward_bench.sh
say "stage 1 returned; folds: $(ls results/queue/reward_bench_m1r3/*_scores_rank_001_*.json 2>/dev/null | wc -l)"
say "--- stage 2: MD panel (37 complexes, 2 ns) ---"
NS=2 EQUIL=200 bash scripts/queue_md_gpu3.sh
say "stage 2 returned; MD systems: $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l)"
say "=== RUN TONIGHT COMPLETE ==="
} >> "$LOG" 2>&1
