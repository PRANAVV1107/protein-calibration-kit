#!/bin/bash
# Re-run expert iteration round 1 with valid a3m inputs, after the noise PoC.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/requeue_expert_r1.log
cd "$KIT" || exit 1
unset PYTHONPATH
say(){ echo "[$(date '+%F %T')] $*"; }
busy(){ ps -eo args --width 400 2>/dev/null \
        | grep -cE 'envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch|[p]rotein_mpnn_run|[b]ash scripts/noise_poc\.sh'; }
{
say "=== REQUEUE expert_r1 (valid a3m) ==="
say "waiting for the noise PoC (busy: $(busy))"
while [ "$(busy)" -gt 0 ]; do sleep 120; done
say "GPU free."
bash scripts/expert_iter_round1.sh
say "returned. folds: $(ls results/rl/expert_r1/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l)/120"
cat results/EXPERT_R1_RESULT.txt 2>/dev/null
say "=== DONE ==="
} >> "$LOG" 2>&1
