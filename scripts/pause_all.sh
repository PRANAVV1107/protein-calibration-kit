#!/bin/bash
# Release the GPU. Stops every scheduler and the in-flight fold.
#
# Safe to run at any time: ColabFold writes a .done.txt per query as it
# completes, and skips those on restart, so the only work lost is the single
# fold in flight. Order matters -- the scheduler goes first, otherwise killing
# the fold just makes queue_v10 advance to the next job and take the card back.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT"

c(){ ls results/queue/adaptyv_highacc_folds/*_scores_rank_001_*.json 2>/dev/null | wc -l; }
echo "high-accuracy arm at $(c)/165 before pause"

pkill -f 'bash scripts/queue_v10\.sh'      && echo "stopped queue_v10"      || echo "queue_v10 not running"
pkill -f 'bash scripts/watch_highacc\.sh'  && echo "stopped watch_highacc"  || echo "watch_highacc not running"
pkill -f 'bash scripts/queue_md_gpu2\.sh'  && echo "stopped queue_md_gpu2"  || echo "queue_md_gpu2 not running"
sleep 2
pkill -f 'envs/colabfold/bin/colabfold_batch' && echo "stopped colabfold" || echo "colabfold not running"
pkill -f 'scripts/md_stability\.py'           && echo "stopped MD"        || true
sleep 5

echo
echo "GPU now:"
nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader
free -g | head -2
echo
echo "Resume with: bash scripts/resume_all.sh"
