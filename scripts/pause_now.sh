#!/bin/bash
# Release the GPU for the user. Covers the currently-running process set.
#
# Safe: ColabFold writes a .done.txt per query and skips those on restart, so
# only the single in-flight fold is lost. MD is waiting, not running, so it
# loses nothing. Schedulers go before their children, otherwise killing a fold
# just makes the queue advance and take the card straight back.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
c(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }
echo "before pause: v2_pdl1 $(c results/queue/v2refold_folds_pdl1)/285  v2_rbd $(c results/queue/v2refold_folds_rbd)/278"

for s in queue_v11.sh queue_md_gpu3.sh queue_md_gpu2.sh watch_highacc.sh; do
  pkill -f "bash scripts/$s" && echo "stopped $s" || true
done
sleep 2
pkill -f 'envs/colabfold/bin/colabfold_batch' && echo "stopped in-flight fold" || echo "no fold running"
pkill -f 'scripts/md_stability\.py' && echo "stopped MD" || true
sleep 5
echo
echo "GPU:  $(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader)"
free -g | sed -n 2p
echo
if ps -eo args --width 300 | grep -q '[b]ash scripts/sched_0200.sh'; then
  echo "2am scheduler: STILL ARMED (fires 02:00, starts the reward benchmark)"
else
  echo "2am scheduler: not running"
fi
