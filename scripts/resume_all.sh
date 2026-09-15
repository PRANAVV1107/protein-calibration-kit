#!/bin/bash
# Restart the pipeline after pause_all.sh.
#
# Does NOT re-run queue_v10.sh. That script's first act is to fold the cheap
# adaptyv arm, which is already complete, and its chunked helper would then
# walk the v2 refold -- correct, but it would restart the high-accuracy arm
# from its own beginning. ColabFold skips finished queries via .done.txt, so
# resuming the high-accuracy fold directly is both faster and easier to read.
#
# Order is deliberate: finish the high-accuracy arm (the answer everything
# else is waiting on), then the v2 refold, then the null panel. MD stays
# disarmed until folding is done, as before.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v11.log
PY=/root/miniconda3/envs/colabfold/bin/python

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

cd "$KIT"

if ps -eo args --width 300 | grep -q 'bash scripts/queue_v11\.sh'; then
  echo "already running"; exit 0
fi

bash -n scripts/queue_v11.sh || { echo "syntax error in queue_v11.sh"; exit 1; }
nohup setsid bash scripts/queue_v11.sh > /dev/null 2>&1 < /dev/null &
nohup setsid bash scripts/watch_highacc.sh > /dev/null 2>&1 < /dev/null &
sleep 6
echo "restarted. logs:"
tail -3 "$LOG" 2>/dev/null
tail -2 "$KIT/results/watch_highacc.log" 2>/dev/null
