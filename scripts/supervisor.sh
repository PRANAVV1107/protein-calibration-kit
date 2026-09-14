#!/bin/bash
# Unattended supervisor for the overnight run.
#
# Fixes a scheduling problem: queue_v8 still holds the rescreen and the null
# panel after the Adaptyv validation, and queue_v9 waits for queue_v8 to exit.
# Left alone, that spends ~10 h on work that was deliberately demoted before
# starting the v2 refold, which is the higher-value job. So once Adaptyv is
# complete this stops queue_v8, which lets queue_v9 take over with the intended
# ordering.
#
# It also runs the Adaptyv analysis the moment the data is complete, so the
# result is waiting rather than the raw folds.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/supervisor.log
PY=/root/miniconda3/envs/colabfold/bin/python
TARGET=365

say(){ echo "[$(date '+%F %T')] $*"; }
adaptyv_count(){ ls "$KIT/results/queue/adaptyv_val_folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

{
say "=== SUPERVISOR START (waiting for Adaptyv to reach $TARGET) ==="
last=0
stall=0
while [ "$(adaptyv_count)" -lt "$TARGET" ]; do
  n=$(adaptyv_count)
  if [ "$n" -eq "$last" ]; then
    stall=$((stall+1))
  else
    stall=0; last=$n
    say "  adaptyv $n/$TARGET"
  fi
  # 40 x 5 min = ~3.3 h with no new fold means the job died; stop waiting
  if [ "$stall" -ge 40 ]; then
    say "  no progress in ~3h at $n/$TARGET -- proceeding with what exists"
    break
  fi
  sleep 300
done

n=$(adaptyv_count)
say "Adaptyv at $n/$TARGET. Running analysis."
cd "$KIT"
"$PY" scripts/analyze_adaptyv.py > "$KIT/results/ADAPTYV_RESULT.txt" 2>&1 \
  && say "analysis written to results/ADAPTYV_RESULT.txt" \
  || say "analysis FAILED (see that file)"

# stop queue_v8 so its demoted tail does not block the v2 refold
if pgrep -f 'queue_v8.sh' > /dev/null; then
  say "stopping queue_v8 so queue_v9 can start the v2 refold"
  pkill -f 'queue_v8.sh' || true
  sleep 5
  # the fold it launched is a child; stop that too, its results are on disk
  pkill -f 'colabfold_batch --model' || true
  say "queue_v8 stopped; queue_v9 should pick up within ~60s"
fi

say "=== SUPERVISOR DONE ==="
} >> "$LOG" 2>&1
