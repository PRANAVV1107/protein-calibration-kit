#!/bin/bash
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
pkill -f 'bash scripts/sched_0200.sh' && echo "disarmed 2am scheduler" || echo "2am scheduler not running"
sleep 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/run_tonight.sh'; then
  echo "run_tonight already active"; exit 0
fi
nohup setsid bash scripts/run_tonight.sh >/dev/null 2>&1 </dev/null &
sleep 25
echo "--- run_tonight.log ---"
tail -6 results/run_tonight.log 2>/dev/null || echo "(none yet)"
echo "--- processes ---"
ps -eo etime,args --width 400 | grep -E 'run_tonight|colabfold_batch' | grep -v grep | cut -c1-90
