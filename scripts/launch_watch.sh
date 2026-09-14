#!/bin/bash
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT"
if ps -eo args --width 300 | grep -q 'bash scripts/watch_highacc\.sh'; then
  echo "already running"; exit 0
fi
bash -n scripts/watch_highacc.sh && echo "[OK] syntax"
nohup setsid bash scripts/watch_highacc.sh > /dev/null 2>&1 < /dev/null &
sleep 5
cat results/watch_highacc.log
