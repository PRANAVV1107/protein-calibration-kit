#!/bin/bash
# Sleep until 02:00 on 2026-09-16, then run the reward benchmark.
#
# A plain sleep rather than cron/at: neither is reliably running in this WSL
# distro, and this survives as an ordinary background process. It does NOT
# survive a reboot -- Windows Update is paused to Oct 1, which was the cause
# last time, but if the machine goes down this needs relaunching.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/sched_0200.log
TARGET=$(date -d '2026-09-16 02:00:00' +%s 2>/dev/null)
{
now=$(date +%s)
wait_s=$(( TARGET - now ))
echo "[$(date '+%F %T')] scheduled for $(date -d @$TARGET '+%F %T') (in $(( wait_s/60 )) min)"
if [ "$wait_s" -gt 0 ]; then sleep "$wait_s"; fi
echo "[$(date '+%F %T')] firing reward_bench.sh"
bash "$KIT/scripts/reward_bench.sh"
echo "[$(date '+%F %T')] reward_bench returned"
} >> "$LOG" 2>&1
