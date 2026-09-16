KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/resume_bench.sh'; then echo "already queued"; exit 0; fi
nohup setsid bash scripts/resume_bench.sh >/dev/null 2>&1 </dev/null &
sleep 5
tail -3 results/resume_bench.log
