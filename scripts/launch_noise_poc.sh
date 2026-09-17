KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/noise_poc.sh'; then echo "already queued"; exit 0; fi
nohup setsid bash scripts/noise_poc.sh >/dev/null 2>&1 </dev/null &
sleep 5
tail -3 results/noise_poc.log
echo "--- queue ---"
ps -eo etime,args --width 300 | grep -E 'chain_remaining|resume_bench|expert_iter|noise_poc|colabfold_batch' | grep -v grep | cut -c1-78
