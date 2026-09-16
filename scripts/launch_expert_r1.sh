KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/expert_iter_round1.sh'; then echo "already queued"; exit 0; fi
nohup setsid bash scripts/expert_iter_round1.sh >/dev/null 2>&1 </dev/null &
sleep 5
tail -3 results/expert_r1.log
echo "--- queue order now ---"
ps -eo etime,args --width 400 | grep -E 'md_stability|resume_bench|expert_iter_round1|colabfold_batch' | grep -v grep | cut -c1-88
