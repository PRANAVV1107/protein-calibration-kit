KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/queue_mcts_pilot.sh'; then echo "already queued"; exit 0; fi
nohup setsid bash scripts/queue_mcts_pilot.sh >/dev/null 2>&1 </dev/null &
sleep 4
tail -3 results/mcts_pilot.log
echo "--- queue ---"
ps -eo etime,args --width 400 | grep -E 'md_stability|resume_bench|expert_iter|queue_mcts|colabfold_batch' | grep -v grep | cut -c1-80
