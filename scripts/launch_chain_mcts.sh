KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
for s in queue_mcts_pilot.sh resume_bench.sh expert_iter_round1.sh; do
  pkill -f "bash scripts/$s" && echo "unqueued $s (re-chained)" || true
done
sleep 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/chain_mcts_first.sh'; then echo "chain already running"; exit 0; fi
nohup setsid bash scripts/chain_mcts_first.sh >/dev/null 2>&1 </dev/null &
sleep 4
tail -3 results/chain_mcts.log
echo "--- now running ---"
ps -eo etime,args --width 400 | grep -E 'chain_mcts|md_stability|mcts_design' | grep -v grep | cut -c1-80
