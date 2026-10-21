#!/bin/bash
# Launch the pre-registered search. ONE instance, guarded.
#
# Two concurrent mcts_design.py processes were left running earlier, writing to
# the same output directory and each spawning its own ColabFold -- they halved
# the GPU and would have corrupted state.json. That happened because a timeout
# on the wrapper killed the outer shell while the nohup'd search kept running,
# and it was then relaunched on top. The guard below makes that impossible.
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
cd "$KIT" || exit 1

n=$(ps -eo args --width 300 | grep -c '[m]cts_design\.py')
if [ "$n" -gt 0 ]; then
  echo "REFUSING: $n mcts_design.py already running"
  ps -eo pid,etime,args --width 200 | grep '[m]cts_design' | cut -c1-90
  exit 1
fi

nohup setsid /root/miniconda3/envs/colabfold/bin/python scripts/mcts_design.py \
  --root-pdb "$ROOT/real_test_results/rbd_200_batch/batch_0/rfdiffusion/binder_47.pdb" \
  --target rbd --out-dir "$KIT/results/rl/mcts_run1" \
  --budget 1000 --children 4 --seqs 4 --partial-t 10 --max-depth 6 \
  --metric ipTM --chunk 8 \
  > "$KIT/results/mcts_run1.log" 2>&1 </dev/null &
sleep 8
echo "launched. pids:"
ps -eo pid,etime,args --width 200 | grep '[m]cts_design' | cut -c1-90
tail -3 "$KIT/results/mcts_run1.log"
