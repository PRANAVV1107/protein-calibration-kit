#!/bin/bash
# Run 2: identical to the pre-registered search except nodes are valued by MEAN.
#
# Run 1 valued nodes by max-of-4 and lost to blind sampling (0.860 vs 0.870) at
# 17x the compute. The noise PoC then showed why: max is a tail statistic, and
# on RBD under this scoring config (skew +3.38) the error in a node's value
# GROWS with the number of draws, 0.120 at K=2 to 0.214 at K=8, while mean-of-K
# falls 0.071 to 0.044.
#
# Everything else is held fixed -- same root, same budget, same children, same
# sequences per node, same partial_T, same metric, same hardware -- so the only
# difference is the aggregation statistic. The best single sequence is still
# tracked and reported, because that is the design one would actually
# synthesise; only the value used to STEER the search changes.
#
# Pre-registered comparison, unchanged from run 1:
#   blind sampling, 60 folds   0.870   <- the bar
#   run 1, max-valued          0.860 in 1012 folds
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
cd "$KIT" || exit 1
unset PYTHONPATH
n=$(ps -eo args --width 300 | grep -c '[m]cts_design\.py')
if [ "$n" -gt 0 ]; then echo "REFUSING: $n mcts_design.py already running"; exit 1; fi
nohup setsid /root/miniconda3/envs/colabfold/bin/python scripts/mcts_design.py \
  --root-pdb "$ROOT/real_test_results/rbd_200_batch/batch_0/rfdiffusion/binder_47.pdb" \
  --target rbd --out-dir "$KIT/results/rl/mcts_run2_mean" \
  --budget 1000 --children 4 --seqs 4 --partial-t 10 --max-depth 6 \
  --metric ipTM --node-value mean --chunk 8 \
  > "$KIT/results/mcts_run2.log" 2>&1 </dev/null &
sleep 8
echo "launched:"; ps -eo pid,etime,args --width 200 | grep '[m]cts_design' | cut -c1-92
tail -3 "$KIT/results/mcts_run2.log"
