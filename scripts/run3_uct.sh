#!/bin/bash
# Run 3: mean-valued nodes + corrected UCT selection, RBD, budget-matched.
# Pre-registered in reports/RUN3_PREREGISTRATION.md before launch.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/real_test_results/rbd_200_batch/batch_0/rfdiffusion/binder_47.pdb
OUT=$KIT/results/rl/mcts_run3_uct
LOG=$KIT/results/mcts_run3.log

# /root/rf_deps is pandas for the rfdiffusion env's Python 3.9 and breaks
# ColabFold's 3.10. Cleared here, not trusted to the caller -- this runs
# unattended and that contamination has killed jobs twice.
unset PYTHONPATH

# The benchmark's proven JAX settings. A speculative change to MEM_FRACTION
# once stalled every fold on this machine and was misdiagnosed three times.
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

# Refuse to start alongside an existing driver. Matches the interpreter path AND
# the script path rather than a bare -f fragment, which has self-matched and
# deadlocked this queue three times.
n=$(ps -eo args --width 400 | grep -cE 'envs/colabfold/bin/python .*scripts/mcts_design\.py')
if [ "$n" -gt 0 ]; then
  echo "[abort] an mcts_design.py driver is already running"; exit 1
fi

cd "$KIT" || exit 1
mkdir -p "$OUT"
echo "[$(date '+%F %T')] === RUN 3: mean value + corrected UCT ===" >> "$LOG"
exec /root/miniconda3/envs/colabfold/bin/python scripts/mcts_design.py \
  --root-pdb "$ROOT" --target rbd --out-dir "$OUT" \
  --budget 1000 --children 4 --seqs 4 --partial-t 10 --max-depth 6 \
  --metric ipTM --node-value mean --chunk 8 >> "$LOG" 2>&1
