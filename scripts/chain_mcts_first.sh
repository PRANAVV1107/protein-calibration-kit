#!/bin/bash
# MD finishes, then MCTS, then the remaining paper folds.
#
# Reordered at the user's request so the search runs tonight. MD is left to
# finish rather than killed -- it is ~1 h from done and is the stability data
# the paper needs; killing it would discard that hour for nothing.
#
# MCTS sizing, from measured primitives rather than estimates:
#   partial diffusion 0.5 min/design (n=12), fold 63 s/design (n=85)
#   1000 folds at 4 children x 4 seqs = ~62 expansions, ~250 backbones, ~21 h
#
# The search saves state every iteration and is fully resumable, so it can be
# stopped at any point and the tree so far is intact and analysable.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
LOG=$KIT/results/chain_mcts.log
cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }
md_busy(){ ps -eo args --width 400 2>/dev/null | grep -c '[m]d_stability\.py'; }

{
say "=== CHAIN: MD -> MCTS -> benchmark -> expert_r1 ==="
say "waiting for MD (running: $(md_busy))"
while [ "$(md_busy)" -gt 0 ]; do sleep 120; done
say "MD done: $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l)/37 systems"

say "--- MCTS (budget 1000 folds, ~21 h) ---"
export PYTHONPATH=/root/rf_deps:${PYTHONPATH:-}
/root/miniconda3/envs/colabfold/bin/python scripts/mcts_design.py \
  --root-pdb "$ROOT/real_test_results/rbd_200_batch/batch_0/rfdiffusion/binder_15.pdb" \
  --target rbd --out-dir results/rl/mcts_rbd \
  --budget 1000 --children 4 --seqs 4 --partial-t 10 --max-depth 6 \
  --metric ipTM --chunk 16
say "MCTS returned"
tail -15 results/rl/mcts_rbd/RESULT.txt 2>/dev/null

say "--- remaining paper folds ---"
bash scripts/resume_bench.sh
bash scripts/expert_iter_round1.sh
say "=== CHAIN COMPLETE ==="
} >> "$LOG" 2>&1
