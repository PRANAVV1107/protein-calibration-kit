#!/bin/bash
# 52-fold MCTS pilot from binder_15, queued behind everything else.
#
# What 52 folds CAN establish:
#   - the loop works end to end (RFdiffusion -> MPNN -> fold -> backprop -> resume)
#   - true cost per node, which is what sizes any cloud run
#   - whether partial_T neighbours beat their parent, the premise of the method
#
# What it CANNOT establish: that search beats best-of-N. A tree of ~13 nodes is
# far too small for that comparison and it should not be reported as one.
#
# Root is binder_15, chosen because its baseline is already paid for: the
# variance study folded 6 sequences on it at high accuracy, mean ipTM 0.333,
# best 0.71. Two free reference points follow --
#   root best-of-6                            0.71
#   60 folds blind across 10 RBD backbones    0.87   <- the real bar, equal budget
# A middling backbone was picked deliberately over a known-good one, so there is
# headroom for improvement to be visible.
#
# Metric is ipTM so it compares directly against those numbers.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
LOG=$KIT/results/mcts_pilot.log
cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }
busy(){ ps -eo args --width 400 2>/dev/null \
        | grep -cE '[m]d_stability\.py|envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch|[p]rotein_mpnn_run'; }
{
say "=== MCTS PILOT queued ==="
say "waiting for MD / benchmark / expert-r1 (busy: $(busy))"
while [ "$(busy)" -gt 0 ]; do sleep 120; done
say "GPU free. starting."
export PYTHONPATH=/root/rf_deps:${PYTHONPATH:-}
/root/miniconda3/envs/colabfold/bin/python scripts/mcts_design.py \
  --root-pdb "$ROOT/real_test_results/rbd_200_batch/batch_0/rfdiffusion/binder_15.pdb" \
  --target rbd --out-dir results/rl/mcts_pilot \
  --budget 52 --children 3 --seqs 4 --partial-t 10 --max-depth 3 \
  --metric ipTM --chunk 12
say "=== MCTS PILOT returned ==="
tail -20 results/rl/mcts_pilot/RESULT.txt 2>/dev/null
} >> "$LOG" 2>&1
