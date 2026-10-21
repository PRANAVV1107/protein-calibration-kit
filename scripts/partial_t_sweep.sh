#!/bin/bash
# Does the search space have any locality at all?
#
# Run 3 measured parent->best-child Spearman rho = 0.168 (p = 0.384) at
# partial_T = 10, with every record-setting design descending from a
# bottom-decile parent. If expanding a node gives an uncorrelated draw, the tree
# has no edges worth following and best-first search IS random search.
#
# But partial_T = 10 was never chosen on evidence. This sweeps it. The question
# is whether the score landscape has a correlation length at all, and if so
# where it is.
#
#   rho rises sharply at low partial_T -> locality exists, the hyperparameter was
#                                          simply wrong, and every negative search
#                                          result here is conditional on it
#   rho null at every partial_T         -> tree search over RFdiffusion backbones
#                                          is structurally dead for this pipeline,
#                                          with a measured mechanism
#
# Both outcomes are publishable; the second is the stronger paper.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
ROOTDIR=/mnt/c/Users/prana/Downloads/proteinfoldingexp/real_test_results/rbd_200_batch/batch_0/rfdiffusion
OUT=$KIT/results/rl/partial_t
LOG=$KIT/results/partial_t.log
TS=${TS:-"2 4 10 20"}          # partial_T values to compare
PARENTS=${PARENTS:-6}          # distinct parent backbones
KIDS=${KIDS:-4}                # children per parent per T
SEQS=${SEQS:-4}                # MPNN draws per node

unset PYTHONPATH
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

say(){ echo "[$(date '+%F %T')] $*"; }
n=$(ps -eo args --width 400 | grep -cE '[e]nvs/colabfold/bin/python .*scripts/mcts_design\.py|[e]nvs/colabfold/bin/colabfold_batch|[m]d_stability\.py')
if [ "$n" -gt 0 ]; then echo "[abort] GPU busy ($n processes)"; exit 1; fi

cd "$KIT" || exit 1
mkdir -p "$OUT"
{
say "=== PARTIAL_T SWEEP ==="
say "T values: $TS | $PARENTS parents x $KIDS children x $SEQS seqs"
say "budget: $(( $(echo $TS | wc -w) * PARENTS * (KIDS + 1) * SEQS )) folds"

# Parents are drawn WITHOUT reference to score -- selecting good parents would
# bias the correlation by restricting range, which is the classic way to make a
# real signal vanish.
for T in $TS; do
  for i in $(seq 0 $((PARENTS-1))); do
    P=$(printf "binder_%d" $((i*7+3)))     # spread across the batch, fixed rule
    TAG="t${T}_p${i}"
    [ -d "$OUT/$TAG" ] && { say "skip $TAG (exists)"; continue; }
    say "--- $TAG: parent $P at partial_T=$T ---"
    /root/miniconda3/envs/colabfold/bin/python scripts/mcts_design.py \
      --root-pdb "$ROOTDIR/$P.pdb" --target rbd \
      --out-dir "$OUT/$TAG" \
      --budget $(( (KIDS + 1) * SEQS )) \
      --children "$KIDS" --seqs "$SEQS" --partial-t "$T" \
      --max-depth 1 --metric ipTM --node-value mean --chunk 8 \
      || say "  $TAG returned nonzero"
  done
done

say "folds on disk: $(find "$OUT" -name '*_scores_rank_001_*.json' | wc -l)"
/root/miniconda3/envs/colabfold/bin/python scripts/partial_t_analyze.py --dir "$OUT" \
  > "$KIT/results/PARTIAL_T_RESULT.txt" 2>&1 \
  && say "analysed -> results/PARTIAL_T_RESULT.txt" || say "analysis FAILED"
say "=== SWEEP COMPLETE ==="
} >> "$LOG" 2>&1
