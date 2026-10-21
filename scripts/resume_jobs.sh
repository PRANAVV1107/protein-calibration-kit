#!/bin/bash
# Clear the pause flag and restart the overnight jobs where they stopped.
# Safe to run at any time: the runner refuses to start if folding is already
# in progress, and rebuilds its to-do list from what is on disk, so nothing is
# ever folded twice.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
PAUSE=$KIT/results/.PAUSE

if [ -f "$PAUSE" ]; then rm -f "$PAUSE"; echo "[resume] pause flag cleared"; fi

busy=$(ps -eo args --width 400 | grep -cE '[c]olabfold_batch|[b]ash .*scripts/overnight\.sh')
if [ "$busy" -gt 0 ]; then
  echo "[resume] already running ($busy processes) -- nothing to do"
  exit 0
fi

cd "$KIT" || exit 1
echo "[resume] state before restart:"
for d in rog_prospective calibration; do
  n=$(find "results/$d/folds" -name '*_scores_rank_001_*.json' 2>/dev/null | wc -l)
  t=$(grep -c '^>' "results/$d/to_fold.fasta" 2>/dev/null || echo "not yet built")
  echo "  $d: $n / $t"
done

setsid nohup bash scripts/overnight.sh > /dev/null 2>&1 < /dev/null &
sleep 20
echo "[resume] runner alive: $(ps -eo args --width 400 | grep -cE '[b]ash .*scripts/overnight\.sh')"
echo "[resume] tail:"
tail -4 results/overnight.log 2>/dev/null
