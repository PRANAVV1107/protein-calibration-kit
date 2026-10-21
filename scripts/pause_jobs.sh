#!/bin/bash
# Pause the overnight jobs. Default is graceful: the runner checks the flag
# between chunks, finishes the chunk in flight, and exits cleanly with nothing
# lost. Pass --now to also kill the fold in progress, which discards only that
# chunk (it is redone on resume).
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
PAUSE=$KIT/results/.PAUSE
CHUNK_EST=5

touch "$PAUSE"
echo "[pause] flag set: $PAUSE"

pids_of(){ ps -eo pid,args --width 500 | grep -E "$1" | awk '{print $1}'; }

if [ "${1:-}" = "--now" ]; then
  echo "[pause] --now: killing the fold in progress"
  for p in $(pids_of '[b]ash .*scripts/overnight\.sh'); do kill "$p" 2>/dev/null && echo "  TERM runner $p"; done
  sleep 3
  for p in $(pids_of '[c]olabfold_batch'); do kill "$p" 2>/dev/null && echo "  TERM colabfold $p"; done
  sleep 6
  for p in $(pids_of '[c]olabfold_batch|[b]ash .*scripts/overnight\.sh'); do
    kill -9 "$p" 2>/dev/null && echo "  KILL $p"; done
  sleep 2
else
  echo "[pause] graceful: the runner will stop after the current chunk"
  echo "        (up to ~$CHUNK_EST designs, roughly 20 min). Use --now to stop immediately."
fi

echo
echo "[pause] state on disk:"
for d in rog_prospective calibration; do
  n=$(find "$KIT/results/$d/folds" -name '*_scores_rank_001_*.json' 2>/dev/null | wc -l)
  t=$(grep -c '^>' "$KIT/results/$d/to_fold.fasta" 2>/dev/null || echo 0)
  echo "  $d: $n / $t folded"
done
echo
echo "[pause] still running: $(ps -eo args --width 400 | grep -cE '[c]olabfold_batch|[b]ash .*scripts/overnight\.sh')"
echo "[pause] resume with: bash scripts/resume_jobs.sh"
