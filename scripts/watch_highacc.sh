#!/bin/bash
# Run the controlled comparison the moment the high-accuracy arm completes.
#
# queue_v10.sh already runs analyze_adaptyv.py on that directory, but that
# script scores one arm in isolation. The comparison that answers the actual
# question -- is the cheap screen's null a failure of the screen or of the
# metrics -- is paired across both arms on identical designs, which is
# compare_arms.py. This runs it without waiting for me to be at the keyboard.
#
# Also writes an interim read every 25 folds, so a trend is visible before the
# arm finishes and the run can be cut short if it is clearly going nowhere.
#
# Guard convention: match the interpreter path followed by the script path.
# Bare `pgrep -f <fragment>` self-matches watcher shells and has deadlocked
# this project three times.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/watch_highacc.log
PY=/root/miniconda3/envs/colabfold/bin/python
TARGET=165

cd "$KIT"
say(){ echo "[$(date '+%F %T')] $*"; }
c(){ ls results/queue/adaptyv_highacc_folds/*_scores_rank_001_*.json 2>/dev/null | wc -l; }
folding(){
  ps -eo args --width 400 2>/dev/null \
    | grep -c 'envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch'
}

{
say "=== WATCH HIGHACC START (target $TARGET) ==="
last=0; stall=0; nextmark=0
while [ "$(c)" -lt "$TARGET" ]; do
  n=$(c)
  if [ "$n" -eq "$last" ]; then
    stall=$((stall+1))
  else
    stall=0; last=$n
  fi
  if [ "$n" -ge "$nextmark" ] && [ "$n" -ge 50 ]; then
    say "interim at $n/$TARGET"
    "$PY" scripts/compare_arms.py --boot 2000 \
      > "results/ARM_COMPARISON_interim_${n}.txt" 2>&1 \
      && say "  -> results/ARM_COMPARISON_interim_${n}.txt"
    nextmark=$((n+25))
  fi
  # the fold job died or finished early: 30 x 10 min = 5 h with no new fold
  if [ "$stall" -ge 30 ]; then
    say "no progress in ~5h at $n/$TARGET -- proceeding with what exists"
    break
  fi
  sleep 600
done

n=$(c)
say "high-accuracy arm at $n/$TARGET. Waiting for the fold process to exit."
while [ "$(folding)" -gt 0 ]; do sleep 60; done

say "Running final comparison."
"$PY" scripts/compare_arms.py --boot 10000 \
  > "$KIT/results/ARM_COMPARISON_FINAL.txt" 2>&1 \
  && say "FINAL -> results/ARM_COMPARISON_FINAL.txt" \
  || say "comparison FAILED (see that file)"
say "=== WATCH HIGHACC DONE ==="
} >> "$LOG" 2>&1
