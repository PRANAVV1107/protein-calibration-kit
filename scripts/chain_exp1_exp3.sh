#!/bin/bash
# Chain: wait for Experiment 1 (rescore) to finish, then start Experiment 3
# (partial_T sweep) without waiting for a human.
#
# Why a chain rather than a cron: partial_t_sweep.sh refuses to start while
# anything else holds the GPU, so a fixed-time trigger would either fire early
# and abort, or fire late and waste the card. This waits on the actual
# completion signal.
#
# Completion is judged by BOTH the fold count and the absence of a folding
# process. Either alone is unsafe: the count can reach target while the last
# chunk is still writing, and the process can be briefly absent between chunks.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/chain_exp1_exp3.log
POLL=${POLL:-120}
MAX_WAIT=${MAX_WAIT:-86400}        # 24 h. The 4 h estimate was wrong: measured rate is
                                   # ~10 min/design (remote MSA fetch + 3 models x 8
                                   # recycles), so 60 designs is ~10 h. 36000 would have
                                   # timed out just short of the finish.

cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }

# Bracket the first character so the grep's own argv, which ps prints, cannot
# match. An unbracketed alternative on this exact line aborted the first launch
# of Experiment 1 -- the fourth time this bug has appeared in this project.
busy(){ ps -eo args --width 400 | grep -cE '[e]nvs/colabfold/bin/colabfold_batch|[e]nvs/colabfold/bin/python .*scripts/mcts_design\.py'; }
folded(){ find "$KIT/results/rescore/folds" -name '*_scores_rank_001_*.json' 2>/dev/null | wc -l; }

{
say "=== CHAIN: exp1 (rescore) -> exp3 (partial_T) ==="
TARGET=$(grep -c '^>' "$KIT/results/rescore/to_fold.fasta" 2>/dev/null || echo 60)
say "waiting for $TARGET re-folds"

waited=0
while [ "$waited" -lt "$MAX_WAIT" ]; do
  n=$(folded); b=$(busy)
  if [ "$n" -ge "$TARGET" ] && [ "$b" -eq 0 ]; then
    say "exp 1 complete: $n/$TARGET folds, GPU idle"
    break
  fi
  # Stalled: nothing folding and the count is short. Do not silently start exp 3
  # on top of a half-finished exp 1 -- the whole point of exp 1 is that the
  # earlier comparison was invalid, and a partial one is no better.
  if [ "$b" -eq 0 ] && [ "$waited" -gt 600 ]; then
    say "STALLED: $n/$TARGET folds and no folding process. Not starting exp 3."
    say "  likely causes: mmseqs2 API rate limit, or colabfold_batch crashed."
    say "  inspect results/rescore.log, then rerun scripts/rescore_winners.sh"
    say "  (it skips designs already folded), and re-run this chain."
    exit 1
  fi
  sleep "$POLL"; waited=$((waited + POLL))
  [ $((waited % 1800)) -lt "$POLL" ] && say "  ... $n/$TARGET folds, ${waited}s elapsed"
done

if [ "$(folded)" -lt "$TARGET" ]; then
  say "TIMEOUT after ${MAX_WAIT}s at $(folded)/$TARGET. Not starting exp 3."
  exit 1
fi

say "--- running exp 1 analysis ---"
unset PYTHONPATH
/root/miniconda3/envs/colabfold/bin/python scripts/rescore_analyze.py --dir results/rescore \
  > "$KIT/results/RESCORE_RESULT.txt" 2>&1 \
  && say "  -> results/RESCORE_RESULT.txt" || say "  analysis FAILED (continuing)"

say "--- starting exp 3: partial_T sweep ---"
sleep 30                                   # let the card settle
setsid nohup bash "$KIT/scripts/partial_t_sweep.sh" > /dev/null 2>&1 < /dev/null &
sleep 60
if [ "$(busy)" -gt 0 ]; then
  say "exp 3 running"
else
  say "WARNING: exp 3 does not appear to be running; check results/partial_t.log"
fi
say "=== CHAIN DONE ==="
} >> "$LOG" 2>&1
