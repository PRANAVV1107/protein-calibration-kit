#!/bin/bash
# Run this ON the rented box. Folds everything pending, chunked and resumable.
#
# Chunked because a single colabfold_batch call over many inputs has been
# OOM-killed three times in this project. Resumable because a spot instance can
# be preempted at any moment -- completed folds write a scores json and are
# skipped, so relaunching after preemption costs only the chunk in flight.
#
# Results are tarred per job as they complete, so a preempted box still has
# finished work ready to pull rather than everything or nothing.
set -u
WORK=${WORK:-/workspace}
CHUNK=${CHUNK:-40}
export PATH="$WORK/localcolabfold/colabfold-conda/bin:$PATH"
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
cd "$WORK" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }

# highest value first: a preempted instance should have finished what matters
ORDER="reward_bench expert_r1 v2_pdl1 v2_rbd null_panel"

for job in $ORDER; do
  src="$WORK/pending/$job"; dst="$WORK/out/$job"
  [ -d "$src" ] || continue
  mkdir -p "$dst"
  todo=()
  for f in "$src"/*.a3m; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .a3m)
    ls "$dst/${b}_scores_rank_001_"*.json >/dev/null 2>&1 || todo+=("$f")
  done
  [ ${#todo[@]} -eq 0 ] && { say "SKIP $job (complete)"; continue; }
  say "START $job: ${#todo[@]} folds, chunks of $CHUNK"
  i=0; c=0
  while [ $i -lt ${#todo[@]} ]; do
    rm -rf "$WORK/.chunk"; mkdir -p "$WORK/.chunk"
    j=0
    while [ $j -lt "$CHUNK" ] && [ $i -lt ${#todo[@]} ]; do
      cp "${todo[$i]}" "$WORK/.chunk/"; i=$((i+1)); j=$((j+1))
    done
    c=$((c+1))
    colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
      "$WORK/.chunk" "$dst" >> "$WORK/$job.log" 2>&1 || say "  chunk $c nonzero"
    say "  $job chunk $c -- $(ls "$dst"/*_scores_rank_001_*.json 2>/dev/null | wc -l)/${#todo[@]}"
  done
  # ship only the scores: the pdb/png output is many times larger and unused
  tar -czf "$WORK/results_$job.tar.gz" -C "$dst" $(cd "$dst" && ls *_scores_rank_001_*.json 2>/dev/null | head -100000) 2>/dev/null
  say "DONE $job -> results_$job.tar.gz"
done
rm -rf "$WORK/.chunk"
say "=== ALL DONE ==="
ls -lh "$WORK"/results_*.tar.gz 2>/dev/null
