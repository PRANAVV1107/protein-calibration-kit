#!/bin/bash
# Continuation of queue_v10 after the gaming pause. Same work, resumed.
#
# Why a new file rather than relaunching v10: bash reads a script by byte
# offset while executing it, so editing a running script corrupts the run --
# which happened here once. The convention since is a new numbered file each
# time. v10 also begins by folding the cheap adaptyv arm, which is finished,
# so restarting it would re-walk completed work before reaching the part that
# matters.
#
# Remaining, in priority order:
#   1  adaptyv HIGH ACCURACY   the controlled comparison, blocks everything
#   2  v2 refold               881 folds, the full design space
#   3  null panel              420 folds, constructed negatives
#
# The high-accuracy arm resumes in place: ColabFold writes a .done.txt per
# query and skips those, so re-invoking on the same fasta picks up where it
# stopped. The remaining designs are the long ones (58-250 aa against a ~200
# aa target) because ColabFold sorts queries shortest-first -- so expect
# slower folds and higher memory than the first 70.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v11.log
PY=/root/miniconda3/envs/colabfold/bin/python

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

CHUNK=${CHUNK:-40}
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

chunked(){                      # name, in_dir, out_dir  -- a3m inputs
  local name="$1" src="$2" dst="$3"
  mkdir -p "$dst"
  local todo=()
  for f in "$src"/*.a3m; do
    local b; b=$(basename "$f" .a3m)
    ls "$dst/${b}_scores_rank_001_"*.json >/dev/null 2>&1 || todo+=("$f")
  done
  [ ${#todo[@]} -eq 0 ] && { say "SKIP  $name (complete)"; return; }
  say "START $name: ${#todo[@]} remaining, chunks of $CHUNK"
  local i=0 n=0
  while [ $i -lt ${#todo[@]} ]; do
    rm -rf "$Q/.c_$name"; mkdir -p "$Q/.c_$name"
    local j=0
    while [ $j -lt "$CHUNK" ] && [ $i -lt ${#todo[@]} ]; do
      cp "${todo[$i]}" "$Q/.c_$name/"; i=$((i+1)); j=$((j+1))
    done
    n=$((n+1))
    colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
      "$Q/.c_$name" "$dst" >> "$Q/$name.log" 2>&1 || say "  chunk $n nonzero (continuing)"
    say "  chunk $n -- $(count "$dst") total"
  done
  rm -rf "$Q/.c_$name"
  say "DONE  $name ($(count "$dst"))"
}

{
say "=== QUEUE V11 START (resumed after pause) ==="
say "high-accuracy arm at $(count "$Q/adaptyv_highacc_folds")/165"

# ---- 1. finish the controlled comparison ----
# Retried once: the remaining designs are the long ones, and a 450-residue
# complex at 3 models / 8 recycles is the configuration that OOM-killed an
# earlier run. A second pass picks up anything the first pass lost, since
# completed queries are skipped.
for attempt in 1 2; do
  n=$(count "$Q/adaptyv_highacc_folds")
  [ "$n" -ge 165 ] && break
  say "highacc attempt $attempt (at $n/165)"
  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
    --num-models 3 --num-recycle 8 \
    "$Q/adaptyv_highacc.fasta" "$Q/adaptyv_highacc_folds" \
    >> "$Q/adaptyv_highacc.log" 2>&1 \
    && say "  attempt $attempt clean" || say "  attempt $attempt nonzero"
  say "  now $(count "$Q/adaptyv_highacc_folds")/165"
done

"$PY" scripts/analyze_adaptyv.py --fold-dir "$Q/adaptyv_highacc_folds" \
      --manifest "$Q/adaptyv_highacc_manifest.csv" \
      --out results/adaptyv_highacc_validation.csv \
      > results/ADAPTYV_HIGHACC_RESULT.txt 2>&1 \
  && say "high-accuracy arm analysed -> results/ADAPTYV_HIGHACC_RESULT.txt"

"$PY" scripts/compare_arms.py --boot 10000 \
      > results/ARM_COMPARISON_FINAL.txt 2>&1 \
  && say "controlled comparison -> results/ARM_COMPARISON_FINAL.txt"

# ---- 2. the full design space ----
for t in il7ra pdl1 rbd; do
  [ -d "$Q/v2refold/$t" ] && chunked "v2_$t" "$Q/v2refold/$t" "$Q/v2refold_folds_$t"
done

# ---- 3. constructed negatives ----
chunked null_panel "$Q/null_panel" "$Q/null_panel_folds"

say "=== QUEUE V11 COMPLETE ==="
} >> "$LOG" 2>&1
