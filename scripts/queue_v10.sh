#!/bin/bash
# Overnight run, reprioritised by the preliminary Adaptyv result.
#
# On 272 of 365 measured outcomes, the cheap hybrid screen does NOT separate
# real binders from real non-binders:
#
#     ipTM         AUC 0.512  [0.422, 0.601]
#     ipSAE_min    AUC 0.450  [0.361, 0.535]
#     pDockQ2_min  AUC 0.434  [0.344, 0.522]
#
# Every interval contains 0.5. (pLDDT and pTM look inverted at 0.383, but that
# is a length confound: binders average 127 aa against 78 aa, and length
# anti-correlates with both. The three interface metrics are length-independent,
# so their null is genuine.)
#
# Those folds used 1 model, 3 recycles and a precomputed target alignment -- the
# CHEAP screen. Adaptyv's own reported ipTM, from a fuller pipeline, reached
# 0.574. So the open question is whether this is a failure of the cheap screen
# or of the metrics themselves, and it is answerable: refold a balanced subset
# at full accuracy and compare AUC on identical designs.
#
# That question comes first tonight. Ranking the whole design space (the v2
# refold) by a criterion that may not predict anything would be premature.
#
#   1  adaptyv hybrid, finish   93 remaining   completes the cheap-screen arm
#   2  adaptyv HIGH ACCURACY   165 folds       the controlled comparison
#   3  v2 refold               881 folds       only meaningful once a metric is validated
#   4  null panel              420 folds
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v10.log
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
say "=== QUEUE V10 START ==="
while pgrep -f 'colabfold_batch --model' > /dev/null; do sleep 20; done
say "GPU free."

# ---- 1. finish the cheap-screen arm ----
chunked adaptyv_val "$Q/adaptyv_val" "$Q/adaptyv_val_folds"
cd "$KIT"
"$PY" scripts/analyze_adaptyv.py > results/ADAPTYV_HYBRID_RESULT.txt 2>&1 \
  && say "hybrid arm analysed -> results/ADAPTYV_HYBRID_RESULT.txt"

# ---- 2. the controlled comparison: same designs, full accuracy ----
mkdir -p "$Q/adaptyv_highacc_folds"
say "START adaptyv_highacc (165 folds, 3 models / 8 recycles / full MSA)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 3 --num-recycle 8 \
  "$Q/adaptyv_highacc.fasta" "$Q/adaptyv_highacc_folds" > "$Q/adaptyv_highacc.log" 2>&1 \
  && say "DONE  adaptyv_highacc ($(count "$Q/adaptyv_highacc_folds"))" \
  || say "PARTIAL adaptyv_highacc ($(count "$Q/adaptyv_highacc_folds") done)"

"$PY" scripts/analyze_adaptyv.py --fold-dir "$Q/adaptyv_highacc_folds" \
      --manifest "$Q/adaptyv_highacc_manifest.csv" \
      --out results/adaptyv_highacc_validation.csv \
      > results/ADAPTYV_HIGHACC_RESULT.txt 2>&1 \
  && say "high-accuracy arm analysed -> results/ADAPTYV_HIGHACC_RESULT.txt"

# ---- 3. the full design space ----
for t in il7ra pdl1 rbd; do
  [ -d "$Q/v2refold/$t" ] && chunked "v2_$t" "$Q/v2refold/$t" "$Q/v2refold_folds_$t"
done

# ---- 4. constructed negatives ----
chunked null_panel "$Q/null_panel" "$Q/null_panel_folds"

say "=== QUEUE V10 COMPLETE ==="
} >> "$LOG" 2>&1
