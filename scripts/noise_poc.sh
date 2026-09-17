#!/bin/bash
# Proof of concept: how fast does node-score noise fall with sequence draws?
#
# 10 backbones per target x 16 ProteinMPNN draws, folded at the search's own
# settings (1 model, 3 recycles, donated target MSA). Two targets deliberately:
# RBD is where the search failed and where K=3 noise is 0.192; IL-7Ra sits at
# 0.027 on the same measure. Running both says whether the noise is a property
# of the pipeline or of hard targets specifically, and that contrast is most of
# the value.
#
# 320 folds, about 5.5 h at the measured ~60 s/fold.
set -u
ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
KIT=$ROOT/calibration-kit
MPNN=$ROOT/external_tools/ProteinMPNN
OUT=$KIT/results/rl/noise_poc
LOG=$KIT/results/noise_poc.log
NDRAW=${NDRAW:-16}
CHUNK=${CHUNK:-16}

export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

cd "$KIT" || exit 1
say(){ echo "[$(date '+%F %T')] $*"; }
busy(){ ps -eo args --width 400 2>/dev/null \
        | grep -cE '[m]cts_design\.py|envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch|[p]rotein_mpnn_run'; }

mkdir -p "$OUT"
{
say "=== NOISE POC START ==="
say "waiting for the current chain (busy: $(busy))"
while [ "$(busy)" -gt 0 ]; do sleep 120; done
say "GPU free."

export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH:-}
RFPY=/root/miniconda3/envs/rfdiffusion/bin/python
echo '{"A": -1.0}' > "$OUT/bias.jsonl"

# target | mpnn source dir | designed chain | source a3m for the target MSA
for entry in \
  "rbd|$ROOT/real_test_results/rbd_200_batch/batch_0/proteinmpnn|A|$KIT/results/queue/v2refold/rbd/rbd__binder_0__s1.a3m" \
  "il7ra|$KIT/results/il7ra_real_pilot50/proteinmpnn|A|$KIT/results/queue/v2refold/il7ra/il7ra__il7ra_binder_0__s1.a3m"
do
  IFS='|' read -r T SRC CH A3M <<< "$entry"
  [ -f "$SRC/parsed.jsonl" ] || { say "SKIP $T (no parsed.jsonl)"; continue; }
  [ -f "$A3M" ] || { say "SKIP $T (no source a3m)"; continue; }
  if [ ! -d "$OUT/mpnn_$T/seqs" ]; then
    say "$T: sampling $NDRAW draws/backbone"
    ( cd "$MPNN" && "$RFPY" protein_mpnn_run.py \
        --jsonl_path "$SRC/parsed.jsonl" --chain_id_jsonl "$SRC/assigned.jsonl" \
        --out_folder "$OUT/mpnn_$T" --num_seq_per_target "$NDRAW" \
        --sampling_temp 0.2 --seed 11 --batch_size 1 \
        --bias_AA_jsonl "$OUT/bias.jsonl" ) >> "$LOG" 2>&1 || say "  $T sampling nonzero"
  fi
  say "$T: $(ls "$OUT/mpnn_$T/seqs"/*.fa 2>/dev/null | wc -l) fasta files"
done

say "building fold inputs (10 backbones/target, the variance-study set)"
/root/miniconda3/envs/colabfold/bin/python scripts/noise_poc_prepare.py --out "$OUT" \
  >> "$LOG" 2>&1 || { say "prepare FAILED"; exit 1; }
say "$(ls "$OUT/fold_in"/*.a3m 2>/dev/null | wc -l) fold inputs"

mkdir -p "$OUT/folds"
todo=()
for f in "$OUT/fold_in"/*.a3m; do
  b=$(basename "$f" .a3m)
  ls "$OUT/folds/${b}_scores_rank_001_"*.json >/dev/null 2>&1 || todo+=("$f")
done
say "${#todo[@]} to fold, chunks of $CHUNK"
i=0; c=0
while [ $i -lt ${#todo[@]} ]; do
  rm -rf "$OUT/.chunk"; mkdir -p "$OUT/.chunk"
  j=0
  while [ $j -lt "$CHUNK" ] && [ $i -lt ${#todo[@]} ]; do
    cp "${todo[$i]}" "$OUT/.chunk/"; i=$((i+1)); j=$((j+1))
  done
  c=$((c+1))
  colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
    "$OUT/.chunk" "$OUT/folds" >> "$LOG" 2>&1 || say "  chunk $c nonzero"
  say "  chunk $c -- $(ls "$OUT/folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l)/${#todo[@]}"
done
rm -rf "$OUT/.chunk"

/root/miniconda3/envs/colabfold/bin/python scripts/noise_curve.py --dir "$OUT" \
  > "$KIT/results/NOISE_CURVE_RESULT.txt" 2>&1 \
  && say "analysed -> results/NOISE_CURVE_RESULT.txt" || say "analysis FAILED"
say "=== NOISE POC COMPLETE ==="
} >> "$LOG" 2>&1
