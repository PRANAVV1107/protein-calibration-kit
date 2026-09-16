#!/bin/bash
# Expert iteration, round 1: does the DPO-tuned model SAMPLE better sequences?
#
# Why this and not MCTS
# ---------------------
# MCTS proper is a 5-7 day build, and there is a deeper obstacle: with the ~500
# folds affordable here, a tree over 60 sequence positions and 20 amino acids is
# so sparsely sampled that search has no advantage over random draws. MCTS needs
# many evaluations or a value function, and value_probe.py measured the cheap
# value function failing (within-backbone rho about -0.06).
#
# But round 1 of the expert-iteration loop needs neither, and it answers the
# question everything else depends on. If sequences sampled from the tuned model
# are no better than sequences from the untuned model, then no tree search built
# on top of that policy will help, and the whole direction is dead cheaply
# rather than after a week of engineering.
#
# It is also the prospective test the DPO result is missing. eval_dpo.py showed
# the tuned model RANKS held-out pairs at 0.490, chance. Ranking existing
# candidates and generating good ones are different abilities, and a model can
# lose the first while keeping the second.
#
# Design
# ------
# 10 RBD backbones, held out of DPO training entirely. 6 sequences sampled from
# each of the untuned and tuned models, same backbones, same sampling
# temperature, same alanine bias -- so the only difference is the weights.
# 120 folds at --num-models 1 --num-recycle 3 with the precomputed target MSA,
# the configuration measured at 63 s/design.
#
# The comparison is PAIRED by backbone. Backbone quality varies far more than
# the effect being measured (this project put between-backbone variance at
# 22-41% of total), so an unpaired comparison of 60 against 60 would be mostly
# noise from which backbones landed where.
#
# Chunked at 30. A single colabfold_batch call over many inputs has been
# OOM-killed three times in this project, most recently at 85/165.
set -u

ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
KIT=$ROOT/calibration-kit
MPNN=$ROOT/external_tools/ProteinMPNN
SRC=$ROOT/real_test_results/rbd_200_batch/batch_0/proteinmpnn
OUT=$KIT/results/rl/expert_r1
LOG=$KIT/results/expert_r1.log
TUNED=$KIT/results/rl/dpo_noRBD/dpo_epoch2.pt
NDRAW=${NDRAW:-6}
CHUNK=${CHUNK:-30}

export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

say(){ echo "[$(date '+%F %T')] $*"; }
busy(){ ps -eo args --width 400 2>/dev/null \
        | grep -cE '[m]d_stability\.py|envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch'; }

mkdir -p "$OUT"
{
say "=== EXPERT ITERATION ROUND 1 ==="
say "waiting for MD and the benchmark resume to finish (busy: $(busy))"
while [ "$(busy)" -gt 0 ]; do sleep 120; done
say "GPU free."

# ---- sample from each model ----
cd "$MPNN" || exit 1
export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH:-}
RFPY=/root/miniconda3/envs/rfdiffusion/bin/python
echo '{"A": -1.0}' > "$OUT/bias_A1.jsonl"

for arm in base tuned; do
  if [ -d "$OUT/mpnn_$arm/seqs" ]; then say "skip sampling $arm (exists)"; continue; fi
  mkdir -p "$OUT/mpnn_$arm"
  EXTRA=""
  [ "$arm" = "tuned" ] && EXTRA="--path_to_model_weights $(dirname $TUNED) --model_name $(basename $TUNED .pt)"
  say "sampling $arm ($NDRAW draws/backbone)"
  $RFPY protein_mpnn_run.py \
    --jsonl_path "$SRC/parsed.jsonl" \
    --chain_id_jsonl "$SRC/assigned.jsonl" \
    --out_folder "$OUT/mpnn_$arm" \
    --num_seq_per_target "$NDRAW" --sampling_temp 0.2 --seed 37 --batch_size 1 \
    --bias_AA_jsonl "$OUT/bias_A1.jsonl" $EXTRA >> "$LOG" 2>&1 \
    || say "  $arm sampling returned nonzero"
  say "  $arm: $(ls $OUT/mpnn_$arm/seqs/*.fa 2>/dev/null | wc -l) fasta files"
done

# ---- build fold inputs for the 10 held-out backbones ----
say "building fold inputs"
/root/miniconda3/envs/colabfold/bin/python "$KIT/scripts/expert_r1_prepare.py" \
  --out "$OUT" >> "$LOG" 2>&1 || { say "prepare FAILED"; exit 1; }
n_in=$(ls "$OUT/fold_in"/*.a3m 2>/dev/null | wc -l)
say "$n_in fold inputs"

# ---- fold, chunked ----
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
  say "  chunk $c -- $(ls $OUT/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l) folded"
done
rm -rf "$OUT/.chunk"

# ---- compare ----
/root/miniconda3/envs/colabfold/bin/python "$KIT/scripts/expert_r1_analyze.py" \
  --dir "$OUT" > "$KIT/results/EXPERT_R1_RESULT.txt" 2>&1 \
  && say "analysed -> results/EXPERT_R1_RESULT.txt" || say "analysis FAILED"
say "=== ROUND 1 COMPLETE ==="
} >> "$LOG" 2>&1
