#!/bin/bash
# Sequence-draw variance study across three targets.
#
# Question: on a FIXED backbone, how much does high-accuracy ipTM vary across
# independent ProteinMPNN sequence draws, and how much does best-of-N gain
# over a single draw?
#
# The IL-7Ra pilot saw a ~0.2-0.25 within-backbone spread across only 3 draws
# on 5 backbones (n too small to generalise). This scales that to 3 targets x
# 10 backbones x 6 draws = 180 high-accuracy folds.
#
# Reuses existing RFdiffusion backbones from all three campaigns -- no new
# backbone generation required.
#
# Waits for any in-flight ColabFold job before starting, so it never contends
# for VRAM (6GB card; concurrent jobs have OOM'd before).
set -e

ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
KIT=$ROOT/calibration-kit
MPNN=$ROOT/external_tools/ProteinMPNN
OUT=$KIT/results/variance_study

N_DRAWS=${N_DRAWS:-6}
N_BACKBONES=${N_BACKBONES:-10}

mkdir -p "$OUT"

echo "======================================================================"
echo "SEQUENCE-DRAW VARIANCE STUDY"
echo "======================================================================"
echo "targets: il7ra, rbd, pdl1   backbones/target: $N_BACKBONES   draws/backbone: $N_DRAWS"
echo "Start: $(date)"

# ---- wait for any in-flight fold so we never contend for VRAM ----
# Match the interpreter's full path to the binary, NOT the bare name: any
# monitoring shell running `pgrep -f colabfold_batch` carries that string in
# its own command line, so a bare-name match makes watchers match each other
# and deadlock. (This happened; it idled the GPU for ~2.5 h.)
REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
  echo ""
  echo "[$(date)] waiting for in-flight ColabFold job to finish..."
  while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
  echo "[$(date)] GPU free, proceeding."
fi

export PATH=/root/miniconda3/bin:$PATH
export DGLBACKEND=pytorch

# ---------------- Stage 1: ProteinMPNN, N draws per backbone ----------------
export CONDA_DEFAULT_ENV=rfdiffusion
export CONDA_PREFIX=/root/miniconda3/envs/rfdiffusion
export PATH=/root/miniconda3/envs/rfdiffusion/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH}

echo '{"A": -1.0}' > "$OUT/bias_A1.jsonl"
cd "$MPNN"

# name|parsed_dir|designed_chain
TARGETS=(
  "il7ra|$KIT/results/il7ra_real_pilot50/proteinmpnn|A"
  "rbd|$ROOT/real_test_results/rbd_200_batch/batch_0/proteinmpnn|A"
  "pdl1|$ROOT/real_test_results/pdl1_pilot50/proteinmpnn|B"
)

for entry in "${TARGETS[@]}"; do
  IFS='|' read -r name pdir chain <<< "$entry"
  echo ""
  echo "[$(date)] Stage 1 [$name]: ProteinMPNN $N_DRAWS draws/backbone (temp 0.2, A bias -1.0, designed chain $chain)..."
  mkdir -p "$OUT/$name"
  python protein_mpnn_run.py \
    --jsonl_path "$pdir/parsed.jsonl" \
    --chain_id_jsonl "$pdir/assigned.jsonl" \
    --out_folder "$OUT/$name/mpnn" \
    --num_seq_per_target "$N_DRAWS" --sampling_temp "0.2" --seed 37 --batch_size 1 \
    --path_to_model_weights "$MPNN/vanilla_model_weights" \
    --bias_AA_jsonl "$OUT/bias_A1.jsonl" > "$OUT/$name/mpnn.log" 2>&1
  echo "   backbone files: $(ls "$OUT/$name/mpnn/seqs"/*.fa 2>/dev/null | wc -l)"
done

# ---------------- Stage 2: build fold set (seeded-random backbone subset) ----
export CONDA_DEFAULT_ENV=colabfold
export CONDA_PREFIX=/root/miniconda3/envs/colabfold
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/colabfold/lib:${LD_LIBRARY_PATH}

echo ""
echo "[$(date)] Stage 2: building fold set..."
cd "$KIT"
python scripts/variance_prepare.py \
  --study-dir "$OUT" --n-backbones "$N_BACKBONES" --draws "$N_DRAWS"

# ---------------- Stage 3: fold everything ----------------
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

for name in il7ra rbd pdl1; do
  echo ""
  echo "[$(date)] Stage 3 [$name]: high-accuracy folding..."
  mkdir -p "$OUT/$name/folds"
  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
    --num-models 3 --num-recycle 8 \
    "$OUT/$name/to_fold.fasta" "$OUT/$name/folds" > "$OUT/$name/fold.log" 2>&1 || \
    echo "   [WARN] $name fold exited non-zero; continuing (partial results usable)"
  echo "   folds done: $(ls "$OUT/$name/folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l)"
done

# ---------------- Stage 4: analyse ----------------
echo ""
echo "[$(date)] Stage 4: analysis..."
python scripts/variance_analyze.py --study-dir "$OUT"

echo ""
echo "======================================================================"
echo "VARIANCE STUDY COMPLETE -> $OUT/variance_results.csv"
echo "======================================================================"
echo "Finish: $(date)"
