#!/bin/bash
# Multi-sequence refinement for a set of already-promising backbones.
#
# Rationale: ipTM is sensitive to the specific ProteinMPNN sequence, not just
# the backbone -- resampling ~10-15% of positions on an unchanged backbone
# moved ipTM by up to 0.27 in the IL-7Ra pilot. So instead of trusting one
# sequence per backbone, draw many (cheap), keep only those with clean
# composition (free), and spend fold budget (expensive) on those.
#
# Usage:
#   refine_sequences.sh <pilot_dir> <out_dir> <target_seq_file> [n_samples] [per_backbone]
#
# <pilot_dir> must contain proteinmpnn/parsed.jsonl and proteinmpnn/assigned.jsonl
# from the original pipeline run (so the same backbones/chain assignment are reused).
set -e

PILOT=${1:?need pilot_dir}
OUT=${2:?need out_dir}
TARGET_FILE=${3:?need target_seq_file}
N_SAMPLES=${4:-8}
PER_BACKBONE=${5:-3}
MAX_PENALTY=${MAX_PENALTY:-0.05}

ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
MPNN=$ROOT/external_tools/ProteinMPNN
KIT=$ROOT/calibration-kit

# Resolve to absolute paths: later stages cd into $MPNN / $KIT, which would
# otherwise break any relative path the caller passed in.
mkdir -p "$OUT"
PILOT=$(cd "$PILOT" && pwd)
OUT=$(cd "$OUT" && pwd)
TARGET_FILE=$(cd "$(dirname "$TARGET_FILE")" && pwd)/$(basename "$TARGET_FILE")
echo '{"A": -1.0}' > "$OUT/bias_A1.jsonl"

echo "======================================================================"
echo "SEQUENCE REFINEMENT"
echo "======================================================================"
echo "Backbones from : $PILOT"
echo "Samples/backbone: $N_SAMPLES   keep/backbone: $PER_BACKBONE"
echo "Alanine bias   : -1.0    temp: 0.2    max composition penalty: $MAX_PENALTY"
echo "Start: $(date)"
echo ""

# ---- Stage 1: ProteinMPNN, many samples per backbone, with alanine bias ----
export PATH=/root/miniconda3/bin:$PATH
export CONDA_DEFAULT_ENV=rfdiffusion
export CONDA_PREFIX=/root/miniconda3/envs/rfdiffusion
export PATH=/root/miniconda3/envs/rfdiffusion/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH}

echo "[$(date)] Stage 1: ProteinMPNN ($N_SAMPLES samples/backbone, A bias -1.0, temp 0.2)..."
cd "$MPNN"
python protein_mpnn_run.py \
  --jsonl_path "$PILOT/proteinmpnn/parsed.jsonl" \
  --chain_id_jsonl "$PILOT/proteinmpnn/assigned.jsonl" \
  --out_folder "$OUT/mpnn" \
  --num_seq_per_target "$N_SAMPLES" --sampling_temp "0.2" --seed 37 --batch_size 1 \
  --path_to_model_weights "$MPNN/vanilla_model_weights" \
  --bias_AA_jsonl "$OUT/bias_A1.jsonl" > "$OUT/mpnn.log" 2>&1
echo "  sequences generated: $(ls "$OUT/mpnn/seqs"/*.fa 2>/dev/null | wc -l) backbone files"

# ---- Stage 2: composition filter (no GPU needed) ----
export CONDA_DEFAULT_ENV=colabfold
export CONDA_PREFIX=/root/miniconda3/envs/colabfold
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/colabfold/lib:${LD_LIBRARY_PATH}

echo ""
echo "[$(date)] Stage 2: composition filter..."
cd "$KIT"
BB_ARG=""
if [ -n "${BACKBONES:-}" ]; then BB_ARG="--backbones $BACKBONES"; fi

python scripts/seq_select.py filter \
  --seqs-dir "$OUT/mpnn/seqs" \
  --per-backbone "$PER_BACKBONE" \
  --max-penalty "$MAX_PENALTY" \
  --target-seq-file "$TARGET_FILE" \
  --out-fasta "$OUT/to_fold.fasta" \
  --out-manifest "$OUT/manifest.csv" \
  $BB_ARG

echo ""
echo "======================================================================"
echo "Stages 1-2 done. Fold input ready: $OUT/to_fold.fasta"
echo "======================================================================"
echo "Next (high-accuracy fold, then select):"
echo ""
echo "  export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false \\"
echo "         XLA_PYTHON_CLIENT_MEM_FRACTION=4.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform"
echo "  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \\"
echo "    --num-models 3 --num-recycle 8 $OUT/to_fold.fasta $OUT/folds"
echo ""
echo "  python scripts/seq_select.py select --fold-dir $OUT/folds \\"
echo "    --manifest $OUT/manifest.csv --out-csv $OUT/ranked.csv"
echo ""
echo "Finish: $(date)"
