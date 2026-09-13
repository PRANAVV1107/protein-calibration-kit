#!/bin/bash
# Generic de novo binder campaign pipeline.
#
#   RFdiffusion -> ProteinMPNN (multi-draw, Ala-biased) -> composition filter
#   -> ColabFold -> select on shaped reward
#
# Differences from the original single-draw flow, all based on measurements
# from the IL-7Ra pilot (see IL7RA_CAMPAIGN.md):
#
#   * ProteinMPNN draws N sequences per backbone instead of 1. ipTM proved
#     highly sensitive to the specific sequence: resampling ~10-15% of
#     positions on an UNCHANGED backbone swung ipTM by up to 0.27, and the
#     within-backbone spread across 3 draws was ~0.2-0.25. One draw is a
#     lottery ticket.
#
#   * An alanine logit bias (-1.0) is applied. Default settings collapsed one
#     design to 40% alanine (another to 59% with a 16-residue run). Raising
#     sampling temperature alone did NOT fix this (mean Ala 20.6% -> 20.4%
#     going 0.1 -> 0.2); the bias did (-> 8.8%, i.e. natural abundance).
#
#   * Sequences are filtered on composition BEFORE folding, because design is
#     cheap (~1.5 s) and folding is not (~4.5 min high-accuracy).
#
#   * Final ranking uses shaped reward (ipTM - composition penalty), the same
#     objective the RL refinement loop optimises.
#
# Measured effect on the IL-7Ra top 5 backbones: mean best-per-backbone ipTM
# 0.840 (1 draw) -> 0.860 (best of 8), with composition penalties falling from
# 0.307 to ~0.000. Clean sequence at no cost in score.
#
# Usage:
#   run_campaign_pipeline.sh <name> <target_pdb> <contig> <hotspots> <target_seq_file> [n_designs] [n_samples] [keep] [fold_mode]
#
# Example (IL-7Ra):
#   run_campaign_pipeline.sh il7ra \
#     benchmark/il7ra_target_only.pdb \
#     '[B17-209/0 50-120]' '[B31,B77,B138,B192]' \
#     benchmark/il7ra_target_seq.txt \
#     50 8 2 highacc
#
# fold_mode: "highacc" (real MSA, 3 models, 8 recycles) or "fast"
#            (single-sequence, 1 model, 3 recycles).
#            NOTE: the fast screen showed NO ranking signal within its top band
#            in the IL-7Ra pilot (Spearman rho=0.06, n=10). Prefer highacc
#            unless you are deliberately triaging a very large backbone set.
set -e

NAME=${1:?need campaign name}
TARGET_PDB=${2:?need target pdb}
CONTIG=${3:?need contig, e.g. [B17-209/0 50-120]}
HOTSPOTS=${4:?need hotspots, e.g. [B31,B77,B138,B192]}
TARGET_FILE=${5:?need target seq file}
N_DESIGNS=${6:-50}
N_SAMPLES=${7:-8}
KEEP=${8:-2}
FOLD_MODE=${9:-highacc}
MAX_PENALTY=${MAX_PENALTY:-0.05}

ROOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp
KIT=$ROOT/calibration-kit
RFD=$ROOT/external_tools/RFdiffusion
MPNN=$ROOT/external_tools/ProteinMPNN
OUT=$KIT/results/${NAME}_campaign

mkdir -p "$OUT"/{rfdiffusion,mpnn,folds}
TARGET_PDB=$(cd "$(dirname "$TARGET_PDB")" && pwd)/$(basename "$TARGET_PDB")
TARGET_FILE=$(cd "$(dirname "$TARGET_FILE")" && pwd)/$(basename "$TARGET_FILE")
echo '{"A": -1.0}' > "$OUT/bias_A1.jsonl"

echo "======================================================================"
echo "CAMPAIGN: $NAME"
echo "======================================================================"
echo "target pdb : $TARGET_PDB"
echo "contig     : $CONTIG"
echo "hotspots   : $HOTSPOTS"
echo "backbones  : $N_DESIGNS   samples/backbone: $N_SAMPLES   keep: $KEEP"
echo "fold mode  : $FOLD_MODE   max composition penalty: $MAX_PENALTY"
echo "Start: $(date)"

# ---------------- Stage 1: RFdiffusion ----------------
export PATH=/root/miniconda3/bin:$PATH
export DGLBACKEND=pytorch
export PYTHONPATH=$ROOT/external_tools:$PYTHONPATH
export CONDA_DEFAULT_ENV=rfdiffusion
export CONDA_PREFIX=/root/miniconda3/envs/rfdiffusion
export PATH=/root/miniconda3/envs/rfdiffusion/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH}

echo ""
echo "[$(date)] Stage 1/5: RFdiffusion ($N_DESIGNS backbones)..."
cd "$RFD"
python scripts/run_inference.py \
  inference.output_prefix="$OUT/rfdiffusion/${NAME}_binder" \
  inference.model_directory_path="$RFD/models" \
  inference.input_pdb="$TARGET_PDB" \
  "contigmap.contigs=$CONTIG" \
  "ppi.hotspot_res=$HOTSPOTS" \
  inference.num_designs=$N_DESIGNS \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0 > "$OUT/rfdiffusion.log" 2>&1
echo "  backbones: $(ls "$OUT/rfdiffusion"/${NAME}_binder_*.pdb 2>/dev/null | wc -l)"

# ---------------- Stage 2: ProteinMPNN, multi-draw + Ala bias ----------------
echo ""
echo "[$(date)] Stage 2/5: ProteinMPNN ($N_SAMPLES samples/backbone, temp 0.2, A bias -1.0)..."
mkdir -p "$OUT/mpnn/pdbs"
cp "$OUT/rfdiffusion"/${NAME}_binder_*.pdb "$OUT/mpnn/pdbs/"
cd "$MPNN"

python helper_scripts/parse_multiple_chains.py \
  --input_path "$OUT/mpnn/pdbs" --output_path "$OUT/mpnn/parsed.jsonl" >> "$OUT/mpnn.log" 2>&1
python helper_scripts/assign_fixed_chains.py \
  --input_path "$OUT/mpnn/parsed.jsonl" --output_path "$OUT/mpnn/assigned.jsonl" \
  --chain_list "A" >> "$OUT/mpnn.log" 2>&1
python protein_mpnn_run.py \
  --jsonl_path "$OUT/mpnn/parsed.jsonl" \
  --chain_id_jsonl "$OUT/mpnn/assigned.jsonl" \
  --out_folder "$OUT/mpnn/output" \
  --num_seq_per_target "$N_SAMPLES" --sampling_temp "0.2" --seed 37 --batch_size 1 \
  --path_to_model_weights "$MPNN/vanilla_model_weights" \
  --bias_AA_jsonl "$OUT/bias_A1.jsonl" >> "$OUT/mpnn.log" 2>&1
echo "  sequence files: $(ls "$OUT/mpnn/output/seqs"/*.fa 2>/dev/null | wc -l)"

# ---------------- Stage 3: composition filter ----------------
export CONDA_DEFAULT_ENV=colabfold
export CONDA_PREFIX=/root/miniconda3/envs/colabfold
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/colabfold/lib:${LD_LIBRARY_PATH}

echo ""
echo "[$(date)] Stage 3/5: composition filter (keep $KEEP/backbone, penalty <= $MAX_PENALTY)..."
cd "$KIT"
python scripts/seq_select.py filter \
  --seqs-dir "$OUT/mpnn/output/seqs" \
  --per-backbone "$KEEP" \
  --max-penalty "$MAX_PENALTY" \
  --target-seq-file "$TARGET_FILE" \
  --out-fasta "$OUT/to_fold.fasta" \
  --out-manifest "$OUT/manifest.csv" \
  --minutes-per-fold $([ "$FOLD_MODE" = "fast" ] && echo 2.5 || echo 4.5)

# ---------------- Stage 4: fold ----------------
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

echo ""
echo "[$(date)] Stage 4/5: ColabFold ($FOLD_MODE)..."
if [ "$FOLD_MODE" = "fast" ]; then
  colabfold_batch --msa-mode single_sequence --model-type alphafold2_multimer_v3 \
    --num-models 1 --num-recycle 3 "$OUT/to_fold.fasta" "$OUT/folds" > "$OUT/fold.log" 2>&1
else
  colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
    --num-models 3 --num-recycle 8 "$OUT/to_fold.fasta" "$OUT/folds" > "$OUT/fold.log" 2>&1
fi
echo "  folds completed: $(ls "$OUT/folds"/*_scores_rank_001_*.json 2>/dev/null | wc -l)"

# ---------------- Stage 5: select on shaped reward ----------------
echo ""
echo "[$(date)] Stage 5/5: ranking by shaped reward..."
python scripts/seq_select.py select \
  --fold-dir "$OUT/folds" \
  --manifest "$OUT/manifest.csv" \
  --out-csv "$OUT/ranked.csv"

echo ""
echo "======================================================================"
echo "CAMPAIGN COMPLETE: $OUT/ranked.csv"
echo "======================================================================"
echo "Reminder: solo-verify the top candidate independently before trusting it,"
echo "and run negative controls (scrambled sequence, off-target) as in"
echo "IL7RA_CAMPAIGN.md before calling anything a hit."
echo "Finish: $(date)"
