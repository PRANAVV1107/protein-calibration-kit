#!/bin/bash
# Real IL-7Ra binder design pipeline: RFdiffusion -> ProteinMPNN -> ColabFold fast-screen
# Mirrors the PD-L1 pilot methodology (50-design pilot, hotspot-guided).
set -e

export PATH=/root/miniconda3/bin:$PATH
export DGLBACKEND=pytorch
export PYTHONPATH=/mnt/c/Users/prana/Downloads/proteinfoldingexp/external_tools:$PYTHONPATH

REPO=/mnt/c/Users/prana/Downloads/proteinfoldingexp/external_tools/RFdiffusion
MODELS=$REPO/models
MPNN=/mnt/c/Users/prana/Downloads/proteinfoldingexp/external_tools/ProteinMPNN
INPUT_PDB=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit/benchmark/il7ra_target_only.pdb
OUT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit/results/il7ra_real_pilot50
TARGET_SEQ="DYSFSCYSQLEVNGSQHSLTCAFEDPDVNTTNLEFEICGALVEVKCLNFRKLQEIYFIETKKFLLIGKSNICVKVGEKSLTCKKIDLTTIVKPEAPFDLSVVYREGANDFVVTFNTSHLQKKYVKVLMHDVAYRQEKDENKWTHVNLSSTKLTLLQRKLQPAAMYEIKVRSIPDHYFKGFWSEWSPSYYFRTP"
NUM_DESIGNS=50

mkdir -p "$OUT/rfdiffusion" "$OUT/proteinmpnn/pdbs" "$OUT/proteinmpnn/output" "$OUT/colabfold_fast"

echo "======================================================================"
echo "IL-7Ra BINDER DESIGN PILOT (50 designs, real RFdiffusion + ProteinMPNN)"
echo "======================================================================"
echo "Start: $(date)"
echo "Target: 3DI3 chain B (IL-7Ra ectodomain), residues 17-209, 193 aa"
echo "Hotspots: B31 (Ser), B77 (Lys), B138 (Lys), B192 (Tyr) -- confirmed vs PDB"
echo ""

# ---- Stage 1: RFdiffusion (blind hotspot-guided generation) ----
export CONDA_DEFAULT_ENV=rfdiffusion
export CONDA_PREFIX=/root/miniconda3/envs/rfdiffusion
export PATH=/root/miniconda3/envs/rfdiffusion/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH}

echo "[$(date)] Stage 1: RFdiffusion ($NUM_DESIGNS backbones, hotspot-guided)..."
cd "$REPO"
python scripts/run_inference.py \
  inference.output_prefix="${OUT}/rfdiffusion/il7ra_binder" \
  inference.model_directory_path="${MODELS}" \
  inference.input_pdb="${INPUT_PDB}" \
  'contigmap.contigs=[B17-209/0 50-120]' \
  'ppi.hotspot_res=[B31,B77,B138,B192]' \
  inference.num_designs=$NUM_DESIGNS \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0 \
  2>&1 | tail -50

RF_COUNT=$(ls "$OUT/rfdiffusion"/il7ra_binder_*.pdb 2>/dev/null | wc -l)
echo "RFdiffusion: $RF_COUNT backbones generated"

# ---- Stage 2: ProteinMPNN ----
echo "[$(date)] Stage 2: ProteinMPNN sequence design..."
for f in "$OUT/rfdiffusion"/il7ra_binder_*.pdb; do cp "$f" "$OUT/proteinmpnn/pdbs/"; done
cd "$MPNN"

python helper_scripts/parse_multiple_chains.py \
  --input_path "$OUT/proteinmpnn/pdbs" --output_path "$OUT/proteinmpnn/parsed.jsonl" 2>&1 | tail -3

python helper_scripts/assign_fixed_chains.py \
  --input_path "$OUT/proteinmpnn/parsed.jsonl" --output_path "$OUT/proteinmpnn/assigned.jsonl" \
  --chain_list "A" 2>&1 | tail -1

python protein_mpnn_run.py \
  --jsonl_path "$OUT/proteinmpnn/parsed.jsonl" \
  --chain_id_jsonl "$OUT/proteinmpnn/assigned.jsonl" \
  --out_folder "$OUT/proteinmpnn/output" \
  --num_seq_per_target 1 --sampling_temp "0.1" --seed 37 --batch_size 1 \
  --path_to_model_weights "${MPNN}/vanilla_model_weights" 2>&1 | tail -5

MP_COUNT=$(ls "$OUT/proteinmpnn/output/seqs" 2>/dev/null | wc -l)
echo "ProteinMPNN: $MP_COUNT sequences generated"

# ---- Stage 3: Build benchmark CSV from real designed sequences ----
echo "[$(date)] Stage 3: Building benchmark CSV from real designs..."
conda deactivate 2>/dev/null || true
export CONDA_DEFAULT_ENV=colabfold
export CONDA_PREFIX=/root/miniconda3/envs/colabfold
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/colabfold/lib:${LD_LIBRARY_PATH}

cd /mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
python3 << PYEOF
import glob, re
import pandas as pd

target_seq = "$TARGET_SEQ"
rows = []
for fa in sorted(glob.glob("$OUT/proteinmpnn/output/seqs/*.fa")):
    with open(fa) as f:
        lines = f.readlines()
    # ProteinMPNN .fa: line 0 = original header, line 1 = original seq, line 2 = designed header, line 3 = designed seq
    binder_seq = lines[3].strip()
    rows.append({
        "sequence": binder_seq,
        "target_seq": target_seq,
        "target_length": len(target_seq),
        "sequence_length": len(binder_seq),
        "outcome": "unknown",
        "source_system": "il7ra_real_design",
    })

df = pd.DataFrame(rows)
df.to_csv("benchmark/il7ra_real_benchmark.csv", index=False)
print(f"[OK] Real IL-7Ra benchmark: {len(df)} designed sequences")
PYEOF

echo ""
echo "======================================================================"
echo "STAGE 1-3 COMPLETE. Next: fast-screen ColabFold"
echo "======================================================================"
echo "Finish: $(date)"
