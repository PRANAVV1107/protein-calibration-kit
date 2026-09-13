#!/bin/bash
# Empirically sweep ProteinMPNN sampling settings against the existing IL-7Ra
# backbones, to pick a temperature / alanine-bias combination that avoids
# low-complexity (poly-Ala) collapse without degrading design quality.
set -e

export PATH=/root/miniconda3/bin:$PATH
MPNN=/mnt/c/Users/prana/Downloads/proteinfoldingexp/external_tools/ProteinMPNN
PILOT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit/results/il7ra_real_pilot50
SWEEP=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit/results/mpnn_sweep

export CONDA_DEFAULT_ENV=rfdiffusion
export CONDA_PREFIX=/root/miniconda3/envs/rfdiffusion
export PATH=/root/miniconda3/envs/rfdiffusion/bin:$PATH
export LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH}

mkdir -p "$SWEEP"
cd "$MPNN"

# Bias files: negative logit bias makes that residue less likely.
echo '{"A": -1.0}'            > "$SWEEP/bias_A1.jsonl"
echo '{"A": -2.0}'            > "$SWEEP/bias_A2.jsonl"
echo '{"A": -1.0, "G": -0.5}' > "$SWEEP/bias_A1G.jsonl"

run_cfg () {
  local name="$1" temp="$2" bias="$3"
  local out="$SWEEP/$name"
  mkdir -p "$out"
  local bias_arg=""
  if [ -n "$bias" ]; then bias_arg="--bias_AA_jsonl $bias"; fi

  python protein_mpnn_run.py \
    --jsonl_path "$PILOT/proteinmpnn/parsed.jsonl" \
    --chain_id_jsonl "$PILOT/proteinmpnn/assigned.jsonl" \
    --out_folder "$out" \
    --num_seq_per_target 1 --sampling_temp "$temp" --seed 37 --batch_size 1 \
    --path_to_model_weights "${MPNN}/vanilla_model_weights" \
    $bias_arg > "$out/run.log" 2>&1
  echo "  done: $name (temp=$temp bias=${bias:-none})"
}

echo "Running ProteinMPNN sweep over 50 IL-7Ra backbones..."
run_cfg baseline_t0.1      0.1  ""
run_cfg t0.2               0.2  ""
run_cfg t0.3               0.3  ""
run_cfg t0.1_biasA1        0.1  "$SWEEP/bias_A1.jsonl"
run_cfg t0.2_biasA1        0.2  "$SWEEP/bias_A1.jsonl"
run_cfg t0.2_biasA2        0.2  "$SWEEP/bias_A2.jsonl"
run_cfg t0.2_biasA1G       0.2  "$SWEEP/bias_A1G.jsonl"

echo "Sweep complete. Results in $SWEEP"
