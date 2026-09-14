#!/bin/bash
# Reordered by impact. The reprice: this project now holds 325 experimentally
# verified NON-binders (Adaptyv EGFR round 2), which changes what the queue is
# worth.
#
#  1  adaptyv seed        1 fold    generates the EGFR alignment
#  2  adaptyv validation  365 folds THE only measured outcomes in this project.
#                                   Everything else validates AlphaFold against
#                                   AlphaFold. This asks whether the selection
#                                   criterion predicts real binding calls, and
#                                   supplies real labels for the calibration
#                                   curve, which is currently fit on synthetic
#                                   ones and disclaimed in the README.
#  3  family panel        75 folds  homolog decoys (gammaC, IL4Ra, IL21R) --
#                                   the cross-reactivity that would actually
#                                   sink a therapeutic. Cheap.
#  4  il7ra rescreen      114 folds the 38 backbones the broken fast screen cut.
#                                   Would likely find a better candidate, but
#                                   that candidate is unvalidatable without a
#                                   wet lab. Demoted.
#  5  null panel          420 folds constructed negatives. Was the right idea
#                                   when measured negatives were unavailable.
#                                   325 real ones now supersede it; retains
#                                   second-order value for cross-target
#                                   generality. Demoted hardest.
#
# Jobs 2-5 make no MSA-server calls: every target alignment is precomputed once
# and reused. Sequential throughout -- concurrent ColabFold has OOM'd on this
# 6GB card.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v7.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

{
say "=== QUEUE V7 START (reordered by impact) ==="
if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
  say "waiting for in-flight fold..."
  while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 20; done
fi
say "GPU free."

# ---- 1. seed: one high-accuracy fold to obtain the EGFR alignment ----
mkdir -p "$Q/adaptyv_seed"
say "START adaptyv_seed (1 fold, generates EGFR domain III alignment)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3 \
  "$Q/adaptyv_seed.fasta" "$Q/adaptyv_seed" > "$Q/adaptyv_seed.log" 2>&1 \
  && say "DONE  adaptyv_seed" || say "FAIL  adaptyv_seed (see log)"

# ---- 2. build the 365 hybrid inputs from that alignment ----
cd "$KIT"
SEED=$(ls "$Q/adaptyv_seed"/*.a3m 2>/dev/null | head -1)
if [ -n "$SEED" ]; then
  say "building validation inputs from $(basename "$SEED")"
  python scripts/build_adaptyv_validation.py --seed-a3m "$SEED" >> "$LOG" 2>&1 \
    || say "FAIL  build_adaptyv_validation"
else
  say "FAIL  no alignment produced; skipping validation build"
fi

# ---- 3. the validation run itself ----
if [ -d "$Q/adaptyv_val" ]; then
  mkdir -p "$Q/adaptyv_val_folds"
  say "START adaptyv_validation (365 designs, 55 binders / 310 non-binders)"
  colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
    "$Q/adaptyv_val" "$Q/adaptyv_val_folds" > "$Q/adaptyv_val.log" 2>&1 \
    && say "DONE  adaptyv_validation ($(count "$Q/adaptyv_val_folds"))" \
    || say "FAIL  adaptyv_validation ($(count "$Q/adaptyv_val_folds") done)"
fi

# ---- 4. family decoys (seed + panel) ----
mkdir -p "$Q/family_seed"
say "START family_seed (5 folds, one alignment per decoy)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3 \
  "$Q/family_seed.fasta" "$Q/family_seed" > "$Q/family_seed.log" 2>&1 \
  && say "DONE  family_seed" || say "FAIL  family_seed"
python scripts/build_family_panel.py >> "$LOG" 2>&1 || say "FAIL  build_family_panel"
if [ -d "$Q/family_panel" ]; then
  mkdir -p "$Q/family_panel_folds"
  say "START family_panel"
  colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
    "$Q/family_panel" "$Q/family_panel_folds" > "$Q/family_panel.log" 2>&1 \
    && say "DONE  family_panel ($(count "$Q/family_panel_folds"))" || say "FAIL  family_panel"
fi

# ---- 5. demoted: rescreen (resumes from the 13 already done) ----
mkdir -p "$Q/il7ra_rescreen_folds"
say "START il7ra_rescreen (114, resumes)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 3 --num-recycle 8 \
  "$Q/il7ra_rescreen.fasta" "$Q/il7ra_rescreen_folds" > "$Q/il7ra_rescreen.log" 2>&1 \
  && say "DONE  il7ra_rescreen ($(count "$Q/il7ra_rescreen_folds"))" || say "FAIL  il7ra_rescreen"

# ---- 6. demoted hardest: constructed null ----
mkdir -p "$Q/null_panel_folds"
say "START null_panel (420 constructed negatives)"
colabfold_batch --model-type alphafold2_multimer_v3 --num-models 1 --num-recycle 3 \
  "$Q/null_panel" "$Q/null_panel_folds" > "$Q/null_panel.log" 2>&1 \
  && say "DONE  null_panel ($(count "$Q/null_panel_folds"))" || say "FAIL  null_panel"

say "=== QUEUE V7 COMPLETE ==="
} >> "$LOG" 2>&1
