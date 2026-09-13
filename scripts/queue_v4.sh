#!/bin/bash
# Remaining work, ordered cheapest-and-most-decisive first.
#
#   1. template_rbd    60 folds  ~45 m  completes the 4-way screen comparison
#                                       (retry: hhsearch was missing, now installed)
#   2. decoys_full     30 folds  ~2.2 h raises the specificity panel from n=5 to n=15,
#                                       enough to test whether interface properties
#                                       predict specificity (n=5 could not)
#   3. il7ra_rescreen 114 folds  ~8.5 h the 38 IL-7Ra backbones never evaluated,
#                                       because the fast screen that cut them does
#                                       not rank -- on this target its bottom decile
#                                       held pair_033 at 0.870, the joint best design
#
# Sequential: concurrent ColabFold on a 6GB card has OOM'd in this project.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
LOG=$KIT/results/queue_v4.log

export PATH=/root/miniconda3/envs/colabfold/bin:/root/miniconda3/bin:/usr/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
export XLA_PYTHON_CLIENT_ALLOCATOR=platform

REAL_FOLD_RE='envs/colabfold/bin/colabfold_batch'
say(){ echo "[$(date '+%F %T')] $*"; }
count(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

{
say "=== QUEUE V4 START ==="
if pgrep -f "$REAL_FOLD_RE" > /dev/null; then
  say "waiting for in-flight fold..."
  while pgrep -f "$REAL_FOLD_RE" > /dev/null; do sleep 30; done
fi
say "hhsearch: $(command -v hhsearch || echo MISSING)"

# ---- 1. template screen (target structure instead of an alignment) ----
mkdir -p "$Q/template_rbd_folds" "$Q/template_cache"
say "START template_rbd (60)"
colabfold_batch --msa-mode single_sequence --model-type alphafold2_multimer_v3 \
  --num-models 1 --num-recycle 3 \
  --templates --custom-template-path "$Q/templates_rbd" \
  --custom-template-cache-path "$Q/template_cache" \
  "$Q/rewardchannel_rbd.fasta" "$Q/template_rbd_folds" > "$Q/template_rbd.log" 2>&1 \
  && say "DONE  template_rbd ($(count "$Q/template_rbd_folds"))" \
  || say "FAIL  template_rbd ($(count "$Q/template_rbd_folds") done; see log)"

# ---- 2. full decoy panel: 15 candidates x 2 decoys ----
mkdir -p "$Q/il7ra_decoys_full_folds"
say "START decoys_full (30)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 3 --num-recycle 8 \
  "$Q/il7ra_decoys_full.fasta" "$Q/il7ra_decoys_full_folds" > "$Q/decoys_full.log" 2>&1 \
  && say "DONE  decoys_full ($(count "$Q/il7ra_decoys_full_folds"))" \
  || say "FAIL  decoys_full ($(count "$Q/il7ra_decoys_full_folds") done; see log)"

# ---- 3. the 38 backbones the broken fast screen discarded ----
mkdir -p "$Q/il7ra_rescreen_folds"
say "START il7ra_rescreen (114)"
colabfold_batch --msa-mode mmseqs2_uniref_env --model-type alphafold2_multimer_v3 \
  --num-models 3 --num-recycle 8 \
  "$Q/il7ra_rescreen.fasta" "$Q/il7ra_rescreen_folds" > "$Q/il7ra_rescreen.log" 2>&1 \
  && say "DONE  il7ra_rescreen ($(count "$Q/il7ra_rescreen_folds"))" \
  || say "FAIL  il7ra_rescreen ($(count "$Q/il7ra_rescreen_folds") done; see log)"

say "=== QUEUE V4 COMPLETE ==="
} >> "$LOG" 2>&1
