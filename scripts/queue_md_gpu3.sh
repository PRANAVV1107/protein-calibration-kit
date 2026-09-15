#!/bin/bash
# GPU molecular dynamics over the stratified panel. Replaces queue_md_gpu2.sh.
#
# What changed and why:
#
#   panel   queue_md_gpu2 ran 15 il7ra sequence-refinement complexes -- one
#           target, and selected for quality. Both are wrong for a correlation.
#           Single-target generalises to nothing, and selecting good designs
#           truncates the predictor's range, which drives a correlation toward
#           zero however good the metric is. This reads the stratified panel
#           from select_md_panel.py: 37 complexes across il7ra, PD-L1 and RBD,
#           spanning the full ipTM range 0.13-0.92, deliberately including
#           complexes expected to fall apart.
#
#   length  2 ns rather than 5. For ranking metrics by correlation, spread
#           across designs buys more than trajectory length per design: at
#           n=15 a Spearman rho needs |rho| > 0.51 for p < 0.05, at n=37 about
#           0.32. Same wall time either way, near 28 h at a measured ~3 ns/h.
#
# The guard matches the interpreter path followed by the script path, a
# two-token sequence only the real ColabFold process produces. Bare
# `pgrep -f <fragment>` self-matches watcher shells and has deadlocked this
# project three times.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/md_gpu3.log
PY=/root/miniconda3/envs/openmm-gpu/bin/python
PANEL=$KIT/results/md_panel_paths.txt

NS=${NS:-2}
EQUIL=${EQUIL:-200}

cd "$KIT" || exit 1
mkdir -p results/md_gpu

say(){ echo "[$(date '+%F %T')] $*"; }
folding(){
  ps -eo args --width 400 2>/dev/null \
    | grep -c 'envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch'
}

{
say "=== GPU MD v3 START (panel, ${NS} ns) ==="
if [ ! -s "$PANEL" ]; then
  say "FATAL: $PANEL missing or empty -- run scripts/select_md_panel.py first"
  exit 1
fi
mapfile -t PDBS < "$PANEL"
say "${#PDBS[@]} complexes queued"

say "waiting for folding to finish (seen: $(folding))"
while [ "$(folding)" -gt 0 ]; do sleep 60; done
say "GPU free."

i=0
for p in "${PDBS[@]}"; do
  i=$((i+1))
  [ -f "$p" ] || { say "  MISSING $p"; continue; }
  b=$(basename "$p" .pdb)
  if [ -f "results/md_gpu/${b}_md.json" ]; then say "  skip  $b"; continue; fi
  # a fold job can start between complexes; re-check before each
  while [ "$(folding)" -gt 0 ]; do sleep 60; done
  say "  START [$i/${#PDBS[@]}] $b"
  "$PY" scripts/md_stability.py "$p" --ns "$NS" --equil-ps "$EQUIL" \
        --platform CUDA --out-dir results/md_gpu >> "$LOG" 2>&1 \
     || say "  FAILED $b"
done

say "=== GPU MD v3 COMPLETE: $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l) systems ==="
} >> "$LOG" 2>&1
