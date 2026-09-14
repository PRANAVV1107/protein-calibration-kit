#!/bin/bash
# MD stability screen on the CPU platform, run IN PARALLEL with GPU folding.
#
# The two jobs use different resources: ColabFold saturates the GPU and one CPU
# core, while this machine has 20 cores sitting at load ~5. MD on the CPU
# platform therefore costs nothing that folding needs.
#
# Parameters are set for a gross-stability SCREEN, not production simulation:
# 1 ns with 20 ps equilibration is enough to see an interface come apart, and
# designed helical binders that dissociate generally do so quickly. State that
# length plainly in any writeup -- it rejects unstable poses, it does not
# confirm stable ones.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/md_screen.log
PY=/root/miniconda3/envs/colabfold/bin/python

NS=${NS:-1.0}
EQUIL=${EQUIL:-20}
THREADS=${THREADS:-7}      # 2 concurrent x 7 = 14 of 20 cores, leaves room for ColabFold

cd "$KIT"
mkdir -p results/md

# candidates: the IL-7Ra finalists, best structure per design
mapfile -t PDBS < <(ls results/il7ra_seqrefine/folds/*_unrelaxed_rank_001_*.pdb 2>/dev/null)

{
echo "[$(date '+%F %T')] MD screen: ${#PDBS[@]} complexes, ${NS} ns each, ${THREADS} threads x2 concurrent"

i=0
for p in "${PDBS[@]}"; do
  b=$(basename "$p" .pdb)
  if [ -f "results/md/${b}_md.json" ]; then
    echo "  skip $b (already done)"; continue
  fi
  "$PY" scripts/md_stability.py "$p" --ns "$NS" --equil-ps "$EQUIL" \
        --threads "$THREADS" --out-dir results/md &
  i=$((i+1))
  # keep at most 2 running so ColabFold is never starved of CPU
  if [ $((i % 2)) -eq 0 ]; then wait; fi
done
wait

echo "[$(date '+%F %T')] MD screen complete: $(ls results/md/*_md.json 2>/dev/null | wc -l) systems"
} >> "$LOG" 2>&1
