#!/bin/bash
# GPU molecular dynamics, queued behind all folding.
#
# The CPU attempt was abandoned: over 90 minutes per complex for 1 ns, and 15
# complexes to do. OpenMM rates its CUDA platform at speed 100 against the CPU
# platform's 1, so the same work is minutes rather than hours -- which changes
# what the simulation can claim. 1 ns only rejects a pose that falls apart
# immediately. 50 ns starts to say something about stability.
#
# Runs from the openmm-gpu env, built separately so the colabfold env's JAX/CUDA
# stack is untouched. That env's OpenMM has no CUDA platform at all, which is
# also why --use-gpu-relax fails there.
#
# Waits for every folding job first: MD and ColabFold would otherwise contend
# for 6GB of VRAM, and concurrent GPU work has OOM'd in this project before.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/md_gpu.log
PY=/root/miniconda3/envs/openmm-gpu/bin/python

NS=${NS:-50}
EQUIL=${EQUIL:-200}

cd "$KIT"
mkdir -p results/md_gpu

say(){ echo "[$(date '+%F %T')] $*"; }

{
say "=== GPU MD START -- waiting for folding to finish ==="
while pgrep -f 'queue_v8.sh' > /dev/null; do sleep 120; done
say "queue_v8 finished."
# match the interpreter plus its arguments, not a bare path fragment: any
# process whose command line merely MENTIONS the fragment would match, which
# has deadlocked this project twice.
while pgrep -f 'colabfold_batch --model' > /dev/null; do sleep 30; done
say "GPU free."

mapfile -t PDBS < <(ls results/il7ra_seqrefine/folds/*_unrelaxed_rank_001_*.pdb 2>/dev/null)
say "${#PDBS[@]} complexes, ${NS} ns each on CUDA"

for p in "${PDBS[@]}"; do
  b=$(basename "$p" .pdb)
  if [ -f "results/md_gpu/${b}_md.json" ]; then say "  skip $b"; continue; fi
  say "  START $b"
  "$PY" scripts/md_stability.py "$p" --ns "$NS" --equil-ps "$EQUIL" \
        --platform CUDA --out-dir results/md_gpu >> "$LOG" 2>&1 \
     || say "  FAILED $b"
done

say "=== GPU MD COMPLETE: $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l) systems ==="
} >> "$LOG" 2>&1
