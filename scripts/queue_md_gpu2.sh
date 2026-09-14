#!/bin/bash
# GPU molecular dynamics, queued behind all folding. Replaces queue_md_gpu.sh.
#
# Two bugs in the previous version, both in the wait guards:
#
#   1. `pgrep -f 'queue_v8.sh'` matched the stale launcher shell, whose command
#      line contains that filename as a Python string literal. queue_v8 had
#      been dead for hours; the guard never released. Third time a bare pgrep
#      fragment has deadlocked this project. This version does not wait on
#      queue_v8 at all -- it is gone.
#
#   2. `pgrep -f 'colabfold_batch --model'` was written to match the binary
#      plus an argument, but the high-accuracy job puts --msa-mode first, so
#      the pattern missed the very job MD must not run alongside. Matching on
#      an argument is fragile because arguments get reordered. This version
#      matches the interpreter path immediately followed by the script path,
#      a two-token sequence only the real process has -- no watcher shell or
#      admin command produces it, and it does not depend on argument order.
#
# Waits for every folding job first: MD and ColabFold contend for 6GB of VRAM,
# and concurrent GPU work has OOM'd here before.
set -u

KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/md_gpu2.log
PY=/root/miniconda3/envs/openmm-gpu/bin/python

NS=${NS:-50}
EQUIL=${EQUIL:-200}

cd "$KIT"
mkdir -p results/md_gpu

say(){ echo "[$(date '+%F %T')] $*"; }

folding(){
  ps -eo args --width 400 2>/dev/null \
    | grep -c 'envs/colabfold/bin/python3\.10 /root/miniconda3/envs/colabfold/bin/colabfold_batch'
}

{
say "=== GPU MD v2 START ==="
n=$(folding)
say "folding processes seen: $n"
while [ "$(folding)" -gt 0 ]; do sleep 60; done
say "GPU free."

mapfile -t PDBS < <(ls results/il7ra_seqrefine/folds/*_unrelaxed_rank_001_*.pdb 2>/dev/null)
say "${#PDBS[@]} complexes, ${NS} ns each on CUDA"

for p in "${PDBS[@]}"; do
  b=$(basename "$p" .pdb)
  if [ -f "results/md_gpu/${b}_md.json" ]; then say "  skip $b"; continue; fi
  # a fold job can start between complexes; re-check before each one
  while [ "$(folding)" -gt 0 ]; do sleep 60; done
  say "  START $b"
  "$PY" scripts/md_stability.py "$p" --ns "$NS" --equil-ps "$EQUIL" \
        --platform CUDA --out-dir results/md_gpu >> "$LOG" 2>&1 \
     || say "  FAILED $b"
done

say "=== GPU MD COMPLETE: $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l) systems ==="
} >> "$LOG" 2>&1
