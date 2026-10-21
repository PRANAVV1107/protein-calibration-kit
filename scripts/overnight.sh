#!/bin/bash
# Two jobs, run in order, fully checkpointed and pausable:
#   1. rog prospective test   80 folds  ~4.5 h   (registers its prediction first)
#   2. instrument calibration 60 folds  ~3.2 h
#
# CHECKPOINTING
# -------------
# Every fold lands as its own JSON under <job>/folds/. On restart each job
# rebuilds its to-do list by diffing the fasta against what is already on disk,
# so nothing is ever re-folded and no state file can go stale. The unit of loss
# from a hard kill is one chunk (default 5 designs, ~20 min).
#
# PAUSING
# -------
#   bash scripts/pause_jobs.sh     graceful: finishes the current chunk, exits
#   bash scripts/resume_jobs.sh    clears the flag and restarts where it stopped
#
# The flag is checked between chunks, so a graceful pause costs at most one
# chunk of waiting and loses nothing. A hard kill (Ctrl-C, reboot) is also safe
# -- it only discards the chunk in flight, which is redone on resume.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
LOG=$KIT/results/overnight.log
PAUSE=$KIT/results/.PAUSE
CHUNK=${CHUNK:-5}

unset PYTHONPATH                 # /root/rf_deps pandas breaks ColabFold's 3.10
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

say(){ echo "[$(date '+%F %T')] $*"; }

# Every alternative bracketed on its first character so this script's own argv,
# which ps prints, cannot match it. Unbracketed patterns have killed the wrong
# process five times in this project.
busy=$(ps -eo args --width 400 | grep -cE '[e]nvs/colabfold/bin/colabfold_batch|[e]nvs/colabfold/bin/python .*scripts/mcts_design\.py|[m]d_stability\.py')
if [ "$busy" -gt 0 ]; then echo "[abort] GPU busy ($busy processes)"; exit 1; fi

cd "$KIT" || exit 1
mkdir -p results
{
say "=== OVERNIGHT: rog prospective test, then instrument calibration ==="
[ -f "$PAUSE" ] && { say "PAUSE flag present at start -- clear it with scripts/resume_jobs.sh"; exit 0; }

run_job () {          # $1 = dir, $2 = label
  local OUT="$KIT/$1" LABEL="$2"
  mkdir -p "$OUT/folds"
  local total done_n
  total=$(grep -c '^>' "$OUT/to_fold.fasta")
  done_n=$(ls "$OUT"/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l)
  say "--- $LABEL: $done_n/$total already folded ---"
  [ "$done_n" -ge "$total" ] && { say "  $LABEL already complete"; return 0; }

  /root/miniconda3/envs/colabfold/bin/python - "$OUT" "$CHUNK" "$PAUSE" "$LABEL" <<'PYEOF'
import os, sys, glob, subprocess, time
out, chunk, pause, label = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
recs, name = [], None
for line in open(os.path.join(out, "to_fold.fasta")):
    line = line.strip()
    if line.startswith(">"): name = line[1:]
    elif name: recs.append((name, line)); name = None

def remaining():
    done = {os.path.basename(f).split("_scores_rank_001_")[0]
            for f in glob.glob(os.path.join(out, "folds", "*_scores_rank_001_*.json"))}
    return [r for r in recs if r[0] not in done]

t0, n0 = time.time(), len(recs) - len(remaining())
while True:
    todo = remaining()
    if not todo:
        print("[%s] complete: %d/%d" % (label, len(recs), len(recs)), flush=True); break
    if os.path.exists(pause):
        print("[%s] PAUSED at %d/%d -- resume with scripts/resume_jobs.sh"
              % (label, len(recs) - len(todo), len(recs)), flush=True)
        sys.exit(42)
    part = todo[:chunk]
    p = os.path.join(out, ".chunk.fasta")
    with open(p, "w") as fh:
        for nm, sq in part: fh.write(">%s\n%s\n" % (nm, sq))
    subprocess.run(["colabfold_batch", "--msa-mode", "mmseqs2_uniref_env",
                    "--model-type", "alphafold2_multimer_v3",
                    "--num-models", "1", "--num-recycle", "3",
                    p, os.path.join(out, "folds")], check=False)
    n = len(recs) - len(remaining())
    el = time.time() - t0
    rate = el / max(1, n - n0)
    print("[%s] %d/%d folded, %.1f min elapsed, %.0f s/design, ~%.1f h left"
          % (label, n, len(recs), el / 60, rate, rate * (len(recs) - n) / 3600), flush=True)
if os.path.exists(os.path.join(out, ".chunk.fasta")):
    os.remove(os.path.join(out, ".chunk.fasta"))
PYEOF
  return $?
}

# ---------- job 1: rog prospective ----------
# Every batch's ProteinMPNN output holds ONE sequence per backbone. The first
# attempt at this experiment took what existed and silently ran 16 backbones at
# K=1 against a registration specifying 20 at K=4. Re-sample first; the prepare
# script now refuses to proceed if the protocol is not satisfiable.
ROGOUT=$KIT/results/rog_prospective
MPNN=/mnt/c/Users/prana/Downloads/proteinfoldingexp/external_tools/ProteinMPNN
if [ ! -d "$ROGOUT/mpnn/seqs" ]; then
  say "sampling 4 ProteinMPNN draws per backbone for batches 1-3"
  mkdir -p "$ROGOUT/mpnn"
  echo '{"A": -1.0}' > "$ROGOUT/bias_A1.jsonl"
  for b in 1 2 3; do
    SRC=/mnt/c/Users/prana/Downloads/proteinfoldingexp/real_test_results/rbd_200_batch/batch_$b/proteinmpnn
    ( cd "$MPNN" && LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH:-} /root/miniconda3/envs/rfdiffusion/bin/python protein_mpnn_run.py --jsonl_path "$SRC/parsed.jsonl" --chain_id_jsonl "$SRC/assigned.jsonl" --out_folder "$ROGOUT/mpnn_b$b" --num_seq_per_target 4 --sampling_temp 0.2 --seed 37 --batch_size 1 --bias_AA_jsonl "$ROGOUT/bias_A1.jsonl" ) || say "  batch $b MPNN returned nonzero"
  done
  # Backbone names REPEAT across batches -- binder_0 exists in all three, and
  # the selected compact group contains both b2/binder_0 and b3/binder_0. Merging
  # on bare names would silently fold one batch's sequences under the other's
  # identity. Key every file by batch AND backbone.
  mkdir -p "$ROGOUT/mpnn/seqs"
  for b in 1 2 3; do
    for f in "$ROGOUT/mpnn_b$b/seqs"/*.fa; do
      [ -f "$f" ] && cp "$f" "$ROGOUT/mpnn/seqs/b${b}_$(basename "$f")"
    done
  done
  say "  $(ls "$ROGOUT/mpnn/seqs"/*.fa 2>/dev/null | wc -l) fasta files merged"
fi

if [ ! -f "$KIT/results/rog_prospective/to_fold.fasta" ]; then
  say "registering the rog prediction BEFORE any fold exists"
  /root/miniconda3/envs/colabfold/bin/python scripts/rog_prospective_prepare.py \
    --out results/rog_prospective || { say "rog prepare FAILED"; exit 1; }
fi
run_job results/rog_prospective "rog"
rc=$?
[ "$rc" -eq 42 ] && { say "paused during rog test"; exit 0; }

# ---------- job 2: instrument calibration ----------
if [ ! -f "$KIT/results/calibration/to_fold.fasta" ]; then
  say "building calibration fold set"
  /root/miniconda3/envs/colabfold/bin/python scripts/calibration_prepare.py \
    --out results/calibration || { say "calibration prepare FAILED"; exit 1; }
fi
run_job results/calibration "calib"
rc=$?
[ "$rc" -eq 42 ] && { say "paused during calibration"; exit 0; }

say "=== BOTH JOBS COMPLETE ==="
say "  rog:   $(ls results/rog_prospective/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l) folds"
say "  calib: $(ls results/calibration/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l) folds"
} >> "$LOG" 2>&1
