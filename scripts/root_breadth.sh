#!/bin/bash
# Root-breadth sweep: 50 RBD backbones x 4 draws, full MSA at 1 model / 3
# recycles -- the configuration section 3.2 recommends.
#
# Tests the tail-statistic prediction that E[best] climbs with the NUMBER of
# roots faster than with draws per root. Ten roots could not test it; the
# allocation table flattens at 0.869 only because there are ten distributions.
#
# Ten of the fifty are already measured at 3 models / 8 recycles, kept in as an
# overlap set. Running a new configuration without one is the mistake that
# invalidated three search runs.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
OUT=$KIT/results/root_breadth
LOG=$KIT/results/root_breadth.log
CHUNK=${CHUNK:-10}

unset PYTHONPATH                 # /root/rf_deps pandas breaks ColabFold's 3.10
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=2.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

say(){ echo "[$(date '+%F %T')] $*"; }

# Every alternative's first character is bracketed so the grep's own argv,
# which ps prints, cannot match. Unbracketed patterns have self-matched and
# killed the wrong thing five times in this project.
n=$(ps -eo args --width 400 | grep -cE '[e]nvs/colabfold/bin/python .*scripts/mcts_design\.py|[e]nvs/colabfold/bin/colabfold_batch|[m]d_stability\.py')
if [ "$n" -gt 0 ]; then echo "[abort] GPU busy ($n processes)"; exit 1; fi

cd "$KIT" || exit 1
{
say "=== ROOT-BREADTH SWEEP ==="
# The campaign's ProteinMPNN run produced a single sequence per backbone, which
# cannot estimate a per-root distribution. Re-sample 4 draws each. Settings match
# the rest of the project: temp 0.2, alanine bias -1.0.
MPNN=/mnt/c/Users/prana/Downloads/proteinfoldingexp/external_tools/ProteinMPNN
SRC=/mnt/c/Users/prana/Downloads/proteinfoldingexp/real_test_results/rbd_200_batch/batch_0/proteinmpnn
if [ ! -d "$OUT/mpnn/seqs" ]; then
  say "sampling 4 ProteinMPNN draws per backbone"
  mkdir -p "$OUT/mpnn"
  echo '{"A": -1.0}' > "$OUT/bias_A1.jsonl"
  ( cd "$MPNN" && LD_LIBRARY_PATH=/root/miniconda3/envs/rfdiffusion/lib:${LD_LIBRARY_PATH:-}     /root/miniconda3/envs/rfdiffusion/bin/python protein_mpnn_run.py       --jsonl_path "$SRC/parsed.jsonl" --chain_id_jsonl "$SRC/assigned.jsonl"       --out_folder "$OUT/mpnn" --num_seq_per_target 4 --sampling_temp 0.2       --seed 37 --batch_size 1 --bias_AA_jsonl "$OUT/bias_A1.jsonl" )     || say "  MPNN returned nonzero"
  say "  $(ls "$OUT/mpnn/seqs"/*.fa 2>/dev/null | wc -l) fasta files"
fi

if [ ! -f "$OUT/to_fold.fasta" ]; then
  say "building fold set"
  /root/miniconda3/envs/colabfold/bin/python scripts/root_breadth_prepare.py --out "$OUT" \
    || { say "prepare FAILED"; exit 1; }
fi
say "$(grep -c '^>' "$OUT/to_fold.fasta") designs queued"

mkdir -p "$OUT/folds"
# Chunked and resumable: a single colabfold_batch over many inputs has been
# OOM-killed three times in this project, most recently at 85/165.
/root/miniconda3/envs/colabfold/bin/python - "$OUT" "$CHUNK" <<'PYEOF'
import os, sys, glob, subprocess, time
out, chunk = sys.argv[1], int(sys.argv[2])
recs, name = [], None
for line in open(os.path.join(out, "to_fold.fasta")):
    line = line.strip()
    if line.startswith(">"): name = line[1:]
    elif name: recs.append((name, line)); name = None
done = {os.path.basename(f).split("_scores_rank_001_")[0]
        for f in glob.glob(os.path.join(out, "folds", "*_scores_rank_001_*.json"))}
todo = [r for r in recs if r[0] not in done]
print("[breadth] %d done, %d to fold" % (len(done), len(todo)), flush=True)
t0 = time.time()
for i in range(0, len(todo), chunk):
    part = todo[i:i+chunk]
    p = os.path.join(out, ".chunk.fasta")
    with open(p, "w") as fh:
        for nm, sq in part: fh.write(">%s\n%s\n" % (nm, sq))
    subprocess.run(["colabfold_batch", "--msa-mode", "mmseqs2_uniref_env",
                    "--model-type", "alphafold2_multimer_v3",
                    "--num-models", "1", "--num-recycle", "3",
                    p, os.path.join(out, "folds")], check=False)
    n = len(glob.glob(os.path.join(out, "folds", "*_scores_rank_001_*.json")))
    el = time.time() - t0
    rate = el / max(1, n - len(done))
    print("[breadth] %d folded, %.1f min elapsed, %.1f s/design, ~%.1f h left"
          % (n, el/60, rate, rate * max(0, len(recs) - n) / 3600), flush=True)
if os.path.exists(os.path.join(out, ".chunk.fasta")):
    os.remove(os.path.join(out, ".chunk.fasta"))
PYEOF

say "folded: $(ls "$OUT"/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l)"
/root/miniconda3/envs/colabfold/bin/python scripts/root_breadth_analyze.py --dir "$OUT" \
  > "$KIT/results/ROOT_BREADTH_RESULT.txt" 2>&1 \
  && say "analysed -> results/ROOT_BREADTH_RESULT.txt" || say "analysis FAILED"
say "=== SWEEP COMPLETE ==="
} >> "$LOG" 2>&1
