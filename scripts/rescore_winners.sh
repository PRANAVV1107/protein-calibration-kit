#!/bin/bash
# Re-score search designs at HIGH ACCURACY, the same instrument the 0.870 bar
# was measured on. See scripts/rescore_prepare.py for why this is necessary.
#
# Recipe copied deliberately from variance_study.sh stage 3 -- full MSA via
# mmseqs2, 3 models, 8 recycles. Any deviation reintroduces the instrument
# mismatch this experiment exists to remove.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
OUT=$KIT/results/rescore
LOG=$KIT/results/rescore.log

unset PYTHONPATH                 # /root/rf_deps pandas breaks ColabFold's 3.10
export PATH=/root/miniconda3/envs/colabfold/bin:$PATH
export TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0 XLA_PYTHON_CLIENT_ALLOCATOR=platform

say(){ echo "[$(date '+%F %T')] $*"; }

# Refuse to run alongside anything else on the card. Matches interpreter path
# AND script path -- bare -f fragments have self-matched and deadlocked this
# queue three times.
n=$(ps -eo args --width 400 | grep -cE '[e]nvs/colabfold/bin/python .*scripts/mcts_design\.py|[e]nvs/colabfold/bin/colabfold_batch|[m]d_stability\.py')
if [ "$n" -gt 0 ]; then echo "[abort] GPU busy ($n processes)"; exit 1; fi

cd "$KIT" || exit 1
{
say "=== RESCORE SEARCH WINNERS AT HIGH ACCURACY ==="
if [ ! -f "$OUT/to_fold.fasta" ]; then
  say "building fold set"
  /root/miniconda3/envs/colabfold/bin/python scripts/rescore_prepare.py --out "$OUT" \
    || { say "prepare FAILED"; exit 1; }
fi
n_in=$(grep -c '^>' "$OUT/to_fold.fasta")
say "$n_in designs to re-score"

mkdir -p "$OUT/folds"
# Chunked: a single colabfold_batch over many inputs has been OOM-killed three
# times in this project, most recently at 85/165.
CHUNK=${CHUNK:-10}
/root/miniconda3/envs/colabfold/bin/python - "$OUT" "$CHUNK" <<'PYEOF'
import os, sys, glob, subprocess
out, chunk = sys.argv[1], int(sys.argv[2])
recs, name = [], None
for line in open(os.path.join(out, "to_fold.fasta")):
    line = line.strip()
    if line.startswith(">"): name = line[1:]
    elif name: recs.append((name, line)); name = None
done = {os.path.basename(f).split("_scores_rank_001_")[0]
        for f in glob.glob(os.path.join(out, "folds", "*_scores_rank_001_*.json"))}
todo = [r for r in recs if r[0] not in done]
print("[rescore] %d done, %d to fold" % (len(done), len(todo)), flush=True)
for i in range(0, len(todo), chunk):
    part = todo[i:i+chunk]
    p = os.path.join(out, ".chunk.fasta")
    with open(p, "w") as fh:
        for nm, sq in part: fh.write(">%s\n%s\n" % (nm, sq))
    print("[rescore] chunk %d: %d designs" % (i//chunk + 1, len(part)), flush=True)
    subprocess.run(["colabfold_batch", "--msa-mode", "mmseqs2_uniref_env",
                    "--model-type", "alphafold2_multimer_v3",
                    "--num-models", "3", "--num-recycle", "8",
                    p, os.path.join(out, "folds")], check=False)
    n = len(glob.glob(os.path.join(out, "folds", "*_scores_rank_001_*.json")))
    print("[rescore]   total folded: %d" % n, flush=True)
if os.path.exists(os.path.join(out, ".chunk.fasta")):
    os.remove(os.path.join(out, ".chunk.fasta"))
PYEOF

say "folded: $(ls "$OUT"/folds/*_scores_rank_001_*.json 2>/dev/null | wc -l)"
/root/miniconda3/envs/colabfold/bin/python scripts/rescore_analyze.py --dir "$OUT" \
  > "$KIT/results/RESCORE_RESULT.txt" 2>&1 \
  && say "analysed -> results/RESCORE_RESULT.txt" || say "analysis FAILED"
say "=== RESCORE COMPLETE ==="
} >> "$LOG" 2>&1
