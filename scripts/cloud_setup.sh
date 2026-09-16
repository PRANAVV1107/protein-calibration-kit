#!/bin/bash
# Run this ON the rented box, once. Installs ColabFold and unpacks the work.
#
# localcolabfold is used rather than a conda build because it ships a working
# JAX/CUDA pairing; mismatched CUDA is the single most common way this install
# fails. It downloads ~4 GB of AlphaFold weights, so expect 10-20 minutes.
#
# MSAs are NOT searched remotely -- every a3m in the bundle already carries its
# alignment, computed locally. That matters: the MMseqs2 search was the dominant
# cost in earlier runs, and paying it again on rented time would waste most of
# the money.
set -eu
WORK=${WORK:-/workspace}
mkdir -p "$WORK"
cd "$WORK"

echo "=== GPU ==="
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

if [ ! -x "$WORK/localcolabfold/colabfold-conda/bin/colabfold_batch" ]; then
  echo "=== installing ColabFold (10-20 min, ~4GB weights) ==="
  apt-get update -qq && apt-get install -y -qq wget curl git tar >/dev/null
  wget -q https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
  bash install_colabbatch_linux.sh
else
  echo "ColabFold already installed"
fi
export PATH="$WORK/localcolabfold/colabfold-conda/bin:$PATH"
colabfold_batch --help >/dev/null 2>&1 && echo "colabfold_batch OK" || { echo "INSTALL FAILED"; exit 1; }

if [ -f "$WORK/pending.tar.gz" ] && [ ! -d "$WORK/pending" ]; then
  echo "=== unpacking work ==="
  tar -xzf "$WORK/pending.tar.gz" -C "$WORK"
fi
echo "jobs:"; for d in "$WORK"/pending/*/; do
  echo "  $(basename "$d"): $(ls "$d"/*.a3m 2>/dev/null | wc -l)"
done
echo
echo "Now run:  bash cloud_run.sh"
