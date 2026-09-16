#!/bin/bash
# Bundle the pending fold inputs for upload to a rented GPU box.
#
# a3m files are plain text and compress about 10:1, so this ships ~100 MB
# instead of ~900 MB. Only PENDING work is included -- anything already folded
# locally is skipped, so a resumed bundle after a preempted instance is small.
set -u
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
Q=$KIT/results/queue
OUT=$KIT/results/cloud
mkdir -p "$OUT/pending"
cd "$KIT" || exit 1

pending(){   # in_dir  out_dir  label
  local src="$1" done_dir="$2" label="$3" n=0
  mkdir -p "$OUT/pending/$label"
  for f in "$src"/*.a3m; do
    [ -e "$f" ] || continue
    local b; b=$(basename "$f" .a3m)
    if ! ls "$done_dir/${b}_scores_rank_001_"*.json >/dev/null 2>&1; then
      cp "$f" "$OUT/pending/$label/"; n=$((n+1))
    fi
  done
  echo "  $label: $n pending"
}

echo "collecting pending inputs..."
pending "$Q/reward_bench_in"          "$Q/reward_bench_m1r3"        reward_bench
pending "$Q/null_panel"               "$Q/null_panel_folds"         null_panel
pending "$KIT/results/rl/expert_r1/fold_in" "$KIT/results/rl/expert_r1/folds" expert_r1
for t in pdl1 rbd; do
  pending "$Q/v2refold/$t" "$Q/v2refold_folds_$t" "v2_$t"
done

echo "compressing..."
tar -czf "$OUT/pending.tar.gz" -C "$OUT" pending
echo
echo "bundle: $(du -h "$OUT/pending.tar.gz" | cut -f1)  ->  $OUT/pending.tar.gz"
echo "files : $(find "$OUT/pending" -name '*.a3m' | wc -l)"
