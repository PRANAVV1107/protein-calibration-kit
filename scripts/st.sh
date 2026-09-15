#!/bin/bash
# Compact pipeline status. Lives in a file because passing this through
# `wsl -- bash -c "..."` collapses newlines and eats $VAR expansions, which has
# twice produced a status line reading 0 folds when 70 were on disk.
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
c(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }

echo "highacc   $(c results/queue/adaptyv_highacc_folds) / 165"
echo "v2_il7ra  $(c results/queue/v2refold_folds_il7ra) / 318"
echo "v2_pdl1   $(c results/queue/v2refold_folds_pdl1) / 285"
echo "v2_rbd    $(c results/queue/v2refold_folds_rbd) / 278"
echo "null      $(c results/queue/null_panel_folds) / 420"
echo "md        $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l) / 15"
echo
echo "gpu  $(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader)"
echo "ram  $(free -g | sed -n 2p)"
echo
echo "-- alive --"
# --width must exceed the full interpreter+script path (~100 chars) or ps
# truncates the args BEFORE the word colabfold_batch and the grep misses the
# very process being looked for. Truncate for display only, after matching.
ps -eo etime,args --width 400 2>/dev/null \
  | grep -E 'envs/colabfold/bin/colabfold_batch|md_stability\.py|queue_v11\.sh|watch_highacc\.sh|queue_md_gpu2\.sh' \
  | grep -v grep | cut -c1-96
echo
echo "-- queue_v11 --"
tail -5 results/queue_v11.log 2>/dev/null
echo "-- watcher --"
tail -3 results/watch_highacc.log 2>/dev/null
