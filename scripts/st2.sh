#!/bin/bash
K=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$K" || exit 1
echo "bench folds : $(ls results/queue/reward_bench_m1r3/*_scores_rank_001_*.json 2>/dev/null | wc -l) / 165"
echo "md systems  : $(ls results/md_gpu/*_md.json 2>/dev/null | wc -l) / 37"
echo "gpu  $(nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total --format=csv,noheader)"
echo "ram  $(free -g | sed -n 2p)"
echo "-- alive --"
ps -eo etime,args --width 400 | grep -E 'run_tonight|colabfold_batch|md_stability' | grep -v grep | cut -c1-92
echo "-- run_tonight.log --"
tail -6 results/run_tonight.log 2>/dev/null
echo "-- md_gpu3.log --"
tail -3 results/md_gpu3.log 2>/dev/null
