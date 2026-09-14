#!/bin/bash
K=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
c(){ ls "$1"/*_scores_rank_001_*.json 2>/dev/null | wc -l; }
echo "highacc  $(c $K/results/queue/adaptyv_highacc_folds) / 165"
echo "v2_il7ra $(c $K/results/queue/v2refold_folds_il7ra)"
echo "v2_pdl1  $(c $K/results/queue/v2refold_folds_pdl1)"
echo "v2_rbd   $(c $K/results/queue/v2refold_folds_rbd)"
echo "null     $(c $K/results/queue/null_panel_folds) / 420"
echo "-- folding processes --"
ps -eo pid,etime,args | grep 'colabfold_batch --model' | grep -v grep
echo "-- schedulers --"
ps -eo pid,etime,args | grep -E 'queue_v10\.sh|queue_md_gpu\.sh' | grep -v grep
echo "-- queue_v10.log --"
tail -8 "$K/results/queue_v10.log" 2>/dev/null || echo "(no log)"
echo "-- results dir --"
ls -la "$K/results/"ADAPTYV_*.txt 2>/dev/null
