#!/bin/bash
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
pkill -f 'bash scripts/queue_md_gpu2\.sh' && echo "stopped queue_md_gpu2 (old 15-complex il7ra panel)" || echo "queue_md_gpu2 not running"
sleep 2
if ps -eo args --width 300 | grep -q '[b]ash scripts/queue_md_gpu3\.sh'; then echo "v3 already armed"; exit 0; fi
bash -n scripts/queue_md_gpu3.sh && echo "[OK] syntax"
NS=2 EQUIL=200 nohup setsid bash scripts/queue_md_gpu3.sh >/dev/null 2>&1 </dev/null &
sleep 6
tail -5 results/md_gpu3.log
