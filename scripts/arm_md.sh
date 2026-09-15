#!/bin/bash
# Arm the MD scheduler behind the fold queue.
#
# NS=5 rather than the script's 50 default. Measured on CUDA: 0.02 ns took
# 0.8 min wall including setup, so production runs near 3 ns/hour. At 50 ns
# that is ~17 h per complex and ~250 h for the 15, which lands well past the
# deadline and only STARTS after ~40 h of folding. 5 ns x 15 is ~25 h, keeps
# all fifteen complexes rather than trading breadth for trajectory length,
# and is still long enough to see an interface decay rather than merely
# catching a pose that falls apart on contact.
KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
if ps -eo args --width 300 | grep -q '[b]ash scripts/queue_md_gpu2\.sh'; then
  echo "MD scheduler already armed"; exit 0
fi
bash -n scripts/queue_md_gpu2.sh && echo "[OK] syntax"
NS=5 EQUIL=200 nohup setsid bash scripts/queue_md_gpu2.sh >/dev/null 2>&1 </dev/null &
sleep 6
echo "--- md_gpu2.log ---"
tail -4 results/md_gpu2.log
