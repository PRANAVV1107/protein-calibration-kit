KIT=/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit
cd "$KIT" || exit 1
n=$(ps -eo args --width 300 | grep -cE '[m]cts_design\.py|[b]ash scripts/chain_remaining\.sh')
if [ "$n" -gt 0 ]; then echo "REFUSING: $n relevant process already running"; exit 1; fi
nohup setsid bash scripts/chain_remaining.sh >/dev/null 2>&1 </dev/null &
sleep 6
tail -4 results/chain_remaining.log
echo "--- running ---"
ps -eo etime,args --width 300 | grep -E 'chain_remaining|resume_bench|expert_iter|colabfold_batch' | grep -v grep | cut -c1-80
