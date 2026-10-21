cd /mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit || exit 1
unset PYTHONPATH
/root/miniconda3/envs/colabfold/bin/python - <<'PY'
import glob, json, os
import numpy as np
vals=[]
for d in ('results/rl/mcts_run1','results/rl/mcts_run2_mean','results/rl/noise_poc'):
    for f in glob.glob(d+'/**/*_scores_rank_001_*.json', recursive=True):
        if 'il7ra' in os.path.basename(f): continue
        try: vals.append(float(json.load(open(f)).get('iptm') or 0))
        except Exception: pass
v=np.array([x for x in vals if x>0])
print("pooled RBD folds: %d"%len(v))
print("  median %.3f   p90 %.3f   p99 %.3f   max %.3f"%(np.median(v),*np.percentile(v,[90,99]),v.max()))
print("  frac >0.80: %.2f%%   >0.85: %.2f%%   >0.90: %.2f%%"
      %(100*(v>0.80).mean(),100*(v>0.85).mean(),100*(v>0.90).mean()))
print()
print("Expected best-of-N (bootstrap, 2000 reps):")
rng=np.random.default_rng(0)
prev=None
for N in (60,120,250,500,1000,2000,4000,8000):
    if N>len(v): 
        m=np.mean([rng.choice(v,N,replace=True).max() for _ in range(2000)])
    else:
        m=np.mean([rng.choice(v,N,replace=False).max() for _ in range(2000)])
    d="" if prev is None else "  (+%.3f vs prev)"%(m-prev)
    print("  N=%-5d  expected max %.3f%s"%(N,m,d))
    prev=m
PY
