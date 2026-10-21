cd /mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit || exit 1
/root/miniconda3/envs/colabfold/bin/python - <<'PY'
import glob, json, os
import numpy as np, pandas as pd
K='/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit/'
def load(d, tag):
    r={}
    for f in glob.glob(K+d+'/*_scores_rank_001_*.json'):
        t=os.path.basename(f).split('_scores_rank_001_')[0]
        r[t]=float(json.load(open(f)).get('iptm') or 0)
    return pd.Series(r, name=tag)
hi = load('results/queue/adaptyv_highacc_folds','high')   # 3 models, 8 recycles
lo = load('results/queue/reward_bench_m1r3','m1r3')       # 1 model, 3 recycles
d = pd.concat([hi,lo],axis=1).dropna()
print('designs folded BOTH ways: %d'%len(d))
print()
print('              high    m1r3   ratio')
for q in (0.25,0.5,0.75,0.9,1.0):
    h=d.high.quantile(q); l=d.m1r3.quantile(q)
    print('  q%-4.2f     %.3f   %.3f   %.2f'%(q,h,l,(l/h if h else float('nan'))))
print()
print('mean high %.3f   mean m1r3 %.3f   mean ratio %.2f'%(d.high.mean(),d.m1r3.mean(),(d.m1r3/d.high.replace(0,np.nan)).mean()))
print('max  high %.3f   max  m1r3 %.3f'%(d.high.max(),d.m1r3.max()))
print()
print('spearman high vs m1r3: %.3f'%d.high.corr(d.m1r3,method='spearman'))
print()
print('How many designs reach >=0.87 under each setting?')
print('  high accuracy: %d of %d'%((d.high>=0.87).sum(),len(d)))
print('  m1r3        : %d of %d'%((d.m1r3>=0.87).sum(),len(d)))
PY
