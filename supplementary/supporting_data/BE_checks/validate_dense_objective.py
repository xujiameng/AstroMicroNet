from pathlib import Path
import numpy as np,pandas as pd,time,json
from diagnostic_core import kernel_fit,evaluate,direct_score
from independent_be import categorical_fit
D=Path(__file__).resolve().parent;rows=[];labels={};start=time.time()
for ix,path in enumerate(sorted((D/'inputs').glob('*.npz'))):
 z=np.load(path);W=z['R'];truth=z['truth'].astype(bool);ii,jj=np.triu_indices(len(W),1)
 theta=float(np.quantile(W[ii,jj],.05));B=(W>=theta).astype(float);np.fill_diagonal(B,0)
 qtruth=direct_score(B,truth)
 for method in ['checked_10','checked_100','independent_64']:
  if method=='independent_64':c,q,log=categorical_fit(B,64,0)
  else:c,q=kernel_fit(B,'checked',int(method.split('_')[1]),0)
  assert abs(q-direct_score(B,c))<1e-9
  rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1],method=method,theta=theta,Q=q,Q_truth=qtruth,Q_minus_truth=q-qtruth,**evaluate(c,truth)))
  labels[path.stem+'__'+method]=c
 if (ix+1)%10==0:print('Optimization check',ix+1,'/80',round(time.time()-start,1),'s',flush=True)
df=pd.DataFrame(rows);df.to_csv(D/'results/dense_objective_validation.csv',index=False)
np.savez_compressed(D/'results/dense_objective_labels.npz',**labels)
print(df.groupby(['method','condition']).agg(mean=('f1','mean'),sd=('f1','std'),min_q_above_truth=('Q_minus_truth','min'),mean_core=('n_core','mean')).to_string(),flush=True)

