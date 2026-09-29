from pathlib import Path
import numpy as np,pandas as pd,json,time
from scipy.sparse import csr_matrix
from diagnostic_core import kernel_fit,threshold2,evaluate
D=Path(__file__).resolve().parent
qs=[.01,.025,.05,.075,.1,.125,.15,.2]
(D/'source/dense_input_protocol.json').write_text(json.dumps({'quantiles':qs,'retained_edge_fraction':'approximately 1-quantile','num_runs':10,'seeds':[0,1,2],'versions':['checked','original'],'purpose':'Follow-up range diagnosis; all results reported, no change to S09'},indent=2),encoding='utf-8')
rows=[];predictions={};start=time.time()
for ix,path in enumerate(sorted((D/'inputs').glob('*.npz'))):
 z=np.load(path);W=z['R'];t=z['truth'].astype(bool);i,j=np.triu_indices(len(t),1)
 for q in qs:
  theta=float(np.quantile(W[i,j],q));B=(W>=theta).astype(float);np.fill_diagonal(B,0)
  for ver in ['original','checked']:
   for seed in [0,1,2]:
    c,Q=kernel_fit(B,ver,10,seed)
    rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1],quantile=q,retained_fraction=float(B[i,j].mean()),theta=theta,version=ver,seed=seed,num_runs=10,Q=Q,**evaluate(c,t)))
    predictions[f'{path.stem}__q{q}__{ver}_seed{seed}']=c
 if (ix+1)%10==0:
  pd.DataFrame(rows).to_csv(D/'results/dense_input_all.csv',index=False)
  print('Density follow-up',ix+1,'/80',round(time.time()-start,1),'s',flush=True)
df=pd.DataFrame(rows);df.to_csv(D/'results/dense_input_all.csv',index=False)
df.groupby(['quantile','version','seed','condition']).f1.agg(['count','mean','std','min','max']).to_csv(D/'results/dense_input_summary.csv')
np.savez_compressed(D/'results/dense_input_labels.npz',**predictions)
print(df[df.version.eq('checked')&df.seed.eq(0)].groupby(['quantile','condition']).f1.agg(['mean','std','min','max']).to_string(),flush=True)

