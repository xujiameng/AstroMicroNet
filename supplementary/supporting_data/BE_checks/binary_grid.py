from pathlib import Path
import time,json,hashlib,csv
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from diagnostic_core import kernel_fit,threshold2,evaluate,direct_score
D=Path(__file__).resolve().parent
# Diagnostic range fixed before this run; every combination is retained, no selection by desired F1.
configs=[('adaptive',0.)]+[('offset',x) for x in [-.3,-.2,-.15,-.1,-.05,.05,.1,.15,.2,.3]]+[('fixed',x) for x in [.3,.4,.5,.6,.65,.7,.75,.8,.85,.9,.95]]+[('quantile',x) for x in [.1,.25,.5,.75,.9,.95]]+[('raw_weighted',0.),('stored_zeros',0.)]
(D/'source/binary_grid_protocol.json').write_text(json.dumps({'configs':configs,'versions':['original','checked'],'num_runs':[1,10],'seeds':[0],'scope':'exploratory diagnosis; compare all cells/conditions without replacing S09'},indent=2),encoding='utf-8')
rows=[];labels={};start=time.time()
paths=sorted((D/'inputs').glob('*.npz'))
for zidx,path in enumerate(paths):
 z=np.load(path);A=z['R'];truth=z['truth'].astype(bool);n=len(A);ii,jj=np.triu_indices(n,1);auto=threshold2(A)
 assert np.allclose(A,A.T) and np.all(np.diag(A)==0)
 for kind,value in configs:
  theta=auto if kind in ['adaptive','stored_zeros'] else (auto+value if kind=='offset' else value if kind=='fixed' else np.quantile(A[ii,jj],value) if kind=='quantile' else float('nan'))
  B=A.copy() if kind=='raw_weighted' else (A>=theta).astype(float)
  np.fill_diagonal(B,0)
  if kind=='stored_zeros':
   # Deliberately retain zeros in CSR data to test a concrete sparse-input trap.
   B=csr_matrix((B.ravel(),np.tile(np.arange(n),n),np.arange(0,n*n+1,n)),shape=(n,n))
  actual=B.toarray() if kind=='stored_zeros' else B
  stored=csr_matrix(B).nnz
  valid=kind not in ['raw_weighted','stored_zeros'] and actual[ii,jj].std()>1e-14
  for version in ['original','checked']:
   for runs in [1,10]:
    c,q=kernel_fit(B,version,runs,seed=0)
    independent=direct_score(actual,c)
    if valid:assert abs(q-independent)<1e-9,(path,kind,value,version,q,independent)
    key=f'{path.stem}__{kind}_{value:g}__{version}_runs{runs}'
    labels[key]=c
    rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1],kind=kind,value=value,theta=theta,version=version,num_runs=runs,seed=0,valid_binary_use=valid,stored_entries=stored,actual_nonzero=int(np.count_nonzero(actual)),edge_density=float(np.count_nonzero(actual[ii,jj])/len(ii)),Q=q,Q_direct=independent,**evaluate(c,truth)))
 if (zidx+1)%10==0:
  pd.DataFrame(rows).to_csv(D/'results/binary_grid_all.csv',index=False)
  print(zidx+1,'/80 networks',len(rows),'fits',round(time.time()-start,1),'s',flush=True)
df=pd.DataFrame(rows)
df.to_csv(D/'results/binary_grid_all.csv',index=False)
df.groupby(['kind','value','version','num_runs','condition']).f1.agg(['count','mean','std','min','max']).to_csv(D/'results/binary_grid_summary.csv')
np.savez_compressed(D/'results/binary_grid_labels.npz',**labels)
print(df[(df.version=='checked')&df.num_runs.eq(10)&((df.kind=='fixed')|df.kind.isin(['adaptive','raw_weighted','stored_zeros']))].groupby(['kind','value','condition']).f1.agg(['mean','std','min','max']).to_string(),flush=True)

