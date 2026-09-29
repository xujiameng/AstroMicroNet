from pathlib import Path
import numpy as np,pandas as pd,json,time,hashlib
from scipy.optimize import check_grad
from independent_be import categorical_fit,categorical_score,continuous_fit,continuous_loss,hill
from diagnostic_core import threshold2,evaluate,kernel_fit
D=Path(__file__).resolve().parent
rng=np.random.default_rng(91027);audit=[]
# Exact direct Pearson and finite-difference validation on small networks.
for example in range(12):
 n=9;raw=rng.random((n,n));W=(raw+raw.T)/2;np.fill_diagonal(W,0)
 if example%2==0:W=(W>=.5).astype(float);np.fill_diagonal(W,0)
 values=[]
 for bits in range(1,2**n-1):
  c=np.array([(bits>>i)&1 for i in range(n)],bool)
  if c.sum()>n-2:continue
  values.append(categorical_score(W,c))
 c,q,logs=categorical_fit(W,64,example)
 globalq=max(values);audit.append(dict(example=example,global_q=globalq,found_q=q,gap=globalq-q))
 assert abs(q-categorical_score(W,c))<1e-10
 assert q<=globalq+1e-9
 if example>=6:assert abs(globalq-q)<1e-9
 gradient=check_grad(lambda x:continuous_loss(x,W)[0],lambda x:continuous_loss(x,W)[1],rng.uniform(.2,.8,n))
 assert gradient<1e-5,gradient
pd.DataFrame(audit).to_csv(D/'results/independent_small_graph_audit.csv',index=False)
rows=[];labs={};fitlogs={};start=time.time()
for ix,path in enumerate(sorted((D/'inputs').glob('*.npz'))):
 z=np.load(path);W=z['R'];truth=z['truth'].astype(bool)
 for kind,A in [('weighted',W),('binary_adaptive',(W>=threshold2(W)).astype(float))]:
  A=A.copy();np.fill_diagonal(A,0);c,q,log=categorical_fit(A,32,0)
  bestflip=max(categorical_score(A,c^(np.arange(len(c))==j)) for j in range(len(c)) if 1<=(c^(np.arange(len(c))==j)).sum()<=len(c)-2)
  assert bestflip<=q+1e-9
  labs[path.stem+'__categorical_'+kind]=c
  rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1],variant='independent_categorical',input_kind=kind,rule='direct_optimized_labels',parameter='',Q=q,best_one_flip=bestflip,**evaluate(c,truth)))
 c,q,logs=continuous_fit(W,10,0);labs[path.stem+'__continuous_scores']=c;fitlogs[path.stem]=logs
 for cutoff in [.3,.4,.5,.6,.7,.8,.9]:
  pred=c>=cutoff
  rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1],variant='independent_continuous',input_kind='weighted',rule='max_normalized_score_cutoff',parameter=cutoff,Q=q,best_one_flip=np.nan,**evaluate(pred,truth)))
 for fraction in [.1,.2,.28,.4]:
  pred=np.zeros(len(c),bool);pred[np.argsort(-c,kind='stable')[:int(np.floor(len(c)*fraction))]]=True
  rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1],variant='independent_continuous',input_kind='weighted',rule='top_fraction',parameter=fraction,Q=q,best_one_flip=np.nan,**evaluate(pred,truth)))
 if (ix+1)%10==0:
  pd.DataFrame(rows).to_csv(D/'results/independent_variants_all.csv',index=False)
  print('Independent',ix+1,'/80',round(time.time()-start,1),'s',flush=True)
df=pd.DataFrame(rows);df.to_csv(D/'results/independent_variants_all.csv',index=False)
df.groupby(['variant','input_kind','rule','parameter','condition'],dropna=False).f1.agg(['count','mean','std','min','max']).to_csv(D/'results/independent_variants_summary.csv')
np.savez_compressed(D/'results/independent_labels_and_scores.npz',**labs)
(D/'results/continuous_optimizer_log.json').write_text(json.dumps(fitlogs,indent=2),encoding='utf-8')
print(df.groupby(['variant','input_kind','rule','parameter','condition'],dropna=False).f1.agg(['mean','std']).to_string())

