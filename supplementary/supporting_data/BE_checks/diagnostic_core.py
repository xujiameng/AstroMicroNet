"""Portable diagnostic helpers. Truth labels are used only by evaluate()."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.sparse import csr_matrix
import numba
from vendor import be_original_kernels,be_checked_kernels
@numba.njit(cache=True)
def seed_numba(seed):
 np.random.seed(seed)
def threshold2(A):
 v=np.sort(A[np.triu_indices(len(A),1)]);s=v.cumsum();s2=(v*v).cumsum();k=np.arange(1,len(v))
 e=s2[:-1]-s[:-1]**2/k+(s2[-1]-s2[:-1])-(s[-1]-s[:-1])**2/(len(v)-k)
 e[v[:-1]==v[1:]]=np.inf
 if not np.isfinite(e).any():return float(v[0])
 q=e.argmin()+1
 return float((s[q-1]/q+(s[-1]-s[q-1])/(len(v)-q))/2)
def evaluate(c,t):
 c=np.asarray(c,bool);t=np.asarray(t,bool)
 tp=int(np.sum(c&t));fp=int(np.sum(c&~t));fn=int(np.sum(~c&t));tn=int(np.sum(~c&~t))
 return dict(tp=tp,fp=fp,fn=fn,tn=tn,n_core=int(c.sum()),f1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else float('nan'))
def direct_score(A,c):
 i,j=np.triu_indices(len(A),1);x=A[i,j];y=(np.asarray(c,bool)[i]|np.asarray(c,bool)[j]).astype(float)
 if x.std()<1e-14 or y.std()<1e-14:return float('nan')
 return float(np.corrcoef(x,y)[0,1])
def kernel_fit(A,version='checked',num_runs=10,seed=0):
 B=csr_matrix(A,dtype=np.float64)
 mod=be_checked_kernels if version=='checked' else be_original_kernels
 seed_numba(seed);bestq=-np.inf;best=None;n=len(B.indptr)-1
 for _ in range(num_runs):
  c=mod._kernighan_lin_(B.indptr,B.indices,B.data,n).astype(np.int64)
  q,_=mod._score_(B.indptr,B.indices,B.data,np.zeros(n,dtype=np.int64),c,n)
  if best is None or q>bestq:bestq=float(q);best=c.copy()
 return best.astype(bool),bestq
