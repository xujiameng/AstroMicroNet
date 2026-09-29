"""Independent reference objectives from Borgatti-Everett (2000).
No upstream optimization code, no true labels in either optimizer.
Categorical OR and continuous-product variants are deliberately separate.
"""
import numpy as np
from scipy.optimize import minimize
import numba
def categorical_score(W,c):
 i,j=np.triu_indices(len(W),1);a=W[i,j];b=(c[i]|c[j]).astype(float)
 if np.std(a)<1e-14 or np.std(b)<1e-14:return np.nan
 return float(np.corrcoef(a,b)[0,1])
@numba.njit(cache=True)
def hill(W,start,max_iter=1000):
 n=len(W);m=n*(n-1)/2.;c=start.copy();k=int(c.sum());total=0.;square=0.
 for i in range(n):
  for j in range(i+1,n):total+=W[i,j];square+=W[i,j]*W[i,j]
 mean=total/m;sxx=square-total*total/m
 d=W@(1.-c);pp=0.
 for i in range(n):pp+=(1-c[i])*d[i]/2.
 def score(kk,ppsum):
  nb=m-(n-kk)*(n-kk-1)/2.
  vv=nb-nb*nb/m
  if vv<1e-15 or sxx<1e-15:return -1.e99
  return (total-ppsum-mean*nb)/np.sqrt(sxx*vv)
 q=score(k,pp)
 for it in range(max_iter):
  bestq=q;best=-1;newk=k;newpp=pp
  for i in range(n):
   kk=k+1 if c[i]==0 else k-1
   if kk<1 or kk>n-2:continue
   pp2=pp-d[i] if c[i]==0 else pp+d[i]
   qi=score(kk,pp2)
   if qi>bestq+1e-12:bestq=qi;best=i;newk=kk;newpp=pp2
  if best<0:return c,q,it
  sign=-1. if c[best]==0 else 1.
  for j in range(n):d[j]+=sign*W[j,best]
  c[best]=1-c[best];k=newk;pp=newpp;q=bestq
 return c,q,max_iter
def categorical_fit(W,num_starts=32,seed=0):
 W=np.asarray(W,float);n=len(W)
 if not(np.allclose(W,W.T) and np.allclose(np.diag(W),0)):raise ValueError('Expected symmetric matrix with zero diagonal')
 rng=np.random.default_rng(seed);order=np.argsort(-W.sum(1),kind='stable')
 # Choose a data-only initialization by scanning the degree/strength order.
 trials=[]
 for k in range(1,n-1):
  c=np.zeros(n,bool);c[order[:k]]=True;trials.append((categorical_score(W,c),c))
 best_start=max(trials,key=lambda z:z[0])[1]
 starts=[best_start]
 for ix in range(num_starts-1):
  k=int(rng.integers(1,n-1));c=np.zeros(n,bool);c[rng.permutation(n)[:k]]=True;starts.append(c)
 bestq=-np.inf;best=None;logs=[]
 for c in starts:
  out,q,steps=hill(W,c.astype(np.float64));out=out.astype(bool)
  independent=categorical_score(W,out);assert abs(q-independent)<1e-9,(q,independent)
  logs.append((q,steps))
  if q>bestq:bestq=q;best=out
 return best,bestq,logs
def continuous_loss(c,W):
 n=len(c);i,j=np.triu_indices(n,1)
 x=W[i,j];x=x-x.mean()
 y=c[i]*c[j];y=y-y.mean()
 sxx=np.dot(x,x);syy=np.dot(y,y)
 if syy<1e-20 or sxx<1e-20:return 1.,np.zeros_like(c)
 den=np.sqrt(sxx*syy);q=np.dot(x,y)/den
 gy=x/den-q*y/syy
 g=np.bincount(i,weights=gy*c[j],minlength=n)+np.bincount(j,weights=gy*c[i],minlength=n)
 return -q,-g
def continuous_fit(W,num_starts=10,seed=0):
 n=len(W);rng=np.random.default_rng(seed)
 strength=W.sum(1);first=.2+.8*(strength-strength.min())/max(np.ptp(strength),1e-12)
 starts=[first]+[rng.uniform(.1,1.,n) for _ in range(num_starts-1)]
 results=[]
 for start in starts:
  r=minimize(continuous_loss,start,args=(W,),jac=True,method='L-BFGS-B',bounds=[(0.,1.)]*n,options={'maxiter':1000,'ftol':1e-12,'gtol':1e-7,'maxls':50})
  results.append(r)
 best=min(results,key=lambda r:r.fun);c=best.x/max(best.x.max(),1e-15)
 i,j=np.triu_indices(n,1);q=float(np.corrcoef(W[i,j],c[i]*c[j])[0,1])
 assert abs(q+best.fun)<1e-8,(q,float(best.fun),float(np.ptp(c)),str(best.message))
 return c,q,[dict(success=bool(r.success),nit=int(r.nit),q=float(-r.fun),message=str(r.message)) for r in results]
