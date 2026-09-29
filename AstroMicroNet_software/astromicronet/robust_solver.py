"""Same fixed-k objective with matrix-native diverse initialization.
No historical labels/coreness, target overlap or biological outcome enters search.
"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import numpy as np

def deterministic_seed(key):
    return int.from_bytes(hashlib.sha256(repr(key).encode()).digest()[:4],'little')

def objective(A,c):
    ix=np.flatnonzero(c)
    return float(A[np.ix_(ix,ix)].sum()/len(ix))

def polish(A,initial,max_swaps=2000,tolerance=1e-11):
    c=np.asarray(initial,bool).copy();score0=objective(A,c);history=[score0]
    for step in range(max_swaps):
        C=np.flatnonzero(c);P=np.flatnonzero(~c);strength=A[:,C].sum(1)
        gain=2*(strength[P,None]-strength[C][None,:]-A[np.ix_(P,C)])
        row,col=np.unravel_index(np.argmax(gain),gain.shape);remaining=float(gain[row,col]/len(C))
        if gain[row,col]<=tolerance:break
        c[C[col]]=False;c[P[row]]=True;history.append(objective(A,c))
    assert history[-1]>=score0-1e-10
    return c,dict(fitness=history[-1],start_fitness=score0,swaps=len(history)-1,
                 remaining_swap_gain=remaining,at_budget=bool(step==max_swaps-1 and remaining>tolerance),history=history)

def matrix_native_ranks(A,key,num_random=128,num_row_seeds=16):
    N=len(A);strength=A.sum(1);ranks=[strength];names=['node_strength']
    # Farthest-first row profiles cover different connectivity patterns.
    norm=np.linalg.norm(A,axis=1);V=A/np.maximum(norm[:,None],1e-15)
    chosen=[];distance=np.full(N,np.inf);index=int(np.argmax(strength))
    for _ in range(min(num_row_seeds,N)):
        chosen.append(index)
        score=A[index].copy();score[index]=float(A.max()+1)
        ranks.append(score);names.append('row_profile_'+str(index))
        distance=np.minimum(distance,np.sum((V-V[index])**2,axis=1));distance[chosen]=-np.inf
        index=int(np.argmax(distance))
    rng=np.random.default_rng(deterministic_seed(('rho_lag_solver_v1',key)))
    ranks.extend(rng.random((num_random,N)));names.extend('random_'+str(i) for i in range(num_random))
    return np.asarray(ranks),names

def solve(A,rho,key,num_random=128,num_row_seeds=16,ranks=None,names=None):
    A=np.asarray(A,float);N=len(A);k=int(np.floor(rho*N))
    assert 0<k<N and np.allclose(A,A.T,atol=1e-12) and np.all(np.isfinite(A))
    assert np.allclose(np.diag(A),0)
    if ranks is None:ranks,names=matrix_native_ranks(A,key,num_random,num_row_seeds)
    outputs=[];logs=[]
    for name,rank in zip(names,ranks):
        c=np.zeros(N,bool);c[np.argsort(-rank,kind='stable')[:k]]=True
        out,log=polish(A,c);outputs.append(out);logs.append(dict(start=name,**log))
    cores=np.asarray(outputs);fitness=np.asarray([x['fitness'] for x in logs])
    near=np.flatnonzero(fitness.max()-fitness<=max(1e-10,abs(fitness.max())*1e-12))
    best=min(near,key=lambda i:tuple(np.flatnonzero(cores[i])))
    return dict(candidate=cores[best],best_fitness=logs[best]['fitness'],best_start=best,
                all_candidates=cores,all_fitness=np.asarray([x['fitness'] for x in logs]),logs=logs)

def resize_candidate(A,c,k):
    c=c.copy()
    while c.sum()>k:
        C=np.flatnonzero(c);s=A[np.ix_(C,C)].sum(1);c[C[np.argmin(s)]]=False
    while c.sum()<k:
        P=np.flatnonzero(~c);s=A[P][:,c].sum(1);c[P[np.argmax(s)]]=True
    return c

def threshold2(A):
    v=np.sort(A[np.triu_indices(len(A),1)]);s=v.cumsum();s2=(v*v).cumsum();k=np.arange(1,len(v))
    e=s2[:-1]-s[:-1]**2/k+(s2[-1]-s2[:-1])-(s[-1]-s[:-1])**2/(len(v)-k)
    e[v[:-1]==v[1:]]=np.inf
    if not np.isfinite(e).any():return float(v[0])
    q=e.argmin()+1
    return float((s[q-1]/q+(s[-1]-s[q-1])/(len(v)-q))/2)

def refine(A,c,theta,all_node_voters=False,tie_tolerance=1e-12,tie_policy='first'):
    C=np.flatnonzero(c);P=np.arange(len(A)) if all_node_voters else np.flatnonzero(~c)
    x=A[np.ix_(P,C)].copy()
    if all_node_voters:x[P[:,None]==C[None,:]]=-np.inf
    maxima=x.max(axis=1)
    # C is sorted by ROI index; near-equal maxima use the smallest index.
    ix=np.argmax(x>=maxima[:,None]-tie_tolerance,axis=1)
    w=x[np.arange(len(P)),ix];valid=(maxima>=theta)&(maxima>0)
    if tie_policy=='fractional':
        ties=x>=maxima[:,None]-tie_tolerance
        support=ties/np.maximum(ties.sum(1,keepdims=True),1)
        support[~valid]=0
        votes=np.zeros(len(A),float);votes[C]=support.sum(0)
    else:
        assert tie_policy=='first'
        votes=np.zeros(len(A),int);np.add.at(votes,C[ix[valid]],1)
    return c&(votes>0),votes

def near_optimal_consensus(A,result,rho,relative_tolerance=.01):
    """Fitness-weighted ensemble within fixed relative objective tolerance.
    Top-k preserves candidate-size constraint. If fused set loses more than
    tolerance, choose the most consensus-supported qualifying original set.
    The selector never sees default labels or cross-rho overlap.
    """
    C=np.asarray(result['all_candidates'],bool);F=np.asarray(result['all_fitness'],float)
    best=float(F.max());allowed=(best-F<=max(1e-10,abs(best)*relative_tolerance))
    P=C[allowed];f=F[allowed];k=int(np.floor(rho*len(A)))
    if best<=0:
        candidate=result['candidate'].copy();votes=candidate.astype(float);kind='best_nonpositive'
    else:
        weights=np.maximum(f,0);votes=weights@P/weights.sum()
        candidate=np.zeros(len(A),bool);candidate[np.argsort(-votes,kind='stable')[:k]]=True
        value=objective(A,candidate)
        if best-value>max(1e-10,abs(best)*relative_tolerance):
            agreements=P@votes
            ix=np.flatnonzero(agreements.max()-agreements<=1e-12)
            q=max(f[ix]);ix=ix[q-f[ix]<=1e-10]
            index=min(ix,key=lambda i:tuple(np.flatnonzero(P[i])))
            candidate=P[index].copy();kind='supported_near_optimal_set'
        else:kind='fitness_weighted_top_k'
    value=objective(A,candidate)
    return candidate,votes,dict(relative_tolerance=relative_tolerance,qualifying_starts=int(allowed.sum()),
        aggregation=kind,best_fitness=best,selected_fitness=value,
        relative_objective_loss=(best-value)/abs(best) if best else 0.)

def overlap(a,b):
    i=int((a&b).sum());u=int((a|b).sum());m=min(a.sum(),b.sum())
    return dict(jaccard=i/u if u else 1.,default_retention=i/a.sum() if a.any() else np.nan,
                small_retention=i/m if m else np.nan,n_default=int(a.sum()),n_variant=int(b.sum()))
