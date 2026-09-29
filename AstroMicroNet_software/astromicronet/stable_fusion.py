"""Numerically consistent original fitness fusion; no new network objective.
Near-equal fused scores are ordered by immutable original ROI record IDs.
Groups are anchored at their largest score, not decimal rounding buckets.
"""
import numpy as np
FUSION_TIE_TOLERANCE=1e-12
def stable_score_order(scores,original_roi_ids,tolerance=FUSION_TIE_TOLERANCE):
    v=np.asarray(scores,float);ids=np.asarray(original_roi_ids)
    assert v.ndim==ids.ndim==1 and len(v)==len(ids) and len(np.unique(ids))==len(ids)
    order=np.lexsort((ids,-v));result=[];left=0
    while left<len(order):
        right=left+1;anchor=v[order[left]]
        while right<len(order) and anchor-v[order[right]]<=tolerance:right+=1
        group=order[left:right];result.extend(group[np.argsort(ids[group],kind='stable')]);left=right
    return np.asarray(result,dtype=int)
def fuse_stable(A,all_candidates,all_fitness,rho,original_roi_ids=None,aggregation='dot'):
    C=np.asarray(all_candidates,bool);f=np.maximum(np.asarray(all_fitness,float),0);assert f.sum()>0
    if aggregation=='dot':v=f@C/f.sum()
    elif aggregation=='average':v=np.average(C,axis=0,weights=f)
    elif aggregation=='sequential':
        v=np.zeros(C.shape[1])
        for weight,c in zip(f,C):v+=weight*c
        v/=f.sum()
    else:raise ValueError(aggregation)
    ids=np.arange(len(A)) if original_roi_ids is None else np.asarray(original_roi_ids)
    candidate=np.zeros(len(A),bool);candidate[stable_score_order(v,ids)[:int(np.floor(rho*len(A)))]]=True
    return candidate,v
