"""Regenerate all archived synthetic inputs in memory and compare numerical values."""
from pathlib import Path
import sys
import numpy as np
from scipy.signal import lfilter
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'software'))
from astromicronet.robust_solver import deterministic_seed

def main():
    support=ROOT/'supplementary/supporting_data'
    n=100;I,J=np.triu_indices(n,1);count=0
    for i in range(1,21):
        rng=np.random.default_rng(deterministic_seed(('fair_retuning_new_test_20260921_v1',i)))
        truth=np.zeros(n,bool);truth[rng.choice(n,int(rng.integers(24,29)),replace=False)]=True
        kind=truth[I].astype(int)+truth[J].astype(int)
        base=rng.uniform(np.array([.3,.6,.8])[kind],np.array([.6,.8,1.])[kind])
        noise=rng.standard_normal(len(I));background=rng.uniform(.3,1,len(I))
        bank=[base,np.clip(base+.1*noise,0,1),np.clip(base+.2*noise,0,1),.25*base+.75*background]
        for condition,v in zip(['baseline','noise_0.1','noise_0.2','overlap_0.75'],bank):
            z=np.load(support/'BE_checks/inputs'/f'V{i:02d}_{condition}.npz',allow_pickle=False)
            np.testing.assert_array_equal(z['truth'],truth)
            np.testing.assert_array_equal(z['R'][I,J],v);count+=1
    for i in range(1,501):
        rng=np.random.default_rng(deterministic_seed(('review_null_expansion_20260921_v1',i)))
        z=np.load(support/'null_networks/inputs'/f'NULL{i:04d}.npz',allow_pickle=False)
        np.testing.assert_array_equal(z['R'][I,J],rng.uniform(.3,1,len(I)))
        assert not z['truth'].any();count+=1
    for seed in range(100,120):
        rng=np.random.default_rng(92022000+seed);t=1500;k=20
        events=rng.poisson(.05/5,size=(t+100,n+1))*rng.uniform(.6,1.4,size=(t+100,n+1))
        calcium=lfilter([1],[1,-np.exp(-.2/1.2)],events,axis=0)
        common=calcium[:,0];individual=calcium[:,1:]
        truth=np.zeros(n,bool);truth[rng.permutation(n)[:k]]=True
        noise=rng.normal(0,.15,size=(t,n));delays=rng.integers(-15,16,size=n)
        for condition in ['synchronous','delayed']:
            signal=np.column_stack([common[50+(delays[j] if condition=='delayed' else 0):50+(delays[j] if condition=='delayed' else 0)+t] for j in range(n)])
            X=signal*np.where(truth,.9,.25)+individual[50:50+t]*np.where(truth,.45,.95)+noise
            z=np.load(ROOT/'validation/waveform/inputs'/f'{seed:02d}_{condition}.npz',allow_pickle=False)
            np.testing.assert_allclose(z['X'],X,rtol=0,atol=1e-14)
            np.testing.assert_array_equal(z['truth'],truth);np.testing.assert_array_equal(z['delays'],delays);count+=1
    print(f'PASS: regenerated {count} synthetic inputs (80 weighted, 500 null, 40 waveform).')
if __name__=='__main__':main()
