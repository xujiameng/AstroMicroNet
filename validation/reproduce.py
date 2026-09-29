"""Replay archived synthetic inputs with their original settings and seed keys."""
from pathlib import Path
import argparse, sys, json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'software'))
from astromicronet.pipeline import Config,fit

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bank',choices=['method','null','waveform'],default='method')
    ap.add_argument('--all',action='store_true',help='Check every network; may take substantial time.')
    args=ap.parse_args()
    support=ROOT/'supplementary/supporting_data'
    source={'method':support/'BE_checks/inputs','null':support/'null_networks/inputs','waveform':ROOT/'validation/waveform/inputs'}[args.bank]
    files=sorted(source.glob('*.npz'))
    if not args.all:files=files[:1]
    count=0
    for p in files:
        z=np.load(p,allow_pickle=False);name=p.stem
        if args.bank=='method':
            expected=np.load(ROOT/'validation/method_expected'/p.name,allow_pickle=False)
            r=fit(z['R'],Config(rho=.28),seed_key=('fair_retuning_test_solver_20260921',name),gate_key=('fair_retuning_test_gate_20260921',name))
            np.testing.assert_array_equal(r['final'],expected['AstroMicroNet_frozen_labels'])
            np.testing.assert_allclose(r['support'],expected['AstroMicroNet_frozen_scores'],rtol=0,atol=1e-12)
        elif args.bank=='null':
            i=int(name[4:]);expected=np.load(support/'null_networks/outputs'/p.name,allow_pickle=False)
            r=fit(z['R'],Config(rho=.28),seed_key=('review_null_solver_20260921',i),gate_key=('review_null_gate_20260921',i))
            np.testing.assert_array_equal(r['final'],expected['final'])
            np.testing.assert_allclose(r['support'],expected['support'],rtol=0,atol=1e-12)
        else:
            from astromicronet.pipeline import correlation_matrix
            expected=np.load(ROOT/'validation/waveform/outputs'/p.name,allow_pickle=False)
            configs=json.loads((support/'waveform/selected_configs.json').read_text())
            i=int(name.split('_')[0])
            for definition in ['xcorr_full','pearson','spearman']:
                A=correlation_matrix(z['X'],'xcorr' if definition=='xcorr_full' else definition)
                np.testing.assert_allclose(A,expected[definition+'_R'],rtol=0,atol=1e-12)
                r=fit(A,Config(**configs[definition]),np.arange(len(A)),('waveform_20260922',i,definition),('waveform_screen_20260922',i,definition))
                np.testing.assert_array_equal(r['final'],expected[definition+'_labels'])
                np.testing.assert_allclose(r['support'],expected[definition+'_scores'],rtol=0,atol=1e-12)
        print('PASS',args.bank,name,flush=True);count+=1
    assert count>0
    print(f'Checked {count} archived synthetic inputs. No results overwritten.')
if __name__=='__main__':main()
