"""Run the diagnostic on another computer; no repository-specific paths."""
from pathlib import Path
import argparse,csv,json,hashlib,platform,sys
import numpy as np,scipy,numba
from scipy.io import loadmat
from diagnostic_core import threshold2,kernel_fit,direct_score,evaluate
from independent_be import categorical_fit
def sha_array(x,dtype):return hashlib.sha256(np.ascontiguousarray(x,dtype=dtype).tobytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--input',type=Path,default=Path(__file__).resolve().parent/'inputs')
 ap.add_argument('--output',type=Path,default=Path('crosscheck_output'))
 ap.add_argument('--mode',choices=['quantile','adaptive','fixed'],default='quantile')
 ap.add_argument('--value',type=float,default=.05,help='Quantile in [0,1], or fixed threshold; ignored for adaptive')
 ap.add_argument('--offset',type=float,default=0.,help='Added only to adaptive threshold')
 ap.add_argument('--runs',type=int,default=10)
 ap.add_argument('--seed',type=int,default=0)
 ap.add_argument('--implementation',choices=['checked','original','independent-binary','independent-weighted'],default='checked')
 ap.add_argument('--matrix-key',default='R');ap.add_argument('--truth-key',default='truth')
 args=ap.parse_args()
 if args.runs<1:ap.error('--runs must be positive')
 if args.mode=='quantile' and not 0<=args.value<=1:ap.error('--value must be in [0,1] for quantiles')
 paths=sorted(args.input.glob('*.npz')) if args.input.is_dir() else [args.input]
 if not paths:raise ValueError('No NPZ files found; pass one NPZ or MAT file using --input')
 args.output.mkdir(parents=True,exist_ok=True);(args.output/'labels').mkdir(exist_ok=True)
 rows=[]
 for path in paths:
  data=loadmat(path) if path.suffix.lower()=='.mat' else np.load(path,allow_pickle=False)
  W=np.asarray(data[args.matrix_key],dtype=np.float64);raw=np.asarray(data[args.truth_key]).ravel()
  if W.ndim!=2 or W.shape[0]!=W.shape[1] or not np.all(np.isfinite(W)):raise ValueError(f'{path.name}: expected finite square matrix')
  if len(raw)!=len(W) or not np.isin(raw,[0,1]).all():raise ValueError('truth must be a 0/1 vector in the SAME node order as R, not a list of indices')
  if not np.allclose(W,W.T,atol=1e-12) or not np.allclose(np.diag(W),0):raise ValueError('Expected undirected network without self loops; do not silently change it')
  truth=raw.astype(bool);ii,jj=np.triu_indices(len(W),1)
  if args.implementation=='independent-weighted':
   theta=None;B=None
   c,q,_=categorical_fit(W,args.runs,args.seed);qcheck=direct_score(W,c)
  else:
   theta=threshold2(W)+args.offset if args.mode=='adaptive' else float(np.quantile(W[ii,jj],args.value)) if args.mode=='quantile' else args.value
   B=(W>=theta).astype(float);np.fill_diagonal(B,0)
   if args.implementation=='independent-binary':c,q,_=categorical_fit(B,args.runs,args.seed)
   else:c,q=kernel_fit(B,args.implementation,args.runs,args.seed)
   qcheck=direct_score(B,c)
  if not np.isfinite(qcheck):raise ValueError('Constant input/pattern gives undefined Pearson objective; this is not a valid BE comparison')
  if abs(q-qcheck)>1e-9:raise ValueError('Objective check failed')
  rows.append(dict(input=path.stem,condition=path.stem.split('_',1)[1] if '_' in path.stem else '',input_kind='weighted' if B is None else 'binary',mode='not_applied' if B is None else args.mode,value=None if B is None else args.value,offset=None if B is None else args.offset,theta=theta,implementation=args.implementation,num_runs=args.runs,seed=args.seed,N=len(W),matrix_sha256=sha_array(W,'<f8'),truth_sha256=sha_array(truth,'u1'),binary_sha256=None if B is None else sha_array(B,'u1'),edge_density=None if B is None else float(B[ii,jj].mean()),Q=q,Q_direct=qcheck,**evaluate(c,truth)))
  with (args.output/'labels'/f'{path.stem}.csv').open('w',encoding='utf-8-sig',newline='') as f:
   writer=csv.writer(f);writer.writerow(['node_index_0based','node_index_1based','truth_core','predicted_core'])
   writer.writerows((i,i+1,int(truth[i]),int(c[i])) for i in range(len(c)))
 with (args.output/'results.csv').open('w',encoding='utf-8-sig',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 versions={'python':platform.python_version(),'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'numba':numba.__version__,'command':sys.argv,'source_hashes':{str(p.relative_to(Path(__file__).resolve().parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).resolve().parent/'diagnostic_core.py',Path(__file__).resolve().parent/'independent_be.py',*sorted((Path(__file__).resolve().parent/'vendor').glob('*.py'))]}}
 (args.output/'environment.json').write_text(json.dumps(versions,ensure_ascii=False,indent=2),encoding='utf-8')
 for condition in sorted(set(r['condition'] for r in rows)):
  a=np.array([r['f1'] for r in rows if r['condition']==condition])
  print(f'{condition}: n={len(a)}, F1 mean={a.mean():.6f}, SD={a.std(ddof=1) if len(a)>1 else 0:.6f}')
 print('Saved',args.output)
if __name__=='__main__':main()
