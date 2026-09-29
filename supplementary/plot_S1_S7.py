"""Replot archived results without rerunning analysis or changing parameters.
Run next to data/; outputs are genuine vector SVG/PDF (image backgrounds excepted).
"""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from lxml import etree
BASE=Path(__file__).resolve().parent
assert (BASE/'data').is_dir(), 'The data folder must be adjacent to this script.'
D=BASE/'data';F=BASE/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'white','axes.linewidth':.8})
COL=['#287E8E','#C87935','#79619F','#6199A3'];GC={'WT':'#286b98','cKO':'#ce7535'}
def read(n):return pd.read_csv(D/n)
def save(fig,name):
    for ext in ['svg','pdf','png']:fig.savefig(F/f'{name}.{ext}',dpi=250,bbox_inches='tight')
    plt.close(fig)
def letter(ax,s):ax.text(-.15,1.03,s,transform=ax.transAxes,fontweight='bold',fontsize=12,va='bottom')
def sig(p):return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'
def box(ax,vals,labels,rows=None,ylabel='F1 score',ylim=None):
    bp=ax.boxplot(vals,positions=np.arange(len(vals)),widths=.48,showfliers=False,patch_artist=True,medianprops={'color':'#333333','linewidth':1},whiskerprops={'color':'#666666'},capprops={'color':'#666666'})
    for patch,c in zip(bp['boxes'],COL):patch.set_facecolor(c);patch.set_alpha(.22);patch.set_edgecolor('#777777')
    for i,(v,c) in enumerate(zip(vals,COL)):
        # Jitter depends only on point order, never on outcome.
        rng=np.random.default_rng(20260927+i);ax.scatter(i+rng.uniform(-.15,.15,len(v)),v,s=21,c=c,alpha=.85,edgecolors='white',linewidths=.25,zorder=3)
    ax.set(xticks=np.arange(len(labels)),xticklabels=labels,ylabel=ylabel)
    lo=min(np.min(v) for v in vals);hi=max(np.max(v) for v in vals);span=max(hi-lo,.10)
    if rows is not None:
        for j,r in enumerate(rows):
            a,b=labels.index(r['a']),labels.index(r['b']);y=hi+span*(.13+j*.17);h=.035*span
            ax.plot([a,a,b,b],[y,y+h,y+h,y],lw=.7,color='#666666');ax.text((a+b)/2,y+h+.008*span,sig(r.get('adjusted_p',r.get('q_bh'))),ha='center',va='bottom',fontsize=9)
        ax.set_ylim(lo-.12*span,hi+span*(.25+len(rows)*.17))
    elif ylim:ax.set_ylim(*ylim)
    else:ax.set_ylim(lo-.12*span,hi+.16*span)
    # Reserve space for significance brackets without labeling impossible values.
    lim=ax.get_ylim();ticks=[v for v in ax.get_yticks() if lim[0]<=v<=min(1.,lim[1])+1e-10]
    ax.set_yticks(ticks);ax.set_ylim(lim);ax.spines['left'].set_bounds(lim[0],min(1.,lim[1]))
def sensitivity(ax,df,x,values,xlabel):
    for _,d in df.groupby('cell'):
        d=d.set_index(x).reindex(values);ax.plot(values,d.jaccard,c='#a9cbd3',alpha=.55,lw=.8)
    a=df.groupby(x).jaccard.quantile([.25,.5,.75]).unstack().reindex(values)
    ax.errorbar(values,a[.5],yerr=[a[.5]-a[.25],a[.75]-a[.5]],color=COL[0],marker='o',lw=1.5,capsize=3,ms=4)
    ax.set(xlabel=xlabel,ylabel='Core-set Jaccard index',xticks=values,ylim=(0,1.05))
sim=read('method_comparison.csv');conds=['baseline','noise_0.1','noise_0.2','overlap_0.75']
fig,ax=plt.subplots(figsize=(6.0,3.8),layout='constrained');v=[sim[sim.condition.eq(c)&sim.method.eq('AstroMicroNet')].sort_values('network').f1.values for c in conds]
box(ax,v,['Baseline','Noise σ = 0.1','Noise σ = 0.2','Overlap'],ylim=(.68,1.04));save(fig,'Figure_S1')
c=read('cell_metrics.csv');fig,axes=plt.subplots(1,2,figsize=(7.2,3.6),layout='constrained');axes[0].scatter(c.N,c.C,s=30,c=COL[0]);axes[1].scatter(c.N,100*c.core_fraction,s=30,c=COL[0]);axes[0].set(xlabel='Total microdomain number N',ylabel='Core microdomain number C');axes[1].set(xlabel='Total microdomain number N',ylabel='Core fraction (%)',ylim=(7,22));axes[1].axhline(20,c='#888888',ls='--',lw=.8);letter(axes[0],'A');letter(axes[1],'B');save(fig,'Figure_S7')
