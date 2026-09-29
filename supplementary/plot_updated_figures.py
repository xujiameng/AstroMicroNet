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
    if name not in ['Figure_S2','Figure_S3','Figure_S4','Figure_S5','Figure_S6','Figure_S8']:
        plt.close(fig);return
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
    a=df.groupby(x).jaccard.agg(['mean','sem']).reindex(values)
    ax.errorbar(values,a['mean'],yerr=a['sem'],color=COL[0],marker='o',lw=1.5,capsize=3,ms=4)
    ax.set(xlabel=xlabel,ylabel='Core-set Jaccard index',xticks=values,ylim=(0,1.05))
sim=read('method_comparison.csv');conds=['baseline','noise_0.1','noise_0.2','overlap_0.75']
fig,ax=plt.subplots(figsize=(6.0,3.8),layout='constrained');v=[sim[sim.condition.eq(c)&sim.method.eq('AstroMicroNet')].sort_values('network').f1.values for c in conds]
box(ax,v,['Baseline','Noise σ = 0.1','Noise σ = 0.2','Overlap'],ylim=(.68,1.04));save(fig,'Figure_S1')
r=read('rho_membership.csv');r=r.drop_duplicates(['cell','rho']);r=r[r.rho.between(.18,.22)]
fig,ax=plt.subplots(figsize=(5.8,3.7),layout='constrained');sensitivity(ax,r,'rho',[.18,.19,.20,.21,.22],'Initial core fraction ρ');save(fig,'Figure_S2')
lag=read('lag_membership.csv');lag=lag[lag.stage.str.endswith('after_gate')];order=['lag0','lag1','lag2','lag5','lag10','lag30','lag60','lag120','lag300','full']
fig,axes=plt.subplots(1,2,figsize=(7.5,3.7),layout='constrained',gridspec_kw={'width_ratios':[1.3,1]})
for i,c in enumerate(order):
    v=lag[lag.condition.eq(c)].jaccard.values;q=[v.mean(),v.std(ddof=1)/np.sqrt(len(v))];axes[0].scatter(i+np.linspace(-.12,.12,len(v)),v,s=14,c='#91bbc3',alpha=.8);axes[0].errorbar(i,q[0],yerr=q[1],fmt='o',color=COL[0],capsize=3,ms=4)
axes[0].set(xticks=range(10),xticklabels=['0','±1','±2','±5','±10','±30','±60','±120','±300','Full'],xlabel='Permitted lag range (s)',ylabel='Core-set Jaccard index',ylim=(-.03,1.06));axes[0].tick_params(axis='x',labelrotation=45)
wav=read('waveform_results.csv');wav=wav[wav.comparison.eq('definition')&wav.condition.eq('delayed')];defs=['xcorr_full','pearson','spearman'];names=['Normalized\ncross-correlation','Pearson','Spearman'];pairs=read('waveform_pairwise.csv');rows=pairs[pairs.condition.eq('delayed')].replace(dict(zip(defs,names))).to_dict('records')
box(axes[1],[wav[wav.method.eq(m)].sort_values('network').f1.values for m in defs],names,rows);letter(axes[0],'A');letter(axes[1],'B');save(fig,'Figure_S3')
h=read('heldout_values.csv');h=h[h.segments.eq(3)];p=read('heldout_pairwise.csv');groups=['CC','CP','PP']
fig,axes=plt.subplots(1,2,figsize=(7.2,3.8),layout='constrained')
for ax,m,l,yl in zip(axes,['fc_retained','density'],'AB',['Mean retained-edge cross-correlation','Connection density']):
    rows=p[p.segments.eq(3)&p.metric.eq(m)].to_dict('records');box(ax,[h[g+'_'+m].values for g in groups],groups,rows,ylabel=yl);letter(ax,l)
# Preserve all ticks and their positions; extend B's spine to A's physical height.
a_bounds=axes[0].spines['left'].get_bounds();a_lim=axes[0].get_ylim();b_lim=axes[1].get_ylim()
a_frac=[(x-a_lim[0])/(a_lim[1]-a_lim[0]) for x in a_bounds]
axes[1].spines['left'].set_bounds(*[b_lim[0]+f*(b_lim[1]-b_lim[0]) for f in a_frac])
save(fig,'Figure_S4')
a=np.load(D/'average40_network.npz')['R'];b=np.load(D/'generated5_network.npz')['R'];ii,jj=np.triu_indices(len(a),1)
fig,ax=plt.subplots(figsize=(5.0,4.6),layout='constrained');ax.scatter(a[ii,jj],b[ii,jj],s=2,alpha=.18,color=COL[0],rasterized=False);ax.plot([0,1],[0,1],'--',color='#777777',lw=.8);ax.set(xlabel='Cross-correlation at 40 Hz',ylabel='Cross-correlation at 5 Hz',xlim=(0,1),ylim=(0,1));ax.set_aspect('equal');save(fig,'Figure_S5')
t=read('threshold_membership.csv');t=t[t.rho.eq(.2)&t.reference.eq('rho020_theta0')];fig,ax=plt.subplots(figsize=(5.8,3.7),layout='constrained');sensitivity(ax,t,'theta_offset',[-.05,-.025,0,.025,.05],'Offset from the adaptive connection threshold');save(fig,'Figure_S6')
c=read('cell_metrics.csv');fig,axes=plt.subplots(1,2,figsize=(7.2,3.6),layout='constrained');axes[0].scatter(c.N,c.C,s=30,c=COL[0]);axes[1].scatter(c.N,100*c.core_fraction,s=30,c=COL[0]);axes[0].set(xlabel='Total microdomain number N',ylabel='Core microdomain number C');axes[1].set(xlabel='Total microdomain number N',ylabel='Core fraction (%)',ylim=(7,22));axes[1].axhline(20,c='#888888',ls='--',lw=.8);letter(axes[0],'A');letter(axes[1],'B');save(fig,'Figure_S7')
p=read('method_pairwise.csv');methods=['AstroMicroNet','BE','Rombach'];fig,axes=plt.subplots(2,2,figsize=(7.2,6.7),layout='constrained')
for ax,cond,l in zip(axes.flat,conds,'ABCD'):
    d=sim[sim.condition.eq(cond)];box(ax,[d[d.method.eq(m)].sort_values('network').f1.values for m in methods],methods,p[p.condition.eq(cond)].to_dict('records'));letter(ax,l);ax.set_xticklabels(['AstroMicroNet','Borgatti Everett\nmodel','Rombach\nmethod']);ax.tick_params(axis='x',labelsize=8)
save(fig,'Figure_S8')

for p in [F/f'Figure_S{i}.svg' for i in (2,3,4,5,6,8)]:
    root=etree.parse(str(p));assert not root.findall(".//{http://www.w3.org/2000/svg}image")
    assert root.findall(".//{http://www.w3.org/2000/svg}text")
print("Six figures exported as editable SVG, vector PDF and PNG")
