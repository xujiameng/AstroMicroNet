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
assert (BASE/'data').is_dir()
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
# Repository variant: omit microscopy panel A; preserve the archived B-F values and plotting.
stats=read('H1R_statistics.csv').set_index('measure');df=read('H1R_cell_metrics.csv');means=read('H1R_slice_means.csv');traces=read('H1R_slice_traces.csv')
fig=plt.figure(figsize=(9.1,6.6),layout='constrained');grid=fig.add_gridspec(2,3);sub=grid[0,0].subgridspec(1,2,wspace=.06);audit=[]
ax=fig.add_subplot(grid[0,0]); ax.axis('off')
ax.text(.5,.5,'Panel A microscopy images\nare not distributed in this repository.',ha='center',va='center',transform=ax.transAxes,fontsize=10)
ax=fig.add_subplot(grid[0,1]);letter(ax,'B')
for g in ['WT','cKO']:
    t=traces[traces.group.eq(g)]
    for _,d in t.groupby('slice_id'):ax.plot(d.time_relative_s,100*d.mean_local_dff,color=GC[g],alpha=.15,lw=.6)
    q=t.groupby('time_relative_s').mean_local_dff;m=q.mean();se=q.sem();ax.fill_between(m.index,100*(m-se),100*(m+se),color=GC[g],alpha=.2);ax.plot(m.index,100*m.values,c=GC[g],lw=1.8,label=g)
ax.axvline(0,c='#888888',ls='--',lw=.7);ax.axvspan(10,130,color='#aaaaaa',alpha=.15);ax.set(xlabel='Time from histamine arrival (s)',ylabel='Mean local ΔF/F (%)',xlim=(-60,180));ax.legend(frameon=False,fontsize=8);lo,hi=ax.get_ylim();y=hi+.03*(hi-lo);ax.plot([10,10,130,130],[y,y+3,y+3,y],c='#666666',lw=.7);ax.text(70,y+4,sig(stats.loc['mean_local_dff_post','q_bh']),ha='center');ax.set_ylim(lo,y+18)
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from matplotlib.cbook import boxplot_stats
box_audit=[]
point_rng=np.random.default_rng(20260928)

def box_and_points(ax,values,position,group,metric,width,jitter,size):
    values=np.asarray(values,dtype=float)
    assert len(values)==9 and np.isfinite(values).all()
    color=GC[group]
    bp=ax.boxplot([values],positions=[position],widths=width,patch_artist=True,
                  showfliers=False,whis=1.5,manage_ticks=False,
                  boxprops={'facecolor':to_rgba(color,.18),'edgecolor':color,'linewidth':1.0},
                  medianprops={'color':'#222222','linewidth':1.3},
                  whiskerprops={'color':color,'linewidth':.9},
                  capprops={'color':color,'linewidth':.9},zorder=2)
    delta=point_rng.uniform(-jitter,jitter,len(values))
    points=ax.scatter(position+delta,values,s=size,color=color,alpha=.85,
                      edgecolors='white',linewidths=.25,zorder=3)
    # Preserve median visibility when near-coincident points lie on the center line.
    bp['medians'][0].set_zorder(4)
    bp['boxes'][0].set_gid(f'box_{metric}_{group}')
    points.set_gid(f'points_{metric}_{group}')
    b=boxplot_stats(values,whis=1.5)[0]
    box_audit.append({'metric':metric,'group':group,'n':len(values),
                      'values':values.tolist(),'q1':float(b['q1']),'median':float(b['med']),
                      'q3':float(b['q3']),'whisker_low':float(b['whislo']),'whisker_high':float(b['whishi'])})

def paired(ax,m,l,yl,factor=1):
    w=means.pivot(index='slice_id',columns='group',values=m);letter(ax,l)
    for i,g in enumerate(['WT','cKO']):
        box_and_points(ax,w[g].to_numpy()*factor,i,g,m,.42,.09,22)
    vals=w.to_numpy()*factor;lo,hi=vals.min(),vals.max();s=max(hi-lo,.05);y=hi+.13*s
    ax.plot([0,0,1,1],[y,y+.035*s,y+.035*s,y],c='#666666',lw=.7)
    ax.text(.5,y+.05*s,sig(stats.loc[m,'q_bh']),ha='center')
    ax.set(xticks=[0,1],xticklabels=['WT','H1R cKO'],ylabel=yl,xlim=(-.35,1.35),ylim=(lo-.1*s,hi+.36*s))
paired(fig.add_subplot(grid[0,2]),'core_fraction','C','Core fraction (%)',100)
ax=fig.add_subplot(grid[1,0]);letter(ax,'D')
for shift,g in [(-.16,'WT'),(.16,'cKO')]:
    d=means[means.group.eq(g)]
    for j,cc in enumerate(['CC','CP','PP']):
        v=d[cc+'_density'].values
        box_and_points(ax,v,j+shift,g,cc+'_density',.26,.055,14)
for j,cc in enumerate(['CC','CP','PP']):
    y=1.03;ax.plot([j-.16,j-.16,j+.16,j+.16],[y,y+.025,y+.025,y],c='#666666',lw=.7);ax.text(j,y+.035,sig(stats.loc[cc+'_density','q_bh']),ha='center')
ax.set(xticks=[0,1,2],xticklabels=['CC','CP','PP'],ylabel='Connection density',ylim=(0,1.17))
ax.legend(handles=[Patch(facecolor=to_rgba(GC[g],.18),edgecolor=GC[g],label=g) for g in ['WT','cKO']],frameon=False,fontsize=8)
paired(fig.add_subplot(grid[1,1]),'density_CC_minus_PP','E','CC density − PP density');paired(fig.add_subplot(grid[1,2]),'mean_pairwise_FC','F','Mean pairwise FC');save(fig,'Figure_S9_quantitative_only')

(BASE/'S9_boxplot_values.json').write_text(json.dumps(box_audit,ensure_ascii=False,indent=2),encoding='utf-8')
