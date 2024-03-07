# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:17:52 2023

@author: avanikoparkar
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import seaborn as sns

wavpath='D:/analysis/data_for_elife_mMAN/Figure2_figuresupplement4_sourcedata1'
#listbirds=os.listdir(wavpath)
listbirds=list(['bird'+str(a) for a in range(1,8)])
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = "Arial"
# Then, "ALWAYS use sans-serif fonts"
plt.rcParams['font.family'] = "sans-serif"


fig, axes = plt.subplots(7, 5, sharex='col',figsize=(8.3, 11.7))
#fig.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}')) # No decimal places

#fig.suptitle(f'TEST: tight_layout={tl} constrained_layout={cl}')
fig.supxlabel('Pre lesion')
fig.supylabel('Post lesion')
cols = [col for col in ['entS','meanT (ms)','F1 (kHz)','F2 (kHz)','F3 (kHz)']]
rows = ['Bird {}'.format(row) for row in range(1,8)]
markers=["o","v","^","<",">","s","p","*","h","d"]

for numbirds,birds in enumerate(listbirds):

    # first specify the path of the bird
    birdpath='D:/analysis/data_for_elife_mMAN/Figure2_figuresupplement4_sourcedata1/'+birds

    df=pd.read_csv('D:/analysis/data_for_elife_mMAN/Figure2_figuresupplement4_sourcedata1/'+birds+'_new.csv',sep='\t')
    vocSelTableGrouped = df.groupby(['calltype','preorpost'])
    vocSelTableGroupedAgg = vocSelTableGrouped.aggregate('mean', numeric_only=True).reset_index()
    
    #print(vocSelTableGroupedAgg)
    fontsize=1
    
    #now do it for all birds
    # plot multiple violin plots
    XFeatureNames = np.hstack(('entS','meanT','F1','F2','F3'))
    #set seaborn plotting aesthetics as default
    #sns.set()
    
    #define plotting region (1 row, 2 columns)
    
    #plt.subplots_adjust(top=0.936,
#bottom=0.053,
#left=0.092,
#right=0.962,
#hspace=0.295,
#wspace=0.285)
    #fig.tight_layout()
    mins=[0,0.08,750,2000,3000]
    maxs=[1,0.2,3000,5000,6000]
    units=['','ms','kHz','kHz','kHz']
    for num in range(0,5):
        df_grouped = df.groupby(['preorpost','nameofsyl'])
        df_groupedagg = df_grouped.aggregate('mean', numeric_only=True).reset_index()
        max_lim=maxs[num]
        min_lim=mins[num]
        #print(df_groupedagg)
        fontsize=1
        df_cv=df_groupedagg[["preorpost","nameofsyl",XFeatureNames[num]]]
        #make new dataframe with columns nameofsyl, cvfund_pre, cvfund_post
        df_sub_pre=df_cv[df_cv["preorpost"]=="pre"]
        df_sub_pos=df_cv[df_cv["preorpost"]=="pos"]
        xcol=df_sub_pre[XFeatureNames[num]].tolist()
        ycol=df_sub_pos[XFeatureNames[num]].tolist()
        sylnames=df_sub_pre["nameofsyl"].tolist()
        df_new=pd.DataFrame({'nameofsyls':sylnames,'PRE':xcol,'POST':ycol})
        #sns.violinplot(data=df,x="nameofsyl",y=XFeatureNames[num],hue='preorpost',ax=axes[numbirds,num],fontsize=fontsize)

        #axes[numbirds,num].set_title(XFeatureNames[num])
        axes[numbirds,num].spines[['top', 'right']].set_visible(False)
        axes[numbirds,num].axline((0, 0), slope=1,linewidth=1, color='k',linestyle='--')    
        axes[numbirds,num].yaxis.grid(False)
        axes[numbirds,num].set_xlim(min_lim, max_lim)
        axes[numbirds,num].set_ylim(min_lim, max_lim)
        #axes[numbirds,num].yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}e3'))
        ax=sns.scatterplot(data=df_new,x="PRE",y="POST",style='nameofsyls',hue='nameofsyls',markers=markers,ax=axes[numbirds,num],legend=False,s=100)
        #xlabels = ['{:,.2f}'.format(x) + '10e3' for x in g.get_xticks()/1000]
        #g.set_xticklabels(xlabels)
        if num>1:
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1000)))
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{:,.0f}'.format(x/1000)))
        axes[numbirds,num].set(xlabel=None, ylabel=None)
   
    #fig.tight_layout()
for ax, col in zip(axes[0], cols):
    ax.set_title(col)

for ax, row in zip(axes[:,0], rows):
    ax.set_ylabel(row, rotation=90, size='large')

fig.subplots_adjust(top=0.88,
bottom=0.11,
left=0.125,
right=0.915,
hspace=0.2,
wspace=0.48)
#fig.tight_layout()
plt.show()
#fig.savefig('D:/analysis/figures/forpublication_finalfigs/pitch_analysis/allbirds_scatter.svg')








