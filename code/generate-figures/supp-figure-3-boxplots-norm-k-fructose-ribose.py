#!/usr/bin/env python

# import off-the-shelf packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# import custom library
from utils import *

# READ DATA
df_summ = read_csv("../../amiga-biolog/summary/merged_summary_norm_sub.txt")

# summary metric of interest: 
# subtraction-normalized carrying capacity in the untransformed lienar scale
param = 'norm(k_lin)'

# the following are variables needed for this figure
varbs = ['Strain_ID','Ribotype','Clade','SequenceType','Substrate',param]

# get median carrying capacity for each strain on each substrate
df_long = df_summ.loc[:,varbs]
df_long = df_long.groupby(['Strain_ID','Substrate'])
df_long = df_long.median(param).reset_index()

# pivot dataframe
df_wide = df_long.pivot(index='Substrate',columns='Strain_ID',values=param).T

# get strain typing data
df_meta = df_summ.loc[:,['Strain_ID','Ribotype','Clade']].drop_duplicates()
df_meta = df_meta.replace({'RTUnk':'Unknown'}).fillna('Unknown')
df_meta = df_meta.set_index('Strain_ID')

# substrates to plot
substrates = ['D-Fructose','D-Ribose']

# define data
data = df_wide.loc[:,substrates]
data = data.unstack().to_frame().reset_index()
data = data.rename(columns={'level_0':'Substrate',0:'Norm. K'})
data = data.merge(df_meta,left_on='Strain_ID',right_index=True)
data = data[data.Ribotype != 'Unknown']
data['Clade'] = data['Clade'].replace({1.:'1',2.:'2',3.:'3',4:'4',5:'5'})

# get ribotype colors
dict_ribotype_annot = read_csv(f"{dir_config}/ribotype_annotations.tsv")
dict_ribotype_annot = dict_ribotype_annot.to_dict()['Color'] 

# initialize figure
fig,axes = init_figure(
    figsize = [13,6],
    kwargs = {
        'nrows':1,
        'ncols':2,
        'sharey':False
        }
)

# plot each substrate in a sub-plot
for ax, substrate in zip(np.ravel(axes), substrates):

    # get data
    toplot = data[data.Substrate==substrate]

    sns.boxplot(
        data = toplot,
        x = "Clade",
        y = "Norm. K",
        order = ['1','2','3','4','5'], 
        width = 0.76,
        showfliers = False,
        showcaps = False,
        medianprops = {'lw':3,'color':'black'},
        whiskerprops = {'lw':1.5},
        boxprops={'lw':1.5,'facecolor':(0.9,0.9,0.9)},
        ax = ax
    )

    sns.swarmplot(
        data = toplot,
        x = "Clade",
        y = "Norm. K",
        hue = "Ribotype",
        palette = dict_ribotype_annot,
        order = ['1','2','3','4','5'], 
        size = 9,
        alpha = 1, 
        edgecolor = 'black',
        linewidth = 0.2,
        zorder = 10,
        dodge = False, 
        legend = False, 
        ax = ax
    )

    # set limits
    ax.set_xlim(-1,5)

    # set grid and spines
    ax.add_grid(lw=0.15)
    ax.thicken_spines(lw=2)
    ax.axhline(0,0,1,ls='--',lw=1,color='black')

    # set labels, title, and fontsize
    ylabel = substrate.split('-')[-1]
    ax.set_ylabel(f'OD(MM+{ylabel}) - OD(MM)', fontsize=15)
    ax.set_xlabel('Phylogenetic Clade', fontsize=15)
    ax.set_title(ylabel,fontsize=15,fontweight='bold',y=1.015)
    ax.enlarge_tick_labels(fontsize=15)

# adjust limits
axes[0].set_ylim([-0.025,0.725])
axes[1].set_ylim([-0.075,0.375])

# adjust spines and ticks
for ax in axes:
    ax.thicken_spines(lw=0,which='trb')
    ax.thicken_spines(lw=2,which='l')
    ax.axhline(0,0,1,lw=2,color='k')
    ax.tick_params(width=0)

# adjust whitespaces
plt.subplots_adjust(wspace=0.25,hspace=0.0)

# SAVE FIGURE
plt.savefig(f"{dir_figure}/supp/supp-figure-3-boxplots-norm-k-fructose-ribose.png",dpi=300,bbox_inches='tight')
plt.close()
