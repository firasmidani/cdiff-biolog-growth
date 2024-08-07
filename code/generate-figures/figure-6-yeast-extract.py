#!/usr/bin/env python

# import off-the-shelf packages
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# import custom library
from utils import *
   
## BIOLOG FIGURES

# read biolog data
df_biolog = read_csv("../../amiga-yeast-extract-biolog/summary/merged_summary_norm_div_by_median.txt")

# indicate low vs high yeast extract based on plate ID
get_extract_conc = lambda ii: ii.split('_')[1].replace('YT1','').replace('YT','')
df_conc = [get_extract_conc(ii) for ii in df_biolog.Isolate.values]
df_conc = pd.DataFrame(
    data = df_conc,
    index = df_biolog.index.values,
    columns = ['YeastExtract']
)

# update summary dataframe with extract concentration
df_biolog = df_biolog.join(df_conc).reset_index()

# focus on susbtrates tha yield growth of at least 20% 
ls_subs = df_biolog[df_biolog['norm(k_lin)'] > 1.2].Substrate.unique()
toplot = df_biolog[df_biolog.Substrate.isin(ls_subs)]
toplot = toplot.pivot(index='Substrate',values='norm(k_lin)',columns='YeastExtract')

# identify substrates that may be consumed differntly based on yeast extract cocnentration
ls_subs_diff_a = toplot[(toplot['High']/toplot['Low'])>1.1].index.values
ls_subs_diff_b = toplot[(toplot['High']/toplot['Low'])<0.9].index.values
ls_subs_diff = list(set(ls_subs_diff_a).union(ls_subs_diff_b))
ls_subs_other = list(set(ls_subs).difference(set(ls_subs_diff)))

# cluster substrates based on normalized carrying capacities
ls_subs_diff = cluster_rows(toplot.loc[ls_subs_diff,:])
ls_subs_other = cluster_rows(toplot.loc[ls_subs_other,:])
ratio_subs = float(len(ls_subs_other))/len(ls_subs_diff)

# heatmap for substrates that are impact by yeast extract concentration

fig,ax = init_figure(figsize = [4,4])

# initialize figure grid
fig, axes = init_figure(
    figsize = [3,10],
    kwargs = dict(
        nrows = 2,
        ncols = 1,
        height_ratios = [1.2,2]
    )
)

# heatmap for substrates that are impact by yeast extract concentration
h = sns.heatmap(
    data = toplot.loc[ls_subs_diff,['Low','High']],
    ax = axes[0],
    vmin = 0,
    center = 1,
    vmax = 3,
    annot = True,
    fmt = '.2f'
)

# remaining substrates that can promote growth of C. diff by at least 1.1 fold-change
sns.heatmap(
    data = toplot.loc[ls_subs_other,['Low','High']],
    ax = axes[1],
    vmin = 0,
    center = 1,
    vmax = 4,
    annot = True,
    fmt = '.2f'
)

# modify labels, title, and ticks
for ax in axes:
    ax.set_title('Carrying Capacity Fold-Change',fontsize=12,y=1.04)
    ax.set_xlabel('Concentration of Yeast Extract',fontsize=12,labelpad=8.0)
    ax.set_ylabel(None)
    ax.enlarge_tick_labels(fontsize=12)

# SAVE FIGURE
plt.subplots_adjust(hspace=0.4)
plt.savefig(f"{dir_figure}/main/figure-6-yeast-extract.png",dpi=300,bbox_inches='tight')
plt.close()
