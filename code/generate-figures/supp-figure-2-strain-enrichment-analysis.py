#!/usr/bin/env python

# import off-the-shelf packages
import matplotlib.pyplot as plt
import seaborn as sns

# import off-the-shelf packages
from matplotlib.ticker import MultipleLocator

# import custom library
from utils import *

# READ DATA
df_summ = read_csv("../../amiga-biolog/summary/merged_summary_norm_sub_by_median.txt")

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

# READ RESULTS OF  ENRICHMENT ANALYSIS
df_results = read_csv('../../tables/strain_enrichment_analysis.tsv')

# COUNT ENRICHMENTS FOR EACH STRAIN SET

# distinguish positive vs negative enrichments
df_results['Sign'] = df_results.apply(lambda x: [1 if x['NES'] > 0 else -1][0],axis=1)

df_results_all = df_results.copy()
df_results_sig = df_results[df_results['FDR q-val'] < 0.05]

counts_all = df_results_all.groupby(['Term','Sign'])['Substrate'].count().reset_index()
counts_sig = df_results_sig.groupby(['Term','Sign'])['Substrate'].count().reset_index()

# get order of strain sets to match Figure 3
df_pivot = df_results.pivot(index='Substrate',columns='Term',values ='NES').astype(float)
g = sns.clustermap(data = df_pivot)
_, col_labels = get_clustermap_labels(g)
plt.close()

def add_bar_plot(data,color):
    sns.barplot(
        data = data,
        x = 'Term',
        y = 'Substrate',
        order = col_labels,
        color = color,
        zorder = 2,
        ax = ax
    )

fig, ax = init_figure(figsize=[3,3])

# visualize counts of all positive enrichments
to_plot = counts_all[counts_all.Sign > 0]
add_bar_plot(to_plot,rgb([252,218,200]))
         
# visualize counts of all negative enrichments
to_plot = counts_all[counts_all.Sign < 0]
to_plot.loc[:,'Substrate'] = -1 * to_plot['Substrate']
add_bar_plot(to_plot,rgb([204,226,238]))

# visualize counts of significant positive enrichments
to_plot = counts_sig[counts_sig.Sign > 0]
add_bar_plot(to_plot,rgb([118,12,32]))

# visualize counts of significant negative enrichments
to_plot = counts_sig[counts_sig.Sign < 0]
to_plot.loc[:,'Substrate'] = -1 * to_plot['Substrate']
add_bar_plot(to_plot,rgb([16,56,109]))

xlim = ax.get_xlim()[1]

# explicitly add legend markers (set outside of final axes bounds)
ax.bar(xlim+2,1,color='black',label='Significant')
ax.bar(xlim+2,1,color='black',alpha=0.1,label='NS')

# rotate x-axis tick labels
[ii.set(rotation=90) for ii in ax.get_xticklabels()]

# adjust axis limits
ax.set_xlim([-1,xlim+0.5])
ax.set_ylim([-30,30])

# modify y-axis tick-labels
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_formatter(lambda x, pos: int(abs(x)))

# modify spines and grid
ax.thicken_spines(lw=2)
ax.axhline(0,color='black',lw=2)
ax.yaxis.grid(lw=0.25,zorder=1)

# modify fonts and labels
fontsize = 10
ax.enlarge_tick_labels(fontsize=fontsize)
ax.text(-3.75,15,'Positive',rotation=90,va='center',fontsize=fontsize)
ax.text(-3.,15,'Counts',rotation=90,va='center',fontsize=fontsize)
ax.text(-3.75,-15,'Negative',rotation=90,va='center',fontsize=fontsize)
ax.text(-3.,-15,'Counts',rotation=90,va='center',fontsize=fontsize)
ax.set_title('Enrichment Scores',fontsize=fontsize)
ax.set_xlabel(None)
ax.set_ylabel(None)

# add legend
ax.legend(fontsize=9)

# SAVE FIGURE
plt.savefig(f"{dir_figure}/supp/supp-figure-2-strain-enrichment-heatmap.png",dpi=300,bbox_inches='tight')
plt.close()
