#!/usr/bin/env python

# import off-the-shelf packages
import matplotlib.pyplot as plt
import seaborn as sns

# import off-the-shelf packages
from matplotlib.ticker import MultipleLocator

# import custom library
from utils import *

# define auxiliary function
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

dict_counts_all = {}
dict_counts_sig = {}
dict_df_results = {}
dict_col_labels = {}

fig, axes = plt.subplots(1,3,figsize=[7,3],width_ratios=[3,0.5,3])

for param,suffix,ax,title in zip(
    ['norm(k_lin)','norm(gr)'],
    ['k_lin','gr'],
    [axes[0],axes[2]],
    ['normalized carrying capacity','normalized growth rate']):

    # READ DATA
    df_summ = read_csv(f"../../amiga-biolog/summary/merged_summary_norm_sub_by_median_{suffix}.txt")

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
    df_results = read_csv(f"../../tables/strain_enrichment_analysis_norm_{suffix}.tsv")
    df_results = df_results.replace({'LabAdapted':'Lab-adapted'})
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
    [ax.spines[ii].set(lw=2) for ii in ['top','bottom','right','left']]
    ax.axhline(0,color='black',lw=2)
    ax.yaxis.grid(lw=0.25,zorder=1)

    # modify fonts and labels
    fontsize = 10
    [ii.set_fontsize(fontsize) for ii in ax.get_xticklabels()+ax.get_yticklabels()]
    ax.text(-3.75,15,'Positive',rotation=90,va='center',fontsize=fontsize)
    ax.text(-3.,15,'Counts',rotation=90,va='center',fontsize=fontsize)
    ax.text(-3.75,-15,'Negative',rotation=90,va='center',fontsize=fontsize)
    ax.text(-3.,-15,'Counts',rotation=90,va='center',fontsize=fontsize)
    ax.set_title(f'Enrichment scores\n({title})',fontsize=fontsize,y=1.05)
    ax.set_xlabel(None)
    ax.set_ylabel(None)

    # add legend
    ax.legend(fontsize=9)

# SAVE FIGURE
axes[1].remove()
plt.savefig(f"{dir_figure}/supp/supp-figure-4-strain-enrichmenty-bar-plots.png",dpi=600,bbox_inches='tight')
plt.close()