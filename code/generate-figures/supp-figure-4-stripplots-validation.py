#!/usr/bin/env python

# import off-the-shelf packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# import off-the-shelf packages
from matplotlib.ticker import MultipleLocator

# import custom library
from utils import *

# READ DATA
df_summ = read_csv("../../amiga-validation/summary/merged_summary_norm_sub.txt")

# Parameter of Interest
metric = 'norm(k_lin)'

# get strain typing data
df_meta = df_summ.loc[:,['Isolate','Ribotype']].drop_duplicates()
df_meta = df_meta.replace({'RTUnk':'Unknown'}).fillna('Unknown')
df_meta = df_meta.set_index('Isolate')

# get ribotype colors
df_ribotype_annots = read_csv(f"{dir_config}/ribotype_annotations.tsv")
dict_colors= df_ribotype_annots.Color.to_dict()

# initialize figure grid
fig, axes = init_figure(
    figsize = [7,7],
    kwargs = dict(
        nrows = 2,
        ncols = 2,
        sharey = False
    )
)

substrates = ['Fructose','Mannitol','Salicin','Ribose']
for ax, substrate in zip(np.ravel(axes),substrates):

    # select substrate 
    df_plot = subsetDf(df_summ,{'Carbon Source':substrate})

    # get medians for each condition (isolate x substrate)
    df_plot = df_plot[~df_plot.Ribotype.isna()]
    df_plot = df_plot.groupby(['Isolate'])[metric].median().to_frame()

    # re-add isolate typing data
    df_plot = df_plot.join(df_meta)

    # get substrate concentration
    conc = list(set(df_summ[df_summ['Carbon Source']==substrate].Carbon_mM))[0]

    # order ribotypes by median growth
    order = df_plot.groupby(['Ribotype']).median().sort_values(metric).index.values

    # only plot median bar from boxplot
    sns.boxplot(
        data = df_plot,
        y = metric,
        x = 'Ribotype',
        order = order,
        color = 'black',
        zorder = 10,
        showbox = False,
        showfliers = False,
        showcaps = False,
        medianprops = {'visible':True, 'color':'black', 'ls':'-', 'lw':2},
        whiskerprops = {'visible':False},
        ax = ax
    )

    sns.stripplot(
        data = df_plot,
        y = metric,
        x = 'Ribotype',
        hue = 'Ribotype',
        palette = dict_colors,
        order = order,
        size = 7.5, 
        jitter = False,
        ax = ax
    )

    # adjust tick labels
    ax.enlarge_tick_labels(fontsize=12)
    [ii.set(rotation=90) for ii in ax.get_xticklabels()]

    # adjust spines, grid, legend, and axes
    ax.thicken_spines(lw=2)
    ax.add_grid(lw=0.2)
    ax.set_ylim([0.0,0.4])
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))

    # adjust labels
    ax.set_xlabel(None)
    ax.set_ylabel('Norm. Carrying Capacity', fontsize=12)
    ax.set_title(f"{substrate} ({conc} mM)",fontsize=12, y=1.05)

# adjust figure whitespace
plt.subplots_adjust(hspace=0.7, wspace=0.5)

# SAVE FIGURE
plt.savefig(f"{dir_figure}/supp/supp-figure-4-stripplots-validation.png",dpi=300,bbox_inches='tight')
plt.close()
