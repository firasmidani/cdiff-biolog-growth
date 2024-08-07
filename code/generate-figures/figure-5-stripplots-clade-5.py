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
df_summ = read_csv("../../amiga-clade-5/summary/merged_summary_norm_sub.txt")

# get strain typing data
df_meta = df_summ.loc[:,['Isolate','Ribotype','Clade']].drop_duplicates()
df_meta = df_meta.replace({'RTUnk':'Unknown'}).fillna('Unknown')
df_meta = df_meta.set_index('Isolate')

# convert clades calls to str (whole numbers only) 
def numeric_to_str(ii): 
    if isinstance(ii,float):  return str(int(ii))
    else: return str(ii)
df_meta['Clade'] = df_meta['Clade'].apply(lambda x: numeric_to_str(x))

# Parameter of Interest
metric = 'norm(k_lin)'

# get clade annotations (color and markers)
df_clade_annots = read_csv(f"{dir_config}/clade_annotations.tsv")
dict_markers = df_clade_annots['Marker'].to_dict()

# get ribotype annotations
df_ribotype_annots = read_csv(f"{dir_config}/ribotype_annotations.tsv")
dict_colors = df_ribotype_annots['Color'].to_dict()

# initialize figure grid
fig, axes = init_figure(
    figsize = [10,6],
    kwargs = dict(
        nrows = 2,
        ncols = 4,
        sharey = False
    )
)

# create temporary collections or markers needed to adjust marker shapes in plot
dict_path = {}
for ii in dict_markers.values():
    s = axes[0,0].scatter([1,2],[3,4],marker=ii)
    dict_path[ii] = s.get_paths()[0]
    s.remove()

order_subs = ['Glucose',
              'Fructose',
              'Tagatose',
              'N-acetylneuraminic acid',
              'N-acetylglucosamine',
              'Mannitol',
              'Salicin',
              'Ribose']

for ax, substrate in zip(np.ravel(axes),order_subs):

    # select substrate 
    df_plot = subsetDf(df_summ,{'Carbon Source':substrate})

    # get medians for each condition (isolate x substrate)
    df_plot = df_plot.groupby(['Isolate'])[metric].median().to_frame()

    # re-add isolate typing data
    df_plot = df_plot.join(df_meta)

    # list of ribotypes and ribotypes
    ls_ribotypes = df_plot.Ribotype.values
    ls_clades = df_plot.Clade.values

    # group isolates into 'Clade 5' or 'Other'
    df_plot['Group'] = ['Clade 5' if clade == '5' else 'Other' for clade in ls_clades]

    order = df_plot.groupby(['Ribotype'])[metric].median().sort_values().index.values

    sns.boxplot(
        data = df_plot,
        y = metric,
        x = 'Group',
        order = ['Clade 5','Other'],
        color = 'black',
        medianprops = {'visible':True, 'color':'black', 'ls':'-', 'lw':2},
        whiskerprops = {'visible':False},
        zorder = 10,
        showfliers = False,
        showbox = False, 
        showcaps = False, 
        legend = False,
        ax = ax
    )

    sns.swarmplot(
        data = df_plot,
        y = metric, 
        x = 'Group',
        order = ['Clade 5','Other'],
        hue = 'Ribotype',
        size = 7.5,
        palette = df_ribotype_annots.Color.to_dict(),
        legend = False,
        ax = ax
    )

    # next few steps: replace default circle marker with clade-specific markers

    # split isolates into groups
    g0 = df_plot[df_plot.Group=='Clade 5'].sort_index(ascending=True)
    g1 = df_plot[df_plot.Group=='Other'].sort_index(ascending=True)

    # get marker for isolates
    g0_markers = g0.apply(lambda x: dict_markers[x['Clade']],axis=1)
    g1_markers = g1.apply(lambda x: dict_markers[x['Clade']],axis=1)
    
    g0 = g0.join(g0_markers.to_frame().rename(columns={0:'marker'}))
    g1 = g1.join(g1_markers.to_frame().rename(columns={0:'marker'}))
    
    # modify markers (paths are ordered by row index and so are g0/g1 dataframes)
    ax.collections[0].set_paths([dict_path[ii] for ii in g0.marker.values])
    ax.collections[1].set_paths([dict_path[ii] for ii in g1.marker.values])

    # adjust tick labels
    ax.enlarge_tick_labels(fontsize=12)
    [ii.set(rotation=0) for ii in ax.get_xticklabels()]
    ax.yaxis.set_major_locator(MultipleLocator(0.1))

    # adjust spines and grid
    ax.thicken_spines(lw=2)
    ax.add_grid(lw=0.2)
    ax.set_ylim([-0.05,0.55])

    # adjust labels and title
    ax.set_ylabel('Norm. Carrying Capacity',fontsize=12)
    ax.set_xlabel(None)
    ax.set_title(f'{substrate}',fontsize=12,y=1.025)

# adjust figure whitespace
plt.subplots_adjust(hspace=0.5, wspace=0.7)

# SAVE FIGURE
plt.savefig(f"{dir_figure}/main/figure-5-stripplots-clade-5.png",dpi=300,bbox_inches='tight')
plt.close()
