#!/usr/bin/env python

# import off-the-shelf packages
import pandas as pd
import matplotlib.pyplot as plt

# import off-the-shelf functions
from statsmodels.multivariate.pca import PCA
from matplotlib import projections
from adjustText import adjust_text

# import custom library
from utils import *

# READ DATA
df_summ = read_csv("../../amiga-biolog/summary/merged_summary_norm_sub_by_median.txt")

# summary metric of interest: 
# subtraction-normalized carrying capacity in the untransformed lienar scale
param = 'norm(k_lin)'

# the following are variables needed for this figure
varbs = ['Strain_ID','Ribotype','Clade','Substrate',param]

# get median carrying capacity for each strain on each substrate
df_long = df_summ.loc[:,varbs]
df_long = df_long.groupby(['Strain_ID','Substrate'])
df_long = df_long.median(param).reset_index()

# pivot dataframe
df_wide = df_long.pivot(index='Substrate',columns='Strain_ID',values=param).T

# IDENTIFY SUBSTRATES OF INTEREST

# substrates that support growth
thresh_od_growth = 0.1 # OD must increase by at least 0.1
thresh_od_death = -0.05 # OD must decrease by at least 0.05
thresh_hit = 9 # OD increase must occur for at least 9 isolates

matches = lambda x: [1 if x >= thresh_od_growth else 0][0]
matches = df_wide.map(matches).sum(0).sort_values(ascending=False)
matches = list(matches[matches >= thresh_hit].index)

# PREPARE DATA FOR PLOTTING

# raw data
pca_input = df_wide.loc[:,matches].dropna()

# PCA on covariance matrix
model = PCA(
    pca_input,
    ncomp = len(matches),
    demean = True,
    standardize = False
)

# position of samples on principal components
df_scores = model.scores 

# principal component loading factors
df_loadings = model.loadings

# proportion of variance explained by prinicipal components (%)
prop_var = (100 * (model.eigenvals/model.eigenvals.sum())).values

# ~~~~~~~~~~~~ initialize plot ~~~~~~~~~~~~

# initialize figure
projections.register_projection(customAxes)
fig = plt.figure(figsize = [13.3,15])

# define plot lay-out
grid = plt.GridSpec(
    nrows = 3,
    ncols = 5,
    wspace = 0.1,
    hspace = 0.35,
    width_ratios = [5.0, 1.5, 5.0, 0.2, 3.5],
) 

# map axes to subplot
subplot_kwargs = {'projection':'biolog'}

ax_00 = fig.add_subplot(grid[0,0],**subplot_kwargs)
ax_02 = fig.add_subplot(grid[0,2],**subplot_kwargs)

ax_10 = plt.subplot(grid[1,0],**subplot_kwargs)
ax_12 = plt.subplot(grid[1,2],**subplot_kwargs)

ax_20 = plt.subplot(grid[2,0],**subplot_kwargs)
ax_22 = plt.subplot(grid[2,2],**subplot_kwargs)

ax_legend_clades = plt.subplot(grid[0,4],**subplot_kwargs)
ax_legend_ribotypes = plt.subplot(grid[1:,4],**subplot_kwargs)

# ~~~~~~~~~~~~ PCA highlighting clade differences  ~~~~~~~~~~~~

# get typing information about each strain
df_typing = df_summ.loc[:,['Strain_ID','Ribotype','Clade']].drop_duplicates()

# handle unknowns
#df_typing = df_typing.replace({'RTUnk':'Unknown'})
df_typing['Ribotype'] = df_typing['Ribotype'].fillna('Unknown')
df_typing['Clade'] = df_typing['Clade'].replace({'Unknown':0})
#df_typing['Clade'] = df_typing['Clade'].fillna(0).astype(int).astype(str)
df_typing['Clade'] = df_typing['Clade'].astype(int).astype(str)
df_typing['Clade'] = df_typing['Clade'].replace({0:'Unknown'})
df_typing = df_typing.set_index('Strain_ID')

# merge typing information with PCA results
df_plot = df_scores.join(df_typing,how='left')
df_plot.to_csv('/Users/midani/Downloads/scores.tsv',sep='\t',header=True,index=True)

# get clade annotations (color and markers)
df_clade_annots = read_csv(f"{dir_config}/clade_annotations.tsv")
dict_clade_color  = df_clade_annots.loc[:,'Color'].to_dict()
dict_clade_color.pop('Unknown')

# plot
for ax, pcs in zip([ax_00,ax_02],[[0,1],[2,3]]):

    # define desired column headers for this sub-plot
    x_pos_pc = f"comp_0{pcs[0]}"
    y_pos_pc = f"comp_0{pcs[1]}"

    for ii, (clade, clade_color) in enumerate(dict_clade_color.items()):

        # get positions for isolates belonging to current clade
        df_plot_clade = df_plot[df_plot.Clade==clade]
        df_plot_clade = df_plot_clade.loc[:,[x_pos_pc,y_pos_pc]]

        # get centroid of all positions
        x_cent, y_cent = comp_centeroid(df_plot_clade.values)

        # plot centroid
        ax.scatter(
            x = x_cent,
            y = y_cent,
            color = clade_color,
            edgecolor = 'black',
            s = 50,
            zorder = 5
        )

        # define position of marker patch in legend
        x_pos = 0.00
        y_pos = 0.65 - 0.1 * ii

        # add marker to legend
        ax_legend_clades.scatter(
            x = x_pos,
            y = y_pos,
            s = 85,
            color = clade_color,
            edgecolor = 'black',
            zorder = 5
        )

        # add line to legend
        ax_legend_clades.plot(
            [x_pos - 3, x_pos + 3],
            [y_pos, y_pos],
            color = clade_color, 
            lw = 2,
            zorder = 1
        )

        # add label to legend
        ax_legend_clades.text(
            x = x_pos + 6,
            y = y_pos - 0.0025,
            s = f"Clade {clade}",
            fontsize = fontsize,
            va = 'center',
            ha = 'left'
        )

        # plot individual isolates
        for _, row in df_plot_clade.iterrows():

            x_pos, y_pos = row[x_pos_pc], row[y_pos_pc]

            ax.plot(
                [x_cent, x_pos],
                [y_cent, y_pos],
                color = clade_color,
                zorder = 1
            )

        # adjus ticks, spines, and grid
        ax.enlarge_tick_labels(fontsize=fontsize)
        ax.thicken_spines()
        ax.add_grid()

        # add axis labels
        base_label = "Principal Component {}\n({:.1f}% variance)"
        ax.set_xlabel(base_label.format(pcs[0]+1, prop_var[pcs[0]]),fontsize=fontsize)
        ax.set_ylabel(base_label.format(pcs[1]+1, prop_var[pcs[1]]),fontsize=fontsize)
    
    # adjust clade legend axis limits
    ax_00.set_xlim([-0.35,0.35])
    ax_00.set_ylim([-0.35,0.25])

    ax_02.set_xlim([-0.35,0.35])
    ax_02.set_ylim([-0.35,0.25])

    ax_legend_clades.set_xlim([-5,40])
    ax_legend_clades.set_ylim([0,0.9])

    # rmeove ticks and spines
    plt.setp(ax_legend_clades,xticks=[],yticks=[])
    ax_legend_clades.thicken_spines(lw=0)

# ~~~~~~~~~~~~ PCA highlighting ribotype differences  ~~~~~~~~~~~~

# identify ribotype and clade for each isolate
df_plot = df_scores.join(df_typing,how='left')

# add ribotype color
df_ribotype_annots = read_csv(f"{dir_config}/ribotype_annotations.tsv")
df_plot = pd.merge(
    left = df_plot,
    right = df_ribotype_annots,
    left_on = 'Ribotype',
    right_index = True
)

# add clade marker
df_clade_annots = read_csv(f"{dir_config}/clade_annotations.tsv")
df_plot = pd.merge(
    left = df_plot,
    right = df_clade_annots.loc[:,['Marker']],
    left_on = 'Clade',
    right_index = True
)

for ax, pcs in zip([ax_10,ax_12],[[0,1],[2,3]]):
    
    pc1 = df_plot.columns[pcs[0]]
    pc2 = df_plot.columns[pcs[1]]

    for isolate, row in df_plot.iterrows():

        ax.scatter(
            x = row[pc1],
            y = row[pc2],
            marker = row['Marker'],
            color = row['Color'],
            edgecolor = (1,1,1,0.5),
            alpha = 0.8,
            lw = 0.5,
            s = 200,
            zorder = 10
        )

        # adjus ticks, spines, and grid
        ax.enlarge_tick_labels(fontsize=fontsize)
        ax.thicken_spines()
        ax.add_grid()
        
        # add axis labels
        base_label = "Principal Component {}\n({:.1f}% variance)"
        ax.set_xlabel(base_label.format(pcs[0]+1, prop_var[pcs[0]]),fontsize=fontsize)
        ax.set_ylabel(base_label.format(pcs[1]+1, prop_var[pcs[1]]),fontsize=fontsize)

# adjust clade legend axis limits
ax_10.set_xlim([-0.35,0.35])
ax_10.set_ylim([-0.35,0.25])

ax_12.set_xlim([-0.25,0.35])
ax_12.set_ylim([-0.35,0.25])

# ~~~~~~~~~~~~ Legend for ribotypes grouped by clade  ~~~~~~~~~~~~

# get all necessary data
df_legend = df_plot.loc[:,['Clade','Ribotype','Color','Marker']].drop_duplicates()
df_legend = df_legend.reset_index(drop=True).sort_values(['Ribotype'],ascending=False)

# init parameters
offset = 0

for clade in ['5', '4', '3', '2', '1' ]:

    # get all related ribotypes
    df_ribotypes = df_legend[df_legend.Clade == clade].drop(['Clade'],axis=1)
    df_ribotypes = df_ribotypes.drop_duplicates().sort_values(['Ribotype'],ascending=False)
    df_ribotypes = df_ribotypes.reset_index(drop=True)

    for _, row in df_ribotypes.iterrows():

        # plot marker
        ax_legend_ribotypes.scatter(
            x = 0,
            y = offset,
            color = row['Color'],
            marker = row['Marker'],
            s = 85
        )

        # annotate marker
        ax_legend_ribotypes.annotate(
            text = row['Ribotype'],
            xy = (0, offset),
            xytext = (3, offset), 
            xycoords = 'data',
            ha = 'left',
            va = 'center',
            fontsize = 10.5 
        )

        offset += 1

    # annotate clade group
    ax_legend_ribotypes.text(
        x = 3, 
        y = offset,
        s = f"Clade {clade}",
        fontsize = 10.5,
        fontweight = 'bold'
    )

    offset += 2

ax_legend_ribotypes.set_xlim([-2.5, 40])
ax_legend_ribotypes.set_ylim([-1, offset-0.5])

plt.setp(ax_legend_ribotypes, xticks=[], yticks=[])
ax_legend_ribotypes.thicken_spines(lw=0.5)

# ~~~~~~~~~~~~ Contribution of Substrates to PCA Loading Factors  ~~~~~~~~~~~~

# get color and shorter names for substrates
df_substrate_annots = read_csv(f"{dir_config}/substrate_annotations.tsv")
dict_sub_labels = df_substrate_annots.loc[:,'Tidy_Name'].to_dict()
dict_sub_colors = df_substrate_annots.loc[:,'Color'].to_dict()

for ax, pcs in zip([ax_20,ax_22], [[0,1],[2,3]]):

    # identify substrates 
    #   for the top two principal components, rank substrates by contribution to loadings
    #   then get substrates that are one of the top 6 contributors for either of the components
    top_loadings = df_loadings.iloc[:,[pcs[0],pcs[1]]]
    top_loadings = top_loadings[(top_loadings.abs().rank(ascending=False) < 6).any(axis=1)]
    
    tmp = df_loadings.iloc[:,[0,1,2,3]]
    tmp = tmp.abs().rank(ascending=False)

    # get substrate colors and names
    dict_labels = [dict_sub_labels[substrate] for substrate in top_loadings.index.values]
    dict_colors = [dict_sub_colors[substrate] for substrate in top_loadings.index.values]

    texts, arrows = [], []

    for i in range(top_loadings.shape[0]):

        arrow_dx = top_loadings.iloc[i,0]
        arrow_dy = top_loadings.iloc[i,1]

        arrows.append(
            ax.arrow(
                0, 0,
                arrow_dx, arrow_dy, 
                width = top_loadings.values.max() / 75, 
                head_width = top_loadings.values.max() / 10, 
                color = dict_colors[i],
                alpha = 0.6, 
                linestyle = '-',
                length_includes_head = True, 
                zorder = 10
            )
        )

        texts.append(
            ax.text(
                arrow_dx,
                arrow_dy,
                f" {dict_labels[i]} ",
                color = dict_colors[i],
                ha = 'left',
                va = 'center',
                fontsize = fontsize,
                zorder = 10
            )
        )

    ax.set_xlim([-1.35,1.35])
    ax.set_ylim([-1.35,1.35])

    x = top_loadings.iloc[:,0].values
    y = top_loadings.iloc[:,1].values

    adjust_text(
        texts,
        x, y, 
        arrows,
        only_move = {
            'text':'y',
            'points':'y',
            'objects':'x'
        },
        time_lime = 15, 
        avoid_self = True,
        ax = ax,
        ensure_inside_axes = True
    )

    ax.enlarge_tick_labels(fontsize=fontsize)
    ax.thicken_spines()
    ax.add_grid()
        
    # add axis labels
    base_label = "Principal Component {}\n({:.1f}% variance)"
    ax.set_xlabel(base_label.format(pcs[0]+1, prop_var[pcs[0]]),fontsize=fontsize)
    ax.set_ylabel(base_label.format(pcs[1]+1, prop_var[pcs[1]]),fontsize=fontsize)

# save plot
plt.savefig(f"{dir_figure}/main/figure-2-pca-norm-k-lin.png",dpi=600,bbox_inches='tight')
plt.close()
