#!/usr/bin/env python

# import off-the-shelf packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# import off-the-shelf functions
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Rectangle
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import variation

# import custom library
from utils import *

# READ DATA
df_summ = read_csv("../../amiga-biolog/summary/merged_summary_norm_sub_by_median.txt")

# summary metric of interest: 
# subtraction-normalized carrying capacity in the untransformed lienar scale
param = 'norm(k_lin)'

# the following are variables needed for this figure
varbs = ['Strain_ID','Substrate',param]

# get median carrying capacity for each strain on each substrate
df_long = df_summ.loc[:,varbs]
df_long = df_long.groupby(['Strain_ID','Substrate'])
df_long = df_long.median(param).reset_index()

# pivot dataframe
df_wide = df_long.pivot(index='Substrate',columns='Strain_ID',values=param).T

# DESCRIBE MEDIAN CARRYING CAPACITY WITH HISTOGRAM

# initialize figure object
fig,ax = init_figure(figsize=[5,4])

# compute histogram
bins = np.arange(-0.4,1.0,0.1)
hist = np.ravel(df_long[param].values)

# plot histogram
ax.hist(x = hist,
        bins = bins,
        color = 'black',
        edgecolor = 'white',
        lw = 1.5,
        zorder=2)

ax.set_yscale('log')
ax.set_ylabel('Number of wells',fontsize=fontsize)
ax.set_xlabel('Normalized Carrying Capacity',fontsize=fontsize)
ax.yaxis.set_label_coords(-0.175,.5)
ax.xaxis.set_label_coords(0.5,-0.15)

# save histogram
plt.savefig(f"{dir_figure}/other/histogram_norm_k_lin.png")
plt.close()

# IDENTIFY SUBSTRATES OF INTEREST

# substrates that support growth
thresh_od_growth = 0.1 # OD must increase by at least 0.1
thresh_od_death = -0.05 # OD must decrease by at least 0.05
num_isolates = df_wide.shape[0] # number of isolates
thresh_hit = np.floor(0.1*num_isolates) # OD increase must occur for at least 10% of isolates

matches = lambda x: [1 if x >= thresh_od_growth else 0][0]
matches = df_wide.map(matches).sum(0).sort_values(ascending=False)
matches = list(matches[matches >= thresh_hit].index)

# substrates that inhibit growth
inhibitors = df_wide.T[df_wide.mean() < thresh_od_death].index.values

# CREATE BOXPLOTS OF NORMALIZED CARRYING CAPACITIES FOR SUBSTRATES OF INTEREST

# ~~~~~~~~~~~~ initialize plot ~~~~~~~~~~~~

fig,axes = init_figure(
    figsize = [7,14],
    kwargs = {
        'nrows':5,
        'ncols':3,
        'gridspec_kw':{
            'width_ratios':[13.5,3.0,0.9],
            'height_ratios':[1.5,0.75,10,1.5,3.0]
        }
    }
)

# create data-frame that maps substrates to group and color

# define plot lay-out
axl_pos = axes[2,0] # positive growth box plots
axc_pos = axes[2,1] # positive growth coefficients of variation
axr_pos = axes[2,2] # substrate color map

axl_neg = axes[4,0] # negative growth box plots
axc_neg = axes[4,1] # negative growth coefficients of variation

ax_legend = axes[0,0] # substrate color legend

# ~~~~~~~~~~~~ positive growth box plots ~~~~~~~~~~~~

# define data
toplot = df_wide.loc[:,matches]
sorted_substrates = toplot.median().sort_values(ascending=True).index.values

toplot_long = toplot.unstack().to_frame().reset_index()
toplot_long = toplot_long.rename(columns={0:param})

# define plot parameters
kwargs_stripplot = {
    's':4,
    'jitter':True,
    'facecolor':(0,0,0,0.65),
    'linewidth':0,
    'zorder':1
}

kwargs_boxplot = {
    'zorder':3,
    'widths':0.6,
    'whis':1.5,
    'showmeans':False,
    'showfliers':False,
    'showcaps':False,
    'boxprops':{
        'facecolor':(0,0,0,0.0),
        'edgecolor':'black',
        'linewidth':2
    },
    'whiskerprops':{
        'linewidth':2,
        'color':'black',
        'zorder':2
    },
    'medianprops':{
        'linewidth':5,
        'color':'blue',
        'zorder':2
    }
}


kwargs_plot = {
    'data':toplot_long,
    'x':param,
    'y':'Substrate',
    'order':sorted_substrates
}

# plot data
sns.boxplot(ax=axl_pos,**kwargs_plot,**kwargs_boxplot)
sns.stripplot(ax=axl_pos,**kwargs_plot,**kwargs_stripplot)

# adjust axes limits and ticks
axl_pos.set_xlim([-0.099,0.799])
axl_pos.set_ylim([-0.75,len(sorted_substrates)-0.25])
axl_pos.xaxis.set_major_locator(MultipleLocator(0.1))
axl_pos.tick_params(axis='both',which='both',length=0)

# adjust grid and remove spines
axl_pos.add_grid()
axl_pos.thicken_spines()

# adjust text and text label sizes
axl_pos.enlarge_tick_labels(fontsize=12)
axl_pos.set_xlabel(None)
axl_pos.set_ylabel(None)
axl_pos.text(
    x=axl_pos.get_center_xlim(),
    y=-2.5,
    s='Normalized Carrying Capacity',
    ha='center',
    va='center',
    fontsize=12,
    transform=axl_pos.transData
)

# ~~~~~~~~~~~~ positive growth coefficients of variation ~~~~~~~~~~~~

# compute CoVs
covs = toplot.apply(comp_covariance).loc[sorted_substrates]

# scale CoVs, for visualization purposes
scaler = MinMaxScaler(feature_range=(0.2,1.0))
covs_scaled = np.ravel(scaler.fit_transform(covs.values[:,np.newaxis]))
covs_colors = [(1-ii,1-ii,1-ii) for ii in covs_scaled]

axc_pos.barh(
    y=np.arange(len(covs)),
    width=covs,
    height=0.6,
    color=covs_colors,
    alpha=1,
    zorder=2
)

# adjust axes limits and ticks
axc_pos.set_xlim([0,1.1])
axc_pos.set_ylim([-0.75,len(covs_scaled)-0.25])
axc_pos.yaxis.set_major_locator(MultipleLocator(1))
plt.setp(axc_pos,xticks=[0,0.5,1])
plt.setp(axc_pos,yticklabels=[])
axc_pos.tick_params(axis='x',which='major',length=0)
axc_pos.tick_params(axis='y',which='major',color=(0,0,0,0.25),width=0.25,zorder=1)

# add grid and remove spines
axc_pos.add_grid()
axc_pos.thicken_spines()

# adjust text and text label size
axc_pos.enlarge_tick_labels(fontsize=12)
axc_pos.set_centered_xlabel(y=-2.5,s='Coef. of Variation')

# ~~~~~~~~~~~~ substrate color patches ~~~~~~~~~~~~

# read substrate --> {group, color} mapping
df_sub_mapping = read_csv(f"{dir_config}/substrate_annotations.tsv")
df_group_color = read_csv(f"{dir_config}/substrate_groups_colors.tsv")

# adda rectanglular patch for each substrate with appropriate color
for y_position, substrate in enumerate(sorted_substrates):
    
    box_group = df_sub_mapping.loc[substrate,'Figure_Label_Group']
    box_color = df_group_color.loc[box_group,'color']
    axr_pos.add_patch(
        Rectangle(
            xy = (0.05, y_position + 0.05 - 0.5),
            width = 0.9, 
            height = 0.9,
            facecolor = box_color,
            alpha = 0.8
        )
    )

# adjust liits, remove ticks and speines
axr_pos.set_xlim([0,1])
axr_pos.set_ylim([-0.75,len(sorted_substrates)-0.25])
plt.setp(axr_pos,xticks=[],yticks=[],xticklabels=[],yticklabels=[])
axr_pos.thicken_spines(lw=0)

# ~~~~~~~~~~~~ negative growth substrate colors legend ~~~~~~~~~~~~

# define raw data
toplot = df_wide.loc[:,inhibitors]
sorted_substrates = toplot.median().sort_values(ascending=True).index.values

toplot_long = toplot.unstack().to_frame().reset_index()
toplot_long = toplot_long.rename(columns={0:param})

# plot raw data

# unlike positive growth boxplot, median bars should be colored red
kwargs_boxplot['medianprops']['color'] = 'red'

kwargs_plot = {
    'data':toplot_long,
    'x':param,
    'y':'Substrate',
    'order':sorted_substrates
}

sns.boxplot(ax=axl_neg,**kwargs_plot,**kwargs_boxplot)
sns.stripplot(ax=axl_neg,**kwargs_plot,**kwargs_stripplot)

# adjust axes limits and ticks
axl_neg.set_xlim([-0.399,0.399])
axl_neg.set_ylim([-1,len(sorted_substrates)])
axl_neg.xaxis.set_major_locator(MultipleLocator(0.1))
axl_neg.tick_params(axis='both',which='both',length=0)

# adjust grid and remove spines
axl_neg.add_grid()
axl_neg.thicken_spines()

# adjust text and text label sizes
axl_neg.enlarge_tick_labels(fontsize=12)
axl_neg.set_xlabel(None)
axl_neg.set_ylabel(None)
axl_neg.text(
    x=axl_neg.get_center_xlim(),
    y=-2.5,
    s='Normalized Carrying Capacity',
    ha='center',
    va='center',
    fontsize=12,
    transform=axl_neg.transData
)

# ~~~~~~~~~~~~ negative growth coefficients of variation ~~~~~~~~~~~~

# compute CoVs
covs = toplot.apply(comp_covariance).loc[sorted_substrates]

# scale CoVs, for visualization purposes
scaler = MinMaxScaler(feature_range=(0.2,1.0))
covs_scaled = np.ravel(scaler.fit_transform(covs.values[:,np.newaxis]))

# assign shade to color black based on scaled CoVs
limit = lambda x: np.max([0,x]) # do not return negative numbers
covs_colors = [(limit(1-ii),limit(1-ii),limit(1-ii)) for ii in covs_scaled]

axc_neg.barh(
    y=np.arange(len(covs)),
    width=covs,
    height=0.6,
    color=covs_colors,
    alpha=1,
    zorder=2
)

# adjust axes limits and ticks
axc_neg.set_xlim([0,1.1])
axc_neg.set_ylim([-1,len(covs_scaled)])
axc_neg.yaxis.set_major_locator(MultipleLocator(1))
plt.setp(axc_neg,xticks=[0,0.5,1])
plt.setp(axc_neg,yticklabels=[])
axc_neg.tick_params(axis='x',which='major',length=0)
axc_neg.tick_params(axis='y',which='major',color=(0,0,0,0.25),width=0.25,zorder=1)

# add grid and remove spines
axc_neg.add_grid()
axc_neg.thicken_spines()

# adjust text and text label size
axc_neg.enlarge_tick_labels(fontsize=12)
axc_neg.set_centered_xlabel(y=-2.5,s='Coef. of Variation')

# ~~~~~~~~~~~~ legend for substrate group colors ~~~~~~~~~~~~

# define orders and positions of groups
sorted_subtypes = [
    'Sugar & Amino Alcohols',
    'Double & Triple Sugars',
    'Simple Sugars',
    'Carboxylic Acids',
    'Amino Acids',
    'Sugar Derivatives'
    ]

x_positions = [0,0,0,9,9,9]
y_positions = [0,1,2,0,1,2]

for _, (subtype, xpos, ypos) in enumerate(zip(sorted_subtypes, x_positions, y_positions)):
    
    # get color for substrate group
    color = df_group_color.loc[subtype,'color']
    
    # add color patch
    ax_legend.add_patch(
        Rectangle(
            xy = (xpos, ypos+0.05),
            facecolor = color,
            alpha = 0.8,
            width = 0.9,
            height = 0.9, 
            label = subtype
        )
    )
    
    # add label
    ax_legend.text(
        x = xpos + 1.15,
        y = ypos + 0.5,
        s = subtype,
        va = 'center',
        ha = 'left',
        fontsize = 12
    )

# adjust limits, remove ticks and spines
ax_legend.set_xlim([-0.5,16.5])
ax_legend.set_ylim([-0.5,3.5])
plt.setp(ax_legend,xticks=[],xticklabels=[])
plt.setp(ax_legend,yticks=[],yticklabels=[])
ax.thicken_spines(lw=2)

# ~~~~~~~~~~~~ finalize and save plot ~~~~~~~~~~~~

# create artificial whitespace by removing unused axes
for ax in axes[0,1:]: ax.remove()
for ax in axes[1,:]: ax.remove()
for ax in axes[3,:]: ax.remove()
axes[4,2].remove()

# save plot
plt.subplots_adjust(wspace=0.05,hspace=0.0)
plt.savefig(f"{dir_figure}/main/figure-1-boxplots-norm-k-lin.png",dpi=600,bbox_inches='tight')
plt.close()
