#!/usr/bin/env python

# import off-the-shelf packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# import off-the-shelf functions
from matplotlib.ticker import MultipleLocator

# import custom library
from utils import *

# ~~~~ DATA FOR TOP BIOLOG PANEL ~~~~~

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

# IDENTIFY SUBSTRATES OF INTEREST

# read enrichment analysis
df_results = read_csv('../../tables/strain_enrichment_analysis.tsv')

# focus on significant enrichment for emerging ribotypes
cond_1 = df_results['FDR q-val'] < 0.05
cond_2 = df_results.Term.isin(['RT255','RT023'])
df_results_substrates = df_results.loc[cond_1 & cond_2]
substrates = df_results_substrates.Substrate.unique()

# define data
data_biolog = df_wide.loc[:,substrates]
data_biolog = data_biolog.unstack().to_frame().reset_index()
data_biolog = data_biolog.rename(columns={'level_0':'Substrate',0:'Norm. K'})
data_biolog = data_biolog.merge(df_meta,left_on='Strain_ID',right_index=True)
data_biolog = data_biolog[data_biolog.Ribotype != 'Unknown']

# if it is not an emerging lineage, label it as "Other"
data_biolog['Ribotype'] = [ii if ii in ['RT023','RT255'] else 'Other' for ii in data_biolog.Ribotype.values]

# ~~~~ DATA FOR BOTTOM VALIDATION PANELS  ~~~~

# Read data
df_summ = read_csv("../../amiga-ribotype-255/summary/merged_summary.txt")

# handle NAs and data types
df_summ.Carbohydrate = df_summ.Carbohydrate.fillna('None')
df_summ.Concentration_mM = df_summ.Concentration_mM.fillna(0).astype(int)

# Parameters of interest
varbs = ['Isolate','Ribotype','Carbohydrate','Concentration_mM','Media']
metric = 'k_lin'
substrates = ['None','Fructose','Ribose']
label = 'Carrying Capacity'

# Reduce data
data_validation = df_summ.copy().loc[:, varbs + [metric]].dropna()

# compute medians for each unique condition of isolate x carb x conc. x media
data_validation = data_validation.groupby(varbs).median().reset_index()

# classify each isolate as RT255 or Other
data_validation['Group'] = [ii if ii=='RT255' else 'Other' for ii in data_validation['Ribotype'].values]
data_validation = data_validation.rename(columns={'Ribotype':'RT','Group':'Ribotype'})

# ~~~ INITIALIZE FIGURE ~~~~

# initialize figure
projections.register_projection(customAxes)
fig = plt.figure(figsize = [15,12.5])

# define plot lay-out
grid = plt.GridSpec(
    nrows = 3,
    ncols = 4,
    wspace = 0.1,
    hspace = 0.1,
    width_ratios = [1,1,2,2],
    height_ratios = [4,2,4]
) 

# map axes to subplot
subplot_kwargs = {'projection':'biolog'}

ax_top = fig.add_subplot(grid[0,:],**subplot_kwargs)

ax_b0 = plt.subplot(grid[2,0],**subplot_kwargs)
ax_b1 = plt.subplot(grid[2,1],**subplot_kwargs)
ax_b2 = plt.subplot(grid[2,2],**subplot_kwargs)
ax_b3 = plt.subplot(grid[2,3],**subplot_kwargs)

# ~~~~ ADD BIOLOG FIGURE ~~~~

# get order of strain sets to match Figure 3
df_medians = data_biolog.groupby(['Substrate'])['Norm. K'].median()
sorted_substrates = df_medians.sort_values(ascending=False).index.values

dict_palette = {
    'RT255':(0.5529,0.6275,0.7961),
    'RT023':(0.4,0.7608,0.6471),
    'Other':(0.9882,0.5529,0.3843)
}


sns.boxplot(
    data = data_biolog,
    x = "Substrate",
    y = "Norm. K",
    hue = "Ribotype",
    palette = dict_palette,#"Set2",
    hue_order = ["RT255","Other","RT023"],
    order = sorted_substrates,
    showfliers = False,
    showcaps = False,
    whis = [0,100],
    zorder = 5,
    ax = ax_top
)

sns.swarmplot(
    data = data_biolog,
    x = "Substrate",
    y = "Norm. K",
    hue = "Ribotype",
    hue_order = ["RT255","Other","RT023"],
    order = sorted_substrates,
    dodge = True,
    legend = False,
    size = 1.75,
    palette = "dark:k",
    zorder = 10,
    ax = ax_top
)

# adjust labels
ax_top.enlarge_tick_labels(fontsize=14)
[ii.set(rotation=45,ha='right') for ii in ax_top.get_xticklabels()]

# adjust limits, grid, and spines
ax_top.set_ylim([-0.05,0.75])
ax_top.thicken_spines(lw=2)
ax_top.add_grid(lw=0.1, zorder=1)
ax_top.set_axisbelow(True)

# add labels
ax_top.set_xlabel(None)
ax_top.set_ylabel('Normalized Carrying Capacity',fontsize=14)

# adjust legend
handles, labels = ax_top.get_legend_handles_labels()
ax_top.legend(ncols=3)

# COLOR SUBSTRATES BASED ON SIGNIFICANCE 

# get q-values
df_results_substrates = df_results_substrates.pivot(
    index = 'Substrate',
    columns = 'Term',
    values = 'FDR q-val'
)

# summarize: significance? yes or no ?
df_sig = df_results_substrates.map(lambda x: ~np.isnan(x))

dict_text_colors = {}
for substrate, is_sig in df_sig.iterrows():

    if not is_sig['RT255']:color = 'green'
    elif not is_sig['RT023']: color = 'slateblue'
    else: color = 'black'

    dict_text_colors[substrate] = color

[ii.set(color=dict_text_colors[ii.get_text()]) for ii in ax_top.get_xticklabels()];

# ~~~~ ADD VLIDATION FIGURES ~~~~

# Define criteria that select smples plotted in each subplot/axis
criteria_0 = {'Carbohydrate':['None'],'Media':'BHIS'}
criteria_1 = {'Carbohydrate':['None'],'Media':'CDMM'}
criteria_2 = {'Carbohydrate':['Fructose'],'Concentration_mM':[10,20,40],'Media':'CDMM'}
criteria_3 = {'Carbohydrate':['Ribose'],'Concentration_mM':[12,24,48],'Media':'CDMM'}

# Select samples
dict_toplot = {
    0: subsetDf(data_validation,criteria_0),
    1: subsetDf(data_validation,criteria_1),
    2: subsetDf(data_validation,criteria_2),
    3: subsetDf(data_validation,criteria_3)
}


# modify boxplot arguments to only plot median bar
kwargs_boxplot = {
    'palette':'dark:k',
    'medianprops':{'visible':True,'color':'k','ls':'-','lw':3},
    'whiskerprops':{'visible':False},
    'zorder':12,
    'showfliers':False,
    'showbox':False,
    'showcaps':False,
    'legend':False
}

kwargs_swarmplot = {
    's':6,
    'palette':{
        "RT255":"darkmagenta",
        "Other":"Gray"
    }
}

# # sub-plots for growth on rich and minimal media
for ii,ax in zip([0,1], [ax_b0,ax_b1]):   

    sns.boxplot(
        data = dict_toplot[ii],
        x = 'Ribotype',
        y = metric,
        hue = 'Ribotype',
        order = ['RT255','Other'],
        dodge = False,
        ax = ax,
        **kwargs_boxplot
    )
    
    sns.swarmplot(
        data = dict_toplot[ii],
        x = 'Ribotype',
        y = metric,
        hue = 'Ribotype',
        order = ['RT255','Other'],
        dodge = False,
        legend = False,
        ax = ax,
        **kwargs_swarmplot
    )

    ax.set_xlabel('Ribotype',fontsize=14)

# sub-plots for growth on minimal media supplemented with individual carbohydrates
for ii,ax in zip([2,3], [ax_b2,ax_b3]):
    
    
    sns.boxplot(
        data = dict_toplot[ii],
        x = 'Concentration_mM',
        y = metric,
        hue = 'Ribotype',
        hue_order = ['RT255','Other'],
        order = sorted(dict_toplot[ii]['Concentration_mM'].unique()),
        dodge = True,
        ax = ax,
        **kwargs_boxplot
    )
    
    sns.swarmplot(
        data = dict_toplot[ii],
        x = 'Concentration_mM',
        y = metric,
        hue = 'Ribotype',
        hue_order = ['RT255','Other'],
        order = sorted(dict_toplot[ii]['Concentration_mM'].unique()),
        dodge = True,
        legend = True,
        ax = ax,
        **kwargs_swarmplot
    )
    
    ax.set_xlabel('Concentration (mM)',fontsize=14)
    ax.legend(fontsize=12)

# adjust plot aesthetics
for ax in [ax_b0,ax_b1,ax_b2,ax_b3]:
    ax.enlarge_tick_labels(fontsize=14)
    ax.thicken_spines(lw=2)
    ax.set_ylim([-0.02,1.02])
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.grid(which='both',lw=0.2)
    ax.xaxis.grid(lw=0.2)
    ax.set_ylabel('Carrying Capacity',fontsize=14)

# only keep leftmost y-axis label and tick labels
for ax in [ax_b1,ax_b2,ax_b3]:
    plt.setp(ax,yticklabels=[])
    ax.set_ylabel(None)

# adjust titles
ax_b0.set_title('Rich Media',fontsize=fontsize,y=1.02)
ax_b1.set_title('Minimal Media',fontsize=fontsize,y=1.02)
ax_b2.set_title('Minimal Media + Fructose',fontsize=fontsize,y=1.02)
ax_b3.set_title('Minimal Media + Ribose',fontsize=fontsize,y=1.02)

# save plot
plt.savefig(f"{dir_figure}/main/figure-4-emerging-lineages.png",dpi=600,bbox_inches='tight')
plt.close()


    