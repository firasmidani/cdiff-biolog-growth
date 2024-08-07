#!/usr/bin/env python

# import off-the-shelf packages
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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

# IDENTIFY SUBSTRATES OF INTEREST

# substrates that support growth
thresh_od_growth = 0.1 # OD must increase by at least 0.1
thresh_od_death = -0.05 # OD must decrease by at least 0.05
num_isolates = df_wide.shape[0] # number of isolates
thresh_hit = np.floor(0.1*num_isolates) # OD increase must occur for at least 10% of isolates

matches = lambda x: [1 if x >= thresh_od_growth else 0][0]
matches = df_wide.map(matches).sum(0).sort_values(ascending=False)
matches = list(matches[matches >= thresh_hit].index)

# order substrates and clades
df_medians = df_wide.loc[:,matches].median()
sorted_substrates = df_medians.sort_values(ascending=False).index.values
sorted_clades = [1,2,3,4,5]

# define data
data = df_wide.loc[:,sorted_substrates]
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
    figsize = [19,25],
    kwargs = {
        'nrows':6,
        'ncols':5,
        'sharey':False
        }
)

for ax, substrate in zip(np.ravel(axes), sorted_substrates):

    toplot = data[data.Substrate==substrate]

    sns.boxplot(
        data = toplot,
        x = "Clade",
        y = "Norm. K",
        order = sorted_clades, 
        width = 0.76,
        showfliers = False,
        showcaps = False,
        medianprops = {'lw':3,'color':'black'},
        whiskerprops = {'lw':1.5},
        boxprops={'lw':1.5,'facecolor':(0.9,0.9,0.9)},
        ax = ax
    )

    sns.stripplot(
        data = toplot,
        x = "Clade",
        y = "Norm. K",
        hue = "Ribotype",
        palette = dict_ribotype_annot,
        order = sorted_clades,
        s = 6,
        alpha = 0.8, 
        edgecolor = 'black',
        linewidth = 0.1,
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
    ax.set_title(substrate,fontsize=15,fontweight='bold',y=1.015)
    ax.set_xlabel('Phylogenetic Clade', fontsize=15)
    ax.set_ylabel('Norm. K', fontsize=15)
    ax.enlarge_tick_labels(fontsize=15)

# delete extra axes
[fig.delaxes(axes[-1,ii]) for ii in [1,2,3,4]]

# adjust whitespaces
plt.subplots_adjust(wspace=0.45,hspace=0.7)

# SAVE FIGURE
plt.savefig(f"{dir_figure}/supp/supp-figure-1-boxplots-norm-k.png",dpi=600,bbox_inches='tight')
plt.close()
