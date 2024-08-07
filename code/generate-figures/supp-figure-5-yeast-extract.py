#!/usr/bin/env python

# import off-the-shelf packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# import off-the-shelf packages
from matplotlib.ticker import MultipleLocator

# import custom library
from utils import *
   
# define main plotting function
def grouped_bar_plot(ax,df,metric="k_lin"):

    fontsize = 12
    linewidth = 2

    df = df.replace({'None':'MM','Trehalose':'MM + Trehalose'})    

    kwargs_data = dict(
        ax = ax,
        data = df,
        x = "Substrate",
        y = metric,
        hue = "Dilution",
        order = ["MM","MM + Trehalose"],
    )

    kwargs_barplot = {
        'zorder':1,
        'width':0.8,
        'lw':0.1,
        'edgecolor':'black',
        'capsize':0,
        'estimator':'median',
        'errorbar':None,
        'err_kws':dict(
            linewidth = 3,
            color = 'black',
            zorder = 3
        ),
        'palette':{
            '1:10':'#d94701',
            '1:20':'#fd8d3c',
            '1:40':'#fdbe85',
            '1:80':'#feedde'
        }
    }

    kwargs_swarmplot = {
        'zorder':5,
        'dodge':True,
        'legend':False,
        'size':3.5,
        'linewidth':0.5,
        'edgecolor':'black',
        'palette':'light:white'
    }

    kwargs_stripplot = {
        'zorder':5,
        'dodge':True,
        'legend':False,
        'jitter':True,
        'size':3.5,
        'linewidth':0.5,
        'edgecolor':'black',
        'palette':'light:white'
    }

    sns.barplot(**kwargs_data,**kwargs_barplot)
    sns.stripplot(**kwargs_data,**kwargs_stripplot)
    #sns.swarmplot(**kwargs_data,**kwargs_swarmplot)

    # adjust tick labels
    ax.tick_params(width=linewidth,length=5)
    ax.enlarge_tick_labels(fontsize=fontsize)

    # adjust spines and grid
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.thicken_spines(lw=linewidth)
    ax.yaxis.grid(which='both',lw=0.1)
    ax.set_axisbelow(True)

    # adjust labels
    ax.set_xlabel(None,fontsize=fontsize)

    # adjust legend
    ax.legend(loc='upper left')

    return ax


# ~~~ plot additional validation data ~~~

# initialize figure grid
fig, axes = init_figure(
    figsize = [9,9],
    kwargs = dict(
        nrows = 2,
        ncols = 2,
        sharey = True
    )
)

dir_data = "../../amiga-yeast-extract-validation/data-endpoints-only"
dir_maps = "../../amiga-yeast-extract-validation/mapping"

df_m68 = []
for filename in ['11062021-001','11062021-002']:

    # read files
    tmp_map = pd.read_csv(f"{dir_maps}/{filename}.txt",sep='\t',header=0,index_col=0)
    tmp_dat = pd.read_csv(f"{dir_data}/{filename}.txt",sep='\t',header=0,index_col=0)

    # format plate layout to tidy layout
    tmp_dat = tmp_dat.unstack().dropna().to_frame().reset_index()
    tmp_dat = tmp_dat.rename(columns={'level_0':'Column','<>':'Row',0:'OD600'})

    # create Well ID which is necessary for merging with meta-data
    tmp_dat['Well_ID'] = tmp_dat.apply(lambda x: x['Row']+x['Column'],axis=1)
    tmp_dat = tmp_dat.set_index('Well_ID').drop(['Column','Row'],axis=1)

    # merge and handle missing values
    tmp_dat = tmp_dat.join(tmp_map)
    tmp_dat['Substrate'] = tmp_dat['Substrate'].fillna('None')
    tmp_dat['Concentration'] = tmp_dat['Concentration'].fillna('None')

    df_m68.append(tmp_dat)

# concat plate data
df_m68 = pd.concat(df_m68)
df_m68 = df_m68[df_m68.Substrate.isin(['None','Trehalose'])]

cond_1 = subsetDf(df_m68,{'Substrate':'Trehalose','Concentration':40})
cond_2 = subsetDf(df_m68,{'Substrate':'None'})
df_plot = pd.concat([cond_1,cond_2])

# plot growth using unwashed/washed inocula
grouped_bar_plot(
    ax = axes[0,0],
    df = subsetDf(df_plot,{'Isolate':'PRB827','Washing':'No'}),
    metric = "OD600"
)
grouped_bar_plot(
    ax = axes[1,0],
    df = subsetDf(df_plot,{'Isolate':'PRB827','Washing':'Yes'}),
    metric = "OD600"
)

# ~~~ read CD2015 (RT027) data ~~~
dir_data = "../../amiga-yeast-extract-validation/data"
dir_maps = "../../amiga-yeast-extract-validation/mapping"

df_cd2015 = []
for filename in ['11062021-003','11062021-004']:

    # read files
    tmp_map = pd.read_csv(f"{dir_maps}/{filename}.txt",sep='\t',header=0,index_col=0)
    tmp_dat = pd.read_csv(f"{dir_data}/{filename}.asc",sep='\t',header=None,index_col=None,encoding='UTF-16',engine='python')

    # grab column of interest at ~17 hours
    col = int(17*3600/600)
    tmp_dat = tmp_dat.iloc[:,[0,col]]

    # format plate layout to tidy layout
    tmp_dat = tmp_dat.rename(columns={0:'Well_ID',col:'OD600'})
    tmp_dat = tmp_dat.set_index('Well_ID')

    # merge and handle missing values
    tmp_dat = tmp_dat.join(tmp_map)
    tmp_dat['Substrate'] = tmp_dat['Substrate'].fillna('None')
    tmp_dat['Concentration'] = tmp_dat['Concentration'].fillna('None')

    df_cd2015.append(tmp_dat)

# concat plate data
df_cd2015 = pd.concat(df_cd2015)
df_cd2015 = df_cd2015[df_cd2015.Substrate.isin(['None','Trehalose'])]

cond_1 = subsetDf(df_cd2015,{'Substrate':'Trehalose','Concentration':40})
cond_2 = subsetDf(df_cd2015,{'Substrate':'None'})

df_plot = pd.concat([cond_1,cond_2])

# plot growth using unwashed/washed inocula
grouped_bar_plot(
    ax = axes[0,1],
    df = subsetDf(df_plot,{'Isolate':'PRB268','Washing':'No'}),
    metric = "OD600"
)
grouped_bar_plot(
    ax = axes[1,1],
    df = subsetDf(df_plot,{'Isolate':'PRB268','Washing':'Yes'}),
    metric = "OD600"
)

# adjust axes and legend title
for ax in np.ravel(axes):
    ax.set_ylim([0,1.0])
    ax.legend(title='Dilution',loc=2)

# adjust labels and titles
axes[0,0].set_title(f"M68 (RT017) unwashed inocula",fontsize=12,y=1.02)
axes[1,0].set_title(f"M68 (RT017) washed inocula",fontsize=12,y=1.02)

axes[0,1].set_title(f"CD2015 (RT027) unwashed inocula",fontsize=12,y=1.02)
axes[1,1].set_title(f"CD2015 (RT027) washed inocula",fontsize=12,y=1.02)

axes[0,0].set_ylabel(f"OD600 at 17 hours",fontsize=12)
axes[1,0].set_ylabel(f"OD600 at 17 hours",fontsize=12)

# SAVE FIGURE
plt.subplots_adjust(hspace=0.3,wspace=0.1)
plt.savefig(f"{dir_figure}/supp/supp-figure-5-yeast-extract.png",dpi=300,bbox_inches='tight')
plt.close()
