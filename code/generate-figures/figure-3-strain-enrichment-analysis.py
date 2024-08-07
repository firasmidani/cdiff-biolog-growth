#!/usr/bin/env python

# import off-the-shelf packages
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
df_meta = df_meta.set_index('Strain_ID')

# READ RESULTS OF  ENRICHMENT ANALYSIS
df_results = read_csv('../../tables/strain_enrichment_analysis.tsv')

# PLOT ENRICHMENT ANALYSIS HEATMAP

# define order of substrates in heatmap
substrates = df_results.Substrate.unique()
df_medians = df_wide.loc[:,substrates].median()
substrate_order = df_medians.sort_values(ascending=False).index.values

# tabulate enrichment scores (substrates x strain sets)
df_pivot = df_results.pivot(
    index = 'Substrate',
    columns = 'Term',
    values = 'NES'
)
df_pivot = df_pivot.astype(float)
df_pivot = df_pivot.loc[substrate_order,:]

# define heatmap color scale bounds
vmax = df_pivot.abs().max().max()
vmin = -1 * vmax

g = sns.clustermap(
    data = df_pivot,
    center = 0,
    vmin = vmin, 
    vmax = vmax,
    col_cluster = True,
    row_cluster = False, 
    tree_kws = dict(linewidth=0),
    cbar_kws = dict(orientation='horizontal'),
    cmap = 'RdBu_r',
    figsize=[5.5,10.5]
)

# ADD SIGNIFICANCE MARKERS

# get row and column lables
row_labels, col_labels = get_clustermap_labels(g)

# get q-values
df_qvals = df_results.pivot(
    index = 'Substrate',
    columns = 'Term',
    values = 'FDR q-val'
)
# get significance-based annotations
df_annots = df_qvals.loc[row_labels, col_labels].map(lambda x: get_sig_annot(x))
df_annots = df_annots.reset_index(drop=True).T.reset_index(drop=True).T
df_annots = df_annots.unstack().dropna().to_frame()

for (col,row), annot in df_annots.iterrows():
    text = g.ax_heatmap.text(
        x = col + 0.5, 
        y = row + 0.75,
        s = annot[0],
        ha = 'center',
        va = 'center', 
        color = 'white',
        fontsize =15
    )

# MODIFY HEATMAP AESTHETICS

# adjust fontsize
for ax in [g.ax_heatmap,g.ax_cbar]:
    [ii.set_fontsize(12) for ii in ax.get_xticklabels()+ax.get_yticklabels()]

# lay color bar horizontally on top of the heatmap
pos_col_dend = g.ax_col_dendrogram.get_position() # position of dendrogram, if plotted
x0 = pos_col_dend.x0
y0 = pos_col_dend.y0
w = pos_col_dend.width
h = pos_col_dend.height

# place left corner of color near left corner of dendogram
# grant color bar same width as dendrogram but very small height
g.ax_cbar.set_position([x0, y0*1.05, w, 0.03])

# add titles and labels
g.ax_cbar.set_title('Normalized Enrichment Score', fontsize=12, y=1.05)
g.ax_heatmap.set_xlabel(None)

# adjust spines
kwargs_spine = {'color':'k','lw':4,'zorder':2}
[g.ax_heatmap.axhline(ii,**kwargs_spine) for ii in [0,g.data.shape[0]]]
[g.ax_heatmap.axvline(ii,**kwargs_spine) for ii in [0,g.data.shape[1]]]
[g.ax_cbar.axhline(ii,**kwargs_spine) for ii in list(g.ax_cbar.get_ylim())]
[g.ax_cbar.axvline(ii,**kwargs_spine) for ii in list(g.ax_cbar.get_xlim())]

# SAVE FIGURE
plt.savefig(f"{dir_figure}/main/figure-3-strain-enrichment-heatmap.png",dpi=600,bbox_inches='tight')
plt.close()


