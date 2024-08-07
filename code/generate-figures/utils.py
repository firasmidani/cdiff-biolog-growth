#!/usr/bin/env python

# TABLE OF CONTENTS

# LIST OF CLASSES
# customAxes (matplotlib.axes.Axes)
#   add_grid
#   enlarge_tick_labels
#   get_center_xlim
#   thicken_spines

# LIST OF FUNCTIONS
# get_centeroid
# comp_covariance
# init_figure
# read_csv
# read_table_as_dict
# subsetDf

# import off-the-shelf packages
import operator
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

from matplotlib import rcParams, projections
from matplotlib.axes import Axes
from scipy.stats import variation

# define matplotlib-related parameters
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['savefig.dpi'] = 300
rcParams['savefig.format'] = 'png'
rcParams['savefig.bbox'] = 'tight'
rcParams['savefig.pad_inches'] = 0.4 # ~1 cm
rcParams['savefig.transparent'] = False

fontsize = 15

# define input/output figures
dir_tables = "../../tables"
dir_figure = "../../figures"
dir_config = "../../configs"

# def custom matplotlib axes class
class customAxes(Axes):
    '''
    Custom sub-class of matplotlib.axes.Axes which makes it easer to streamline all plots.
    
    Args:
        Axes (matplotlib.axes.Axes) which encapsulates all elements of a matplotlib sub-plot
    '''
    
    name = 'biolog'
    
    def add_grid(self,lw=0.25,zorder=1): 
        '''
        Adds boilerplate grid lines.
        '''
        
        self.grid(lw=lw,color='black',alpha=0.25,zorder=zorder)
    
    def enlarge_tick_labels(self,fontsize=20):
        '''
        Set fontsize for x-axis and y-axis tick labels. 
        
        Args:
            fontsize (numeric) default is 20
        '''
        
        [ii.set(fontsize=fontsize) for ii in self.get_xticklabels()+self.get_yticklabels()]
        
    def get_center_xlim(self):
        
        x_bounds = self.get_xlim()
        x_center = 0.5 * (x_bounds[1] - x_bounds[0]) + x_bounds[0]
        
        return x_center
    
    def set_centered_xlabel(self,y=-1,s='x-axis placeholder',fontsize=12):

        self.text(
            self.get_center_xlim(),
            y = y,
            s = s,
            ha = 'center',
            va = 'center',
            fontsize = fontsize,
            transform = self.transData
        )
    
    def thicken_spines(self,lw=2,which ='all'):
        '''
        Adjust thickness of spines for a sub-plot.
        
        Args:
            lw (float) line width
            which (str) acceptable option are 'all', 'l', 'r', 't', 'b'
        '''
        
        dict_which = {'l':'left','r':'right','t':'top','b':'bottom'}
        
        if which == 'all':
            which = dict_which.values()
        elif isinstance(which,list) or isinstance(which,str):
            which = [dict_which[ii] for ii in which]
        else:
            raise AttributeError(f'{which} is not a valid argument for "which"')
        
        [self.spines[ii].set(lw=lw) for ii in which]

def cluster_rows(df):
    '''
    Given a table, cluster rows based on values in columns then return 
    order of the rows after clustering. 
    
    Args:
        df (pandas.DataFrame)

    Returns
        list
    '''

    linkage = sch.linkage(df,method='average',optimal_ordering=False)
    dendrogram = sch.dendrogram(linkage,no_plot=True)
    leaves = dendrogram['leaves']
    leaves = df.index.values[leaves]
    
    return leaves 

def comp_centeroid(arr):
    '''
    Compute centroid for multi-dimensional array of numbers. 

    Args:
        arr (numpy.Array) 2-dimensional array
    
    Params:
        2-tuple (float, float) centroids for the first and second dimensions
    '''
    
    if arr.ndim != 2:
    
        print(f"UserError: array must be 2-dimensional")
        return None
    
    elif arr.shape[1] != 2:
    
        print(f"UserError: array must be 2-dimensional")
        return None 

    else:
    
        length = arr.shape[0]
        sum_x = float(np.sum(arr[:, 0]))
        sum_y = float(np.sum(arr[:, 1]))
        x_cent = sum_x/length
        y_cent = sum_y/length

    return (x_cent, y_cent)

def comp_covariance(series):

    values = [ii for ii in series.values if not np.isnan(ii)]
    mu = np.mean(values)

    # if mean is negative, let's invert all values 
    if mu < 0: weight = -1
    else: weight = 1

    weighted_values = [weight*ii for ii in values if not np.isnan(ii)]

    return variation(weighted_values)

def get_clustermap_labels(g):
    '''
    Get the row and column labels after the index and column headers of the clustermap
    input dataframe were re-ordered by the clustering. 

    Args:
        g (seaborn.matrix.ClusterGrid)
    
    Returns:
        row_labels (list)
        col_labels (list)
    '''

    df = g.data

    if g.dendrogram_col:
        col_labels = [df.columns[ii] for ii in g.dendrogram_col.reordered_ind]
    else:
        col_labels = df.columns
    
    if g.dendrogram_row:
        row_labels = [df.index.values[ii] for ii in g.dendrogram_row.reordered_ind]
    else:
        row_labels = df.index.values
    
    return row_labels, col_labels

def get_sig_annot(num):
    '''
    Returns symbols indicating level of significance for input number.

    Args:
        num (float)
    
    Returns:
        (str or None)
    '''

    if num <= 1e-3: return '***'
    elif num <= 1e-2: return '**'
    elif num <= 0.05: return '*'

    return None

def init_figure(figsize=(4,4),kwargs={}):
    '''
    Initializes a subplot with matplotlib.figure.Figure and custom matplotlib.axes.Axes classes
    
    Args:
        figsize (float, float) width, height in inches.
    
    Returns
        fig (matplotlib.figure.Figure)
        ax (matplotlib.axis.Axes)
    '''
    
    projections.register_projection(customAxes)
    
    fig,ax = plt.subplots(figsize=figsize,subplot_kw={'projection':'biolog'},**kwargs)
    
    return fig,ax 

def read_csv(filepath,kwargs={}): 
    '''
    Read tab-delimited filer with header and index column using low memory settings
    
    Args:
        filepath (str) 
        kwargs (dictionary) keyword arguments that can be passed to pandas.read_csv
    
    Returns:
        (pandas.DataFrame)
    '''
    kwargs_read_csv = {'sep':'\t','header':0,'index_col':0,'low_memory':False}
    
    return pd.read_csv(filepath,**kwargs_read_csv,**kwargs)

def read_table_as_dict(filepath,column=None):
    '''
    Reads a column in a tab-delimited text file and creates a dictionary that maps
    the index in each row with the value in the column. 
    
    Args:
        filepath (str)
        
    Returns
        table_dict (dictionary)
    '''
    
    df = read_csv(filepath)

    if isinstance(column,str) and df.shape[1] > 0:

        table_dict = read_csv(filepath).to_dict([column])
    
    elif (column is None) and df.shape[1] == 1:
    
        # reset column headers and simply grab the only one which corresponds to `0`
        table_dict = df.T.reset_index().T.to_dict()[0]
    
    else:
        
        table_dict = None
    
    return table_dict

def rgb(ls): 
    '''
    Convert RGB color code from 0-255 scale to 0-1 scale.
    
    Args:
        ls (list)
    
    Returns:
        (3-tuple)
    '''

    return tuple([float(ii)/255 for ii in ls])

def savefig(fig_path):
    '''
    Saves high-resolution PNG.
    '''

    plt.savefig(fig_path,format="png",dpi=300,bbox_inches='tight')

def subsetDf(df,criteria):
    '''
    Retains only rows in a pandas.DataFrame that match select criteria. 
    
    Args:
        df (pandas.DataFrame)
        criteria (dictionary): keys (str) are column headers in df, and values (list) are respective column values
        
    Returns (pandas.DataFrame)
    
    Pasted here from https://github.com/firas/midani/amiga
    '''
    
    if criteria is None: return df
    
    # if a criteria (dictionary) has a value that is not list format, put it into a list
    for key,value in criteria.items():
        if not isinstance(value,list):
            criteria[key] = [value]
            
    return df[df.isin(criteria).sum(1)==len(criteria)]

def get_pid_varb_medians(df,varb,norm_ref='A1'):
    '''
    Returns dictionary where keys are Plate_IDs and values are medians for the user-specific metric.

    Args:
        df (pandas.DataFrame)
        varb (str) must be a summary metric (i.e. column header) in the dataframe
        norm_ref ('str') must be a Well ID in the dataframe or 'median', default is `A1`

    Returns: 
        dict_pid_medians (dictionary)
    '''

    data = df.loc[:,['Well','Plate_ID',varb]]
    
    if norm_ref=='A1':

        # get dictionary where keys are 
        pid_medians = data[data.Well=='A1'].drop(['Well'],axis=1)
        dict_pid_medians = pid_medians.set_index('Plate_ID').to_dict()[varb]

    elif norm_ref=='median':
        pid_medians = data.loc[:,['Plate_ID',varb]]
        dict_pid_medians = pid_medians.groupby(['Plate_ID']).median().to_dict()[varb]

    return dict_pid_medians

def renormalize_summary_metrics(df,varb='k_lin',norm_ref='A1',method='sub'):
    '''
    Normalizes a specific metric (varb) using a specific reference (norm_ref) and normalization (method).

    Args:
        df (pandas.DataFrame)
        varb (str) must be a summary metric (i.e. column header) in the dataframe
        norm_ref (str) must be a Well ID in the dataframe or 'median', default is `A1`
        'method' (str) must be either 'sub' (for substraction) or 'div' (for division)

    Returns:
        df (pandas.DataFrame)
    '''

    if method == 'sub':
        op = operator.sub
    elif method == 'div':
        op = operator.truediv
    
    # most likely to be used varbs (included here for convenience)
    varbs = ['auc_lin','auc_log','k_lin','k_log','death_lin','death_log',
             'gr','dr','td','lagC','lagP','t_k','t_gr','t_dr']
    
    df_renorm = df.copy()

    # get dictionary with Plate_IDs as keys and corresponding medians as values 
    varb_medians = get_pid_varb_medians(df_renorm,varb=varb,norm_ref=norm_ref)
    varb_renorms = df.apply(lambda row: op(row[varb],varb_medians[row['Plate_ID']]), axis=1)
    df_renorm.loc[:,f"norm({varb})"] = varb_renorms

    return df_renorm

def process_amiga_output(work_dir,quality_check=False):

    # READ SUMMARY TABLE
    df_summ = read_csv(f"{work_dir}/summary/merged_summary.txt")
    
    if quality_check: 

        # read quality calls for each unique Pleate ID (pid)
        df_pid_quality = read_csv(f"{work_dir}/notes/df_pid_quality.txt")
        df_pid_quality.set_index('Plate_ID',inplace=True)

        # read warning flags for each unique well ID (wid)
        df_wid_quality = read_csv(f"{work_dir}/notes/df_well_flags.txt")
        df_wid_quality = df_wid_quality.map(lambda x: x.strip() if isinstance(x,str) else x)

        # add well flags
        df_summ = pd.merge(df_summ,df_wid_quality,on=['Plate_ID','Well'],how='outer')
        df_summ['Exclude'] = df_summ['Exclude'].fillna(0)

        # add plate quality calls
        df_summ = pd.merge(df_summ,df_pid_quality,on=['Plate_ID'],how='outer')

        # remove low-quality plates
        df_summ = df_summ[df_summ.Quality.isin(['Good','Okay'])]

        # remove low-quality wells
        df_summ = df_summ[df_summ.Exclude==0]

        # clean-up
        df_summ = df_summ.drop(['Quality','Exclude'],axis=1)


    # NORMALIZE SUMMARY METRICS BY MEDIAN WELLS

    # apply functions
    df_summ_renorm_sub = renormalize_summary_metrics(df_summ,varb="k_lin",norm_ref="median",method="sub")
    df_summ_renorm_div = renormalize_summary_metrics(df_summ,varb="k_lin",norm_ref="median",method="div")

    # save renormalized dataframes
    df_summ_renorm_sub.to_csv(f"{work_dir}/summary/merged_summary_norm_sub_by_median.txt",sep="\t")
    df_summ_renorm_div.to_csv(f"{work_dir}/summary/merged_summary_norm_div_by_median.txt",sep="\t")

