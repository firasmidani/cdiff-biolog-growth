#!/usr/bin/env python

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

if __name__ == "__main__":

    # import off-the-shelf packages
    import sys
    import operator
    import pandas as pd

    # import custom library
    sys.path.append('../generate-figures/')
    from utils import *

    # READ ARGUMENTS
    work_dir = sys.argv[1]
    qc = str(sys.argv[2])

    # PARSE QUALITY CHECK (QC) ARGUMENT
    if (qc == "True") or (qc == "1"): qc = True
    else: qc = False
    
    # READ SUMMARY TABLE
    df_summ = read_csv(f"{work_dir}/summary/merged_summary.txt")
    
    if qc: 

        # read quality calls for each unique Pleate ID (pid)
        df_pid_quality = read_csv(f"{work_dir}/notes/df_pid_quality.txt")

        # read warning flags for each unique well ID (wid)
        df_wid_quality = read_csv(f"{work_dir}/notes/df_well_flags.txt")
        df_wid_quality = df_wid_quality.applymap(lambda x: x.strip() if isinstance(x,str) else x)

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
        df_summ = df_summ.drop(['Exclude','Quality'],axis=1)

    # NORMALIZE SUMMARY METRICS BY MEDIAN WELLS

    # apply functions
    df_summ_renorm_sub = renormalize_summary_metrics(df_summ,varb="k_lin",norm_ref="median",method="sub")
    df_summ_renorm_div = renormalize_summary_metrics(df_summ,varb="k_lin",norm_ref="median",method="div")

    # save renormalized dataframes
    df_summ_renorm_sub.to_csv(f"{work_dir}/summary/merged_summary_norm_sub_by_median.txt",sep="\t")
    df_summ_renorm_div.to_csv(f"{work_dir}/summary/merged_summary_norm_div_by_median.txt",sep="\t")
