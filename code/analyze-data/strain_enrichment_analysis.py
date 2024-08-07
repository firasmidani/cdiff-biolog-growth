#!/usr/bin/env python

# import off-the-shelf packages
import gseapy as gp
import numpy as np
import pandas as pd

# READ DATA
df_summ = pd.read_csv(
    "../../amiga-biolog/summary/merged_summary_norm_sub_by_median.txt",
    sep='\t',header=0,index_col=0,low_memory=False
    )

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
df_meta = df_meta.fillna('Unknown')
df_meta = df_meta.set_index('Strain_ID')

# IDENTIFY SUBSTRATES OF INTEREST

# substrates that support growth
thresh_od_growth = 0.1 # OD must increase by at least 0.1
thresh_od_death = -0.05 # OD must decrease by at least 0.05
num_isolates = df_wide.shape[0] # number of isolates
thresh_hit = np.floor(0.1*num_isolates) # OD increase must occur for at least 10% of isolates

matches = lambda x: [1 if x >= thresh_od_growth else 0][0]
matches = df_wide.map(matches).sum(0).sort_values(ascending=False)
substrates = list(matches[matches >= thresh_hit].index)

# DEFINE STRAIN SETS FOR ENRICHMENT ANALYSIS 

strain_sets = {}

# sets that are simply isolates belonging to same ribotype 
ls_sets_ribotype = ['RT017','RT023','RT255','RT078']
for ribotype in ls_sets_ribotype: 
    strain_sets[ribotype] = list(df_meta[df_meta.Ribotype==ribotype].index.values)

# sets where strains certainly or possibly belong to one of these specific robptypes
# for example, isolates that are typed as RT027;RT036 are included in the RT027+ set
for ribotype in ['RT027','RT106','RT014']:
    set_name = f"{ribotype}+"
    ls_strains = [ii for ii in df_meta.Ribotype.values if ii.startswith(ribotype)]
    strain_sets[f"{ribotype}+"] = df_meta[df_meta.Ribotype.isin(ls_strains)].index.values

# set for Clade 5 isolates excluding RT078 isolates
ls_strains = [ii for ii in df_meta.Ribotype.values if ii in ['RT033','RT126','RT288']]
strain_sets['Clade5+'] = df_meta[df_meta.Ribotype.isin(ls_strains)].index.values

# set for lab-adapted strains
strain_sets['LabAdapted'] = ['R20291','M68','VPI 10463','CD630-Britton','CD630-ermS','CD630-Savidge']

# RUN ENRICHMENT ANALYSIS
isolate_ranks = df_wide.loc[:,substrates].join(df_meta)

ls_results = []
for substrate in substrates: 
    
    ranks = isolate_ranks.loc[:,[substrate]].dropna()
    result = gp.prerank(
        rnk = ranks,
        gene_sets = strain_sets, 
        min_size = 4, 
        permutation_num = 100000, 
        weight = 0.5
    )
    result = result.res2d
    result.loc[:,'Substrate'] = substrate
    ls_results.append(result)

df_results = pd.concat(ls_results)

# save results
df_results.to_csv('../../tables/strain_enrichment_analysis.tsv',sep="\t")

