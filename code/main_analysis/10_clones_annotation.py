"""
Clones annotation.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')


# Read data
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
meta = adata.obs.copy()


##


# Number of clones per timepoint
df_ = (
    meta
    .groupby('sample')
    ['GBC'].value_counts()
    .reset_index(name='n')
    .assign(
        f=lambda x: x['n'] / x.groupby('sample')['n'].transform('sum')
    )
    .assign(dataset=lambda x: 'M' + x['sample'].map(lambda s: s.split('_')[1]))
    .assign(timepoint=lambda x: x['sample'].map(lambda s: s.split('_')[0]))
)

# Long clones
L = []
for dataset in df_['dataset'].unique():
    df_sample = (
        df_.query('dataset == @dataset')
        .pivot_table(index='GBC', columns='timepoint', values='n', fill_value=0)
    )
    L.append(
        df_sample
        .loc[(np.sum(df_sample>0, axis=1)==3) & (np.sum(df_sample>=10, axis=1)==3)]
        [['PT', 'lung', 'CTC']]
        .sort_values('CTC', ascending=False)
        .assign(dataset=dataset)
    )
df_long = pd.concat(L, axis=0)
df_long = df_long.sort_values('CTC', ascending=False)
df_long.shape

# Non clones
L = []
for dataset in df_['dataset'].unique():
    df_sample = (
        df_.query('dataset == @dataset')
        .pivot_table(index='GBC', columns='timepoint', values='n', fill_value=0)
    )
    L.append(
        df_sample
        .loc[lambda x: (x['CTC']==0) & (x['lung']==0) & (x['PT']>=50)]
        [['PT', 'lung', 'CTC']]
        .sort_values('CTC', ascending=False)
        .assign(dataset=dataset)
    )
df_non = pd.concat(L, axis=0)
df_non = df_non.sort_values('PT', ascending=False)
df_non.shape


##


# Annotate GBCs in adata
adata.obs['metastatic_status'] = np.select(
    [adata.obs['GBC'].isin(df_long.index) & (adata.obs['timepoint'] == 'PT'),
     adata.obs['GBC'].isin(df_non.index) & (adata.obs['timepoint'] == 'PT')], 
    ['Pro-metastatic', 'Non-metastatic'],
    default='other'
)   
adata.obs['metastatic_status'] = adata.obs['metastatic_status'].to_list()
adata.obs['metastatic_status'].value_counts()


##


# Adata
# adata.write(os.path.join(path_data, 'adata.h5ad'))


##

