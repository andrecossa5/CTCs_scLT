"""
Control clones.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib.pyplot as plt
plu.set_rcParams()


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')


# Read old and new meta
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
meta = adata.obs.copy()

# Number of clones per timepoint
meta['GBC'] = meta['GBC'].astype('str')
meta['sample'] = meta['sample'].astype('str')
meta.groupby('timepoint')['GBC'].nunique()

df_ = (
    meta.groupby('sample')
    ['GBC'].nunique()
    .reset_index()
    .assign(timepoint=lambda x: x['sample'].map(lambda s: s.split('_')[0]))
)

fig, ax = plt.subplots(figsize=(3,3))
order = ['PT', 'lung', 'CTC']
plu.strip(df_, 'timepoint', 'GBC', ax=ax, color='k', x_order=order)#, add_stats=True, pairs=[('PT', 'lung'), ('PT', 'CTC'), ('lung', 'CTC')])
ax.set_yscale('log')
fig.tight_layout()
plt.show()

df_ = (
    meta
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='f')
    .groupby('sample')
    .apply(lambda x: np.sum(-np.log10(x['f'])*x['f']))
    .reset_index(name='SH')
    .assign(timepoint=lambda x: x['sample'].map(lambda s: s.split('_')[0]))
)

fig, ax = plt.subplots(figsize=(3,3))
order = ['PT', 'lung', 'CTC']
plu.strip(df_, 'timepoint', 'SH', ax=ax, color='k', x_order=order)#, add_stats=True, pairs=[('PT', 'lung'), ('PT', 'CTC'), ('lung', 'CTC')])
fig.tight_layout()
plt.show()



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


# Log clones
L = []
for dataset in df_['dataset'].unique():
    df_sample = df_.query('dataset == @dataset').pivot_table(index='GBC', columns='timepoint', values='n', fill_value=0)
    L.append(
        df_sample
        .loc[(np.sum(df_sample>0, axis=1)==1) & (np.sum(df_sample, axis=1)>=50)]
        [['PT', 'lung', 'CTC']]
        .sort_values('CTC', ascending=False)
        .assign(dataset=dataset)
    )

df_long = pd.concat(L, axis=0)

# 105 longitudinal clones (11 + 16 + 28 + 50 across datasets)



