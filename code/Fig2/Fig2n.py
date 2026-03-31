"""
Fig2 m. Longitudinal vs PT-restricted clones.
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import plotting_utils as plu
from adjustText import adjust_text
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Prep df
adata = adata[adata.obs['metastatic_status'] != 'other'].copy()
modules = adata.obs.columns[adata.obs.columns.str.startswith('H')]
adata.obs[modules] = (
    adata.obs[modules]
    .apply(lambda x: (x - x.mean()) / x.std(), axis=0)
)


##


# Summarize per metastatic_status
fig, ax = plt.subplots(figsize=(2.2, 1.6))

df_ = (
    adata.obs[['metastatic_status'] + ['H1', 'H20', 'H14']]
    .melt(id_vars='metastatic_status', var_name='module', value_name='zscore')
)
df_['metastatic_status'] = df_['metastatic_status'].astype('str')
df_.loc[df_['metastatic_status']=='Pro-metastatic', 'metastatic_status'] = 'Longitudinal'
df_.loc[df_['metastatic_status']=='Non-metastatic', 'metastatic_status'] = 'PT-restricted'

colors = {'Longitudinal': '#bd2a09', 'PT-restricted': '#0aad69'}

plu.violin(df_, x='module', y='zscore', by='metastatic_status', 
    categorical_cmap=colors,
    by_order=['Longitudinal', 'PT-restricted'], ax=ax,
    # add_stats=True, pairs=[['Longitudinal', 'PT-restricted']],
)
plu.format_ax(ax, xlabel='', ylabel='z-score', reduced_spines=True)
plu.add_legend(
    ax=ax, colors=colors, label='', loc='center', bbox_to_anchor=(0.5, 1.1),
    ncols=2, ticks_size=6, artists_size=7
)

fig.subplots_adjust(left=.25, top=.8, bottom=.2, right=.8)
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2n.pdf'))


##