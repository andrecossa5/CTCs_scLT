"""
Control clones.
"""

import os
import numpy as np
import pandas as pd
import plotting_utils as plu
import matplotlib.pyplot as plt
plu.set_rcParams()


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')


# Read old and new meta
meta = pd.read_csv(os.path.join(path_data, 'cells_meta_new.csv'), index_col=0)

meta.groupby('timepoint')['GBC'].nunique()

df_ = (
    meta.groupby('sample')
    ['GBC'].nunique()
    .reset_index()
    .assign(timepoint=lambda x: x['sample'].str.split('_').str[0])
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
    .assign(timepoint=lambda x: x['sample'].str.split('_').str[0])
)

fig, ax = plt.subplots(figsize=(3,3))
order = ['PT', 'lung', 'CTC']
plu.strip(df_, 'timepoint', 'SH', ax=ax, color='k', x_order=order)#, add_stats=True, pairs=[('PT', 'lung'), ('PT', 'CTC'), ('lung', 'CTC')])
fig.tight_layout()
plt.show()

df_ = (
    meta
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='f')
    .assign(timepoint=lambda x: x['sample'].str.split('_').str[0])
)
fig, ax = plt.subplots(figsize=(3,3))
