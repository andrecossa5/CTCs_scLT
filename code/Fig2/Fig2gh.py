"""
Fig2 g-h. Hotspot modules.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib.pyplot as plt
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Prep df
modules = adata.obs.columns[adata.obs.columns.str.startswith('H')]
adata.obs[modules] = (
    adata.obs[modules]
    .apply(lambda x: (x - x.mean()) / x.std(), axis=0)
)

# Group by timepoint and average
df_ = (
    adata.obs
    .groupby('timepoint')
    [modules].mean()
    .loc[['PT', 'lung', 'CTC']]
    .apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    .T
)
_dom = df_.idxmax(axis=1)
order = [
    idx for tp in ['PT', 'lung', 'CTC']
    for idx in df_[_dom == tp].sort_values(tp, ascending=False).index
]
df_ = df_.loc[order]


##

fig, ax = plt.subplots(figsize=(1, 3.5))
plu.plot_heatmap(df_, ax=ax, palette='mako', cb=False)
plu.format_ax(ax, xlabel='', ylabel='', rotx=90)
plu.add_cbar(
    df_.values.flatten(), ax=ax, palette='mako',
    ticks_size=6,
    label='z-score', label_size=7,
    layout=( (0.1,1.02,.8,.02), 'top', 'horizontal' ),
)
fig.tight_layout()
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2g.pdf'))


##


fig, axs = plt.subplots(3,1,figsize=(1,3))

sc.pl.umap(adata, color='H1', ax=axs[0],
           show=False, frameon=False, cmap='mako', title='', colorbar_loc=None)
sc.pl.umap(adata, color='H12', ax=axs[1],
           show=False, frameon=False, cmap='mako', title='', colorbar_loc=None)
sc.pl.umap(adata, color='H20', ax=axs[2],
           show=False, frameon=False, cmap='mako', title='', colorbar_loc=None)
# plu.add_cbar([0,1], palette='magma', label='DPT', ax=axs[2], 
#              layout=( (0.1,.05,.25,.03), 'top', 'horizontal' ),
#              label_size=6, ticks_size=4
#              )

fig.subplots_adjust(left=0.05, right=0.95, 
                    top=0.95, bottom=0.05, hspace=0.05)
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2h.pdf'))


##