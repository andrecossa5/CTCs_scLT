"""
Fig2 a-d. Preprocessed dataset visualzation
"""

import os
import pickle
import numpy as np
import scanpy as sc
import plotting_utils as plu
import matplotlib.pyplot as plt
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Read colors
with open(os.path.join(path_data, 'colors.pkl'), 'rb') as f:
    COLORS = pickle.load(f)


##


fig, axs = plt.subplots(1,3,figsize=(6,2), width_ratios=[1,1,1])

sc.pl.umap(adata, color='timepoint', ax=axs[0],
           show=False, frameon=False, legend_loc='on data', size=2.5, title='')
sc.pl.umap(adata, color='cell_state', ax=axs[1],
           show=False, frameon=False, legend_loc='on data', size=2.5, title='')
sc.pl.umap(adata, color='dpt_pseudotime', ax=axs[2],
           show=False, frameon=False, size=2.5, cmap='magma', title='', colorbar_loc=None)
plu.add_cbar([0,1], palette='magma', label='DPT', ax=axs[2], 
             layout=( (0.1,.05,.25,.03), 'top', 'horizontal' ),
             label_size=6, ticks_size=4
             )


fig.subplots_adjust(top=0.8, bottom=0.1, left=0.1, right=0.9, wspace=0.2)
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2ac.pdf'))


##


fig, ax = plt.subplots(figsize=(2,1.2))
plu.bb_plot(adata.obs, cov1='timepoint', cov2='cell_state', ax=ax, categorical_cmap=COLORS['cell_state'])
fig.tight_layout()
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2d.pdf'))


##


fig, ax = plt.subplots(figsize=(2,1.1))
plu.dist(adata.obs, x='dpt_pseudotime', by='timepoint', ax=ax, categorical_cmap=COLORS['timepoint'])
plu.format_ax(ax, xlabel='DPT', ylabel='Density', title='', reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2e.pdf'))


##