"""
Fig2. Clone trajectories.
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib.pyplot as plt
plu.set_rcParams()


##


# Source utils: insert your <code/main_analysis> path here
sys.path.insert(0, '/Users/cossa/Desktop/projects/CTCs_scLT/code/main_analysis')
from utils import *

# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata


##


# Plot changing modules vs pseudotime for representative clones
GBCs = adata.obs.loc[adata.obs['metastatic_status']!='other', 'GBC'].astype(str).unique()
up = ['H1']
down = ['H20', 'H14']

i = 7
gbc_ = GBCs[i]
cells_ = adata.obs.loc[lambda x: x['GBC']==GBCs[i]].index
df_ = adata.obs.loc[cells_, ['timepoint', 'dpt_pseudotime']+up+down]
df_[up+down] = df_[up+down].apply(lambda x: (x-x.mean()) / x.std(), axis=0)
df_['timepoint'] = pd.Categorical(df_['timepoint'], categories=['PT', 'lung', 'CTC'], ordered=True)

# Statistical test
cols = up + down
gam_results = {col: fit_gam_association(df_, x_col='dpt_pseudotime', y_col=col, n_splines=5) for col in cols}

fig, axs = plt.subplots(len(cols), 1, figsize=(2, 3.5), sharex=True)
for ax, col in zip(axs, cols):
    r = gam_results[col]
    ax.scatter(df_['dpt_pseudotime'], df_[col], c='lightgray', s=2, alpha=0.4, zorder=1)
    ax.plot(r['x_test'], r['y_pred'], color='steelblue', lw=2, zorder=3)
    ax.fill_between(r['x_test'], r['ci'][:, 0], r['ci'][:, 1], color='steelblue', alpha=0.2, zorder=2)
    ax.set_ylabel(col)
    ax.set_title(f'R²={r["pseudo_r2"]:.2f}  p={r["pval"]:.1e}')
axs[2].set_xlabel('DPT')
fig.tight_layout()
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', f'{gbc_}_Fig2i_clone_trajectories.pdf'))


##
