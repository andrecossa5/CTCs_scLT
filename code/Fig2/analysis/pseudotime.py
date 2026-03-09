"""
Trajectory inference: pseudotime 
"""

import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import plotting_utils as plu
from scipy.sparse import csr_matrix, save_npz
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))


##


# 1. Pseudotime inference with scanpy Diffusion Pseudotime (DPT)
np.random.seed(1234)
test = (adata.obsm['X_umap'][:,0]>8) & \
       (adata.obsm['X_umap'][:,0]<12) & \
       (adata.obsm['X_umap'][:,1]>3.2) & \
       (adata.obsm['X_umap'][:,1]<3.5)
root = np.random.choice(adata.obs_names[test])
adata.uns["iroot"] = adata.obs_names.get_loc(root)
sc.tl.dpt(adata, n_dcs=10)

# Viz
fig, ax = plt.subplots(figsize=(3,3))
sc.pl.umap(adata, color='dpt_pseudotime', ax=ax,
           show=False, frameon=False, legend_loc=None, size=2.5, title='')
fig.tight_layout()
plt.show()

# FA-layout visualization
sc.tl.draw_graph(adata, layout='fa',
                 init_pos='paga', root=root, n_jobs=-1, random_state=1234)
fig, ax = plt.subplots(figsize=(3,3))
sc.pl.draw_graph(adata, color='timepoint', ax=ax,
           show=False, frameon=False, legend_loc=None, size=2.5, title='')
fig.tight_layout()
plt.show()

# Write anndata, and its components for CosPar
# anndata.settings.allow_write_nullable_strings = True
# adata.write(os.path.join(path_data, 'adata.h5ad'))

# Cell meta
adata.obs.to_csv(os.path.join(path_data, 'cells_meta.csv'))
# Connectivities
save_npz(os.path.join(path_data, 'connectivities.npz'), adata.obsp['connectivities'])


##