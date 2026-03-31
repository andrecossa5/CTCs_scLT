"""
Dataset preprocessing.
"""

import os
import anndata
import pickle
import colorsys 
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import save_npz
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Dimred I: batch correction
X_scVI = pd.read_csv(os.path.join(path_data, 'X_scVI.csv'), index_col=0)
adata.obsm['X_scVI'] = X_scVI.loc[adata.obs_names].values

# Dimred II: initial visualization and neighborhood graph
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_scVI', metric='cosine')
sc.tl.diffmap(adata, n_comps=10)
sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_diffmap', metric='cosine')
# sc.tl.umap(adata)

# fig, axs =  plt.subplots(1,4, figsize=(8.5,2))
# sc.pl.umap(adata, color='n_UMIs', ax=axs[0], show=False)
# sc.pl.umap(adata, color='n_genes', ax=axs[1], show=False)
# sc.pl.umap(adata, color='mouse', ax=axs[2], show=False)
# sc.pl.umap(adata, color='timepoint', ax=axs[3], show=False)
# fig.tight_layout()
# plt.show()

# Clustering
r = 0.1 # 0.1-0.3, chosen 0.1 for cell state annotation, 0.3 for PAGA and visualization
sc.tl.leiden(adata, resolution=r, key_added=f'leiden_{r}')

# Pseudotime 

# fig, axs =  plt.subplots(1,4,figsize=(8.5,2))
# sc.pl.diffmap(adata, components=['2,4'], color='n_UMIs', ax=axs[0], show=False)
# sc.pl.diffmap(adata, components=['2,4'], color='n_genes', ax=axs[1], show=False)
# sc.pl.diffmap(adata, components=['2,4'], color='percent_mt', ax=axs[2], show=False)
# sc.pl.diffmap(adata, components=['2,4'], color='timepoint', ax=axs[3], show=False)
# fig.tight_layout()
# plt.show()

adata.uns["iroot"] = adata.obsm['X_diffmap'][:,1].argmax()
sc.tl.dpt(adata, n_dcs=10)

# Dimred III: fine-tuned visualization
sc.tl.paga(adata, groups='leiden_0.3')
sc.pl.paga(adata, color='leiden_0.3', show=False)
sc.tl.umap(adata, init_pos='paga')
sc.tl.draw_graph(adata, layout='fa', init_pos='paga', root=adata.uns["iroot"], n_jobs=-1, random_state=1234)

# Set colors for visualization

# Colors
timepoint_colors = plu.create_palette(adata.obs, 'timepoint', plu.darjeeling)
leiden_colors = plu.create_palette(adata.obs, 'leiden_0.1', sc.pl.palettes.vega_10_scanpy)                                                                                                                                                                                                                     
n_clones = adata.obs['GBC'].nunique()                                                                                                                           
colors = [                                                                                                                                                      
    '#%02x%02x%02x' % tuple(int(c * 255) for c in colorsys.hsv_to_rgb(h, 0.7, 0.9))                                                                             
    for h in np.linspace(0, 0.85, n_clones)                                                                                                                     
]
colors = list(np.random.default_rng(seed=10).permutation(colors))
colors = list(np.random.permutation(colors))
gbc_colors = plu.create_palette(adata.obs, 'GBC', colors, saturation=0.5)
adata.uns['timepoint_colors'] = list(timepoint_colors.values())
adata.uns['leiden_0.1_colors'] = list(leiden_colors.values())
adata.uns['GBC_colors'] = list(gbc_colors.values())
with open(os.path.join(path_data, 'colors.pkl'), 'wb') as f:
    pickle.dump(
        {'timepoint': timepoint_colors, 'leiden_0.1': leiden_colors, 'GBC': gbc_colors}, 
        f
    )

# fig, axs = plt.subplots(1,3,figsize=(7,2.5))
# sc.pl.umap(adata, color='timepoint', ax=axs[0],
#            show=False, frameon=False, legend_loc=None, size=2.5, title='')
# sc.pl.umap(adata, color='leiden_0.1', ax=axs[1],
#            show=False, frameon=False, legend_loc=None, size=2.5, title='')
# sc.pl.umap(adata, color='GBC', ax=axs[2],
#            show=False, frameon=False, legend_loc=None, size=2.5, title='')
# fig.tight_layout()
# plt.show()

# Clean up anndata here
# adata.obs.drop(columns=['leiden_0.3'], inplace=True)
# del adata.uns['mouse_colors']
# ...

# Write anndata and its components
# anndata.settings.allow_write_nullable_strings = True
# adata.write(os.path.join(path_data, 'adata.h5ad'))

# Cell meta
adata.obs.to_csv(os.path.join(path_data, 'cells_meta.csv'))

# Connectivities
save_npz(os.path.join(path_data, 'connectivities.npz'), adata.obsp['connectivities'])


##

