"""
Visualize preprocessed dataset.
"""

import os
import numpy as np
import pandas as pd
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

# Colors
timepoint_colors = plu.create_palette(adata.obs, 'timepoint', plu.darjeeling)
mouse_colors = plu.create_palette(adata.obs, 'mouse', sc.pl.palettes.vega_10_scanpy)
gbc_colors = plu.create_palette(adata.obs, 'GBC', sc.pl.palettes.zeileis_28)
adata.uns['timepoint_colors'] = timepoint_colors.values()
adata.uns['mouse_colors'] = mouse_colors.values()
adata.uns['GBC_colors'] = gbc_colors.values()

# Viz
fig, ax = plt.subplots(1,5,figsize=(10,1.5))
sc.pl.umap(adata, color='n_UMIs', ax=ax[0],
           show=False, frameon=False, legend_loc=None, size=1, title='')
sc.pl.umap(adata, color='n_genes', ax=ax[1],
           show=False, frameon=False, legend_loc=None, size=1, title='')
sc.pl.umap(adata, color='timepoint', ax=ax[2],
           show=False, frameon=False, legend_loc=None, size=1, title='')
plu.add_legend(ax=ax[2], colors=timepoint_colors,
               label='Timepoint', loc='upper left', bbox_to_anchor=(0,1))
sc.pl.umap(adata, color='mouse', ax=ax[3],
           show=False, frameon=False, legend_loc=None, size=1, title='')
plu.add_legend(ax=ax[3], colors=mouse_colors,
               label='Mouse', loc='upper left', bbox_to_anchor=(0,1))
sc.pl.umap(adata, color='GBC', ax=ax[4],
           show=False, frameon=False, legend_loc=None, size=1, title='')

fig.tight_layout()
plt.show()

# Viz
fig, ax = plt.subplots(figsize=(3,3))
sc.pl.umap(adata, color='timepoint', ax=ax,
           show=False, frameon=False, legend_loc=None, size=2.5, title='')
plu.add_legend(ax=ax, colors=timepoint_colors,
               label='Timepoint', loc='upper left', bbox_to_anchor=(0,1))
fig.tight_layout()
plt.show()


##
