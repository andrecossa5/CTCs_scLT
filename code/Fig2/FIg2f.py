"""
Fig2 f. Hotspot modules.
"""

import os
import colorsys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors as mcolors
from scipy.cluster.hierarchy import linkage, leaves_list
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Refactor column
adata.var['hotspot'] = adata.var['hotspot'].map(lambda x: f'H{x}')

# Read colors
with open(os.path.join(path_data, 'colors.pkl'), 'rb') as f:
    COLORS = pickle.load(f)

# Read hotspot local correlations
hs = pd.read_csv(
    os.path.join(path_data, 'hotspot_local_correlation_z.csv'),
    index_col=0
)

# Read hotspot modules annotations
modules = pd.read_csv(os.path.join(path_data, 'modules_labels.csv'), index_col=0)
modules.loc[modules['Label'].str.startswith('TN'), 'Label'] = 'TNF-alpha'
modules.loc[modules['Label']=='EMT ', 'Label'] = 'EMT'
modules_annotation = modules['Label'].to_dict()
order = leaves_list(linkage(hs.values, method='average', metric='cosine'))

# --- Module colors ---

def _moderate(color, vmin=0.5, vmax=0.85):
    r, g, b = mcolors.to_rgb(color)
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    v = np.clip(v, vmin, vmax)
    return mcolors.to_hex(colorsys.hsv_to_rgb(h, s, v))

genes_ordered = hs.index[order].tolist()
gene_modules = adata.var.loc[genes_ordered, 'hotspot'].map(modules_annotation).values
unique_modules = sorted(gene_modules)
n_mods = len(unique_modules)
palette = [_moderate(c) for c in sns.color_palette("husl", n_mods)]
module_colors = {m: c for m, c in zip(unique_modules, palette)}
row_colors = [module_colors[m] for m in gene_modules]


##


fig, ax = plt.subplots(figsize=(5, 3))

# Heatmap
ax.imshow(hs.iloc[order, order], cmap='Spectral_r', vmin=-25, vmax=25)
plu.format_ax(ax=ax, rotx=90, xticks=[], yticks=[], xticks_size=.2)

# Row annotation strip (left)
axins_row = ax.inset_axes((1.005, 0, 0.05, 1))
cb_row = plt.colorbar(
    matplotlib.cm.ScalarMappable(
    cmap=matplotlib.colors.ListedColormap(row_colors[::-1])),
    cax=axins_row, orientation='vertical'
)
cb_row.ax.set(xticks=[], yticks=[])
cb_row.outline.set_linewidth(0.1)

# Column annotation strip (top)
axins_col = ax.inset_axes((0, 1.005, 1, 0.05))
cb_col = plt.colorbar(
    matplotlib.cm.ScalarMappable(
    cmap=matplotlib.colors.ListedColormap(row_colors)),
    cax=axins_col, orientation='horizontal'
)
cb_col.ax.set(xticks=[], yticks=[])
cb_col.outline.set_linewidth(0.1)

# Labels
ax.set(xlabel='Genes', ylabel='Genes')

# Cbar modules
plu.add_legend(module_colors, label='Modules', ax=ax, 
               artists_size=7, ticks_size=6, ncols=1, bbox_to_anchor=(1.1, 1.1))

# Label
plu.add_cbar(
    hs.values.flatten(),
    ax=ax, palette='Spectral_r', vmin=-25, vmax=25, label='Local correlation', 
    layout=( (1.25,-.01,.3,.03), 'top', 'horizontal' )
)

fig.subplots_adjust(left=0.05, right=0.7, top=0.9, bottom=0.1)
plt.show()
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2f.pdf'))
