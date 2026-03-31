"""
Cell state annotation.
"""

import os
import anndata
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


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata_orig = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata = adata_orig.copy()

# Read colors
with open(os.path.join(path_data, 'colors.pkl'), 'rb') as f:
    COLORS = pickle.load(f)


##


# Annotate cell states

# Add for compatibility with scanpy's rank_genes_groups
adata.uns['log1p']['base'] = np.exp(1) 

# Filter too lowly expressed genes for markers genes computation
chosen = 'leiden_0.1'
min_n_cells = adata.obs[chosen].value_counts().min() * .5
test = (~adata.var_names.str.contains('.', regex=False)) & \
       (np.sum(adata.layers['raw'].toarray()>0, axis=0) >= min_n_cells)
adata = adata[:,list(test)].copy()

# Marker genes
sc.tl.rank_genes_groups(adata, groupby=chosen, method='wilcoxon', pts=True)

# Filter and visualize marker genes
res_de = format_rank_genes_groups(adata, key='cell_state_markers')
order = order_groups(adata, groupby=chosen, obsm_key='X_diffmap', n_dims=10)
res_de_filtered = format_rank_genes_groups(adata, key='cell_state_markers', filter_results=True, min_pct_group=.5, max_pct_rest=.5, min_log2FC=1.5, max_pval_adj=.05)
top_markers = get_top_markers(res_de_filtered, order_groups=order, ntop=3)
df_plot = res_de.query('gene in @top_markers')

# fig, ax = plt.subplots(figsize=(6,3.5))
# plu.dotplot(df_plot, 'gene', 'group', 
#             order_x=top_markers, order_y=order,
#             color='log2FC', size='pct_group', ax=ax, vmin=-5, vmax=5)
# plu.format_ax(ax=ax, rotx=90, xlabel='', ylabel='Clusters')
# ax.get_legend().set_bbox_to_anchor((1,1.2))
# fig.tight_layout()
# plt.show()

# fig, ax = plt.subplots(figsize=(5,3))
# plu.bb_plot(adata.obs, cov1='timepoint', cov2=chosen, categorical_cmap=COLORS[chosen], ax=ax)
# fig.tight_layout()
# plt.show()

# Annotate individual clusters: e.g. cluster 4 --> IFN response
mapping = {
    '0' : 'NFkB',
    '1' : 'OXPHOS',
    '2' : 'Myc/E2F',
    '3' : 'Hypoxia/Glycolisis',
    '4' : 'IFN'
}
# gsea_path = os.path.join(path_data, 'CELL_STATE_MARKERS.xlsx')
# with pd.ExcelWriter(gsea_path) as writer:
#     for cluster in adata.obs['leiden_0.1'].unique():
#         res_gsea = run_GSEA(
#             res_de.query('group==@cluster').set_index('gene')['log2FC'], 
#             max_size_set=1000, 
#             collections='MSigDB_Hallmark_2020'
#         )[1]
#         res_gsea['Term'] = res_gsea['Term'].map(lambda x: x.replace('MSigDB_Hallmark_2020__', ''))
#         res_gsea
#         res_gsea.to_excel(writer, sheet_name=mapping[cluster], index=False)

# Add slots to original adata
adata_orig.obs['cell_state'] = adata_orig.obs[chosen].map(mapping)
adata_orig.obs['cell_state'] = adata_orig.obs['cell_state'].astype('category')
adata_orig.uns['cell_state_markers'] = adata.uns['rank_genes_groups'].copy()

# Update colors
colors = plu.create_palette(adata_orig.obs, 'cell_state', sc.pl.palettes.vega_10_scanpy)                                                                                                                                                                                                                     
adata_orig.uns['cell_state_colors'] = list(colors.values())
with open(os.path.join(path_data, 'colors.pkl'), 'wb') as f:
    pickle.dump({**COLORS, 'cell_state': colors}, f)


##


# Hotspot modules
modules = pd.read_csv(os.path.join(path_data, 'hotspot_modules.csv'), index_col=0)
adata_orig.var['hotspot'] = adata_orig.var_names.map(modules['Module'].to_dict()).fillna(-1).astype(int)

# Score modules
# ora_path = os.path.join(path_data, 'GENE_MODULES.xlsx')
# with pd.ExcelWriter(ora_path) as writer:
#     for module in modules['Module'].unique():
#         if module == -1:
#             continue
#         gene_list = list(adata_orig.var.query('hotspot==@module').index)
#         # sc.tl.score_genes(adata_orig, gene_list=gene_list, score_name=f'H{module}', n_bins=40)
#         df_,_ = run_ORA(gene_list)
#         df_.to_excel(writer, sheet_name=f'H{module}', index=False)


# Confirm cell state annotations
# cols = adata_orig.obs.columns[adata_orig.obs.columns.str.startswith('H')]
# adata_orig.obs.groupby('timepoint')[cols].median().T.plot.bar(figsize=(5,3))
# plt.show()

# Individual gene programs: e.g. H2, a marker of the IFN response module
# g = adata_orig.var.query('hotspot == 2').index.to_list()
# run_ORA(g)

# Write final anndata and its components
# anndata.settings.allow_write_nullable_strings = True
# adata_orig.write(os.path.join(path_data, 'adata.h5ad'))


##