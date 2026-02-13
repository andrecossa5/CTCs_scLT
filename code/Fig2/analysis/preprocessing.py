"""
Dataset preprocessing.
"""

import os
import random
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import plotting_utils as plu
from sklearn.metrics import pairwise_distances
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')
path_matrices = os.path.join(path_data, 'STARSolo')

# Curated gene list to expand HVGs selection
curated_genes = pd.read_csv(
    os.path.join(path_data, 'curated_genes.txt'), header=None, sep='\t'
)[0]

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Normalization and HVGs selection
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, n_bins=40, batch_key='sample')

# Expand HVGs to curated subsets
adata.var['curated'] = adata.var_names.isin(curated_genes)
adata.var['highly_variable'] = (adata.var['highly_variable']) | (adata.var['curated'])
ki67_idx_hvgs = adata[:,adata.var['highly_variable']].var_names.get_loc('MKI67')
corr = pairwise_distances(
    adata[:,adata.var['highly_variable']].X.toarray().T,
    metric='correlation',
    n_jobs=-1
)
cc_corr = 1-corr[ki67_idx_hvgs,:]
del corr

# Visualize gene-gene correlation distribution (with MKI67 gene)
# sns.kdeplot(cc_corr, fill=True)
# plt.show()

# Remove cell-cycle correlated genes
adata.var['cell_cycle'] = False
adata.var.loc[adata.var['highly_variable'], 'cell_cycle'] = cc_corr > 0.25
adata.var['highly_variable'] = (adata.var['highly_variable']) & (~adata.var['cell_cycle'])
adata.var['highly_variable'].sum()

# Reduce dimensions


# Load scVI batch-corrected latent representation
# ...

# sc.tl.pca(adata, n_comps=100)
# sc.pl.pca_variance_ratio(adata, n_pcs=100)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=25, metric='cosine')
sc.tl.diffmap(adata, n_comps=10)
sc.pp.neighbors(adata, n_neighbors=20, use_rep='X_diffmap', metric='cosine')
random.seed(0)  # Set seed for PAGA reproducibility
sc.tl.paga(adata, groups='timepoint')
sc.pl.paga(adata)
sc.tl.umap(adata, init_pos='paga')

# Vizualize QC metrics in PCA and UMAP space
# fig, axs = plt.subplots(1,4,figsize=(7.5,2))
# sc.pl.pca(adata, color='n_UMIs', ax=axs[0], show=False)
# sc.pl.pca(adata, color='n_genes', ax=axs[1], show=False)
# sc.pl.pca(adata, color='percent_mt', ax=axs[2], show=False)
# sc.pl.pca(adata, color='timepoint', ax=axs[3], show=False)
# fig.tight_layout()
# plt.show()

# fig, axs =  plt.subplots(1,4, figsize=(8.5,2))
# sc.pl.umap(adata, color='n_UMIs', ax=axs[0], show=False)
# sc.pl.umap(adata, color='n_genes', ax=axs[1], show=False)
# sc.pl.umap(adata, color='percent_mt', ax=axs[2], show=False)
# sc.pl.umap(adata, color='timepoint', ax=axs[3], show=False)
# fig.tight_layout()
# plt.show()

# fig, axs =  plt.subplots(1,4,figsize=(8.5,2))
# sc.pl.diffmap(adata, color='n_UMIs', ax=axs[0], show=False)
# sc.pl.diffmap(adata, color='n_genes', ax=axs[1], show=False)
# sc.pl.diffmap(adata, color='percent_mt', ax=axs[2], show=False)
# sc.pl.diffmap(adata, color='timepoint', ax=axs[3], show=False)
# fig.tight_layout()
# plt.show()


##


# Save preprocessed adata
# import anndata
# anndata.settings.allow_write_nullable_strings = True
# adata.write(os.path.join(path_data, 'adata.h5ad'))


##
