"""
Dataset preprocessing.
"""

import os
import pandas as pd
import scanpy as sc
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

# Write anndata with expanded HVGs selection and cell cycle annotation
# import anndata
# anndata.settings.allow_write_nullable_strings = True
# adata.write(os.path.join(path_data, 'adata.h5ad'))


##
