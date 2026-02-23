"""
Batch-correction with scVI.
"""

import sys
import os
import scvi
import scanpy as sc


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Subset HVGs
adata = adata[:, adata.var['highly_variable']].copy()
scvi.model.SCVI.setup_anndata(
    adata,
    layer="raw",
    categorical_covariate_keys=["mouse"],
    continuous_covariate_keys=["n_UMIs", "n_genes"],
)

# Train scVI model
model = scvi.model.SCVI(adata)
model.train(max_epochs=250, early_stopping=True)

# Save model
model_dir = os.path.join(path_data, "scvi_model")
model.save(model_dir, overwrite=True)


##
