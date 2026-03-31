"""
Gene modules with Hotspot.
"""

import os
import pickle
import scanpy as sc
import pandas as pd
import numpy as np
import hotspot as hsc


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')


if __name__ == '__main__':

    # Read adata
    adata = sc.read(os.path.join(path_data, 'adata.h5ad'))

    # Hotspot object
    hs = hsc.Hotspot(
        adata[:,adata.var['highly_variable']].copy(),
        layer_key="raw",
        model='danb',
        distances_obsp_key='distances',
        umi_counts_obs_key="n_UMIs"
    )

    # Run
    hs.create_knn_graph(weighted_graph=False, n_neighbors=30)
    hs_results = hs.compute_autocorrelations(jobs=16)
    hs_genes = hs_results.loc[hs_results.FDR < 0.05].index  # Select genes
    local_correlations = hs.compute_local_correlations(hs_genes, jobs=16)
    modules = hs.create_modules(
        min_gene_threshold=30, core_only=True, fdr_threshold=0.05
    )

    # Save
    modules.to_csv(os.path.join(path_data, 'hotspot_modules.csv'))
    with open(os.path.join(path_data, 'hotspot.pkl'), 'wb') as f:
        pickle.dump(hs, f)


##



