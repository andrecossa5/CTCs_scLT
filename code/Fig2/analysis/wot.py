"""
Trajectory inference: CoSpar
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import wot as wot


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')
path_wot = os.path.join(path_data, 'wot')

# Re-create adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))

# Fix time info
adata.obs['time_info'] = np.select(
    [adata.obs['timepoint'] == 'PT', 
     adata.obs['timepoint'] == 'lung', 
     adata.obs['timepoint'] == 'CTC'],
    [0, 1, 2]
)

# create OTModel
ot_model = wot.ot.OTModel(adata[:,adata.var['highly_variable']].copy(), day_field='time_info', epsilon = 0.05, lambda1 = 1,lambda2 = 50) 

# OT maps
ot_model.compute_all_transport_maps(tmap_out=path_wot)


##