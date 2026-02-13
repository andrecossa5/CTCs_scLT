"""
Define cells states.
"""

import os
import random
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import plotting_utils as plu
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
