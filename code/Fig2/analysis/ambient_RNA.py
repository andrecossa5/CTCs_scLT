"""
Prepare mtx files for ambient RNA analysis.
"""

import os
import anndata
import pandas as pd
import matplotlib.pyplot as plt


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data', 'STARSolo')

# Write .h5ad files
for sample in os.listdir(path_data):
    
    print(f'Processing sample: {sample}')

    # adata = anndata._io.read_mtx(os.path.join(path_data, sample, 'raw', 'matrix.mtx.gz'))
    # adata = adata.T.copy()
    # adata.var_names = (
    #     pd.read_csv(
    #     os.path.join(path_data, sample, 'raw', 'features.tsv.gz'), 
    #     sep='\t', header=None
    #     )[0]
    # )
    # adata.obs_names = (
    #     pd.read_csv(
    #     os.path.join(path_data, sample, 'raw', 'barcodes.tsv.gz'), 
    #     sep='\t', header=None
    #     )[0]
    # )
    # adata.write(os.path.join(path_data, sample, 'raw.h5ad'))
    # 
    # ##
    # 
    # fig, ax = plt.subplots(figsize=(3,3))
    # order = adata.X.sum(axis=1).A1.argsort()[::-1]
    # sums = adata[order,:].X.sum(axis=1).A1
    # ax.plot(sums)
    # ax.set(xlabel='Droplets, ordered', ylabel='UMI counts per droplet')
    # fig.tight_layout()
    # ax.set_yscale('log')
    # ax.set_xscale('log')
    # fig.savefig(
    #     os.path.join(
    #         path_data, sample, 'Elbow_10x.pdf'
    #     ), dpi=1000
    # )

    path_input = os.path.join(path_data, sample, 'raw.h5ad')
    path_output = os.path.join(path_data, sample, 'cellbender.h5ad')
    os.system(
        f"cellbender remove-background --input {path_input} --output {path_output} --expected-cells 5000 --epochs 100"
    )

##