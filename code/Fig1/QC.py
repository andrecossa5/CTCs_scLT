"""
QC single-cell dataset.
"""

import os
import pandas as pd
import scanpy as sc
import pyranges as pr
import matplotlib.pyplot as plt
import seaborn as sns
import plotting_utils as plu
plu.set_rcParams()


##


# Utils
def filter_genes(adata, gtf):
    """
    Filter non-coding genes and pseudogenes, and add meta information.
    """
    genes = (
        gtf.df
        .query("Feature=='gene' and gene_type=='protein_coding'")
        [["gene_id", "gene_name", "gene_type"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # Annotate and filter genes
    genes['gene_id'] = genes['gene_id'].map(lambda x: x.split('.')[0])
    adata = adata[:, adata.var['gene_ids'].isin(genes['gene_id'])].copy()
    test = (
        (adata.var_names.str.startswith('LC0')) | \
        (adata.var_names.str.startswith('LC1')) | \
        (adata.var_names.str.startswith('AC0')) | \
        (adata.var_names.str.startswith('AC1')) | \
        (adata.var_names.str.startswith('RPS')) | \
        (adata.var_names.str.startswith('RPL'))
    )
    adata = adata[:,~test].copy()

    return adata


##


def filter_cells(adata, thr=.001,
                 min_umis=1000, max_umis=50000,
                 min_genes=1000, max_genes=10000,
                 max_percent_mt=0.15):
    """
    Filter cells based on QC metrics.
    """
    adata_ = adata.copy()
    test = (
        (adata_.obs['n_UMIs'] >= min_umis) &
        (adata_.obs['n_UMIs'] <= max_umis) &
        (adata_.obs['n_genes'] >= min_genes) &
        (adata_.obs['n_genes'] <= max_genes) &
        (adata_.obs['percent_mt'] <= max_percent_mt)
    )
    adata_ = adata_[test, :].copy()

    # Refilter genes, at least expressed in <thr>% of cells
    cell_threshold = thr * adata_.n_obs
    test = (adata_.X > 0).sum(axis=0).A1 >= cell_threshold
    adata_ = adata_[:, test].copy()

    return adata_


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')
path_matrices = os.path.join(path_data, 'STARSolo')

# Gene annotation
gtf = pr.read_gtf(os.path.join(path_data, 'gencode.v45.annotation.gtf'))
# Cell assignment to GBC clones
cell_assignment = pd.read_csv(
    os.path.join(path_data,'CBC_GBC_assignment.csv'), index_col=0
)

# Read matrices, filter genes, add cell QC covariates
SAMPLES = os.listdir(path_matrices)
L = []
for i,sample in enumerate(SAMPLES):

    # Read matrices
    if sample.startswith('.'):
        continue
    adata = sc.read_10x_h5(os.path.join(path_matrices, sample, 'cb_filtered.h5'))
    adata = filter_genes(adata, gtf)

    # Add cell QC covariates
    adata.obs['n_UMIs'] = adata.X.sum(axis=1).A1
    adata.obs['n_genes'] = (adata.X>0).sum(axis=1).A1
    mt_counts = adata.X[:, adata.var_names.str.startswith('MT-')].sum(axis=1).A1
    adata.obs['percent_mt'] = mt_counts / adata.obs['n_UMIs']
    adata.obs['sample'] = sample

    # Add clonal assignments and subset to only cells with assignments
    adata.obs_names = adata.obs_names.map(lambda x: f'{x}_{sample}')
    cell_assignment_ = cell_assignment.query('sample==@sample')
    cells_ = list(set(adata.obs_names) & set(cell_assignment_.index))
    adata = adata[cells_].copy()
    adata.obs['GBC'] = cell_assignment_.loc[cells_, 'GBC']
    L.append(adata)

# Concatenate
adata = sc.concat(L)
del L

# Remove genes with 0 counts
adata = adata[:,adata.X.sum(axis=0).A1>0].copy()
adata.var_names[adata.X.sum(axis=0).A1.argsort()[:-1]][:100]

# Add timepoint and mouse info
adata.obs['timepoint'] = adata.obs['sample'].map(lambda x: x.split('_')[0])
adata.obs['mouse'] = adata.obs['sample'].map(lambda x: x.split('_')[1])

# Store raw counts before filtering cells
adata.layers['raw'] = adata.X.copy()
adata


##


# Cell quality, exploratory
qc_metrics = ['n_UMIs', 'n_genes', 'percent_mt']
adata.obs[qc_metrics+['timepoint']].groupby('timepoint').describe().T

# fig, axs = plt.subplots(1,3,figsize=(8,2.5))
# sns.kdeplot(data=adata.obs, x='n_UMIs', ax=axs[0], hue='timepoint', fill=True)
# axs[0].axvline(1000, color='red', linestyle='--')
# axs[0].axvline(75000, color='red', linestyle='--')
# sns.kdeplot(data=adata.obs, x='n_genes', ax=axs[1], hue='timepoint', fill=True, legend=False)
# axs[1].axvline(1000, color='red', linestyle='--')
# sns.kdeplot(data=adata.obs, x='percent_mt', ax=axs[2], hue='timepoint', fill=True, legend=False)
# axs[2].axvline(.25, color='red', linestyle='--')
# fig.tight_layout()
# plt.show()
# 
# fig, ax = plt.subplots(figsize=(3,3))
# sns.kdeplot(data=adata.obs.loc[adata.obs.query('n_UMIs<100000').index],
#             x='n_UMIs', ax=ax, hue='timepoint', fill=True)
# ax.axvline(1000, color='red', linestyle='--')
# ax.axvline(75000, color='red', linestyle='--')
# fig.tight_layout()
# plt.show()


##


# Cell QC
adata = filter_cells(
    adata,
    min_umis=1000, max_umis=75000,
    min_genes=500, max_genes=20000,
    max_percent_mt=0.25,
    thr=.0005
)
adata
adata.obs['sample'].value_counts()
adata.obs['mouse'].value_counts()
adata.obs['timepoint'].value_counts()
adata.obs[qc_metrics+['timepoint']].groupby('timepoint').describe().T
adata.obs[qc_metrics+['sample']].groupby('sample').median()

# Save QC adata
import anndata
anndata.settings.allow_write_nullable_strings = True
adata.write(os.path.join(path_data, 'adata.h5ad'))


##
