"""
Trajectory inference: pseudotime and OT-based (Moslin) analyses.
"""

import os
import anndata
from click import clear
import numpy as np
import pandas as pd
import scanpy as sc
import plotting_utils as plu
import matplotlib.pyplot as plt
from pygam import LinearGAM
from moscot.problems.time import LineageProblem
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))


##


# 1. Pseudotime inference with Diffusion Pseudotime (DPT)
np.random.seed(1234)
test = (adata.obsm['X_umap'][:,0]>8) & \
       (adata.obsm['X_umap'][:,0]<12) & \
       (adata.obsm['X_umap'][:,1]>3.2) & \
       (adata.obsm['X_umap'][:,1]<3.5)
root = np.random.choice(adata.obs_names[test])
adata.uns["iroot"] = adata.obs_names.get_loc(root)
sc.tl.dpt(adata, n_dcs=10)

# Viz
fig, ax = plt.subplots(figsize=(3,3))
sc.pl.umap(adata, color='dpt_pseudotime', ax=ax,
           show=False, frameon=False, legend_loc=None, size=2.5, title='')
fig.tight_layout()
plt.show()

# FA-layout visualization
sc.tl.draw_graph(adata, layout='fa',
                 init_pos='paga', root=root, n_jobs=-1, random_state=1234)
fig, ax = plt.subplots(figsize=(3,3))
sc.pl.draw_graph(adata, color='timepoint', ax=ax,
           show=False, frameon=False, legend_loc=None, size=2.5, title='')
fig.tight_layout()
plt.show()

# Write anndata
anndata.settings.allow_write_nullable_strings = True
adata.write(os.path.join(path_data, 'adata.h5ad'))


##


# Gene expression along pseudotime with GAMs
pseudotime = adata.obs['dpt_pseudotime'].values
genes = adata[:,adata.var['highly_variable']].var_names

res = []
for i, gene in enumerate(genes):
    gam = LinearGAM().fit(pseudotime, adata[:,gene].X.toarray().flatten())
    r2 = gam.statistics_['pseudo_r2']['explained_deviance']
    p = gam.statistics_['p_values'][1]
    res.append([gene, r2, p])
    if i%100==0:
        print(f'Processed {i} genes')

res = pd.DataFrame(res, columns=['gene', 'effect_size', 'pvalue'])

res.sort_values('effect_size', ascending=False).head(50)



sc.tl.rank_genes_groups(adata, groupby='timepoint', method='wilcoxon')

pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges']).head(10)
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)


sc.tl.score_genes(
    adata,
    gene_list=adata.var_names[adata.var_names.str.startswith('NDUF')],
    score_name='score'
)

adata.obs.groupby('timepoint')['score'].describe().T


adata.obs['score'] = (adata.obs['score'] - adata.obs['score'].mean()) / adata.obs['score'].std()
sc.pl.umap(adata, color='score')





