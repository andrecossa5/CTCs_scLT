"""
Trajectory inference: pseudotime and OT-based (Moslin) analyses.
"""

import os
import anndata
import numpy as np
import scanpy as sc
import cellrank as cr
import plotting_utils as plu
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances
from scipy.sparse import csr_matrix
from moscot.problems.time import LineageProblem
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))


##


# 1. Pseudotime inference with scanpy Diffusion Pseudotime (DPT)
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


# 2. Apply Moslin (Lineage problem) to compute OT couplings

# Compute lineage cost matrix based on GBC (Jaccard distance)
GBC_bin = (
    adata.obs
    .assign(value=1)
    .reset_index()
    .pivot_table(index='index', values='value', columns='GBC', fill_value=0)
)
D = pairwise_distances(GBC_bin.values, metric='jaccard', n_jobs=-1)
adata.obsp['gbc_distances'] = csr_matrix(D)

# Run Moslin
adata.obs['time'] = np.select(
    [adata.obs['timepoint']=='PT',
     adata.obs['timepoint']=='lung',
     adata.obs['timepoint']=='CTC'],
     [0, 1, 2]
)
adata.obs['time'] = adata.obs['time'].astype('category')
lp = LineageProblem(adata)
lp = lp.prepare(
    time_key="time",
    joint_attr="X_scVI",
    lineage_attr={"attr": "obsp", "key": "gbc_distances", "cost": "custom"},
)
lp = lp.solve(alpha=0.75, epsilon=1e-3, tau_a=0.99, batch_size=1000)


##


# 3. Use cellrank GPCCA to build and analyze the transition matrix from Moslin
rtk = cr.kernels.RealTimeKernel.from_moscot(lp)
rtk = rtk.compute_transition_matrix()

# GPCCA
g = cr.estimators.GPCCA(rtk)
g.compute_schur(n_components=20)

# Write out
g.write(os.path.join(path_data, 'gpcca.pkl'))


##


# Choose number of macrostates
# fig = g.plot_spectrum(real_only=True, figsize=(4,4))  # look for the 'elbow' to pick n_states
# fig = g._plot_real_spectrum(10)
# fig.gca().set_title('')
# fig.tight_layout()
# plt.show()

# Compute macrostates
# g = cr.estimators.GPCCA.read(os.path.join(path_data, 'gpcca.pkl'))
# g.compute_macrostates(n_states=4)
# g.adata

# Viz
# sc.pl.umap(adata_, color='macrostates_fwd', legend_loc=False)
# g.plot_macrostates(which='all', discrete=True, legend_loc='right')
# g.plot_coarse_T(annotate=True)   # coarse-grained transition matrix between macrostates


##


# Optional: Gene expression along pseudotime with GAMs
# pseudotime = adata.obs['dpt_pseudotime'].values
# genes = adata[:,adata.var['highly_variable']].var_names
# 
# res = []
# for i, gene in enumerate(genes):
#     gam = LinearGAM().fit(pseudotime, adata[:,gene].X.toarray().flatten())
#     r2 = gam.statistics_['pseudo_r2']['explained_deviance']
#     p = gam.statistics_['p_values'][1]
#     res.append([gene, r2, p])
#     if i%100==0:
#         print(f'Processed {i} genes')
# 
# res = pd.DataFrame(res, columns=['gene', 'effect_size', 'pvalue'])
# 
# res.sort_values('effect_size', ascending=False).head(50)
# 
# 
# 
# sc.tl.rank_genes_groups(adata, groupby='timepoint', method='wilcoxon')
# 
# pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges']).head(10)
# pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
# 
# 
# sc.tl.score_genes(
#     adata,
#     gene_list=adata.var_names[adata.var_names.str.startswith('NDUF')],
#     score_name='score'
# )
# 
# adata.obs.groupby('timepoint')['score'].describe().T
# 
# 
# adata.obs['score'] = (adata.obs['score'] - adata.obs['score'].mean()) / adata.obs['score'].std()
# sc.pl.umap(adata, color='score')


##
