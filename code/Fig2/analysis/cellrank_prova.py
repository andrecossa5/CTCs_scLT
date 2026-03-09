"""
Trajectory inference: CellRank
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import cellrank as cr
import plotting_utils as plu
import matplotlib.pyplot as plt
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')
path_cospar = os.path.join(path_data, 'cospar')

# Load adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))

# Time info
adata.obs['time_info'] = np.select(
    [adata.obs['timepoint'] == 'PT', 
     adata.obs['timepoint'] == 'lung', 
     adata.obs['timepoint'] == 'CTC'],
    [0, 1, 2]
)
adata.obs['time_info'] = pd.Categorical(
    adata.obs['time_info'], categories=[0, 1, 2], ordered=True
)

# Create transition kernel
rtk = cr.kernels.RealTimeKernel.from_wot(adata, path=path_cospar, time_key='time_info')
rtk.compute_transition_matrix()

# Analyze with GPCCA
g = cr.estimators.GPCCA(rtk)
g.compute_schur()

# Choose number of macrostates
fig = g.plot_spectrum(real_only=True, figsize=(8,4))  # look for the 'elbow' to pick n_states
fig = g._plot_real_spectrum(20)
fig.gca().set_title('')
fig.tight_layout()
fig.savefig(os.path.join(path_cospar, 'schur_spectrum.pdf'))

# Save class
with open(os.path.join(path_cospar, 'gpcca_cospar.pkl'), 'wb') as f:
    pickle.dump(g, f)


##


# Compute macrostates
# g.compute_macrostates(n_states=6, n_cells=100)
# g.adata.obs['macrostates_fwd'].value_counts()
# classes = []
# for i in range(g.adata.obsm['macrostates_fwd_memberships'].shape[0]):
#     memberships = np.array(g.adata.obsm['macrostates_fwd_memberships'][i,:].tolist())
#     maximum = memberships.max()
#     argmax = memberships.argmax()
#     if maximum > 0.5:
#         classes.append(argmax)
#     else:
#         classes.append(np.nan)

# adata.obs['classes'] = classes
# adata.obs['classes'] = pd.Categorical(adata.obs['classes'], categories=range(6), ordered=True)
# adata.obs['classes'].value_counts()
#.value_counts()  # check if any cell is assigned to multiple macrostates

# sc.pl.draw_graph(adata, color='classes', legend_loc=False)

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


##
