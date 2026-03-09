"""
Trajectory inference: CoSpar
"""

import os
import numpy as np
import pandas as pd
import cospar as cs
import cospar.tmap._utils as cs_utils
from anndata import AnnData
from scipy.sparse import csr_matrix, load_npz


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')
path_cospar = os.path.join(path_data, 'cospar')

# Re-create adata
conn = load_npz(os.path.join(path_data, 'connectivities.npz'))
obs = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
adata = AnnData(obs=obs, obsp={'connectivities': conn})

# Re-shape clonal info
adata.obs['GBC'] = adata.obs['GBC'].astype(str) 
adata.obsm['X_clone'] = csr_matrix(
    pd.get_dummies(adata.obs['GBC']).astype(np.int8)
)

# Fix time info
adata.obs['time_info'] = np.select(
    [adata.obs['timepoint'] == 'PT', 
     adata.obs['timepoint'] == 'lung', 
     adata.obs['timepoint'] == 'CTC'],
    [0, 1, 2]
)
order = [0, 1, 2]
cs.hf.update_time_ordering(adata, updated_ordering=order)

# Set params
cs.settings.verbosity = 2
cs.settings.data_path = path_cospar  
adata.uns['data_des'] = ['Cospar']      
trunca_threshold = 0.001    
temp_str = "0" + str(trunca_threshold)[2:]  # → "0001"
CoSpar_KNN = 30
fname = os.path.join(
    cs.settings.data_path,
    f"Cospar_Similarity_matrix_with_all_cell_states_kNN{CoSpar_KNN}_Truncate{temp_str}"
)

# Generate CosPar similarity matrix (S) from precomputed connectivities
S = cs_utils.generate_similarity_matrix(
    adata,
    fname,
    round_of_smooth=15,
    neighbor_N=CoSpar_KNN,
    truncation_threshold=trunca_threshold,
    use_existing_KNN_graph=True,   # <-- uses adata.obsp['connectivities'] as-is
    compute_new_Smatrix=True
)

# Infer Cospar transition map T from multi-time clonal data, extend to all cells
adata_ = cs.tmap.infer_Tmap_from_multitime_clones(
    adata,
    clonal_time_points=order,
    smooth_array=[15, 10, 5],    
    CoSpar_KNN=CoSpar_KNN,  
    compute_new=False,
    extend_Tmap_space=True
)

##

# Write out separate coupling matrices
Tmap_t1 = adata_.uns['Tmap_cell_id_t1']
Tmap_t2 = adata_.uns['Tmap_cell_id_t2']

for tp_src, tp_tgt in zip(order[:-1], order[1:]):
    src_cells = adata_.obs_names[adata_.obs['time_info'] == tp_src]
    tgt_cells = adata_.obs_names[adata_.obs['time_info'] == tp_tgt]
    row_pos = np.where(np.isin(adata_.obs_names[Tmap_t1], src_cells))[0]
    col_pos = np.where(np.isin(adata_.obs_names[Tmap_t2], tgt_cells))[0]
    tmap = AnnData(
        X=adata_.uns['transition_map'][np.ix_(row_pos, col_pos)],
        obs=pd.DataFrame(index=adata_.obs_names[Tmap_t1[row_pos]]),
        var=pd.DataFrame(index=adata_.obs_names[Tmap_t2[col_pos]]),
    )
    tmap.write(os.path.join(path_cospar, f'tmap_{tp_src}_{tp_tgt}.h5ad'))


##