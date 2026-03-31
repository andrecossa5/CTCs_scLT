"""
Fig2 m. Longitudinal vs PT-restricted clones.
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import plotting_utils as plu
from adjustText import adjust_text
plu.set_rcParams()


##


# Source utils: insert your <code/main_analysis> path here
sys.path.insert(0, '/Users/cossa/Desktop/projects/CTCs_scLT/code/main_analysis')
from utils import *


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Subset cells
adata_ = adata[adata.obs['metastatic_status']!='other'].copy()
adata_.obs['metastatic_status'].value_counts()

# Subset: genes
adata_ = adata_[:, np.sum(adata_.layers['raw']>0, axis=0)>5000].copy()

# DE
sc.tl.rank_genes_groups(adata_, groupby='metastatic_status', method='wilcoxon', pts=True, key_added='met_markers')
res_de = format_rank_genes_groups(adata_, key='met_markers')
res_de_filtered = format_rank_genes_groups(adata_, key='met_markers', 
                                           filter_results=True, min_pct_group=.5, max_pct_rest=1, min_log2FC=1, max_pval_adj=.1)

# Annotate
df_ = res_de.query('group=="Pro-metastatic"').set_index('gene')
df_['group'] = np.select(
    [( df_['pval_adj']<.1) & (df_['log2FC']>1),
     ( df_['pval_adj']<.1) & (df_['log2FC']<-1)],
    ['up', 'down'],
    default='other'
)


##


# Viz
fig, ax = plt.subplots(figsize=(2,2))

plu.scatter(df=df_.query('group=="other"'), x='log2FC', y='score', 
    color='grey', alpha=0.1, size=1, ax=ax
)
plu.scatter(df=df_.query('group=="up"'), x='log2FC', y='score', 
    color='red', alpha=0.8, size=15, ax=ax, kwargs={'edgecolor':'black', 'linewidth':0.5}
)
plu.scatter(df=df_.query('group=="down"'), x='log2FC', y='score', 
    color='blue', alpha=0.8, size=15, ax=ax, kwargs={'edgecolor':'black', 'linewidth':0.5}
)
plu.format_ax(ax, xlabel='log2FC', ylabel='Re-scaled adj-pval', reduced_spines=True)

_top_up   = df_.query('group=="up"').nlargest(5, 'log2FC')
_top_down = df_.query('group=="down"').nsmallest(5, 'log2FC')
_to_label = pd.concat([_top_up, _top_down])

texts = []
for gene, row in _to_label.iterrows():
    texts.append(ax.text(row['log2FC'], row['score'], gene, fontsize=5))
    
adjust_text(
    texts, ax=ax,
    arrowprops=dict(arrowstyle='-', color='black', lw=0.4),
    expand=(1.2, 1.4),
)

fig.tight_layout()
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2m.pdf'))


##