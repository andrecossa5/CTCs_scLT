"""
Fig2 i. Clones annotation.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import plotting_utils as plu
plu.set_rcParams()


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Read adata
adata = sc.read(os.path.join(path_data, 'adata.h5ad'))
adata

# Calculate frequencies
df_ = (
    adata.obs
    .groupby('sample')
    ['GBC'].value_counts()
    .reset_index(name='n')
    .assign(
        f=lambda x: x['n'] / x.groupby('sample')['n'].transform('sum')
    )
    .assign(dataset=lambda x: 'M' + x['sample'].map(lambda s: s.split('_')[1]))
    .assign(timepoint=lambda x: x['sample'].map(lambda s: s.split('_')[0]))
)

# Only pro e non-met clones
df_ = df_.merge(
    (
        adata.obs[['GBC', 'metastatic_status']]
        .drop_duplicates().reset_index(drop=True)
        .query('metastatic_status!="other"')
    ),
    on='GBC',
    how='left'
)
df_ = df_.dropna()

# Re-arrange
df_wide = (
    df_.pivot(index=['GBC', 'dataset'], columns='timepoint', values='f').fillna(0)
    .loc[lambda x: x.sum(axis=1)>0]
    .reset_index()
    .melt(
        id_vars=['GBC', 'dataset'], 
        value_vars=['PT', 'lung', 'CTC'], 
        var_name='timepoint', 
        value_name='f'
    )
    .assign(GBC_dataset=lambda x: x['GBC'].astype(str) + '_' + x['dataset'])
)

# Load colors
with open(os.path.join(path_data, 'colors.pkl'), 'rb') as f:
    COLORS = pickle.load(f)

# Order clones: highest CTC frequency on top
tp_order = ['PT', 'lung', 'CTC']
ctc_freq = (
    df_wide[df_wide['timepoint'] == 'CTC']
    .set_index('GBC_dataset')['f']
    .reindex(df_wide['GBC_dataset'].unique(), fill_value=0)
    .sort_values(ascending=True)   # ascending → lowest CTC at bottom (y=0)
)
clone_order = ctc_freq.index.tolist()
n_clones = len(clone_order)
df_wide['GBC_dataset'] = pd.Categorical(df_wide['GBC_dataset'], categories=clone_order, ordered=True)
df_wide['timepoint']   = pd.Categorical(df_wide['timepoint'],   categories=tp_order,    ordered=True)


##


# Figure
fig, ax = plt.subplots(figsize=(1.5, 4))

for tp in tp_order:
    sub = df_wide[df_wide['timepoint'] == tp]
    ax.scatter(
        x=sub['timepoint'].cat.codes,
        y=sub['GBC_dataset'].cat.codes,
        s=sub['f'] * 100,
        color=COLORS['timepoint'][tp],
        alpha=0.7,
        edgecolors='w',
        linewidths=0.3,
        zorder=3,
        clip_on=False,
    )

# Axes formatting
n_long = df_.query('metastatic_status=="Pro-metastatic"')['GBC'].unique().size
n_PT = df_.query('metastatic_status=="Non-metastatic"')['GBC'].unique().size
plu.format_ax(
    ax=ax, ylabel='Lentiviral clones', yticks=[], xticks=tp_order, rotx=90, 
    title=f'n longitudinal: {n_long}\nn PT-restricted: {n_PT}'
)
ax.set_xlim(-0.5, len(tp_order) - 0.5)
ax.set_ylim(-2.5, n_clones + 5.5)

# Size legend
for f_val, label in [(0.05, '5%'), (0.15, '15%'), (0.30, '30%')]:
    ax.scatter([], [], s=f_val * 100, color='grey', alpha=0.85,
               edgecolors='w', linewidths=0.3, label=label)
ax.legend(
    title='Prevalence', fontsize=6, title_fontsize=6,
    bbox_to_anchor=(1, 0), loc='lower right', frameon=False,
)

fig.subplots_adjust(left=0.3, right=0.8, top=0.85, bottom=0.2)
fig.savefig(os.path.join(path_main, 'figures', 'Fig2', 'Fig2i.pdf'))


##

