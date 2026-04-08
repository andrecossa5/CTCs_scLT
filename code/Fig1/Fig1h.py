"""
Fig1 h. Fishplot longitudinal clones
"""

import os
import numpy as np
import scanpy as sc
import plotting_utils as plu
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import norm
matplotlib.use('macOSX')


# Paths
path_main = "/Users/ieo7295/Desktop/data_ctc"
path_results= os.path.join(path_main, "results")

# Utils
def gaussian_smooth(x, y, grid, sd):
    weights = np.transpose([norm.pdf(grid, m, sd) for m in x])
    weights = weights / weights.sum(0)
    return (weights * y).sum(1)

#Read adata
adata = sc.read(os.path.join(path_main, 'adata.h5ad'))
meta = adata.obs.copy()

##

# Calculate clonal frequencies
df_ = (
    meta.groupby('sample')
    .apply(lambda x: x['GBC'].value_counts(normalize=True))
    .reset_index()
    .pivot_table(columns='sample', index='GBC', values='GBC', fill_value=0)
)


d = {
    'Dataset 1, late' : ['PT_1_late', 'lung_1_late', 'CTC_1_late'],
    'Dataset 2, late' : ['PT_2_late', 'lung_2_late', 'CTC_2_late'],
    'Dataset 3, late' : ['PT_3_late', 'lung_3_late', 'CTC_3_late'],
    'Dataset 4, late' : ['PT_4_late', 'lung_4_late', 'CTC_4_late']
}

# Palettes for datasets 
palettes = [
    sc.pl.palettes.default_102[::-1],  
    sc.pl.palettes.default_102,        
    sc.pl.palettes.default_102[75:] , 
    sc.pl.palettes.default_102[60:]   
]
palette_idx = 0

#Fishplot
for dataset, sample_list in d.items():
    df_long = df_[sample_list]
    df_long = df_long.sort_values(sample_list[-1], ascending=False)

    color_palette = palettes[palette_idx]
    offset = (palette_idx * 25) % len(color_palette)  
    selected_colors = [color_palette[(offset + i) % len(color_palette)] for i in range(len(df_long))]
    colors = {k: v for k, v in zip(df_long.index, selected_colors)}
    palette_idx += 1
    x = np.arange(df_long.shape[1])
    y = [ np.array(x) for x in df_long.values.tolist() ]
    grid = np.linspace(-1.3, 4, num=500)
    y_smoothed = [ gaussian_smooth(x, y_, grid, .35) for y_ in y ]

    fig, ax = plt.subplots(figsize=(10,5))
    ax.stackplot(grid, y_smoothed, baseline="sym", colors=colors.values())
    plu.format_ax(ax, xticks=['PT', 'lung', 'CTC'], yticks=[])
    ax.spines[['right', 'bottom', 'top', 'left']].set_visible(False)
    for l in x:    
        ax.axvline(x=l, color='k', linewidth=.5, linestyle='dashed')
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, f'Fig_1h_{dataset}.pdf'), dpi=400)
