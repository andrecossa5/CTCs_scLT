"""
Cell assignment to lentiviral clones (GBC).
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


##


def filter_and_pivot(df_combos, umi_treshold,
                     max_ratio_treshold, normalized_abundance_treshold):
    """
    Filter a CBC-GBC-UMI table and return it in large format.
    """

    test = (df_combos['umi'] >= umi_treshold) & \
        (df_combos['max_ratio'] >= max_ratio_treshold) & \
        (df_combos['normalized_abundance'] >= normalized_abundance_treshold)
    df_combos['status'] = np.where(test, 'supported', 'unsupported')
    M = (
        df_combos
        .query('status=="supported"')
        .pivot_table(index='CBC', columns='GBC', values='umi')
    )
    M[M.isna()] = 0

    return M, df_combos


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Args
umi_treshold=3
max_ratio_treshold=.5
normalized_abundance_treshold=.5

# Read CBC-GBC combos
df = pd.read_csv(
    os.path.join(path_data, 'CBC_GBC_combos_merged.tsv.gz'),
    sep='\t', index_col=0
)
L = []
for sample in df['sample'].unique():
    df_sample = df.query('sample==@sample')
    M, _ = filter_and_pivot(df_sample,
        umi_treshold=umi_treshold,
        max_ratio_treshold=max_ratio_treshold,
        normalized_abundance_treshold=normalized_abundance_treshold
    )
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    filtered_M = M.loc[unique_cells]
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC')
        .assign(sample=sample)
    )
    cells_df.index = cells_df.index.map(lambda x: f'{x}_{sample}')
    L.append(cells_df)

df_assignment = pd.concat(L)
df_assignment.to_csv(os.path.join(path_data, 'CBC_GBC_assignment.csv'))


##
