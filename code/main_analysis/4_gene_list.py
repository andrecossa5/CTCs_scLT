"""
Curated gene list for HVGs selection.
"""

import os
import pandas as pd


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data')

# Load genes
GENES = []
for subfolder in os.listdir(os.path.join(path_data, 'gene_lists')):
    for gene_list in os.listdir(os.path.join(path_data, 'gene_lists', subfolder)):
        GENES += (
            pd.read_csv(
                os.path.join(path_data, 'gene_lists', subfolder, gene_list),
                sep='\t', header=None
        )
        [0].to_list()
    )

# Write curated gene list
(
    pd.Series(
        list(set([ x.strip(',') for x in GENES ]))
    )
    .to_csv(os.path.join(path_data, 'curated_genes.txt'), index=False, header=False)
)


##
