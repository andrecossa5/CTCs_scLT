"""
Create loom file for pyscenic (raw)
"""

import os 
import pandas as pd
import numpy as np
import loompy as lp
import scanpy as sc

#Files required by pyscenic 
# "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
# "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
# "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
# "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
# "allTFs_hg38.txt"
# "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

#Paths
path_main = "/Users/ieo7295/Desktop/data_CTC"
path_tfs=os.path.join(path_main, "allTFs_hg38.txt")

#output path
f_loom_path_scenic= os.path.join(path_main,"CTC_scenic.loom")

#Data
adata= sc.read(os.path.join(path_main, "adata.h5ad" ))
tfs= [tf.strip() for tf in open(path_tfs)]

# inspect that our anndata has transcription factors listed in our main annotations.
print(
    f"{np.sum(adata.var.index.isin(tfs))} out of {len(tfs)} TFs are found in the object"
)

#select raw counts
adata.X=adata.layers['raw'].copy()

#create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


##
#!/bin/bash
#SBATCH --job-name=scenic_p
#SBATCH --output=scenic_p.%j.out 
#SBATCH --error=scenic_p.%j.err 
#SBATCH --mail-user=noemi.bulla@ieo.it
#SBATCH --nodes=1                         # PBS: select=1
#SBATCH --cpus-per-task=5                 # PBS: ncpus=1
#SBATCH --mem=56G
#SBATCH --mail-type=END,FAIL 
#SBATCH --time=48:00:00 

# cd /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/code || exit

# source ~/.bashrc

# mamba activate pyscenic

# pyscenic grn \
#  /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/resources/CTC_scenic.loom \
#  /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/old_scenic/resources/allTFs_hg38.txt \
#  -o /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/results/matr_adjacencies_raw.csv \
#  --num_workers 5

# mamba deactivate


##
#!/bin/bash
#SBATCH --job-name=scenic_ctx
#SBATCH --output=scenic_ctx.%j.out 
#SBATCH --error=scenic_ctx.%j.err 
#SBATCH --mail-user=noemi.bulla@ieo.it
#SBATCH --nodes=1                         # PBS: select=1
#SBATCH --cpus-per-task=4                 # PBS: ncpus=1
#SBATCH --mem=42G
#SBATCH --mail-type=END,FAIL 
#SBATCH --time=48:00:00  

# cd /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/code || exit

# source ~/.bashrc

# mamba activate pyscenic


# pyscenic ctx /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/results/matr_adjacencies_raw.csv \
#     /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/old_scenic/resources/hg38_*_full_tx_v10_clust.genes_vs_motifs.*.feather \
#     --annotations_fname /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/old_scenic/resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
#     --expression_mtx_fname /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/resources/CTC_scenic.loom \
#     --output /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/results/ctc_regulons_masked.csv \
#     --mask_dropouts \
#     --num_workers 4

# mamba deactivate

##
#!/bin/bash
#SBATCH --job-name=scenic_auc_mask
#SBATCH --output=scenic_auc_mask.%j.out 
#SBATCH --error=scenic_auc_mask.%j.err 
#SBATCH --mail-user=noemi.bulla@ieo.it
#SBATCH --nodes=1                         # PBS: select=1
#SBATCH --cpus-per-task=4                 # PBS: ncpus=1
#SBATCH --mem=52G
#SBATCH --mail-type=END,FAIL 
#SBATCH --time=48:00:00 

# cd /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/code || exit

# source ~/.bashrc

# mamba activate pyscenic

# pyscenic aucell \
#     /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/resources/CTC_scenic.loom \
#     /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/results/ctc_regulons_masked.csv \
#     --output /hpcnfs/scratch/PGP/nbulla/albe_breast_ctc/CTC_scenic/scenic_paper/results/ctc_aucell_masked.loom \
#     --num_workers 4

# mamba deactivate