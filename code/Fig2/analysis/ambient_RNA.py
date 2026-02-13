"""
Prepare mtx files for ambient RNA analysis.
"""

import os
import cellbender
import scanpy as sc
import matplotlib.pyplot as plt


##


# Paths
path_main = '/Users/cossa/Desktop/projects/CTCs_scLT'
path_data = os.path.join(path_main, 'data', 'STARSolo')
path_h5ad = os.path.join(path_main, 'data', 'raw_h5ad')

# Write raw cellranger matrices as .h5ad files, run cellbender
for sample in os.listdir(path_data):

    if sample.startswith('.'):
     continue
    print(f'Processing sample: {sample}')
    adata = sc.read_10x_mtx(os.path.join(path_data, sample, 'raw'))
    adata.write(os.path.join(path_h5ad, f'{sample}_raw.h5ad'))
    
    ##
    
    fig, ax = plt.subplots(figsize=(3,3))
    order = adata.X.sum(axis=1).A1.argsort()[::-1]
    sums = adata[order,:].X.sum(axis=1).A1
    ax.plot(sums)
    ax.set(xlabel='Droplets, ordered', ylabel='UMI counts per droplet')
    ax.set_yscale('log')
    ax.set_xscale('log')
    fig.tight_layout()
    fig.savefig(
        os.path.join(
            path_data, sample, 'Elbow_10x.pdf'
        ), dpi=1000
    )

   # path_input = os.path.join(path_h5ad, f'{sample}_raw.h5ad')
   # path_output = os.path.join(path_data, sample, 'cellbender.h5')
   # os.system(
   #     f"cellbender remove-background --input {path_input} --output {path_output} --expected-cells 5000 --epochs 100"
   # )


##


# #!/bin/bash
# #SBATCH --time=3:00:00
# #SBATCH --ntasks=1
# #SBATCH --nodes=1
# #SBATCH --mem=15G
# #SBATCH -J PT_1_late
# #SBATCH -o PT_1_late.stdout
# #SBATCH -e PT_1_late.stderr
# #SBATCH --mail-user=cossa@ebi.ac.uk
# #SBATCH --mail-type=END
# #SBATCH —gres=gpu:1
# 
# mamba activate ambient_RNA
# 
# path_wd=/homes/cossa/preprocess_CTCs/work/cellbender
# path_input=/homes/cossa/preprocess_CTCs/raw_h5ad/PT_1_late_raw.h5ad
# path_output=/homes/cossa/preprocess_CTCs/raw_h5ad/PT_1_late_cb.h5ad
# 
# cd $path_wd
# cellbender remove-background \
#     --cuda \
#     --input $path_input \
#     --output $path_output \
#     --expected-cells 5000 \
#     --epochs 100
