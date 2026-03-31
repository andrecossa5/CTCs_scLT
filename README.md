# Single-cell lineage tracing of Breast Cancer circulating tumor cells (CTCs)

This repo holds all necessary code to reproduce the analysis of the CTCs paper.
This is the template we will serve until the final submission.

```
.
├── code
│   ├── Fig1
│   │   ├── Fig1c.py
│   │   └── Fig1d.py
│   ├── Fig2
│   │   ├── Fig2ae.py
│   │   ├── FIg2f.py
│   │   ├── Fig2gh.py
│   │   ├── Fig2i.py
│   │   ├── Fig2l.py
│   │   ├── Fig2m.py
│   │   └── Fig2n.py
│   └── main_analysis
│       ├── __pycache__
│       │   └── utils.cpython-311.pyc
│       ├── 1_ambient_RNA.py
│       ├── 10_clones_annotation.py
│       ├── 2_clone_calling.py
│       ├── 3_QC.py
│       ├── 4_gene_list.py
│       ├── 5_preprocessing_I.py
│       ├── 6_scVI.py
│       ├── 7_preprocessing_II.py
│       ├── 8_Hotspot.py
│       ├── 9_cell_state_annotation.py
│       └── utils.py
├── data
│   ├── adata.h5ad
│   ├── CBC_GBC_assignment.csv
│   ├── CBC_GBC_combos_merged.tsv.gz
│   ├── cell_state_markers_filtered.csv
│   ├── cell_state_markers_full.csv
│   ├── cells_meta.csv
│   ├── colors.pkl
│   ├── connectivities.npz
│   ├── curated_genes.txt
│   ├── gencode.v45.annotation.gtf
│   ├── gene_lists
│   │   ├── Breast
│   │   │   ├── breast_atlas_top_markers.txt
│   │   │   ├── CTCs_11_gene_list_markers.txt
│   │   │   ├── CTCs_11_gene_list_top50.txt
│   │   │   ├── CTCs_7_cluster3_list.txt
│   │   │   ├── CTCs_7_list_markers_all.txt
│   │   │   ├── CTCs_8_gene_list_oxphos.txt
│   │   │   └── CTCs_8_top50_gene_list.txt
│   │   └── Other
│   │       ├── CTCs_1_gene_list_pdac.txt
│   │       ├── CTCs_10_gene_list_mel.txt
│   │       ├── CTCs_15_gene_list_hcc.txt
│   │       └── Peer.txt
│   ├── hotspot_local_correlation_z.csv
│   ├── hotspot_modules.csv
│   ├── hotspot.pkl
│   ├── modules_labels.csv
│   ├── modules_labels.numbers
│   ├── STARSolo
│   │   ├── CTC_1_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── CTC_2_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── CTC_3_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── CTC_4_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── lung_1_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── lung_2_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── lung_3_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── lung_4_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── PT_1_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── PT_2_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   ├── PT_3_late
│   │   │   ├── cb_filtered.h5
│   │   │   ├── Elbow_10x.pdf
│   │   │   ├── filtered
│   │   │   │   ├── barcodes.tsv.gz
│   │   │   │   ├── features.tsv.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   └── raw
│   │   │       ├── barcodes.tsv.gz
│   │   │       ├── features.tsv.gz
│   │   │       └── matrix.mtx.gz
│   │   └── PT_4_late
│   │       ├── cb_filtered.h5
│   │       ├── Elbow_10x.pdf
│   │       ├── filtered
│   │       │   ├── barcodes.tsv.gz
│   │       │   ├── features.tsv.gz
│   │       │   └── matrix.mtx.gz
│   │       └── raw
│   │           ├── barcodes.tsv.gz
│   │           ├── features.tsv.gz
│   │           └── matrix.mtx.gz
│   └── X_scVI.csv
├── envs
│   ├── cospar.yml
│   ├── ctcs.yml
│   └── hotspot.yml
├── figures
│   ├── ~$Fig2_andre.pptx
│   ├── Fig_schemes.pptx
│   ├── Fig.1_clones.pdf
│   ├── fig.3_metabolomics.pdf
│   ├── fig.4_seahorse.pdf
│   ├── fig.5_ATO.pdf
│   ├── Fig2
│   │   ├── CGCGTCACACTGTCGGGC_Fig2i_clone_trajectories.pdf
│   │   ├── CTGCGGTTTCGTTAACGC_Fig2i_clone_trajectories.pdf
│   │   ├── Fig2ac.pdf
│   │   ├── Fig2d.pdf
│   │   ├── Fig2e.pdf
│   │   ├── Fig2f.pdf
│   │   ├── Fig2g.pdf
│   │   ├── Fig2h.pdf
│   │   ├── Fig2i.pdf
│   │   ├── Fig2m.pdf
│   │   └── Fig2n.pdf
│   ├── Fig2_andre.pdf
│   ├── Fig2_andre.pptx
│   └── Fig2.zip
├── LICENSE
└── README.md
```

`data`, and `figures` are not synchronized with the repo (they are in the `.gitignore` file), but should **always** be the same that we are sharing on Google Drive.
The `code` folder is divided per Figure. If a figure requires complex analyses and visualizations, the necessary code would be further splitted into multiple scripts and panels (e.g. `Fig1c.py` script reproduce Fig1c panel).
To ensure no pain at revision, each of these script needs to be fully independent from the others, and must read and write from the shared `data` and `results` folders only.
To avoid pain now, each **small** analysis task should require a separate branch, to be merged with master and removed at completion.
The exception to this would be major updates in the whole folder structures, that may be committed *directly* on the master branch.


