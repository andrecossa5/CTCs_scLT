# Single-cell lineage tracing of Breast Cancer circulating tumor cells (CTCs)

This repo holds all necessary code to reproduce the analysis of the CTCs paper.
This is the template we will serve until the final submission.

```
.
├── code
│   ├── Fig1
│   │   ├── clone_calling.py
│   │   ├── Fig1c.py
│   │   ├── Fig1d.py
│   │   ├── QC.py
│   │   └── viz_clones.py
│   ├── Fig2
│   │   ├── analysis
│   │   │   ├── __pycache__
│   │   │   │   ├── clone_calling.cpython-311.pyc
│   │   │   │   └── QC.cpython-312.pyc
│   │   │   ├── ambient_RNA.py
│   │   │   ├── cell_states.py
│   │   │   ├── gene_list.py
│   │   │   ├── GRN_inference.py
│   │   │   ├── preprocessing.py
│   │   │   ├── scVI.py
│   │   │   └── trajectories.py
│   │   └── umaps.py
│   ├── Fig3
│   ├── Fig4
│   ├── Fig5
│   └── Supplementary
├── data
│   ├── adata.h5ad
│   ├── CBC_GBC_assignment.csv
│   ├── CBC_GBC_combos_merged.tsv.gz
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
│   ├── old
│   │   ├── clustered.h5ad
│   │   └── ctc_regulon.csv
│   └── STARSolo
│       ├── CTC_1_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── CTC_2_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── CTC_3_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── CTC_4_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── lung_1_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── lung_2_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── lung_3_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── lung_4_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── PT_1_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── PT_2_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       ├── PT_3_late
│       │   ├── cb_filtered.h5
│       │   ├── Elbow_10x.pdf
│       │   ├── filtered
│       │   │   ├── barcodes.tsv.gz
│       │   │   ├── features.tsv.gz
│       │   │   └── matrix.mtx.gz
│       │   └── raw
│       │       ├── barcodes.tsv.gz
│       │       ├── features.tsv.gz
│       │       └── matrix.mtx.gz
│       └── PT_4_late
│           ├── cb_filtered.h5
│           ├── Elbow_10x.pdf
│           ├── filtered
│           │   ├── barcodes.tsv.gz
│           │   ├── features.tsv.gz
│           │   └── matrix.mtx.gz
│           └── raw
│               ├── barcodes.tsv.gz
│               ├── features.tsv.gz
│               └── matrix.mtx.gz
├── envs
│   └── env_1.yml
├── figures
│   ├── Fig_schemes.pptx
│   ├── Fig.1_clones.pdf
│   ├── Fig.2_transcription.pdf
│   ├── fig.3_metabolomics.pdf
│   ├── fig.4_seahorse.pdf
│   └── fig.5_ATO.pdf
├── LICENSE
├── README.md
└── results
```

`data`, `results`, and `figures` are not synchronized with the repo (they are in the `.gitignore` file), but should **always** be the same that we are sharing on Google Drive.
The `code` folder is divided per Figure. If a figure requires complex analyses and visualizations, the necessary code would be further splitted into multiple scripts and panels (e.g. `Fig1c.py` script reproduce Fig1c panel).
To ensure no pain at revision, each of these script needs to be fully independent from the others, and must read and write from the shared `data` and `results` folders only.
To avoid pain now, each **small** analysis task should require a separate branch, to be merged with master and removed at completion.
The exception to this would be major updates in the whole folder structures, that may be committed *directly* on the master branch.


