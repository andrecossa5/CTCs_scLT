# Single-cell lineage tracing of Breast Cancer circulating tumor cells (CTCs)

This repo holds all necessary code to reproduce the analysis of the CTCs paper.
This is the template we will serve until the final submission.

```
.
├── code
│   ├── Fig1
│   ├── Fig2
│   │   └── analysis
│   │       ├── ambient_RNA.py
│   │       ├── GRN_inference.py
│   │       ├── preprocessing.py
│   │       ├── QC.py
│   │       └── trajectories.py
│   ├── Fig3
│   ├── Fig4
│   └── Fig5
├── data
├── envs
├── figures
│   ├── ~$Fig_schemes.pptx
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

`data`, `results`, and `figures` are not synchronized with the repo (they are in the `.gitignore` file), but should be the same we are sharing on Google Drive. The `code` is divided per Figure. If a final figure requires complex analysis and visualizations, the necessary code is further splitted into analysis scripts and panels (e.g. `Fig1a.py` script reproduce Fig1a panel). 

To ensure no pain at revision, each figure panel visualization script needs to be fully independent from the others, and must read and write from the shared `data` and `results` folders only.


