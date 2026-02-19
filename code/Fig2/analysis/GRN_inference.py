"""
Script to analyze pyscenic results 
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import plotting_utils as plu
from gseapy import enrichr
from sklearn.preprocessing import scale
matplotlib.use('macOSX')

#Path
path_main = "/Users/ieo7295/Desktop/data_ctc"
path_results= os.path.join(path_main, "results")

#Data
adata= sc.read_h5ad(os.path.join(path_main, "adata.h5ad"))
regulons = pd.read_csv(os.path.join(path_main, "ctc_regulon.csv"), index_col=0)

#extract gene sets
regulons['gene_set'] = regulons['gene_set'].map(lambda x: x.split(', '))

# Score regulons in anndata with sc.tl.score_genes 
scores = np.zeros((adata.shape[0], regulons.shape[0]))
for i,tf in enumerate(regulons.index):
    genes_set = regulons.loc[tf,'gene_set']
    genes_present = [ gene in adata.var_names for gene in genes_set ]
    test = np.sum(genes_present)>= len(genes_set) * .9
    if test:
        sc.tl.score_genes(adata, gene_list=genes_set)
        scores[:,i] = adata.obs['score'].values

# filter weakly expressed regulons, and z-score their expression score
test = np.sum(scores, axis=0)!=0
scores = pd.DataFrame(scale(scores[:,test]), index=adata.obs_names, columns=regulons.index[test])

# Effect size analysis
adata_reg = sc.AnnData(
    X=scores.values,
    obs=adata.obs[['timepoint']].copy(),
    var=pd.DataFrame(index=scores.columns)
)

def wilcoxon_effects(adata_reg, scores_df, group="CTC", reference=None,
                          groupby="timepoint", padj_max=0.1,
                          min_delta_median=0.2, min_delta_prev=None,
                          prev_threshold=0.0):
    """
    reference=None -> tests group vs rest (scanpy behavior)
    reference="PT" -> group vs PT
    Adds effect sizes:
      - delta_median
      - delta_mean
      - delta_prev (fraction > prev_threshold)
    """
    kwargs = dict(groupby=groupby, method="wilcoxon")
    if reference is None:
        sc.tl.rank_genes_groups(adata_reg, **kwargs)
    else:
        sc.tl.rank_genes_groups(adata_reg, groups=[group], reference=reference, **kwargs)

    res = sc.get.rank_genes_groups_df(adata_reg, group=group).copy()
    res = res.drop(columns=['logfoldchanges'])

    g = adata_reg.obs[groupby].astype(str).values
    mask_g = (g == group)
    if reference is None:
        mask_ref = ~mask_g
        ref_name = "rest"
    else:
        mask_ref = (g == str(reference))
        ref_name = str(reference)

    g_vals   = scores_df.loc[adata_reg.obs_names[mask_g]]
    ref_vals = scores_df.loc[adata_reg.obs_names[mask_ref]]

    delta_median = g_vals.median(axis=0) - ref_vals.median(axis=0)
    delta_mean   = g_vals.mean(axis=0)   - ref_vals.mean(axis=0)
    delta_prev   = (g_vals > prev_threshold).mean(axis=0) - (ref_vals > prev_threshold).mean(axis=0)

    res["delta_median"] = res["names"].map(delta_median)
    res["delta_mean"]   = res["names"].map(delta_mean)
    res["delta_prev"]   = res["names"].map(delta_prev)
    res["reference"]    = ref_name

    # Apply thresholds
    filt = (res["pvals_adj"] < padj_max) & (res["delta_median"] > min_delta_median)
    if min_delta_prev is not None:
        filt &= (res["delta_prev"] > min_delta_prev)

    res_filt = res.loc[filt].sort_values("delta_median", ascending=False)
    return res, res_filt

# thresholds 
padj_max = 0.05
min_delta_median = 0.2   
min_delta_prev = 0.05    

# CTC vs rest
res_all, res_all_filt = wilcoxon_effects(
    adata_reg, scores, group="CTC", reference=None,
    padj_max=padj_max, min_delta_median=min_delta_median,
    min_delta_prev=min_delta_prev, prev_threshold=0.0
)
 
# CTC vs PT
res_pt, res_pt_filt = wilcoxon_effects(
    adata_reg, scores, group="CTC", reference="PT",
    padj_max=padj_max, min_delta_median=min_delta_median,
    min_delta_prev=min_delta_prev, prev_threshold=0.0
)

# CTC vs lung
res_lung, res_lung_filt = wilcoxon_effects(
    adata_reg, scores, group="CTC", reference="lung",
    padj_max=padj_max, min_delta_median=min_delta_median,
    min_delta_prev=min_delta_prev, prev_threshold=0.0
)

# Top regulons by robust effect
top_regulons_all  = res_all_filt.head(30)
top_regulons_pt   = res_pt_filt.head(30)
top_regulons_lung = res_lung_filt.head(30)

top_all  = set(top_regulons_all["names"])
top_pt   = set(top_regulons_pt["names"])
top_lung = set(top_regulons_lung["names"])


## Pathway enrichment for each TF gene set
def run_ora_per_tf(res_de_top, comparison_name):
    """Over-Representation Analysis for each TF individually"""
    results_all = []
    
    for tf in res_de_top['names']:
        if tf not in regulons.index:
            print(f"TF {tf} not found in regulons, skipping.")
            continue
        
        gene_list = regulons.loc[tf, 'gene_set']
        
        # Run enrichment
        enr = enrichr(gene_list=gene_list, gene_sets='MSigDB_Hallmark_2020') #GO_Biological_Process_2025
        ora_df = enr.results.copy()

        #Filter and reformat
        ora_df = ora_df[ora_df['Adjusted P-value'] < 0.1].copy()
        ora_df['TF'] = tf
        ora_df['n_tf_genes'] = len(gene_list)
        
        results_all.append(ora_df)
    
    if results_all:
        final_df = pd.concat(results_all, ignore_index=True)
        final_df = final_df.sort_values(['TF', 'Adjusted P-value'])
        final_df.to_csv(os.path.join(path_results, f'ORA_per_TF_{comparison_name}_hallmark.csv'), index=False)

# Run 
run_ora_per_tf(top_regulons_all, 'all_conditions')
run_ora_per_tf(top_regulons_pt, 'CTC_vs_PT')
run_ora_per_tf(top_regulons_lung, 'CTC_vs_lung')

## Manual annotation
manual_annot = {   

    'RARG(+)': '''
    Gene set : 46

    Top pathways Hallmark: Hypoxia	3/200
                           Glycolysis	3/200
                           Inflammatory Response	3/200

    scores: 65.64017486572266
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.5721968740988068
    delta_mean: 0.550289157088761
    delta_prev: 0.24911671069030106
    reference: rest

    scores: 58.55464553833008
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.5389667064553096
    delta_mean: 0.5116750942339101
    delta_prev: 0.23386466362663383
    reference: PT

    scores: 63.08513641357422
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.6294702732856918
    delta_mean: 0.6283528907042908
    delta_prev: 0.27995085895243005
    reference: lung

    Description: Retinoic acid receptors bind as heterodimers to their 
                 target response elements in response to their ligands, 
                 all-trans or 9-cis retinoic acid, and regulate gene expression 
                 in various biological processes. 

                 Transcriptome analysis revealed that after doxorubicin treatment, 
                 cells with the RARG variant showed significantly more activation of 
                 oxidative phosphorylation and other genes associated with metabolic 
                 stress than control cells(Magdy et al., 2021) 
    ''',

    'CEBPG(+)': '''
    Gene set : 367

    scores: 43.88399887084961
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.35395624569997763
    delta_mean: 0.3949288182233196
    delta_prev: 0.16431417086748695
    reference: rest

    scores: 47.621795654296875
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.40389508034412114
    delta_mean: 0.4467780886647092
    delta_prev: 0.1856658700525382
    reference: PT

    scores: 26.97256088256836
    pvals: 3.102209120943079e-160
    pvals_adj: 1.0857731923300777e-159
    delta_median: 0.2510117792475668
    delta_mean: 0.2901082596845612
    delta_prev: 0.12114872195613596
    reference: lung
    
    Top pathways Hallmark: mTORC1 Signaling	20/200
                           Unfolded Protein Response	15/113
                           Estrogen Response Early	14/200
    Top pathways GO: Sulfur Amino Acid Transport (GO:0000101)	4/7
                     L-alpha-amino Acid Transmembrane Transport (GO:1902475)	7/38
                     L-leucine Transport (GO:0015820)	4/9

    Description: CEBPG transcription factor correlates with antioxidant and DNA repair 
                 genes in normal bronchial epithelial cells but not in individuals with 
                 bronchogenic carcinoma (D'Anna et al., 2005)

                 C/EBP-γ is essential for the induction of the integrated stress response.
                 C/EBP-γ/ATF4 appears to regulate the transcription of genes necessary for glutathione synthesis
                 (e.g., Cth, Slc7a1, Gpx7, Mthfd2, Slc1a5, Gpt2). Importantly, C/EBP-γ was determined to be a key 
                 regulator of glutathione synthesis in both normal and lung adenocarcinoma cells (Renfro et al., 2022)
    ''',

    'TFDP1(+)': '''
    Gene set : 1617
    
    scores: 61.670799255371094
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.846390539110212
    delta_mean: 0.5281153330343921
    delta_prev: 0.21438561789017535
    reference: PT

    scores: 39.315895080566406
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.490764909304703
    delta_mean: 0.32238957604061735
    delta_prev: 0.12199400020963402
    reference: rest

    Top pathways Hallmark: Oxidative Phosphorylation	25/200
                           Pperoxisome	15/104
                           PI3K/AKT/mTOR  Signaling	15/105
    Top pathways GO: DNA Metabolic Process (GO:0006259)	113/302
                     DNA Repair (GO:0006281)	108/302
                     DNA Damage Response (GO:0006974)	108/421

    Description : This gene encodes a member of a family of transcription factors that 
                  heterodimerize with E2F proteins to enhance their DNA-binding activity and
                  promote transcription from E2F target genes. The encoded protein functions as 
                  part of this complex to control the transcriptional activity of numerous genes involved 
                  in cell cycle progression from G1 to S phase.

                  In TNBC, transcription factors TFDP1, CEBPB, and ETS2 target enzymes involved in 
                  oxidative phosphorylation. In endometrial cancer, transcription factors XBP1, TFDP1, E2F1, and FOXM1 
                  primarily regulate enzymes involved in the TCA cycle and one-carbon amino acid metabolism (Srivastava et al., 2025)
    ''',

    'IRF9(+)': '''
    Gene set : 1283

    scores: 74.29396057128906
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.66931895087795
    delta_mean: 0.7160222762946582
    delta_prev: 0.25163250321849484
    reference: rest

    scores: 62.13469696044922
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.5644182839427228
    delta_mean: 0.6150728192360557
    delta_prev: 0.208437311523189
    reference: PT
    
    scores: 78.82821655273438
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.9078932133397356
    delta_mean: 0.9201057355793658
    delta_prev: 0.33895763009388374
    reference: lung

    Top pathways GO: Negative Regulation of Viral Genome Replication (GO:0045071)	16/48
                     Defense Response to Virus (GO:0051607)	33/198
                     Regulation of Viral Entry Into Host Cell (GO:0046596)	9/18 
    Top pathways Hallmark: Interferon Alpha Response	45/97
                     Interferon Gamma Response	57/200
                     Estrogen Response Late	30/200
    
    Description : On the one hand, IRF9 directly regulates the expression of PHB1 and also directly 
                  interacts with AKT to inhibit the phosphorylation of AKT at the Thr308 site but not at the Ser473 site;
                  on the other hand, IRF9 overexpression results in mitochondrial dysfunction by increasing opening of mPTP and 
                  mitochondrial permeability. Given the important regulatory role of PHB1 and AKT in mitochondrial function, our results 
                  indicated that IRF9 may regulate mitochondrial function and energy metabolism to promote the proliferation of PASMCs in 
                  PAH via the IRF9- PHB1-AKT axis.
                  IRF9 is frequently overexpressed in human lung cancer and is associated with decreased patient survival. (Chen et al., 2021)
    ''',

    'ILF2(+)' : '''
    Gene set : 2367

    scores: 67.8245849609375
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.607255429893069
    delta_mean: 0.6150015787462548
    delta_prev: 0.23693422498747374
    reference: rest

    scores: 87.48237609863281
    pvals: 0.0
    pvals_adj: 0.0
    delta_median: 0.8513189815927997
    delta_mean: 0.8233100135341631
    delta_prev: 0.32994454206018015
    reference: PT

    Top pathways Hallmark: Myc Targets V1	95/200
                  E2F Targets	79/200
                  Oxidative Phosphorylation	74/200
    Top pathways GO: Mitochondrial Translation (GO:0032543)	37/94
                  Mitochondrial Respiratory Chain Complex Assembly (GO:0033108)	36/90
                  Mitochondrial Gene Expression (GO:0140053)	38/100

    Description : Chromatin-interacting protein that forms a stable heterodimer with interleukin enhancer-binding factor 
                  3/ILF3 and plays a role in several biological processes including transcription, innate immunity or cell growth.
                  Essential for the efficient reshuttling of ILF3 (isoform 1 and isoform 2) into the nucleus. Together with ILF3, 
                  forms an RNA-binding complex that is required for mitotic progression and cytokinesis by regulating the expression of a cluster of mitotic genes.

                  ILF2 promotes SCLC tumor growth in vitro and in vivo. ILF2 elevates oxidative phosphorylation expression and declines glucose intake 
                  and lactate production(Zhao et al., 2019)

                  ILF2 gene and/or protein are involved in cancer development and progression through multiple mechanisms, such as regulating cell cycle and 
                  apoptosis, participating in tumor metabolism, and maintaining tumor mitochondrial homeostasis. (Liu et al.,2024)
    ''',

}


#Viz 
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
tfs = ['RARG(+)', 'CEBPG(+)', 'TFDP1(+)', 'IRF9(+)', 'ILF2(+)']

categorical_cmap = {'PT': '#1f77b4', 'lung': '#2ca02c', 'CTC': '#ff7f0e'}

for idx, tf in enumerate(tfs):
    data_list = []
    for tp in ['PT', 'lung', 'CTC']:
        tp_scores = scores.loc[adata.obs['timepoint'] == tp, tf].values
        data_list.append(pd.DataFrame({
            'timepoint': [tp] * len(tp_scores),
            'score': tp_scores
        }))
    df_tf = pd.concat(data_list, ignore_index=True)
    
    # Violin
    ax_violin = axes[0, idx]
    plu.violin(
        df_tf,
        x='timepoint',
        y='score',
        x_order=['PT', 'lung', 'CTC'],
        categorical_cmap=categorical_cmap,
        ax=ax_violin,
        kwargs={'inner': 'box'}
    )
    plu.format_ax(
        ax=ax_violin,
        title=tf,
        xlabel='',
        ylabel='Score (z-scored)',
        reduced_spines=True
    )
    
    # Histograms overlay
    pt_scores = scores.loc[adata.obs['timepoint'] == 'PT', tf].values
    lung_scores = scores.loc[adata.obs['timepoint'] == 'lung', tf].values
    ctc_scores = scores.loc[adata.obs['timepoint'] == 'CTC', tf].values
    
    axes[1, idx].hist(pt_scores, bins=50, alpha=0.5, label='PT', density=True, color=categorical_cmap['PT'])
    axes[1, idx].hist(lung_scores, bins=50, alpha=0.5, label='lung', density=True, color=categorical_cmap['lung'])
    axes[1, idx].hist(ctc_scores, bins=50, alpha=0.5, label='CTC', density=True, color=categorical_cmap['CTC'])
    axes[1, idx].set_xlabel('Score (z-scored)')
    axes[1, idx].set_ylabel('Density')
    axes[1, idx].set_title(tf)
    axes[1, idx].legend(fontsize=8)
    axes[1, idx].grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig(os.path.join(path_results, 'regulon_scores_distributions_all.pdf'), dpi=300)

