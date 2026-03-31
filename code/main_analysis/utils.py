"""
Utils for main analysis.
"""

import numpy as np
import pandas as pd
import gseapy
from typing import List
from anndata import AnnData
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse
from scipy.cluster.hierarchy import linkage, leaves_list
from pygam import LinearGAM, s as gam_s


##


# Utils
def filter_genes(adata, gtf):
    """
    Filter non-coding genes and pseudogenes, and add meta information.
    """
    genes = (
        gtf.df
        .query("Feature=='gene' and gene_type=='protein_coding'")
        [["gene_id", "gene_name", "gene_type"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # Annotate and filter genes
    genes['gene_id'] = genes['gene_id'].map(lambda x: x.split('.')[0])
    adata = adata[:, adata.var['gene_ids'].isin(genes['gene_id'])].copy()
    test = (
        (adata.var_names.str.startswith('LC0')) | \
        (adata.var_names.str.startswith('LC1')) | \
        (adata.var_names.str.startswith('AC0')) | \
        (adata.var_names.str.startswith('AC1')) | \
        (adata.var_names.str.startswith('RPS')) | \
        (adata.var_names.str.startswith('RPL'))
    )
    adata = adata[:,~test].copy()

    return adata


##


def filter_cells(adata, thr=.001,
                 min_umis=1000, max_umis=50000,
                 min_genes=1000, max_genes=10000,
                 max_percent_mt=0.15):
    """
    Filter cells based on QC metrics.
    """
    adata_ = adata.copy()
    test = (
        (adata_.obs['n_UMIs'] >= min_umis) &
        (adata_.obs['n_UMIs'] <= max_umis) &
        (adata_.obs['n_genes'] >= min_genes) &
        (adata_.obs['n_genes'] <= max_genes) &
        (adata_.obs['percent_mt'] <= max_percent_mt)
    )
    adata_ = adata_[test, :].copy()

    # Refilter genes, at least expressed in <thr>% of cells
    cell_threshold = thr * adata_.n_obs
    test = np.asarray((adata_.X > 0).sum(axis=0)).ravel() >= cell_threshold
    adata_ = adata_[:, test].copy()

    return adata_


##


def format_rank_genes_groups(
    adata: AnnData,
    key: str = 'rank_genes_groups',
    filter_results: bool = False,
    rank_by: str = 'log2FC',
    max_pval_adj: float = .05, 
    min_log2FC: float = 1, 
    min_pct_group: float = .5, 
    max_pct_rest: float = .5
    ) -> pd.DataFrame :

    L = []
    cats = adata.uns[key]['names'].dtype.names
    group_col = adata.uns[key]['params']['groupby']
    
    for cat in cats:
        genes = adata.uns[key]['names'][cat]
        df_ = pd.DataFrame({
            'gene' : genes,
            'score' : adata.uns[key]['scores'][cat],
            'log2FC' : adata.uns[key]['logfoldchanges'][cat],
            'pval_adj' : adata.uns[key]['pvals_adj'][cat]
        })
        df_ = df_.dropna(subset=['gene'])
        df_ = df_[df_['gene'].isin(adata.var_names)]
        df_['pct_group'] = adata.uns[key]['pts'][cat].loc[df_['gene']].values
        df_['pct_rest'] = adata.uns[key]['pts_rest'][cat].loc[df_['gene']].values
        df_['group'] = cat
        cell_group = adata.obs_names[adata.obs[group_col] == cat].to_list()
        cell_rest = adata.obs_names[adata.obs[group_col] != cat].to_list()
        def _to_dense(X):
            return X.toarray() if issparse(X) else np.asarray(X)
        df_['mean_exp_group'] = _to_dense(adata[cell_group, df_['gene'].to_list()].X).mean(axis=0)
        df_['mean_exp_rest'] = _to_dense(adata[cell_rest, df_['gene'].to_list()].X).mean(axis=0)
        
        if filter_results:
            df_ = df_.query('log2FC>=@min_log2FC and pct_group>@min_pct_group and pct_rest<=@max_pct_rest')
            df_ = df_.query('pval_adj<=@max_pval_adj')
            df_ = df_.sort_values(rank_by, ascending=False)

        L.append(df_)

    df_results = pd.concat(L)

    return df_results


##


def order_groups(adata, groupby=None, obsm_key='X_pca', n_dims=15):
    """
    Sort groups in adata.obs[groupby].
    """

    X_pca_agg = (
        pd.DataFrame(adata.obsm[obsm_key][:,:n_dims])
        .assign(group=adata.obs[groupby].astype('str').values)
        .groupby('group').median()
    )
    sorted_groups = (
        X_pca_agg.index.values
        [leaves_list(linkage(pairwise_distances(X_pca_agg)))]
    )

    return sorted_groups


##


def get_top_markers(df, groupby='group', sort_by='log2FC', ascending=False, order_groups=None, ntop=3):
    """
    Get top markers for groups in adata.obs[groupby], according to order_groups.
    """

    top_markers = []
    for cat in order_groups:
        markers = (
            df.loc[df[groupby]==cat]
            .sort_values(sort_by, ascending=ascending)
            .head(ntop)
            ['gene'].to_list()
        )
        top_markers.extend(markers)

    return top_markers



##


def run_GSEA(
    ranked_list: pd.Series, 
    collections: str|List[str] = 'MSigDB_Hallmark_2020',
    max_pval_adj: float = .01,
    min_size_set: int = 15,
    max_size_set: int =  1000
    ) -> pd.DataFrame :
    """
    Fast GSEA.
    """

    # names = pd.Series(gseapy.get_library_name())
    # names[names.str.contains('Hall')]

    L = collections if isinstance(collections, list) else [collections]
    results = gseapy.prerank(
        rnk=ranked_list,
        gene_sets=L,
        threads=-1,
        min_size=min_size_set,
        max_size=max_size_set,
        permutation_num=200, 
        outdir=None, 
        seed=1234,
        verbose=True,
    )
    df = (
        results.res2d
        [[ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]]
        .rename(columns={'FDR q-val' : 'pval_adj'})
        .query('pval_adj<=@max_pval_adj')
        .sort_values('NES', ascending=False)
    )

    return results, df


##


def run_ORA(
    gene_list: List[str], 
    collections: str|List[str] = 'MSigDB_Hallmark_2020',
    max_pval_adj: float = .01
    ) -> pd.DataFrame :

    L = collections if isinstance(collections, list) else [collections]
    results = gseapy.enrichr(
        gene_list=gene_list,
        gene_sets=L,
        cutoff=.1,
        no_plot=True,
        outdir=None, 
    ).results
    df = (
        results
        [[ 'Term', 'Overlap', 'Odds Ratio', 'Adjusted P-value', 'Genes' ]]
        .rename(columns={'Adjusted P-value' : 'pval_adj'})
        .query('pval_adj<=@max_pval_adj')
        .sort_values('Odds Ratio', ascending=False)
    )

    return results, df


##


def fit_gam_association(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    n_splines: int = 6,
    spline_order: int = 3,
    lam: float = None,
    n_test: int = 200,
) -> dict:
    """
    Test the association between a continuous response (y_col) and a continuous
    predictor (x_col) using a GAM with normal distribution and identity link.

    Appropriate for z-scored covariates (e.g. scaled gene module scores vs pseudotime).
    Mirrors CellRank's GAM backend: cubic B-splines, derivative+L2 penalties, grid search.

    Effect size:
      - pseudo_r2 : explained deviance (0-1), primary effect size
      - delta     : max-min of fitted spline in units of the (z-scored) response

    Significance:
      - pval : chi-squared approximation for the spline term (pyGAM)
               Note: can be anti-conservative with small n; interpret with caution.
    """

    x = df[x_col].dropna().values
    y = df.loc[df[x_col].notna(), y_col].values

    term = gam_s(0, n_splines=n_splines, spline_order=spline_order,
                 penalties=['derivative', 'l2'])
    gam = LinearGAM(term)

    if lam is None:
        gam.gridsearch(x[:, None], y, progress=False)
    else:
        gam.set_params(lam=lam).fit(x[:, None], y)

    x_test = np.linspace(x.min(), x.max(), n_test)
    y_pred = gam.predict(x_test)
    ci = gam.confidence_intervals(x_test, width=0.95)

    return {
        'gam'      : gam,
        'x_test'   : x_test,
        'y_pred'   : y_pred,
        'ci'       : ci,
        'pval'     : gam.statistics_['p_values'][0],
        'pseudo_r2': gam.statistics_['pseudo_r2']['explained_deviance'],
        'delta'    : float(y_pred.max() - y_pred.min()),
    }


##