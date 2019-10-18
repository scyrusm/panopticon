import pandas as pd
from pyensembl import EnsemblRelease
from panopticon.utilities import get_valid_gene_info
from tqdm import tqdm
import numpy as np
from scipy import stats

def segmentation_to_copy_ratio_dict(genes,
                                             segmentation,
                                             chrom_col='chrom',
                                             start_col='chromStart',
                                             end_col='chromEnd',
                                             score_col='copyRatio',
                                             log2=False):
    """


    Parameters
    ----------
    segmentation : File containing the segmentation.  Should be readable with pandas.read_table.
        
    chrom_col : Column name indicating chromosomes (bed file column name is default)
        (Default value = 'chrom')
    start_col : Column name indicating region start (bed file column name is default)
        (Default value = 'chromStart')
    end_col : Column name indicating region end (bed file column name is default)
        (Default value = 'chromEnd')
    score_col : Column name indicating the copy number of that region (or log_2 of the copy number, if log2 == True
        (Default value = 'score')
    log2 :  Indicates whether the score_col column indicates the copy ratio, or the base-2 log of the copy ratio
         (Default value = False)

    Returns
    -------
    A dictionary of genes to their copy ratios
    
    """

    gene_names, gene_contigs, gene_starts, gene_ends = get_valid_gene_info(
        genes)
    segmentation_cp = segmentation.copy(
    )  # So intermediate calculations don't get saved
    if log2:
        segmentation_cp[score_col] = 2**segmentation_cp[score_col]

    chromosome_scale = 10e12
    gene_df = pd.DataFrame(gene_names)
    gene_df.columns = ['gene']
    gene_df[chrom_col] = gene_contigs
    gene_df[start_col] = gene_starts
    gene_df[end_col] = gene_ends
    gene_df['pseudoposition'] = gene_df[chrom_col].apply(
        lambda x: 23
        if x == 'X' else int(x)) * chromosome_scale + gene_df[start_col]

    segmentation_cp['pseudoposition'] = segmentation_cp[chrom_col].apply(
        lambda x: 23
        if x == 'X' else int(x)) * chromosome_scale + segmentation[start_col]
    return pd.merge_asof(
        gene_df.sort_values('pseudoposition'),
        segmentation_cp[[
            chrom_col, start_col, end_col, score_col, 'pseudoposition'
        ]],
        on='pseudoposition',
        suffixes=('_gene_df', '_segmentation'),
        direction='backward').set_index('gene')[score_col].to_dict()


def multireference_dna_correspondence(loom,
                                           loomquery, 
                                           *segmentations):

    list_of_correlations = []
    for segmentation in tqdm(segmentations, desc='Processing segmentations'):
        crd = loom.ra[segmentation]
        expressions = loom[:, loomquery]

        genemask = ~np.isnan(crd)#np.isin(loom.ra['gene'], list(crd.keys()))
        expressions = expressions[genemask, :]
        validgenes = loom.ra['gene'][genemask]
#        tcrs = np.array([crd[gene] for gene in validgenes])  #[low_dropout]
        tcrs = crd[genemask]

        correlations = []
        for icell in range(expressions.shape[1]):
            correlations.append(
                stats.kendalltau(tcrs, expressions[:, icell]).correlation)
        list_of_correlations.append(correlations)
    return list_of_correlations

