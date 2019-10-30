from typing import List, Tuple
import binascii
import sys
from pyensembl import EnsemblRelease
import numpy as np
import pandas as pd

valid_chromosomes = [str(x) for x in range(1, 23)] + ['X']

def get_valid_gene_info(genes: List[str]) -> Tuple[List[str], List[int], List[int], List[int]]:
    """

    Parameters
    ----------
    genes : A list of genes
        

    Returns
    -------
    gene_names : The names of those genes
    gene_contigs : The chromosomes (contigs) associated with those genes
    gene_starts : The starting positions of those genes
    gene_ends : The ending positions of those genes

    """
    assembly = EnsemblRelease(75)
    gene_names = []
    gene_contigs = []
    gene_starts = []
    gene_ends = []
    for gene in np.intersect1d(genes, [
            gene.gene_name
            for gene in assembly.genes() if gene.contig in valid_chromosomes
    ]):  # Toss genes not in hg19
        gene_info = assembly.genes_by_name(gene)
        gene_info = gene_info[0]
        gene_names.append(gene)
        gene_contigs.append(gene_info.contig)
        gene_starts.append(gene_info.start)
        gene_ends.append(gene_info.end)
    return gene_names, gene_contigs, gene_starts, gene_ends

def get_module_score(loom, signature_name, querymask=None, nbins=100, ncontrol=5):
    """

    Parameters
    ----------
    loom : loom object on which to calculate score
        
    signature_name : Name of signature (pre-loaded into loom object) over which to calculate score
        
    nbins : Number of quantile bins to use
         (Default value = 100)
    ncontrol : Number of genes in each matched quantile
         (Default value = 5)

    Returns
    -------
    Module score

    Intended to be identical to Seurat's AddModuleScore (See Tirosh 2016)

    """
    if querymask is None:
        querymask = np.array([True]*loom.shape[1])
    alldata = loom[:,querymask]
    nonsigdata = alldata[loom.ra[signature_name]==0]
    sigdata = alldata[loom.ra[signature_name]==1]

    gene_quantiles = pd.qcut(alldata.mean(axis=1), nbins, duplicates='drop',labels=False)
    sigdata_quantiles = gene_quantiles[loom.ra[signature_name]==1]
    nonsigdata_quantiles = gene_quantiles[loom.ra[signature_name]==0]
    signature = sigdata.mean(axis=0)
    control_group = []
    for quantile in np.unique(sigdata_quantiles):
        noccurrences = np.sum(sigdata_quantiles==quantile)
        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
        control_group += list(np.random.choice(np.where(nonsigdata_quantiles==quantile)[0], size=ncontrol*noccurrences, replace=False))

    control_group = np.array(control_group)
    control = nonsigdata[control_group].mean(axis=0)
    return signature - control


def intify(df_init):
    df = df_init.copy()
    for col in df.columns:
        if col.endswith('_ad'):
            raise Exception("Don't append you column names with _ad! -- Samuel")
        df[col] = df[col].apply(lambda x:  int(binascii.hexlify(x.encode()), 16))
    while np.sum(df.max()>sys.maxsize)>0:
        for col in df.columns:
            if df[col].max() > sys.maxsize:
                df[col+'_ad'] = df[col] // sys.maxsize
                df[col] = df[col] % sys.maxsize
    return df.astype(np.int64)

def deintify(df_init):
    df = df_init.copy()
    while np.sum([x.endswith('_ad') for x in df.columns]) > 0:
        for col in df.columns:
            if col.endswith('_ad') and col+'_ad' not in df.columns:
                df[col[0:-3]] = df[col[0:-3]].astype(object)
                df[col] = df[col].astype(object)
                df[col[0:-3]] = df[col[0:-3]] + sys.maxsize * df[col]

                df.drop(col, axis=1, inplace=True)

    for col in df.columns:
        df[col] = df[col].apply(lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode())
    return df
