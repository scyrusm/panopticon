from pyensembl import EnsemblRelease
import numpy as np
import pandas as pd

valid_chromosomes = [str(x) for x in range(1, 23)] + ['X']

def get_valid_gene_info(genes):
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
