from typing import List, Tuple
import binascii
import sys
from pyensembl import EnsemblRelease
import numpy as np
import pandas as pd
import os
import loompy
import binascii
from tqdm import tqdm
import umap

valid_chromosomes = [str(x) for x in range(1, 23)] + ['X']


def get_valid_gene_info(
        genes: List[str]) -> Tuple[List[str], List[int], List[int], List[int]]:
    """

    Parameters
    ----------
    genes : A list of genes
        
    genes: List[str] :
        

    Returns
    -------

    
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


def get_module_score_loom(loom,
                     signature_name,
                     querymask=None,
                     nbins=100,
                     ncontrol=5):
    """

    Parameters
    ----------
    loom : loom object on which to calculate score
        
    signature_name : Name of signature (pre-loaded into loom object) over which to calculate score
        
    nbins : Number of quantile bins to use
        (Default value = 100)
    ncontrol : Number of genes in each matched quantile
        (Default value = 5)
    querymask :
         (Default value = None)

    Returns
    -------

    
    """
    if querymask is None:
        querymask = np.array([True] * loom.shape[1])
    alldata = loom[:, querymask]
    nonsigdata = alldata[loom.ra[signature_name] == 0]
    sigdata = alldata[loom.ra[signature_name] == 1]

    gene_quantiles = pd.qcut(
        alldata.mean(axis=1), nbins, duplicates='drop', labels=False)
    sigdata_quantiles = gene_quantiles[loom.ra[signature_name] == 1]
    nonsigdata_quantiles = gene_quantiles[loom.ra[signature_name] == 0]
    signature = sigdata.mean(axis=0)
    control_group = []
    for quantile in np.unique(sigdata_quantiles):
        noccurrences = np.sum(sigdata_quantiles == quantile)
        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
        control_group += list(
            np.random.choice(
                np.where(nonsigdata_quantiles == quantile)[0],
                size=ncontrol * noccurrences,
                replace=False))

    control_group = np.array(control_group)
    control = nonsigdata[control_group].mean(axis=0)
    return signature - control

def get_module_score_matrix(alldata,
                     signature_mask,
                     nbins=100,
                     ncontrol=5):
    """

    Parameters
    ----------
    alldata : matrix
        
    signature_mask : indices corresponding to signature
        
    nbins : Number of quantile bins to use
        (Default value = 100)
    ncontrol : Number of genes in each matched quantile
        (Default value = 5)

    Returns
    -------

    
    """
    assert len(signature_mask) == alldata.shape[0]
    nonsigdata = alldata[~signature_mask,:]
    sigdata = alldata[signature_mask,:]

    gene_quantiles = pd.qcut(
        alldata.mean(axis=1), nbins, duplicates='drop', labels=False)
    sigdata_quantiles = gene_quantiles[signature_mask]
    nonsigdata_quantiles = gene_quantiles[~signature_mask]
    signature = sigdata.mean(axis=0)
    control_group = []
    for quantile in np.unique(sigdata_quantiles):
        noccurrences = np.sum(sigdata_quantiles == quantile)
        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
        control_group += list(
            np.random.choice(
                np.where(nonsigdata_quantiles == quantile)[0],
                size=ncontrol * noccurrences,
                replace=False))

    control_group = np.array(control_group)
    control = nonsigdata[control_group].mean(axis=0)
    return signature - control

def intify(df_init):
    """

    Parameters
    ----------
    df_init :
        

    Returns
    -------

    """
    df = df_init.copy()
    for col in df.columns:
        if col.endswith('_ad'):
            raise Exception(
                "Don't append you column names with _ad! -- Samuel")
        df[col] = df[col].apply(lambda x: int(
            binascii.hexlify(x.encode()), 16))
    while np.sum(df.max() > sys.maxsize) > 0:
        for col in df.columns:
            if df[col].max() > sys.maxsize:
                df[col + '_ad'] = df[col] // sys.maxsize
                df[col] = df[col] % sys.maxsize
    return df.astype(np.int64)


def deintify(df_init):
    """

    Parameters
    ----------
    df_init :
        

    Returns
    -------

    """
    df = df_init.copy()
    while np.sum([x.endswith('_ad') for x in df.columns]) > 0:
        for col in df.columns:
            if col.endswith('_ad') and col + '_ad' not in df.columns:
                df[col[0:-3]] = df[col[0:-3]].astype(object)
                df[col] = df[col].astype(object)
                df[col[0:-3]] = df[col[0:-3]] + sys.maxsize * df[col]

                df.drop(col, axis=1, inplace=True)

    for col in df.columns:
        df[col] = df[col].apply(lambda x: binascii.unhexlify(
            hex(x)[2::].encode()).decode())
    return df


def seurat_to_loom(seuratrds, patient_id_column, celltype_column,
                   complexity_column, loomfile):
    """

    Parameters
    ----------
    seuratrds :
        
    patient_id_column :
        
    celltype_column :
        
    complexity_column :
        
    loomfile :
        

    Returns
    -------

    """
    import rpy2.robjects as robjects
    import rpy2.robjects as ro
    from scipy import sparse
    from rpy2.robjects import pandas2ri
    import loompy
    robjects.r('''
    library(Seurat)
    seurat2rawandmeta <- function(seuratrds) {
        seuratobj <- readRDS(seuratrds)
        return(list(genes=rownames(seuratobj@data), metadata=seuratobj@meta.data, data=as.data.frame(summary(seuratobj@data))))
    }
    ''')
    seurat_grab = robjects.r['seurat2rawandmeta'](seuratrds)
    genes = pd.DataFrame(np.array(seurat_grab.rx2('genes')))
    genes.columns = ['gene']
    metadata = pandas2ri.rpy2py_dataframe(seurat_grab.rx2('metadata'))

    if patient_id_column != 'patient_ID':
        metadata['patient_ID'] = metadata[patient_id_column]
        metadata.drop(patient_id_column, inplace=True)
    if celltype_column != 'cell_type':
        metadata['cell_type'] = metadata[celltype_column]
        metadata.drop(celltype_column, inplace=True)
    if complexity_column != 'complexity':
        metadata['complexity'] = metadata[complexity_column]
        metadata.drop(complexity_column, inplace=True)
    data_df = pandas2ri.rpy2py_dataframe(seurat_grab.rx2('data'))
    sparsedata = sparse.coo_matrix((data_df['x'], (data_df['i'] - 1,
                                                   data_df['j'] - 1))).tocsc()
    sparsedata.resize((genes.shape[0], metadata.shape[0]))

    loompy.create(loomfile, sparsedata, genes.to_dict("list"),
                  metadata.to_dict("list"))



def intify(df_init):
    """

    Parameters
    ----------
    df_init :
        

    Returns
    -------

    """
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
    """

    Parameters
    ----------
    df_init :
        

    Returns
    -------

    """
    df = df_init.copy()
    while np.sum([x.endswith('_ad') for x in df.columns]) > 0:
        for col in df.columns:
            if col.endswith('_ad') and col+'_ad' not in df.columns:
                df[col[0:-3]] = df[col[0:-3]].astype(object)
                df[col] = df[col].astype(object)
                df[col[0:-3]] = df[col[0:-3]] + sys.maxsize * df[col]

                df.drop(col, axis=1, inplace=True)

    for col in df.columns:
        try:
            df[col] = df[col].apply(lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode())
        except:
            print(df[col].apply(lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode()))
            raise Exception("whoops")
    return df


def recover_meta(db, do_deint=False):
    """

    Parameters
    ----------
    db :
        
    do_deint :
         (Default value = False)

    Returns
    -------

    """
    colmeta = None
    for key in db.ca.keys():
        if colmeta is None:
            colmeta = pd.DataFrame(db.ca[key])
            colmeta.columns = [key]
        else:
            colmeta[key] = db.ca[key]
    if do_deint:
        colmeta = deintify(colmeta.astype(np.int64))
    rowmeta = None
    for key in db.ra.keys():
        if rowmeta is None:
            rowmeta = pd.DataFrame(db.ra[key])
            rowmeta.columns = [key]
        else:
            rowmeta[key] = db.ra[key]
    if do_deint:
        rowmeta = deintify(rowmeta.astype(np.int64))    
    return rowmeta, colmeta

def generate_pca_loadings(loom, pca_type='log_tpm'):
    if pca_type == 'log_tpm':
        if [x for x in loom.ra.keys() if 'Log(TPM+1) PC' in x] == []:
            raise Exception("It seems that Log(TPM+1) PCs have not yet been calculated; ergo no loadings can be calculated")
        n_pca_cols = np.max([int(x.split(' PC ')[1]) for x in loom.ra.keys() if 'Log(TPM+1) PC' in x])
    else:
        raise Exception("Only Log(TPM+1) PC accepted currently")
#    elif pca_tpye == 'rank':
#        n_pca_cols = np.max([int(x.split(' PC ')[1]) for x in loom.ra.keys() if 'rank PC' in x])
    pcas = []
    for col in ['Log(TPM+1) PC {}'.format(x) for x in range(1,n_pca_cols+1)]:
        pcas.append(loom.ra[col])
    cellpca = np.vstack(pcas).T
    sparselogtpm = loom['log(TPM+1)'].sparse().tocsr()
    compressedlogtpm = (sparselogtpm.transpose() @ cellpca)
    for iloading in range(compressedlogtpm.shape[1]):
        loom.ca['Log(TPM+1) PC {} Loading'.format(iloading)] = compressedlogtpm[:,iloading]

def generate_embedding(loom, pca_type='log_tpm'):
    if pca_type == 'log_tpm':
        if len([x for x in loom.ra.keys() if 'Log(TPM+1) PC' in x]) == 0:
            raise Exception("It seems that Log(TPM+1) PC Loadings have not yet been calculated; ergo no UMAP embedding can be calculated")
        n_pca_cols = np.max([int(x.split(' PC ')[1].replace('Loading','')) for x in loom.ca.keys() if 'Log(TPM+1) PC ' in x and 'Loading' in x])
    else:
        raise Exception("Only Log(TPM+1) PC accepted currently")
    pca_loadings = []
    for col in ['Log(TPM+1) PC {} Loading'.format(x) for x in range(1,n_pca_cols+1)]:
        pca_loadings.append(loom.ca[col])
    compressedlogtpm = np.vstack(pca_loadings).T
    reducer = umap.UMAP(random_state=17,min_dist=0.0001,n_neighbors=30,metric='correlation',verbose=True, n_epochs=1000)
    reducer.fit(compressedlogtpm,)
    embedding = reducer.transform(compressedlogtpm)
    loom.ca['UMAP embedding 1'] = embedding[:,0]
    loom.ca['UMAP embedding 2'] = embedding[:,1]
#    plt.scatter(embedding[:,0], embedding[:,1])
#    plt.show()
