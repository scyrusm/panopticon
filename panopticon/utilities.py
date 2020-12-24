import sys
from typing import List, Tuple

import numpy as np
import pandas as pd


def get_valid_gene_info(
        genes: List[str]) -> Tuple[List[str], List[int], List[int], List[int]]:
    """Returns gene locations for all genes in ensemble 93  --S Markson 3 June 2020

    Parameters
    ----------
    genes : A list of genes
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes: List[str] :
        

    Returns
    -------

    
    """
    from pyensembl import EnsemblRelease
    assembly = EnsemblRelease(93)
    valid_chromosomes = [str(x) for x in range(1, 23)] + ['X']
    gene_names = []
    gene_contigs = []
    gene_starts = []
    gene_ends = []
    for gene in np.intersect1d(genes, [
            gene.gene_name
            for gene in assembly.genes() if gene.contig in valid_chromosomes
    ]):  # Toss genes not in hg38 release 93
        gene_info = assembly.genes_by_name(gene)
        gene_info = gene_info[0]
        gene_names.append(gene)
        gene_contigs.append(gene_info.contig)
        gene_starts.append(gene_info.start)
        gene_ends.append(gene_info.end)
    return gene_names, gene_contigs, gene_starts, gene_ends


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
    sparsedata = sparse.coo_matrix(
        (data_df['x'], (data_df['i'] - 1, data_df['j'] - 1))).tocsc()
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
    import binascii
    df = df_init.copy()
    for col in df.columns:
        if col.endswith('_ad'):
            raise Exception(
                "Don't append you column names with _ad! -- Samuel")
        df[col] = df[col].apply(
            lambda x: int(binascii.hexlify(x.encode()), 16))
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
    import binascii
    df = df_init.copy()
    while np.sum([x.endswith('_ad') for x in df.columns]) > 0:
        for col in df.columns:
            if col.endswith('_ad') and col + '_ad' not in df.columns:
                df[col[0:-3]] = df[col[0:-3]].astype(object)
                df[col] = df[col].astype(object)
                df[col[0:-3]] = df[col[0:-3]] + sys.maxsize * df[col]

                df.drop(col, axis=1, inplace=True)

    for col in df.columns:
        try:
            df[col] = df[col].apply(
                lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode())
        except:
            print(df[col].apply(
                lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode()))
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


def we_can_pickle_it(thing, thingname: str):
    """

    Parameters
    ----------
    thing :
        
    thingname : str :
        
    thingname : str :
        
    thingname : str :
        
    thingname: str :
        

    Returns
    -------

    
    """
    import pickle
    with open(thingname, 'wb') as f:
        pickle.dump(thing, f, pickle.HIGHEST_PROTOCOL)


def we_can_unpickle_it(thingname: str):
    """

    Parameters
    ----------
    thingname : str :
        
    thingname : str :
        
    thingname : str :
        
    thingname: str :
        

    Returns
    -------

    
    """
    import pickle
    with open(thingname, 'rb') as f:
        thing = pickle.load(f)
    return thing
