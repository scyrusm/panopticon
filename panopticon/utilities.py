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
    """Returns gene locations for all genes in ensemble 75--this is outdated!  --S Markson 3 June 2020

    Parameters
    ----------
    genes : A list of genes
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
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
    """Calculates a module score over a loom file.  This routine is deprecated--use generate masked module score (S Markson 3 June 2020).

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


def get_module_score_matrix(alldata, signature_mask, nbins=100, ncontrol=5):
    """generates a module score (a la Seurat's AddModuleScore, see Tirosh 2016) on a matrix, with a mask.  I don't call this directly (S Markson 3 June 2020).

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
    nonsigdata = alldata[~signature_mask, :]
    sigdata = alldata[signature_mask, :]

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
        try:
            df[col] = df[col].apply(lambda x: binascii.unhexlify(
                hex(x)[2::].encode()).decode())
        except:
            print(df[col].apply(lambda x: binascii.unhexlify(
                hex(x)[2::].encode()).decode()))
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


def generate_incremental_pca(loom, layername, batch_size=512, n_components=50):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    batch_size :
        (Default value = 512)
    n_components :
        (Default value = 50)

    Returns
    -------

    
    """

    from tqdm import tqdm
    from sklearn.decomposition import IncrementalPCA
    pca = IncrementalPCA(n_components=n_components)
    for (ix, selection, view) in tqdm(
            loom.scan(axis=1, batch_size=batch_size),
            total=loom.shape[1] // batch_size):
        #pca.partial_fit(view[:, :].transpose())
        pca.partial_fit(view[layername][:, :].T)
    for i in range(50):
        loom.ra['{} PC {}'.format(layername, i + 1)] = pca.components_[i]
    loom.attrs['NumberPrincipalComponents_{}'.format(layername)] = n_components
    loom.attrs['PCAExplainedVarianceRatio_{}'.format(
        layername)] = pca.explained_variance_ratio_


def generate_pca_loadings(loom, layername, dosparse=False, batch_size=512):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    dosparse :
        (Default value = False)
    batch_size :
        (Default value = 512)

    Returns
    -------

    
    """
    if len([x for x in loom.ra.keys() if '{} PC'.format(layername) in x]) == 0:
        raise Exception(
            "It seems that {} PCs have not yet been calculated; ergo no UMAP embedding can be calculated"
            .format(layername))


#    n_pca_cols = np.max([
#        int(x.split(' PC ')[1]) for x in loom.ra.keys()
#        if '{} PC '.format(layername) in x
#    ])
    n_pca_cols = loom.attrs['NumberPrincipalComponents_{}'.format(layername)]

    #    elif pca_tpye == 'rank':
    #        n_pca_cols = np.max([int(x.split(' PC ')[1]) for x in loom.ra.keys() if 'rank PC' in x])
    pcas = []
    for col in [
            '{} PC {}'.format(layername, x) for x in range(1, n_pca_cols + 1)
    ]:
        pcas.append(loom.ra[col])
    cellpca = np.vstack(pcas).T
    if dosparse:
        sparsedata = loom[layername].sparse().tocsr()
        compresseddata = (sparsedata.transpose() @ cellpca)
    else:
        compresseddatalist = []
        for (ix, selection, view) in tqdm(
                loom.scan(axis=1, batch_size=batch_size),
                total=loom.shape[1] // batch_size):
            compresseddatalist.append(view[layername][:, :].T @ cellpca)
    compresseddata = np.vstack(compresseddatalist)
    for iloading in range(compresseddata.shape[1]):
        loom.ca['{} PC {} Loading'.format(
            layername, iloading + 1)] = compresseddata[:, iloading]


def generate_embedding(loom,
                       layername,
                       min_dist=0.0001,
                       n_neighbors=30,
                       n_epochs=1000,
                       metric='correlation',
                       random_state=None,
                       components_to_use=None,
                       mode='pca'):
    """

    Parameters
    ----------
    loom :
        
    pca_type :
        (Default value = 'log_tpm')
    layername :
        
    min_dist :
        (Default value = 0.0001)
    n_neighbors :
        (Default value = 30)
    n_epochs :
        (Default value = 1000)
    metric :
        (Default value = 'correlation')
    random_state :
        (Default value = None)
    pca_cols_to_use :
        (Default value = None)
    components_to_use :
        (Default value = None)
    mode :
        (Default value = 'nmf')

    Returns
    -------

    
    """
    if mode not in ['pca', 'nmf']:
        raise Exception("Currently only two modes implemented:  nmf and pca")
    if mode == 'pca':
        n_pca_cols = loom.attrs['NumberPrincipalComponents_{}'.format(
            layername)]
        pca_loadings = []
        if components_to_use != None:
            for col in [
                    '{} PC {} Loading'.format(layername, x)
                    for x in components_to_use
            ]:
                pca_loadings.append(loom.ca[col])
        else:
            for col in [
                    '{} PC {} Loading'.format(layername, x)
                    for x in range(1, n_pca_cols + 1)
            ]:
                pca_loadings.append(loom.ca[col])
        compressed = np.vstack(pca_loadings).T
    elif mode == 'nmf':
        n_nmf_cols = loom.attrs['NumberNMFComponents']
        print(n_nmf_cols)
        nmf_loadings = []
        if components_to_use != None:
            for col in [
                    '{} NMF Loading Component {}'.format(layername, x)
                    for x in components_to_use
            ]:
                nmf_loadings.append(loom.ca[col])
        else:
            for col in [
                    '{} NMF Loading Component {}'.format(layername, x)
                    for x in range(1, n_nmf_cols + 1)
            ]:
                nmf_loadings.append(loom.ca[col])
        print(len(nmf_loadings))
        compressed = np.vstack(nmf_loadings).T
    reducer = umap.UMAP(
        random_state=None,
        min_dist=min_dist,
        n_neighbors=n_neighbors,
        metric=metric,
        verbose=True,
        n_epochs=n_epochs)
    embedding = reducer.fit_transform(compressed)
    loom.ca['{} UMAP embedding 1'.format(mode.upper())] = embedding[:, 0]
    loom.ca['{} UMAP embedding 2'.format(mode.upper())] = embedding[:, 1]


def get_subclustering(X,
                      score_threshold,
                      max_clusters=50,
                      min_input_size=10,
                      silhouette_threshold=0.1):
    """

    Parameters
    ----------
    embedding :
        
    score_threshold :
        
    max_clusters :
        (Default value = 10)
    X :
        
    min_input_size :
        (Default value = 5)

    Returns
    -------

    
    """
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import AgglomerativeClustering
    from tqdm import tqdm

    if X.shape[0] < min_input_size:
        return np.array([0] * X.shape[0])
    else:
        clustering = AgglomerativeClustering(
            n_clusters=2,
            memory='clusteringcachedir/',
            affinity='cosine',
            compute_full_tree=True,
            linkage='average')
        scores = []
        minnk = 2
        for nk in tqdm(range(minnk, np.min([max_clusters, X.shape[0]]), 1)):
            clustering.set_params(n_clusters=nk)
            clustering.fit(X)

            score = silhouette_score(
                X,
                clustering.labels_,
                metric='cosine',
                sample_size=np.min([5000, X.shape[0]]))
            scores.append(score)
            #break
        print(np.array(scores))
        if np.max(scores) >= score_threshold:
            print("Number of clusters:", np.argmax(scores) + minnk)
            clustering.set_params(n_clusters=np.argmax(scores) + minnk)
            clustering.fit(X)
            return clustering.labels_
        else:
            return np.array([0] * X.shape[0])


def generate_clustering(loom,
                        layername,
                        clustering_depth=3,
                        starting_clustering_depth=0,
                        max_clusters='sqrt_rule',
                        mode='pca',
                        silhouette_threshold=0.1):
    """

    Parameters
    ----------
    loom :
        
    clustering_depth :
        (Default value = 3)
    starting_clustering_depth :
        (Default value = 0)
    max_clusters :
        (Default value = 200)
    layername :
        
    mode :
        (Default value = 'pca')

    Returns
    -------

    
    """
    if type(clustering_depth) != int or clustering_depth < 1 or type(
            starting_clustering_depth) != int:
        raise Exception(
            "clustering_depth and starting_clustering_depth must be natural numbers."
        )
    if (starting_clustering_depth > 0) and ('ClusteringIteration{}'.format(
            starting_clustering_depth - 1) not in loom.ca.keys()):
        raise Exception(
            "starting_clustering_depth not yet computed; please run with lower starting_clustering depth, or 0"
        )
    if mode not in ['pca', 'nmf']:
        raise Exception("Currently only implemented for modes:  pca and nmf")
    from time import time
    from sklearn.decomposition import IncrementalPCA
    from sklearn.preprocessing import StandardScaler
    from tqdm import tqdm
    from panopticon.utilities import get_subclustering
    if mode == 'pca':
        from sklearn.decomposition import PCA
    elif mode == 'nmf':
        from sklearn.decomposition import NMF

    if starting_clustering_depth == 0:
        if mode == 'nmf':
            n_nmf_cols = loom.attrs['NumberNMFComponents']
            nmf_loadings = []
            for col in [
                    '{} NMF Loading Component {}'.format(layername, x)
                    for x in range(1, n_nmf_cols + 1)
            ]:
                nmf_loadings.append(loom.ca[col])
            X = np.vstack(nmf_loadings).T
        elif mode == 'pca':
            #            if layername.endswith('_standardized'):
            #                raise Exception(
            #                    "Standardization automatically performed, use non-standardized layer (e.g. log(TPM+1)"
            #                )
            #            if layername + '_standardized' not in loom.layers.keys():
            #                raise Exception(
            #                    "Must pre-compute PCA on standardized version of desired layer"
            #                )
            #            n_pca_cols = loom.attrs[
            #                'NumberPrincipalComponents_{}_standardized'.format(layername)]
            n_pca_cols = loom.attrs['NumberPrincipalComponents_{}'.format(
                layername)]
            pca_loadings = []
            #            for col in [
            #                    '{} PC {} Loading'.format(layername + '_standardized', x)
            #                    for x in range(1, n_pca_cols + 1)
            #            ]:
            for col in [
                    '{} PC {} Loading'.format(layername, x)
                    for x in range(1, n_pca_cols + 1)
            ]:
                pca_loadings.append(loom.ca[col])
            X = np.vstack(pca_loadings).T
        if max_clusters == 'sqrt_rule':
            clustering = get_subclustering(
                X,
                silhouette_threshold,
                max_clusters=int(np.floor(np.sqrt(X.shape[0])))
            )  # This shouldn't be hard-coded S Markson 9 June 2020
        else:
            clustering = get_subclustering(
                X, silhouette_threshold, max_clusters=max_clusters
            )  # This shouldn't be hard-coded S Markson 9 June 2020

        loom.ca['ClusteringIteration0'] = clustering
        starting_clustering_depth = 1

    for subi in range(starting_clustering_depth, clustering_depth):

        loom.ca['ClusteringIteration{}'.format(subi)] = ['U'] * len(
            loom.ca['ClusteringIteration{}'.format(subi - 1)])

        for cluster in set([
                x for x in loom.ca['ClusteringIteration{}'.format(subi - 1)]
                if x != 'U'
        ]):  #will need to fix
            mask = loom.ca['ClusteringIteration{}'.format(
                subi -
                1)] == cluster  #first mask, check for top level clustering
            #break
            start = time()
            data_c = loom[layername][:, mask]
            print("processing cluster", cluster, "; time to load: ",
                  time() - start, ", mask size: ", np.sum(mask))
            if mode == 'nmf':
                model = NMF(
                    n_components=np.min([50, data_c.shape[1]]),
                    init='random',
                    random_state=0)
                X = model.fit_transform(data_c.T)
            elif mode == 'pca':

                #                data_c = StandardScaler().fit_transform(
                #                    data_c.T
                #                )  # This transposition is confusing.  S Markson 7 June 2020
                data_c = data_c.T
                if data_c.shape[0] > 5000:
                    model = IncrementalPCA(n_components=10)
                    for chunk in tqdm(
                            np.array_split(
                                data_c, data_c.shape[0] // 512, axis=0),
                            desc='partial fitting over chunks of masked data'):
                        model.partial_fit(chunk)
                    X = model.transform(data_c)
                    print("EV", model.explained_variance_)
                    print("EVR", model.explained_variance_ratio_)
                else:
                    model = PCA(
                        n_components=np.min([10, data_c.shape[0]]),
                        random_state=0)

                    X = model.fit_transform(data_c)
                    print("EV", model.explained_variance_)
                    print("EVR", model.explained_variance_ratio_)

            if max_clusters == 'sqrt_rule':
                print("xshape", X.shape)
                nopath_clustering = get_subclustering(
                    X,
                    silhouette_threshold,
                    max_clusters=int(np.floor(np.sqrt(X.shape[0])))
                )  # This shouldn't be hard-coded S Markson 9 June 2020
            else:
                nopath_clustering = get_subclustering(
                    X, silhouette_threshold, max_clusters=max_clusters
                )  # This shouldn't be hard-coded S Markson 9 June 2020


#            nopath_clustering = get_subclustering(X, score_threshold=silhouette_threshold)  #Really shouldn't be hard-coded S Markson 9 June 2020
            fullpath_clustering = [
                '{}-{}'.format(cluster, x) for x in nopath_clustering
            ]
            loom.ca['ClusteringIteration{}'.format(
                subi)][mask] = fullpath_clustering
        loom.ca['ClusteringIteration{}'.format(subi)] = loom.ca[
            'ClusteringIteration{}'.format(
                subi)]  #This is to force the changes to save to disk

import numpy as np
from tqdm import tqdm


def cluster_differential_expression(loom,
                                    cluster_level,
                                    layername,
                                    ident1,
                                    ident2=None,
                                    mask1=None,
                                    mask2=None,
                                    verbose=False,
                                    ident1_downsample_size=None,
                                    ident2_downsample_size=None):
    """

    Parameters
    ----------
    loom :
        
    cluster_level :
        
    layername :
        
    ident1 :
        
    ident2 :
        (Default value = None)
    verbose :
        (Default value = False)
    ident1_downsample_size :
        (Default value = None)
    ident2_downsample_size :
        (Default value = None)
    mask1 :
        (Default value = None)
    mask2 :
        (Default value = None)

    Returns
    -------

    
    """
    from scipy.stats import mannwhitneyu
    from time import time
    import pandas as pd

    if (mask1 is not None) and (mask2 is not None):
        print("ignoring ident1, ident2")

    elif (mask1 is not None) or (mask2 is not None):
        raise Exception(
            "Either both or neither of mask1, mask2 must be specified")
    else:
        if type(ident1) != list:
            ident1 = [ident1]
        mask1 = np.isin(loom.ca[cluster_level], ident1)
        if ident2 == None:
            print(
                "Automatic complement: Cells in same subcluster except lowest subcluster"
            )
            if cluster_level == 'ClusteringIteration0':
                mask2 = ~np.isin(loom.ca[cluster_level], ident1)
                clusterset = np.unique(loom.ca[cluster_level])
                ident2 = list(np.setdiff1d(clusterset, ident1))
            else:
                cluster_level_number = int(cluster_level[-1])
                prefices = [
                    '-'.join(x.split('-')[0:(cluster_level_number)]) +
                    '-'  # 13 Apr 2020--check that this works
                    for x in ident1
                ]
                if len(np.unique(prefices)) > 1:
                    raise Exception(
                        "Cluster Differential expression with automatic complement must use ident1 only from cells from same n-1th subcluster"
                    )
                prefix = prefices[0]
                clusterset = [
                    x for x in np.unique(loom.ca[cluster_level])
                    if x.startswith(prefix)
                ]
                ident2 = list(np.setdiff1d(clusterset, ident1))
                mask2 = np.isin(loom.ca[cluster_level], ident2)

        else:
            if type(ident2) != list:
                ident2 = [ident2]
            mask2 = np.isin(loom.ca[cluster_level], ident2)
        print("Comparison of", ident1, "against", ident2)
    if (np.sum(mask1) == 0) or (np.sum(mask2) == 0):
        return np.nan
    if ident1_downsample_size:
        p = np.min([ident1_downsample_size, np.sum(mask1)]) / np.sum(mask1)
        mask1 *= np.random.choice([True, False],
                                  p=[p, 1 - p],
                                  size=mask1.shape[0])
    if ident2_downsample_size:
        p = np.min([ident2_downsample_size, np.sum(mask2)]) / np.sum(mask2)
        mask2 *= np.random.choice([True, False],
                                  p=[p, 1 - p],
                                  size=mask2.shape[0])
    print(np.sum(mask1), np.sum(mask2))
    pvalues = []
    genes = []
    meanexpr1 = []
    meanexpr2 = []
    start = time()
    data1 = loom[layername][:, mask1]
    if verbose:
        print('First matrix extracted in', time() - start, 'seconds')
    start = time()
    data2 = loom[layername][:, mask2]
    if verbose:
        print('Second matrix extracted', time() - start, 'seconds')
    for igene, gene in enumerate(
            tqdm(loom.ra['gene'], desc='Computing Mann-Whitney p-values')):
        genes.append(gene)
        if np.std(data1[igene, :]) + np.std(data2[igene, :]) < 1e-14:
            pvalues.append(1)
        else:
            pvalues.append(
                mannwhitneyu(data1[igene, :], data2[igene, :]).pvalue)
        meanexpr1.append(data1[igene, :].mean())
        meanexpr2.append(data2[igene, :].mean())
    output = pd.DataFrame(genes)
    output.columns = ['gene']
    output['pvalue'] = pvalues
    output['MeanExpr1'] = meanexpr1
    output['MeanExpr2'] = meanexpr2
    return output.sort_values('pvalue', ascending=True)


def get_cluster_markers(loom, layername, cluster_level):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    cluster_level :
        

    Returns
    -------

    
    """
    from panopticon.utilities import cluster_differential_expression
    diffex = {}
    for cluster in np.unique(loom.ca[cluster_level]):
        diffex[cluster] = cluster_differential_expression(
            loom,
            cluster_level,
            layername,
            ident1=cluster,
            verbose=True,
            ident1_downsample_size=500,
            ident2_downsample_size=500).query('MeanExpr1 > MeanExpr2')
    return diffex


def we_can_pickle_it(thing, thingname: str):
    """

    Parameters
    ----------
    thing :
        
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
        
    thingname: str :
        

    Returns
    -------

    
    """
    import pickle
    with open('diffex_level1.pkl', 'rb') as f:
        thing = pickle.load(f)
    return thing
    with open(thingname, 'wb') as f:
        pickle.dump(thing, f, pickle.HIGHEST_PROTOCOL)


def get_gsea_with_selenium(diffi):
    """If you aren't Sam, probably don't use this.

    Parameters
    ----------
    diffi :
        

    Returns
    -------

    
    """
    from selenium import webdriver
    from selenium.webdriver.common.keys import Keys
    from selenium.webdriver.support.ui import Select
    import os
    import time
    import sys
    import pandas as pd

    driver = webdriver.Firefox()
    genelisthitdict = {}
    flag = True
    for diffex in diffi:
        for key in diffex.keys():
            driver.get("https://www.gsea-msigdb.org/gsea/msigdb/annotate.jsp")
            if flag:
                elem = driver.find_element_by_name("j_username")
                elem.send_keys("smarkson@broadinstitute.org")
                elem = driver.find_element_by_name('login')
                elem.click()
                flag = False
            elem = driver.find_element_by_id("geneList")
            elem.send_keys('\n'.join(diffex[key]['gene'][0:10].values))
            for genelistcheckbox in ["H", "C2", "C5", "C6", "C7"]:
                elem = driver.find_element_by_id(
                    "{}checkbox".format(genelistcheckbox))
                elem.click()
            #driver.close()
            elem = driver.find_element_by_name("viewOverlaps")
            elem.click()
            try:
                elem = driver.find_element_by_link_text("Text")
                elem.click()
                #driver.close()
                footersize = 0
                while True:
                    try:
                        genelisthits = pd.read_table(
                            "/home/smarkson/Downloads/overlap.tsv",
                            skiprows=10,
                            skipfooter=footersize,
                            header=None)
                        break
                    except:
                        footersize += 1
                genelisthits = pd.read_table(
                    "/home/smarkson/Downloads/overlap.tsv",
                    skiprows=10,
                    skipfooter=footersize + 2,
                    header=None)
                genelisthitdict[key] = genelisthits
                os.system("rm /home/smarkson/Downloads/overlap.tsv")
            except:
                genelisthitdict[key] = 'no hits'
            #print(genelisthitdict[key])
    return genelisthitdict


def get_cluster_embedding(loom,
                          layername,
                          cluster,
                          min_dist=0.01,
                          n_neighbors=None,
                          verbose=False,
                          mask=None):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    cluster :
        
    min_dist :
        (Default value = 0.01)
    n_neighbors :
        (Default value = None)
    verbose :
        (Default value = False)
    mask :
        (Default value = None)

    Returns
    -------

    
    """
    from sklearn.decomposition import PCA
    import umap
    import numpy as np
    clustering_level = len(str(cluster).split('-')) - 1
    if mask is None:
        mask = loom.ca['ClusteringIteration{}'.format(
            clustering_level)] == cluster
    else:
        print("Clustering over custom mask--ignoring cluster argument")
    if n_neighbors is None:
        n_neighbors = int(np.sqrt(np.sum(mask)))
    data = loom[layername][:, mask]

    pca = PCA(n_components=np.min([50, data.shape[1]]))
    pca.fit(data[:, :].transpose())

    cellpca = (data.T @ pca.components_.T)
    reducer = umap.UMAP(
        random_state=17,
        verbose=verbose,
        min_dist=min_dist,
        n_neighbors=n_neighbors,
        metric='correlation')
    reducer.fit(cellpca, )
    embedding = reducer.transform(cellpca)
    return embedding


def plot_subclusters(loom,
                     layername,
                     cluster,
                     sublayers=1,
                     plot_output=None,
                     label_clusters=True,
                     complexity_cutoff=0):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    cluster :
        
    sublayers :
        (Default value = 1)
    plot_output :
        (Default value = None)
    label_clusters :
        (Default value = True)
    complexity_cutoff :
        (Default value = 0)

    Returns
    -------

    
    """
    from panopticon.utilities import get_cluster_embedding
    import matplotlib.pyplot as plt
    current_layer = len(str(cluster).split('-')) - 1
    if complexity_cutoff > 0:
        print("Using complexity cutoff")
        supermask = (loom.ca['ClusteringIteration{}'.format(current_layer)] ==
                     cluster) & (loom.ca['nGene'] >= complexity_cutoff
                                 )  # I kinda hate this (30 Apr 2020 smarkson)
        embedding = get_cluster_embedding(
            loom, layername, cluster, mask=supermask)
    else:
        embedding = get_cluster_embedding(loom, layername, cluster)
        supermask = loom.ca['ClusteringIteration{}'.format(
            current_layer)] == cluster
    subclusters = [
        x for x in np.unique(loom.ca['ClusteringIteration{}'.format(
            current_layer + sublayers)])
        if str(x).startswith(str(cluster) + '-')
    ]
    for subcluster in subclusters:
        mask = loom.ca['ClusteringIteration{}'.format(
            current_layer + sublayers)][supermask] == subcluster
        plt.scatter(embedding[mask, 0], embedding[mask, 1], label=subcluster)
        if label_clusters:
            plt.annotate(
                subcluster,
                (np.mean(embedding[mask, 0]), np.mean(embedding[mask, 1])))
    plt.legend()
    if plot_output is not None:
        plt.savefig(plot_output)
    else:
        plt.show()


def get_metafield_breakdown(loom,
                            cluster,
                            field,
                            complexity_cutoff=0,
                            mask=None):
    """

    Parameters
    ----------
    loom :
        
    cluster :
        
    field :
        
    complexity_cutoff :
        (Default value = 0)
    mask :
        (Default value = None)

    Returns
    -------

    
    """
    cluster_level = len(str(cluster).split('-')) - 1
    if mask is None:
        mask = (loom.ca['ClusteringIteration{}'.format(cluster_level)] ==
                cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    else:
        print("ignoring cluster, using custom mask")
        mask = (mask) & (loom.ca['nGene'] >= complexity_cutoff)
    return pd.DataFrame(loom.ca[field][mask])[0].value_counts()


def generate_masked_module_score(loom, layername, mask, genelist, ca_name):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    mask :
        
    genelist :
        
    ca_name :
        

    Returns
    -------

    
    """
    from panopticon.utilities import get_module_score_matrix
    matrix = loom[layername][:, mask]
    sigmask = np.isin(loom.ra['gene'], genelist)
    sig_score = get_module_score_matrix(matrix, sigmask)
    maskedscores = []
    counter = 0
    for flag in mask:
        if flag:
            maskedscores.append(sig_score[counter])
            counter += 1
        else:
            maskedscores.append(np.nan)
    loom.ca[ca_name] = maskedscores


def generate_gene_variances(loom, layername):
    """

    Parameters
    ----------
    loom :
        
    layername :
        

    Returns
    -------

    
    """
    loom.ra['GeneVar'] = loom[layername].map([np.var])[0]


def generate_nmf_and_loadings(loom,
                              layername,
                              nvargenes=2000,
                              n_components=100,
                              verbose=False):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    nvargenes :
        (Default value = 2000)
    n_components :
        (Default value = 100)
    verbose :
        (Default value = False)

    Returns
    -------

    
    """
    from sklearn.decomposition import NMF

    if 'GeneVar' not in loom.ra.keys():
        raise Exception(
            "Necessary to have already generated gene expression variances")
    vargenemask = loom.ra['GeneVar'] > np.sort(
        loom.ra['GeneVar'])[::-1][nvargenes]
    X = loom[layername][vargenemask, :]
    model = NMF(
        n_components=n_components,
        init='random',
        random_state=0,
        verbose=verbose)
    W = model.fit_transform(X)
    H = model.components_

    # record NMF basis
    nmf_factors = []
    counter = 0
    for isvargene in vargenemask:
        if isvargene:
            nmf_factors.append(W[counter, :])
            counter += 1
        else:
            nmf_factors.append(np.zeros(W.shape[1]))
    nmf_factors = np.vstack(nmf_factors)
    factor_sums = []
    for i in range(nmf_factors.shape[1]):
        loom.ra['{} NMF Component {}'.format(
            layername, i + 1)] = nmf_factors[:, i] / np.sum(nmf_factors[:, i])
        factor_sums.append(np.sum(nmf_factors[:, i]))
    factor_sums = np.array(factor_sums)
    # record NMF loadings
    for i in range(H.shape[0]):
        loom.ca['{} NMF Loading Component {}'.format(
            layername, i + 1)] = H[i, :] * factor_sums[i]

    loom.attrs['NumberNMFComponents'] = n_components


def generate_cell_and_gene_quality_metrics(loom, layername):
    """

    Parameters
    ----------
    loom :
        
    layername :
        

    Returns
    -------

    
    """
    import numpy as np
    from statsmodels import robust

    rpgenemask = np.array([x.startswith('RP') for x in loom.ra['gene']])
    madrp, meanrp, maxrp = loom[layername].map([robust.mad, np.mean, np.max],
                                               axis=1,
                                               selection=rpgenemask)
    madall, meanall = loom[layername].map([robust.mad, np.mean], axis=1)
    loom.ca['RibosomalRelativeMeanAbsoluteDeviation'] = madrp / meanrp
    loom.ca['RibosomalMaxOverMean'] = maxrp / meanrp
    loom.ca['AllGeneRelativeMeanAbsoluteDeviation'] = madall / meanall

    def complexity(vec):
        """

        Parameters
        ----------
        vec :
            

        Returns
        -------

        
        """
        return np.sum(vec > 0)

    loom.ca['nGene'] = loom[layername].map([complexity], axis=1)[0]
    loom.ra['nCell'] = loom[layername].map([complexity], axis=0)[0]


def get_cluster_specific_greater_than_cutoff_mask(loom,
                                                  metric,
                                                  cluster_level,
                                                  default_cutoff,
                                                  exception_dict={}):
    """

    Parameters
    ----------
    loom :
        
    metric :
        
    cluster_level :
        
    default_cutoff :
        
    exception_dict :
        (Default value = {})

    Returns
    -------

    
    """
    mask = []
    for metric_val, cluster in zip(loom.ca[metric], loom.ca[cluster_level]):

        if cluster in exception_dict.keys():
            mask.append(metric_val > exception_dict[cluster])
        else:
            mask.append(metric_val > default_cutoff)
    return np.array(mask)


def create_subsetted_loom(loom, output_loom, cellmask, genemask):
    """

    Parameters
    ----------
    loom :
        
    output_loom :
        
    cellmask :
        
    genemask :
        

    Returns
    -------

    
    """
    from panopticon.utilities import recover_meta
    if '' not in loom.layers.keys():
        raise Exception("Expecting '' layer, yet none found")
    rowmeta, colmeta = recover_meta(loom)
    loompy.create(
        output_loom, loom[''].sparse().tocsr()[:, cellmask][genemask, :],
        rowmeta[genemask].to_dict("list"), colmeta[cellmask].to_dict("list"))
    with loompy.connect(output_loom) as smallerloom:
        for layer in [x for x in loom.layer.keys() if x != '']:
            smallerloom[layer] = loom[layer][:, cellmask][genemask, :]


def get_patient_averaged_table(loom,
                               patient_key='patient_ID',
                               column_attributes=[],
                               n_cell_cutoff=0):
    """

    Parameters
    ----------
    loom :
        
    patient_key :
        (Default value = 'patient_ID')
    column_attributes :
        (Default value = [])
    n_cell_cutoff :
        (Default value = 0)

    Returns
    -------

    
    """
    unfiltered = pd.DataFrame(loom.ca[patient_key])
    unfiltered.columns = [patient_key]
    for column_attribute in column_attributes:
        unfiltered[column_attribute] = loom.ca[column_attribute]
    unfiltered.groupby(patient_key).mean()
    threshold_filter = lambda x: np.sum(np.isnan(x)) > n_cell_cutoff
    filtered = (
        unfiltered.groupby(patient_key).apply(threshold_filter)).replace(
            to_replace=False,
            value=np.nan) * unfiltered.groupby(patient_key).mean()
    return filtered


def generate_standardized_layer(loom, layername, variance_axis='cell'):
    """

    Parameters
    ----------
    loom :
        
    layername :
        

    Returns
    -------

    
    """
    if layername.endswith('_standardized'):
        raise Exception(
            "It appears that layer is already standardized; note that _standardized suffix is reserved"
        )
    if variance_axis not in ['gene', 'cell']:
        raise Exception(
            "Can only set variance to be 1 along cell axis or gene axis")
    from sklearn.preprocessing import StandardScaler
    if variance_axis == 'gene':
        loom[layername +
             '_gene_standardized'] = StandardScaler().fit_transform(
                 loom[layername][:, :].T).T
    elif variance_axis == 'cell':
        loom[layername +
             '_cell_standardized'] = StandardScaler().fit_transform(
                 loom[layername][:, :])


def generate_count_normalization(loom,
                                 raw_count_layername,
                                 output_layername,
                                 denominator=10**5):
    """

    Parameters
    ----------
    loom :
        
    raw_count_layername :
        
    output_layername :
        
    denominator :
         (Default value = 10**5)

    Returns
    -------

    """
    import numpy as np
    colsums = loom[raw_count_layername].map([np.sum], axis=1)[0]
    sparselayer = loom[raw_count_layername].sparse()
    sparselayer = sparselayer.multiply(denominator / colsums).log1p().multiply(
        1 / np.log(2))
    loom[output_layername] = sparselayer


def generate_malignancy_score(loom,
                              layername,
                              cell_sort_key='CellSort',
                              patient_id_key='patient_ID',
                              malignant_sort_label='45neg',
                              cell_name_key='cellname'):
    """

    For calculating malignancy scores for cells based on inferred CNV.  This subroutine isn't terribly future proof.  S Markson 6 June 2020.  

    Parameters
    ----------
    loom :
        
    layername :
        
    cell_sort_key :
         (Default value = 'CellSort')
    patient_id_key :
         (Default value = 'patient_ID')
    malignant_sort_label :
         (Default value = '45neg')
    cellname_key :
         (Default value = 'cellname')from panopticon.wme import get_list_of_gene_windows)
    robust_mean_windowed_expressionsfrom sklearn.decomposition import PCAfrom tqdm import tqdmcnv_scores_dict :
         (Default value = {}for patient in tqdm(np.unique(bm.ca[patient_id_key]))
    desc :
         (Default value = 'Computing per-patient)
    per-cell malignancy scores' :
        

    Returns
    -------

    """
    from panopticon.wme import get_list_of_gene_windows, robust_mean_windowed_expressions
    from sklearn.decomposition import PCA
    from tqdm import tqdm
    cnv_scores_dict = {}
    cnv_quantiles_dict = {}
    for patient in tqdm(
            np.unique(loom.ca[patient_id_key]),
            desc='Computing per-patient, per-cell malignancy scores'):
        mask = loom.ca['patient_ID'] == patient
        if malignant_sort_label in loom.ca[cell_sort_key][mask]:
            gene_windows = get_list_of_gene_windows(loom.ra['gene'])
            single_patient_expression = loom[layername][:, mask]
            mwe = robust_mean_windowed_expressions(
                loom.ra['gene'],
                gene_windows,
                single_patient_expression,
                upper_cut=2)
            pca = PCA(n_components=1)
            pca1 = pca.fit_transform(mwe.T)[:, 0]
            mask1 = loom.ca[cell_sort_key][mask] == malignant_sort_label
            mask2 = loom.ca[cell_sort_key][mask] != malignant_sort_label
            #if loom.ca[cell_sort_key][mask]
            if pca1[mask1].mean() > pca1[mask2].mean():
                scores45neg = np.sum(
                    np.greater.outer(pca1[mask1], pca1[mask2]),
                    axis=1) / np.sum(mask2)
                cnv_quantiles = (np.argsort(pca1) / np.sum(mask))
            elif pca1[mask1].mean() < pca1[mask2].mean():
                scores45neg = np.sum(
                    np.less.outer(pca1[mask1], pca1[mask2]),
                    axis=1) / np.sum(mask2)
                cnv_quantiles = (np.argsort(pca1) / np.sum(mask))[::-1]
            else:
                raise Exception(
                    "Unlikely event that CNV same for 45+/- has occurred")

        else:
            scores45neg = []
            cnv_quantiles = [np.nan] * np.sum(mask)
        counter45neg = 0
        for i, (cell_name, cell_sort) in enumerate(
                zip(loom.ca[cell_name_key][mask],
                    loom.ca[cell_sort_key][mask])):
            cnv_quantiles_dict[cell_name] = cnv_quantiles[i]
            if cell_sort == malignant_sort_label:
                cnv_scores_dict[cell_name] = scores45neg[counter45neg]
                counter45neg += 1
            else:
                cnv_scores_dict[cell_name] = 0
        if counter45neg != len(scores45neg):
            raise Exception("Some 45- cells unaccounted for")
    loom.ca['MalignantCNVScore'] = [
        cnv_scores_dict[x] for x in loom.ca[cell_name_key]
    ]
    loom.ca['MalignantCNVQuantile'] = [
        cnv_quantiles_dict[x] for x in loom.ca[cell_name_key]
    ]


def get_cosine_self_similarity(loom, layername, cluster, self_mean=None):
    from sklearn.metrics.pairwise import cosine_similarity
    clustering_level = len(str(cluster).split('-')) - 1
    mask = loom.ca['ClusteringIteration{}'.format(clustering_level)] == cluster
    if mask.sum() == 0:
        raise Exception("Mask is empty")
    if self_mean is None:
        self_mean = loom[layername][:, mask].mean(axis=1)
    return cosine_similarity(
        loom[layername][:, mask].T, Y=np.array([self_mean]))


def get_dictionary_of_cluster_means(loom, layername, clustering_level):
    from tqdm import tqdm
    mean_dict = {}
    for cluster in tqdm(
            np.unique(loom.ca[clustering_level]),
            desc='looping over clusters'):
        mask = loom.ca[clustering_level] == cluster
        if mask.sum() < 5000:
            mean_dict[cluster] = loom[layername][:, mask].mean(axis=1)
        else:
            mean_dict[cluster] = loom[layername].map([np.mean],
                                                     selection=mask)[0]

    return mean_dict


def get_differential_expression_dict(loom,
                                     layername,
                                     output=None,
                                     downsample_size=500):
    from panopticon.utilities import cluster_differential_expression
    from panopticon.utilities import we_can_pickle_it
    diffex = {}
    for i in range(5):
        for cluster in np.unique(loom.ca['ClusteringIteration{}'.format(i)]):
            print(cluster)
            diffex[cluster] = cluster_differential_expression(
                loom,
                'ClusteringIteration{}'.format(i),
                layername,
                ident1=cluster,
                ident1_downsample_size=downsample_size,
                ident2_downsample_size=downsample_size)
            if type(diffex[cluster]) != float:
                diffex[cluster] = diffex[cluster].query(
                    'MeanExpr1 >MeanExpr2').head(500)
                print(diffex[cluster].head(20))
            print('')
    if output is not None:
        we_can_pickle_it(diffex, output)
    return diffex
