import numpy as np
import pandas as pd


def get_module_score_matrix(loom,
                            layername,
                            cellmask,
                            signature_mask,
                            nbins=10,
                            ncontrol=5):
    """generates a module score (a la Seurat's AddModuleScore, see Tirosh 2016) on a matrix, with a mask.  I don't call this directly (S Markson 3 June 2020).

    Parameters
    ----------
    alldata : matrix
        
    signature_mask : indices corresponding to signature
        
    nbins : Number of quantile bins to use
        (Default value = 10)
    ncontrol : Number of genes in each matched quantile
        (Default value = 5)
    loom :
        
    layername :
        
    cellmask :
        

    Returns
    -------

    
    """
    assert len(signature_mask) == loom.shape[0]
    #nonsigdata = alldata[~signature_mask, :]
    #sigdata = alldata[signature_mask, :]
    if cellmask is None:
        gene_quantiles = pd.qcut(loom[layername].map(
            [np.mean],
            axis=0,
        )[0],
                                 nbins,
                                 duplicates='drop',
                                 labels=False)
    else:
        if len(cellmask) != loom.shape[1]:
            raise Exception(
                "cellmask must be boolean mask with length equal to the number of columns of loom"
            )
        gene_quantiles = pd.qcut(loom[layername].map(
            [np.mean], axis=0, selection=cellmask.nonzero()[0])[0],
                                 nbins,
                                 duplicates='drop',
                                 labels=False)
    sigdata_quantiles = gene_quantiles[signature_mask]
    nonsigdata_quantiles = gene_quantiles[~signature_mask]
    if len(signature_mask) != loom.shape[0]:
        raise Exception(
            "signature_mask must be boolean mask with length equal to the number of rows of loom"
        )
    signature = loom[layername].map([np.mean],
                                    axis=1,
                                    selection=signature_mask.nonzero()[0])[0]
    sigdata_quantiles = gene_quantiles[signature_mask]
    control_group = []
    for quantile in np.unique(sigdata_quantiles):
        noccurrences = np.sum(sigdata_quantiles == quantile)
        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
        control_group += list(
            np.random.choice(np.where(nonsigdata_quantiles == quantile)[0],
                             size=ncontrol * noccurrences,
                             replace=False))

    control_group = np.array(control_group)
    #control = nonsigdata[control_group].mean(axis=0)
    control = loom[layername].map([np.mean], axis=1,
                                  selection=control_group)[0]
    return signature - control


def generate_masked_module_score(loom,
                                 layername,
                                 cellmask,
                                 genelist,
                                 ca_name,
                                 nbins=10,
                                 ncontrol=5,
                                 gene_ra='gene'):
    """

    Parameters
    ----------
    loom : Name of loom object of interest.
        
    layername : Layername on which the module score will be calculated.
        
    cellmask : Mask over cells over which the score will be calculated ("None" for all cells)
        
    genelist : list of gene names in signature
        
    ca_name : Desired name of signature to be made into a column attribute.
        
    nbins :
        (Default value = 10)
    ncontrol :
        (Default value = 5)
    gene_ra :
        (Default value = 'gene')

    Returns
    -------

    
    """
    from panopticon.analysis import get_module_score_matrix
    #    if cellmask is None:
    #        cellmask = np.array([True] * loom.shape[1])
    #matrix = loom[layername][:, mask]
    sigmask = np.isin(loom.ra[gene_ra], genelist)
    sig_score = get_module_score_matrix(loom,
                                        layername,
                                        cellmask,
                                        sigmask,
                                        nbins=nbins,
                                        ncontrol=ncontrol)
    if cellmask is not None:
        maskedscores = []
        counter = 0
        for flag in cellmask:
            if flag:
                maskedscores.append(sig_score[counter])
                counter += 1
            else:
                maskedscores.append(np.nan)
    else:
        maskedscores = sig_score
    loom.ca[ca_name] = maskedscores


#def get_module_score_matrix(alldata, signature_mask, nbins=10, ncontrol=5):
#    """generates a module score (a la Seurat's AddModuleScore, see Tirosh 2016) on a matrix, with a mask.  I don't call this directly (S Markson 3 June 2020).
#
#    Parameters
#    ----------
#    alldata : matrix
#
#    signature_mask : indices corresponding to signature
#
#    nbins : Number of quantile bins to use
#        (Default value = 10)
#    ncontrol : Number of genes in each matched quantile
#        (Default value = 5)
#
#    Returns
#    -------
#
#
#    """
#    assert len(signature_mask) == alldata.shape[0]
#    nonsigdata = alldata[~signature_mask, :]
#    sigdata = alldata[signature_mask, :]
#
#    gene_quantiles = pd.qcut(alldata.mean(axis=1),
#                             nbins,
#                             duplicates='drop',
#                             labels=False)
#    sigdata_quantiles = gene_quantiles[signature_mask]
#    nonsigdata_quantiles = gene_quantiles[~signature_mask]
#    signature = sigdata.mean(axis=0)
#    control_group = []
#    for quantile in np.unique(sigdata_quantiles):
#        noccurrences = np.sum(sigdata_quantiles == quantile)
#        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
#        control_group += list(
#            np.random.choice(np.where(nonsigdata_quantiles == quantile)[0],
#                             size=ncontrol * noccurrences,
#                             replace=False))
#
#    control_group = np.array(control_group)
#    control = nonsigdata[control_group].mean(axis=0)
#    return signature - control
#
#
#def generate_masked_module_score(loom,
#                                 layername,
#                                 mask,
#                                 genelist,
#                                 ca_name,
#                                 nbins=10,
#                                 ncontrol=5,
#                                 gene_ra='gene'):
#    """
#
#    Parameters
#    ----------
#    loom : Name of loom object of interest.
#
#    layername : Layername on which the module score will be calculated.
#
#    mask : Mask over cells over which the score will be calculated ("None" for all cells)
#
#    genelist : list of gene names in signature
#
#    ca_name : Desired name of signature to be made into a column attribute.
#
#    nbins :
#        (Default value = 10)
#    ncontrol :
#        (Default value = 5)
#    gene_ra :
#        (Default value = 'gene')
#
#    Returns
#    -------
#
#
#    """
#    from panopticon.analysis import get_module_score_matrix
#    if mask is None:
#        mask = np.array([True] * loom.shape[1])
#    matrix = loom[layername][:, mask]
#    sigmask = np.isin(loom.ra[gene_ra], genelist)
#    sig_score = get_module_score_matrix(matrix,
#                                        sigmask,
#                                        nbins=nbins,
#                                        ncontrol=ncontrol)
#    maskedscores = []
#    counter = 0
#    for flag in mask:
#        if flag:
#            maskedscores.append(sig_score[counter])
#            counter += 1
#        else:
#            maskedscores.append(np.nan)
#    loom.ca[ca_name] = maskedscores


def generate_nmf_and_loadings(loom,
                              layername,
                              nvargenes=2000,
                              n_components=100,
                              verbose=False):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
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
    if len(vargenemask) != loom.shape[0]:
        raise Exception(
            "vargenemask must be boolean mask with length equal to the number of rows of loom"
        )
    X = loom[layername][vargenemask.nonzero()[0], :]
    model = NMF(n_components=n_components,
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


def generate_incremental_pca(loom,
                             layername,
                             batch_size=512,
                             n_components=50,
                             min_size_for_incrementalization=5000):
    """Computes a principal component analysis (PCA) over a layer of interest.  Defaults to incremental PCA (using IncrementalPCA from sklearn.decomposition) but will switch to conventional PCA for LoomConnections with cell
    numbers below a min_size_for_incrementalization.  Will write the n_components principal components as row attributes:
    - (layer) PC (PC number, 1-indexed)
    
    The following are written as attributes:
    - NumberPrincipalComponents_(layername).  This is simply n_components.
    - PCExplainedVariancedRatio_(layername).  This is explained_variance_ratio_ from the PCA model.
    
    Will also run panopticon.analysis.generate_pca_loadings.

    Parameters
    ----------
    loom : The LoomConnection instance upon which PCA will be calculated.
        
    layername : The layer of the loom file over which the PCs will be computed.
        
    batch_size :
        (Default value = 512)
    n_components :
        (Default value = 50)
    min_size_for_incrementalization :
        (Default value = 5000)

    Returns
    -------

    
    """

    from tqdm import tqdm
    from sklearn.decomposition import IncrementalPCA, PCA
    from panopticon.analysis import generate_pca_loadings
    batch_size_altered = False
    while loom.shape[1] % batch_size < n_components:
        batch_size += 1
        batch_size_altered = True
    if batch_size_altered:
        print(
            "Batch size increased to {} so that smallest batch will be greater than n_components"
            .format(batch_size))
    if loom.shape[1] < min_size_for_incrementalization:
        print(
            "Loom size below threshold for incremental PCA; running conventional PCA"
        )
        pca = PCA(n_components=n_components)
        pca.fit(loom[layername][:, :].T)
    else:
        pca = IncrementalPCA(n_components=n_components)
        for (ix, selection, view) in tqdm(loom.scan(axis=1,
                                                    batch_size=batch_size),
                                          total=loom.shape[1] // batch_size):
            #pca.partial_fit(view[:, :].transpose())

            pca.partial_fit(view[layername][:, :].T)
    for i in range(50):
        loom.ra['{} PC {}'.format(layername, i + 1)] = pca.components_[i]
    loom.attrs['NumberPrincipalComponents_{}'.format(layername)] = n_components
    loom.attrs['PCAExplainedVarianceRatio_{}'.format(
        layername)] = pca.explained_variance_ratio_
    generate_pca_loadings(loom, layername, batch_size=batch_size)


def generate_pca_loadings(loom, layername, dosparse=False, batch_size=1024):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    layername :
        
    dosparse :
        (Default value = False)
    batch_size :
        (Default value = 512)

    Returns
    -------

    
    """
    from tqdm import tqdm

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
        for (ix, selection, view) in tqdm(loom.scan(axis=1,
                                                    batch_size=batch_size),
                                          total=loom.shape[1] // batch_size):
            compresseddatalist.append(view[layername][:, :].T @ cellpca)
    compresseddata = np.vstack(compresseddatalist)
    for iloading in range(compresseddata.shape[1]):
        loom.ca['{} PC {} Loading'.format(layername, iloading +
                                          1)] = compresseddata[:, iloading]


def get_pca_loadings_matrix(loom, layername, n_components=None):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    layername : corresponding layer from which to retrieve PCA loadings matrix
        
    components_to_use :
        (Default value = None)
    n_components :
        (Default value = None)

    Returns
    -------

    
    """
    pca_loadings = []
    if n_components != None:
        for col in [
                '{} PC {} Loading'.format(layername, x)
                for x in range(1, n_components + 1)
        ]:
            pca_loadings.append(loom.ca[col])
    else:
        n_components = loom.attrs['NumberPrincipalComponents_{}'.format(
            layername)]
        for col in [
                '{} PC {} Loading'.format(layername, x)
                for x in range(1, n_components + 1)
        ]:
            pca_loadings.append(loom.ca[col])
    return np.vstack(pca_loadings).T


def generate_embedding(loom,
                       layername,
                       min_dist=0.0001,
                       n_neighbors=30,
                       n_epochs=1000,
                       metric='correlation',
                       random_state=None,
                       n_pca_components=None,
                       mode='pca'):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
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
    n_pca_components :
        (Default value = None)

    Returns
    -------

    
    """

    import umap

    if mode not in ['pca', 'nmf']:
        raise Exception("Currently only two modes implemented:  nmf and pca")
    if mode == 'pca':
        from panopticon.analysis import get_pca_loadings_matrix
        compressed = get_pca_loadings_matrix(loom,
                                             layername,
                                             n_components=n_pca_components)
    elif mode == 'nmf':
        n_nmf_cols = loom.attrs['NumberNMFComponents']
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
        compressed = np.vstack(nmf_loadings).T
    reducer = umap.UMAP(random_state=None,
                        min_dist=min_dist,
                        n_neighbors=n_neighbors,
                        metric=metric,
                        verbose=True,
                        n_epochs=n_epochs)
    embedding = reducer.fit_transform(compressed)
    loom.ca['{} {} UMAP embedding 1'.format(layername,
                                            mode.upper())] = embedding[:, 0]
    loom.ca['{} {} UMAP embedding 2'.format(layername,
                                            mode.upper())] = embedding[:, 1]


def get_subclustering(X,
                      score_threshold,
                      max_clusters=50,
                      min_input_size=10,
                      silhouette_threshold=0.2,
                      regularization_factor=0.01,
                      clusteringcachedir='clusteringcachedir/',
                      show_dendrogram=False,
                      linkage='average',
                      silhouette_score_sample_size=None):
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
    silhouette_threshold :
        (Default value = 0.2)
    regularization_factor :
        (Default value = 0.01)
    clusteringcachedir :
        (Default value = 'clusteringcachedir/')

    Returns
    -------

    
    """
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import AgglomerativeClustering
    from tqdm import tqdm

    if X.shape[0] < min_input_size:
        return np.array([0] * X.shape[0])
    else:
        print('Computing agglomerative clustering with cosine affinity, {} linkage'.format(linkage))
        clustering = AgglomerativeClustering(n_clusters=2,
                                             memory=clusteringcachedir,
                                             affinity='cosine',
                                             compute_full_tree=True,
                                             linkage=linkage)
        scores = []
        minnk = 2
        for nk in tqdm(range(minnk, np.min([max_clusters, X.shape[0]]), 1), desc='Computing silhouette scores'):
            clustering.set_params(n_clusters=nk)
            clustering.fit(X)
            if silhouette_score_sample_size is not None:
                silhouette_score_sample_size = np.min([X.shape[0],silhouette_score_sample_size])
            score = silhouette_score(X,
                                     clustering.labels_,
                                     metric='cosine',
                                     sample_size=silhouette_score_sample_size)
            # sample_size=np.min([5000, X.shape[0]]))
            scores.append(score)
            #break
        print("scores", np.array(scores))
        if show_dendrogram and hasattr(clustering, 'children_'):
            if linkage == 'average':
                method = 'weighted'
            else:
                method = linkage
            from scipy.cluster.hierarchy import dendrogram, linkage
            import matplotlib.pyplot as plt
            Z = linkage(X, metric='cosine', method=method)
            d = dendrogram(Z, truncate_mode='level', p=10)
            plt.show()

        if np.max(scores) >= score_threshold:
            clustering.set_params(n_clusters=np.argmax(scores) + minnk)
            clustering.fit(X)
            print(
                np.argmax(scores) + minnk, "clusters, with",
                list(
                    pd.DataFrame(clustering.labels_)[0].value_counts().values),
                "cells")
            return clustering.labels_
        else:
            print("No score exceeded score threshold")
            return np.array([0] * X.shape[0])


def generate_clustering(loom,
                        layername,
                        starting_clustering_depth=0,
                        n_clustering_iterations=3,
                        max_clusters='cbrt_rule',
                        mode='pca',
                        n_components=10,
                        silhouette_threshold=0.1,
                        clusteringcachedir='clusteringcachedir/',
                        out_of_core_batch_size=1024,
                        min_subclustering_size=100,
                        first_round_leiden=False,
                        optimized_leiden=True,
                        leiden_nneighbors=100,
                        leiden_iterations=10,
                        incremental_pca_threshold=10000,
                        show_dendrogram=False,
                        linkage='average'):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    final_clustering_depth : The clustering iteration on which to terminate; final_clustering_depth=3 will assign values through column attribute ClusteringIteration3
        (Default value = 3)
    starting_clustering_depth : The clustering iteration on which to begin; starting_clustering_depth=0 will assign values to column attribute ClusteringIteration0
        (Default value = 0)
    max_clusters :
        (Default value = 200)
    layername :
        
    mode :
        (Default value = 'pca')
    silhouette_threshold :
        (Default value = 0.1)
    clusteringcachedir :
        (Default value = 'clusteringcachedir/')
    n_components :
        (Default value = 10)
    out_of_core_batch_size :
        (Default value = 512)
    n_clustering_iterations :
        (Default value = 3)
    min_subclustering_size :
        (Default value = 50)
    first_round_leiden :
        (Default value = False)
    leiden_nneighbors :
        (Default value = 20)
    leiden_iterations :
        (Default value = 10)
    incremental_pca_threshold :
        (Default value = 10000)

    Returns
    -------

    
    """

    if type(n_clustering_iterations
            ) != int or n_clustering_iterations < 1 or type(
                starting_clustering_depth) != int:
        raise Exception(
            "final_clustering_depth and starting_clustering_depth must be natural numbers."
        )
    if (starting_clustering_depth > 0) and (
            'ClusteringIteration{}'.format(starting_clustering_depth - 1)
            not in loom.ca.keys()):
        raise Exception(
            "starting_clustering_depth not yet computed; please run with lower starting_clustering depth, or 0"
        )
    if mode not in ['pca', 'nmf']:
        raise Exception("Currently only implemented for modes:  pca and nmf")

    from time import time
    from sklearn.decomposition import IncrementalPCA
    from tqdm import tqdm
    from panopticon.analysis import get_subclustering
    if mode == 'pca':
        from sklearn.decomposition import PCA
    elif mode == 'nmf':
        from sklearn.decomposition import NMF

    final_clustering_depth = starting_clustering_depth + n_clustering_iterations - 1
    if starting_clustering_depth == 0:
        if first_round_leiden:
            from sklearn.neighbors import kneighbors_graph
            from panopticon.analysis import get_pca_loadings_matrix
            from panopticon.utilities import get_igraph_from_adjacency
            from panopticon.utilities import import_check

            X = get_pca_loadings_matrix(loom, layername, n_components=n_components)
            if optimized_leiden:
                from panopticon.clustering import silhouette_optimized_leiden
                leiden_output = silhouette_optimized_leiden(X)
            else:
                from panopticon.clustering import leiden_with_silhouette_score
                leiden_output = leiden_with_silhouette_score(X, leiden_nneighbors, leiden_iterations=leiden_iterations)
            clustering = leiden_output.clustering 

        else:
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
                from panopticon.analysis import get_pca_loadings_matrix
                X = get_pca_loadings_matrix(loom,
                                            layername,
                                            n_components=n_components)
            if max_clusters == 'sqrt_rule':
                mc = int(np.floor(np.sqrt(X.shape[0])))
            elif max_clusters == 'cbrt_rule':
                mc = int(np.floor(np.cbrt(X.shape[0])))
            else:
                mc = max_clusters
            clustering = get_subclustering(
                X,
                silhouette_threshold,
                max_clusters=mc,
                clusteringcachedir=clusteringcachedir,
                show_dendrogram=show_dendrogram,
                linkage=linkage)

        loom.ca['ClusteringIteration0'] = clustering
        starting_clustering_depth = 1

    for subi in range(starting_clustering_depth, final_clustering_depth + 1):

        loom.ca['ClusteringIteration{}'.format(subi)] = ['U'] * len(
            loom.ca['ClusteringIteration{}'.format(subi - 1)])

        for cluster in set([
                x for x in loom.ca['ClusteringIteration{}'.format(subi - 1)]
                if x != 'U'
        ]):  #
            mask = loom.ca['ClusteringIteration{}'.format(
                subi -
                1)] == cluster  #first mask, check for top level clustering
            #break
            #            mask = np.nonzero(mask) # workaround to accommodate h5py 3.* bug -- S. Markson 5 Aug 2021
            print("No. cells in cluster", mask.sum())
            if mode == 'nmf':
                print(
                    "NMF operation currently requires loading full layer into memory"
                )
                start = time()
                data_c = loom[layername][:, mask.nonzero()[0]]
                print("processing cluster", cluster, "; time to load: ",
                      time() - start, ", mask size: ", np.sum(mask))
                model = NMF(n_components=np.min([50, data_c.shape[1]]),
                            init='random',
                            random_state=0)
                X = model.fit_transform(data_c.T)
            elif mode == 'pca':

                if mask.sum() > incremental_pca_threshold:
#                    print(
#                        "Warning: running incremental PCA with loom.scan with masked items; this can lead to batches smaller than the number of principal components.  If computation fails, try adjusting out_of_core_batch_size, or raising incremental_pca_threshold."
#                    )
                    pca = IncrementalPCA(n_components=n_components)
                    slice_to_merge = None
                    for (ix, selection, view) in tqdm(
                            loom.scan(axis=1,
                                      batch_size=out_of_core_batch_size,
                                      items=mask.nonzero()[0],
                                      layers=[layername]),
                            total=loom.shape[1] // out_of_core_batch_size, desc='calculating masked incremental pca'):
                        if slice_to_merge is None:
                            #if view.shape[1]<n_components:
                            if view.shape[1]<out_of_core_batch_size:
                                slice_to_merge = view[layername][:,:].T
                            else:
                                pca.partial_fit(view[layername][:, :].T)
                        else:
#                            print("merging small slices: ",slice_to_merge.shape[0],view.shape[1])
                            slice_to_merge = np.vstack([slice_to_merge, view[layername][:,:].T])
#                            if slice_to_merge.shape[0]>=n_components:
                            if slice_to_merge.shape[0]>=out_of_core_batch_size:
#                                print('Merged slice: ',slice_to_merge.shape[0])
                                pca.partial_fit(slice_to_merge)
                                slice_to_merge = None
                    if slice_to_merge is not None:
                        if slice_to_merge.shape[0]>=n_components:
#                             print('Merged slice: ',slice_to_merge.shape[0])
                            pca.partial_fit(slice_to_merge)
                            slice_to_merge = None
                        else:
                            print("Warning: {} cells neglected.  Consider re-running with different out_of_core_batch_size".format(slice_to_merge.shape[0]))

    
                    compresseddatalist = []
                    for (ix, selection, view) in tqdm(
                            loom.scan(axis=1,
                                      batch_size=out_of_core_batch_size,
                                      items=mask.nonzero()[0],
                                      layers=[layername]),
                            total=loom.shape[1] // out_of_core_batch_size, desc='calculating masked incremental pca loadings'):
                        compresseddatalist.append(
                            view[layername][:, :].T @ pca.components_.T)
                    X = np.vstack(compresseddatalist)
                elif mask.sum() < min_subclustering_size:
                    X = np.array(
                        [None] * mask.sum()
                    )  # This is a hack to avoid computing PCA in cases where no clustering will be performed
                else:
                    data_c = loom[layername][:, mask.nonzero()[0]].T
                    model = PCA(n_components=np.min([10, data_c.shape[0]]),
                                random_state=0)

                    X = model.fit_transform(data_c)
            print(X.shape)
            if max_clusters == 'sqrt_rule':
                mc = int(np.floor(np.sqrt(X.shape[0])))
            elif max_clusters == 'cbrt_rule':
                mc = int(np.floor(np.cbrt(X.shape[0])))
            else:
                mc = max_clusters
            nopath_clustering = get_subclustering(
                X,
                silhouette_threshold,
                max_clusters=mc,
                clusteringcachedir=clusteringcachedir,
                min_input_size=min_subclustering_size,
                show_dendrogram=show_dendrogram,
                linkage=linkage)

            fullpath_clustering = [
                '{}-{}'.format(cluster, x) for x in nopath_clustering
            ]
            loom.ca['ClusteringIteration{}'.format(subi)][
                mask.nonzero()[0]] = fullpath_clustering
        loom.ca['ClusteringIteration{}'.format(subi)] = loom.ca[
            'ClusteringIteration{}'.format(
                subi)]  #This is to force the changes to save to disk


def get_cluster_markers(loom, layername, cluster_level):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    layername :
        
    cluster_level :
        

    Returns
    -------

    
    """
    from panopticon.analysis import get_cluster_differential_expression
    diffex = {}
    for cluster in np.unique(loom.ca[cluster_level]):
        try:
            diffex[cluster] = get_cluster_differential_expression(
                loom,
                layername,
                cluster_level=cluster_level,
                ident1=cluster,
                verbose=True,
                ident1_downsample_size=500,
                ident2_downsample_size=500).query('MeanExpr1 > MeanExpr2')
        except:
            print("Some issue processing cluster {}".format(cluster))
    return diffex


def get_cluster_embedding(loom,
                          layername,
                          cluster,
                          min_dist=0.01,
                          n_neighbors=None,
                          verbose=False,
                          mask=None,
                          genemask=None,
                          n_components_pca=50):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
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
    genemask :
        (Default value = None)
    n_components_pca :
        (Default value = 50)

    Returns
    -------

    
    """
    from sklearn.decomposition import PCA
    import umap

    clustering_level = len(str(cluster).split('-')) - 1
    if mask is None:
        mask = loom.ca['ClusteringIteration{}'.format(
            clustering_level)] == cluster
    else:
        print("Clustering over custom mask--ignoring cluster argument")
    if len(mask) != loom.shape[1]:
        raise Exception(
            "mask must be boolean mask with length equal to the number of columns of loom"
        )
    if n_neighbors is None:
        n_neighbors = int(np.sqrt(np.sum(mask)))
    if genemask is None:
        data = loom[layername][:, mask.nonzero()[0]]
    else:
        if len(genemask) != loom.shape[0]:
            raise Exception(
                "genemask must be boolean mask with length equal to the number of rows of loom"
            )
        data = loom[layername][genemask.nonzero()[0], :][:, mask.nonzero()[0]]

    pca = PCA(n_components=np.min([n_components_pca, data.shape[1]]))
    pca.fit(data[:, :].transpose())

    cellpca = (data.T @ pca.components_.T)
    reducer = umap.UMAP(random_state=17,
                        verbose=verbose,
                        min_dist=min_dist,
                        n_neighbors=n_neighbors,
                        metric='correlation')
    #    reducer.fit(cellpca, )
    embedding = reducer.fit_transform(cellpca)
    return embedding


def get_metafield_breakdown(loom,
                            cluster,
                            field,
                            complexity_cutoff=0,
                            mask=None):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
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
        mask = (loom.ca['ClusteringIteration{}'.format(cluster_level)]
                == cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    else:
        print("ignoring cluster, using custom mask")
        if len(mask) != loom.shape[1]:
            raise Exception(
                "mask must be boolean mask with length equal to the number of columns of loom"
            )
        mask = (mask) & (loom.ca['nGene'] >= complexity_cutoff)
    return pd.DataFrame(loom.ca[field][mask.nonzero()[0]])[0].value_counts()


def get_patient_averaged_table(loom,
                               patient_key='patient_ID',
                               column_attributes=[],
                               n_cell_cutoff=0):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
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


def generate_malignancy_score(loom,
                              layername,
                              cell_sort_key='CellSort',
                              patient_id_key='patient_ID',
                              malignant_sort_label='45neg',
                              cell_name_key='cellname'):
    """For calculating malignancy scores for cells based on inferred CNV.  This subroutine isn't terribly future proof.  S Markson 6 June 2020.

    Parameters
    ----------
    loom : LoomConnection object
        
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
        
    cell_name_key :
        (Default value = 'cellname')

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
        if malignant_sort_label in loom.ca[cell_sort_key][mask.nonzero()[0]]:
            gene_windows = get_list_of_gene_windows(loom.ra['gene'])
            single_patient_expression = loom[layername][:, mask.nonzero()[0]]
            mwe = robust_mean_windowed_expressions(loom.ra['gene'],
                                                   gene_windows,
                                                   single_patient_expression,
                                                   upper_cut=2)
            pca = PCA(n_components=1)
            pca1 = pca.fit_transform(mwe.T)[:, 0]
            mask1 = loom.ca[cell_sort_key][mask.nonzero()
                                           [0]] == malignant_sort_label
            mask2 = loom.ca[cell_sort_key][mask.nonzero()
                                           [0]] != malignant_sort_label
            #if loom.ca[cell_sort_key][mask]
            if pca1[mask1].mean() > pca1[mask2].mean():
                scores45neg = np.sum(np.greater.outer(pca1[mask1],
                                                      pca1[mask2]),
                                     axis=1) / np.sum(mask2)
                cnv_quantiles = (np.argsort(pca1) / np.sum(mask))
            elif pca1[mask1].mean() < pca1[mask2].mean():
                scores45neg = np.sum(np.less.outer(pca1[mask1], pca1[mask2]),
                                     axis=1) / np.sum(mask2)
                cnv_quantiles = (np.argsort(pca1) / np.sum(mask))[::-1]
            elif np.sum(
                    mask2
            ) == 0:  #  Case where we have no non-malignant reference cells
                scores45neg = [np.nan] * np.sum(mask)
                cnv_quantiles = [np.nan] * np.sum(mask)

            else:
                raise Exception(
                    "Unlikely event that CNV same for 45+/- has occurred")

        else:
            scores45neg = []
            cnv_quantiles = [np.nan] * np.sum(mask)
        counter45neg = 0
        for i, (cell_name, cell_sort) in enumerate(
                zip(loom.ca[cell_name_key][mask.nonzero()[0]],
                    loom.ca[cell_sort_key][mask.nonzero()[0]])):
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
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    layername :
        
    cluster :
        
    self_mean :
        (Default value = None)

    Returns
    -------

    
    """
    from sklearn.metrics.pairwise import cosine_similarity
    clustering_level = len(str(cluster).split('-')) - 1
    mask = loom.ca['ClusteringIteration{}'.format(clustering_level)] == cluster
    if mask.sum() == 0:
        raise Exception("Mask is empty")
    if self_mean is None:
        self_mean = loom[layername][:, mask.nonzero()[0]].mean(axis=1)
    return cosine_similarity(loom[layername][:, mask.nonzero()[0]].T,
                             Y=np.array([self_mean]))


def get_dictionary_of_cluster_means(loom, layername, clustering_level):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    layername :
        
    clustering_level :
        

    Returns
    -------

    
    """
    from tqdm import tqdm
    mean_dict = {}
    for cluster in tqdm(np.unique(loom.ca[clustering_level]),
                        desc='looping over clusters'):
        mask = loom.ca[clustering_level] == cluster
        if mask.sum() < 5000:
            mean_dict[cluster] = loom[layername][:, mask.nonzero()[0]].mean(
                axis=1)
        else:
            mean_dict[cluster] = loom[layername].map(
                [np.mean], selection=mask.nonzero()[0])[0]

    return mean_dict


def get_differential_expression_dict(loom,
                                     layername,
                                     output=None,
                                     downsample_size=500,
                                     starting_iteration=0,
                                     final_iteration=3,
                                     min_cluster_size=50,
                                     gene_alternate_name=None,
                                     verbose=True):
    """Runs get_cluster_differential_expression over multiple clustering iterations (From ClusteringIteration(x) to ClusteringIteration(y), inclusive, where x = starting_iteration, and y = final_iteration), where ident1 is a cluster, and ident2 is the set of all other clusters which differ only in the terminal iteration (e.g. if there are clusters 0-0, 0-1, and 0-2, 1-0, and 1-1, differential expression will compare 0-0 with 0-1 and 0-2, 0-1 with 0-0 and 0-2, etc).  Outputs a dictionary with each of these differential expression result, with key equal to ident1.

    Parameters
    ----------
    loom : LoomConnection object
        
    layername : layer key of loom, over which differential expression will be computed
        
    output : Optional filename whereto a .pkl object will be written with dictionary output, or an xlsx, with each key assigned to a separate sheet
        (Default value = None)
    downsample_size : Number of cells from each cluster to downsample to prior to running differential expression
        (Default value = 500)
    starting_iteration : if 0, will start with ClusteringIteration0, for example
        (Default value = 0)
    final_iteration : if 3, will continue to ClusteringIteration3, for example
        (Default value = 3)
    min_cluster_size : minimum size of clusters to consider (if one of clusters if below this threshold, will output nan instead of a differential expression dataframe for that particular key)
        (Default value = 50)
    gene_alternate_name :
        (Default value = None)
    verbose :
        (Default value = True)

    Returns
    -------

    
    """
    from panopticon.analysis import get_cluster_differential_expression
    from panopticon.utilities import we_can_pickle_it

    if gene_alternate_name is None and 'gene_common_name' in loom.ra.keys():
        gene_alternate_name = 'gene_common_name'

    diffex = {}
    for i in range(starting_iteration, final_iteration + 1):
        for cluster in np.unique(loom.ca['ClusteringIteration{}'.format(i)]):
            if verbose:
                print(cluster)
            diffex[cluster] = get_cluster_differential_expression(
                loom,
                layername,
                cluster_level='ClusteringIteration{}'.format(i),
                ident1=cluster,
                ident1_downsample_size=downsample_size,
                ident2_downsample_size=downsample_size,
                min_cluster_size=min_cluster_size,
                gene_alternate_name=gene_alternate_name)
            if type(diffex[cluster]) != float:
                diffex[cluster] = diffex[cluster].query(
                    'MeanExpr1 >MeanExpr2').head(500)
                if verbose:
                    print(diffex[cluster].head(20))
            if verbose:
                print('')
    if output is not None:
        if output.endswith('.xlsx'):
            try:
                import xlsxwriter
            except ImportError as e:
                print(
                    "xlsxwriter not installed; returning output without writer to excel file"
                )
                return diffex

            relkeys = [
                x for x in diffex.keys()
                if type(diffex[x]) == pd.core.frame.DataFrame
            ]
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                for key in relkeys:
                    prefix = '-'.join(str(key).split('-')[0:-1])
                    complement = ', '.join([
                        x for x in np.unique(loom.ca[
                            'ClusteringIteration{}'.format(
                                len(prefix.split('-')))])
                        if x != key and x.startswith(prefix)
                    ])
                    sheet_name = '{} up vs. {}'.format(key, complement)
                    if len(
                            sheet_name
                    ) >= 32:  # excel doesn't like it when sheets have names with more than 31 characters
                        complement = complement.replace(
                            '{}-'.format(prefix), '-')
                    sheet_name = '{} up vs. {}'.format(key, complement)
                    if len(sheet_name) >= 32:
                        sheet_name = '{} up'.format(key)
                    diffex[key].to_excel(writer,
                                         sheet_name=sheet_name,
                                         index=False)
                writer.save()

        else:
            we_can_pickle_it(diffex, output)
    return diffex


def scrna2tracer_mapping(scrna_cellnames, tracer_cellnames):
    """

    Parameters
    ----------
    scrna_cellnames :
        
    tracer_cellnames :
        

    Returns
    -------

    
    """
    # I hate everything about this--S Markson 7 September 2020
    tracer2scrna_name = {}
    for tracer_cellname in tracer_cellnames:
        if tracer_cellname in scrna_cellnames:
            tracer2scrna_name[tracer_cellname] = tracer_cellname
        else:
            for pattern in ['_Lym', '_nucseq', '_45pos', '_45neg']:
                if tracer_cellname + pattern in scrna_cellnames:
                    tracer2scrna_name[
                        tracer_cellname] = tracer_cellname + pattern
        if tracer_cellname not in tracer2scrna_name.keys(
        ) and tracer_cellname.startswith('M'):
            #print(tracer_cellname)
            pass
    return tracer2scrna_name


def get_cluster_differential_expression(loom,
                                        layername,
                                        cluster_level=None,
                                        ident1=None,
                                        ident2=None,
                                        mask1=None,
                                        mask2=None,
                                        verbose=False,
                                        ident1_downsample_size=None,
                                        ident2_downsample_size=None,
                                        min_cluster_size=0,
                                        gene_alternate_name=None):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    cluster_level :
        (Default value = None)
    layername :
        
    ident1 :
        (Default value = None)
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
    min_cluster_size :
        (Default value = 0)
    gene_alternate_name :
        (Default value = None)

    Returns
    -------

    
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import fdrcorrection

    from time import time
    from tqdm import tqdm

    if gene_alternate_name is None and 'gene_common_name' in loom.ra.keys():
        gene_alternate_name = 'gene_common_name'

    if (mask1 is not None) and (mask2 is not None):
        if verbose:
            print("ignoring ident1, ident2")

    elif (mask1 is not None) or (mask2 is not None):
        raise Exception(
            "Either both or neither of mask1, mask2 must be specified")
    else:
        if cluster_level is None:
            raise Exception(
                "cluster_level must be specified when running with cluster identities, i.e. without specifying an explicit mask"
            )
        if ident1 is None:
            raise Exception(
                "ident1 must be specified when running with cluster identities, i.e. without specifying an explicit mask"
            )
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
    if (np.sum(mask1) < min_cluster_size) or (np.sum(mask2) <
                                              min_cluster_size):
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
    if verbose:
        print('Group 1 size: ',np.sum(mask1), ', group 2 size: ',np.sum(mask2))
    pvalues = []
    uvalues = []
    genes = []
    meanexpr1 = []
    meanexpr2 = []
    meanexpexpr1 = []
    meanexpexpr2 = []
    fracexpr1 = []
    fracexpr2 = []
    start = time()
    data1 = loom[layername][:, mask1.nonzero()[0]]
    if verbose:
        print('First matrix extracted in', time() - start, 'seconds')
    start = time()
    data2 = loom[layername][:, mask2.nonzero()[0]]
    if verbose:
        print('Second matrix extracted', time() - start, 'seconds')
    for igene, gene in enumerate(
            tqdm(loom.ra['gene'], desc='Computing Mann-Whitney p-values')):
        genes.append(gene)
        if np.std(data1[igene, :]) + np.std(data2[igene, :]) < 1e-14:
            pvalues.append(1)
            uvalues.append(np.nan)
        else:
            mw = mannwhitneyu(data1[igene, :],
                              data2[igene, :],
                              alternative='two-sided')
            pvalues.append(mw.pvalue)
            uvalues.append(mw.statistic)
        meanexpr1.append(data1[igene, :].mean())
        meanexpr2.append(data2[igene, :].mean())
        meanexpexpr1.append(np.mean(2**data1[igene, :]))
        meanexpexpr2.append(np.mean(2**data2[igene, :]))
        fracexpr1.append((data1[igene, :] > 0).mean())
        fracexpr2.append((data2[igene, :] > 0).mean())
    output = pd.DataFrame(genes)
    output.columns = ['gene']
    output['pvalue'] = pvalues
    output['CommonLanguageEffectSize'] = np.array(uvalues) / (data1.shape[1] *
                                                              data2.shape[1])
    output['MeanExpr1'] = meanexpr1
    output['MeanExpr2'] = meanexpr2
    output['MeanExpExpr1'] = meanexpexpr1
    output['MeanExpExpr2'] = meanexpexpr2
    output['FracExpr1'] = fracexpr1
    output['FracExpr2'] = fracexpr2
    if gene_alternate_name is not None:
        gene2altname = {
            gene: altname
            for gene, altname in zip(loom.ra['gene'],
                                     loom.ra[gene_alternate_name])
        }
        altnames = [gene2altname[x] for x in genes]
        output['GeneAlternateName'] = altnames

    output = output.sort_values('CommonLanguageEffectSize', ascending=False)
    output['BenjaminiHochbergQ'] = fdrcorrection(output['pvalue'],
                                                 is_sorted=False)[1]
    return output


def get_differential_expression_over_continuum(loom,
                                               layer,
                                               mask,
                                               covariate,
                                               method='spearman',
                                               gene_alternate_name=None):
    """

    Parameters
    ----------
    loom : LoomConnection object
        
    layer :
        
    mask :
        
    covariate :
        
    method :
        (Default value = 'spearman')

    Returns
    -------

    
    """
    if method not in ['spearman', 'kendall']:
        raise Exception(
            "Requested method not implemented.  Only `kendall` or `spearman` currently available."
        )
    if np.sum(mask) != len(covariate):
        raise Exception(
            "Length of covariate vector does not match mask length.")
    if gene_alternate_name is None and 'gene_common_name' in loom.ra.keys():
        gene_alternate_name='gene_common_name'


    from tqdm import tqdm

    if method == 'spearman':
        from scipy.stats import spearmanr
        method_module = spearmanr
    elif method == 'kendall':
        from scipy.stats import kendalltau
        method_module = kendalltau
    corrs = []
    pvals = []
    genes = []
    X = loom[layer][:, mask.nonzero()[0]]
    for i, gene in enumerate(tqdm(loom.ra['gene'], desc='looping over genes')):
        #try:
        if np.std(X[i, :]) < 1e-14:
            pvals.append(1)
            corrs.append(0)
        else:

            result = method_module(X[i, :], covariate, nan_policy='omit')
            pvals.append(result.pvalue)
            corrs.append(result.correlation)
        genes.append(gene)
        #except:
        #    pass
    df = pd.DataFrame(genes)
    df.columns = ['gene']
    df['GeneAlternateName'] = loom.ra[gene_alternate_name]
    df['pval'] = pvals
    df['corr'] = corrs

    return df.sort_values('pval')


def get_differential_expression_custom(X1, X2, genes, axis=0):
    """

    Parameters
    ----------
    X1 :
        
    X2 :
        
    genes :
        
    axis :
        (Default value = 0)

    Returns
    -------

    
    """
    from tqdm import tqdm
    from scipy.stats import mannwhitneyu
    if (len(genes) != X1.shape[axis]) or ((len(genes) != X1.shape[axis])):
        raise Exception(
            "(len(genes)!=X1.shape[axis]) or ((len(genes)!=X1.shape[axis]))")
    if axis not in [0, 1]:
        raise Exception("axis must be 0 or 1")
    if axis == 1:
        X1 = X1.T
        X2 = X2.T
    pvalues = []
    uvalues = []
    meanexpr1 = []
    meanexpr2 = []
    meanexpexpr1 = []
    meanexpexpr2 = []
    fracexpr1 = []
    fracexpr2 = []

    for igene, gene in enumerate(
            tqdm(genes, desc='Computing Mann-Whitney p-values')):
        if np.std(X1[igene, :]) + np.std(X2[igene, :]) < 1e-14:
            pvalues.append(1)
            uvalues.append(np.nan)
        else:
            mw = mannwhitneyu(X1[igene, :],
                              X2[igene, :],
                              alternative='two-sided').pvalue
            pvalues.append(mw.pvalues)
            uvalues.append(mw.statistic)
        meanexpr1.append(X1[igene, :].mean())
        meanexpr2.append(X2[igene, :].mean())
        meanexpexpr1.append(np.mean(2**X1[igene, :]))
        meanexpexpr2.append(np.mean(2**X2[igene, :]))
        fracexpr1.append((X1[igene, :] > 0).mean())
        fracexpr2.append((X2[igene, :] > 0).mean())
    output = pd.DataFrame(genes)
    output.columns = ['gene']
    output['pvalue'] = pvalues
    output['CommonLanguageEffectSize'] = np.array(uvalues) / (X1.shape[1] *
                                                              X2.shape[1])
    output['MeanExpr1'] = meanexpr1
    output['MeanExpr2'] = meanexpr2
    output['MeanExpExpr1'] = meanexpexpr1
    output['MeanExpExpr2'] = meanexpexpr2
    output['FracExpr1'] = fracexpr1
    output['FracExpr2'] = fracexpr2
    return output.sort_values('pvalue', ascending=True)


def simpson(x, with_replacement=False):
    """For computing simpson index directly from counts (or frequencies, if with_replacement=True)

    Parameters
    ----------
    x :
        
    with_replacement :
        (Default value = False)

    Returns
    -------

    
    """
    total = np.sum(x)
    if with_replacement:
        return np.sum([(y / total) * (y / total) for y in x])
    else:
        return np.sum([(y / total) * ((y - 1) / (total - 1)) for y in x])


def generate_ca_simpson_index(loom,
                              ca,
                              blacklisted_ca_values=[],
                              second_ca=None,
                              output_name=None,
                              overwrite=False,
                              with_replacement=False):

    if output_name is None:
        raise Exception("output_name must be specified")
    if output_name in loom.ca.keys() and overwrite is False:
        raise Exception(
            "overwrite must be True to write over existing ca ({})".format(
                output_name))

    from panopticon.analysis import simpson
    ca2counts = pd.DataFrame(loom.ca[ca])[0].value_counts().to_dict()
    df = pd.DataFrame(loom.ca[ca], columns=[ca])
    if second_ca is None:
        df = df[~df[ca].isin(blacklisted_ca_values)]
        s = df[ca].value_counts().agg(simpson,
                                      with_replacement=with_replacement)
        loom.ca[output_name] = s
    else:
        df[second_ca] = loom.ca[second_ca]
        df = df[~df[ca].isin(blacklisted_ca_values)]
        s_dict = df.groupby(second_ca)[ca].value_counts().groupby(
            second_ca).agg(simpson,
                           with_replacement=with_replacement).to_dict()
        loom.ca[output_name] = [s_dict[x] for x in loom.ca[second_ca]]


def get_enrichment_score(genes,
                         geneset,
                         scores=None,
                         presorted=False,
                         return_es_curve=False,
                         return_pvalue=False,
                         n_pvalue_permutations=1000):
    """Returns an enrichment score (ES) in the manner of Subramanian et al. 2005 (https://doi.org/10.1073/pnas.0506580102).

    Parameters
    ----------
    genes :
        
    geneset :
        
    scores :
        (Default value = None)
    presorted :
        (Default value = False)
    return_es_curve :
        (Default value = True)
    return_pvalue :
        (Default value = False)
    n_pvalue_permutations :
        (Default value = 1000)

    Returns
    -------

    
    """
    from collections import namedtuple
    from tqdm import tqdm

    enrichment_score_output = namedtuple(
        "EnrichmentScoreOutput",
        "enrichment_score enrichment_score_curve p_value")

    if presorted is False and scores is None:
        raise Exception("Scores must be specified when sorted==False")
    elif presorted is False:
        genes = np.array(genes)[np.argsort(scores)[::-1]]
    running_es = []
    phit = 0
    pmiss = 0
    n_genes = len(genes)
    n_genes_in_geneset = len(np.intersect1d(genes, geneset))
    if n_genes_in_geneset == 0:
        raise Exception("Overlap of geneset with gene list is zero")

    for gene in genes:
        if gene in geneset:
            phit += 1 / n_genes_in_geneset
        else:
            pmiss += 1 / (n_genes - n_genes_in_geneset)

        running_es.append(phit - pmiss)
    es = running_es[np.argmax(np.abs(running_es))]

    if return_pvalue is False:
        if return_es_curve:
            return enrichment_score_output(es, running_es, None)
        else:
            return enrichment_score_output(es, None, None)
    else:
        null_enrichments = []
        for i in tqdm(range(n_pvalue_permutations)):
            null_enrichments.append(
                get_enrichment_score(np.random.permutation(genes),
                                     geneset,
                                     return_es_curve=False,
                                     return_pvalue=False,
                                     presorted=True).enrichment_score)

        pval = np.mean(np.abs(es) < np.abs(np.array(null_enrichments)))

        if return_es_curve:
            return enrichment_score_output(es, running_es, pval)
        else:
            return enrichment_score_output(es, None, pval)


def hutcheson_t(x, y):
    """

    Parameters
    ----------
    x :
        
    y :
        

    Returns
    -------

    """
    from collections import namedtuple
    from panopticon.utilities import import_check
    exit_code = import_check(
        "rpy2", 'pip install rpy2 (will also require installation of R)')
    if exit_code != 0:
        return
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    robjects.r('''                                      
    library(ecolTest)                                     
    hutcheson <- function(x,                              
                    y){
    
    Hutcheson_t_test(x,y)
    }
    ''')
    out = namedtuple('Hutcheson_t_test', [
        'statistic', 'parameter', 'p_value', 'estimate', 'null_value',
        'method', 'alternative', 'data_name'
    ])
    test = robjects.r['hutcheson'](np.array(x), np.array(y))
    return out(*[x[0] for x in test])


def generate_diffusion_coordinates(loom,
                                   layername,
                                   sigma,
                                   n_coordinates=10,
                                   verbose=False,
                                   metric='euclidean'):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    sigma :
        
    n_coordinates :
         (Default value = 10)
    verbose :
         (Default value = False)
    metric :
         (Default value = 'euclidean')

    Returns
    -------

    """
    from sklearn.metrics import pairwise_distances
    from numpy.linalg import eig
    if verbose:
        print("Calculating transition matrix...")
    distances = pairwise_distances(loom[layername][:, :].T, metric=metric)
    k = np.exp(-distances / sigma)
    p = np.sum(k, axis=1)
    transition_matrix = np.divide(
        k, p)  # transition_matrix.sum(axis=0) is all ones
    if verbose:
        print("Transition matrix calculation complete. Diagonalizing...")
    vals, vecs = eig(transition_matrix)
    #    loom.attrs['diffusion_sigma'] = sigma # May implement later
    for i in range(1, n_coordinates + 1):
        loom.ca['{} DC {}'.format(
            layername, i)] = vals[i] * np.matmul(transition_matrix, vecs[:, i])


def conditional_simpson(x, x_conditional, x_total, with_replacement=False):
    """For computing simpson index directly from counts (or frequencies, if with_replacement=True), where the first selected element is conditional on some feature

    Parameters
    ----------
    x :
        
    with_replacement :
         (Default value = False)
    x_conditional :
        
    x_total :
        

    Returns
    -------

    
    """
    #total = np.sum(x)
    total_conditional = np.sum(x_conditional)

    if with_replacement:
        return np.sum([(y / x_total) * (y_conditional / total_conditional)
                       for y, y_conditional in zip(x, x_conditional)])
    else:
        return np.sum([
            ((y - 1) / (x_total - 1)) * (y_conditional / total_conditional)
            for y, y_conditional in zip(x, x_conditional)
        ])


def get_cluster_enrichment_dataframes(x, y, data):
    """

    Parameters
    ----------
    x :
        
    y :
        
    data :
        

    Returns
    -------

    """
    from panopticon.utilities import phi_coefficient
    from collections import namedtuple
    from scipy.stats import fisher_exact
    fishers_exact_p_dict = {}
    phi_coefficient_dict = {}
    counts_dict = {}
    cluster_fraction_incluster_dict = {}
    cluster_fraction_ingroup_dict = {}
    for group in data[x].unique():
        fishers_exact_p_dict[group] = {}
        phi_coefficient_dict[group] = {}
        counts_dict[group] = {}
        cluster_fraction_incluster_dict[group] = {}
        cluster_fraction_ingroup_dict[group] = {}
        for cluster in data[y].unique():
            t11 = data[(data[x] == group)
                       & (data[y] == cluster)].shape[0]  #.sum()
            t12 = data[(data[x] == group)
                       & (data[y] != cluster)].shape[0]  #.sum()
            t21 = data[(data[x] != group) & (data[y] == cluster)].shape[0]
            t22 = data[(data[x] != group) & (data[y] != cluster)].shape[0]
            table = np.array([[t11, t12], [t21, t22]])

            fishers_exact_p_dict[group][cluster] = fisher_exact(table)[1]
            phi_coefficient_dict[group][cluster] = phi_coefficient(table)
            counts_dict[group][cluster] = t11
            cluster_fraction_incluster_dict[group][cluster] = t11 / (
                t11 + t21
            )  # this represents the fraction of cells in a cluster that are in a group
            cluster_fraction_ingroup_dict[group][cluster] = t11 / (
                t11 + t12
            )  # this represents the fraction of cells in a group that are in a cluster

    fishers_exact_p_df = pd.DataFrame.from_dict(fishers_exact_p_dict)
    phi_coefficient_df = pd.DataFrame.from_dict(phi_coefficient_dict)
    counts_df = pd.DataFrame.from_dict(counts_dict)
    cluster_fraction_incluster_df = pd.DataFrame.from_dict(
        cluster_fraction_incluster_dict)
    cluster_fraction_ingroup_df = pd.DataFrame.from_dict(
        cluster_fraction_ingroup_dict)

    ClusterEnrichment = namedtuple('ClusterEnrichment', [
        'FishersExactP', 'PhiCoefficient', 'Counts', 'FractionOfCluster',
        'FractionOfGroup'
    ])

    return ClusterEnrichment(fishers_exact_p_df, phi_coefficient_df, counts_df,
                             cluster_fraction_incluster_df,
                             cluster_fraction_ingroup_df)
