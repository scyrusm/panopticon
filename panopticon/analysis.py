import numpy as np
import pandas as pd


def get_module_score_matrix(alldata, signature_mask, nbins=10, ncontrol=5):
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

    gene_quantiles = pd.qcut(alldata.mean(axis=1),
                             nbins,
                             duplicates='drop',
                             labels=False)
    sigdata_quantiles = gene_quantiles[signature_mask]
    nonsigdata_quantiles = gene_quantiles[~signature_mask]
    signature = sigdata.mean(axis=0)
    control_group = []
    for quantile in np.unique(sigdata_quantiles):
        noccurrences = np.sum(sigdata_quantiles == quantile)
        # there's an edge case wherein if size is greater than the number of genes to be taken without replacement, this will generate an error.  Will be an issue in the few-gene regime
        control_group += list(
            np.random.choice(np.where(nonsigdata_quantiles == quantile)[0],
                             size=ncontrol * noccurrences,
                             replace=False))

    control_group = np.array(control_group)
    control = nonsigdata[control_group].mean(axis=0)
    return signature - control


def generate_masked_module_score(loom, layername, mask, genelist, ca_name):
    """

    Parameters
    ----------
    loom : Name of loom object of interest.
        
    layername : Layername on which the module score will be calculated.
        
    mask : Mask over cells over which the score will be calculated ("None" for all cells)
        
    genelist : list of gene names in signature
        
    ca_name : Desired name of signature to be made into a column attribute.
        

    Returns
    -------

    
    """
    from panopticon.analysis import get_module_score_matrix
    if mask is None:
        mask = np.array([True] * loom.shape[1])
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
    """

    Parameters
    ----------
    loom :
        
    layername :
        
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
    generate_pca_loadings(loom, layername)


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

    import umap

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
    reducer = umap.UMAP(random_state=None,
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
                      silhouette_threshold=0.2,
                      regularization_factor=0.01,
                      clusteringcachedir='clusteringcachedir/'):
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

    Returns
    -------

    
    """
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import AgglomerativeClustering
    from tqdm import tqdm

    if X.shape[0] < min_input_size:
        return np.array([0] * X.shape[0])
    else:
        clustering = AgglomerativeClustering(n_clusters=2,
                                             memory=clusteringcachedir,
                                             affinity='cosine',
                                             compute_full_tree=True,
                                             linkage='average')
        scores = []
        minnk = 2
        for nk in tqdm(range(minnk, np.min([max_clusters, X.shape[0]]), 1)):
            clustering.set_params(n_clusters=nk)
            clustering.fit(X)

            score = silhouette_score(X,
                                     clustering.labels_,
                                     metric='cosine',
                                     sample_size=None)
            # sample_size=np.min([5000, X.shape[0]]))
            scores.append(score)
            #break
        print("scores", np.array(scores))
        print("ignoring regularization factor")
        #        scores = scores - np.arange(len(scores))*regularization_factor
        #        print("corrected scores",np.array(scores))
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
                        silhouette_threshold=0.1,
                        clusteringcachedir='clusteringcachedir/'):
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
    silhouette_threshold :
         (Default value = 0.1)

    Returns
    -------

    
    """
    if type(clustering_depth) != int or clustering_depth < 1 or type(
            starting_clustering_depth) != int:
        raise Exception(
            "clustering_depth and starting_clustering_depth must be natural numbers."
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
            n_pca_cols = loom.attrs['NumberPrincipalComponents_{}'.format(
                layername)]
            pca_loadings = []
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
                max_clusters=int(np.floor(np.sqrt(X.shape[0]))),
                                 clusteringcachedir=clusteringcachedir
            )  # This shouldn't be hard-coded S Markson 9 June 2020
        else:
            clustering = get_subclustering(
                X,
                silhouette_threshold,
                max_clusters=max_clusters,
                clusteringcachedir=clusteringcachedir
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
                model = NMF(n_components=np.min([50, data_c.shape[1]]),
                            init='random',
                            random_state=0)
                X = model.fit_transform(data_c.T)
            elif mode == 'pca':

                data_c = data_c.T
                if data_c.shape[0] > 5000:
                    model = IncrementalPCA(n_components=10)
                    for chunk in tqdm(
                            np.array_split(data_c,
                                           data_c.shape[0] // 512,
                                           axis=0),
                            desc='partial fitting over chunks of masked data'):
                        model.partial_fit(chunk)
                    X = model.transform(data_c)
                    print("EV", model.explained_variance_)
                    print("EVR", model.explained_variance_ratio_)
                else:
                    model = PCA(n_components=np.min([10, data_c.shape[0]]),
                                random_state=0)

                    X = model.fit_transform(data_c)
                    print("EV", model.explained_variance_)
                    print("EVR", model.explained_variance_ratio_)

            if max_clusters == 'sqrt_rule':
                print("xshape", X.shape)
                nopath_clustering = get_subclustering(
                    X,
                    silhouette_threshold,
                    max_clusters=int(np.floor(np.sqrt(X.shape[0]))),
                    clusteringcachedir=clusteringcachedir
                )  # This shouldn't be hard-coded S Markson 9 June 2020
            else:
                nopath_clustering = get_subclustering(
                    X,
                    silhouette_threshold,
                    max_clusters=max_clusters,
                    clusteringcachedir=clusteringcachedir
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


def cluster_differential_expression(loom,
                                    cluster_level,
                                    layername,
                                    ident1,
                                    ident2=None,
                                    mask1=None,
                                    mask2=None,
                                    verbose=False,
                                    ident1_downsample_size=None,
                                    ident2_downsample_size=None,
                                    min_cluster_size=0):
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
    min_cluster_size :
         (Default value = 0)

    Returns
    -------

    
    """
    from scipy.stats import mannwhitneyu
    from time import time
    from tqdm import tqdm

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
    print(np.sum(mask1), np.sum(mask2))
    pvalues = []
    genes = []
    meanexpr1 = []
    meanexpr2 = []
    fracexpr1 = []
    fracexpr2 = []
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
        fracexpr1.append((data1[igene, :] > 0).mean())
        fracexpr2.append((data2[igene, :] > 0).mean())
    output = pd.DataFrame(genes)
    output.columns = ['gene']
    output['pvalue'] = pvalues
    output['MeanExpr1'] = meanexpr1
    output['MeanExpr2'] = meanexpr2
    output['FracExpr1'] = fracexpr1
    output['FracExpr2'] = fracexpr2
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
    from panopticon.analysis import cluster_differential_expression
    diffex = {}
    for cluster in np.unique(loom.ca[cluster_level]):
        try:
            diffex[cluster] = cluster_differential_expression(
                loom,
                cluster_level,
                layername,
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
    if n_neighbors is None:
        n_neighbors = int(np.sqrt(np.sum(mask)))
    if genemask is None:
        data = loom[layername][:, mask]
    else:
        data = loom[layername][genemask, :][:, mask]

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
        mask = (loom.ca['ClusteringIteration{}'.format(cluster_level)]
                == cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    else:
        print("ignoring cluster, using custom mask")
        mask = (mask) & (loom.ca['nGene'] >= complexity_cutoff)
    return pd.DataFrame(loom.ca[field][mask])[0].value_counts()


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


def generate_malignancy_score(loom,
                              layername,
                              cell_sort_key='CellSort',
                              patient_id_key='patient_ID',
                              malignant_sort_label='45neg',
                              cell_name_key='cellname'):
    """For calculating malignancy scores for cells based on inferred CNV.  This subroutine isn't terribly future proof.  S Markson 6 June 2020.

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
        if malignant_sort_label in loom.ca[cell_sort_key][mask]:
            gene_windows = get_list_of_gene_windows(loom.ra['gene'])
            single_patient_expression = loom[layername][:, mask]
            mwe = robust_mean_windowed_expressions(loom.ra['gene'],
                                                   gene_windows,
                                                   single_patient_expression,
                                                   upper_cut=2)
            pca = PCA(n_components=1)
            pca1 = pca.fit_transform(mwe.T)[:, 0]
            mask1 = loom.ca[cell_sort_key][mask] == malignant_sort_label
            mask2 = loom.ca[cell_sort_key][mask] != malignant_sort_label
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
                from IPython.core.debugger import set_trace
                set_trace()
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
    """

    Parameters
    ----------
    loom :
        
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
        self_mean = loom[layername][:, mask].mean(axis=1)
    return cosine_similarity(loom[layername][:, mask].T,
                             Y=np.array([self_mean]))


def get_dictionary_of_cluster_means(loom, layername, clustering_level):
    """

    Parameters
    ----------
    loom :
        
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
            mean_dict[cluster] = loom[layername][:, mask].mean(axis=1)
        else:
            mean_dict[cluster] = loom[layername].map([np.mean],
                                                     selection=mask)[0]

    return mean_dict


def get_differential_expression_dict(loom,
                                     layername,
                                     output=None,
                                     downsample_size=500,
                                     starting_iteration=0,
                                     final_iteration=3,
                                     min_cluster_size=50):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    output :
         (Default value = None)
    downsample_size :
         (Default value = 500)
    starting_iteration :
         (Default value = 0)
    final_iteration :
         (Default value = 3)
    min_cluster_size :
         (Default value = 50)

    Returns
    -------

    """
    from panopticon.analysis import cluster_differential_expression
    from panopticon.utilities import we_can_pickle_it
    diffex = {}
    for i in range(starting_iteration, final_iteration + 1):
        for cluster in np.unique(loom.ca['ClusteringIteration{}'.format(i)]):
            print(cluster)
            diffex[cluster] = cluster_differential_expression(
                loom,
                'ClusteringIteration{}'.format(i),
                layername,
                ident1=cluster,
                ident1_downsample_size=downsample_size,
                ident2_downsample_size=downsample_size,
                min_cluster_size=min_cluster_size)
            if type(diffex[cluster]) != float:
                diffex[cluster] = diffex[cluster].query(
                    'MeanExpr1 >MeanExpr2').head(500)
                print(diffex[cluster].head(20))
            print('')
    if output is not None:
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
