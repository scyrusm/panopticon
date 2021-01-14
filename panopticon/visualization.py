import numpy as np
import matplotlib.pyplot as plt


def plot_subclusters(loom,
                     layername,
                     cluster,
                     sublayers=1,
                     plot_output=None,
                     label_clusters=True,
                     complexity_cutoff=0,
                     downsample_to=500,
                     blacklist=[]):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    cluster : Cluster, or list of clusters
        
    sublayers :
        (Default value = 1)
    plot_output :
        (Default value = None)
    label_clusters :
        (Default value = True)
    complexity_cutoff :
        (Default value = 0)
    downsample_to :
         (Default value = 500)
    blacklist :
         (Default value = [])

    Returns
    -------

    
    """
    from panopticon.analysis import get_cluster_embedding
    import matplotlib.pyplot as plt
    if type(cluster) == str or type(cluster) == int:
        current_layer = len(str(cluster).split('-')) - 1
        supermask = (loom.ca['ClusteringIteration{}'.format(current_layer)]
                     == cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    elif type(cluster) == list or type(cluster) == np.ndarray:
        current_layer_options = [len(str(x).split('-')) - 1 for x in cluster]
        if len(np.unique(current_layer_options)) != 1:
            raise Exception("All clusters of interest must be from same layer")
        else:
            current_layer = np.unique(current_layer_options)[0]
        supermask = np.isin(
            loom.ca['ClusteringIteration{}'.format(current_layer)],
            cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    print(np.sum(supermask))
    keep_probability = np.min([1, downsample_to / np.sum(supermask)])
    supermask *= np.random.choice([False, True],
                                  size=len(supermask),
                                  p=[1 - keep_probability, keep_probability])
    # I kinda hate this (30 Apr 2020 smarkson)
    embedding = get_cluster_embedding(loom, layername, cluster, mask=supermask)
    if type(cluster) == str or type(cluster) == int:
        subclusters = [
            x for x in np.unique(loom.ca['ClusteringIteration{}'.format(
                current_layer + sublayers)])
            if str(x).startswith(str(cluster) + '-') and (x not in blacklist)
        ]
    elif (type(cluster) == list) or (type(cluster) == np.ndarray):
        subclusters = []
        for cluster_part in cluster:
            subclusters += [
                x for x in np.unique(loom.ca['ClusteringIteration{}'.format(
                    current_layer + sublayers)])
                if str(x).startswith(str(cluster_part) +
                                     '-') and (x not in blacklist)
            ]
    for subcluster in subclusters:
        mask = loom.ca['ClusteringIteration{}'.format(
            current_layer + sublayers)][supermask] == subcluster
        plt.scatter(embedding[mask, 0], embedding[mask, 1], label=subcluster)
        if label_clusters:
            plt.annotate(
                subcluster,
                (np.mean(embedding[mask, 0]), np.mean(embedding[mask, 1])))
    plt.legend(ncol=len(subclusters) // 14 + 1, bbox_to_anchor=(1.05, 0.95))
    if plot_output is not None:
        plt.savefig(plot_output, bbox_inches='tight')
    else:
        plt.show()
    return embedding


def plot_cluster_umap(loom,
                      layername,
                      cluster,
                      mask=None,
                      plot_output=None,
                      label_clusters=True,
                      complexity_cutoff=0,
                      downsample_to=500,
                      blacklist=[]):
    """

    Parameters
    ----------
    loom :
        
    layername :
        
    cluster : Cluster, or list of clusters
        
    sublayers :
        (Default value = 1)
    plot_output :
        (Default value = None)
    label_clusters :
        (Default value = True)
    complexity_cutoff :
        (Default value = 0)
    mask :
         (Default value = None)
    downsample_to :
         (Default value = 500)
    blacklist :
         (Default value = [])

    Returns
    -------

    
    """
    from panopticon.analysis import get_cluster_embedding
    import matplotlib.pyplot as plt
    if mask is not None:
        supermask = mask
    elif type(cluster) == str or type(cluster) == int:
        current_layer = len(str(cluster).split('-')) - 1
        supermask = (loom.ca['ClusteringIteration{}'.format(current_layer)]
                     == cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    elif type(cluster) == list or type(cluster) == np.ndarray:
        current_layer_options = [len(str(x).split('-')) - 1 for x in cluster]
        if len(np.unique(current_layer_options)) != 1:
            raise Exception("All clusters of interest must be from same layer")
        else:
            current_layer = np.unique(current_layer_options)[0]
        supermask = np.isin(
            loom.ca['ClusteringIteration{}'.format(current_layer)],
            cluster) & (loom.ca['nGene'] >= complexity_cutoff)
    print(np.sum(supermask))
    keep_probability = np.min([1, downsample_to / np.sum(supermask)])
    supermask *= np.random.choice([False, True],
                                  size=len(supermask),
                                  p=[1 - keep_probability, keep_probability])
    # I kinda hate this (30 Apr 2020 smarkson)
    embedding = get_cluster_embedding(loom, layername, cluster, mask=supermask)
    if plot_output is not None:
        plt.savefig(plot_output, bbox_inches='tight')
    else:
        from IPython.core.debugger import set_trace
        set_trace()
        plt.show()
    return embedding


def get_cluster_differential_expression_heatmap(loom,
                                                layer,
                                                clusteringlevel,
                                                diffex={},
                                                output=None,
                                                bracketwidth=3.5,
                                                ff_offset=0,
                                                ff_scale=7,
                                                bracket_width=3.5):
    """

    Parameters
    ----------
    loom :
        
    layer :
        
    clusteringlevel :
        
    diffex :
         (Default value = {})
    output :
         (Default value = None)
    bracketwidth :
         (Default value = 3.5)
    ff_offset :
         (Default value = 0)
    ff_scale :
         (Default value = 7)
    bracket_width :
         (Default value = 3.5)

    Returns
    -------

    """
    from panopticon.analysis import cluster_differential_expression
    import seaborn as sns
    import pandas as pd

    clusteredmask = []
    for cluster in np.unique(loom.ca[clusteringlevel]):
        mask = loom.ca[clusteringlevel] == cluster
        if mask.sum() > 2:
            clusteredmask.append(np.where(mask)[0])
    clusteredmask = np.hstack(clusteredmask)
    allgenes = []
    rawX = []
    clusters = np.unique(loom.ca[clusteringlevel])
    for cluster in clusters:
        mask = loom.ca[clusteringlevel] == cluster
        if mask.sum() > 2:
            if cluster not in diffex.keys():
                diffex[cluster] = cluster_differential_expression(
                    loom, clusteringlevel, layer, cluster)
            if np.sum(loom.ca[clusteringlevel] == cluster) > 1:
                genes = diffex[cluster][~diffex[cluster]['gene'].isin(
                    allgenes)].query('MeanExpr1 > MeanExpr2').query(
                        'FracExpr2<.9').head(10)['gene'].values
                genemask = np.isin(loom.ra['gene'], genes)
                rawX.append(loom[layer][genemask, :][:, clusteredmask])
                allgenes.append(genes)

    clusteredmask = np.hstack(clusteredmask)

    hmdf = pd.DataFrame(np.vstack(rawX))
    hmdf.index = np.hstack(allgenes)

    #fig, ax = plt.subplots(figsize=(5,12))
    fig = plt.figure(figsize=(5, 12))
    ax = fig.add_axes([0.4, 0.1, 0.4, .8])
    sns.heatmap(hmdf,
                cmap='coolwarm',
                yticklabels=1,
                xticklabels=False,
                ax=ax,
                cbar_kws={'label': r'log${}_2$(TP100k+1) Expression'})
    plt.ylabel('Gene')
    plt.xlabel('Cell')

    for i, cluster in enumerate(clusters):
        ax.annotate(
            cluster,
            xy=(-hmdf.shape[1] * 0.8,
                hmdf.shape[0] - ff_offset - 3.5 - i * ff_scale),
            xytext=(-hmdf.shape[1] * 0.9,
                    hmdf.shape[0] - ff_offset - 3.5 - i * ff_scale),
            ha='right',
            va='center',
            bbox=dict(boxstyle='square', fc='white'),
            arrowprops=dict(
                arrowstyle='-[, widthB={}, lengthB=1.5'.format(bracket_width),
                lw=2.0),
            annotation_clip=False)
    #plt.savefig("./figures/AllCD8TCellClustersHeatmap14Sept2020.pdf")
    if output is not None:
        plt.savefig(output)
    plt.show()
    return diffex


def plot_dotmap(loom,
                diffex,
                clusterlevel,
                topn=10,
                scale=1,
                output=None,
                title=None,
                keyblacklist=[],
                geneblacklist=[]):
    """

    Parameters
    ----------
    loom :
        
    diffex :
        
    clusterlevel :
        
    topn :
         (Default value = 10)
    scale :
         (Default value = 1)
    output :
         (Default value = None)
    title :
         (Default value = None)
    keyblacklist :
         (Default value = [])
    geneblacklist :
         (Default value = [])

    Returns
    -------

    """
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(figsize=(25, 5))
    print(diffex.keys())
    markers = [[]]
    for key in diffex.keys():
        markers.append(
            diffex[key][(~np.isin(diffex[key]['gene'], geneblacklist))
                        & (~np.isin(diffex[key]['gene'], np.hstack(markers)))
                        & (~diffex[key]['gene'].apply(lambda x: '.' in x))].
            query('MeanExpr1 > MeanExpr2')['gene'].head(topn).values)
    markers = np.hstack(markers)
    print(markers)
    markernames = markers
    key2x = {
        key: x
        for key, x in zip(diffex.keys(), range(len(diffex.keys())))
    }
    marker2y = {marker: x for marker, x in zip(markers, range(len(markers)))}
    exprs = []
    for marker in markers:
        markerindex = np.where(loom.ra['gene'] == marker)[0][0]
        for key in [key for key in diffex.keys() if key not in keyblacklist]:

            expr = loom['log2(TP100k+1)'][markerindex, ][loom.ca[clusterlevel]
                                                         == key].mean()
            exprs.append(expr)
            #plt.scatter(key2x[key],marker2y[marker],s=expr*scale,cmap='coolwarm',c=expr, )
            #print(expr)
    ## This is very dumb -- S. Markson 8 October 2020
    for marker in markers:
        markerindex = np.where(loom.ra['gene'] == marker)[0][0]
        #for key in [key for key in diffex.keys() if key!='1-0']:
        for key in diffex.keys():

            expr = loom['log2(TP100k+1)'][markerindex, ][loom.ca[clusterlevel]
                                                         == key].mean()
            sc = plt.scatter(marker2y[marker],
                             key2x[key],
                             s=expr * scale,
                             cmap='coolwarm',
                             c=expr,
                             vmin=np.min(exprs),
                             vmax=np.max(exprs))
            #print(expr)

    ax.set_xticklabels(markernames)
    ax.set_xticks(range(len(markers)))
    ax.set_yticklabels(
        [key for key in key2x.keys() if key not in keyblacklist])
    #print(np.array([val for val in diffex.values() if key not in keyblacklist]))
    ax.set_yticks(
        np.array([val for val in key2x.values() if key not in keyblacklist]))
    ax.set_xlim([-0.5, len(markers) - 0.5])
    ax.set_ylim([-0.5, len(key2x) - 0.5])
    for i in np.arange(-0.5, len(markers) - 0.5, topn)[1::]:
        plt.axvline(ls='--', x=i, color='k')
        #if flag:
    cbar = plt.colorbar(sc, )
    cbar.set_label(
        'Mean log(TP100k+1) Expression',
        rotation=90,
    )
    plt.xticks(rotation=90)
    if output is not None:
        plt.savefig(output)
    if title is not None:
        plt.title(title)
    plt.show()
