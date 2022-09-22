import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


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
        param layername:
    cluster :
        type cluster: Cluster, or list of clusters
    sublayers :
        Default value = 1)
    plot_output :
        Default value = None)
    label_clusters :
        Default value = True)
    complexity_cutoff :
        Default value = 0)
    downsample_to :
        Default value = 500)
    blacklist :
        Default value = [])
    layername :
        

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
            current_layer + sublayers)][supermask.nonzero()[0]] == subcluster
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
        param layername:
    cluster :
        type cluster: Cluster, or list of clusters
    sublayers :
        Default value = 1)
    plot_output :
        Default value = None)
    label_clusters :
        Default value = True)
    complexity_cutoff :
        Default value = 0)
    mask :
        Default value = None)
    downsample_to :
        Default value = 500)
    blacklist :
        Default value = [])
    layername :
        

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
    keep_probability = np.min([1, downsample_to / np.sum(supermask)])
    supermask *= np.random.choice([False, True],
                                  size=len(supermask),
                                  p=[1 - keep_probability, keep_probability])
    # I kinda hate this (30 Apr 2020 smarkson)
    embedding = get_cluster_embedding(loom, layername, cluster, mask=supermask)
    if plot_output is not None:
        plt.savefig(plot_output, bbox_inches='tight')
    else:
        plt.show()
    return embedding


import numpy as np
import matplotlib.pyplot as plt


def cluster_differential_expression_heatmap(
        loom,
        layer,
        clusteringlevel,
        diffex={},
        output=None,
        min_cluster_size=2,
        figsize=(5, 5),
        cbar_label=None,
        n_top_genes=10,
        gene_sort_criterion='CommonLanguageEffectSize',
        vmin=None,
        vmax=None,
        cmap='coolwarm',
        average_over_clusters=True,
        return_fig_ax=False,
        custom_gene_list=None,
        gene_ra='gene',
        rotate=False,
        cluster_blacklist=[]):
    """
    Generates a heatmap, with expression of marker genes displayed in heatmap form. Can also be used with hand-picked genes using the `custom_gene_list` argument, as well as with custom labels by setting `clusteringlevel` to be a column attribute representing the clusters of interest.
    When using a custom set of genes, this command will automatically cluster those genes using the seaborn `clustermap` command. 


    Parameters
    ----------
    loom : LoomConnection
        LoomConnection object for which the differential expression over clusters will be computed
        param layer:
    layer : str
        Specified loom layer to be used as values for heatmap.
    clusteringlevel : int or str
        Specifies column attribute of `loom` to be uses as clusters. Typical values might be `ClusteringIteration0`.  However, any column attribute can be used if the diffex dictionary object is pre-computed. (Default value = {})
    diffex : dict
        Specifies pre-computed output of get_cluster_differential expression.  Differential expression for clusters that don't exist as keys in `diffex` will be computed on-the-fly. (Default value = {})
    output : NoneType or str
        If str, then will write heatmap to file with filename `output`. (Default value = None)
    min_cluster_size : int
         Minimum number of cells in cluster for cluster to be included in plot. (Default value = 2)
    figsize : tuple
         Size of figure, in inches. (Default value = (5, 5)) :
    cbar_label : str
         Specifies the label of colorbar axis. (Default value = None)
    n_top_genes : int
         Specifies the number of top marker genes for each cluster to be displayed in heatmap. (Default value = 10)
    gene_sort_criterion : str
         (Default value = 'CommonLanguageEffectSize')
    vmin : float
         Sets the minimum value for heatmap/clustermap. None indicates no minimum value. (Default value = None)
    vmax : float
         Sets the maximum value for heatmap/clustermap. None indicates no maximum value. (Default value = None)
    cmap : str
         Sets the colormap to be used for heatmap/clustermap. Must be a valid matplotlib colormap. (Default value = 'coolwarm')
    average_over_clusters : bool
         If set to `True`, will average expression over all cells of a given cluster. If `False`, will group clusters together, without averaging. (Default value = True)
    return_fig_ax : bool
         If set to `True` will return a tuple (`fig`,`ax`) for the figure and axis respectively with the heatmap of interest. If a seaborn clustermap was generated, will instead return the full seaborn.clustermap output, wherefrom figure and axis objects can be accessed. (Default value = False)
    custom_gene_list : list, ndarry or NoneType
         If not None, specifies a custome gene list to use for heatmap/clustermap. (Default value = None)
    gene_ra : str
         Specifies the row attribute of LoomConnection `loom` indicating the gene name. (Default value = 'gene')

    Returns
    -------

    
    """
    from panopticon.analysis import get_cluster_differential_expression
    import seaborn as sns
    import pandas as pd

    if custom_gene_list is not None:
        print("Using custom gene list: ")
        print(
            "Removing the following genes:\n",
            np.unique(
                [x for x in custom_gene_list if x not in loom.ra[gene_ra]]))
        custom_gene_list = np.unique(
            [x for x in custom_gene_list if x in loom.ra[gene_ra]])
        print("Using the following genes:\n", custom_gene_list)

    if cbar_label is None:
        cbar_label = layer

    clusteredmask = []
    clusters = []
    cluster_cols = []
    for cluster in np.unique(loom.ca[clusteringlevel]):
        mask = loom.ca[clusteringlevel] == cluster
        if mask.sum() >= min_cluster_size:
            clusteredmask.append(np.where(mask)[0])
            clusters.append(cluster)
            cluster_cols += mask.sum() * [cluster]
    clusteredmask = np.hstack(clusteredmask)
    all_gene_common_names = []
    allgenes = []
    rawX = []
    clusters_above_min_size_threshold = []
    for cluster in clusters:
        mask = loom.ca[clusteringlevel] == cluster
        if mask.sum() >= min_cluster_size:
            if custom_gene_list is None:
                if cluster not in diffex.keys():
                    diffex[cluster] = get_cluster_differential_expression(
                        loom,
                        layer,
                        cluster_level=clusteringlevel,
                        ident1=cluster)
                if type(diffex[cluster]
                        ) != float and cluster not in cluster_blacklist:
                    genes = diffex[cluster][
                        ~diffex[cluster]['gene'].isin(  # need to fix this...
                            allgenes)].query(
                                'MeanExpr1 > MeanExpr2').sort_values(
                                    gene_sort_criterion,
                                    ascending=False).head(n_top_genes)[
                                        'gene'].values  # need to fix this...
                    genemask = np.isin(loom.ra[gene_ra], genes)
                    rawX.append(
                        loom[layer][genemask.nonzero()[0], :][:,
                                                              clusteredmask])
                    allgenes += list(genes)
                    all_gene_common_names += list(
                        loom.ra['gene_common_name']
                        [genemask])  # this shouldn't be hard-coded...
                    clusters_above_min_size_threshold.append(cluster)


#                else:
#                    clusters.remove(cluster)
            else:
                if type(custom_gene_list) != list and type(
                        custom_gene_list) != np.ndarray:
                    raise Exception("custom_gene_list must have type List")
                genemask = np.isin(loom.ra[gene_ra], custom_gene_list)
                rawX.append(
                    loom[layer][genemask.nonzero()[0], :][:, clusteredmask])
                allgenes.append(custom_gene_list)
                break
    if len(clusters_above_min_size_threshold) > 0:
        clusters = clusters_above_min_size_threshold
    clusteredmask = np.hstack(clusteredmask)

    hmdf = pd.DataFrame(np.vstack(rawX))
    #return hmdf, allgenes
    hmdf.index = all_gene_common_names
    hmdf.columns = cluster_cols
    if average_over_clusters:
        hmdf = hmdf.T.reset_index().groupby('index').mean().T
        hmdf = hmdf[clusters]
    if rotate:
        hmdf = hmdf.T
    if vmax is None and vmin is None:
        annotations = None
    else:
        annotations = np.round(hmdf.values, decimals=1)
        if vmax is not None and vmin is not None:
            mask = (annotations <= vmax) & (annotations >= vmin)
        elif vmax is not None:
            mask = (annotations <= vmax)
        elif vmin is not None:
            mask = (annotations >= vmin)
        annotations = annotations.astype(str)
        annotations[mask] = ""

    if custom_gene_list is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.4, 0.1, 0.4, .8])
        sns.heatmap(hmdf,
                    cmap=cmap,
                    yticklabels=1,
                    ax=ax,
                    cbar_kws={'label': cbar_label},
                    annot_kws={"fontsize": 8},
                    annot=annotations,
                    fmt='s',
                    vmin=vmin,
                    vmax=vmax)
        for i, cluster in enumerate(clusters):
            if rotate:
                xy = ((i + 0.5) / len(clusters), -.3)
                xytext = (
                    (i + 0.5) / len(clusters),
                    -.6,
                )
                rotation = 90
                # ff = 2.54*4
                scale = 2.54 * fig.get_size_inches()[0] * .4
            else:
                xy = (-.7, (len(clusters) - i - 0.5) / len(clusters))
                xytext = (-.9, (len(clusters) - i - 0.5) / len(clusters))
                rotation = 0
                # ff = 2.54
                scale = 2.54 * fig.get_size_inches()[1]
            ax.annotate(cluster,
                        xy=xy,
                        xytext=xytext,
                        ha='center',
                        va='center',
                        bbox=dict(boxstyle='square', fc='white'),
                        arrowprops=dict(
                            arrowstyle='-[, widthB={}, lengthB=1'.format(
                                scale / len(clusters)),
                            lw=2.0),
                        annotation_clip=False,
                        xycoords='axes fraction',
                        rotation=rotation)
        if rotate:
            ax.set_xlabel('Gene')
            if average_over_clusters:
                ax.set_ylabel('Cluster')
            else:
                ax.set_ylabel('Cell')
                ax.set_yticklabels([])
        else:
            ax.set_ylabel('Gene')
            if average_over_clusters:
                ax.set_xlabel('Cluster')
            else:
                ax.set_xlabel('Cell')
                ax.set_xticklabels([])
    else:
        g = sns.clustermap(hmdf,
                           cmap=cmap,
                           yticklabels=1,
                           cbar_kws={'label': cbar_label},
                           annot_kws={"fontsize": 8},
                           vmin=vmin,
                           vmax=vmax,
                           annot=annotations,
                           fmt='s',
                           col_cluster=False,
                           figsize=figsize)

    plt.tight_layout()

    if output is not None:
        plt.savefig(output)
    if return_fig_ax:
        if custom_gene_list is None:
            return fig, ax
        else:
            return g
    else:
        plt.show()


def swarmviolin(data,
                x,
                y,
                hue=None,
                ax=None,
                split=False,
                dodge=True,
                alpha=0.2,
                noswarm=False,
                violinplot_kwargs={},
                swarmplot_kwargs={},
                swarm_downsample_percentage=None,
                annotate_hue_pvalues=False,
                annotate_hue_effect_size=False,
                annotate_hue_n=False,
                annotate_hue_pvalue_fmt_str='p: {0:.2f}',
                annotate_hue_effect_size_fmt_str='es: {0:.2f}',
                annotate_hue_n_fmt_str='n: {}, {}',
                effect_size='cohensd',
                pvalue='mannwhitney',
                custom_annotation_dict={},
                custom_annotation_fontsize=6):
    """

    Parameters
    ----------
    data :
        type data: pandas dataframe in form that would be acceptable input to seaborn violin, swarmplots
    diffex :
        param x:
    y :
        type y: column of data to be used for violin, swarmplot y argument
    hue : column of data to be used for violin, swarmplot hue argument
        Default value = 10)
    ax : matplotlib axis
        Default value = None)
    split : split: split argument to be passed to violin, swarmplot
        Default value = True)
    alpha : alpha: alpha to be used for seaborn violinplot
        Default value = 0.2)
    violinplot_kwargs :
        Default value = {})
    swarmplot_kwargs :
        Default value = {})
    swarm_downsample_percentage :
        Default value = None)
    annotate_hue_pvalues :
        Default value = False)
    annotate_hue_pvalue_fmt_str :
        Default value = 'p: {0:.2f}')
    x :
        
    noswarm :
        (Default value = False)
    annotate_hue_effect_size :
        (Default value = False)
    annotate_hue_effect_size_fmt_str :
        (Default value = 'es: {0:.2f}')
    effect_size : If annotating the effect size between hues, this will set the relevant means of calculation.  Must be one of 'cohensd' or 'cles' for Cohen's d, or common language effect size, respectively.
        (Default value = 'cohensd')

    Returns
    -------

    
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    if hue is None:
        hue_order = None
    else:
        hue_order = data[hue].unique()
    sns.violinplot(data=data,
                   x=x,
                   hue=hue,
                   hue_order=hue_order,
                   y=y,
                   split=split,
                   inner='quartile',
                   cut=0,
                   alpha=alpha,
                   ax=ax,
                   width=.8,
                   **violinplot_kwargs)
    for violin in ax.collections:
        violin.set_alpha(alpha)
    if not noswarm:
        if swarm_downsample_percentage is None:
            sns.swarmplot(
                data=data,
                x=x,
                hue=hue,
                hue_order=hue_order,
                y=y,
                alpha=1,
                ax=ax,
                dodge=dodge,
                #          size=2.5,
                **swarmplot_kwargs)
        else:
            downsample_mask = np.random.choice(
                [True, False],
                size=data.shape[0],
                p=[
                    swarm_downsample_percentage,
                    1 - swarm_downsample_percentage
                ],
            )
            sns.swarmplot(data=data[downsample_mask],
                          x=x,
                          hue=hue,
                          hue_order=hue_order,
                          y=y,
                          dodge=dodge,
                          alpha=1,
                          ax=ax,
                          **swarmplot_kwargs)

    def legend_without_duplicate_labels(ax):
        """

        Parameters
        ----------
        ax :
            

        Returns
        -------

        
        """
        handles, labels = ax.get_legend_handles_labels()
        unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels))
                  if l not in labels[:i]]
        ax.legend(*zip(*unique), title='', loc=(1.05, .8))

    legend_without_duplicate_labels(ax)
    if annotate_hue_pvalues or annotate_hue_effect_size or annotate_hue_n and len(
            custom_annotation_dict.keys()) == 0:
        if len(data[hue].unique()) != 2:
            raise Exception(
                "hue must be a categorical variable with 2 unique values")
        from scipy.stats import mannwhitneyu, ttest_ind
        if np.issubdtype(data[y].dtype, np.number):
            category_col = x
            continuous_col = y
            vertical_violins = True
            ticklabels = ax.get_xmajorticklabels()
        else:
            category_col = y
            continuous_col = x
            vertical_violins = False
            ticklabels = ax.get_ymajorticklabels()
        for ticklabel in ticklabels:
            category = ticklabel.get_text()
            hue1, hue2 = data[hue].unique()

            a = data[(data[hue] == hue1)
                     & (data[category_col] == category)][continuous_col].values
            b = data[(data[hue] == hue2)
                     & (data[category_col] == category)][continuous_col].values
            a = np.array([x for x in a if not np.isnan(x)])
            b = np.array([x for x in b if not np.isnan(x)])
            mw = mannwhitneyu(a, b, alternative='two-sided')
            tt = ttest_ind(a, b)
            if pvalue == 'mannwhitney':
                pval = mw.pvalue
            elif pvalue == 'ttest':
                pval = tt[1]
            else:
                raise Exception(
                    "`pvalue` must be either \'mannwhitney\' or \'ttest\'")

            if effect_size == 'cohensd':
                from panopticon.utilities import cohensd
                es = cohensd(a, b)
            elif effect_size == 'cles':
                es = mw.statistic / len(a) / len(b)
            else:
                raise Exception(
                    "effect_size must be either \'cohensd\' or \'cles\'")
            annotation_string = ''
            if vertical_violins:
                if annotate_hue_pvalues:
                    annotation_string += annotate_hue_pvalue_fmt_str.format(
                        pval) + '\n'
                if annotate_hue_effect_size:
                    annotation_string += annotate_hue_effect_size_fmt_str.format(
                        es) + '\n'
                if annotate_hue_n:
                    annotation_string += annotate_hue_n_fmt_str.format(
                        len(a), len(b)) + '\n'
                ax.annotate(
                    annotation_string,
                    (ticklabel.get_position()[0], np.max(np.hstack((a, b)))),
                    ha='center',
                    va='bottom')
            else:
                if annotate_hue_pvalues:
                    annotation_string += ' ' + annotate_hue_pvalue_fmt_str.format(
                        pval)
                if annotate_hue_effect_size:
                    annotation_string += '\n ' + annotate_hue_effect_size_fmt_str.format(
                        es)
                if annotate_hue_n:
                    annotation_string += '\n' + annotate_hue_n_fmt_str.format(
                        len(a), len(b))

                ax.annotate(annotation_string, (
                    np.max(np.hstack((a, b))),
                    ticklabel.get_position()[1],
                ),
                            ha='left',
                            va='center')
    if len(custom_annotation_dict.keys()) > 0:
        if np.issubdtype(data[y].dtype, np.number):
            category_col = x
            continuous_col = y
            vertical_violins = True
            ticklabels = ax.get_xmajorticklabels()
        else:
            category_col = y
            continuous_col = x
            vertical_violins = False
            ticklabels = ax.get_ymajorticklabels()

        for ticklabel in ticklabels:
            annotation_string = ''
            category = ticklabel.get_text()
            annotation_pos = np.max(data[data[category_col]==category][continuous_col].values)

            if vertical_violins:
                if ticklabel.get_text() in custom_annotation_dict.keys():
                    annotation_string += '\n' + custom_annotation_dict[
                        ticklabel.get_text()]
                ax.annotate(annotation_string,
                            (ticklabel.get_position()[0], annotation_pos),
                            ha='center',
                            va='bottom',
                            fontsize=custom_annotation_fontsize)
            else:
                if ticklabel in custom_annotation_dict.keys():
                    annotation_string += '\n' + custom_annotation_dict[
                        ticklabel]

                ax.annotate(annotation_string, (
                    annotation_pos,
                    ticklabel.get_position()[1],
                ),
                            ha='left',
                            va='center',
                            fontsize=custom_annotation_fontsize)

    return ax


def volcano(diffex,
            ax=None,
            gene_column='gene',
            pval_column='pvalue',
            mean_expr_left='MeanExpr1',
            mean_expr_right='MeanExpr2',
            left_name='',
            right_name='',
            genemarklist=[],
            logfoldchange_importance_threshold=0.5,
            neglogpval_importance_threshold=5,
            title='',
            output=None,
            positions=None,
            show=True,
            gene_label_offset_scale=1):
    """

    Parameters
    ----------
    diffex :
        param ax:  (Default value = None)
    gene_column :
        Default value = 'gene')
    pval_column :
        Default value = 'pvalue')
    mean_expr_left :
        Default value = 'MeanExpr1')
    mean_expr_right :
        Default value = 'MeanExpr2')
    left_name : Genes toward the left in the volcano plot, which are upregulated in group2 of 'diffex'
        Default value = '')
    right_name : Genes toward the right in the volcano plot, which are upregulated in group1 of 'diffex'
        Default value = '')
    genemarklist :
        Default value = [])
    logfoldchange_importance_threshold :
        Default value = 0.5)
    neglogpval_importance_threshold :
        Default value = 5)
    title :
        Default value = '')
    output :
        Default value = None)
    positions :
        Default value = None)
    show :
        Default value = True)
    gene_label_offset_scale :
        Default value = 1)
    ax :
        (Default value = None)

    Returns
    -------

    
    """

    import matplotlib.pyplot as plt
    import matplotlib
    import matplotlib.patheffects as pe

    matplotlib.rcParams['axes.linewidth'] = 3
    if ax is None:
        fig, ax = plt.subplots(figsize=(
            5,
            5,
        ))
    neglogpvalues = -np.log(diffex[pval_column].values) / np.log(10)
    logfoldchange = diffex[mean_expr_left].values - diffex[
        mean_expr_right].values
    important_mask = (np.abs(logfoldchange) >
                      logfoldchange_importance_threshold)
    important_mask = important_mask & (neglogpvalues >
                                       neglogpval_importance_threshold)
    ax.scatter(logfoldchange[important_mask],
               neglogpvalues[important_mask],
               alpha=1,
               marker='.',
               s=3,
               c='b')
    ax.scatter(logfoldchange[~important_mask],
               neglogpvalues[~important_mask],
               alpha=.1,
               marker='.',
               s=2,
               c='b')
    maxx = np.nanmax(np.abs(logfoldchange))
    maxy = np.nanmax(neglogpvalues)
    xoffset = .1
    yoffset = .1
    offsetcounter = 0
    offsetcounter2 = 0
    #noteable_genes = diffex[diffex[gene_column].isin(genemarklist)].sort_values(pval_column)[gene_column].values
    if positions is None:
        positions = ['t'] * len(genemarklist)
    elif positions == 'auto':
        from sklearn.metrics import euclidean_distances
        annomask = np.isin(diffex[gene_column], genemarklist)
        distances = euclidean_distances(
            np.vstack((logfoldchange[annomask], neglogpvalues[annomask])).T)
        positions = []
        for i in range(1, len(distances)):
            mindist = np.min(distances[i, 0:i])
            if mindist < 0.2:
                positions.append('l')
            else:
                positions.append('r')

    for gene, position in zip(
            genemarklist,
            positions,
    ):
        genedf = diffex[diffex[gene_column] == gene]
        negpval = -np.log(genedf.iloc[0][pval_column]) / np.log(10)
        logfoldchange = genedf.iloc[0].MeanExpr1 - genedf.iloc[0].MeanExpr2
        ax.scatter(logfoldchange, negpval, marker='.', color='k')
        if position == 'b':
            ax.annotate(
                gene, (logfoldchange, negpval),
                (logfoldchange,
                 negpval + .015 * maxy * gene_label_offset_scale),
                va='bottom',
                ha='center',
                path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        elif position == 't':
            ax.annotate(
                gene, (logfoldchange, negpval),
                (logfoldchange,
                 negpval - .015 * maxy * gene_label_offset_scale),
                va='top',
                ha='center',
                path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        elif position == 'l':
            ax.annotate(
                gene, (logfoldchange, negpval),
                (logfoldchange + .03 * maxx * gene_label_offset_scale,
                 negpval),
                va='center',
                ha='left',
                path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        elif position == 'r':
            ax.annotate(
                gene, (logfoldchange, negpval),
                (logfoldchange - .03 * maxx * gene_label_offset_scale,
                 negpval),
                va='center',
                ha='right',
                path_effects=[pe.withStroke(linewidth=1, foreground="white")])
        else:
            raise Exception("invalid position character selection")
    plt.axvline(0, ls='--')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    ax.set_xlabel('log fold change\n' + '$\leftarrow$' + left_name + '\n' +
                  right_name + r'$\rightarrow$',
                  fontsize=14)

    ax.set_ylabel('-log' + r'${}_{10}$' + '(p-value)', fontsize=14)

    ax.set_xlim([-maxx * 1.04, maxx * 1.04])
    ax.set_ylim([0, maxy * 1.01])
    plt.tight_layout()
    ax.set_title(title)
    if output is not None:
        plt.savefig(output, bbox_inches='tight', transparent='true', dpi=300)
    if show:
        plt.show()


def repertoire_plot(x=None,
                    y=None,
                    data=None,
                    hue=None,
                    ax=None,
                    fig=None,
                    show=False,
                    output=None,
                    pre_normalize_by_cohort=False,
                    normalize=False,
                    piechart=False,
                    ylabel='',
                    color_palette=None,
                    stack_order='agnostic',
                    smear=False,
                    annotate_simpson=False):
    """Repertoire plot, designed for plotting cluster compositions or TCR repertoires as stacked bar plots or pies, with stack height indicating the size of a given TCR clone.  See https://doi.org/10.1101/2021.08.25.456956, Fig. 3e.  In this context, input should consist of a dataframe ('data'), with each row representing a cell.  Argument 'y' should be a column of 'data' representing the cell's clone or other grouping of cells.  Argument 'x' should be a column of 'data' representing the sample whence the cell came.

    Parameters
    ----------
    x : str
       Column of data indicating sample. (Default value = None)
    y :  str
        Column of data indicate clone. (Default value = None)
    data : pandas.DataFrame object, with necessary columns specified by arguments x, y.
        (Default value = None)
    hue : Optional column of data indicating additional grouping of samples
        Default value = None)
    ax : matplotlib matplotlib.axes._subplots.AxesSubplot object or array thereof, optional
        Default value = None)
    fig : matplotlib.figure.Figure object, optional
        Default value = None)
    show : if 'True', will run matplotlib.show() upon completion
        Default value = False)
    output : argument to matplotlib.savefig
        Default value = None)
    normalize :
        (Default value = False)
    piechart :
        (Default value = False)
    ylabel :
        (Default value = '')
    legend :
        (Default value = False)
    color_palette :
        (Default value = None)
    stack_order :
        (Default value = 'agnostic')

    Returns
    -------

    
    """
    if x is None:
        raise Exception("x is a required argument")
    if y is None:
        raise Exception("y is a required argument")
    if data is None:
        raise Exception("data is a required argument")
    if output is not None:
        if type(output) != str:
            raise Exception(
                "argument \"output\" must be a string, representing a desired filename"
            )
    #if piechar

    import matplotlib.pyplot as plt
    from tqdm import tqdm
    import matplotlib.transforms as transforms

    if ax is None and fig is not None:
        raise Exception(
            "argument \"ax\" must be included when argument \"fig\" is included"
        )
    elif fig is None and ax is not None:
        raise Exception(
            "argument \"ax\" must be included when argument \"fig\" is included"
        )
    if piechart and ax is not None:
        if type(ax) not in [np.ndarray, list]:
            raise Exception(
                "ax must be a list or array of matplotlib.axes._subplots.AxesSubplot objects when argument piechart==True"
            )
    if stack_order not in ['matched', 'agnostic']:
        raise Exception("stack_order must be one of \'matched\', \'agnostic\'")
    if color_palette is None:
        if stack_order == 'matched':
            import seaborn as sns
            color_palette = sns.palettes.color_palette('colorblind')
        elif stack_order == 'agnostic':
            color_palette = ['w']
    if (not smear and type(color_palette)
            == str) or (smear and type(color_palette) != str):
        raise Exception(
            "type of `color_palette` must be str if and only if `smear`==True")
    if smear and normalize:
        raise Exception(
            "Color smear only permissible in without normalization")

    all_heights = []

    if hue is None:
        if stack_order == 'agnostic':
            grouped_data = data.groupby(x)[y].value_counts()
        elif stack_order == 'matched':
            grouped_data = data.groupby(x)[y].value_counts(sort=False)

            # this code ensures that all categories are represented (even with 0) in all bars
            newname = "{}_{} count".format(x, y)  # to cover an edge case
            original_order = data[y].unique()
            grouped_data_df = grouped_data.rename(newname).reset_index()
            all_y = grouped_data_df[y].unique()
            grouped_data_df_dict = grouped_data_df.groupby(
                x)[y].unique().to_dict()
            for key in grouped_data_df_dict:
                for missing_element in np.setdiff1d(
                        all_y, grouped_data_df_dict[key]
                ):  # may fail if cluster name are of different types
                    zero_count_insert = pd.DataFrame(
                        np.array([key, missing_element, 0]).reshape(1, -1),
                        columns=grouped_data_df.columns).astype(
                            grouped_data_df.dtypes)
                    grouped_data_df = pd.concat(
                        [grouped_data_df, zero_count_insert], )
            y_to_original_order = {x: i for i, x in enumerate(original_order)}
            grouped_data_df = grouped_data_df.sort_values(
                y, key=np.vectorize(lambda x: y_to_original_order[x]))
            grouped_data = grouped_data_df.set_index([x, y])[newname]

        total = data.groupby(x)[y].unique().apply(lambda x: len(x)).max()
        groupings = data[x].unique()
        ind = np.arange(len(groupings))

    else:
        if stack_order == 'agnostic':
            grouped_data = data.groupby([x, hue])[y].value_counts()
        elif stack_order == 'matched':
            grouped_data = data.groupby([x, hue])[y].value_counts(sort=False)
            original_order = data[y].unique()
            newname = "{}_{}_{} count".format(x, y,
                                              hue)  # to cover an edge case
            grouped_data_df = grouped_data.rename(newname, ).reset_index()
            all_y = grouped_data_df[y].unique()
            grouped_data_df_dict = grouped_data_df.groupby(
                [x, hue])[y].unique().to_dict()

            for key in grouped_data_df_dict:
                for missing_element in np.setdiff1d(
                        all_y, grouped_data_df_dict[key]
                ):  # may fail if cluster name are of different types
                    zero_count_insert = pd.DataFrame(
                        np.array([key[0], key[1], missing_element,
                                  0]).reshape(1, -1),
                        columns=grouped_data_df.columns).astype(
                            grouped_data_df.dtypes)
                    grouped_data_df = pd.concat(
                        [grouped_data_df, zero_count_insert], )
            y_to_original_order = {x: i for i, x in enumerate(original_order)}
            grouped_data_df = grouped_data_df.sort_values(
                y, key=np.vectorize(lambda x: y_to_original_order[x]))
            grouped_data = grouped_data_df.set_index([x, hue, y])[newname]

        total = data.groupby([x,
                              hue])[y].unique().apply(lambda x: len(x)).max()
        groupings = list(data.groupby([x, hue]).groups.keys())
        ind = []
        counter = -1
        for i, grouping in enumerate(groupings):
            if grouping[0] != groupings[i - 1][0]:
                counter += 1
            ind.append(counter)
            counter += 1
        ind = np.array(ind)


# return grouped_data
    for grouping in groupings:
        heights = []
        for n in grouped_data[grouping].values:
            heights.append(n)
        heights = heights + [0] * (total - len(heights))
        heights = np.array(heights)[::-1]
        all_heights.append(heights)
    all_heights = np.vstack(all_heights)
    all_heights = all_heights[:, ::-1]
    if pre_normalize_by_cohort:
        all_heights = np.divide(all_heights, all_heights.sum(axis=0))
    if normalize:
        all_heights = np.divide(all_heights.T, all_heights.sum(axis=1)).T

    bottoms = np.array([0.0] * all_heights.shape[0])

    if (not piechart) and (ax is None) and (fig is None):
        if smear:
            fig, ax = plt.subplots(1,
                                   2,
                                   figsize=(5.5, 5),
                                   gridspec_kw={'width_ratios': [10, 1]})
        else:
            fig, ax = plt.subplots(figsize=(5, 5))
    elif (piechart) and (ax is None) and (fig is None):
        import math
        ngroups = len(groupings)

        nrows = int(np.round(ngroups / np.sqrt(ngroups)))
        ncols = math.ceil(ngroups / np.sqrt(ngroups))
        if smear:
            print("Constructing subplots with {} rows, {} columns".format(
                nrows, 2 * ncols))
            fig, ax = plt.subplots(
                nrows,
                2 * ncols,
                figsize=(5.5, 5),
                gridspec_kw={'width_ratios': [10, 1] * ncols})
        else:
            print("Constructing subplots with {} rows, {} columns".format(
                nrows, ncols))
            fig, ax = plt.subplots(nrows, ncols, figsize=(5, 5))

    if piechart:
        for i in tqdm(range(all_heights.shape[0])):
            if stack_order == 'matched':
                labels = grouped_data[grouping].index.values
            else:
                labels = None

            if smear is True:
                import seaborn as sns
                countrange = range(1, np.max(all_heights[i, :] + 1))
                if type(color_palette) == str:
                    count2color = {
                        i: sns.color_palette(color_palette,
                                             len(countrange))[i - 1]
                        for i in countrange
                    }
                    count2color[0] = count2color[1]
                else:
                    raise Exception("`color_palette` must be of type str")

                colors = [count2color[x] for x in all_heights[i, :]]
                wedgeprops = {}
            else:
                colors = color_palette
                wedgeprops = {"edgecolor": "k"}
            if len(np.shape(ax)) == 1:
                if smear:
                    subax = ax[i * 2]
                    cbarax = ax[i * 2 + 1]
                else:
                    subax = ax[i]
                    cbarax = None
            elif len(np.shape(ax)) == 2:
                maxrows, maxcols = np.shape(ax)
                if smear:
                    if 2 * all_heights.shape[0] // maxrows > maxcols:
                        raise Exception(
                            "Insufficient number of subplots for number of grouping (pies)"
                        )
                    subax = ax[2 * i // maxcols, 2 * i % maxcols]
                    cbarax = ax[(2 * i + 1) // maxcols, (2 * i + 1) % maxcols]
                    #edgecolor='k')
                else:
                    if all_heights.shape[0] // maxrows > maxcols:
                        raise Exception(
                            "Insufficient number of subplots for number of grouping (pies)"
                        )
                    subax = ax[i // maxcols, i % maxcols]
                    cbarax = None
                    #edgecolor='k')
            title = groupings[i]
            if annotate_simpson:
                from panopticon.analysis import simpson
                si = simpson(all_heights[i, :])
                title += '\nSI: {0:.5f}'.format(si)
            subax.set_title(title)

            subax.pie(all_heights[i, :],
                      colors=colors,
                      wedgeprops=wedgeprops,
                      labels=labels)
            if cbarax is not None:
                import matplotlib as mpl
                from matplotlib.colors import LinearSegmentedColormap

                norm = mpl.colors.Normalize(vmin=0, vmax=np.max(countrange))
                colorlist = [
                    count2color[x]
                    for x in range(1, 1 + np.max(list(count2color.keys())))
                ]
                cm = LinearSegmentedColormap.from_list('custom',
                                                       colorlist,
                                                       N=len(colorlist))

                cb1 = mpl.colorbar.ColorbarBase(cbarax,
                                                cmap=cm,
                                                norm=norm,
                                                orientation='vertical')
                cb1.set_label('Clone size (No. cells)')
                if len(colorlist) + 1 >= 5:
                    ticklabels = np.arange(1,
                                           len(colorlist) + 1,
                                           (len(colorlist) + 1) // 5)
                else:
                    ticklabels = np.arange(1, len(colorlist) + 1, 1)
                ticks = [x - 0.5 for x in ticklabels]
                cb1.set_ticks(ticks)
                cb1.set_ticklabels(ticklabels)

    else:  # conventional stacked bar plot mode
        if smear:
            import seaborn as sns
            countrange = range(1, np.max(all_heights) + 1)
            if type(color_palette) == str:
                count2color = {
                    i: sns.color_palette(color_palette, len(countrange))[i - 1]
                    for i in countrange
                }
                count2color[0] = count2color[1]
            else:
                raise Exception("`color_palette` must be of type str")
        if smear:
            subax = ax[0]
            cbarax = ax[1]
        else:
            subax = ax
        for i in tqdm(range(all_heights.shape[1])):

            if stack_order == 'matched':
                label = grouped_data[grouping].index[i]
            else:
                label = None

            if smear:
                subax.bar(ind,
                          all_heights[:, i],
                          0.8,
                          bottom=bottoms,
                          color=[count2color[h] for h in all_heights[:, i]],
                          label=label)
                if cbarax is not None:
                    import matplotlib as mpl
                    from matplotlib.colors import LinearSegmentedColormap
                    norm = mpl.colors.Normalize(vmin=0,
                                                vmax=np.max(countrange))
                    colorlist = [
                        count2color[x]
                        for x in range(1, 1 + np.max(list(count2color.keys())))
                    ]
                    cm = LinearSegmentedColormap.from_list('custom',
                                                           colorlist,
                                                           N=len(colorlist))
                    cb1 = mpl.colorbar.ColorbarBase(cbarax,
                                                    cmap=cm,
                                                    norm=norm,
                                                    orientation='vertical')
                    cb1.set_label('Clone size (No. cells)')
                    if len(colorlist) + 1 >= 5:
                        ticklabels = np.arange(1,
                                               len(colorlist) + 1,
                                               (len(colorlist) + 1) // 5)
                    else:
                        ticklabels = np.arange(1, len(colorlist) + 1, 1)
                    ticks = [x - 0.5 for x in ticklabels]
                    cb1.set_ticks(ticks)
                    cb1.set_ticklabels(ticklabels)
            else:
                subax.bar(ind,
                          all_heights[:, i],
                          0.8,
                          bottom=bottoms,
                          color=color_palette[i % len(color_palette)],
                          edgecolor='k',
                          label=label)
                subax.set_xticks(ind)

            bottoms += all_heights[:, i]
            subax.set_xticks(ind)
        if hue is None:
            subax.set_xticklabels(groupings, rotation=90)
        else:
            subax.set_xticklabels([grouping[1] for grouping in groupings],
                                  rotation=90)
            plt.tight_layout()

            xpositions = []
            labels = []
            xpos = 0
            for i, grouping in enumerate(groupings):
                if grouping[0] != groupings[i - 1][0]:
                    xpositions.append(i)
                    labels.append(grouping[0])
            xpositions.append(len(groupings))
            xpositions = np.array(xpositions)[0:-1] + np.diff(xpositions) * 0.5
            xpositions += np.arange(len(xpositions))
            xpositions -= 0.5
            for xposition, label in zip(xpositions, labels):
                trans = transforms.blended_transform_factory(
                    subax.transData, fig.transFigure)

                subax.annotate(label, (xposition, 0),
                               annotation_clip=False,
                               ha='center',
                               va='bottom',
                               xycoords=trans,
                               fontsize=13)
        if ylabel == '':
            if normalize:
                ylabel = 'cell fraction'
            else:
                ylabel = 'cell count'
        subax.set_ylabel(ylabel)
    if stack_order == 'matched':
        plt.legend(bbox_to_anchor=(1, 1))
    plt.tight_layout()
    if output is not None:

        if output.endswith('.png'):
            plt.savefig(output, dpi=300)
        else:
            plt.savefig(output)
    if show:
        plt.show()


def data_to_grid_kde(x, y, xmin=None, xmax=None, ymin=None, ymax=None, px=100):
    """

    Parameters
    ----------
    x : Vector
        
    y : Vector
        
    px :
        (Default value = 100)
    xmin :
        (Default value = None)
    xmax :
        (Default value = None)
    ymin :
        (Default value = None)
    ymax :
        (Default value = None)

    Returns
    -------

    
    """
    from scipy.stats import gaussian_kde
    xy = np.vstack((x, y))
    kde = gaussian_kde(xy)
    if np.any([z is None for z in [xmin, xmax, ymin, ymax]]):
        if np.all([z is None for z in [xmin, xmax, ymin, ymax]]):
            xmin, xmax = (x.min(), x.max())
            ymin, ymax = (y.min(), y.max())
        else:
            raise Exception(
                "Either all of xmin, xmax, ymin, ymax must be \"None\", or none may be."
            )
    xgrid = np.arange(xmin, xmax, (xmax - xmin) / px)
    ygrid = np.arange(ymin, ymax, (ymax - ymin) / px)
    xx, yy = np.meshgrid(xgrid, ygrid)
    zz = kde.evaluate(
        np.array([xx.reshape(-1, 1)[:, 0],
                  yy.reshape(-1, 1)[:, 0]])).reshape(xx.shape, )
    return zz


def plot_differential_density(x, y, mask1, mask2, ax=None, cmap=plt.cm.RdBu_r):
    """

    Parameters
    ----------
    x :
        
    y :
        
    mask1 :
        
    mask2 :
        
    ax :
        (Default value = None)
    cmap :
        (Default value = plt.cm.RdBu_r)

    Returns
    -------

    
    """
    xmin, xmax = (x.min(), x.max())
    ymin, ymax = (y.min(), y.max())

    import numpy as np
    from panopticon.visualization import data_to_grid_kde
    zz1 = data_to_grid_kde(x[mask1],
                           y[mask1],
                           xmin=xmin,
                           ymin=ymin,
                           xmax=xmax,
                           ymax=ymax)
    from panopticon.visualization import data_to_grid_kde
    zz2 = data_to_grid_kde(x[mask2],
                           y[mask2],
                           xmin=xmin,
                           ymin=ymin,
                           xmax=xmax,
                           ymax=ymax)
    vmin = np.min(zz1 - zz2)
    vmax = np.max(zz1 - zz2)
    if np.abs(vmax) > np.abs(vmin):
        vmin = -vmax
    else:
        vmax = -vmin
    if ax is None:
        plt.imshow(zz1 - zz2, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        plt.show()
    else:
        ax.imshow(zz1 - zz2, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)


def plot_density(x, y, ax=None, cmap=plt.cm.twilight_r):
    """

    Parameters
    ----------
    x :
        
    y :
        
    ax :
        (Default value = None)
    cmap :
        (Default value = plt.cm.twilight_r)

    Returns
    -------

    
    """

    from panopticon.visualization import data_to_grid_kde
    zz = data_to_grid_kde(x, y)
    if ax is None:
        plt.imshow(zz, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        plt.show()
    else:
        ax.imshow(zz, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)


def cluster_enrichment_heatmap(x,
                               y,
                               data,
                               show=True,
                               output=None,
                               fig=None,
                               cax=None,
                               ax=None,
                               side_annotation=True,
                               heatmap_shading_key='FractionOfCluster',
                               annotation_key='Counts',
                               annotation_fmt='.5g',
                               figsize=(5, 5)):
    """
    Produces a heatmap indicating the fraction of cell clusters across groups.  For example, if there are `m` experimental groups and `n` clusters of cells, will produce a heatmap with
    `n` rows and `m` columns. 
    `heatmap_shading_key` can be any field of the named tuple output of `panopticon.analysis.get_cluster_enrichment_dataframes`. These include:

    - If `heatmap_shading_key` = "FractionOfCluster", heatmap color will be row-normalized; that is, it will indicate the fraction of cells in a cluster that are in groups.
    - If `heatmap_shading_key` = "FractionOfGroup", heatmap color will be column-normalized; that is, it will indicate the fraction of cells in a group that are in a given cluster. 
    - If `heatmap_shading_key` = "Counts", heatmap color will depict raw counts.
    - If `heatmap_shading_key` = "PhiCoefficient", heatmap color will depict phi-coefficients (described below).
    - If `heatmap_shading_key` = "FishersExactP", heatmap color will depict Fisher's exact test p-values (described below).

    P-values and phi-coefficients are computed by constructing the contigency matrices as follows:

    .. list-table:: 
       :widths: 25 25
    
       * - a
         - b
       * - c
         - d

    where `a` represents counts in cluster (**not** normalized) in group, `b` counts not in cluster in group, `c` counts in group, not in cluster, and `d` counts not in group, not in cluster. This is most intuitive for two groups, but can be computed in all cases (margins of the contingency matrix will be unchanged).
    P-values are computed via `scipy.stats.fisher_exact`, and effect sizes by phi coefficient (`panopticon.utilities.phi_coefficient`).

    If there are only two groups, side annotation can also be used in order to display counts, normalized counts, fisher's exact p-values and phi-coefficients all on one plot (see https://doi.org/10.1101/2021.08.25.456956, Fig. 4c, f and Fig. 5c). This will only work in the two-group case however. 

    Parameters
    ----------
    x : str
        column of `data` indicating the group (e.g. experimental group)
        
    y :  str
        column of `data` indicating the cell cluster
        
    data : pandas.DataFrame
       `pandas.DataFrame` object with 
        
    show : bool
        (Default value = True)
    output : NoneType or str
        (Default value = None)
    fig : matplotlib.figure.Figure or None
        (Default value = None)
    cax : matplotlib.axes._subplots.AxesSubplot or None
        (Default value = None)
    ax : matplotlib.axes._subplots.AxesSubplot
        (Default value = None)
    side_annotation : bool
        (Default value = True)
    heatmap_shading_key : str
        (Default value = 'ClusterFraction')
    annotation_key : str
        (Default value = 'Counts')
    annotation_fmt : str
        (Default value = '.5g')

    Returns
    -------
    None.  

    See Also
    --------

    panopticon.analysis.get_cluster_enrichment_dataframes : routine for generating dataframes used in this visualization
    scipy.stats.fisher_exact : Fisher's exact test
    panopticon.utilities.phi_coefficient : phi coefficient

    
    """
    from panopticon.analysis import get_cluster_enrichment_dataframes
    import matplotlib.pyplot as plt
    import seaborn as sns

    cluster_enrichment_dataframes = get_cluster_enrichment_dataframes(
        x, y, data)
    if fig is None and cax is None and ax is None:
        fig, (cax,
              ax) = plt.subplots(nrows=2,
                                 figsize=figsize,
                                 gridspec_kw={"height_ratios": [0.025, 1]})
    elif fig is not None and cax is not None and ax is not None:
        pass
    else:
        raise Exception(
            'Either none of "fig", "ax", and "cax" must be None, or all must be'
        )

    if len(data[x].unique()) != 2 and side_annotation:
        raise Exception(
            "Side annotation only supported for data with two unique groups")

    if (heatmap_shading_key not in cluster_enrichment_dataframes._fields):
        raise Exception("heatmap_shading key must be one of {}".format(
            cluster_enrichment_dataframes._fields))
    if (annotation_key not in cluster_enrichment_dataframes._fields):
        raise Exception("annotation_key key must be one of {}".format(
            cluster_enrichment_dataframes._fields))

    sns.heatmap(cluster_enrichment_dataframes._asdict()[heatmap_shading_key],
                cmap='Blues',
                annot=cluster_enrichment_dataframes._asdict()[annotation_key],
                cbar=False,
                vmin=0,
                vmax=1,
                fmt=annotation_fmt)
    fig.colorbar(ax.get_children()[0],
                 cax=cax,
                 orientation="horizontal",
                 label='proportion of cells (in-cluster)')
    #      annotation_clip=False)    #scipy.stats.chisquare
    cax.xaxis.tick_top()
    cax.xaxis.set_label_position('top')
    if side_annotation:
        for (i, fep_row), (j, phi_row) in zip(
                cluster_enrichment_dataframes.FishersExactP.reset_index(
                    drop=True).iterrows(),
                cluster_enrichment_dataframes.PhiCoefficient.reset_index(
                    drop=True).iterrows()):
            phi = phi_row.values[0]
            pval = fep_row.values[0]
            ax.annotate('p = {0:.2e}, '.format(pval) + r'$\phi$' +
                        ' = {0:.4f}'.format(phi),
                        xy=(2, i + .5),
                        xytext=(2.3, i + .5),
                        ha='left',
                        va='center',
                        bbox=dict(boxstyle='square', fc='white', ec='black'),
                        arrowprops=dict(
                            arrowstyle='-[, widthB={}, lengthB=0'.format(1),
                            lw=2.0,
                            color='k'),
                        annotation_clip=False)  #scipy.stats.chisquare
    if output is not None:
        plt.tight_layout()
        if output.endswith('.png'):
            plt.savefig(output, dpi=600)
        else:
            plt.savefig(output)
    plt.show()


def expand_value_counts(df, counts_col):
    return df.iloc[np.hstack([[i] * df[counts_col].values[i]
                              for i in range(len(df_filtered))])]
