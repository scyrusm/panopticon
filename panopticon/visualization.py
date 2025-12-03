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


#    plt.legend(ncol=len(subclusters) // 14 + 1, bbox_to_anchor=(1.05, 0.95))
#    if plot_output is not None:
#        plt.savefig(plot_output, bbox_inches='tight')
#    else:
#        plt.show()
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
    return hmdf, allgenes
    if len(all_gene_common_names) > 0:
        hmdf.index = all_gene_common_names
    else:
        hmdf.index = allgenes
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

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def legend_without_duplicate_labels(ax, loc=(1.05, .8)):
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
    ax.legend(*zip(*unique), title='', loc=loc)


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
                annotate_pvalue_vs_all=False,
                annotate_effect_size_vs_all=False,
                annotate_n=False,
                annotate_pvalue_vs_all_fmt_str='p: {0:.2f}',
                annotate_effect_size_vs_all_fmt_str='es: {0:.2f}',
                annotate_n_fmt_str='n: {}',
                annotate_hue_single_line=False,
                paired_hue_matching_col=None,
                effect_size='cohensd',
                pvalue='mannwhitney',
                custom_annotation_dict={},
                custom_annotation_fontsize=6,
                pairing_column=None,
                only_plot_hue_pairing=False,
                pairing_line_alpha=0.1,
                censor_inf_in_statistics=False):
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
    if dodge:
        sns.violinplot(data=data,
                       x=x,
                       hue=hue,
                       hue_order=hue_order,
                       y=y,
                       split=split,
                       scale='width',
                       cut=0,
                       alpha=alpha,
                       ax=ax,
                       width=.8,
                       **violinplot_kwargs)
    else:
        sns.violinplot(data=data,
                       x=x,
                       y=y,
                       split=split,
                       scale='width',
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

    from panopticon.visualization import legend_without_duplicate_labels
    legend_without_duplicate_labels(ax)
    if annotate_hue_pvalues or annotate_hue_effect_size or annotate_hue_n and len(
            custom_annotation_dict.keys()) == 0:
        if len(data[hue].unique()) != 2:
            raise Exception(
                "hue must be a categorical variable with 2 unique values")
        from scipy.stats import mannwhitneyu, ttest_ind, ttest_rel, wilcoxon
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
            if pvalue in ['wilcoxon_signed_rank', 'ttest_rel']:
                if paired_hue_matching_col is None:
                    raise Exception(
                        "paired_hue_matching_col must be specified")
                a = data[(data[hue] == hue1)
                         & (data[category_col].astype(str)
                            == category)]  #[continuous_col].values
                b = data[(data[hue] == hue2)
                         & (data[category_col].astype(str) == category)]  #
                ab = pd.merge(a,
                              b,
                              on=paired_hue_matching_col,
                              suffixes=('_a', '_b'))
                a = ab[continuous_col + '_a'].values
                b = ab[continuous_col + '_b'].values
                #  [continuous_col].values
                # print(data[(data[hue] == hue1)
                #         & (data[category_col] == category)])
                #print(category, category_col, continuous_col)
                if pvalue == 'wilcoxon_signed_rank':
                    f = wilcoxon
                elif pvalue == 'ttest_rel':
                    f = ttest_rel
                pval = f(a, b).pvalue

            elif pvalue in ['mannwhitney', 'ttest']:
                a = data[(data[hue] == hue1)
                         & (data[category_col].astype(str) == category
                            )][continuous_col].values
                b = data[(data[hue] == hue2)
                         & (data[category_col].astype(str) == category
                            )][continuous_col].values
                a = np.array([x for x in a if not np.isnan(x)])
                b = np.array([x for x in b if not np.isnan(x)])
                if censor_inf_in_statistics:
                    a = np.array([x for x in a if not np.isinf(x)])
                    b = np.array([x for x in b if not np.isinf(x)])

                if len(a) == 0 or len(b) == 0:
                    pval = np.nan
                else:
                    mw = mannwhitneyu(a, b, alternative='two-sided')
                    tt = ttest_ind(a, b)
                    if pvalue == 'mannwhitney':
                        pval = mw.pvalue
                    elif pvalue == 'ttest':
                        pval = tt[1]
            else:
                raise Exception(
                    "`pvalue` must be either \'mannwhitney\' or \'ttest\' or \'wilcoxon_signed_rank\' or \'ttest_rel\'"
                )
            if len(a) == 0 or len(b) == 0:
                es = np.nan
            else:
                if effect_size == 'cohensd':
                    from panopticon.utilities import cohensd
                    es = cohensd(a, b)
                elif effect_size == 'cles':
                    es = mw.statistic / len(a) / len(b)
                else:
                    raise Exception(
                        "effect_size must be either \'cohensd\' or \'cles\'")
            annotation_string = ''
            if len(a) == 0 or len(b) == 0:
                pass
            elif vertical_violins:
                if annotate_hue_pvalues:
                    annotation_string += annotate_hue_pvalue_fmt_str.format(
                        pval) + '\n'
                if annotate_hue_effect_size:
                    annotation_string += annotate_hue_effect_size_fmt_str.format(
                        es) + '\n'
                if annotate_hue_n:
                    annotation_string += annotate_hue_n_fmt_str.format(
                        len(a), len(b)) + '\n'
                if annotate_hue_single_line:
                    anno_x, anno_y = ticklabel.get_position(
                    )[0], data[continuous_col].max()
                else:
                    anno_x, anno_y = ticklabel.get_position()[0], np.max(
                        np.hstack((a, b)))

                ax.annotate(annotation_string, (anno_x, anno_y),
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
                if annotate_hue_single_line:
                    anno_x, annoy_y = data[continuous_col].max(
                    ), ticklabel.get_position()[1]
                else:
                    anno_x, anno_y = np.max(np.hstack(
                        (a, b))), ticklabel.get_position()[1]

                ax.annotate(annotation_string, (anno_x, anno_y),
                            ha='left',
                            va='center')
    if annotate_pvalue_vs_all or annotate_effect_size_vs_all or annotate_n and len(
            custom_annotation_dict.keys()) == 0:
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
            if pvalue in ['mannwhitney', 'ttest']:
                a = data[data[category_col] == category][continuous_col].values
                b = data[data[category_col] != category][continuous_col].values
                a = np.array([x for x in a if not np.isnan(x)])
                b = np.array([x for x in b if not np.isnan(x)])
                if censor_inf_in_statistics:
                    a = np.array([x for x in a if not np.isinf(x)])
                    b = np.array([x for x in b if not np.isinf(x)])
                mw = mannwhitneyu(a, b, alternative='two-sided')
                tt = ttest_ind(a, b)
                if pvalue == 'mannwhitney':
                    pval = mw.pvalue
                elif pvalue == 'ttest':
                    pval = tt[1]
            else:
                raise Exception(
                    "`pvalue` must be either \'mannwhitney\' or \'ttest\' ")

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
                if annotate_pvalue_vs_all:
                    annotation_string += annotate_pvalue_vs_all_fmt_str.format(
                        pval) + '\n'
                if annotate_effect_size_vs_all:
                    annotation_string += annotate_effect_size_vs_all_fmt_str.format(
                        es) + '\n'
                if annotate_n:
                    annotation_string += annotate_n_fmt_str.format(
                        len(a)) + '\n'
                ax.annotate(annotation_string,
                            (ticklabel.get_position()[0], np.max(a)),
                            ha='center',
                            va='bottom')
            else:
                if annotate_pvalue_vs_all:
                    annotation_string += ' ' + annotate_pvalue_vs_all_fmt_str.format(
                        pval)
                if annotate_effect_size_vs_all:
                    annotation_string += '\n ' + annotate_effect_size_vs_all_fmt_str.format(
                        es)
                if annotate_n:
                    annotation_string += '\n' + annotate_n_fmt_str.format(
                        len(a))

                ax.annotate(annotation_string, (
                    np.max(a),
                    ticklabel.get_position()[1],
                ),
                            ha='left',
                            va='center',
                            annotation_clip=False)
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
            annotation_pos = np.max(
                data[data[category_col] == category][continuous_col].values)

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
    if pairing_column is not None:
        xoffsets = []
        yoffsets = []
        if hue is None:
            groups = [x.get_text() for x in ax.get_xticklabels()]
        else:
            groups = list(
                np.hstack([[x.get_text()] * 2 for x in ax.get_xticklabels()]))
            group_notdoubled = [x.get_text() for x in ax.get_xticklabels()]
        n_empty_collections = 0
        for collection in ax.collections:
            if collection.get_offsets().shape[0] <= 1:
                n_empty_collections += 1
        for collection in ax.collections[n_empty_collections::]:
            if type(collection.get_offsets()) == np.ma.core.MaskedArray:

                xoffsets.append(collection.get_offsets().data[:, 0])
                yoffsets.append(collection.get_offsets().data[:, 1])
            elif type(collection.get_offsets()) == np.ndarray:
                xoffsets.append(collection.get_offsets()[:, 0])
                yoffsets.append(collection.get_offsets()[:, 1])
            else:
                raise Exception('Collection offsets unfamiliar type')


#        print(xoffsets)
        xoffsets = np.vstack(xoffsets)
        yoffsets = np.vstack(yoffsets)
        offset_df = pd.DataFrame(groups, columns=['group'])
        for i in range(xoffsets.shape[1]):
            offset_df['{}_x'.format(i)] = xoffsets[:, i]
            offset_df['{}_y'.format(i)] = yoffsets[:, i]
            if (hue is not None) and only_plot_hue_pairing:
                for j in range(0, xoffsets.shape[0], 2):
                    ax.plot([xoffsets[j, i], xoffsets[j + 1, i]],
                            [yoffsets[j, i], yoffsets[j + 1, i]],
                            color='k',
                            alpha=pairing_line_alpha,
                            ls='--')
            else:
                ax.plot(xoffsets[:, i],
                        yoffsets[:, i],
                        color='k',
                        alpha=pairing_line_alpha,
                        ls='--')
        replicate_matches = []
        for col in [y for y in offset_df.columns if y.endswith('_y')]:
            for replicate in data[pairing_column].unique():
                if hue is None:
                    if offset_df.set_index('group')[col].equals(
                            data[data[pairing_column] == replicate].set_index(
                                x).loc[groups][y]):
                        replicate_matches.append(replicate)
                else:
                    if offset_df.set_index('group')[col].reset_index(
                            drop=True).equals(data[
                                data[pairing_column] == replicate].set_index(
                                    x).loc[group_notdoubled][y].reset_index(
                                        drop=True)):

                        replicate_matches.append(replicate)

        if len(np.unique(replicate_matches)) != len(
                np.unique(data[pairing_column])):
            raise Exception(
                "Problem with pairing column--it is possible that pairing annotation was done incorrectly"
            )
    return ax


def volcano(diffex,
            ax=None,
            gene_column='gene',
            pval_column='pvalue',
            effect_size_col='CommonLanguageEffectSize',
            left_name='',
            right_name='',
            genemarklist=[],
            neglogpval_importance_threshold=5,
            title='',
            output=None,
            positions=None,
            show=True,
            gene_label_offset_scale=1,
            side_annotation_gene_label_offset_scale=1,
            no_effect_line=0,
            counterscale=1,
            lcounter_init=0,
            rcounter_init=0,
            verbose=False,
            draggable_annotations=False,
            gene_position_dict_for_side_annotations={},
            boldgenes=[],
            specialcolorgenelist=[],
            specialcolor='r',
            value_for_significance_truncation=0):
    """

    Parameters
    ----------
    diffex :
        param ax:  (Default value = None)
    gene_column :
        Default value = 'gene')
    pval_column :
        Default value = 'pvalue')
    left_name : Genes toward the left in the volcano plot, which are upregulated in group2 of 'diffex'
        Default value = '')
    right_name : Genes toward the right in the volcano plot, which are upregulated in group1 of 'diffex'
        Default value = '')
    genemarklist :
        Default value = [])
    effect_size_importance_threshold :
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

    #    matplotlib.rcParams['axes.linewidth'] = 3
    if ax is None:
        fig, ax = plt.subplots(figsize=(
            5,
            5,
        ))
    diffex[pval_column] = diffex[pval_column].apply(lambda x: x if x>value_for_significance_truncation else value_for_significance_truncation)
    neglogpvalues = -np.log(diffex[pval_column].values) / np.log(10)
    effect_size = diffex[effect_size_col].values
    important_mask = (neglogpvalues > neglogpval_importance_threshold)
    if len(specialcolorgenelist)>0:
        color=diffex[gene_column].isin(specialcolorgenelist).apply(lambda x: specialcolor if x else 'b').values
    else:
        color='b'
    ax.scatter(effect_size[important_mask],
               neglogpvalues[important_mask],
               alpha=1,
               marker='.',
               s=3,
               c=color[important_mask])
    ax.scatter(effect_size[~important_mask],
               neglogpvalues[~important_mask],
               alpha=.1,
               marker='.',
               s=2,
               c=color[~important_mask])
    maxx = np.nanmax(effect_size)
    minx = np.nanmin(effect_size)

    maxy = np.nanmax(neglogpvalues)
    left_edge = minx + .4 * (minx - maxx)
    right_edge = maxx + .4 * (maxx - minx)
    top_edge = maxy * 1.01
    bottom_edge = 0
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
            np.vstack((effect_size[annomask], neglogpvalues[annomask])).T)
        positions = []
        for i in range(1, len(distances)):
            mindist = np.min(distances[i, 0:i])
            if mindist < 0.2:
                positions.append('l')
            else:
                positions.append('r')

    def position_to_xytext_habt_va(position, effect_size, negpval, maxx, maxy,
                                   gene_label_offset_scale):
        if position in ['br', 'bl', 'tr', 'tl']:
            if position[1] == 'r':
                habt = 'right'
            elif position[1] == 'l':
                habt = 'left'
            else:
                raise Exception(
                    "\'{}\':  invalid position character selection".format(
                        position))
            position = position[0]
        else:
            habt = 'center'
        if position == 'b':
            xytext = (effect_size,
                      negpval + .03 * maxy * gene_label_offset_scale)
            va = 'bottom'
        elif position == 't':
            xytext = (effect_size,
                      negpval - .03 * maxy * gene_label_offset_scale)
            va = 'top'
        elif position == 'l':
            xytext = (effect_size + .03 * maxx * gene_label_offset_scale,
                      negpval)
            va = 'center'
        elif position == 'r':
            xytext = (effect_size - .03 * maxx * gene_label_offset_scale,
                      negpval)
            va = 'center'
        else:
            raise Exception(
                "\'{}\':  invalid position character selection".format(
                    position))
        return xytext, habt, va

    negpvals = []
    for gene in genemarklist:
        genedf = diffex[diffex[gene_column] == gene]
        negpval = -np.log(genedf.iloc[0][pval_column]) / np.log(10)
        negpvals.append(negpval)
    genemarklist = list(np.array(genemarklist)[np.argsort(negpvals)][::-1])
    if positions != 'side':
        if type(positions) == dict:
            positions = [positions[key] for key in genemarklist]
        for gene, position in zip(
                genemarklist,
                positions,
        ):
            if verbose:
                print(gene, position)
            genedf = diffex[diffex[gene_column] == gene]
            negpval = -np.log(genedf.iloc[0][pval_column]) / np.log(10)
            effect_size = genedf.iloc[0][effect_size_col]
            #if len(specialcolorgenelist)>0:
            #    color=diffex[gene_column].isin(specialcolorgenelist).apply(lambda x: specialcolor if x else 'k').values
            #else:
            #color='k'
            ax.scatter(effect_size, negpval, marker='.', color='k')
            xytext, habt, va = position_to_xytext_habt_va(
                position, effect_size, negpval, maxx, maxy,
                gene_label_offset_scale)
            anno = ax.annotate(
                gene, (effect_size, negpval),
                xytext,
                va=va,
                ha=habt,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])
            if draggable_annotations:
                anno.draggable()
    else:
        lcounter = lcounter_init
        rcounter = rcounter_init
        if type(positions) == dict:
            positions = [positions[key] for key in genemarklist]
        for gene in genemarklist:
            genedf = diffex[diffex[gene_column] == gene]
            negpval = -np.log(genedf.iloc[0][pval_column]) / np.log(10)
            effect_size = genedf.iloc[0][effect_size_col]
            fontweight='normal'
            if gene in boldgenes:
                fontweight='bold'
            if gene in gene_position_dict_for_side_annotations.keys():
                position = gene_position_dict_for_side_annotations[gene]
                xytext, habt, va = position_to_xytext_habt_va(
                    position, effect_size, negpval, maxx, maxy,
                    gene_label_offset_scale)
                anno = ax.annotate(gene, (effect_size, negpval),
                                   xytext,
                                   va=va,
                                   ha=habt,
                                   fontweight=fontweight,
                                   path_effects=[
                                       pe.withStroke(linewidth=2,
                                                     foreground="white")
                                   ])
            else:

                if effect_size < no_effect_line:
                    xytext = (
                        left_edge +
                        .03 * maxx * side_annotation_gene_label_offset_scale,
                        top_edge - lcounter)
                    habt = 'right'
                    va = 'center'
                    lcounter += counterscale
                else:
                    xytext = (
                        right_edge -
                        .03 * maxx * side_annotation_gene_label_offset_scale,
                        top_edge - rcounter)
                    habt = 'left'
                    va = 'center'
                    rcounter += counterscale
                anno = ax.annotate(gene, (effect_size, negpval),
                                   xytext=xytext,
                                   va=va,
                                   ha=habt,
                                   fontweight=fontweight,
                                   arrowprops=dict(facecolor='black',
                                                   width=0.1,
                                                   headwidth=0,
                                                   alpha=0.25),
                                   path_effects=[
                                       pe.withStroke(linewidth=2,
                                                     foreground="white")
                                   ])
            if draggable_annotations:
                anno.draggable()
            ax.scatter(effect_size, negpval, marker='.', color='k')

    ax.axvline(no_effect_line, ls='--', color='k', alpha=0.25)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)
    ax.set_xlabel(effect_size_col + '\n' + '$\leftarrow$' + left_name + '\n' +
                  right_name + r'$\rightarrow$',
                  fontsize=14)

    ax.set_ylabel('-log' + r'${}_{10}$' + '(p-value)', fontsize=14)
    ax.set_xlim([left_edge, right_edge])
    ax.set_ylim([bottom_edge, top_edge])
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
                    annotate_simpson=False,
                    weights=None,
                    colorkey_col=None):
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
    if colorkey_col is not None and (stack_order != 'agnostic'
                                     or piechart == True or weights is None):
        raise Exception(
            "Color key currently implemented only for agnostic stack order stacked bar plot, with weights"
        )

    all_heights = []
    if colorkey_col is not None:
        all_colorkeys = []

    if hue is None:
        if stack_order == 'agnostic':
            if weights is None:
                grouped_data = data.groupby(x)[y].value_counts()
            else:
                if data[weights].dtypes != int:
                    raise Exception(
                        "dtype of column \'{}\' must be int".format(weights))
                if colorkey_col is None:
                    grouped_data = data.groupby(
                        [x, y])[weights].sum().sort_values(ascending=False)
                else:
                    grouped_data = data.groupby([
                        x, y, colorkey_col
                    ])[weights].sum().sort_values(ascending=False)
        elif stack_order == 'matched':
            if weights is None:
                grouped_data = data.groupby(x)[y].value_counts(sort=False)
            else:
                if data[weights].dtypes != int:
                    raise Exception(
                        "dtype of column \'{}\' must be int".format(weights))
                grouped_data = data.groupby([x, y])[weights].sum()

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
            if weights is None:
                grouped_data = data.groupby([x, hue])[y].value_counts()
            else:
                #                from IPython.core.debugger import set_trace; set_trace()
                if data[weights].dtypes != int:
                    raise Exception(
                        "dtype of column \'{}\' must be int".format(weights))
                if colorkey_col is None:
                    grouped_data = data.groupby(
                        [x, hue,
                         y])[weights].sum().sort_values(ascending=False)
                else:
                    grouped_data = data.groupby([
                        x, hue, y, colorkey_col
                    ])[weights].sum().sort_values(ascending=False)

        elif stack_order == 'matched':
            if weights is None:
                grouped_data = data.groupby([x,
                                             hue])[y].value_counts(sort=False)
            else:
                if data[weights].dtypes != int:
                    raise Exception(
                        "dtype of column \'{}\' must be int".format(weights))
                grouped_data = data.groupby([x, hue, y])[weights].sum()

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
#    from IPython.core.debugger import set_trace; set_trace()
    for grouping in groupings:
        heights = []
        for n in grouped_data[grouping].values:
            heights.append(n)
        heights = heights + [0] * (total - len(heights))
        heights = np.array(heights)[::-1]
        all_heights.append(heights)

        if colorkey_col is not None:
            colors = []
            for index in grouped_data[grouping].index.values:
                colors.append(index[-1])
            all_colorkeys.append(colors)

    all_heights = np.vstack(all_heights)
    all_heights = all_heights[:, ::-1]

    if colorkey_col is not None:
        #        from IPython.core.debugger import set_trace; set_trace()
        all_colorkeys_padded = []
        max_colorkey_length = np.max([len(x) for x in all_colorkeys])
        for colorkey in all_colorkeys:
            all_colorkeys_padded.append(
                colorkey + (max_colorkey_length - len(colorkey)) *
                [colorkey[-1]])  # padds the last value as a dummy
#        from IPython.core.debugger import set_trace; set_trace()
        all_colorkeys = np.vstack(all_colorkeys_padded)


#        all_colorkeys = all_colorkeys[:,::-1]

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
            title = str(groupings[i])
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
                if colorkey_col is not None:
                    #                    return np.unique(all_colorkeys)
                    colorkey2color = {
                        key: color_palette[i % len(color_palette)]
                        for i, key in enumerate(np.unique(all_colorkeys))
                    }
                    subax.bar(
                        ind,
                        all_heights[:, i],
                        0.8,
                        bottom=bottoms,
                        color=[colorkey2color[x] for x in all_colorkeys[:, i]],
                        edgecolor='k',
                        lw=0.1,
                        label=label)
                    subax.set_xticks(ind)
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
        plt.legend(bbox_to_anchor=(1, 1), reverse=True)
    elif colorkey_col is not None:
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor=colorkey2color[colorkey],
                  edgecolor='k',
                  label=colorkey) for colorkey in colorkey2color.keys()
        ]
        plt.legend(handles=legend_elements,
                   bbox_to_anchor=(1, 1),
                   title=colorkey_col,
                   reverse=True)

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


def plot_differential_density(x,
                              y,
                              mask1,
                              mask2,
                              ax=None,
                              cmap=plt.cm.RdBu_r,
                              buffer_factor=0):
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
    xbuffer = (x.max() - x.min()) * buffer_factor
    ybuffer = (y.max() - y.min()) * buffer_factor
    xmin, xmax = (x.min() - xbuffer, x.max() + xbuffer)
    ymin, ymax = (y.min() - ybuffer, y.max() + ybuffer)

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
                               figsize=(5, 5),
                               weights=None,
                               sort_by_phi_col=None):
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
        x, y, data, weights=weights, sort_by_phi_col=sort_by_phi_col)
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
    if (annotation_key
            not in cluster_enrichment_dataframes._fields) and (annotation_key
                                                               is not None):
        raise Exception("annotation_key key must be one of {}".format(
            cluster_enrichment_dataframes._fields))
    if annotation_key is not None:
        annot = cluster_enrichment_dataframes._asdict()[annotation_key]
    else:
        annot = None
    sns.heatmap(cluster_enrichment_dataframes._asdict()[heatmap_shading_key],
                cmap='Blues',
                annot=annot,
                cbar=False,
                vmin=0,
                vmax=1,
                fmt=annotation_fmt)
    if heatmap_shading_key == 'FractionOfCluster':
        cbar_label = 'proportion of cells (row-normalized)'
    elif heatmap_shading_key == 'FractionOfGroup':
        cbar_label = 'proportion of cells (column-normalized)'
    fig.colorbar(ax.get_children()[0],
                 cax=cax,
                 orientation="horizontal",
                 label=cbar_label)
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
    if show:
        plt.show()
    return fig


def expand_value_counts(df, counts_col, scale=500):
    return df.iloc[np.hstack([[i] * df[counts_col].values[i]
                              for i in range(len(df_filtered))])]

    import matplotlib.pyplot as plt
    fig, (ax, cax) = plt.subplots(1,
                                  2,
                                  figsize=(4, 20),
                                  gridspec_kw={'width_ratios': [20, 1]})
    key2x = {key: x for key, x in zip(df.index.values, range(df.shape[0]))}
    marker2y = {
        marker: x
        for marker, x in zip(df.columns, range(len(df.columns)))
    }
    exprs = []
    for col in df.columns:
        for key in [key for key in df.index.values]:
            exprs.append(df.loc[key][col])

    for col in df.columns:
        exprs = []
        for key in [key for key in df.index.values]:
            exprs.append(df.loc[key][col])
        for key in [key for key in df.index.values]:
            expr = df.loc[key][col]
            sc = ax.scatter(key2x[key],
                            marker2y[col],
                            s=(expr - np.min(exprs)) * scale,
                            cmap='coolwarm',
                            c=expr,
                            vmin=np.min(exprs),
                            vmax=np.max(exprs))

    ax.set_yticklabels([x.split('_With')[0] for x in df.columns])
    ax.set_yticks(range(len(df.columns)))

    ax.set_xticks(range(len(df.index.values)))
    ax.set_xticklabels(list(range(len(df.index.values))))
    ax.set_xlabel('Cluster')
    ax.set_ylim([-0.5, len(df.columns)])
    cbar = plt.colorbar(
        sc,
        aspect=30,
        cax=cax,
    )
    cbar.set_label('module score (min-max normalized)',
                   rotation=90,
                   fontsize=17)
    yticklabels = cbar.ax.get_yticklabels()
    cax_ylims = cax.get_ylim()
    cbar.ax.set_yticks(cax_ylims)
    cbar.ax.set_yticklabels([0, 1])
    plt.tight_layout()


def plot_differential_expression_barplot(
        diffex,
        gene_col='gene',
        effect_size='CommonLanguageEffectSize',
        n_genes=10,
        gene_blacklist=[],
        fig=None,
        ax=None,
        output=None,
        show=True):
    import seaborn as sns
    import matplotlib.pyplot as plt
    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=(4, 4))
    sns.barplot(data=pd.concat([
        diffex[~diffex[gene_col].isin(gene_blacklist)].sort_values(
            effect_size, ascending=False).head(n_genes),
        diffex[~diffex[gene_col].isin(gene_blacklist)].sort_values(
            effect_size, ascending=True).head(n_genes)[::-1]
    ]),
                x=effect_size,
                y=gene_col,
                color='k',
                ax=ax)
    if output is not None:
        plt.tight_layout()
        plt.savefig(output)
    if show:
        plt.show()


def plot_dot_plot(loom,
                  x_column_attribute=None,
                  y_genes=None,
                  y_column_attributes=None,
                  x_column_blacklist=[],
                  layername='log2(TP10k+1)',
                  gene_ra_name='gene',
                  scale=100,
                  figsize=(4, 15),
                  dot_size_mode='FractionNonzero',
                  cbar_label='expression/module score (min-max normalized)',
                  orientation='horizontal',
                  cmap='coolwarm',
                  cbar_range='zero_centered',
                  x_column_attribute_sortkey=None,
                  legend_bbox_to_anchor=(-.05, 1),
                  z_score=False,
                  minmax_normalize=False):
    if x_column_attribute not in loom.ca.keys():
        raise Exception(
            "x_column_attribute not a column attribute of loomfile")
    if y_genes is None and y_column_attributes is None:
        raise Exception(
            "One of genes or y_column_attributes must be an iterable")
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    import matplotlib.pyplot as plt
    if dot_size_mode not in ['FractionNonzero', 'MinMaxNormalization']:
        raise Exception(
            'dot_size_mode must be one of "FractionNonzero" or "MinMaxNormalization"'
        )

    df = pd.DataFrame(loom.ca[x_column_attribute],
                      columns=[x_column_attribute])
    if y_genes is not None:
        for gene in y_genes:
            if gene not in loom.ra[gene_ra_name]:
                raise Exception("{} not in loom.ra[{}]".format(
                    gene, gene_ra_name))
            else:
                igene = np.where(loom.ra[gene_ra_name] == gene)[0]
            df[gene] = loom[layername][igene, :].sum(axis=0)
            if dot_size_mode == 'FractionNonzero':
                if gene + '_fraction' in df.columns:
                    raise Exception(
                        '"{}" already one of the designated column attributes'.
                        format(gene + '_fraction'))
                df[gene + '_fraction'] = df[gene] > 0
    if y_column_attributes is not None:
        for y_column_attribute in y_column_attributes:
            if y_column_attribute not in loom.ca.keys():
                raise Exception(
                    "{} not in loom.ca.keys()".format(y_column_attribute))
            df[y_column_attribute] = loom.ca[y_column_attribute]
            if dot_size_mode == 'FractionNonzero':
                if y_column_attribute + '_fraction' in df.columns:
                    raise Exception(
                        '"{}" already one of the designated column attributes'.
                        format(y_column_attribute + '_fraction'))
                df[y_column_attribute +
                   '_fraction'] = df[y_column_attribute] > 0
    df = df[~df[x_column_attribute].isin(x_column_blacklist)]
    df = df.groupby(x_column_attribute).mean()
    df = df.sort_index(ascending=False)
    if x_column_attribute_sortkey is not None:
        df = df.sort_index(key=np.vectorize(x_column_attribute_sortkey),
                           ascending=True)

    if z_score and minmax_normalize:
        raise Exception("Only one of z_score or minmax_normalize may be true")
    if z_score:
        for col in df.columns:
            df[col] = (df[col] - df[col].mean()) / df[col].std()
    elif minmax_normalize:
        for col in df.columns:
            df[col] = (df[col] - df[col].min()) / (df[col].max()-df[col].min())
    fig, (ax, cax) = plt.subplots(
        1,
        2,
        figsize=figsize,
        gridspec_kw={
            'width_ratios': [100, 1],
            'wspace': 0.05
        },
    )
    #if orientation

    if dot_size_mode == 'FractionNonzero':

        relevant_columns = [
            x for x in df.columns if not x.endswith('_fraction')
        ]
        #if orientation=='vertical':
        key2x = {key: x for key, x in zip(df.index.values, range(df.shape[0]))}
        marker2y = {
            marker: x
            for marker, x in zip(relevant_columns, range(len(
                relevant_columns)))
        }

        exprs = [
        ]  # this loops over all expression values in order to create a consistent color
        for col in tqdm(relevant_columns):
            for key in [key for key in df.index.values]:
                exprs.append(df.loc[key][col])
        if cbar_range == 'zero_centered':
            vmin = np.min([np.min(exprs), -np.max(exprs)])
            vmax = np.max([-np.min(exprs), np.max(exprs)])
        elif cbar_range == 'zero_based':
            vmin = 0
            vmax = np.max(exprs)
        elif cbar_range == 'data_range':
            vmin = np.min(exprs)
            vmax = np.max(exprs)
        else:
            raise Exception(
                'cbar_range must be in ["zero_centered","zero_based","data_range"]'
            )
        for col in relevant_columns:
            for key in [key for key in df.index.values]:
                expr = df.loc[key][col]
                frac = df.loc[key][col + '_fraction']
                if orientation == 'vertical':
                    x = key2x[key]
                    y = marker2y[col]
                elif orientation == 'horizontal':
                    y = key2x[key]
                    x = marker2y[col]
                sc = ax.scatter(x,
                                y,
                                s=frac * scale + 1e-9,
                                cmap=cmap,
                                c=expr,
                                vmin=vmin,
                                vmax=vmax,
                                edgecolors='black',
                                linewidth=0.1)
    elif dot_size_mode == 'MinMaxNormalization':
        relevant_columns = df.columns
        key2x = {key: x for key, x in zip(df.index.values, range(df.shape[0]))}
        marker2y = {
            marker: x
            for marker, x in zip(relevant_columns, range(len(
                relevant_columns)))
        }

        exprs = [
        ]  # this loops over all expression values in order to create a consistent color
        for col in tqdm(relevant_columns):
            for key in [key for key in df.index.values]:
                exprs.append(df.loc[key][col])
        if cbar_range == 'zero_centered':
            vmin = np.min([np.min(exprs), -np.max(exprs)])
            vmax = np.max([-np.min(exprs), np.max(exprs)])
        elif cbar_range == 'zero_based':
            vmin = 0
            vmax = np.max(exprs)
        elif cbar_range == 'data_range':
            vmin = np.min(exprs)
            vmax = np.max(exprs)
        else:
            raise Exception(
                'cbar_range must be in ["zero_centered","zero_based","data_range"]'
            )
        for col in relevant_columns:
            for key in [key for key in df.index.values]:

                expr = df.loc[key][col]
                if orientation == 'vertical':
                    x = key2x[key]
                    y = marker2y[col]
                elif orientation == 'horizontal':
                    y = key2x[key]
                    x = marker2y[col]
                sc = ax.scatter(x,
                                y,
                                s=(expr - np.min(exprs)) * scale + 1e-9,
                                cmap=cmap,
                                c=expr,
                                vmin=vmin,
                                vmax=vmax,
                                edgecolors='black',
                                linewidth=0.1)
    if orientation == 'vertical':
        ax.set_yticklabels([x.split('_With')[0] for x in relevant_columns])
        ax.set_yticks(range(len(relevant_columns)))

        ax.set_xticks(np.array(range(len(df.index.values))))
        ax.set_xticklabels(list(df.index.values))
        ax.set_xlabel('Cluster')
        ax.set_ylim([-0.5, len(relevant_columns) - 0.5])
        ax.set_xlim([-0.5, len(df.index) - 0.5])
        #for i in np.arange(-0.5,len(markers)-0.5,topn)[1::]:
        #    plt.axhline(ls='--',y=i,color='k')
        #if flag:
        cbar = plt.colorbar(
            sc,
            aspect=30,
            cax=cax,
        )
        cbar.set_label(cbar_label, rotation=90, fontsize=17)
        yticklabels = cbar.ax.get_yticklabels()
        cax_ylims = cax.get_ylim()
        cbar.ax.set_yticks(cax_ylims)
        cbar.ax.set_yticklabels([0, 1])
    elif orientation == 'horizontal':
        major_tick_columns = relevant_columns[0::2]
        minor_tick_columns = relevant_columns[1::2]
        major_tick_positions = list(np.arange(0,len(relevant_columns),2))
        minor_tick_positions = list(np.arange(1,len(relevant_columns),2))
        ax.set_xlim([-0.5, len(relevant_columns) - 0.5])
        ax.set_ylim([-0.5, len(df.index) - 0.5])
        ax.set_xticks(major_tick_positions)
        ax.set_xticks(minor_tick_positions,minor=True)
        ax.set_xticklabels([x.split('_With')[0] for x in major_tick_columns])
        ax.set_xticklabels([x.split('_With')[0] for x in minor_tick_columns],minor=True)

        ax.set_yticks(np.array(range(len(df.index.values))))
        ax.set_yticklabels(list(df.index.values))
        ax.set_ylabel('Cluster')

        #for i in np.arange(-0.5,len(markers)-0.5,topn)[1::]:
        #    plt.axhline(ls='--',y=i,color='k')
        #if flag:
        cbar = plt.colorbar(
            sc,
            aspect=30,
            cax=cax,
        )
        cbar.set_label(cbar_label, rotation=0, fontsize=14, va='center')
        yticklabels = cbar.ax.get_xticklabels()
        cax_ylims = cax.get_xlim()
        cbar.ax.set_xticks(cax_ylims)
        cbar.ax.set_xticklabels([])
        from matplotlib import ticker

#        print(relevant_columns[1:(len(relevant_columns)+1):2])
#        print(ax
#        ax.xaxis.set_minor_formatter(
#            ticker.FixedFormatter(
#                                  ['']*0+relevant_columns[1:(len(relevant_columns) +
#                                                      1):2]+['']*0))
#
        ax.tick_params(axis='x', which='minor', length=15)
        ax.tick_params(axis='x', which='both', color='lightgrey')
        ax.autoscale(enable=True, axis='x', tight=True)
        ax.set_xlim([-0.5, len(relevant_columns) - 0.5])
        ax.set_ylim([-0.5, len(df.index) - 0.5])
        ax.set_ylabel('')
        #plt.ylabel('y_description', loc='right')

    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D
    handles, labels = ax.get_legend_handles_labels()

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if dot_size_mode == 'FractionNonzero':
        import seaborn as sns
        points = []
        for sizereference in [0.05, 0.5, 0.9]:
            points.append(
                Line2D(
                    [0],
                    [0],
                    label='{:.0%}'.format(sizereference),
                    marker='o',
                    markersize=np.sqrt(sizereference * scale + 1e-9),
                    markeredgewidth=0.1,
                    markeredgecolor='k',
                    markerfacecolor=sns.color_palette('coolwarm')[-2],
                    linestyle='',
                ))
        # add manual symbols to auto legend
        handles.extend(points)

        ax.legend(
            handles=handles,
            bbox_to_anchor=legend_bbox_to_anchor,
            title='Fraction of cells with\npositive expression',
            borderpad=1,
            labelspacing=1.5,
        )
    elif dot_size_mode == 'MinMaxNormalization':
        point1 = Line2D([0], [0],
                        label='manual point',
                        marker='o',
                        markersize=10,
                        markeredgecolor='r',
                        markerfacecolor='k',
                        linestyle='')

        # add manual symbols to auto legend
        handles.extend([point1])

        ax.legend(handles=handles,
                  bbox_to_anchor=legend_bbox_to_anchor,
                  title='')
    plt.tight_layout()

    return fig


def _reconstruct_iterative_clustering_tree(
        loom,
        diffdict=None,
        anno_n_genes=5,
        marker_anno_delimiter='\n',
        marker_detail_level=1,
        markers_to_highlight=[],
        marker_gene_gene_ra='gene_common_name',
        marker_gene_highlight_layer='log2(TP10k+1)'):
    import numpy as np
    from panopticon.utilities import import_check
    exit_code = import_check("treelib", 'pip install treelib')
    if exit_code != 0:
        return
    from treelib import Tree
    import pandas as pd

    clustering_iteration_ca = [
        x for x in loom.ca.keys() if x.startswith('ClusteringIteration')
    ]
    n_clustering_iterations = len(clustering_iteration_ca)
    tree = Tree()
    tree.create_node("(n={})".format(loom.shape[0]), '')

    for ica in range(n_clustering_iterations):
        if 'U' in list(loom.ca['ClusteringIteration{}'.format(ica)]):
            break
        else:
            vc = pd.DataFrame(loom.ca['ClusteringIteration{}'.format(
                ica)])[0].value_counts().sort_index()
            for cluster, row in vc.items():
                if ica == 0:
                    parent = ''
                else:
                    parent = '-'.join(cluster.split('-')[0:-1])
                if len(
                        np.unique([
                            x for x in loom.ca['ClusteringIteration{}'.format(
                                ica)] if str(x).startswith(parent)
                        ])) > 1:
                    anno = "{} (n={})".format(cluster, row)
                    if diffdict is not None:
                        if type(diffdict[cluster]) != float:
                            if marker_detail_level == 1:
                                anno += '\n' + marker_anno_delimiter.join(
                                    diffdict[cluster]['GeneAlternateName'].
                                    head(anno_n_genes).values)
                            if marker_detail_level == 2:
                                #      anno +='\n'
                                for ianno, rowanno in diffdict[cluster].head(
                                        anno_n_genes).iterrows():
                                    anno += marker_anno_delimiter + '{0} (ME: {1:.2f}, FE: {2:.2f})'.format(
                                        rowanno['GeneAlternateName'],
                                        rowanno['MeanExpr1'],
                                        rowanno['FracExpr1']
                                    ) + marker_anno_delimiter
                            if marker_detail_level == 3:
                                #     anno +='\n'
                                for ianno, rowanno in diffdict[cluster].head(
                                        anno_n_genes).iterrows():
                                    anno += marker_anno_delimiter + '{0} (ME: {1:.2f}, FE: {2:.2f} vs.ME: {3:.2f}, FE: {4:.2f}'.format(
                                        rowanno['GeneAlternateName'],
                                        rowanno['MeanExpr1'],
                                        rowanno['FracExpr1'],
                                        rowanno['MeanExpr2'],
                                        rowanno['FracExpr2'])
                        for marker in markers_to_highlight:
                            cluster_level = len(str(cluster).split('-')) - 1
                            mask1 = loom.ca['ClusteringIteration{}'.format(
                                cluster_level)] == cluster
                            if marker in loom.ra[marker_gene_gene_ra]:
                                igene = np.where(loom.ra[marker_gene_gene_ra]
                                                 == marker)[0][0]
                                anno += '\n' + 'mean {0}: {1:.2f}'.format(
                                    marker, loom[marker_gene_highlight_layer][
                                        igene, :][mask1].mean())
                            #print(loom[marker_gene_highlight_layer][igene,:][mask1].mean(), loom[marker_gene_highlight_layer][igene,:][mask2].mean())
                            elif marker in loom.ca.keys():
                                if loom.ca[marker].dtype == float:
                                    anno += '\n' + 'mean {0}: {1:.2f}'.format(
                                        marker, loom.ca[marker][mask1].mean())
                                else:
                                    proportions = '{}'.format(
                                        100 *
                                        pd.DataFrame(loom.ca[marker][mask1],
                                                     columns=[marker])
                                        [marker].value_counts(normalize=True))
                                    if len(proportions.split('\n')[1:-1]) > 1:
                                        proportions = '%\n'.join(
                                            proportions.split('\n')
                                            [1:-1]) + '%'
                                        anno += '\n\n' + '{} proportions:\n{}'.format(
                                            marker, proportions)
                                    else:
                                        proportions = '%\n'.join(
                                            proportions.split('\n')
                                            [1:-1]) + '%'
                                        anno += '\n\n' + '{} :\n{}'.format(
                                            marker, proportions)

                            else:
                                raise Exception(
                                    "{} not a valid gene or ca".format(marker))
                    tree.create_node(anno, '{}'.format(cluster), parent=parent)
    return tree


def _plot_tree(t, mode='treelib', gv_filename='test.gv', shape='box'):
    if mode == 'treelib':
        t.show()
    elif mode == 'graphviz':
        t.to_graphviz(gv_filename, shape=shape)
        from panopticon.utilities import import_check
        exit_code = import_check("graphviz", 'pip install graphviz')
        if exit_code != 0:
            return
        from graphviz import Source
        s = Source.from_file(gv_filename)
        s.view()
    else:
        raise Exception("mode must be one of ('treelib','graphviz')")


def plot_iterative_clustering_tree(loom,
                                   mode='treelib',
                                   diffdict=None,
                                   graphviz_node_shape='box',
                                   marker_detail_level=1,
                                   anno_n_genes=5,
                                   markers_to_highlight=[]):
    #from panopticon.visualization import _reconstruct_iterative_clustering_tree, _plot_tree
    t = _reconstruct_iterative_clustering_tree(
        loom,
        diffdict=diffdict,
        anno_n_genes=anno_n_genes,
        marker_detail_level=marker_detail_level,
        markers_to_highlight=markers_to_highlight)
    _plot_tree(t,
               mode=mode,
               gv_filename=loom.filename + '.gv',
               shape=graphviz_node_shape)


def plot_color_coded_embedding(loom,
                               x_ca,
                               y_ca,
                               category_ca=None,
                               category_as_continuum=False,
                               use_gex_as_ca=False,
                               gene_ra='gene_common_name',
                               gex_layer='log2(TP10k+1)',
                               fig=None,
                               ax=None,
                               color_palette='colorblind',
                               legend=True,
                               on_figure_annotation=False,
                               s=2):
    import numpy as np
    import seaborn as sns
    if fig is not None:
        if ax is None:
            raise Exception("Both or neither of fig, ax may be None")
    if ax is None:
        if fig is not None:
            raise Exception("Both or neither of fig, ax may be None")
        fig, ax = plt.subplots(figsize=(4, 4))
    if category_as_continuum:
        if use_gex_as_ca:
            igene = np.where(loom.ra[gene_ra] == category_ca)[0][0]
            c = loom[gex_layer][igene, :]
        else:
            c = loom.ca[category_ca]
        g = ax.scatter(loom.ca[x_ca],
                       loom.ca[y_ca],
                       s=s,
                       c=c,
                       cmap=color_palette)
        plt.colorbar(g, label=category_ca)
    else:
        shuffle = np.arange(loom.shape[1])
        np.random.shuffle(shuffle)
        if type(color_palette) == str:
            color_palette = sns.color_palette(color_palette)
        elif type(
                sns.color_palette('colorblind')) == sns.palettes._ColorPalette:
            pass
        else:
            raise Exception(
                "color palette must be str or seaborn color palette")
        category2color = {
            category: sns.color_palette(color_palette)[i]
            for i, category in enumerate(np.unique(loom.ca[category_ca]))
        }
        ax.scatter(
            loom.ca[x_ca][shuffle],
            loom.ca[y_ca][shuffle],
            s=s,
            c=[category2color[x] for x in loom.ca[category_ca][shuffle]])
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0],
                   marker='o',
                   color='w',
                   label=category,
                   markerfacecolor=category2color[category],
                   markersize=10)
            for category in np.unique(loom.ca[category_ca])
        ]
        if on_figure_annotation:
            import matplotlib.patheffects as pe
            for category in np.unique(loom.ca[category_ca]):
                mask = loom.ca[category_ca] == category
                x, y = loom.ca[x_ca][mask].mean(), loom.ca[y_ca][mask].mean()
                ax.annotate(
                    str(category), (x, y),
                    path_effects=[pe.withStroke(foreground='w', linewidth=2)],
                    ha='center',
                    va='center')

        # Create the figure
        if legend:
            ax.legend(handles=legend_elements,
                      bbox_to_anchor=(1, 1),
                      title=category_ca)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    left, right = ax.get_xlim()
    bottom, top = ax.get_ylim()
    scalefactor = 0.2
    newright = (right - left) * scalefactor + left
    newtop = (top - bottom) * scalefactor + bottom
    ax.spines['bottom'].set_bounds(
        left,
        newright,
    )
    ax.spines['left'].set_bounds(
        bottom,
        newtop,
    )
    ax.plot(newright, bottom, ">k", scalex=False, scaley=False, clip_on=False)
    ax.plot(left, newtop, "^k", scalex=False, scaley=False, clip_on=False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(x_ca, loc='left', fontsize=14)
    ax.set_ylabel(y_ca, loc='bottom', fontsize=14)

    return fig, ax


def gsea_plot(ranking,
              pathway2genelist_dict,
              left_label='Enriched for genes at beginning of ranking',
              right_label='Enriched for genes at end of ranking',
              figsize=(9 / 2, 7 / 2)):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from panopticon.analysis import get_enrichment_score
    import matplotlib

    palette = sns.color_palette('colorblind', 12)
    fig, axes = plt.subplots(len(pathway2genelist_dict.keys()) + 1,
                             1,
                             figsize=figsize,
                             height_ratios=[20] +
                             [1] * len(pathway2genelist_dict.keys()),
                             sharex=True)
    mins = []
    maxs = []
    for ikey, key in enumerate([y for y in pathway2genelist_dict.keys()]):
        es = get_enrichment_score(ranking,
                                  pathway2genelist_dict[key],
                                  presorted=True,
                                  return_es_curve=True,
                                  return_pvalue=True,
                                  use_fgsea=True)
        axes[0].plot(es.enrichment_score_curve,
                     label=key,
                     lw=3,
                     color=palette[ikey])
        maxs.append(np.max(es.enrichment_score_curve))
        mins.append(np.min(es.enrichment_score_curve))

        for x in np.where(np.isin(ranking, pathway2genelist_dict[key]))[0]:
            axes[ikey + 1].axvline(x)
        axes[ikey + 1].set_ylabel(key + '\n' * 5,
                                  rotation=0,
                                  ha='left',
                                  va='bottom')
        axes[ikey + 1].set_yticks([])
        axes[ikey + 1].yaxis.set_label_position("right")
        axes[ikey + 1].yaxis.tick_right()
        for side in ['top', 'bottom', 'right', 'left']:
            axes[ikey + 1].spines[side].set_visible(False)
        axes[ikey + 1].set_xticks([])
        if es.p_value == 0:
            key_with_pval = key + ', p<{0:.5g}'.format(1 / 10000)
        else:
            key_with_pval = key + ', p={0:.5g}'.format(es.p_value)
        axes[ikey + 1].set_ylabel(key_with_pval,
                                  rotation=0,
                                  ha='left',
                                  va='center')

    axes[0].legend(bbox_to_anchor=(1, 1))
    axes[0].set_ylabel('running enrichment score')
    axes[-1].set_xlabel('rank in gene list\n' + r'$\leftarrow$ ' + left_label +
                        ' ' * 20 + right_label + r' $\rightarrow$')
    axes[0].set_ylim([np.min(mins), np.max(maxs)])
    axes[0].spines['top'].set_position(('data', 0))
    axes[0].spines['bottom'].set_position(('data', 0))
    axes[0].spines['right'].set_visible(False)
    plt.tight_layout()
    return fig, axes


def boxenplot_with_statistics(
        data,
        x,
        y,
        hue=None,
        ax=None,
        annotate_hue_pvalues=False,
        annotate_hue_effect_size=False,
        annotate_hue_n=False,
        annotate_hue_pvalue_fmt_str='p: {0:.2f}',
        annotate_hue_effect_size_fmt_str='es: {0:.2f}',
        annotate_hue_n_fmt_str='n: {}, {}',
        annotate_pvalue_vs_all=False,
        annotate_effect_size_vs_all=False,
        annotate_n=False,
        annotate_pvalue_vs_all_fmt_str='p: {0:.2f}',
        annotate_effect_size_vs_all_fmt_str='es: {0:.2f}',
        annotate_n_fmt_str='n: {}',
        annotate_hue_single_line=False,
        effect_size='cohensd',
        pvalue='mannwhitney',
        custom_annotation_dict={},
        custom_annotation_fontsize=6,
        censor_inf_in_statistics=False):
    """

    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    if hue is None:
        hue_order = None
    else:
        hue_order = data[hue].unique()
    sns.boxenplot(data=data,
                  x=x,
                  y=y,
                  hue=hue,
                  width_method='exponential',
                  ax=ax,
                  outlier_prop=0.05,
                  saturation=1)

    from panopticon.visualization import legend_without_duplicate_labels
    legend_without_duplicate_labels(ax)
    if annotate_hue_pvalues or annotate_hue_effect_size or annotate_hue_n and len(
            custom_annotation_dict.keys()) == 0:
        if len(data[hue].unique()) != 2:
            raise Exception(
                "hue must be a categorical variable with 2 unique values")
        from scipy.stats import mannwhitneyu, ttest_ind, ttest_rel, wilcoxon
        if np.issubdtype(data[y].dtype, np.number):
            category_col = x
            continuous_col = y
            vertical_boxes = True
            ticklabels = ax.get_xmajorticklabels()
        else:
            category_col = y
            continuous_col = x
            vertical_boxes = False
            ticklabels = ax.get_ymajorticklabels()
        for ticklabel in ticklabels:
            category = ticklabel.get_text()
            hue1, hue2 = data[hue].unique()

            if pvalue in ['mannwhitney', 'ttest']:
                a = data[(data[hue] == hue1)
                         & (data[category_col].astype(str) == category
                            )][continuous_col].values
                b = data[(data[hue] == hue2)
                         & (data[category_col].astype(str) == category
                            )][continuous_col].values
                a = np.array([x for x in a if not np.isnan(x)])
                b = np.array([x for x in b if not np.isnan(x)])
                if censor_inf_in_statistics:
                    a = np.array([x for x in a if not np.isinf(x)])
                    b = np.array([x for x in b if not np.isinf(x)])

                if len(a) == 0 or len(b) == 0:
                    pval = np.nan
                else:
                    mw = mannwhitneyu(a, b, alternative='two-sided')
                    tt = ttest_ind(a, b)
                    if pvalue == 'mannwhitney':
                        pval = mw.pvalue
                    elif pvalue == 'ttest':
                        pval = tt[1]
            else:
                raise Exception(
                    "`pvalue` must be either \'mannwhitney\' or \'ttest\' ")
            if len(a) == 0 or len(b) == 0:
                es = np.nan
            else:
                if effect_size == 'cohensd':
                    from panopticon.utilities import cohensd
                    es = cohensd(a, b)
                elif effect_size == 'cles':
                    es = mw.statistic / len(a) / len(b)
                else:
                    raise Exception(
                        "effect_size must be either \'cohensd\' or \'cles\'")
            annotation_string = ''
            if len(a) == 0 or len(b) == 0:
                pass
            elif vertical_boxes:
                if annotate_hue_pvalues:
                    annotation_string += annotate_hue_pvalue_fmt_str.format(
                        pval) + '\n'
                if annotate_hue_effect_size:
                    annotation_string += annotate_hue_effect_size_fmt_str.format(
                        es) + '\n'
                if annotate_hue_n:
                    annotation_string += annotate_hue_n_fmt_str.format(
                        len(a), len(b)) + '\n'
                if annotate_hue_single_line:
                    anno_x, anno_y = ticklabel.get_position(
                    )[0], data[continuous_col].max()
                else:
                    anno_x, anno_y = ticklabel.get_position()[0], np.max(
                        np.hstack((a, b)))

                ax.annotate(annotation_string, (anno_x, anno_y),
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
                if annotate_hue_single_line:
                    anno_x, annoy_y = data[continuous_col].max(
                    ), ticklabel.get_position()[1]
                else:
                    anno_x, anno_y = np.max(np.hstack(
                        (a, b))), ticklabel.get_position()[1]

                ax.annotate(annotation_string, (anno_x, anno_y),
                            ha='left',
                            va='center')
    if annotate_pvalue_vs_all or annotate_effect_size_vs_all or annotate_n and len(
            custom_annotation_dict.keys()) == 0:
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
            if pvalue in ['mannwhitney', 'ttest']:
                a = data[data[category_col] == category][continuous_col].values
                b = data[data[category_col] != category][continuous_col].values
                a = np.array([x for x in a if not np.isnan(x)])
                b = np.array([x for x in b if not np.isnan(x)])
                if censor_inf_in_statistics:
                    a = np.array([x for x in a if not np.isinf(x)])
                    b = np.array([x for x in b if not np.isinf(x)])
                mw = mannwhitneyu(a, b, alternative='two-sided')
                tt = ttest_ind(a, b)
                if pvalue == 'mannwhitney':
                    pval = mw.pvalue
                elif pvalue == 'ttest':
                    pval = tt[1]
            else:
                raise Exception(
                    "`pvalue` must be either \'mannwhitney\' or \'ttest\' ")

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
                if annotate_pvalue_vs_all:
                    annotation_string += annotate_pvalue_vs_all_fmt_str.format(
                        pval) + '\n'
                if annotate_effect_size_vs_all:
                    annotation_string += annotate_effect_size_vs_all_fmt_str.format(
                        es) + '\n'
                if annotate_n:
                    annotation_string += annotate_n_fmt_str.format(
                        len(a)) + '\n'
                ax.annotate(annotation_string,
                            (ticklabel.get_position()[0], np.max(a)),
                            ha='center',
                            va='bottom')
            else:
                if annotate_pvalue_vs_all:
                    annotation_string += ' ' + annotate_pvalue_vs_all_fmt_str.format(
                        pval)
                if annotate_effect_size_vs_all:
                    annotation_string += '\n ' + annotate_effect_size_vs_all_fmt_str.format(
                        es)
                if annotate_n:
                    annotation_string += '\n' + annotate_n_fmt_str.format(
                        len(a))

                ax.annotate(annotation_string, (
                    np.max(a),
                    ticklabel.get_position()[1],
                ),
                            ha='left',
                            va='center',
                            annotation_clip=False)
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
            annotation_pos = np.max(
                data[data[category_col] == category][continuous_col].values)

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
