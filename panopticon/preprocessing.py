import numpy as np


def generate_gene_variances(loom, layername):
    """Computes the variance (across cells) of genes, and assigns this to the row attribute `GeneVar`

    Parameters
    ----------
    loom : LoomConnection
        The LoomConnection instance upon which gene variance will be calculated
        
    layername : str
        The name of the layer of the LoomConnection upon which the gene variance will be calculated.
        

    Returns
    -------
    Writes a new variable `Genevar` as a row attribute to LoomConnection `loom`.Returns NoneType object.

    
    """
    loom.ra['GeneVar'] = loom[layername].map([np.var])[0]


def generate_cell_and_gene_quality_metrics(loom,
                                           layername,
                                           gene_ra='gene',
                                           ribosomal_qc=False,
                                           mitochondrial_qc=False,
                                           ribosomal_gene_mask=None,
                                           mitochondrial_gene_mask=None,
                                           verbose=False):
    """Calculates multiple QC-related quantities and writes them to the LoomConnection instance specified in loom:
   
    - nCountsForCell: absolute number of non-normalized counts for a given cell (across all genes)
    - nCountsForGene: absolute number of non-normalized counts for a give gene (across all cells)
    - RibosomalRelativeMeanAbsolutionDeviation:  this is madrp / meanrp, where madrp is the mean absolute deviation of ribosomal protein genes (RP*) and meanrp is the mean expression of ribosomal protein genes (RP*) (means, mads compute on a cell-by-cell basis).
    - RibosomalMaxOverMean: this is maxrp / meanrp, where maxrp is the maximum RP* expression over all RP* genes; meanrp is as above (max, mean on a cell-by-cell basis).
    - MitochondrialMean: this is the mean expression of genes from the mitochondrial genome.
    - MitochondrialRelativeMeanAbsoluteDeviation: this is madmt / meanmt, where madmt is the mean absolute deviation of mitochondrial genes (MT*) and meanmt is the mean expression of mitochondrial protein genes (MT*) (means, mads computed on a cell-by-cell basis). 
    - AllGeneRelativeMeanAbsoluteDeviation: this is madall / meanall, where madall is the MAD over all genes, and meanall is the mean over all genes (mad, mean on a cell-by-cell basis).
    - nGene: number of unique genes expressed (>0) on a cell-by-cell basis (also known as complexity).
    - nCell: number of unique cells expressing a gene (>0), on a gene-by-gene basis.

    Parameters
    ----------
    loom : LoomConnection
        The LoomConnection instance upon which cell and gene quality metrics will be annotated.
        
    layername : str
        The name of the layer upon which cell and gene quality metrics will be calculated.
        
    gene_ra : str
         The row attribute used for genes (recommended to use the HUGO names for genes, as this is used by default to determine which genes are mitochondrial or ribosomal) (Default value = 'gene')
    ribosomal_qc : bool
         If True, will compute ribosomal-based QC metrics (RibosomalRelativeMeanAbsolutionDeviation:, RibosomalMaxOverMean) (Default value = False)
    mitochondrial_qc : bool
         If True, will compute mitochondrial-based QC metric ((Default value = False)
    ribosomal_gene_mask : boolean mask with length equal to loom.shape[0]
         Specifies the indices (via Boolean mask) of the rows corresponding to ribosomal genes.  If None, will generate mask from all genes whose name starts with 'RP' or 'Rp' (Default value = None)
    mitochondrial_gene_mask : boolean mask with length equal to loom.shape[0]
         Specifies the indices (via Boolean mask) of the rows corresponding to mitochrondrial genes.  If None, will generate mask from all genes whose name starts with 'mt.' or 'mt.' (Default value = None)
    verbose : bool
         If True, will print notification whenever a QC calculation is completed.  (Default value = False)

    Returns
    -------
    Calculated qualities are written to LoomConnection `loom`. Returns NoneType object.

    
    """
    from statsmodels import robust
    ncounts_for_cell = loom[layername].map([np.sum], axis=1)[0]
    if verbose:
        print("Completed ncounts_for_cell")
    ncounts_for_gene = loom[layername].map([np.sum], axis=0)[0]
    if verbose:
        print("Completed ncounts_for_gene")

    if (ncounts_for_cell == (
            ncounts_for_cell).astype(int)).sum() == len(ncounts_for_cell):
        loom.ca['nCountsForCell'] = ncounts_for_cell
        loom.ra['nCountsForGene'] = ncounts_for_gene
        israw = True
    else:
        print(
            "layer does not appear to be raw counts; skipping raw count sums per cell, gene, as well as statistics involving percentages of overall counts"
        )
        israw = False

    if ribosomal_qc:
        if ribosomal_gene_mask is None:
            ribosomal_gene_mask = np.array([
                x.startswith('RP') or x.startswith('Rp')
                for x in loom.ra[gene_ra]
            ])
        if len(ribosomal_gene_mask)!=loom.shape[0]:
            raise Exception("ribosomal_gene_mask must be boolean mask with length equal to the number of rows of loom")
        madrp, meanrp, maxrp = loom[layername].map(
            [robust.mad, np.mean, np.max],
            axis=1,
            selection=ribosomal_gene_mask.nonzero()[0])
        loom.ca['RibosomalRelativeMeanAbsoluteDeviation'] = madrp / meanrp
        loom.ca['RibosomalMaxOverMean'] = maxrp / meanrp
        if israw:
            loom.ca['RibosomalCountFraction'] = loom[layername].map(
                [np.sum], axis=1, selection=ribosomal_gene_mask.nonzero()
                [0])[0] / ncounts_for_cell
        if verbose:
            print("Completed ribosomal QC")
    if mitochondrial_qc:
        if mitochondrial_gene_mask is None:
            mitochondrial_gene_mask = np.array([
                x.lower().startswith('mt-') or x.lower().startswith('mt.')
                for x in loom.ra[gene_ra]
            ])
            print(loom.ra[gene_ra][mitochondrial_gene_mask])
        if len(mitochondrial_gene_mask)!=loom.shape[0]:
            raise Exception("mitochondrial_gene_mask must be boolean mask with length equal to the number of rows of loom")
        madmt, meanmt, maxmt = loom[layername].map(
            [robust.mad, np.mean, np.max],
            axis=1,
            selection=mitochondrial_gene_mask.nonzero()[0])
        loom.ca['MitochondrialMean'] = meanmt
        loom.ca['MitochondrialRelativeMeanAbsoluteDeviation'] = madmt / meanmt
        if israw:
            loom.ca['MitochondrialCountFraction'] = loom[layername].map(
                [np.sum],
                axis=1,
                selection=mitochondrial_gene_mask.nonzero()
                [0])[0] / ncounts_for_cell
        if verbose:
            print("Completed mitochondrial QC")

    madall, meanall = loom[layername].map([robust.mad, np.mean], axis=1)
    loom.ca['AllGeneRelativeMeanAbsoluteDeviation'] = madall / meanall
    if verbose:
        print("Completed madall/meanall")

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
    if verbose:
        print("Completed nGene, nCell")


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


def generate_count_normalization(loom,
                                 raw_count_layername,
                                 output_layername,
                                 denominator=10**4):
    """Generates a new layer with log2(transcripts per denominator).  By default this will produce a log2(TP100k+1) layer; adjusting the denominator will permit e.g. a log2(TP10k+1) or log2(TPM+1), etc.

    Parameters
    ----------
    loom : The LoomConnection instance upon which a count normalized layer will be computed.
        
    raw_count_layername : The layername corresponding to raw counts, or normalized log counts (counts of TPM).
        
    output_layername : The desired output layername.
        
    denominator : The denominator for transcript count normalization (e.g. transcripts per denominator).
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


def generate_standardized_layer(loom,
                                layername,
                                variance_axis='cell',
                                batch_size=512,
                                out_of_core_cell_threshold=20000):
    """

    Parameters
    ----------
    loom : The LoomConnection instance upon which the standardized layer will be added.
        
    layername : The name of the layer which will be standardized.
        
    variance_axis : The axis over which the standardization will proceed (i.e. over cells or over genes).
        (Default value = 'cell')
    batch_size :
         (Default value = 512)
    out_of_core_cell_threshold :
         (Default value = 20000)

    Returns
    -------

    
    """
    from tqdm import tqdm

    if layername.endswith('_standardized'):
        raise Exception(
            "It appears that layer is already standardized; note that _standardized suffix is reserved"
        )
    if variance_axis not in ['gene', 'cell']:
        raise Exception(
            "Can only set variance to be 1 along cell axis or gene axis")
    from sklearn.preprocessing import StandardScaler
    if variance_axis == 'gene':
        if loom.shape[1] < out_of_core_cell_threshold:
            loom[layername +
                 '_gene_standardized'] = StandardScaler().fit_transform(
                     loom[layername][:, :].T).T
        else:
            scaler = StandardScaler()
            loom[layername + '_gene_standardized'] = "float64"
            for (ix, selection,
                 view) in tqdm(loom.scan(axis=0, batch_size=batch_size),
                               total=loom.shape[1] // batch_size):
                loom[layername + '_gene_standardized'][
                    selection, :] = scaler.fit_transform(
                        view[layername][:, :].T).T

    elif variance_axis == 'cell':
        if loom.shape[1] < out_of_core_cell_threshold:
            loom[layername +
                 '_cell_standardized'] = StandardScaler().fit_transform(
                     loom[layername][:, :])
        else:
            scaler = StandardScaler()
            loom[layername + '_cell_standardized'] = "float64"
            for (ix, selection,
                 view) in tqdm(loom.scan(axis=1, batch_size=batch_size),
                               total=loom.shape[1] // batch_size):
                loom[layername +
                     '_cell_standardized'][:,
                                           selection] = scaler.fit_transform(
                                               view[layername][:, :])


def generate_antibody_prediction(loom,
                                 raw_antibody_counts_df=None,
                                 antibodies=None,
                                 pseudocount=1,
                                 overwrite=False,
                                 only_generate_zscore=False):
    """
    This approach takes some inspiration from the dsb approach: https://doi.org/10.1101/2020.02.24.963603.
    However, there is no use of isotypes.  Therefore, is amounts only to "step 1" of that procedure. This routine can also take into
    acccount the distribution of antibody/feature barcode counts from empty droplet using the optional argument `raw_antibody_counts_df`. 
    This routine uses a 2-component Gaussian-mixture model on z-scored log1p antibody counts (with or without background correction) 
    to predict whether a given cells is "positive" or "negative" for the feature barcode in question.

    Parameters
    ----------
    loom : LoomConnection
        The LoomConnection object on which to make antibody predictions. 
        
    raw_antibody_counts_df : pandas.DataFrame object, or NoneType
        If `None`(Default value = None)
    antibodies : str, or iterable of strings
        interable of strings, each of which should be a column attribute of `loom`. Each should represent the column attributes indicating the raw counts of a particular feature barcode. (Default value = None)
    pseudocount : int
        Indicates the pseudocount of the feature barcode to use when taking the log(count+pseudocount) (Default value = 1)
    overwrite : bool
        If True, will overwrite existing prediction. If False, with raise an Exception if prediction had previously been generated. (Default value = False)
    only_generate_zscore : bool
        If True, will only generate log1pseudo(counts) and z-score thereof, but will not compute a Gaussian mixture model nor make subsequent prediction. (Default value = False)

    Returns
    -------
    None (all output are written as column attributes to LoomConnection `loom`)

    """
    if raw_antibody_counts_df is None and antibodies is None:
        raise Exception(
            "one of antibodies, raw_antibody_counts_df must be other than None"
        )
    if raw_antibody_counts_df is not None and antibodies is not None:
        raise Exception(
            "only one of antibodies, raw_antibody_counts_df must be other than None"
        )

    if raw_antibody_counts_df is not None:
        antibodies = raw_antibody_counts_df.columns

    for antibody in antibodies:
        if antibody not in loom.ca.keys():
            raise Exception(
                "raw_antibody_count_df must be prepared such that columns match column attributes in loom corresponding to raw antibody conjugate counts"
            )

        new_ca_name = antibody + '_zscore_log{}p'.format(pseudocount)
        if new_ca_name in loom.ca.keys() and overwrite == False:
            raise Exception(
                "{} already in loom.ca.keys(); rename antibody column attribute and re-run, or set overwrite argument to True"
                .format(new_ca_name))
        if raw_antibody_counts_df is not None:
            if pseudocount == 1:
                logcounts_cells = np.log1p(loom.ca[antibody])
                logcounts_empty_droplets = np.log1p(
                    raw_antibody_counts_df[antibody])
            else:
                logcounts_cells = np.log(loom.ca[antibody] + pseudocount)
                logcounts_empty_droplets = np.log(
                    raw_antibody_counts_df[antibody].values + pseudocount)
            logcounts_empty_droplets_mean = np.nanmean(
                logcounts_empty_droplets)
            logcounts_empty_droplets_std = np.nanstd(logcounts_empty_droplets)
            loom.ca[new_ca_name] = (
                logcounts_cells -
                logcounts_empty_droplets_mean) / logcounts_empty_droplets_std
        else:
            if pseudocount == 1:
                logcounts_cells = np.log1p(loom.ca[antibody])
            else:
                logcounts_cells = np.log(loom.ca[antibody] + pseudocount)
            logcounts_cells_mean = np.nanmean(logcounts_cells)
            logcounts_cells_std = np.nanstd(logcounts_cells)
            loom.ca[new_ca_name] = (logcounts_cells -
                                    logcounts_cells_mean) / logcounts_cells_std

        if not only_generate_zscore:
            from sklearn import mixture

            model = mixture.GaussianMixture(n_components=2)
            prediction_ca_name = antibody + '_prediction'
            if prediction_ca_name in loom.ca.keys() and overwrite == False:
                raise Exception(
                    "{} already in loom.ca.keys(); rename antibody column attribute and re-run, or set overwrite argument to True"
                    .format(prediction_ca_name))
            if np.isnan(loom.ca[new_ca_name]).sum() > 0:
                cellmask = ~np.isnan(loom.ca[new_ca_name])
                model.fit(loom.ca[new_ca_name][cellmask].reshape(-1, 1))
                predictions = []
                for val in loom.ca[new_ca_name]:
                    if np.isnan(val):
                        predictions.append(np.nan)
                    else:
                        predictions.append(
                            model.predict(np.array(val).reshape(-1, 1))[0])
                predictions = np.array(predictions)

            else:
                predictions = model.fit_predict(loom.ca[new_ca_name].reshape(
                    -1, 1))

            if model.means_[0][0] > model.means_[1][0]:
                predictions = 1 - predictions

            loom.ca[prediction_ca_name] = np.array(predictions)


def generate_guide_rna_prediction(
        loom,
        guide_rnas,
        nguide_ca='nGuide',
        nguide_reads_ca='nGuideReads',
        cell_prediction_summary_ca='CellGuidePrediction',
        overwrite=False,
        only_generate_log2=False,
        ncell_threshold_for_guide=10,
        nguide_threshold_for_cell=10):
    """
    This approach is inspired by Replogle et a. 2018 (https://doi.org/10.1038/s41587-020-0470-y). However, instead of a Gaussian/Poisson mixture, this routine uses a Poisson/Poisson mixture. 
    This routine uses the pomegranate package (https://github.com/jmschrei/pomegranate). 

    Parameters
    ----------
    loom : LoomConnection
        A LoomConnection object upon which guide rna predictions will be made
        
    guide_rnas : iterable of strings
        a list or other iterable of the strings, each corresponding to a column attribute of `loom` indicate the raw counts of a given guide RNA over cells 
    nguide_ca : str
         QC metric, indicating the name of the column attribute to use to indicate the number of predicted guide RNAs for a cell (Default value = 'nGuide')
    nguide_reads_ca :
         QC metric, indicating the name of the column attribute to use to indicate the total number of guide RNA reads for a cell(Default value = 'nGuideReads')
    cell_prediction_summary_ca : str
         Indicates the name of the column attribute to use to indicate a summary of positively-predicted guide RNAs for a cell(Default value = 'CellGuidePrediction')
    overwrite : bool
         If False, will raise exception if requested column attributes have already been written.  If True, will overwrite existing column attributes. (Default value = False)
    only_generate_log2 : bool
         If true, will generate log2 guide RNA counts, but will not apply any mixture model prediction. (Default value = False)
    ncell_threshold_for_guide : int
         Threshold for the number of cells wherein guide should have nonzero counts for mixture model to attempt prediction. (Default value = 10)
    nguide_threshold_for_cell : int
         Threshold for the number of guides to be detected in a given cell to attempt to make a prediction for that particular cell. (Default value = 10)

    Returns
    -------

    """
    from panopticon.utilities import import_check
    exit_code = import_check("pomegranate",
                             'conda install -c anaconda pomegranate')
    if exit_code != 0:
        return
    import pandas as pd

    if nguide_reads_ca in loom.ca.keys() and overwrite == False:
        raise Exception(
            "{} already in loom.ca.keys(); if intended, set overwrite argument to True"
            .format(nguide_reads_ca))

    guide_rna_dfs = []
    for guide_rna in guide_rnas:
        guide_rna_dfs.append(
            pd.DataFrame(loom.ca[guide_rna], columns=[guide_rna], copy=True))
    guide_rna_dfs = pd.concat(guide_rna_dfs, axis=1)
    loom.ca[nguide_reads_ca] = guide_rna_dfs.sum(axis=1).values
    threshold_for_cell_mask = loom.ca[
        nguide_reads_ca] >= nguide_threshold_for_cell
    prediction_ca_names = []
    for guide_rna in guide_rnas:
        if guide_rna not in loom.ca.keys():
            raise Exception(
                "raw_antibody_count_df must be prepared such that columns match column attributes in loom corresponding to raw antibody conjugate counts"
            )

        new_ca_name = guide_rna + '_log2'
        if new_ca_name in loom.ca.keys() and overwrite == False:
            raise Exception(
                "{} already in loom.ca.keys(); rename guide column attribute and re-run, or set overwrite argument to True"
                .format(new_ca_name))

        loom.ca[new_ca_name] = np.log2(loom.ca[guide_rna])
        if not only_generate_log2:
            from pomegranate import GeneralMixtureModel, PoissonDistribution

            prediction_ca_name = guide_rna + '_prediction'
            prediction_ca_names.append(prediction_ca_name)
            if prediction_ca_name in loom.ca.keys() and overwrite == False:
                raise Exception(
                    "{} already in loom.ca.keys(); rename guide rna column attribute and re-run, or set overwrite argument to True"
                    .format(prediction_ca_name))
            if (~np.isfinite(loom.ca[new_ca_name])).sum() > 0:
                cellmask = np.isfinite(loom.ca[new_ca_name])
                if cellmask.sum(
                ) >= ncell_threshold_for_guide:  # have minimum cells for guide
                    model = GeneralMixtureModel.from_samples(
                        [PoissonDistribution, PoissonDistribution],
                        n_components=2,
                        X=loom.ca[new_ca_name][cellmask.nonzero()[0]].reshape(
                            -1, 1))
                    predictions = []
                    for val in loom.ca[new_ca_name]:
                        if not np.isfinite(val):
                            predictions.append(np.nan)
                        else:
                            predictions.append(
                                model.predict(np.array(val).reshape(-1, 1))[0])
                else:
                    predictions = [0] * loom.shape[1]
                predictions = np.array(predictions)

            else:
                #print(guide_rna, loom.ca[guide_rna].sum())
                # print('Warning:  pomegrante Poisson/Normal mixture model has predicted a Poisson component with greater log(UMI+1) counts than normal component.  This is unusual behavior!')
                model = GeneralMixtureModel.from_samples(
                    [PoissonDistribution, PoissonDistribution],
                    n_components=2,
                    X=loom.ca[new_ca_name].reshape(-1, 1))
                #model.fit(loom.ca[new_ca_name].reshape(-1, 1))
                predictions = model.predict(loom.ca[new_ca_name].reshape(
                    -1, 1))

            if loom.ca[new_ca_name][np.array(predictions) == 0].mean(
            ) > loom.ca[new_ca_name][np.array(predictions) == 1].mean():
                predictions = 1 - predictions
            predictions = np.array(predictions)
            predictions = np.nan_to_num(predictions, nan=0.0)
            predictions *= threshold_for_cell_mask
            loom.ca[prediction_ca_name] = predictions

    guide_prediction_dfs = []
    for prediction_ca_name in prediction_ca_names:
        guide_prediction_dfs.append(
            pd.DataFrame(loom.ca[prediction_ca_name],
                         columns=[prediction_ca_name],
                         copy=True))
    guide_prediction_dfs = pd.concat(guide_prediction_dfs, axis=1)
    loom.ca[nguide_ca] = guide_prediction_dfs.sum(axis=1).values

    loom.ca[cell_prediction_summary_ca] = guide_prediction_dfs.apply(
        lambda x: '+'.join(guide_prediction_dfs.columns[np.where(x == 1)[0]]),
        axis=1).values


def get_clustering_based_outlier_prediction(
        loom,
        max_cluster_fraction_break_threshold=0.99,
        cluster_proportion_minimum=0.01):
    """

    Parameters
    ----------
    loom :
        
    max_cluster_fraction_break_threshold :
         (Default value = 0.99)
    cluster_proportion_minimum :
         (Default value = 0.01)

    Returns
    -------

    """
    import pandas as pd

    levels = [x for x in loom.ca.keys() if x.startswith('ClusteringIteration')]
    max_clustering_level = np.max(
        [int(x.replace('ClusteringIteration', '')) for x in levels])
    stub_cells = []
    mask = np.array([True] * loom.shape[1])
    for level in range(max_clustering_level + 1):
        #if level == 0:
        counts = pd.DataFrame(loom.ca['ClusteringIteration{}'.format(level)]
                              [mask])[0].value_counts()
        fractions = counts / counts.sum()

        if fractions.values[0] < max_cluster_fraction_break_threshold:
            break
        outlier_clusters = fractions[
            fractions < cluster_proportion_minimum].index.values
        mask *= ~np.isin(loom.ca['ClusteringIteration{}'.format(level)],
                         outlier_clusters)
    return mask
