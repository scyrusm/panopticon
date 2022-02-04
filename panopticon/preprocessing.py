import numpy as np


def generate_gene_variances(loom, layername):
    """
    Computes the variance (across cells) of genes, and assigns this to the row attribute "GeneVar" 

    Parameters
    ----------
    loom : The LoomConnection instance upon which gene variance will be calculated
    
    layername : The name of the layer of the LoomConnection upon which the gene variance will be calculated.
        

    Returns
    -------
    None
    
    """
    loom.ra['GeneVar'] = loom[layername].map([np.var])[0]


def generate_cell_and_gene_quality_metrics(loom,
                                           layername,
                                           gene_ra='gene',
                                           ribosomal_qc=False,
                                           mitochondrial_qc=False,
                                           ribosomal_gene_mask=None,
                                           mitochondrial_gene_mask=None):
    """
    Calculates five quantities and writes them to the LoomConnection instance specified in loom:

    - RibosomalRelativeMeanAbsolutionDeviation:  this is madrp / meanrp, where madrp is the mean absolute deviation of ribosomal protein genes (RP*) and meanrp is the mean expression of ribosomal protein genes (RP*) (means, mads compute on a cell-by-cell basis).
    - RibosomalMaxOverMean: this is maxrp / meanrp, where maxrp is the maximum RP* expression over all RP* genes; meanrp is as above (max, mean on a cell-by-cell basis).
    - AllGeneRelativeMeanAbsoluteDeviation: this is madall / meanall, where madall is the MAD over all genes, and meanall is the mean over all genes (mad, mean on a cell-by-cell basis).  
    - nGene: number of unique genes expressed (>0) on a cell-by-cell basis.
    - nCell: number of unique cells expressing a gene (>0), on a gene-by-gene basis.  

    Parameters
    ----------
    loom : The LoomConnection instance upon which cell and gene quality metrics will be annotated.
        
    layername : The name of the layer upon which cell and gene quality metrics will be calculated.
        

    Returns
    -------
    None

    
    """
    from statsmodels import robust
    ncounts_for_cell = loom[layername].map([np.sum], axis=1)[0]
    ncounts_for_gene = loom[layername].map([np.sum], axis=0)[0]

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
        madrp, meanrp, maxrp = loom[layername].map(
            [robust.mad, np.mean, np.max],
            axis=1,
            selection=ribosomal_gene_mask.nonzero()[0])
        loom.ca['RibosomalRelativeMeanAbsoluteDeviation'] = madrp / meanrp
        loom.ca['RibosomalMaxOverMean'] = maxrp / meanrp
        if israw:
            loom.ca['RibosomalCountFraction'] = loom[layername].map(
                [np.sum], axis=1,
                selection=ribosomal_gene_mask.nonzero()[0])[0] / ncounts_for_cell
    if mitochondrial_qc:
        if mitochondrial_gene_mask is None:
            mitochondrial_gene_mask = np.array([
                x.lower().startswith('mt-') or x.lower().startswith('mt.')
                for x in loom.ra[gene_ra]
            ])
        madmt, meanmt, maxmt = loom[layername].map(
            [robust.mad, np.mean, np.max],
            axis=1,
            selection=mitochondrial_gene_mask.nonzero()[0])
        loom.ca['MitochondrialMean'] = meanmt
        loom.ca['MitochondrialRelativeMeanAbsoluteDeviation'] = madmt / meanmt
        if israw:
            loom.ca['MitochondrialCountFraction'] = loom[layername].map(
                [np.sum], axis=1,
                selection=mitochondrial_gene_mask.nonzero()[0])[0] / ncounts_for_cell

    madall, meanall = loom[layername].map([robust.mad, np.mean], axis=1)
    loom.ca['AllGeneRelativeMeanAbsoluteDeviation'] = madall / meanall

    def complexity(vec):
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
    None

    
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
    output_layername :  The desired output layername.                                                                                                                                                                                                                         
    denominator : The denominator for transcript count normalization (e.g. transcripts per denominator).                                                                                                                                                                                                                           
        (Default value = 10**5)                                                                                                                                                                                                                
                                                                                                                                                                                                                                               
    Returns                                                                                                                                                                                                                                    
    -------                                                                                                                                                                                                                                    
    None                                                                                                                                                                                                                                               
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
                                                                                                                                                                                                                                               
    Returns                                                                                                                                                                                                                                    
    -------                                                                                                                                                                                                                                    
    None 
    
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


def generate_guide_rna_prediction(loom,
                                  guide_rnas,
                                  overwrite=False,
                                  only_generate_log2=False):
    from panopticon.utilities import import_check
    import_check("pomegranate", 'conda install -c anaconda pomegranate')

    for guide_rna in guide_rnas:
        if guide_rna not in loom.ca.keys():
            raise Exception(
                "raw_antibody_count_df must be prepared such that columns match column attributes in loom corresponding to raw antibody conjugate counts"
            )

        new_ca_name = guide_rna + '_log2'
        if new_ca_name in loom.ca.keys() and overwrite == False:
            raise Exception(
                "{} already in loom.ca.keys(); rename antibody column attribute and re-run, or set overwrite argument to True"
                .format(new_ca_name))

        loom.ca[new_ca_name] = np.log2(loom.ca[guide_rna])
        if not only_generate_log2:
            from pomegranate import GeneralMixtureModel, NormalDistribution, PoissonDistribution, LogNormalDistribution

#            model = GeneralMixtureModel.from_samples(
#                [PoissonDistribution, PoissonDistribution],
#                n_components=2,
#                X=counts.reshape(-1, 1))

            prediction_ca_name = guide_rna + '_prediction'
            if prediction_ca_name in loom.ca.keys() and overwrite == False:
                raise Exception(
                    "{} already in loom.ca.keys(); rename antibody column attribute and re-run, or set overwrite argument to True"
                    .format(prediction_ca_name))
            if (~np.isfinite(loom.ca[new_ca_name])).sum() > 0:
                cellmask = np.isfinite(loom.ca[new_ca_name])
                print(cellmask.sum(), guide_rna)
                if cellmask.sum()>10:
                    model = GeneralMixtureModel.from_samples(
                        [PoissonDistribution, PoissonDistribution],
                        n_components=2,
                        X=loom.ca[new_ca_name][cellmask.nonzero()[0]].reshape(-1, 1))
                    predictions = []
                    for val in loom.ca[new_ca_name]:
                        if not np.isfinite(val):
                            predictions.append(np.nan)
    
                        else:
                            predictions.append(
                                model.predict(np.array(val).reshape(-1, 1))[0])
                else:
                    predictions = [0]*loom.shape[1]
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
            loom.ca[prediction_ca_name] = predictions


def get_clustering_based_outlier_prediction(
        loom,
        max_cluster_fraction_break_threshold=0.99,
        cluster_proportion_minimum=0.01):
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
