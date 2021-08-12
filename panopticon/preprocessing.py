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


def generate_cell_and_gene_quality_metrics(loom, layername):
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
                                 denominator=10**5):
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


def generate_standardized_layer(loom, layername, variance_axis='cell', batch_size=512, out_of_core_cell_threshold=20000):
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
            loom[layername+'_gene_standardized'] = "float64"
            for (ix, selection, view) in tqdm(loom.scan(axis=0,
                                                        batch_size=batch_size),
                                              total=loom.shape[1] // batch_size):
                loom[layername+'_gene_standardized'][selection,:] = scaler.fit_transform(view[layername][:,:].T).T

    elif variance_axis == 'cell':
        if loom.shape[1] < out_of_core_cell_threshold:
            loom[layername +
                 '_cell_standardized'] = StandardScaler().fit_transform(
                     loom[layername][:, :])
        else:
            scaler = StandardScaler()
            loom[layername+'_cell_standardized'] = "float64"
            for (ix, selection, view) in tqdm(loom.scan(axis=1,
                                                        batch_size=batch_size),
                                              total=loom.shape[1] // batch_size):
                loom[layername+'_cell_standardized'][:,selection] = scaler.fit_transform(view[layername][:,:])
    
