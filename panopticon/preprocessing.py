import numpy as np

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


def generate_cell_and_gene_quality_metrics(loom, layername):
    """

    Parameters
    ----------
    loom :
        
    layername :
        

    Returns
    -------

    
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
    """Generates a new layer with log2 (TP_something)                                                                                                                                                                                          
                                                                                                                                                                                                                                               
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
