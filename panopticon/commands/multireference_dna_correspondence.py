from tqdm import tqdm
import re

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from panopticon.dna import multireference_dna_correspondence
from panopticon.legacy import get_module_score_loom

# import loompy after the main packages, because sometimes it breaks packages that are imported further:
import loompy

def multireference_dna_correspondence_main(loomfile,
                                           queries,
                                           segmentations,
                                           stratifyby=None,
                                           output=None,
                                           modulescore_signature=None):
    """

    Parameters
    ----------
    loomfile :
        
    queries :
        
    segmentations :
        
    stratifyby :
         (Default value = None)
    output :
         (Default value = None)
    modulescore_signature :
         (Default value = None)

    Returns
    -------

    """
    if len(segmentations) < 2:
        raise Exception(
            "Must name two difference CNV segmentations to perform multireference comparison"
        )
    with loompy.connect(loomfile, validate=False) as loom:
        querymask = np.array([True] * loom.shape[1])
        for query in queries:
            field, threshold = re.split("<|>|<=|>=|==|!=", query)
            operation = query.replace(field, '').replace(threshold, '')
            field = field.strip(' ')
            threshold = threshold.strip(' ')
            field_type = loom.ca[field].dtype.type if loom.ca[field].dtype != 'O' else str
            if operation == '==':
                querymask = (querymask &
                             (loom.ca[field].astype(field_type) == field_type(threshold)))
            elif operation == '!=':
                querymask = querymask & (loom.ca[field].astype(field_type) !=
                                         field_type(threshold))
            elif operation == '<':
                querymask = querymask & (loom.ca[field].astype(field_type) <
                                         field_type(threshold))
            elif operation == '>':
                querymask = querymask & (loom.ca[field].astype(field_type) >
                                         field_type(threshold))
            elif operation == '>=':
                querymask = querymask & (loom.ca[field].astype(field_type) >=
                                         field_type(threshold))
            elif operation == '<=':
                querymask = querymask & (loom.ca[field].astype(field_type) <=
                                         field_type(threshold))
            

        list_of_correlations = multireference_dna_correspondence(
            loom, querymask, *segmentations)
        stratum_masks = []
        stratum_names = []
        if stratifyby is not None:
            for stratum in np.unique(loom.ca[querymask][stratifyby]):
                stratum_names.append(stratum)
                stratum_masks.append(loom.ca[querymask][stratifyby] == stratum)

        else:
            stratum_names.append("All")
            stratum_masks.append(np.array([True] * loom.shape[1])[querymask])

        if modulescore_signature is not None:
            module_score = get_module_score_loom(
                loom, modulescore_signature, querymask=querymask)


#    fig, ax = plt.subplots(len(segmentations) - 1, len(segmentations) - 1, figsize = (len(segmentations)*4, len(segmentations)*4))
    if modulescore_signature is None:
        module_score_tick = 0
    else:
        module_score_tick = 1
    if stratifyby is None:
        fig, ax = plt.subplots(
            2 + module_score_tick,
            (len(segmentations) - 1) * len(segmentations) // 2,
            figsize=(len(segmentations) * (len(segmentations) * 2 - 1), 8))
    else:
        fig, ax = plt.subplots(
            (2 + module_score_tick) * len(stratum_masks),
            (len(segmentations) - 1) * len(segmentations) // 2,
            figsize=(len(segmentations) * (len(segmentations) * 2 - 1),
                     8 * len(stratum_masks)))
    ## Getting the figsize right is annoying, perhaps there's a better solution

    counter = 0
    for iseg, xseg in enumerate(segmentations):
        np.savetxt(xseg + "_corr.txt", list_of_correlations[iseg])
        for jseg, yseg in enumerate(segmentations):
            for imask, stratum_mask in enumerate(stratum_masks):
                if len(segmentations) <= 2:
                    relax1 = ax[0 + imask * (2 + module_score_tick)]
                    relax2 = ax[1 + imask * (2 + module_score_tick)]
                else:
                    relax1 = ax[0 + imask * (2 + module_score_tick), counter]
                    relax2 = ax[1 + imask * (2 + module_score_tick), counter]
                if iseg < jseg:
                    if stratifyby is not None:
                        stratify_title_addendum = ", {} == {}".format(
                            stratifyby, stratum_names[imask])
                    else:
                        stratify_title_addendum = ""
                    relax1.scatter(
                        np.array(list_of_correlations[iseg])[stratum_mask],
                        np.array(list_of_correlations[jseg])[stratum_mask],
                        alpha=0.5)
                    relax1.set_xlabel(xseg)
                    relax1.set_ylabel(yseg)
                    relax1.set_title(
                        "KT Correlations for Different Segmentations" +
                        stratify_title_addendum)
                    sns.distplot(
                        np.array(list_of_correlations[iseg])[stratum_mask] -
                        np.array(list_of_correlations[jseg])[stratum_mask],
                        ax=relax2)

                    relax2.set_title(
                        "Distribution of Differences in KT Correlation" +
                        stratify_title_addendum)
                    relax2.set_xlabel("{} - {}".format(xseg, yseg))
                    relax2.set_ylabel("Probability")

                    if modulescore_signature is not None:
                        if len(segmentations) <= 2:
                            relax3 = ax[2 + imask * (2 + module_score_tick)]
                        else:
                            relax3 = ax[2 + imask * (2 + module_score_tick), counter]
                        relax3.scatter(
                        np.array(list_of_correlations[iseg])[stratum_mask] -
                        np.array(list_of_correlations[jseg])[stratum_mask],
                                module_score[stratum_mask])


            counter += 1
    if output:
        plt.savefig(output)
    else:
        plt.show()
