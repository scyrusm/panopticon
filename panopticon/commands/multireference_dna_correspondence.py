import numpy as np
import re
import loompy
import matplotlib.pyplot as plt
from panopticon.dna import multireference_dna_correspondence
from tqdm import tqdm


def multireference_dna_correspondence_main(loomfile,
                                           queries,
                                           segmentations,
                                           output=None):
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
            if operation == '==':
                querymask = (querymask &
                             (loom.ca[field].astype(str) == threshold))
            elif operation == '!=':
                querymask = querymask & (loom.ca[field].astype(str) !=
                                         threshold)
            elif operation == '<':
                querymask = querymask & (loom.ca[field].astype(str) <
                                         threshold)
            elif operation == '>':
                querymask = querymask & (loom.ca[field].astype(str) >
                                         threshold)
            elif operation == '>=':
                querymask = querymask & (loom.ca[field].astype(str) >=
                                         threshold)
            elif operation == '<=':
                querymask = querymask & (loom.ca[field].astype(str) <=
                                         threshold)

        list_of_correlations = multireference_dna_correspondence(
            loom, querymask, *segmentations)
    
    fig, ax = plt.subplots(len(segmentations) - 1, len(segmentations) - 1, figsize = (len(segmentations)*4, len(segmentations)*4))
    for iseg, xseg in enumerate(segmentations):
        np.savetxt(xseg+"_corr.txt", list_of_correlations[iseg]) 
        for jseg, yseg in enumerate(segmentations):
            if len(segmentations)>2:
                if iseg < jseg:
                    ax[jseg - 1, iseg].scatter(
                        list_of_correlations[iseg],
                        list_of_correlations[jseg],
                        alpha=0.5)
                    ax[jseg - 1, iseg].set_xlabel(xseg)
                    ax[jseg - 1, iseg].set_ylabel(yseg)
                    if iseg < jseg - 1:
                        ax[iseg, jseg - 1].axis('off')
            else:
                if iseg < jseg:
                    ax.scatter(
                        list_of_correlations[iseg],
                        list_of_correlations[jseg],
                        alpha=0.5)
                    ax.set_xlabel(xseg)
                    ax.set_ylabel(yseg)
    if output:
        plt.savefig(output)
    else:
        plt.show()
