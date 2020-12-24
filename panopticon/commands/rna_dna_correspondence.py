#! /usr/bin/env python

from tqdm import tqdm
import argparse
from warnings import warn

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats, sparse
from scipy.stats.mstats import theilslopes
import pymannkendall as mk

from panopticon import wme, dna

# import loompy after the main packages, because sometimes it breaks packages that are imported further:
import loompy

def rna_dna_correspondence_main(loomfile, args_seg_file, args_patient_column,
                                args_patient, args_time_point, output=None):
    """

    Parameters
    ----------
    loomfile :
        
    args_seg_file :
        
    args_patient_column :
        
    args_patient :
        
    args_time_point :
        
    output :
         (Default value = None)

    Returns
    -------

    """
    with loompy.connect(loomfile, validate=False) as loom:
        genes = loom.ra['gene']
        metadata = pd.DataFrame(loom.ca['patient_ID']).astype(str)
        metadata.columns = ['patient_ID']
        metadata['complexity'] = loom.ca['complexity']
        metadata['cell_type'] = loom.ca['cell_type']
        metadata['time_point'] = loom.ca['time_point']

        list_of_gene_windows = wme.get_list_of_gene_windows(
            genes, window_size=800, window_step=300)

        segmentation = pd.read_table(args_seg_file)

        copy_ratio_dict = dna.segmentation_to_copy_ratio_dict(
            genes,
            segmentation,
            chrom_col='Chromosome',
            start_col='Start.bp',
            end_col='End.bp',
            score_col='tau',
            log2=True)
        gotten_cr = [
            np.mean([copy_ratio_dict[gene] for gene in window])
            for window in list_of_gene_windows
        ]
        gotten_wme, gotten_wme_metadata = wme.get_windowed_mean_expression(loom,
            list_of_gene_windows,
            patient=args_patient,
            patient_column='patient_ID',
            upper_cut=0)

        ps = []
        ps_tumor = []
        ps_tumor_low_complexity = []

        ps_nontumor = []
        ps_nontumor_low_complexity = []
        sorted_cr = np.argsort(gotten_cr)
        #    sorted_cr = np.hstack((sorted_cr[0:len(sorted_cr) // 5],
        #                           sorted_cr[4 * len(sorted_cr) // 5:len(sorted_cr)]))

        celltypes = metadata[metadata['patient_ID'] == args_patient]['cell_type'].values
        complexities = metadata[metadata['patient_ID'] == args_patient]['complexity'].values
        tps = metadata[metadata['patient_ID'] == args_patient]['time_point'].values
        for i in tqdm(range(gotten_wme.shape[1])):
#            bah = mk.original_test(gotten_wme[:, i][sorted_cr])
#            p = (1 - stats.norm.cdf(bah.z))
            p, a, b, c = theilslopes(gotten_wme[:,i], gotten_cr)
            if complexities[i] > 500:
                if celltypes[i] == 'tumor':
                    ps_tumor.append(p)
                else:
                    ps_nontumor.append(p)
            else:
                if celltypes[i] == 'tumor':
                    ps_tumor_low_complexity.append(p)
                else:
                    ps_nontumor_low_complexity.append(p)
#                        if complexities[i] > 500:
#                            ps_tumor.append(p)
#                        else:
#                            ps_tumor_low_complexity.append(p)
#        celltypes = metadata.query(
#            'patient_ID == {}'.format(args_patient))['cell_type'].values
#        complexities = metadata.query(
#            'patient_ID == {}'.format(args_patient))['complexity'].values
#        tps = metadata.query(
#            'patient_ID == {}'.format(args_patient))['time_point'].values
#        iterator = 0
#        for i in tqdm(range(gotten_wme.shape[1])):
#            if tps[i] == args_time_point:
#                if celltypes[i] != 'tumor':
#                    bah = mk.original_test(gotten_wme[:, i][sorted_cr])
#                    p = (1 - stats.norm.cdf(bah.z))
#                    p = bah.z
#                    #if p > 1:
#                    #    plt.plot(gotten_wme[:,i][sorted_cr])
#                    #    plt.show()
#                    ps.append(p)
#                    if celltypes[i] == 'tumor':
#                        ps_tumor.append(p)
#                    else:
#                        ps_nontumor.append(p)
#                    if iterator == 500:
#                        break
#                    iterator += 1
#        iterator = 0
#        sorted_cr = np.argsort(gotten_cr)
##        celltypes = metadata.query(
##            'patient_ID == {}'.format(args_patient))['cell_type'].values
#        for i in tqdm(range(gotten_wme.shape[1])):
#            if tps[i] == args_time_point:
#                if celltypes[i] == 'tumor':
#                    bah = mk.original_test(gotten_wme[:, i][sorted_cr])
#                    p = (1 - stats.norm.cdf(bah.z))
#                    p = bah.z
#                    ps.append(p)
#                    if celltypes[i] == 'tumor':
#                        #    if p < 3:
#                        #       plt.plot(gotten_wme[:,i][sorted_cr])
#                        #       plt.show()
#                        if complexities[i] > 1000:
#                            ps_tumor.append(p)
#                        else:
#                            ps_tumor_low_complexity.append(p)
#                    else:
#                        ps_nontumor.append(p)
#                    if iterator > 500 and len(ps_tumor_low_complexity) > 0:
#                        break
#                    iterator += 1
        sns.distplot(ps_nontumor, label='non-malignant', color='r')
        sns.distplot(ps_nontumor_low_complexity, label='low complexity non-malignant', color='y')
        sns.distplot(ps_tumor, label='high complexity_tumor', color='k')
        sns.distplot(ps_tumor_low_complexity, label='low complexity tumor', color='b')
        plt.legend()
        #from IPython.core.debugger import set_trace; set_trace()
        # The IPython is not in the dependencies and the output is disabled from Feb 2020.
        # Ipython debug mode should be removed.
        #if output:

        warn(f"Output is disabled from 20 Feb 2020. Saving to {output}")

        plt.savefig(output)

        #if False:
        #    plt.savefig(output)
        #else:
        #    plt.show()
