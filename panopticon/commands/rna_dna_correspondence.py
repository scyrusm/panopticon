#! /usr/bin/env python
from tqdm import tqdm
import argparse
import pandas as pd
from panopticon import wme, dna
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import sparse
import pymannkendall as mk
import loompy


def rna_dna_correspondence_main(loomfile, args_seg_file, args_patient_column,
                                args_patient, args_time_point, output=None):
    with loompy.connect(loomfile, validate=False) as loom:
        genes = loom.ra['gene']
        metadata = pd.DataFrame(loom.ca['patient_ID'])
        metadata.columns = ['patient_ID']
        metadata['complexity'] = loom.ca['complexity']
        metadata['cell_type'] = loom.ca['cell_type']
        metadata['time_point'] = loom.ca['time_point']

        list_of_gene_windows = wme.get_list_of_gene_windows(
            genes, window_size=200, window_step=50)

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
        gotten_wme = wme.get_windowed_mean_expression(loom,
            list_of_gene_windows,
            patient=args_patient,
            patient_column='patient_ID',
            upper_cut=0)

        ps = []
        ps_tumor = []
        ps_tumor_low_complexity = []

        ps_nontumor = []
        sorted_cr = np.argsort(gotten_cr)
        #    sorted_cr = np.hstack((sorted_cr[0:len(sorted_cr) // 5],
        #                           sorted_cr[4 * len(sorted_cr) // 5:len(sorted_cr)]))
        celltypes = metadata.query(
            'patient_ID == {}'.format(args_patient))['cell_type'].values
        complexities = metadata.query(
            'patient_ID == {}'.format(args_patient))['complexity'].values
        tps = metadata.query(
            'patient_ID == {}'.format(args_patient))['time_point'].values
        iterator = 0
        for i in tqdm(range(gotten_wme.shape[1])):
            if tps[i] == args_time_point:
                if celltypes[i] != 'tumor':
                    bah = mk.original_test(gotten_wme[:, i][sorted_cr])
                    p = (1 - stats.norm.cdf(bah.z))
                    p = bah.z
                    #if p > 1:
                    #    plt.plot(gotten_wme[:,i][sorted_cr])
                    #    plt.show()
                    ps.append(p)
                    if celltypes[i] == 'tumor':
                        ps_tumor.append(p)
                    else:
                        ps_nontumor.append(p)
                    if iterator == 200:
                        break
                    iterator += 1
        iterator = 0
        sorted_cr = np.argsort(gotten_cr)
        celltypes = metadata.query(
            'patient_ID == {}'.format(args_patient))['cell_type'].values
        for i in tqdm(range(gotten_wme.shape[1])):
            if tps[i] == args_time_point:
                if celltypes[i] == 'tumor':
                    bah = mk.original_test(gotten_wme[:, i][sorted_cr])
                    p = (1 - stats.norm.cdf(bah.z))
                    p = bah.z
                    ps.append(p)
                    if celltypes[i] == 'tumor':
                        #    if p < 3:
                        #       plt.plot(gotten_wme[:,i][sorted_cr])
                        #       plt.show()
                        if complexities[i] > 1000:
                            ps_tumor.append(p)
                        else:
                            ps_tumor_low_complexity.append(p)
                    else:
                        ps_nontumor.append(p)
                    if iterator > 200 and len(ps_tumor_low_complexity) > 0:
                        break
                    iterator += 1
        sns.distplot(ps_nontumor)
        sns.distplot(ps_tumor)
        sns.distplot(ps_tumor_low_complexity)
        if output:
            plt.savefig(output)
        else:
            plt.show()
