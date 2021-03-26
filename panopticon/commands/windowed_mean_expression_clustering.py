#! /usr/bin/env python

import click
import os
import argparse
from warnings import warn

import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib.pyplot as plt
import umap

from panopticon import wme
from panopticon.clustering import kt_cluster

# import loompy after the main packages, because sometimes it breaks packages that are imported further:
import loompy


def windowed_mean_expression_clustering_main(loomfile, patient, cell_type, complexity_cutoff, n_clusters, figure_output, raw_data_output):
    """

    Parameters
    ----------
    loomfile :
        
    patient :
        
    cell_type :
        
    complexity_cutoff :
        
    n_clusters :
        
    figure_output :
        

    Returns
    -------

    """
    print(patient)
    with loompy.connect(loomfile, validate=False) as loom:
    
        genes = loom.ra['gene']
        metadata = pd.DataFrame(loom.ca['patient_ID']).astype(str)
        metadata.columns = ['patient_ID']
        metadata['complexity'] = loom.ca['complexity']
        metadata['cell_type'] = loom.ca['cell_type']
        list_of_gene_windows = wme.get_list_of_gene_windows(genes)
        mean_window_expressions, mean_window_metadata = wme.get_windowed_mean_expression(loom,
            list_of_gene_windows,
            patient_column='patient_ID',
            patient=patient,
            cell_type_column='cell_type',
            cell_type=cell_type,
            complexity_cutoff=complexity_cutoff,
            complexity_column='complexity')
        mean_window_expression_ranks = wme.get_ranks(mean_window_expressions)
        reducer = umap.UMAP(random_state=17)
        reducer.fit(mean_window_expression_ranks.T)
        embedding = reducer.transform(mean_window_expression_ranks.T)
        if n_clusters > 1:
            labels, Z = kt_cluster(mean_window_expression_ranks, t=n_clusters)
        else:
            labels = np.ones(len(embedding))
        print(set(labels))
        fig, ax = plt.subplots()
        plt.rcParams.update({'font.size': 22})
        for label in set(labels):
            mask = labels == label
            plt.scatter(embedding[mask, 0], embedding[mask, 1], label=label)
        plt.title(" ".join(patient).replace('CSF','P').replace('DFCI','P'))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(top=False)
        ax.tick_params(labeltop=False)
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        ax.set_xlabel('UMAP 1', fontsize=22)
        ax.set_ylabel('UMAP 2', fontsize=22)
#        plt.legend()
        if raw_data_output:
            raw_data = pd.DataFrame(embedding)
            raw_data.columns = ['UMAP 1','UMAP 2']
            if raw_data_output.endswith('tsv'):
                raw_data.to_csv(raw_data_output,sep='\t')
            else:
                raw_data.to_csv(raw_data_output)
        if figure_output:
            plt.savefig(figure_output)
        else:
            plt.show()
        warn("Save cluster functionality temporarily eliminated")
    
#        clusters_save_choice = click.prompt(
#            'Would you like to save these clusters?',
#            type=click.Choice({'y', 'n'}, case_sensitive=False),
#            default='n')
#        if clusters_save_choice == 'y':
#            default_cluster_file = '/'.join(
#                data.split('/')[0:-1]) + '/clusters.txt'
#            cluster_file = click.prompt(
#                'Where should the cluster labels be saved?',
#                default=default_cluster_file)
#            np.savetxt(cluster_file, labels, fmt='%i')
