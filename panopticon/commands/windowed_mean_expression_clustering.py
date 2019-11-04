#! /usr/bin/env python

import argparse
from panopticon import wme
import umap
import numpy as np
import os
import pandas as pd
from scipy import sparse
from warnings import warn
import matplotlib.pyplot as plt
import click
from panopticon.clustering import kt_cluster
import loompy


def windowed_mean_expression_clustering_main(loomfile, patient, cell_type, complexity_cutoff, n_clusters, figure_output):
    print(patient)
    with loompy.connect(loomfile, validate=False) as loom:
    
        genes = loom.ra['gene']
        metadata = pd.DataFrame(loom.ca['patient_ID'])
        metadata.columns = ['patient_ID']
        metadata['complexity'] = loom.ca['complexity']
        metadata['cell_type'] = loom.ca['cell_type']
        list_of_gene_windows = wme.get_list_of_gene_windows(genes)
        mean_window_expressions = wme.get_windowed_mean_expression(loom,
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
        for label in set(labels):
            mask = labels == label
            plt.scatter(embedding[mask, 0], embedding[mask, 1], label=label)
        plt.legend()
        plt.title("UMAP Visualization")
        if figure_output:
            plt.savefig(figure_output)
        else:
            plt.show()
    
        clusters_save_choice = click.prompt(
            'Would you like to save these clusters?',
            type=click.Choice({'y', 'n'}, case_sensitive=False),
            default='n')
        if clusters_save_choice == 'y':
            default_cluster_file = '/'.join(
                data.split('/')[0:-1]) + '/clusters.txt'
            cluster_file = click.prompt(
                'Where should the cluster labels be saved?',
                default=default_cluster_file)
            np.savetxt(cluster_file, labels, fmt='%i')
