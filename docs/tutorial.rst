===============
Getting started
===============

Making a panopticon-friendly loom file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To get started, find some scRNA data, or download some public data e.g.
from the `Broad Institute Single Cell
Portal <https://singlecell.broadinstitute.org/single_cell>`__. One good
starting place is `Tirosh et al., Science
2016 <https://singlecell.broadinstitute.org/single_cell/study/SCP11/melanoma-intra-tumor-heterogeneity>`__.
Most of the time, scRNA data is shared as a big dense ``.txt``, ``.tsv``
or ``.csv`` matrix with each column corresponding to a cell, and each
row corresponding to a gene, and the entries representing gene
expression, normalized in some way. Of course, sometimes the rows are
cells and columns are genes, so keep your head on a swivel. Usually
there will be some metadata wherein each row gives some extra info
corresponding to each cell.

The basic data object used by panopticon is a
`loom <http://loompy.org/>`__ file, with some expectations for what
certain columns attributes and row attributes are named. Cell names
belong in the ``cellname`` column attribute, gene names belong in the
``gene`` row attribute; complexity (the number of unique genes expressed
by a cell) is the column attribute ``nGene``---the number of unique
cells expressing a gene is ``nCell``. Beyond that, you're pretty
vogelfrei.

For help in getting from whatever data you have to what you can run
other panopticon commands on, run

::

    panopticon book scrna-wizard

This tool should guide you through a series of prompts to make a fully
functioning panopticon-friendly loom file. If it seems to be failing,
it's usually because you have some input data with a formatting edge
case that I haven't yet thought of. Try to resolve that, or, if it
persists, message me.

Data normalization, clustering, and exploration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

9 times out of 10 I'll start out with something like the following,
where the ``''`` layer represents raw transcript counts.

::

    import loompy
    from panopticon.preprocessing import generate_count_normalization
    from panopticon.analysis import generate_incremental_pca, generate_embedding, generate_clustering, generate_masked_module_score
    import matplotlib.pyplot as plt
    import numpy as np

    db = loompy.connect("WhateverYouNamedIt.loom")

    layername = 'log2(TP100k+1)'
    # Generates a new layer with log2(transcripts per
    # denominator).
    generate_count_normalization(db, '', layername)

    generate_incremental_pca(db, layername)
    generate_embedding(db, layername)
    generate_clustering(db, layername, clustering_depth=3)

    fig, ax = plt.subplots(figsize=(6,6))
    for cluster in np.unique(db.ca['ClusteringIteration0']):
        mask = db.ca['ClusteringIteration0'] == cluster
        ax.scatter(db.ca['PCA UMAP embedding 1'][mask], db.ca['PCA UMAP embedding 2'][mask], label=cluster)
    plt.legend()
    plt.show()

This will plot your basic UMAP of all your cells, with cells clusters
and colored based on the first iteration of the panopticon's
agglomerative iterative subclustering procedure.

Making a split-exon gtf file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose that you have a ``.gtf`` file (hereafter ``nameofyourgtf.gtf``). If you want to create a different ``.gtf`` file where a particular gene (say, ``GeneOfInterest``) has been replaced with separate "genes" corresponding to different exons of that gene, you run the following command on the command line:
::
    panopticon create-split-exon-gtf nameofyourgtf.gtf nameofyouroutputgtf.gtf GeneOfInterest

If there are multiple such genes (``GeneOfInterest1``, ``GeneOfInterest2``...), you may "split" them with the command:
::
    panopticon create-split-exon-gtf nameofyourgtf.gtf nameofyouroutputgtf.gtf GeneOfInterest1 GeneOfInterest2

and so on. 
