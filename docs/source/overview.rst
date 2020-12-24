
========
Overview
========

.. figure:: https://upload.wikimedia.org/wikipedia/en/e/e1/Panopticon_Willey_Reveley_1791.png
   :alt: Utopian data structure

   Utopian data structure

Panopticon is a set of tools named after Jeremy Bentham's model for a
prison, whereby multiple cells (!) could be observed with ease by a
single individual.

This package is designed to give me some functionality I didn't find in other general single-cell RNA analysis tools, such as
`Seurat <https://satijalab.org/seurat/>`__ and
`scanpy <https://scanpy.readthedocs.io/en/stable/>`__. 

A few unique features include

* functions for joint analysis of bulk DNA and single-cell RNA

* convenient iterative subclustering (by default agglomerative clustering with a Pearson correlation or cosine metric, with number of clusters selected by silhouette score)

Panopticon relies on `loom <http://loompy.org/>`__ files with certain demands on column and row attributes.  
