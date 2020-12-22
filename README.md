# panopticon

![Utopian data structure](https://upload.wikimedia.org/wikipedia/en/e/e1/Panopticon_Willey_Reveley_1791.png )

# Overview

Panopticon is a set of tools named after Jeremy Bentham's model for a prison, whereby multiple cells (!) could be observed with ease by a single individual.

These tools are highly modular, and at the moment the most up-to-date functionality is bunged up in `panopticon/utilities.py`.  This [smells](https://en.wikipedia.org/wiki/Code_smell), I know.  

This package is designed to give me better performance and flexibility than other general single-cell RNA analysis tools, such as [Seurat](https://satijalab.org/seurat/) and [scanpy](https://scanpy.readthedocs.io/en/stable/).  I don't pretend that this tool is better than those, and if you're just starting out with single-cell RNA analysis, I recommend that you become familiar with those tools first.  

Design-wise, panopticon is a bit similar to scanpy, in that it relies heavily on loom files, and is python-based.  

# Installation

You can install the panopticon on the command line by cloning this git repository.  Because this package is in a rather alpha mode, I recommend installing with the `--editable` flag.
```
    git clone https://github.com/scyrusm/panopticon.git
    cd panopticon
    pip install --editable .
    pyensembl install --release 75 --species homo_sapiens
```

Now try exploring the exciting panopticon programming with

```
    panopticon --help
```

If you would like autocompletion of panopticon commands (though I don't do this myself), put the following line in your `.bashrc`
```
eval "$(_PANOPTICON_COMPLETE=source panopticon)"
```

# Getting started

### Making a panopticon-friendly loom file

To get started, find some scRNA data, or download some public data e.g. from the [Broad Institute Single Cell Portal](https://singlecell.broadinstitute.org/single_cell).  One good starting place is [Tirosh et al., Science 2016](https://singlecell.broadinstitute.org/single_cell/study/SCP11/melanoma-intra-tumor-heterogeneity).  Most of the time, scRNA data is shared as a big dense `.txt`, `.tsv` or `.csv` matrix with each column corresponding to a cell, and each row corresponding to a gene, and the entries representing gene expression, normalized in some way.  Of course, sometimes the rows are cells and columns are genes, so keep your head on a swivel.  Usually there will be some metadata wherein each row gives some extra info corresponding to each cell.   

The basic data object used by panopticon is a [loom](http://loompy.org/) file, with some expectations for what certain columns attributes and row attributes are named.  Cell names belong in the `cellname` column attribute, gene names belong in the `gene` row attribute; complexity (the number of unique genes expressed by a cell) is the column attribute `nGene`---the number of unique cells expressing a gene is `nCell`.  Beyond that, you're pretty vogelfrei.  

For help in getting from whatever data you have to what you can run other panopticon commands on, run
```
panopticon book scrna-wizard
```
This tool should guide you through a series of prompts to make a fully functioning panopticon-friendly loom file.  If it seems to be failing, it's usually because you have some input data with a formatting edge case that I haven't yet thought of.  Try to resolve that, or, if it persists, message me.  

### Data normalization, clustering, and exploration

9 times out of 10 I'll start out with something like the following, where the `''` layer represents raw transcript counts.  
```
import loompy
from panopticon.utilities import generate_incremental_pca, generate_embedding, generate_clustering, generate_count_normalization, generate_masked_module_score
import matplotlib.pyplot as plt
import numpy as np

db = loompy.connect("WhateverYouNamedIt.loom")

generate_count_normalization(db, '','log2(TP100k+1)')
generate_incremental_pca(db, 'log2(TP100k+1)')
generate_embedding(db, 'log2(TP100k+1)' )
generate_clustering(db, 'log2(TP100k+1)', clustering_depth=3)

fig, ax = plt.subplots(figsize=(6,6))
for cluster in np.unique(db.ca['ClusteringIteration0']):
    mask = db.ca['ClusteringIteration0'] == cluster
    ax.scatter(db.ca['PCA UMAP embedding 1'][mask], db.ca['PCA UMAP embedding 2'][mask], label=cluster)
plt.legend()
plt.show()

```
This will plot your basic UMAP of all your cells, with cells clusters and colored based on the first iteration of the panopticon's agglomerative iterative subclustering procedure.

# Additional Documentation

For more in-depth documentation, see panopticon's [Read the Docs](https://panopticon-single-cell.readthedocs.io) page.

# Development

## Development

If you would like to report a bug or request a feature, please do so using the [issues page](https://github.com/scyrusm/panopticon/issues)

## Dependencies
All code has been tested on python 3.7.7.  Package dependencies below.  

```
numpy
pandas
scipy
matplotlib
umap-learn
click
pyensembl
tqdm
seaborn
loompy
pymannkendall

```



# To do
- Clean up docstrings and make sure that they are in numpydoc format
- add in tests, including test data
