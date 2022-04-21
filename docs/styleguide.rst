
======================
Panopticon Style Guide
======================

**Overview**

Rows are genes, columns are cells.  **Row** attributes therefore represent genes and by-gene metadata, including

- gene (can be HUGO Gene Nomenclature Committee names, or the Ensembl ID).
    - It is preferable to use the Ensembl ID as "gene." See `here <https://www.biostars.org/p/344244/>`_ for some discussion. When doing so, use an auxiliary row attribute for the common, HUGO names. We recommend "GeneCommonName".
- Components from dimensionality reduction; for example, "log2(TP10k+1) PC 1" would represent the loadings of genes of the first principal component of the layer "log2(TP10k+1)"
- chromosome


**Columns** represent cells, and by-cell metadata, including

- cellname (can be the cell barcode, or the cellbarcode + sample info)
- Loadings of dimensional reductions. For example, "log2(TP10k+1) PC 1 Loading" would represent the loadings of the first principal component of layer "log2(TP10k+1)" on cells. 
- Use column attributes for any other cell-based annotations, e.g. module scores, quality scores, TCR information, and by-hand cell annotations.  
- When adding an annotation, add a date in DD<latin letter month>YYYY format, e.g. "CellAssignment1Apr2022".

Loom files do not support matrix formats as column attributes, unlike anndata files, which allow "obsm" and "varm" attributes for columns and rows, respectively.  This can be addressed by breaking your matrix into separate vectors (e.g. "PCA UMAP embedding 1", "PCA UMAP embedding 2").  **It is the opinion of the author that this limitation is a good discipline anyway, and that discipline is freedom.** When breaking a matrix into vectors, use 1-indexing.

**Layers** are counts matrices.  Use them to store:

- raw counts
- normalized counts
- standardized counts
- spliced/unspliced counts

**Global attributes** should be used sparingly. 

The author **highly** recommends using file locking when working working with loom (or any hdf5) files. Therefore, consider adding the following to your .bashrc:

::

    export HDF5_USE_FILE_LOCKING=TRUE


