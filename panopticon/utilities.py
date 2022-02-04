import sys
from typing import List, Tuple

import numpy as np
import pandas as pd


def get_valid_gene_info(
    genes: List[str],
    release=102,
    species='homo sapiens'
) -> Tuple[List[str], List[int], List[int], List[int]]:
    """Returns gene locations for all genes in ensembl release 93  --S Markson 3 June 2020

    Parameters
    ----------
    genes : A list of genes
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes : List[str] :
        
    genes: List[str] :
        

    Returns
    -------

    
    """
    from panopticon.utilities import import_check
    exit_code = import_check("pyensembl", 'pip install pyensembl')
    if exit_code != 0:
        return
    from pyensembl import EnsemblRelease
    assembly = EnsemblRelease(release, species=species)
    gene_names = []
    gene_contigs = []
    gene_starts = []
    gene_ends = []
    for gene in np.intersect1d(genes, [
            gene.gene_name for gene in assembly.genes()
            if gene.contig.isnumeric() or gene.contig == 'X'
    ]):  # Toss genes not in hg38 release 93
        gene_info = assembly.genes_by_name(gene)
        gene_info = gene_info[0]
        gene_names.append(gene)
        gene_contigs.append(gene_info.contig)
        gene_starts.append(gene_info.start)
        gene_ends.append(gene_info.end)
    return gene_names, gene_contigs, gene_starts, gene_ends


def seurat_to_loom(seuratrds, patient_id_column, celltype_column,
                   complexity_column, loomfile):
    """

    Parameters
    ----------
    seuratrds :
        
    patient_id_column :
        
    celltype_column :
        
    complexity_column :
        
    loomfile :
        

    Returns
    -------

    
    """
    import rpy2.robjects as robjects
    from scipy import sparse
    from rpy2.robjects import pandas2ri
    import loompy
    robjects.r('''
    library(Seurat)
    seurat2rawandmeta <- function(seuratrds) {
        seuratobj <- readRDS(seuratrds)
        return(list(genes=rownames(seuratobj@data), metadata=seuratobj@meta.data, data=as.data.frame(summary(seuratobj@data))))
    }
    ''')
    seurat_grab = robjects.r['seurat2rawandmeta'](seuratrds)
    genes = pd.DataFrame(np.array(seurat_grab.rx2('genes')))
    genes.columns = ['gene']
    metadata = pandas2ri.rpy2py_dataframe(seurat_grab.rx2('metadata'))

    if patient_id_column != 'patient_ID':
        metadata['patient_ID'] = metadata[patient_id_column]
        metadata.drop(patient_id_column, inplace=True)
    if celltype_column != 'cell_type':
        metadata['cell_type'] = metadata[celltype_column]
        metadata.drop(celltype_column, inplace=True)
    if complexity_column != 'complexity':
        metadata['complexity'] = metadata[complexity_column]
        metadata.drop(complexity_column, inplace=True)
    data_df = pandas2ri.rpy2py_dataframe(seurat_grab.rx2('data'))
    sparsedata = sparse.coo_matrix(
        (data_df['x'], (data_df['i'] - 1, data_df['j'] - 1))).tocsc()
    sparsedata.resize((genes.shape[0], metadata.shape[0]))

    loompy.create(loomfile, sparsedata, genes.to_dict("list"),
                  metadata.to_dict("list"))


def intify(df_init):
    """

    Parameters
    ----------
    df_init :
        

    Returns
    -------

    
    """
    import binascii
    df = df_init.copy()
    for col in df.columns:
        if col.endswith('_ad'):
            raise Exception(
                "Don't append you column names with _ad! -- Samuel")
        df[col] = df[col].apply(
            lambda x: int(binascii.hexlify(x.encode()), 16))
    while np.sum(df.max() > sys.maxsize) > 0:
        for col in df.columns:
            if df[col].max() > sys.maxsize:
                df[col + '_ad'] = df[col] // sys.maxsize
                df[col] = df[col] % sys.maxsize
    return df.astype(np.int64)


def deintify(df_init):
    """

    Parameters
    ----------
    df_init :
        

    Returns
    -------

    
    """
    import binascii
    df = df_init.copy()
    while np.sum([x.endswith('_ad') for x in df.columns]) > 0:
        for col in df.columns:
            if col.endswith('_ad') and col + '_ad' not in df.columns:
                df[col[0:-3]] = df[col[0:-3]].astype(object)
                df[col] = df[col].astype(object)
                df[col[0:-3]] = df[col[0:-3]] + sys.maxsize * df[col]

                df.drop(col, axis=1, inplace=True)

    for col in df.columns:
        try:
            df[col] = df[col].apply(
                lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode())
        except:
            print(df[col].apply(
                lambda x: binascii.unhexlify(hex(x)[2::].encode()).decode()))
            raise Exception("whoops")
    return df


def recover_meta(db, do_deint=False):
    """

    Parameters
    ----------
    db :
        
    do_deint :
        (Default value = False)

    Returns
    -------

    
    """
    colmeta = None
    for key in db.ca.keys():
        if colmeta is None:
            colmeta = pd.DataFrame(db.ca[key])
            colmeta.columns = [key]
        else:
            colmeta[key] = db.ca[key]
    if do_deint:
        colmeta = deintify(colmeta.astype(np.int64))
    rowmeta = None
    for key in db.ra.keys():
        if rowmeta is None:
            rowmeta = pd.DataFrame(db.ra[key])
            rowmeta.columns = [key]
        else:
            rowmeta[key] = db.ra[key]
    if do_deint:
        rowmeta = deintify(rowmeta.astype(np.int64))
    return rowmeta, colmeta


def we_can_pickle_it(thing, thingname: str):
    """

    Parameters
    ----------
    thing :
        
    thingname : str :
        
    thingname : str :
        
    thingname : str :
        
    thingname : str :
        
    thingname: str :
        

    Returns
    -------

    
    """
    import pickle
    with open(thingname, 'wb') as f:
        pickle.dump(thing, f, pickle.HIGHEST_PROTOCOL)


def we_can_unpickle_it(thingname: str):
    """

    Parameters
    ----------
    thingname : str :
        
    thingname : str :
        
    thingname : str :
        
    thingname : str :
        
    thingname: str :
        

    Returns
    -------

    
    """
    import pickle
    with open(thingname, 'rb') as f:
        thing = pickle.load(f)
    return thing


def get_alpha_concave_hull_polygon(xcoords, ycoords, alpha=0.1, buffer=1):
    """Much credit to https://thehumangeo.wordpress.com/2014/05/12/drawing-boundaries-in-python/

    Parameters
    ----------
    xcoords :
        
    ycoords :
        
    alpha :
         (Default value = 0.1)
    buffer :
         (Default value = 1)

    Returns
    -------

    """

    from shapely.ops import cascaded_union, polygonize
    import shapely.geometry as geometry
    from scipy.spatial import Delaunay
    import numpy as np
    import math

    def alpha_shape(points, alpha):
        """Compute the alpha shape (concave hull) of a set
        of points.

        Parameters
        ----------
        points :
            Iterable container of points.
        alpha :
            alpha value to influence the
            gooeyness of the border. Smaller numbers
            don't fall inward as much as larger numbers.
            Too large, and you lose everything!

        Returns
        -------

        """
        if len(points) < 4:
            # When you have a triangle, there is no sense
            # in computing an alpha shape.
            return geometry.MultiPoint(list(points)).convex_hull

        def add_edge(edges, edge_points, coords, i, j):
            """Add a line between the i-th and j-th points,
            if not in the list already

            Parameters
            ----------
            edges :
                
            edge_points :
                
            coords :
                
            i :
                
            j :
                

            Returns
            -------

            """
            if (i, j) in edges or (j, i) in edges:
                # already added
                return
            edges.add((i, j))
            edge_points.append(coords[[i, j]])

        coords = np.array([point.coords[0] for point in points])

        tri = Delaunay(coords)
        edges = set()
        edge_points = []
        # loop over triangles:
        # ia, ib, ic = indices of corner points of the
        # triangle
        for ia, ib, ic in tri.vertices:
            pa = coords[ia]
            pb = coords[ib]
            pc = coords[ic]

            # Lengths of sides of triangle
            a = math.sqrt((pa[0] - pb[0])**2 + (pa[1] - pb[1])**2)
            b = math.sqrt((pb[0] - pc[0])**2 + (pb[1] - pc[1])**2)
            c = math.sqrt((pc[0] - pa[0])**2 + (pc[1] - pa[1])**2)

            # Semiperimeter of triangle
            s = (a + b + c) / 2.0

            # Area of triangle by Heron's formula
            area = math.sqrt(s * (s - a) * (s - b) * (s - c))
            circum_r = a * b * c / (4.0 * area)

            # Here's the radius filter.
            #print circum_r
            if circum_r < 1.0 / alpha:
                add_edge(edges, edge_points, coords, ia, ib)
                add_edge(edges, edge_points, coords, ib, ic)
                add_edge(edges, edge_points, coords, ic, ia)

        m = geometry.MultiLineString(edge_points)
        triangles = list(polygonize(m))
        return cascaded_union(triangles), edge_points

    points = []
    for x, y in zip(xcoords, ycoords):
        points.append(geometry.shape({'type': 'Point', 'coordinates': [x, y]}))

    concave_hull, edge_points = alpha_shape(points, alpha=alpha)
    return concave_hull.buffer(buffer)


def get_outlier_removal_mask(xcoords, ycoords, nth_neighbor=10, quantile=.9):
    """

    Parameters
    ----------
    xcoords :
        
    ycoords :
        
    nth_neighbor :
         (Default value = 10)
    quantile :
         (Default value = .9)

    Returns
    -------

    """
    from scipy.spatial.distance import pdist, squareform
    D = squareform(pdist(np.vstack((xcoords, ycoords)).T))
    distances = D[np.argsort(D, axis=0)[nth_neighbor - 1, :], 0]
    return distances <= np.quantile(distances, quantile)


def cohensd(g1, g2):
    """
    Returns Cohen's D for the effect size of group 1 values (g1) over group 2 values (g2).  

    Parameters
    ----------
    g1 : group 1 values (list or numpy vector)
        
    g2 : group 2 values (list or numpy vector)
        

    Returns
    -------
    (mean(g1) - mean(g2) )/s, where s is the pooled standard deviation of the two groups with Bessel's correction

    """
    n1 = len(g1)
    n2 = len(g2)
    s1 = np.std(g1, ddof=1)
    s2 = np.std(g2, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 * s1 + (n2 - 1) * s2 * s2) / (n1 + n2 - 2))

    return (np.mean(g1) - np.mean(g2)) / s


def phi_coefficient(contingency_table):
    """
    Returns the phi-coefficient for a contingency table.

    Paramenters
    -----------
    contingency_table : contingency table, identical in format to scipy.stats.fisher_exact


    Returns
    -------
    phi coefficient

    """
    table1 = contingency_table[0]
    table2 = contingency_table[1]
    table = np.vstack([table1, table2])
    phitop = (table1[0] * table2[1] - table1[1] * table2[0])
    phibottom = np.sqrt((table2[1]+table2[0])*\
                        (table1[1]+table1[0])*\
                        (table1[0]+table2[0])*\
                        (table2[1]+table1[1]))
    phi = phitop / phibottom
    return phi


def get_igraph_from_adjacency(adjacency, directed=None):
    """This is taken from scanpy._utils.__init__.py as of 12 August 2021
    
    
    Get igraph graph from adjacency matrix."""
    import igraph as ig

    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=directed)
    g.add_vertices(adjacency.shape[0])  # this adds adjacency.shape[0] vertices
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es['weight'] = weights
    except KeyError:
        pass
    if g.vcount() != adjacency.shape[0]:
        logg.warning(f'The constructed graph has only {g.vcount()} nodes. '
                     'Your adjacency matrix contained redundant nodes.')
    return g


def convert_10x_h5(path_10x_h5,
                   output_file,
                   labelkey=None,
                   label='',
                   genes_as_ca=[],
                   gene_whitelist=None,
                   output_type='loom'):

    import loompy
    from panopticon.utilities import import_check
    exit_code = import_check("cellranger", 'conda install -c hcc cellranger')
    if exit_code != 0:
        return
    import cellranger.matrix as cr_matrix
    output_type = output_file.split('.')[-1]
    if output_type not in ['loom', 'pkl']:
        raise Exception(
            "output_file must be have suffix loom or pkl, denoting an output type of loom of pickle respectively"
        )

    filtered_feature_bc_matrix = cr_matrix.CountMatrix.load_h5_file(
        path_10x_h5)
    id2feature = {
        val: key
        for key, val in filtered_feature_bc_matrix.feature_ids_map.items()
    }
    features = [
        id2feature[x].decode("utf-8")
        for x in range(filtered_feature_bc_matrix.features_dim)
    ]
    features_common_names = filtered_feature_bc_matrix.feature_ref.get_feature_names(
    )

    barcodes = filtered_feature_bc_matrix.bcs.astype(str)
    ca = {'cellname': barcodes}
    if labelkey is not None:
        ca[labelkey] = [label] * len(barcodes)

    m = filtered_feature_bc_matrix.m
    if gene_whitelist is not None:
        if len(gene_whitelist) > 0:
            mask = np.isin(features, gene_whitelist)
            m = m[mask, :]
            features = list(np.array(features)[mask])
            features_common_names = list(np.array(features_common_names)[mask])
    if type(genes_as_ca) == str:
        genes_as_ca = [genes_as_ca]
    else:
        genes_as_ca = list(genes_as_ca)
    if len(genes_as_ca) > 0:
        mask = np.isin(features, genes_as_ca)
        if len(genes_as_ca) != mask.sum():
            raise Exception(
                "Improper mapping of row attributes; perhaps gene of interest not in loom.ra[\'gene\']?"
            )
        for gene in genes_as_ca:
            submask = np.array(features) == gene
            gene_common_name = np.array(features_common_names)[submask][0]
            if np.sum(submask) > 1:
                raise Exception("Two or more features with this name")
            elif np.sum(submask) == 0:
                raise Exception("No features with this name")
            ca[gene] = list(m[submask, :].toarray()[0])
            ca[gene_common_name] = ca[gene]
        m = m[~mask, :]
        features = list(np.array(features)[~mask])
        features_common_names = list(np.array(features_common_names)[~mask])

    ra = {'gene': features, 'gene_common_name': features_common_names}
    if output_type == 'loom':
        loompy.create(output_file, m, ra, ca)
    if output_type == 'pkl':
        if gene_whitelist is None:
            raise Exception(
                "pkl output intended only for saving a small subsetted geneset of interest.  Please select a whitelist before saving as dataframe pkl."
            )
        mask = np.isin(features, gene_whitelist)
        features = np.array(features)[mask]
        features_common_names = np.array(features_common_names)[mask]
        df = pd.DataFrame(m[mask, :].toarray())
        df.index = features
        if labelkey is not None:
            df.columns = [labelkey + '_' + x for x in barcodes]
        else:
            df.columns = barcodes
        df.to_pickle(output_file)


def create_split_exon_gtf(input_gtf, output_gtf, gene):
    gtf = pd.read_table(input_gtf, header=None, comment='#')
    gtf.columns = [
        'seqname', 'source', 'feature', 'start', 'end', 'score', 'strand',
        'frame', 'attribute'
    ]
    gtf = gtf[gtf['feature'] == 'exon']
    if type(gene) == str:
        mask = gtf['attribute'].apply(
            lambda x: 'gene_name "{}"'.format(gene) in x)
    elif type(gene) in [list, tuple, np.array]:
        mask = np.array([False] * len(gtf))
        for g in gene:
            mask = mask | gtf['attribute'].apply(
                lambda x: 'gene_name "{}"'.format(g) in x)
    gtf_unchanged = gtf[~mask]
    gtf_changed = gtf[mask]

    def append_exon_number_to_id_and_name(attribute):
        exon_number = attribute.split('exon_number')[1].split(';')[0].split(
            '\"')[-2]

        old_gene_id_str = 'gene_id' + attribute.split('gene_id')[1].split(
            ';')[0]
        new_gene_id_str = '\"'.join(
            old_gene_id_str.split('\"')[0:-1]) + '-exon' + exon_number + '\"'
        attribute = attribute.replace(old_gene_id_str, new_gene_id_str)

        old_gene_name_str = 'gene_name' + attribute.split(
            'gene_name')[1].split(';')[0]
        new_gene_name_str = '\"'.join(
            old_gene_name_str.split('\"')[0:-1]) + '-exon' + exon_number + '\"'
        attribute = attribute.replace(old_gene_name_str, new_gene_name_str)

        old_transcript_id_str = 'transcript_id' + attribute.split(
            'transcript_id')[1].split(';')[0]
        new_transcript_id_str = '\"'.join(
            old_transcript_id_str.split('\"')
            [0:-1]) + '-exon' + exon_number + '\"'
        attribute = attribute.replace(old_transcript_id_str,
                                      new_transcript_id_str)

        old_transcript_name_str = 'transcript_name' + attribute.split(
            'transcript_name')[1].split(';')[0]
        new_transcript_name_str = '\"'.join(
            old_transcript_name_str.split('\"')
            [0:-1]) + '-exon' + exon_number + '\"'
        attribute = attribute.replace(old_transcript_name_str,
                                      new_transcript_name_str)

        if 'ccds_id' in attribute:
            old_ccds_id_str = 'ccds_id' + attribute.split('ccds_id')[1].split(
                ';')[0]
            new_ccds_id_str = '\"'.join(old_ccds_id_str.split('\"')
                                        [0:-1]) + '-exon' + exon_number + '\"'
            attribute = attribute.replace(old_ccds_id_str, new_ccds_id_str)

        return attribute

    gtf_changed['attribute'] = gtf_changed['attribute'].apply(
        append_exon_number_to_id_and_name)
    gtf = pd.concat([gtf_changed, gtf_unchanged])
    gtf.to_csv(output_gtf, sep='\t', index=False, header=None)


def get_umap_from_matrix(X,
                         random_state=17,
                         verbose=True,
                         min_dist=0.001,
                         n_neighbors=20,
                         metric='correlation'):

    import umap

    reducer = umap.UMAP(random_state=random_state,
                        verbose=verbose,
                        min_dist=min_dist,
                        n_neighbors=n_neighbors,
                        metric=metric)
    return reducer.fit_transform(X)


def convert_h5ad(h5ad,
                 output_loom,
                 convert_obsm=True,
                 convert_varm=True,
                 convert_uns=True,
                 convert_layers=True):
    import scanpy
    import loompy

    h5ad = scanpy.read_h5ad(h5ad)
    ra = {'gene': np.array(h5ad.var.index)}
    for col in h5ad.var.columns:
        if col == 'gene':
            raise Exception(
                "var column of h5ad is \"gene\".  This conflicts with panopticon loom format.  You must rename before converting."
            )
        else:
            ra[col] = np.array(h5ad.var[col].values)

    ca = {'cellname': np.array(h5ad.obs.index)}
    for col in h5ad.obs.columns:
        if col == 'cellname':
            raise Exception(
                "obs column of h5ad is \"cellname\".  This conflicts with panopticon loom format.  You must rename before converting."
            )
        else:
            ca[col] = np.array(h5ad.obs[col].values)
    if convert_obsm:
        for obsm_key in h5ad.obsm.keys():
            for i in range(h5ad.obsm[obsm_key].shape[1]):
                ca_key = "{}_{}".format(
                    obsm_key,
                    i + 1)  # one added so that these are 1-indexed by default
                if ca_key in ca.keys():
                    raise Exception(
                        "key\"{}\" already present as column attribute key.  Please rename to avoid."
                    )
                else:
                    ca[ca_key] = h5ad.obsm[obsm_key][:, i]
    if convert_varm:
        for varm_key in h5ad.varm.keys():
            for i in range(h5ad.varm[varm_key].shape[1]):
                ra_key = "{}_{}".format(
                    varm_key,
                    i + 1)  # one added so that these are 1-indexed by default
                if ra_key in ra.keys():
                    raise Exception(
                        "key\"{}\" already present as row attribute key.  Please rename to avoid."
                    )
                else:
                    ra[ra_key] = h5ad.varm[varm_key][:, i]
    loompy.create(output_loom, h5ad.X.T, ra, ca)

    if convert_uns:
        loom = loompy.connect(output_loom)
        for uns_key in h5ad.uns.keys():
            loom.attrs[uns_key] = h5ad.uns[uns_key]
        loom.close()

    if convert_layers:
        loom = loompy.connect(output_loom)
        for layer_key in h5ad.layers.keys():
            loom.layers[layer_key] = h5ad.layers[key].T
        loom.close()


def get_UMI_curve_from_10x_h5(path_10x_h5, save_to_file=None):
    from panopticon.utilities import import_check
    exit_code = import_check("cellranger", 'conda install -c hcc cellranger')
    if exit_code != 0:
        return
    import cellranger.matrix as cr_matrix
    import matplotlib.pyplot as plt

    bc_matrix = cr_matrix.CountMatrix.load_h5_file(path_10x_h5)
    fig, ax = plt.subplots(figsize=(5, 5))

    ax.plot(np.sort(bc_matrix.get_counts_per_bc())[::-1])
    ax.set_title('UMI counts per barcode, sorted')
    ax.set_ylabel('UMI counts')
    ax.set_xlabel('cell rank, UMI counts (most to fewest)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    if save_to_file is None:
        plt.show()
    else:
        plt.savefig(save_to_file)
        plt.cla()


def get_dsb_normalization(cell_antibody_counts,
                          empty_droplet_antibody_counts,
                          use_isotype_control=True,
                          denoise_counts=True,
                          isotype_control_name_vec=None,
                          define_pseudocount=False,
                          pseudocount_use=10,
                          quantile_clipping=False,
                          quantile_clip=[0.001, 0.9995],
                          return_stats=False):

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    if isotype_control_name_vec is None:
        isotype_control_name_vec = robjects.r("NULL")
    if (pseudocount_use != 10) and (not define_pseudocount):
        raise Exception(
            "\"define_pseudocount\" must be set to True to use pseudocount_use"
        )

    rpy2.robjects.numpy2ri.activate()

    robjects.r('''                                      
    library(mclust)                                     
    library(dsb)                                        
                                                        
    dsb <- function(cells,                              
                    empty,                              
                    use.isotype.control=TRUE,           
                    denoise.counts=TRUE,                
                    isotype.control.name.vec = NULL,    
                    define.pseudocount = FALSE,         
                    pseudocount.use = 10,               
                    quantile.clipping = FALSE,          
                    quantile.clip = c(0.001, 0.9995),   
                    return.stats = FALSE){
    
    DSBNormalizeProtein(cells, empty, use.isotype.control=use.isotype.control,
                    isotype.control.name.vec = isotype.control.name.vec,
                    denoise.counts=denoise.counts,
                    define.pseudocount = define.pseudocount, 
                    pseudocount.use = pseudocount.use, 
                    quantile.clipping = quantile.clipping,
                    quantile.clip = quantile.clip, 
                    return.stats = return.stats)
    }
    ''')
    dsb = robjects.r['dsb']
    return dsb(cell_antibody_counts,
               empty_droplet_antibody_counts,
               use_isotype_control=use_isotype_control,
               denoise_counts=denoise_counts,
               isotype_control_name_vec=isotype_control_name_vec,
               define_pseudocount=define_pseudocount,
               pseudocount_use=pseudocount_use,
               quantile_clipping=quantile_clipping,
               quantile_clip=quantile_clip,
               return_stats=return_stats)


def get_cellphonedb_compatible_counts_and_meta(loom,
                                               layername,
                                               celltype_ca,
                                               gene_ra='gene',
                                               cellname_ca='cellname',
                                               return_df=False,
                                               output_prefix=None,
                                               mouse_to_human=False):
    if output_prefix is None and not return_df:
        raise Exception(
            "either output_prefix must be specified, or return_df must be True"
        )

    counts = pd.DataFrame(loom[layername][:, :])
    counts.columns = loom.ca[cellname_ca]

    #counts.insert(0, 'Gene', np.array([x.upper() for x in loom.ra[gene_ra]]))
    genes = loom.ra[gene_ra]
    if mouse_to_human:
        from pybiomart import Server
        server = Server(host="http://www.ensembl.org")
        mouse_dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].
                         datasets['mmusculus_gene_ensembl'])
        mouse_data = mouse_dataset.query(
            attributes=['ensembl_gene_id', 'external_gene_name'])
        mouse_data['Gene upper'] = mouse_data['Gene name'].apply(
            lambda x: str(x).upper())
        human_dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].
                         datasets['hsapiens_gene_ensembl'])
        human_data = human_dataset.query(
            attributes=['ensembl_gene_id', 'external_gene_name'])
        conversion_dict = pd.merge(
            mouse_data, human_data, left_on='Gene upper',
            right_on='Gene name').set_index(
                'Gene stable ID_x')['Gene stable ID_y'].to_dict()

        convertible_mask = np.array(
            [x in conversion_dict.keys() for x in genes])
        genes = [
            conversion_dict[x] if x in conversion_dict.keys() else np.nan
            for x in genes
        ]
    counts.insert(0, 'Gene', genes)
    if mouse_to_human:
        counts = counts.iloc[convertible_mask, :]
        counts = counts.groupby('Gene').first().reset_index()
    meta = pd.DataFrame(loom.ca[cellname_ca])
    meta.columns = ['Cell']
    meta['cell_type'] = loom.ca[celltype_ca]

    if output_prefix is not None:
        counts.to_csv(output_prefix + '_counts.txt', sep='\t', index=False)
        meta.to_csv(output_prefix + '_meta.txt', sep='\t', index=False)
        command = 'cellphonedb method statistical_analysis {0}_meta.txt {0}_counts.txt'.format(
            output_prefix)
        print("Run cellphonedb on command line with \"{}\"".format(command))
    elif return_df:
        return meta, counts


def create_gsea_txt_and_cls(loom,
                            layername,
                            output_prefix,
                            phenotypes,
                            cellmask=None,
                            gene_ra='gene',
                            cellname_ca='cellname'):
    import os
    if cellmask is None:
        cellmask = np.array([True] * loom.shape[1])
    if type(phenotypes) == str:
        phenotypes = loom.ca[phenotypes]
    if len(phenotypes) != cellmask.sum():
        raise Exception(
            "length of phenotypes vector must be equal to number of samples (cells)"
        )

    txt = pd.DataFrame(loom.ra[gene_ra])
    txt.columns = ['NAME']
    txt['DESCRIPTION'] = 'na'
    #txt = pd.concat([txt,pd.DataFrame(loom[layername][:,cellmask])],axis=1)
    #txt.columns = ['NAME','DESCRIPTION'] + list(loom.ca[cellname_ca][cellmask])
    #txt.to_csv(output_prefix+'.txt',index=False,sep='\t')

    total = cellmask.sum()
    nphenotypes = len(np.unique(phenotypes))
    outcls = output_prefix + '.cls'
    if os.path.exists(outcls):
        os.system("rm {}".format(outcls))
        #raise Exception("cls file already present--cannot overwrite")
    line1 = "{} {} 1".format(total, nphenotypes)
    line2 = '# ' + ' '.join(np.unique(phenotypes))
    phenotype2index = {
        phenotype: i
        for i, phenotype in enumerate(np.unique(phenotypes))
    }
    #print(phenotype2index)
    #print([phenotype2index[x] for x in phenotypes])
    line3 = ' '.join([str(phenotype2index[x]) for x in phenotypes])
    for line in [line1, line2, line3]:
        os.system('echo \"{}\">>{}'.format(line, outcls))


def get_cross_column_attribute_heatmap(loom,
                                       ca1,
                                       ca2,
                                       normalization_axis=None):
    #if type(normalization_axis) == list:
    #    outdfs = []
    #    for axis in normalization_axis:
    #        outdfs.append(get_cross_column_attribute_heatmap(loom, ca1, ca2, normalization_axis=axis))
    #    return outdfs

    df = pd.DataFrame(loom.ca[ca1], copy=True)
    df.columns = [ca1]
    df[ca2] = loom.ca[ca2]

    df = pd.DataFrame(df.groupby(ca1, )[ca2].value_counts())
    df.columns = ['counts']

    dfs = []

    for i, df_group in df.reset_index().groupby(ca1):
        dfs.append(
            df_group.rename(columns={
                'counts': 'counts_' + i
            }).set_index(ca2)['counts_' + i])
    outdf = pd.concat(dfs, axis=1)
    if normalization_axis is None:
        return outdf
    elif normalization_axis == 0:
        return np.divide(outdf, outdf.sum(axis=0).values)
    elif normalization_axis == 1:
        return np.divide(outdf.T, outdf.sum(axis=1).values).T
    else:
        raise Exception("normalization axis must be one of \"None\", 0, or 1")


def get_complement_contigency_tables(df):
    if type(df) != pd.core.frame.DataFrame:
        raise Exception("pandas dataframe expected input")
    complement_contigency_table_dict = {}

    for col in df.columns:
        complement_contigency_table_dict[col] = {}

        for index in df.index.values:
            a = df.loc[index][col].sum()
            b = df.loc[index][[x for x in df.columns if x != col]].sum()
            c = df.loc[[x for x in df.index if x != index]][col].sum()
            d = np.sum(df.loc[[x for x in df.index if x != index
                               ]][[x for x in df.columns if x != col]].sum())

            complement_contigency_table_dict[col][index] = [[a, b], [c, d]]
    return complement_contigency_table_dict


def get_cluster_differential_expression_heatmap_df(loom,
                                                   layer,
                                                   clusteringlevel,
                                                   diffex={},
                                                   gene_name='gene',
                                                   cell_name='cellname'):
    """

    Returns
    -------

    """
    from panopticon.analysis import get_cluster_differential_expression
    import seaborn as sns
    import pandas as pd

    clusteredmask = []
    for cluster in np.unique(loom.ca[clusteringlevel]):
        mask = loom.ca[clusteringlevel] == cluster
        if mask.sum() > 2:
            clusteredmask.append(np.where(mask)[0])
    clusteredmask = np.hstack(clusteredmask)
    allgenes = []
    allgeneindices = []
    rawX = []
    clusters = [
        x for x in np.unique(loom.ca[clusteringlevel]) if x in diffex.keys()
    ]
    for cluster in clusters:
        mask = loom.ca[clusteringlevel] == cluster
        genes = diffex[cluster][~diffex[cluster]['gene'].isin(allgenes)].query(
            'MeanExpr1 > MeanExpr2').query('FracExpr2<.9').head(
                10)['gene'].values
        genemask = np.isin(loom.ra['gene'], genes)
        rawX.append(loom[layer][genemask, :][:, clusteredmask.nonzero()[0]])
        allgenes.append(genes)
        allgeneindices.append(np.where(genemask)[0])
    clusteredmask = np.hstack(clusteredmask)
    allgeneindices = np.hstack(allgeneindices)
    hmdf = pd.DataFrame(np.vstack(rawX))
    hmdf.index = np.hstack(loom.ra[gene_name][allgeneindices])
    hmdf.columns = loom.ca[cell_name][clusteredmask]
    return hmdf


def generate_ca_frequency(loom,
                          ca,
                          blacklisted_ca_values=[],
                          second_ca=None,
                          output_name=None,
                          overwrite=False):
    if output_name is None:
        raise Exception("output_name must be specified")
    if output_name in loom.ca.keys() and overwrite is False:
        raise Exception(
            "overwrite must be True to write over existing ca ({})".format(
                output_name))
    ca2counts = pd.DataFrame(loom.ca[ca])[0].value_counts().to_dict()
    for ca_value in blacklisted_ca_values:
        ca2counts[ca_value] = np.nan
    if second_ca is None:
        denominator = loom.shape[1]
        frequencies = [ca2counts[x] / denominator for x in loom.ca[ca]]
    else:
        denominator = pd.DataFrame(
            loom.ca[second_ca])[0].value_counts().to_dict()
        frequencies = [
            ca2counts[x] / denominator[y]
            for x, y in zip(loom.ca[ca], loom.ca[second_ca])
        ]

    loom.ca[output_name] = frequencies

def import_check(package, statement_upon_failure):
    try:
        exec("import {}".format(package))
        return 0
    except:
        print("Import of package \'{}\' failed. We recommend installing with the command \'{}\'".format(package, statement_upon_failure))
        return 1
