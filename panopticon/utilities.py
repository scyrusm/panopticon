import sys
from typing import List, Tuple

import numpy as np
import pandas as pd


def get_valid_gene_info(
        genes: List[str]) -> Tuple[List[str], List[int], List[int], List[int]]:
    """Returns gene locations for all genes in ensemble 93  --S Markson 3 June 2020

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
    from pyensembl import EnsemblRelease
    assembly = EnsemblRelease(93)
    valid_chromosomes = [str(x) for x in range(1, 23)] + ['X']
    gene_names = []
    gene_contigs = []
    gene_starts = []
    gene_ends = []
    for gene in np.intersect1d(genes, [
            gene.gene_name
            for gene in assembly.genes() if gene.contig in valid_chromosomes
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


def convert_10x_h5(path_10x_h5, output_loom, labelkey=None,label='', genes_as_ca=[]):
    import cellranger.matrix as cr_matrix
    import loompy
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
        ca[labelkey] = [label]*len(barcodes)


    m = filtered_feature_bc_matrix.m
    if type(genes_as_ca)==str:
        genes_as_ca = [genes_as_ca]
    if len(genes_as_ca)>0:
        mask = np.isin(features, genes_as_ca)
        if len(genes_as_ca) != mask.sum():
            raise Exception("Improper mapping of row attributes; perhaps gene of interest not in loom.ra[\'gene\']?")
        for gene in genes_as_ca:
            submask = features == gene 
            ca[gene] = list(m[submask,:].toarray()[0])
        m = m[~mask,:]
        features = list(np.array(features)[~mask])
        features_common_names = list(np.array(features_common_names)[~mask])

    ra = {'gene': features, 'gene_common_name': features_common_names}
    loompy.create(output_loom, m, ra, ca)

def create_split_exon_gtf(input_gtf,output_gtf,gene):
    gtf = pd.read_table(input_gtf,header=None, comment='#')
    gtf.columns = ['seqname','source','feature','start','end','score','strand','frame','attribute']
    gtf = gtf[gtf['feature']=='exon']
    if type(gene)==str:
        mask = gtf['attribute'].apply(lambda x: 'gene_name "{}"'.format(gene) in x)
    elif type(gene) in [list, tuple, np.array]:
        mask = np.array([False]*len(gtf))
        for g in gene:
            mask = mask | gtf['attribute'].apply(lambda x: 'gene_name "{}"'.format(g) in x)
    gtf_unchanged = gtf[~mask]
    gtf_changed = gtf[mask]

    def append_exon_number_to_id_and_name(attribute):
        exon_number = attribute.split('exon_number')[1].split(';')[0].split('\"')[-2]

        old_gene_id_str = 'gene_id'+attribute.split('gene_id')[1].split(';')[0]
        new_gene_id_str = '\"'.join(old_gene_id_str.split('\"')[0:-1])+'-exon'+exon_number+'\"'

        old_gene_name_str = 'gene_name'+attribute.split('gene_name')[1].split(';')[0]
        new_gene_name_str = '\"'.join(old_gene_name_str.split('\"')[0:-1])+'-exon'+exon_number+'\"'

        old_transcript_id_str = 'transcript_id'+attribute.split('transcript_id')[1].split(';')[0]
        new_transcript_id_str = '\"'.join(old_transcript_id_str.split('\"')[0:-1])+'-exon'+exon_number+'\"'

        old_transcript_name_str = 'transcript_name'+attribute.split('transcript_name')[1].split(';')[0]
        new_transcript_name_str = '\"'.join(old_transcript_name_str.split('\"')[0:-1])+'-exon'+exon_number+'\"'

        old_ccds_id_str = 'ccds_id'+attribute.split('ccds_id')[1].split(';')[0]
        new_ccds_id_str = '\"'.join(old_ccds_id_str.split('\"')[0:-1])+'-exon'+exon_number+'\"'

        attribute = attribute.replace(old_gene_id_str,new_gene_id_str)
        attribute = attribute.replace(old_gene_name_str,new_gene_name_str)
        attribute = attribute.replace(old_transcript_id_str,new_transcript_id_str)
        attribute = attribute.replace(old_transcript_name_str,new_transcript_name_str)
        attribute = attribute.replace(old_ccds_id_str,new_ccds_id_str)

        return attribute

    gtf_changed['attribute'] = gtf_changed['attribute'].apply(append_exon_number_to_id_and_name)
    gtf = pd.concat([gtf_changed, gtf_unchanged])
    gtf.to_csv(output_gtf,sep='\t',index=False,header=None)    
