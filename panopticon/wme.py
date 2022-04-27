"""
wme.py
====================================
wme
"""

# second version
import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy import stats
from itertools import islice
from scipy.sparse import coo_matrix, save_npz
from panopticon.utilities import get_valid_gene_info


def get_list_of_gene_windows(genes,
                             window_size=200,
                             window_step=50,
                             release=106,
                             species='homo sapiens'):
    """
    This function will, given a set of genes, return a list of lists, where 

    Parameters
    ----------
    genes : list of str
        The list of genes that will be used to generate a list of gene windows (list of lists).
    window_step : int
        How many genes over each window will be "shifted" from the previous. (Default value = 50)
    window_size : int
        The size of the windows. (Default value = 200)
    release : int
        The ensembl release which will be used to sort the genes into windows of contiguous genes along the genome. (Default value = 106)

    Returns
    -------
    A list of lists of strings. Each element of this list will have length window_size. 

    """
    gene_names, gene_contigs, gene_starts, gene_ends = get_valid_gene_info(
        genes, release=release, species=species)

    gene_df = pd.DataFrame(gene_names)
    gene_df.columns = ['name']
    gene_df['contig'] = gene_contigs
    gene_df['start'] = gene_starts
    gene_df['end'] = gene_ends
    gene_df_groupby = gene_df.set_index('name').sort_values('start').groupby(
        'contig')
    list_of_gene_windows = []
    for chromosome in gene_df['contig'].unique():
        list_of_gene_windows += [
            list(gene_df_groupby.groups[chromosome])[i:(i + window_size)]
            for i in np.arange(
                0,
                len(gene_df_groupby.groups[chromosome]) - window_size +
                1, window_step)
        ]

    return list_of_gene_windows


def robust_mean_windowed_expressions(genes,
                                     list_of_gene_windows,
                                     expression_data,
                                     upper_cut=5,
                                     windsor=False,
                                     use_tqdm=True,
                                     tqdm_desc='computing WME'):
    """
    Produces an arithmetic mean over expression in windows determined by list_of_gene_windows.  Highest-expression genes in each window are discarded.  
    Can be made more memory-friendly, by implementing a map function over expression_data--I still haven't done this.  S Markson 4 June 2020.  

    Parameters
    ----------
    genes :
        param list_of_gene_windows:
    expression_data :
        param upper_cut:  (Default value = 0)
    windsor :
        Default value = False)
    tqdm_desc :
        Default value = '')
    list_of_gene_windows :
        
    upper_cut :
        (Default value = 5)

    Returns
    -------

    
    """
    gene_to_index = {gene: i for i, gene in enumerate(genes)}
    mean_window_expressions = np.zeros(
        (len(list_of_gene_windows), expression_data.shape[1]))
    with tqdm(total=len(list_of_gene_windows), desc=tqdm_desc,disable=~use_tqdm) as pbar:
        for i, window in enumerate(list_of_gene_windows):
            window_expression_indices = np.array(
                [gene_to_index[gene] for gene in window])
            exprs = expression_data[window_expression_indices, :]
            robust_cell_means = np.zeros(exprs.shape[1])
            for icell in range(exprs.shape[1]):
                cell_exprs = exprs[:, icell]
                truncated = np.sort(cell_exprs)[::-1][upper_cut::]
                if windsor:
                    robust_cell_means[icell] = np.hstack(
                        ([truncated[0]] * upper_cut, truncated)).mean()
                else:
                    robust_cell_means[icell] = truncated.mean()

            mean_window_expressions[i, :] = robust_cell_means
            pbar.update(1)
    return mean_window_expressions


def get_windowed_mean_expression(loom,
                                 list_of_gene_windows,
                                 patient_column='Patient_ID',
                                 patient=0,
                                 cell_type_column=None,
                                 cell_type=None,
                                 complexity_column='nGene',
                                 complexity_cutoff=0,
                                 upper_cut=5,
                                 log2=False):
    """

    THIS IS DEPRECATED--S. Markson 4 June 2020

    Parameters
    ----------
    genes :
        param metadata:
    expression_data :
        param list_of_gene_windows:
    patient :
        param cell_type:  (Default value = 'tumor')
    complexity_cutoff :
        Default value = 1000)
    cell_type_col_name :
        Default value = 'cell.type')
    patient_col_name :
        Default value = 'patient_ID')
    complexity_col_name :
        Default value = 'nGene')
    metadata :
        
    list_of_gene_windows :
        
    cell_type :
        (Default value = 'tumor')
    patient_columns :
         (Default value = 'Patient_ID')
    cell_type_column :
         (Default value = 'cell.type')

    Returns
    -------

    
    """
    # Nota bene:  patient id gets cast to string below

    genes = loom.ra['gene']
    #  This is very inefficient--make a general function for loom copy-over
    metadata = pd.DataFrame(loom.ca['patient_ID'])
    metadata.columns = ['patient_ID']
    metadata['complexity'] = loom.ca['complexity']
    metadata['cell_type'] = loom.ca['cell_type']
    #    metadata['cell_name'] = loom.ca['cell_names'] # I hate this

    if complexity_cutoff > 0:
        metadata = metadata[metadata[complexity_column] > complexity_cutoff]
    if type(patient) not in [tuple, list]:
        patient = [str(patient)]
    else:
        patient = list(patient)
        patient = [str(x) for x in patient]
    print("debug", patient)
    if cell_type_column == None and cell_type == None:
        relevant_indices = metadata[(
            metadata[patient_column].astype(str).isin(patient))].index.values
    else:
        relevant_indices = metadata[
            (metadata[cell_type_column].astype(str) == str(cell_type))
            &
            (metadata[patient_column].astype(str).isin(patient))].index.values

    if log2:
        relevant_expression_data = 2**loom[:, relevant_indices] - 1
    else:
        relevant_expression_data = loom[:, relevant_indices]

    mean_window_expressions = robust_mean_windowed_expressions(
        genes,
        list_of_gene_windows,
        relevant_expression_data,
        tqdm_desc='Calculating Mean Window Expressions, with "Robustification"',
        upper_cut=upper_cut)
    return mean_window_expressions, metadata.loc[relevant_indices]


def get_ranks(mean_window_expressions):
    """

    Parameters
    ----------
    mean_window_expressions :
        

    Returns
    -------

    
    """
    mean_window_expression_ranks = np.zeros(mean_window_expressions.shape)
    for icell in range(mean_window_expressions.shape[1]):
        mean_window_expression_ranks[:, icell] = stats.rankdata(
            mean_window_expressions[:, icell])
    return mean_window_expression_ranks


def convert_to_sparse(dense_file,
                      sparse_file=None,
                      genes_not_present=False,
                      genelist_file=None,
                      delimiter='\t'):
    """

    Parameters
    ----------
    dense_file :
        
    sparse_file :
         (Default value = None)
    genelist_file :
         (Default value = None)
    delimiter :
         (Default value = '\t')

    Returns
    -------

    """

    N = 20
    iterator = 0
    row = []
    col = []
    data = []
    genes = []
    with open(dense_file, 'r') as infile:
        firstline = islice(infile, 1)
        headings = np.genfromtxt(firstline, dtype=None)
        with tqdm(unit=' rows completed',
                  unit_scale=True,
                  unit_divisor=1024,
                  desc='Converting dense matrix to sparse: ') as pbar:

            while True:
                gen = islice(infile, N)
                chunk = np.genfromtxt(gen, dtype=str, delimiter=delimiter)
                if genes_not_present:
                    expressions = chunk.astype(float)
                else:
                    genes += list(chunk[:, 0])
                    expressions = chunk[:, 1::].astype(float)
                #print(chunk)
                x, y = np.where(expressions > 0)

                for i, j in zip(x, y):
                    row.append(i + iterator)
                    col.append(j)
                    data.append(expressions[i, j])
                if chunk.shape[0] < N:
                    iterator += chunk.shape[0]
                    break
                else:
                    iterator += N
                pbar.update(N)
    expr_mat = coo_matrix((data, (row, col)), shape=(iterator, len(headings)))
    if sparse_file:
        save_npz(sparse_file, expr_mat)
    if genelist_file and not genes_not_present:
        np.savetxt(genelist_file, np.array(genes), delimiter=',', fmt='%s')
    return expr_mat, genes


def get_masked_wme(loom,
                   layername,
                   mask=None,
                   gene_ra='gene',
                   species='homo sapiens',
                   release=106,
                   window_step=50,
                   window_size=50,
                   return_principal_components=None,
                   upper_cut=0,
                   mask_option='load_full'):

    from panopticon.wme import get_list_of_gene_windows, robust_mean_windowed_expressions
    from tqdm import tqdm
    gene_windows = get_list_of_gene_windows(loom.ra[gene_ra],
                                            species=species,
                                            window_step=window_step,
                                            window_size=window_size,
                                            release=release)
    if mask_option == 'scan':
        mwe_parts = []
        with tqdm(total=loom.shape[1]//512, desc='Computing WME in scan mode') as pbar:
            for (ix, selection, view) in loom.scan(items=mask, axis=1):
                mwe_parts.append(
                    robust_mean_windowed_expressions(
                        view.ra[gene_ra],
                        gene_windows,
                        view[layername][:, :],
                        upper_cut=upper_cut,
                        use_tqdm=False
                    ).T)
                pbar.update(1)
            mwe = np.vstack(mwe_parts).T

    else:
        if mask is None:
            X = loom[layername][:, :]
        else:
            if len(mask) != loom.shape[1]:
                raise Exception(
                    "mask must be boolean mask with length equal to the number of columns of loom"
                )

            if mask_option == 'load_full':  # this is to address an h5py performance bog
                X = loom[layername][:, :][:, mask.nonzero()[0]]
            elif mask_option == 'mask_first':
                X = loom[layername][:, mask.nonzero()[0]]
            #if mask_option not in ['load_full','mask_first','scan']:
            else:
                raise Exception(
                    "mask_option must be one of: load_full, mask_first, scan")
        mwe = robust_mean_windowed_expressions(
            loom.ra[gene_ra],
            gene_windows,
            X,
            upper_cut=upper_cut,
        )
    if return_principal_components is not None:
        if type(return_principal_components) != int:
            raise Exception(
                "type of return_principal_components must be None or int")

        from sklearn.decomposition import PCA
        pca = PCA(n_components=return_principal_components)
        return pca.fit_transform(mwe.T)
    else:
        return mwe.T


def generate_wme_pca(loom,
                     layername,
                     gene_ra='gene',
                     species='homo sapiens',
                     release=106,
                     window_step=50,
                     window_size=50,
                     n_principal_components=5,
                     mask_option='scan',
                     overwrite=False):
    from panopticon.wme import get_masked_wme
    pca = get_masked_wme(loom,
                         layername,
                         gene_ra=gene_ra,
                         species=species,
                         release=release,
                         window_step=window_step,
                         window_size=window_size,
                         return_principal_components=n_principal_components,
                         mask_option=mask_option)
    for ipc in range(pca.shape[1]):
        new_ca = '{} WME PC {}'.format(layername, ipc+1)
        if new_ca in loom.ca.keys() and overwrite==False:
            raise Exception("{} in loom.ca.keys(); cannot overwrite unless overwrite is True".format(new_ca))
        loom.ca[new_ca] = pca[:,ipc]

