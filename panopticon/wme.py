"""
wme.py
====================================
wme
"""

## second version
import numpy as np
from tqdm import tqdm
import pandas as pd
from scipy import sparse
from scipy import stats
import matplotlib.pyplot as plt
from itertools import islice
from scipy.sparse import coo_matrix, save_npz
from panopticon.utilities import get_valid_gene_info




def get_list_of_gene_windows(genes, window_size=200, window_step=1):
    """

    Parameters
    ----------
    genes :
        param window_size:  (Default value = 200)
    window_step :
        Default value = 1)
    window_size :
        (Default value = 200)

    Returns
    -------

    
    """
    gene_names, gene_contigs, gene_starts, gene_ends = get_valid_gene_info(genes)

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


def robustify(genes,
              list_of_gene_windows,
              expression_data,
              upper_cut=5,
              windsor=False,
              tqdm_desc=''):
    """

    Parameters
    ----------
    genes :
        param list_of_gene_windows:
    expression_data :
        param upper_cut:  (Default value = 5)
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
    mean_window_expressions = np.zeros((len(list_of_gene_windows),
                                        expression_data.shape[1]))
    with tqdm(total=len(list_of_gene_windows), desc=tqdm_desc) as pbar:
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
    metadata = pd.DataFrame(loom.ca['patient_ID'])
    metadata.columns = ['patient_ID']
    metadata['complexity'] = loom.ca['complexity']
    metadata['cell_type'] = loom.ca['cell_type']

    
    if complexity_cutoff > 0:
        metadata = metadata[metadata[complexity_column]>complexity_cutoff]
    if type(patient) not in [tuple, list]:
        patient = [str(patient)]
    else:
        patient = list(patient)
        patient = [str(x) for x in patient]
    print("debug", patient)
    if cell_type_column==None and cell_type == None:
        relevant_indices = metadata[(metadata[patient_column].astype(str).isin(patient)) ].index.values
    else:
        relevant_indices = metadata[(metadata[cell_type_column].astype(str) == str(cell_type))
                                    & (metadata[patient_column].astype(str).isin(patient))].index.values

    from IPython.core.debugger import set_trace; set_trace()
    if log2:
        relevant_expression_data = 2**loom[:, relevant_indices] - 1
    else:
        relevant_expression_data = loom[:, relevant_indices]

    mean_window_expressions = robustify(
        genes,
        list_of_gene_windows,
        relevant_expression_data,
        tqdm_desc='Calculating Mean Window Expressions, with "Robustification"',
        upper_cut=upper_cut
    )
    return mean_window_expressions


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


def convert_to_sparse(dense_file, sparse_file=None, genes_not_present=False, genelist_file=None, delimiter='\t'):
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
        with tqdm(
                unit=' rows completed',
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
        np.savetxt(genelist_file,np.array(genes),delimiter=',',fmt='%s')
    return expr_mat, genes


