import click
import os
import loompy
import numpy as np
import pandas as pd
from scipy import sparse
from panopticon.dna import segmentation_to_copy_ratio_dict
from panopticon.utilities import get_valid_gene_info

def scrna_wizard_main():
    
    filepath = click.prompt("Location of .loom file", type=click.Path('wb'))
    filename = click.prompt("Name of .loom file")
    if not filename.endswith('.loom'):
        filename+='.loom'
    matrixpath = click.prompt("Data/counts matrix (sparse npz or dense txt)", type=click.File('rb'))
    if str(matrixpath).endswith('.npz'):
        matrix = sparse.load_npz(matrixpath)
    elif str(matrixpath).endswith('.csv'):
        hasheader = click.prompt('Does that file have a header?', type=click.Choice(['n','y']), default='n')
        if hasheader == 'n': 
            matrix = pd.read_table(matrixpath, header=None, sep=',')
        elif hasheader == 'y':
            matrix = pd.read_table(matrixpath, sep=',')
        matrix = matrix.iloc[:,(matrix.dtypes != object).values]
        matrix = matrix.values
    else:
        hasheader = click.prompt('Does that file have a header?', type=click.Choice(['n','y']), default='n')
        if hasheader == 'n': 
            matrix = pd.read_table(matrixpath, header=None)
        elif hasheader == 'y':
            matrix = pd.read_table(matrixpath)
        matrix = matrix.iloc[:,(matrix.dtypes != object).values]
        matrix = matrix.values

    metadatapath = click.prompt("Cell metadata (pandas-loadable)", type=click.File('rb'))
    metadata = pd.read_table(metadatapath)
    iscomplexity = click.prompt("Does this cell metadata file have a column corresponding to complexity?", type=click.Choice(['n','y']), default='n')
    if iscomplexity == 'y':
        complexity_col = click.prompt("Which of these columns corresponds to the cell complexity?", type=click.Choice(metadata.columns))
        if complexity_col != 'complexity':
            metadata['complexity'] = metadata[complexity_col]
            metadata.drop(complexity_col, inplace=True, axis=1)
    ispatientid = click.prompt("Does this cell metadata file have a column corresponding to patient identity?", type=click.Choice(['n','y']), default='n')
    if ispatientid == 'y':
        patient_col = click.prompt("Which of these columns corresponds to the patient identity?", type=click.Choice(metadata.columns))
        if patient_col != 'patient_ID':
            metadata['patient_ID'] = metadata[patient_col]
            metadata.drop(patient_col, inplace=True, axis=1)
    iscelltype = click.prompt("Does this cell metadata file have a column corresponding to cell type?", type=click.Choice(['n','y']), default='n')
    if iscelltype == 'y':
        cell_type_col = click.prompt("Which of these columns corresponds to the cell type?", type=click.Choice(metadata.columns))
        if cell_type_col != 'cell_type':
            metadata['cell_type'] = metadata[cell_type_col]
            metadata.drop(cell_type_col, inplace=True, axis=1)
    genepath = click.prompt("Gene metadata (or simply genelist)", type=click.File('rb'))
    isheader = click.prompt("Does this file have a header?", type=click.Choice(['n','y']), default='n')
    if isheader=='n':
        genes = pd.read_table(genepath, header=None)
        genes.columns = ['gene']
    else:
        genes = pd.read_table(genepath)
        gene_col = click.prompt("Which of these columns corresponds to the gene name?", type=click.Choice(genes.columns))
        if gene_col != 'gene':
            genes['gene'] = gene[gene_col]
            genes.drop(gene_col, inplace=True, axis=1)
    
    loompy.create(filepath+'/'+filename, matrix,  genes.to_dict("list"),metadata.to_dict("list"))
    print("Loom file creation complete.")

def cnv_wizard_main():
    loomfile = click.prompt("Loom file that you would like to augment with cnv/segmentation data: ")
    while not (os.path.isfile(loomfile) and loomfile.endswith('.loom')):
        loomfile = click.prompt("Not a loom file.  Please select loom file that you would like to augment with cnv/segmentation data: ")
    segmentation = None
    segmentationfile = click.prompt("Segmentation file that you would like to add to loom file: " )
    while not (os.path.isfile(segmentationfile)):
        segmentationfile = click.prompt("Not a valid file.  Please select segmentation file that you would like to the loom file: ")
    segmentation = pd.read_table(segmentationfile)
    with loompy.connect(loomfile, validate=False) as loom:
        chromosome = click.prompt("Which column of this segmentation corresponds to the chromosome?", type=click.Choice(segmentation.columns))
        chromosome_start = click.prompt("Which column of this segmentation corresponds to the chromosome start?", type=click.Choice(segmentation.columns))
        chromosome_end = click.prompt("Which column of this segmentation corresponds to the chromosome end?", type=click.Choice(segmentation.columns))
        chromosome_tcr = click.prompt("Which column of this segmentation corresponds to the copy ratio (or log_2(copy ratio))?", type=click.Choice(segmentation.columns))
        log2 = click.prompt("Was that the log2(copy ratio) (i.e., is the copy ratio 2^(value given))?", type=click.Choice(['n','y']), default='n')
        segmentation['chrom'] = segmentation[chromosome]
        segmentation['chromStart'] = segmentation[chromosome_start]
        segmentation['chromEnd'] = segmentation[chromosome_end]
        segmentation['copyRatio'] = segmentation[chromosome_tcr]
        if log2 == 'y':
            segmentation['copyRatio'] = segmentation['copyRatio'].apply(lambda x: 2**x)

        gene_to_cnv = segmentation_to_copy_ratio_dict(loom.ra['gene'], segmentation)
        segmentation_name = click.prompt("How would you like to label this segmentation?")
        loom.ra[segmentation_name] = [gene_to_cnv[gene] if gene in gene_to_cnv.keys() else np.nan for gene in loom.ra['gene']]
    print("CNV Segmentation addition complete.")

def gene_position_augmentation_main(loomfile):
    with loompy.connect(loomfile, validate=False) as loom:
        gene_names, gene_contigs, gene_starts, gene_ends = get_valid_gene_info(loom.ra['gene'])
        gene_to_contig = {gene:contig for gene, contig in zip(gene_names, gene_contigs)}
        gene_to_start = {gene:start for gene, start in zip(gene_names, gene_starts)}
        gene_to_end = {gene:end for gene, end in zip(gene_names, gene_ends)}
        loom.ra.chromosome = [gene_to_contig[gene] if gene in gene_to_contig.keys() else np.nan for gene in loom.ra['gene']]
        loom.ra.start = [gene_to_start[gene] if gene in gene_to_start.keys() else np.nan for gene in loom.ra['gene']]
        loom.ra.end = [gene_to_end[gene] if gene in gene_to_end.keys() else np.nan for gene in loom.ra['gene']]


