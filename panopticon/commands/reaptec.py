import click
import os

def reaptec_main(fastq_dir, cellranger_output_dir, star_reference_dir):
    # create whitelist
    command = "{0}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed -e 's/-1//g' > {0}/outs/filtered_feature_bc_matrix/barcode_whitelist.txt".format(cellranger_output_dir)
    os.system(command)

    # perform whitelisted star mapping
    for fastq in [x for x in os.listdir(fastq_dir) if 'R1' in fastq]:
        command = "STAR --runThreadN 32 --genomeDIR {} --readFilesIn {} {} --solocCBwhitelist {} --soloBarcodeMate 1 --clip5bNbases 39 0 ".format(star_reference_dir, fastq, fastq.replace('R1','R2',cellranger_output_dir+'/outs/filtered_feature_bc_matrix/barcode_whitelist.txt')
        command += "--readFilesCommand zcat --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloStrand Reverse --outFileNamePrefix {}".format(fastq.replace('.fastq.gz','')
        command += "--outSAMtype BAM SortedByCoordinate --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_Directional_UMItools --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM "
        command += "--limitBAMsortRAM 60000000000"
        print(command)






       
