import click
import os


def reaptec_main(fastq_dir, cellranger_output_dir, star_reference_dir,
                 genome_url="http://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.chr.gtf.gz",
                 reference_url="http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"):
    # create whitelist
    command = "zcat {0}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed -e 's/-1//g' > {0}/outs/filtered_feature_bc_matrix/barcode_whitelist.txt".format(
        cellranger_output_dir)
    os.system(command)

    # create reference
    command = "wget {}".format(reference_url)
    print(command)
    os.system(command)

    command = "wget {}".format(genome_url)
    print(command)
    os.system(command)

    reference_file = reference_url.split('/')[-1]
    command = "gunzip {}".format(reference_file)
    print(command)
    os.system(command)

    genome_file = genome_url.split('/')[-1]
    command = "gunzip {}".format(genome_file)
    print(command)
    os.system(command)

    command = "STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./GRCm39_index --genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang 149".format(genome_file, reference_file)
    print(command)
    os.system(command)

    # perform whitelisted star mapping
    for fastq in [x for x in os.listdir(fastq_dir) if 'R1' in x]:
        command = "STAR --runThreadN 32 --genomeDir {0} --readFilesIn {1}/{2} {1}/{3} --soloCBwhitelist {4} --soloBarcodeMate 1 --clip5pNbases 39 0 ".format(
            star_reference_dir, fastq_dir, fastq, fastq.replace('R1', 'R2'),
            cellranger_output_dir +
            '/outs/filtered_feature_bc_matrix/barcode_whitelist.txt')
        command += "--readFilesCommand zcat --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloStrand Reverse --outFileNamePrefix {}".format(
            fastq.replace('.fastq.gz', ''))
        command += "--outSAMtype BAM SortedByCoordinate --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_Directional_UMItools --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM "
        command += "--limitBAMsortRAM 60000000000"
        print(command)
