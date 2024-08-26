import click
import os


def reaptec_main(
    fastq_dir,
    cellranger_output_dir,
    star_reference_dir='./GRCm39_111_stargenome',
    reference_url="http://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.chr.gtf.gz",
    genome_url="http://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
):
    # create whitelist
    command = "zcat {0}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sed -e 's/-1//g' > {0}/outs/filtered_feature_bc_matrix/barcode_whitelist.txt".format(
        cellranger_output_dir)
    os.system(command)

    # create reference
    reference_file = '{}/{}'.format(star_reference_dir, reference_url.split('/')[-1])
    genome_file = '{}/{}'.format(star_reference_dir, genome_url.split('/')[-1])
    if (not os.path.isfile(reference_file)) and (not os.path.isfile(
            reference_file.replace('.gz', ''))):
        command = "wget {} --directory-prefix {}".format(reference_url, star_reference_dir)
        print(command)
        os.system(command)

    if (not os.path.isfile(genome_file)) and (not os.path.isfile(
            genome_file.replace('.gz', ''))):
        command = "wget {} --directory-prefix {}".format(genome_url, star_reference_dir)
        print(command)
        os.system(command)

    if not os.path.isfile(reference_file.replace('.gz', '')):
        command = "gunzip {}".format(reference_file)
        print(command)
        os.system(command)

    if not os.path.isfile(genome_file.replace('.gz', '')):
        command = "gunzip {}".format(genome_file)
        print(command)
        os.system(command)

    reference_file = reference_file.replace('.gz','')
    genome_file = genome_file.replace('.gz','')
    
    if (not os.path.isdir(star_reference_dir)) or (not os.path.isfile('{}/{}'.format(star_reference_dir, 'genomeParameters.txt'))):
        command = "STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./{} --genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang 149 --limitGenomeGenerateRAM 16000000000".format(
            star_reference_dir, genome_file, reference_file)
        print(command)
        os.system(command)

    # perform whitelisted star mapping
    for fastq in [x for x in os.listdir(fastq_dir) if 'R1' in x]:
        outfile_name_prefix = fastq.replace('.fastq.gz', '')
        command = "STAR --runThreadN 32 --genomeDir {0} --readFilesIn {1}/{2} {1}/{3} --soloCBwhitelist {4} --soloBarcodeMate 1 --clip5pNbases 39 0 ".format(
            star_reference_dir, fastq_dir, fastq, fastq.replace('R1', 'R2'),
            cellranger_output_dir +
            '/outs/filtered_feature_bc_matrix/barcode_whitelist.txt')
        command += "--readFilesCommand zcat --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloStrand Reverse --outFileNamePrefix {}".format(
            outfile_name_prefix)
        command += " --outSAMtype BAM SortedByCoordinate --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_Directional_UMItools --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM "
        command += "--limitBAMsortRAM 60000000000"
        print(command)
        os.system(command)

#        # extract only read 1 reads 
#        command = "samtools view -@ 12 -hbf 64 -q 255 {0}_Aligned.sortedByCoord.out.bam > {0}_unique_R1.bam".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#
#        # read 1 5' scRNA-seq reads with an unencoded G extracted from BAM files
#        command = "samtools view -@ 12 -H {0}_unique_R1.bam > {0}_header.sam".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#        command = "samtools view -@ 12 -F 16 {0}_unique_R1.bam".format(outfile_name_prefix)+" | awk -F \'\\t\' 'BEGIN {OFS=\"\\t\"} {BASE = substr($10, 40,1); if ($6 ~ /^40S[0-9]/ && BASE == \"G\") {print $0}}' > "+"{0}_SoftclipG_F.sam".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#
#        command = "samtools view -@ 12 -f 16 {0}_unique_R1.bam".format(outfile_name_prefix)+" | awk -F \'\\t\' 'BEGIN {OFS=\"\\t\"} {ALT = substr($10, length($10)-39,1); if ($6 ~ /[0-9]M40S$/ && ALT == \"C\") {print $0}}' > "+"{0}_SoftclipG_R.sam".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#
#        command = "cat {0}_header.sam {0}_SoftclipG_F.sam {0}_SoftclipG_R.sam ".format(outfile_name_prefix)+"| samtools sort -A 12 -O bam -o SoftclipG_{0}.bam".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#
#        # Duplicated reads remove using umi-tools
#        command = "umi_tools dedup --per-cell -I softclipG_{0}.bam --extract-umi-method=tag --umi-tag=UR --cell-tag=CR -S SoftclipG_{0}_deduplicated.bam".format(outfile_name_prefix)
#        print(command)
#        os.system(command)

#       # reads with cell barcodes corresponding to the list in step 2 are extracted
#       command = "awk \'{print\"CB:Z:\"$1}\' "+"{0}_whitelist.txt > {0}_cell_barcode.txt".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#       command = "samtools view -A 12 -H SoftclipG_{0}_deduplicated.bam > SAM_header".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#        print(command)
#        os.system(command)

#        # a count file was generated for each TSS using the bamToBed function in BEDTools
#        command = "samtools view -A 12 SoftClipG_{}_filtered.bam | ".format(outfile_name_prefix)+"awk \'BEGIN{OFS=\"\\t\"}{print $23}\' > "+"{}_cell_barcode_tmp.txt".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#        command = "bamToBed -i SoftclipG_{0}_filtered.bam | paste - {0}_cell_barcode_tmp.txt |".format(outfile_name_prefix)+" awk \'BEGIN {OFS=\"\\t\'}{if($6==\"+\"){print $1,$2, $2+1,\".\",$7,$6} else {print $1, $4-1,$3,\".\",$7,$6}}\' | sort -k1,1 -k2,2n -k6,6 | bedtools groupby -g 1,2,3,4,5,6 -c 1 -o count | awk \'BEGIN{OFS=\"\\t\"}{if($1 ~/chr/) print $1,$2,$3,$5,$7,$6}\' >"+"{}.CTSS.fwd.rev.cell.barcode.bed".format(outfile_name_prefix)
#        print(command)
#        os.system(command)
#
#        # bidirectionally transcribed candidate enhancers (btcEnhs) were identified.  See https://github.com/anderssonrobin/enhancers/blob/master/scripts/bidir_enhancers.
