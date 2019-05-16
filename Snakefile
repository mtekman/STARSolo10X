from os import path

star_bin="STAR"
# General directories
data_dir="sources"

# You have to find this whitelist files inside the CellRanger distribution
# 
whitelist=path.join(data_dir, "celseq2.barcodes.tsv")

# First read is cDNAFragment sequence
# Second read is the CellBarcodeUMI sequence (Cell+UMI) read
readFilesIn1=path.join(data_dir, "41737_R2.fastq.gz") #R2 is sequences
readFilesIn2=path.join(data_dir, "41737_R1.fastq.gz") #R1 is barcodes

output_dir="outputs"
fasta_dir=path.join(data_dir, "fastas")

# The genome index has to be re-generated with the latest 2.7.0x release
homo_version="GRCh38"
homo_release=96

final_matrix=path.join(output_dir, "matrix.mtx")

rule all:
    input:
        final_matrix


fasta_all=path.join(data_dir, "%s_r%d.all.fa" % (homo_version, homo_release))
annotation=path.join(data_dir, "%s_r%d.gtf" % (homo_version, homo_release))


rule getFasta:
    '''Downloads a Fasta file and puts in the right place'''
    output:
        fasta_all
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-{homo_release}/fasta/homo_sapiens/dna/Homo_sapiens.{homo_version}.dna.toplevel.fa.gz
        gunzip -c Homo_sapiens.{homo_version}.dna.toplevel.fa.gz > {output}
        """
        

rule getAnnotation:
    output:
        gtf=annotation
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-{homo_release}/gtf/homo_sapiens/Homo_sapiens.{homo_version}.{homo_release}.gtf.gz
        gunzip -c Homo_sapiens.{homo_version}.{homo_release}.gtf.gz > {output.gtf}
        """


star_index_directory="star_indices"
star_indices=path.join(star_index_directory, homo_version)
star_index=path.join(star_indices, "SAindex")

star_threads=4
star_overhang=100
        
rule generateStarIndex:
    input:
        fasta=fasta_all,
        gtf=annotation
    output:
        sindex=star_index
    shell:
        "{star_bin} --runThreadN {star_threads} --runMode genomeGenerate --genomeDir {star_indices} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} --sjdbOverhang {star_overhang}"



chromium_protocol="V2"
# 10x Chromium v2 protocol
# Change this for V3
CBstart=1
CBlen=16
UMIstart=17
UMIlen=10
Strand="Forward"
Features="Gene"
UMIdedup="1MM_All"
OutFileNames="%s/ genes.tsv barcodes.tsv matrix.mtx matrixSJ.mtx" % output_dir
        
rule doStarSolo:
    input:
        files1=readFilesIn1,
        files2=readFilesIn2,
        index=star_index,
        gtf=annotation
    output:
        matrix=final_matrix
    shell:
        "{star_bin} --soloType Droplet --soloCBwhitelist {input.gtf} --readFilesIn {input.files1} {input.files2} --soloCBstart {CBstart} --soloCBlen {CBlen} --soloUMIstart {UMIstart} --soloUMIlen {UMIlen} --soloStrand {Strand} --soloFeatures {Features} --soloUMIdedup {UMIdedup} --soloOutFileNames {OutFileNames}"


