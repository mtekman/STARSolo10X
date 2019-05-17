from os import path

star_bin="STAR"
# General directories
data_dir="sources"

# You have to find this whitelist files inside the CellRanger distribution
#
whitelist=path.join(data_dir, "737K-august-2016.txt")

# First read is cDNAFragment sequence
# Second read is the CellBarcodeUMI sequence (Cell+UMI) read
readFilesIn1=path.join(data_dir, "41737_R2.fastq.gz") #R2 is sequences
readFilesIn2=path.join(data_dir, "41737_R1.fastq.gz") #R1 is barcodes

output_dir="outputs"
fasta_dir=path.join(data_dir, "fastas")

# The genome index has to be re-generated with the latest 2.7.0x release
chromium_version="V3"
homo_version="GRCh38"
homo_release=96

final_matrix=path.join(output_dir, chromium_version, "matrix.mtx")

fasta_all=path.join(data_dir, "%s_r%d.all.fa" % (homo_version, homo_release))
annotation=path.join(data_dir, "%s_r%d.gtf" % (homo_version, homo_release))

star_threads=15
star_overhang=100
star_index_directory="star_indices"
star_indices=path.join(star_index_directory, homo_version)
star_index=path.join(star_indices, "SAindex")

OutFileNames="genes.tsv barcodes.tsv matrix.mtx matrixSJ.mtx"

rule all:
    input:
        final_matrix


rule getFasta:
    '''Downloads a Fasta file and puts in the right place'''
    output:
        fasta_all
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-{homo_release}/fasta/homo_sapiens/dna/Homo_sapiens.{homo_version}.dna.primary_assembly.fa.gz
        gunzip -c Homo_sapiens.{homo_version}.dna.primary_assembly.fa.gz > {output}
        """


rule getAnnotation:
    output:
        gtf=annotation
    shell:
        """
        wget ftp://ftp.ensembl.org/pub/release-{homo_release}/gtf/homo_sapiens/Homo_sapiens.{homo_version}.{homo_release}.gtf.gz
        gunzip -c Homo_sapiens.{homo_version}.{homo_release}.gtf.gz > {output.gtf}
        """

rule generateStarIndex:
    input:
        fasta=fasta_all,
        gtf=annotation
    output:
        sindex=star_index
    log:
        "logs/generateStarIndex.log"
    shell:
        "{star_bin} --runThreadN {star_threads} --runMode genomeGenerate --genomeDir {star_indices} --genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} --sjdbOverhang {star_overhang}"


rule doStarSoloV3:
    input:
        files1=readFilesIn1,
        files2=readFilesIn2,
        index=star_index,
    params:
        # See: https://assets.ctfassets.net/an68im79xiti/51xGuiJhVKOeIIceW88gsQ/1db2c9b5c9283d183ff4599fb489a720/CG000183_ChromiumSingleCell3__v3_UG_Rev-A.pdf
        # Page 14
        CBstart=1,
        CBlen=16,
        UMIstart=17,
        UMIlen=12,
        Strand="Forward",
        Features="Gene",
        UMIdedup="1MM_All",
        outdir=path.join(output_dir, chromium_version)
    output:
        matrix=final_matrix
    log:
        "logs/doStarSoloV3.log"
    shell:
        "{star_bin} --soloType Droplet --soloCBwhitelist {whitelist} --readFilesIn {input.files1} {input.files2} --readFilesCommand zcat --genomeDir {star_indices} --soloCBstart {params.CBstart} --soloCBlen {params.CBlen} --soloUMIstart {params.UMIstart} --soloUMIlen {params.UMIlen} --soloStrand {params.Strand} --soloFeatures {params.Features} --soloUMIdedup {params.UMIdedup} --soloOutFileNames {params.outdir} {OutFileNames}"


# rule doStarSoloV2:
#         # See: https://assets.ctfassets.net/an68im79xiti/UhAMGmlaEMmYMaA4A4Uwa/d65ff7b9bb5e88c2bb9e15e58f280e18/CG00052_SingleCell3_ReagentKitv2UserGuide_RevE.pdf
#         # Page 15
#         CBstart=1,
#         CBlen=16,
#         UMIstart=17,
#         UMIlen=10,
#         Strand="Forward",
#         Features="Gene",
#         UMIdedup="1MM_All",
#         outdir=path.join(output_dir, "V2")
