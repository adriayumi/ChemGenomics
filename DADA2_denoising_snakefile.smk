# Snakemake file for running DADA2 Denoising
# DADA2 v.1.9.1

# Set work directory
workdir: 'DADA2'

# Set wildcards for sample names, from adapter trimmed read fastq files
SAMPLES, = glob_wildcards('/chemical_genomics/02-trimming_adapters/{sample}_trimmed_1.fastq.gz')

# Set rule all
rule all:
    input:
        TRACK = 'track_barcodes.txt',
        ASVTAB = 'asvtab_barcodes.txt',
        ASV_FASTA = 'asv_barcodes.fasta'

# Set DADA2 rule with adapter trimmed read fastq files as input
rule dada2:
    input:
        FASTQ_R1 = expand('/chemical_genomics/02-trimming_adapters/{sample}_trimmed_1.fastq.gz', sample=SAMPLES),
        FASTQ_R2 = expand('/chemical_genomics/02-trimming_adapters/{sample}_trimmed_2.fastq.gz', sample=SAMPLES)
    output:
        TRACK = 'track_barcodes.txt',
        ASVTAB = 'asvtab_barcodes.txt',
        ASV_FASTA = 'asv_barcodes.fasta'
    script:
        'DADA2_barcodes_denoising.R'