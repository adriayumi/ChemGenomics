# conda environment: qc-processing

workdir: '/opt_hd/adrielle/chemical_genomics/02-trimming_adapters'

SAMPLES, = glob_wildcards('/opt_hd/adrielle/chemical_genomics/00-data/{sample}_L006_R1_001.fastq.gz')

rule all:
    input:
        'paired_samples_primer_trimming_report.html'

rule cutadapt_5prime:
    input:
        R1 = '/opt_hd/adrielle/chemical_genomics/00-data/{sample}_L006_R1_001.fastq.gz',
        R2 = '/opt_hd/adrielle/chemical_genomics/00-data/{sample}_L006_R2_001.fastq.gz'
    output:
        R1 = temp('/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{sample}_trimmed_5prime_1.fastq.gz'),
        R2 = temp('/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{sample}_trimmed_5prime_2.fastq.gz')
    log:
        'cutadapt_logs/{sample}_trimmed_primers_5prime.log'
    threads: 6
    conda:
        "../env/qc-processing.yml"
    shell:
        'cutadapt \
        -j {threads} \
        --no-indels \
        --discard-untrimmed \
        -g GATGTCCACGAGGTCTCT \
        -G GTCGACCTGCAGCGTACG \
        -o {output.R1}\
        -p {output.R2} \
        {input.R1} {input.R2} > {log}'


rule cutadapt_3prime:
    input:
        R1 = '/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{sample}_trimmed_5prime_1.fastq.gz',
        R2 = '/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{sample}_trimmed_5prime_2.fastq.gz'
    output:
        R1 = '/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{sample}_trimmed_1.fastq.gz',
        R2 = '/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{sample}_trimmed_2.fastq.gz'
    log:
        'cutadapt_logs/{sample}_trimmed_primers_3prime.log'
    threads: 6
    conda:
        "../env/qc-processing.yml"
    shell:
        'cutadapt \
        -j {threads} \
        --no-indels \
        --discard-untrimmed \
        -a CGTACGCTGCAGGTCGAC \
        -A AGAGACCTCGTGGACATC \
        -o {output.R1} \
        -p {output.R2} \q
        {input.R1} {input.R2} > {log}'

rule fastqc:
    input:q
        expand('/opt_hd/adrielle/chemical_genomics/02-trimming_adapters/{{sample}}_trimmed_{n}.fastq.gz', n=['1','2'])
    output:
        temp(expand('{{sample}}_trimmed_{n}_fastqc.html', n=['1', '2'])),
        temp(expand('{{sample}}_trimmed_{n}_fastqc.zip', n=['1', '2']))
    params:
        output_directory = './'
    conda:
        "../env/qc-processing.yml"
    shell:
        'fastqc -o {params.output_directory} {input}'

rule multiqc:
    input:
        expand('{sample}_trimmed_{n}_fastqc.html', sample=SAMPLES, n=['1', '2']),
        expand('{sample}_trimmed_{n}_fastqc.zip', sample=SAMPLES, n=['1', '2']),
        expand('cutadapt_logs/{sample}_trimmed_primers_5prime.log', sample=SAMPLES),
        expand('cutadapt_logs/{sample}_trimmed_primers_3prime.log', sample=SAMPLES)
    output:
        'paired_samples_primer_trimming_report.html'
    params:
        title = 'paired barcodes samples primer trimming report',
        filename = 'paired_samples_primer_trimming_report'
    conda:
        "../env/qc-processing.yml"
    shell:
        'multiqc --interactive -o ./ -i "{params.title}" -n {params.filename} {input}'
