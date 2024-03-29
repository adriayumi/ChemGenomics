# conda environment: qc-processing
# run fastqc and multiqc on raw fastq files

workdir: '/opt_hd/adrielle/chemical_genomics/QC_reports'

SAMPLES, = glob_wildcards('/opt_hd/adrielle/chemical_genomics/raw_data/{sample}_R1_001.fastq.gz')

rule all:
    input:
        'multiqc/barcodes_samples_fastqc_report.html'

rule fastqc:
    input:
        expand('/opt_hd/adrielle/chemical_genomics/raw_data/{{sample}}_R{n}_001.fastq.gz', n=['1','2'])
    output:
        temp(expand('{{sample}}_R{n}_001_fastqc.html', n=['1', '2'])),
        temp(expand('{{sample}}_R{n}_001_fastqc.zip', n=['1', '2']))
    params:
        output_directory = './'
    conda:
        "../env/qc-processing.yml"
    shell:
        'fastqc -o {params.output_directory} {input}'

rule multiqc:
    input:
        expand('{sample}_R{n}_001_fastqc.html', sample=SAMPLES, n=['1', '2']),
        expand('{sample}_R{n}_001_fastqc.zip', sample=SAMPLES, n=['1', '2'])
    output:
        'multiqc/barcodes_samples_fastqc_report.html'
    params:
        title = 'barcodes samples FastQC report',
        filename = 'barcodes_samples_fastqc_report',
        output_directory = 'multiqc'
    conda:
        "../env/qc-processing.yml"
    shell:
        'multiqc --interactive -o {params.output_directory} -i "{params.title}" -n {params.filename} {input}'

