[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4443837-blue.svg)](https://doi.org/10.5281/zenodo.4443837)

# ChemGenomics
Data, Code and Methods used for analysis of Chemical Genomics Profiling data associated with the work "Violacein-induced chaperone system collapse underlies multi-stage antiplasmodial activity" by Tavella et al., accepted for publication by ACS Infectious Diseases. 

Tavella et al. (2021). "Violacein-induced chaperone system collapse underlies multi-stage antiplasmodial activity". ACS Infectious Diseases.
https://dx.doi.org/10.1021/acsinfecdis.0c0045

## Raw data

Raw data (sequencing reads) used in this work is available at SRA with Accession number: PRJNA689872

## File Description

Files are listed in the order the analysis should be executed after quality control of reads.

* DADA2_denoising_snakefile.smk = Snakemake file for running denoising script with DADA2 on multiple sequencing read files.
* DADA2_barcodes_denoising.R = R script for denoising of sequencing reads with DADA2.
* Levenshtein_distance_filtering.ipynb = Jupyter Notebook (.ipynb) with Python code for filtering denoised sequenced barcodes using Levenshtein distance. Notebook also available in html version (.html) to view in web browser. 
* asvtab_barcodes.txt = Count values for ASVs (denoised  sequenced barcodes) resulting from DADA2, this file is used in the "Levenshtein_distance_filtering" Python notebook.
* yeast_pool_barcodes_info.tsv = Table containing barcode sequences and which mutated ORF they represent. File used in Levenshtein_distance_filtering Python notebook.
* Violacein_DESeq2_analysis = RStudio notebook (.Rmd) with R code for differential abundance analysis of barcoded mutant yeast strains treated with violacein or DMSO, using DESeq2.  Notebook also available in html version (.html) to view in web browser. 

## Methods

### Raw data pre-processing
The initial quality of the raw reads was evaluated using the FastQC (version 0.11.7) [1] and MultiQC (version 1.6) [2] tools. Primers and adapters (U1 foward: GATGTCCACGAGGTCTCT, U2 reverse: GTCGACCTGCAGCGTACG) were removed using the Cutadapt version 1.16 tool [3]. In this step, primer sequences with insertions and deletions were not allowed and unprocessed read pairs were discarded (options --no-indels, --discard-untrimmed). After removing primers and adapters, most reads should be around 20 base pairs in size, which is the expected size of the barcode sequence. 

### Denoising of Amplicon sequence variants
The identification of amplicon sequence variants (ASV) was performed using the DADA2 denoising algorithm (version 1.9.1) [4]. First, reads R1 and R2 were filtered to discard those that contained more than one expected error (maxEE = 1), quality score less than 2 (minQ = 2) and size less than 16 bp or greater than 21 bp (minLen = 16, maxLen = 21). Then, the parameters of the error models were obtained by alternating the inference of samples with the estimate of parameters until the convergence was reached. The error models and the duplicate readings grouped from all the samples were used as input for the given function to perform the denoising of reads R1 and R2 (options OMEGA_A = 1e-40, BAND_SIZE = 10, USE_KMERS = TRUE, VECTORIZED_ALIGNMENT = TRUE, GREEDY = FALSE, GAPLESS = FALSE). After this step, the pairs of reads R1 and R2 with a minimum overlap of 10 bp and no mismatch were merged to obtain the ASVs. 
Since it is known that several barcodes can have different sequences from the sequences originally described, we allowed ASVs with a Levenshtein distance of up to 2 to barcodes with greater similarity. Thus, 5447 ASVs, which correspond to 5,405 barcodes, were taken for subsequent analyzes. Finally, the matrix containing the count of the number of reads per barcode has been updated to replace each barcode sequence with the respective mutated ORF name it represents. 

### Differential Abundance Analysis 
The DESeq2 package (version 1.20.0) [5] was used to normalize the counts and estimate the differential abundance of barcodes between samples treated against untreated controls. Initially, a filter was made to remove barcodes with low counts throughout the samples (at least 10 observations in at least 6 samples). Normalization between samples was performed using the library size factor method [6], using only ASVs with a number of observations greater than zero.
Usually in chemical genomics using S. cerevisiae, heterozygous yeast pools are treated with EC20 concentration of a compound of interest for 20 generations. In this work, we opted to treat the yeast pools with violacein EC20 concentration for 5 and 10 generations only, in a way to avoid false positives hits. This decision came with the onus that statistics could be compromised, once the growth gap between yeasts with growth defects and nonaffected mutants were not as dramatic as they would be if they were set for 20 generations. For this reason, to find whether the drug treatment induces a change in yeast abundance at any time point, we used an analysis design that included the generation as a regression variable and then tested using the likelihood ratio test where the generation factor was removed in the reduced formula (test = "LRT"). We considered significant hits strains that presented an adjusted p-value < 0.01 and log2 fold change < 0 (growth defect).

#### References
[1] ANDREWS, Simon et al. FastQC: a quality control tool for high throughput sequence data. 2010.

[2] EWELS, Philip et al. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, v. 32, n. 19, p. 3047-3048, 2016.

[3] MARTIN, Marcel. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, v. 17, n. 1, p. pp. 10-12, 2011.

[4] CALLAHAN, Benjamin J. et al. DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, v. 13, n. 7, p. 581, 2016.

[5] Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome biology 15.12 (2014): 550.

[6] Anders, Simon, and Wolfgang Huber. "Differential expression analysis for sequence count data." Genome biology 11.10 (2010): R106.
