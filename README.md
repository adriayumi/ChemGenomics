[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4443837-blue.svg)](https://doi.org/10.5281/zenodo.4443837)

# ChemGenomics
Data and Code used for analysis of Chemical Genomics Profiling data associated with the work "Violacein-induced chaperone system collapse underlies multi-stage antiplasmodial activity" by Tavella et al. (under review).


## File Description

Files are listed in the order the analysis should be executed after quality control of reads.

* DADA2_denoising_snakefile.smk = Snakemake file for running denoising script with DADA2 on multiple sequencing read files.
* DADA2_barcodes_denoising.R = R script for denoising of sequencing reads with DADA2.
* Levenshtein_distance_filtering.ipynb = Jupyter Notebook with Python code for filtering denoised sequenced barcodes using Levenshtein distance.
* Levenshtein_distance_filtering.html = Html version of Python Notebook with code for filtering denoised sequenced barcodes using Levenshtein distance.
* asvtab_barcodes.txt = Count values for ASVs (denoised  sequenced barcodes) resulting from DADA2, this file is used in the "Levenshtein_distance_filtering" Python notebook.
* yeast_pool_barcodes_info.tsv = Table containing barcode sequences and which mutated ORF they represent. File used in Levenshtein_distance_filtering Python notebook.
* Violacein_DESeq2_analysis.Rmd = RStudio notebook with R code for differential abundance analysis of barcoded mutant yeast strains treated with violacein or DMSO, using DESeq2.
* Violacein_DESeq2_analysis.nb.html = Html version of RStudio notebook with R code for differential abundance analysis of barcoded mutant yeast strains treated with violacein or DMSO, using DESeq2


## Raw data

Raw data (sequencing reads) used in this work is available at SRA with Accession number: PRJNA689872
