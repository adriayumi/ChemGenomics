suppressMessages(library(dada2))

setDadaOpt(OMEGA_A = 1e-40, OMEGA_P = 1e-30, BAND_SIZE = 10, USE_KMERS = TRUE, VECTORIZED_ALIGNMENT = TRUE, GREEDY = FALSE, GAPLESS = FALSE)

getN <- function(x) sum(getUniques(x))

fnFs <- sort(snakemake@input$FASTQ_R1)
fnRs <- sort(snakemake@input$FASTQ_R2)
sample.names <- sapply(strsplit(basename(fnFs), "_trimmed_[1|2].fastq.gz"), `[`, 1) 

filtFs <- file.path("DADA2", paste0(sample.names, "_R1_filtered.fastq.gz"))
filtRs <- file.path("DADA2", paste0(sample.names, "_R2_filtered.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 16, maxLen = 21, maxN = 0, maxEE = 1, minQ = 2, compress = TRUE, multithread = snakemake@threads)

errF <- learnErrors(filtFs, MAX_CONSIST = 25, randomize = FALSE, multithread = snakemake@threads)
errR <- learnErrors(filtRs, MAX_CONSIST = 25, randomize = FALSE, multithread = snakemake@threads)

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, pool = TRUE, BAND_SIZE = 8, multithread = snakemake@threads)
dadaRs <- dada(derepRs, err = errR, pool = TRUE, BAND_SIZE = 8, multithread = snakemake@threads)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, trimOverhang = TRUE)

seqtab <- makeSequenceTable(mergers)

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
track <- cbind("sample" = sample.names, track)
rownames(track) <- c()
write.table(track, file = snakemake@output$TRACK, sep = "\t", row.names = FALSE, quote = FALSE)

count.table <- t(seqtab)
count.table <- cbind("BARCODE" = rownames(count.table), count.table)

rownames(count.table) <- c()
write.table(count.table, file = snakemake@output$ASVTAB, sep = "\t", row.names = FALSE, quote = FALSE)

uniquesToFasta(seqtab, fout = snakemake@output$ASV_FASTA, ids = colnames(seqtab), width = 80)
