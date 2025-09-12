### Processing amplicon data using dada2 ###
# I just pre-process/trim the reads here, then after they are trimmed I rename them to original naming convention and import into QIIME2 to assign taxonomy 

library(dada2)

setwd("/Users/syrenawhitner/Desktop/sediment/5077_Illumina")

# Set file path
path = "/Users/syrenawhitner/Desktop/sediment/5077_Illumina"
fnFs = sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# select first part for sample name 
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# select second part for sample name 
sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 2)

# Visualize quality profiles for the first few samples
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Define file paths for filtered files
filt_path = file.path(path, "filtered")
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# filter samples! --> pick values based on your region, these are for SSU V9 region 
# more stringent filtering, based on EE values instead of length 
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=F, verbose = TRUE)
head(out)

# checking error rates
errF = learnErrors(filtFs, multithread=TRUE)
errR = learnErrors(filtRs, multithread=TRUE)

# plot error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplication
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)
names(derepFs) = sample.names
names(derepRs) = sample.names

# Denoise sequences
dadaFs = dada(derepFs, err=errF, multithread=TRUE)
dadaRs = dada(derepRs, err=errR, multithread=TRUE)