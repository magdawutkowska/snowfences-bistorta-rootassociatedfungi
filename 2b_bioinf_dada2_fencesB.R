#A script by Magda Wutkowska, 20190228
# •	demultiplexed using sabre (barcode+primer)
# •	primers and Ns were removed from them in cutadapt prior to dada2 (here) (grep a lot all the primers in all orientations to check what is going on with the data)
# •	Pairs were matched using fastqCombinePairedEnd.py

rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dada2) #packageVersion("dada2") [1] ‘1.11.1

getwd()
setwd("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test")

# File parsing, IMPORTANT! forward files were treated as reverse, and reverse and forward
pathR <- "/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test/fencesB_FWD"
pathF <- "/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test/fencesB_REV"

filtpathF <- file.path(pathF, "filter")
filtpathF
filtpathR <- file.path(pathR, "filter")
filtpathR

fastqFs <- sort(list.files(pathF, pattern="fq"))
fastqRs <- sort(list.files(pathR, pattern="fq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs), 
              truncQ=0, verbose=TRUE, multithread=TRUE, matchIDs = TRUE)

# File parsing
filtFs <- list.files(filtpathF, pattern="fq", full.names = TRUE) #full.names = TRUE adds a whole path to the file name
filtRs <- list.files(filtpathR, pattern="fq", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_R"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_R"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Learn forward error rates
set.seed(100)
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) #123696999 total bases in 468855 reads from 4 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE) #124025357 total bases in 468855 reads from 4 samples will be used for learning the error rates.

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample inference and merger of paired-end reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.namesR

dada_forward <- dada(derep_forward, err=errF, multithread=TRUE, pool="pseudo", verbose = TRUE)
dada_reverse <- dada(derep_reverse, err=errR, multithread=TRUE, pool="pseudo", verbose = TRUE)

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, verbose=TRUE) #minOverlap = 20, maxMismatch = 10

# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

# Construct sequence table
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab)
dim(seqtab) 
View(t(seqtab))
summary.matrix(t(seqtab))

saveRDS(seqtab, "./seqtab_fencesB.rds")
#to read it back to R: fencesB <- readRDS("./seqtab_fencesB.rds")


# step-testing:
# #chimera removal
# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# dim(seqtab.nochim)
# sum(seqtab.nochim)/sum(seqtab) #frequency of chimera seqs 0.9791554 ;)
# View(t(seqtab.nochim))
# 
# table(nchar(getSequences(seqtab.nochim)))
# 
# 
# #assign taxonomy
# unite.ref <- "./sh_general_release_dynamic_02.02.2019.fasta"
# taxa <- assignTaxonomy(seqtab.nochim, unite.ref, taxLevels = c("Kingdom", "Phylum", "Class",
#                                                                "Order", "Family", "Genus", "Species"), multithread = TRUE, tryRC = TRUE, verbose=TRUE)
# summary(taxa)