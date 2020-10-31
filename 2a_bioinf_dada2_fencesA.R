#A script by Magda Wutkowska, 20190227
# •	demultiplexed using sabre (barcode+primer)
# •	primers and Ns were removed from them in cutadapt prior to dada2 (here) (grep a lot all the primers in all orientations to check what is going on with the data)
# •	Pairs were matched using fastqCombinePairedEnd.py

rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(dada2) #packageVersion("dada2") [1] ‘1.11.1

getwd()
setwd("/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test")
# File parsing
pathF <- "/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test/fencesA_FWD"
pathR <- "/Users/magdalenawutkowska/Dropbox/PhD/PROJECT_3/r-dada2_test/fencesA_REV"

filtpathF <- file.path(pathF, "filter")
filtpathF
filtpathR <- file.path(pathR, "filter")
filtpathR

fastqFs <- sort(list.files(pathF, pattern="fq"))
fastqRs <- sort(list.files(pathR, pattern="fq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs), 
              truncQ=0, verbose=TRUE, multithread=TRUE, matchIDs = TRUE, minLen = 50)


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
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE) #119736528 total bases in 452803 reads from 4 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE) #120394697 total bases in 452803 reads from 4 samples will be used for learning the error rates.


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

saveRDS(seqtab, "./seqtab_fencesA.rds")
#to read it back to R: fencesA <- readRDS("./seqtab_fencesA.rds")


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
# # Kingdom                      Phylum                   Class                  Order    
# # k__Fungi:1603   p__Ascomycota       :928   c__Agaricomycetes :637   o__Helotiales   :474  
# # p__Basidiomycota    :663   c__Leotiomycetes  :506   o__Agaricales   :298  
# # p__Mortierellomycota:  6   c__Dothideomycetes: 90   o__Thelephorales:155  
# # p__Mucoromycota     :  1   c__Eurotiomycetes : 40   o__Sebacinales  : 79  
# # NA's                :  5   c__Sordariomycetes: 40   o__Pleosporales : 61  
# # (Other)           :123   (Other)         :292  
# # NA's              :167   NA's            :244  
# # Family                 Genus               Species    
# # f__Thelephoraceae               :155   g__Tomentella   :123   s__bicolor    :  34  
# # f__Helotiales_fam_Incertae_sedis:116   g__Inocybe      : 93   s__finlandica :  24  
# # f__Helotiaceae                  :108   g__Cortinarius  : 68   s__ericae     :  20  
# # f__Inocybaceae                  : 94   g__Pezoloma     : 43   s__trachyspora:  19  
# # f__Cortinariaceae               : 68   g__Phialocephala: 43   s__badia      :  16  
# # (Other)                         :646   (Other)         :689   (Other)       : 434  
# # NA's                            :416   NA's            :544   NA's          :1056 