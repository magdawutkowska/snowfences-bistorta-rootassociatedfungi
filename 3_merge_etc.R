#script written by M.Wutkowska 20190429

R.version #R version 3.5.2 (2018-12-20)

#clean the environment
rm(list = ls(all = TRUE))
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)


#set dir
getwd()
setwd("/Users/magdalenawutkowska/Dropbox/+fences_bistorta/baseon_mydata/")


#read multiple runs
fencesA <- readRDS("./seqtab_fencesA.rds")
fencesB <- readRDS("./seqtab_fencesB.rds")


#merge all the ASV tables
library(dada2) #packageVersion("dada2") [1] ‘1.11.1
ASVtable <- mergeSequenceTables(fencesA, fencesB)
dim(ASVtable) #[1]   98 4085
rownames(ASVtable)
colnames(ASVtable)


#remove chimeras
ASVtable_nochim <- removeBimeraDenovo(ASVtable, method="consensus", multithread=TRUE, verbose = TRUE) #Identified 1058 bimeras out of 4085 input sequences.
dim(ASVtable_nochim) #98 3027
class(ASVtable_nochim)


##inspect the data
rownames(ASVtable_nochim)
boxplot(rowSums(ASVtable_nochim))
hist(rowSums(ASVtable_nochim))
boxplot(colSums(ASVtable_nochim))
colnames(ASVtable_nochim)

#now let's check seq lengths and remove these that are shorter than 200bp
table(nchar(getSequences(ASVtable_nochim)))
seqlens <- nchar(getSequences(ASVtable_nochim))
ASVtable_rightlen <- ASVtable_nochim[,seqlens >= 200] #9 sequences were removed
dim(ASVtable_rightlen) #98 3018


#let's sum sequences from forward(A) and reverse (B) orientations into samples, to do that the first step would be to:
#(1)convert row names into a first column

df_asvtable <- as.data.frame(ASVtable_rightlen)
dim(df_asvtable)


library(dplyr)
library(tibble)
df_asvtable_tbl <- as_tibble(rownames_to_column(df_asvtable))
dim(df_asvtable_tbl) #98 3019
remove_rownames(df_asvtable_tbl)
rownames(df_asvtable_tbl)
colnames(df_asvtable_tbl)

class(df_asvtable_tbl)
class(df_asvtable_tbl$rowname)
class(df_asvtable_tbl$ACGCAAGTTGCGCCCGAAGCCTTCGGGCCGAGGGCACGTCTGCATGGGCGTCACGCACAGCGTCGCCCCCACCCCACTCGTGGGGCGTGGGGCGGATTCTGGCCCCCCGTGTGCTCCCGCGCGCGGTCGGCCTAAAATCAGACCCCGTGGCCGCGAAATGCCGCGACGATTGGTGGTGTACGTGGCGGCCTCGAGCCTCCGAACATCGCGTCGCGCCTTCCGTGGCCCTCTGGAGTCAAGAGGACCCTCGAGAGCCCTCCGCCGGTGCGGAGGGGCCTCTCAACCGTTGCGACCCCATGTCAGGCGGGACTACCCGCTGAGTTTAA)


#(2)remove whatever is after "_" in sample column
df_asvtable_tbl$rowname <- sapply(strsplit(basename(df_asvtable_tbl$rowname), "_"), `[`, 1)
as.factor(df_asvtable_tbl$rowname)


#(3)sum values for samples, first convert wide format data frame into long and aggregate based on sample name and variable
library(reshape2)
df_asvtable_tbl_long <- melt(df_asvtable_tbl, id.vars = c("rowname"))
df_asvtable_long_sum <- aggregate(. ~ rowname + variable, df_asvtable_tbl_long, sum)
head(df_asvtable_long_sum)

#(4)go back to wide format
#library(reshape2)
df_asvtable_wide <- dcast(df_asvtable_long_sum, rowname ~ variable)
dim(df_asvtable_wide) #49 3019
class(df_asvtable_wide)


#putting back samples names to the rownames
rownames(df_asvtable_wide) <- df_asvtable_wide$rowname
dim(df_asvtable_wide) #49 3019
df_asvtable_wide <- df_asvtable_wide[,-1]
df_asvtable_wide <- df_asvtable_wide[, -colSums(df_asvtable_wide)<=1]
dim(df_asvtable_wide) #49 3018


#identifying fungi in blank sample contains #55ASV)
df_asvtable_wide_t <- t(df_asvtable_wide)
df_asvtable_wide_t <- as.data.frame(df_asvtable_wide_t)
blank_table <- df_asvtable_wide_t[df_asvtable_wide_t$Blank > 0,]
df_asvtable_wide_t <- subset(df_asvtable_wide_t, select=-c(Blank))

saveRDS(blank_table, "./blank_table_fencesroots.rds") #negative control table / blank table
write.table(blank_table, file="blank_table_fencesroots.txt", row.names=TRUE, col.names=TRUE, sep = '\t')


#subsampling to the same number of reads

library(vegan) #This is vegan 2.5-4

sort(colSums(df_asvtable_wide_t)) #2 samples with <17k, next one 40k and up to 240k, therefore rarefy to 42k and 2 samples lost
asvtable_48k <- df_asvtable_wide_t[, !colSums(df_asvtable_wide_t)<48000] #2 samples were removed D12e (5542 reads)  D12f (17431 reads)
dim(asvtable_48k) #3018   46
asvtable_48k <- rrarefy(t(asvtable_48k), 48000)
dim(asvtable_48k) 
rowSums(asvtable_48k) #perfecto

#remove empty ASVs
asvtable_48k <- as.data.frame(t(asvtable_48k[,!(colSums(asvtable_48k)<1)]))
dim(asvtable_48k) #2850   46


# Save the output: abundance table
saveRDS(asvtable_48k, "./df_asvtable_abu_fencesroots.rds") #ASV table
write.table(asvtable_48k, file="df_asvtable_abu_fencesroots.txt", row.names=TRUE, col.names=TRUE, sep = '\t')


#removing singletones
asvtable_48k_pa <- asvtable_48k
asvtable_48k_pa[asvtable_48k_pa > 0] <- 1

asvtable_pa_nosingl <- asvtable_48k_pa[rowSums(asvtable_48k_pa) > 1,]
dim(asvtable_pa_nosingl) #1069   46
colnames(asvtable_pa_nosingl)


matrix_nosingl <- as.matrix(t(asvtable_pa_nosingl))
tax_fences_bistorta_nosinglet <- assignTaxonomy(matrix_nosingl, "./sh_general_release_dynamic_02.02.2019.fasta", multithread=TRUE, verbose = TRUE)
summary(tax_fences_bistorta_nosinglet)

# Kingdom                      Phylum                   Class    
# k__Fungi:1069               p__Ascomycota       :669   c__Leotiomycetes  :405  
#                             p__Basidiomycota    :396   c__Agaricomycetes :389  
#                             p__Mortierellomycota:  1   c__Dothideomycetes: 76  
#                             NA's                :  3   c__GS37           : 21  
#                                                        c__Lecanoromycetes: 21  
#                                                        (Other)           : 67  
#                                                        NA's              : 90  
# Order                                  Family                Genus    
# o__Helotiales   :375   f__Helotiales_fam_Incertae_sedis:103   g__Tomentella  : 84  
# o__Agaricales   :192   f__Thelephoraceae               : 93   g__Cortinarius : 56  
# o__Thelephorales: 93   f__Helotiaceae                  : 79   g__Inocybe     : 44  
# o__Pleosporales : 54   f__Cortinariaceae               : 56   g__Meliniomyces: 38  
# o__Sebacinales  : 49   f__Tricholomataceae             : 53   g__Cadophora   : 33  
# (Other)         :165   (Other)                         :430   (Other)        :477  
# NA's            :141   NA's                            :255   NA's           :337  
# Species   
# s__bicolor       : 27  
# s__finlandica    : 23  
# s__ciliifera     : 18  
# s__diasemospermus: 18  
# s__geophilum     : 14  
# (Other)          :260  
# NA's             :709  


# Save the output
saveRDS(asvtable_pa_nosingl, "./df_asvtable_pa_fencesroots.rds") #ASV table
write.table(asvtable_pa_nosingl, file="df_asvtable_pa_fencesroots.txt", row.names=TRUE, col.names=TRUE, sep = '\t')

saveRDS(tax_fences_bistorta_nosinglet, "./tax_fences_pa_bistorta.rds") #taxonomy table
write.table(tax_fences_bistorta_nosinglet, file="tax_fences_pa_bistorta.txt", row.names=TRUE, col.names=TRUE, sep = '\t')