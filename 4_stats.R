#this script was writen to analyse Bistorta vivipara root-associated fungal data from snow fence study

R.version.string #[1] "R version 3.5.2 (2018-12-20)"
setwd("/Users/magdalenawutkowska/Dropbox/+fences_bistorta/baseon_mydata/")

library(vegan) #his is vegan 2.5-4
library(tidyverse)
# ── Attaching packages ─────────────────── tidyverse 1.2.1 ──
# ✓ ggplot2 3.2.1     ✓ purrr   0.3.4
# ✓ tibble  3.0.1     ✓ dplyr   0.8.5
# ✓ tidyr   0.8.3     ✓ stringr 1.4.0
# ✓ readr   1.3.1     ✓ forcats 0.5.0
# ── Conflicts ────────────────────── tidyverse_conflicts() ──
# x dplyr::filter() masks stats::filter()
# x dplyr::lag()    masks stats::lag()



########################################################################
####Read OTU table, taxonomy assigment table and edaphic parameters ####
########################################################################

asv <- readRDS("df_asvtable_pa_fencesroots.rds")
dim(asv) #1069   46
summary(asv)

tax <- read_delim(file = "tax_fences_pa_bistorta.txt", delim = "\t", col_types = 
                          cols(
                            Identity = col_factor(),
                            ECM_annotation = col_factor(),
                            general_annotation = col_factor(),
                            Kingdom = col_factor(),
                            Phylum = col_factor(),
                            Class = col_factor(),
                            Order = col_factor(),
                            Family = col_factor(),
                            Genus = col_factor(),
                            Species = col_factor()
                          ))
dim(tax) #[1] 1069    7
summary(tax)


env <- read_delim(file = "env.txt", delim = "\t", col_types = 
                    cols(
                      Individual = col_character(),
                      Vegetation = col_factor(),
                      Block = col_factor(),
                      Fence = col_factor(),
                      Treatment = col_factor(),
                      Code = col_factor(),
                      pH = col_number(),
                      Moisture = col_number(),
                      Conductivity = col_number(),
                      OM = col_number(),
                      N = col_number(),
                      C = col_number(),	
                      C_N_ratio	= col_number(),
                      Rh_Length	= col_number(), 
                      Root_weight	= col_number(), 
                      Root_length	= col_number(),
                      Northing = col_number(),
                      Easting = col_number()
                    ))


#number of all detected AVSs, Shannon diversity
env$ASV <- specnumber(t(asv))
asv_t <- t(asv)
asv_t <- as.data.frame(asv_t)

env$shannon <- diversity(asv_t, index = "shannon")

boxplot(env$ASV)
boxplot(env$ASV ~ env$Fence)
boxplot(env$ASV ~ env$Vegetation)
boxplot(env$ASV ~ env$Block)
boxplot(env$ASV ~ env$Treatment)

boxplot(env$shannon)
boxplot(env$shannon ~ env$Fence)
boxplot(env$shannon ~ env$Vegetation)
boxplot(env$shannon ~ env$Block)
boxplot(env$shannon ~ env$Treatment)



#######################################################################################################################################
####GNMDS was based on script https://datadryad.org/bitstream/handle/10255/dryad.74203/Script%20for%20NMDS%20analysis.R?sequence=1 ####
#######################################################################################################################################
####FOR ALL OTUS pa

library(vegan)
library(MASS)
library(stats)

#otu table
attach(asv)
names(asv)
str(asv)

#environmental data
attach(env)
names(env)
str(env)

#making Bray-Curtis dissimilarity matrix:
dist.all<-vegdist(t(asv), method="jaccard", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
dist.all 

mds.all<-NULL
for(i in 1:100)
{mds.all[[i]]<-isoMDS(dist.all,initMDS(dist.all, k=2), k=2, maxit=1000,tol=1e-7)}
mds.stress.all<-unlist(lapply(mds.all,function(v){v[[2]]})) 

dist.all.clust<-hclust(dist.all,"single")      # hierarchical clustering 
plot(dist.all.clust,ylab="Dissimilarity",xlab="",sub="")  # plot dendrogram

# looking at the stress values for 100 mds:
mds.stress.all
#ordering the stress values for the 100 mds:
order(mds.stress.all)
#Saving the order in a vector
ordered.all<-order(mds.stress.all)
ordered.all

#find the stress of the solutions with the lowest and second lowest stress:
mds.stress.all[ordered.all[1]]
mds.stress.all[ordered.all[2]]

#scaling of axes to half change units and varimax rotation
mds.best.all<-postMDS(mds.all[[ordered.all[1]]],dist.all)
mds.best.all
mds.secbest.all<-postMDS(mds.all[[ordered.all[2]]],dist.all)
mds.secbest.all

#Procrustes comparisons
procrustes(mds.best.all,mds.secbest.all,permutations=999)
protest(mds.best.all,mds.secbest.all,permutations=999)
plot(procrustes(mds.best.all,mds.secbest.all,permutations=999))

#making variables from GNMDS axes 1 and 2 for plotting
gnmds1_all<-mds.best.all$points[,1]
gnmds2_all<-mds.best.all$points[,2]

##Biplot ##both using origin=TRUE here
mds.df.all<-data.frame(gnmds1_all,gnmds2_all)
fit.all<-envfit(mds.df.all, env[,c(2:3,5,7:16)], 999)
fit.all 

par(mfrow=c(1,1))
plot(main="GNMDS_all", gnmds1_all,gnmds2_all,xlab="GNMDS1",ylab="GNMDS2",type="n", xlim = c(-0.5, 0.5), ylim= c(-0.5, 0.5))
lines(c(-1,1),c(0,0),lty=2,col=8)
lines(c(0,0),c(-1,1),lty=2,col=8)

points(gnmds1_all[Code=="MD"], gnmds2_all[Code=="MD"], pch=19, col='grey')
points(gnmds1_all[Code=="MC"], gnmds2_all[Code=="MC"], pch=21, col='grey')

points(gnmds1_all[Code=="HD"], gnmds2_all[Code=="HD"], pch=15, col='grey')
points(gnmds1_all[Code=="HC"], gnmds2_all[Code=="HC"], pch=22, col='grey')

ordiellipse(mds.df.all, env$Treatment, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="black", border = "black")
ordiellipse(mds.df.all, env$Vegetation, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "black")
ordiellipse(mds.df.all, env$Block, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "black")
ordiellipse(mds.df.all, env$Code, display = "sites", kind = c("se"), conf = 0.95, draw = c("lines"), alpha = 2, label = TRUE, lwd=1, col="darkgreen", border = "black")

legend("bottomleft", legend = c("MD", "MC", "HD", "HC") , pch=c(19,21,15,22), col=c("grey", "grey", "grey","grey"), bty="n")
plot(fit.all)

detach(asv)
detach(fungi_pa)



###################
#### PERMANOVA ####
###################

adonis(asv_t ~ Vegetation, data = env, method = "jaccard", permutations = 999)
adonis(asv_t ~ Treatment, data = env, method = "jaccard", permutations = 999)
adonis(asv_t ~ Block, data = env, method = "jaccard", permutations = 999)

adonis(asv_t ~ Vegetation * Treatment, data = env, method = "jaccard", permutations = 999)
adonis(asv_t ~ Vegetation * Treatment, data = env, method = "jaccard", permutations = 999, strata = env$Block)
adonis(asv_t ~ Vegetation * Treatment, data = env, method = "jaccard", permutations = 999, strata = env$Fence)

adonis(asv_t ~ Vegetation + Treatment, data = env, method = "jaccard", permutations = 999)
adonis(asv_t ~ Vegetation + Treatment, data = env, method = "jaccard", permutations = 999, strata = env$Block)
adonis(asv_t ~ Vegetation + Treatment, data = env, method = "jaccard", permutations = 999, strata = env$Fence)




###################
#### DIVERSITY ####
###################

table(factor(env$Code))

env %>%
  ggplot(aes(x=Treatment, y=shannon)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Vegetation) +
  theme_bw() +
  geom_jitter()

env %>%
  ggplot(aes(x=Treatment, y=ASV)) + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Vegetation) +
  theme_bw() +
  geom_jitter()

cor.test(env$shannon, env$ASV, method = "pearson") #0.9842 CI 0.9714587 0.9912785
