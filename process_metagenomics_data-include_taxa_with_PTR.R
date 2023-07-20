setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggplot2)
library(dplyr)
library(matrixStats)
library(ggrepel)
library(DistributionUtils)
library(compositions)

load("repclr.RData")
#name B.ovatus subspecies
repclr$species_name[repclr$species_id == "OTU-04707"] <- "Bacteroides ovatus_1"
repclr$species_name[repclr$species_id == "OTU-04709"] <- "Bacteroides ovatus_2"
repclr$species_name[repclr$species_name == "Bacteroides vulgatus"] <- "Phocaeicola vulgatus"

head(repclr)

meta = read.csv("metadata.csv")


sub = data.frame(sample_id = meta$sample_id, Donor = meta$Donor,  time = meta$time)

ind = match(repclr$sample_id, sub$sample_id)
out = sub[ind,]


ordered = subset(repclr, select = c("sample_id","species_name", "reads", "clr", "log2_ptr"))
ordered = bind_cols(ordered, out)

subsamples = c("ae","an","ao","am") #these samples have more than 50 collections

shotgun = list()
for (i in 1:length(subsamples)) {
  temp = ordered[ordered$Donor == subsamples[i],]
  temp = temp[order(temp$time),]
  shotgun[[i]] = temp
  
  shotgun[[i]] = shotgun[[i]][!shotgun[[i]]$species_name == "",]
  
}
names(shotgun) = subsamples
save(shotgun, file = "shotgun.RData")