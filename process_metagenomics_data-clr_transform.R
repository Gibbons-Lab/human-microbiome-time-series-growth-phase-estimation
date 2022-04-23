setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggplot2)
library(dplyr)
library(matrixStats)
library(ggrepel)
library(DistributionUtils)
library(compositions)

metdf = read.csv("reads_and_rates.csv", stringsAsFactors = F)

#only those with PTR
n = length(levels(factor(metdf$sample_id)))
result = list()
for (i in 1:n) {
  temp = metdf[metdf$sample_id == levels(factor(metdf$sample_id))[i],]
  temp = temp[temp$reads>0.1,] #remove 0 reads to avoid making geometric mean to be 0
  temp = temp[complete.cases(temp$log2_ptr),] #include only taxa with replication rates
  geo = geometricmean(as.numeric(temp$reads))
  temp$clr = temp$reads/geo
  temp$clr = log(temp$clr)
  
  
  result[[i]] = temp
}
repclr = bind_rows(result)
save(repclr, file = "repclr.RData")

#all taxa even without PTR
n = length(levels(factor(metdf$sample_id)))
result = list()
for (i in 1:n) {
  temp = metdf[metdf$sample_id == levels(factor(metdf$sample_id))[i],]
  temp = temp[temp$reads>0.1,] #remove 0 reads to avoid making geometric mean to be 0
  geo = geometricmean(as.numeric(temp$reads))
  temp$clr = temp$reads/geo
  temp$clr = log(temp$clr)
  
  
  result[[i]] = temp
}
all_clr = bind_rows(result)
save(all_clr, file = "all_clr.RData")
