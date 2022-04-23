setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggpmisc)
library(broom)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

load("shotgun.RData")

#calculate regression coefficient of ptr vs clr
betas = list()
for (i in 1:length(subsamples)) {
  mf = data.frame(clr = shotgun[[i]]$clr, log2_ptr = shotgun[[i]]$log2_ptr, species = shotgun[[i]]$species_name, time = shotgun[[i]]$time)
  
  #remove taxa that appears less than 5 times for all time
  sdf = mf[mf$species %in% names(which(table(mf$species) > 5)), ]
  len = levels(factor(sdf$species))

  #for consistency 
  len_new = len
  len = len[order(len)]
  beta_donor = list()
  
  #calculate delta for each taxon
  for (j in 1:length(len)) {
    temp = mf[mf$species == len[j],]
    temp = temp[!duplicated(temp$time),]
    
    del = list()
    time_del = list()
    
    for (k in 1:(length(temp$time)-1)) {
      del[[k]] = temp$clr[k+1] - temp$clr[k]
      time_del[[k]] = temp$time[k+1] - temp$time[k]
    }
    
    del = as.data.frame(t(bind_cols(del)))
    time_del = as.data.frame(t(bind_cols(time_del)))
    
    index = time_del$V1>3
    del = del[!index,]
    
    temp = temp[-nrow(temp),]
    temp = temp[!index,]
    temp$delta = unlist(del) #combine in "temp" data frame (contains one taxon within sample)
    
    if (nrow(temp)>1) {

      
      #calculate regression coeff
      fit = lm(log2_ptr ~ clr, data = temp)
      coeff = fit$coefficients
      b = as.numeric(coeff[2])
      summarized = summary(fit)
      p = as.numeric(summarized$coefficients[,4][2])
      beta_donor[[j]] = data.frame(Sample = subsamples[i], Taxon = len[j], Slope = b, R2 = as.numeric(summarized$r.squared), p_val = p)

    }

    
    
  }
  
  
  betas[[i]] = as.data.frame(bind_rows(beta_donor))
  
}






#for color coding trendline
#library(openxlsx)
ae = betas[[1]]
an = betas[[2]]
ao = betas[[3]]
am = betas[[4]]

cv_list = list(ae, an, ao, am)
save(cv_list, file = "cv_list.RData")

pd = position_dodge(0.1)
meta_k = list()
d = list()
for (i in 1:length(subsamples)) {
  mf = data.frame(clr = shotgun[[i]]$clr, log2_ptr = shotgun[[i]]$log2_ptr, species = shotgun[[i]]$species_name, time = shotgun[[i]]$time)
  
  #remove taxa that appears less than 5 times for all time
  sdf = mf[mf$species %in% names(which(table(mf$species) > 5)), ]
  len = levels(factor(sdf$species))
  
  #for consistency 
  len_new = len
  len = len[order(len)]
  beta_donor = list()
  
  #for consistency 
  len_new = len
  len = len[order(len)]
  
  #adjust regression p values
  coef = cv_list[[i]]
  coef$adj_p = p.adjust(coef$p_val, method = "fdr")
  
  coef = coef[order(coef$Taxon),]
  rownames(coef) = NULL
  
  #calculate delta for each taxon
  meta_del = list()
  mean_clr = list()
  mean_ptr = list()
  deltas = list()
  temp2 = list()
  for (j in 1:length(len)) {
    temp = mf[mf$species == len[j],]
    temp = temp[!duplicated(temp$time),]
    
    del = list()
    time_del = list()
    DAC = list()
    
    for (k in 1:(length(temp$time)-1)) {
      del[[k]] = temp$clr[k+1] - temp$clr[k]
      time_del[[k]] = temp$time[k+1] - temp$time[k]
    }
    del = as.data.frame(t(bind_cols(del)))
    time_del = as.data.frame(t(bind_cols(time_del)))
    
    index = time_del$V1>3
    del = del[!index,]
    
    temp = temp[-nrow(temp),]
    temp = temp[!index,]
    temp$delta = unlist(del) #combine in "temp" data frame (contains one taxon within sample)
    
    temp2[[j]] = temp
  }
  
  temp2 = Filter(function(x) nrow(x)>1, temp2)
  taxon_name = rep(0, length(temp2))
  for (j in 1:length(temp2)) {
    temp = temp2[[j]]
    if (nrow(temp)>1) {
      
      taxon_name[j] = temp$species[1]
      #x int
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
      deltas[[j]] = temp$delta
      
      
      
      #calculate x intercept
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
    } else{
      len_new[j] <- NA
      next
    }
    
  }
  
  
  k = data.frame(k = unlist(meta_del))
  c = data.frame(mean_clr = unlist(mean_clr), mean_ptr = unlist(mean_ptr))
  
  k = bind_cols(k,c)
  rownames(k) = taxon_name
  meta_k[[i]] = k
  
  d[[i]] = unlist(deltas)
  
  
  
}
    

###plot as heatmap
library(ComplexHeatmap)
library(circlize)

for (i in 1:length(betas)) {
  temp = betas[[i]]
  temp$adj_p = p.adjust(temp$p_val, method = "fdr")
  temp$sig = ifelse(temp$Slope>0 & temp$adj_p < 0.05, "Sig. pos",
                    ifelse(temp$Slope<0 & temp$adj_p < 0.05, "Sig. neg", "Not sig."))
  
  betas[[i]] = temp
  
}

#combine as one data frame
bind_rows(betas)

ae = betas[[1]]
an = betas[[2]]
ao = betas[[3]]
am = betas[[4]]


taxa = union(betas[[1]]$Taxon, betas[[2]]$Taxon)
taxa = union(taxa, betas[[3]]$Taxon)
taxa = union(taxa, betas[[4]]$Taxon)
taxa = taxa[order(taxa)]


df_for_hm = function(df, taxa){
  
  df = df[,c(2,7)]
  diff = setdiff(taxa, df$Taxon)
  diff = data.frame(Taxon = diff, sig = rep("Not detected", length(diff)))
  df = bind_rows(df, diff)
  df = df[order(df$Taxon),]
  df = df[,2]
  
  return(df)
}

ae = df_for_hm(ae, taxa)  
am = df_for_hm(am, taxa)  
an = df_for_hm(an, taxa)  
ao = df_for_hm(ao, taxa)  

for_heatmap = data.frame(taxa = taxa, ae = ae, am = am, an = an, ao = ao)
rownames(for_heatmap) = for_heatmap$taxa
for_heatmap = for_heatmap[,-1]

for_heatmap = as.matrix(for_heatmap)


(p1 = Heatmap(for_heatmap, show_row_names = T, show_column_names = T, name = "Significant", c("floralwhite", "grey90", "royalblue", "orange"),
        rect_gp = gpar(col = "white", lwd = 0.5), row_names_gp = gpar(fontsize = 12, fontface = "bold"),
        column_names_gp = gpar(fontsize = 18, fontface = "bold"), column_names_rot = 0, column_names_centered = T))

png(filename = "significant_correlations_PTR_CLR.png", width = 6, height = 6, units = "in", res = 300)
print(p1)
dev.off()
