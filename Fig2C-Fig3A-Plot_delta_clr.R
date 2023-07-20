setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggpmisc)
library(broom)
library(ggplot2)
library(dplyr)

load("shotgun.RData")
#calculate delta (abundance at t+delta t -abundance at t)
meta_k = list()
d = list()
pd = position_dodge(0.1)
subsamples = c("ae","an","ao","am") #these samples have more than 50 collections

df_list = list()
for (i in 1:length(subsamples)) {
  mf = data.frame(clr = shotgun[[i]]$clr, log2_ptr = shotgun[[i]]$log2_ptr, species = shotgun[[i]]$species_name, time = shotgun[[i]]$time)
  
  #remove taxa that appears less than 5 times for all time
  mf = mf[mf$species %in% names(which(table(mf$species) > 5)), ]
  
  len = levels(factor(mf$species))
  len_new = len
  
  #calculate delta for each taxon
  meta_del = list()
  mean_clr = list()
  mean_ptr = list()
  deltas = list()
  dfs = list()
  
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
    temp$donor = rep(subsamples[i], nrow(temp))
    
    if (nrow(temp)>1) {
      #estimate carrying capacity
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
      deltas[[j]] = temp$delta
      
      #estimate average
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
      dfs[[j]] = temp
    }
    
    else{
      len_new[j] <- NA
      next
    }
    
    
  }
  
  
  
  
  k = data.frame(k = unlist(meta_del))
  c = data.frame(mean_clr = unlist(mean_clr), mean_ptr = unlist(mean_ptr))
  
  k = bind_cols(k,c)
  rownames(k) = len_new[complete.cases(len_new)]
  meta_k[[i]] = k
  
  d[[i]] = unlist(deltas)
  df_list[[i]] = dfs
  
}



df_list = bind_rows(df_list)
bacs = unique(df_list$species)

xsmallest = list()
xlargest = list()
ysmallest = list()
ylargest = list()
ptrsmallest = list()
ptrlargest = list()
for (i in 1:length(bacs)) {
  ind <- df_list$species == bacs[i]
  
  xsmallest[[i]] = min(df_list$clr[ind])
  xlargest[[i]] = max(df_list$clr[ind])
  ysmallest[[i]] = min(df_list$clr[ind])
  ylargest[[i]] = max(df_list$clr[ind])
  ptrsmallest[[i]] = min(df_list$log2_ptr[ind])
  ptrlargest[[i]] = max(df_list$log2_ptr[ind])
}



ranges = data.frame(taxa = bacs, xmin = unlist(xsmallest), xmax = unlist(xlargest), ymin = unlist(ysmallest), ymax = unlist(ylargest),
                    ptrmin = unlist(ptrsmallest), ptrmax = unlist(ptrlargest))

load("cv_list.RData")
for (i in 1:length(subsamples)) {
  mf = data.frame(clr = shotgun[[i]]$clr, log2_ptr = shotgun[[i]]$log2_ptr, species = shotgun[[i]]$species_name, 
                  time = shotgun[[i]]$time, donor = shotgun[[i]]$Donor)
  
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
 
    
    
    
    if (nrow(temp)>1) {
      #estimate carrying capacity
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
      deltas[[j]] = temp$delta
      
      
      
      
      #clr abundance vs delta for each taxon
      formula = y ~ x
      
      p1=ggplot(data=temp, aes(x=clr,y=`delta`))+ 
        geom_point(position=pd, size=3, shape=21, fill="black") + 
        xlab("CLR Abundance (t1)") +
        ylab("CLR Delta (t2-t1)") +
        theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
              legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
              legend.text=element_text(size=0),axis.text.x = element_text(size=10.5),
              axis.text.y=element_text(size=11),legend.title=element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.title = element_text(size = 14, face = "bold"))+
        labs(title = len[j])+
        geom_smooth(method = "lm", formula = formula, se = F) #+
        #xlim(ranges$xmin[ranges$taxa == len[j]], ranges$xmax[ranges$taxa == len[j]])+
        #ylim(ranges$ymin[ranges$taxa == len[j]], ranges$ymax[ranges$taxa == len[j]])
      #stat_poly_eq(formula = formula, 
      #            aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
      #           parse = TRUE)+
      #stat_fit_glance(method = "lm",
      #               label.y = "bottom",
      #              method.args = list(formula = y ~ x),
      #             mapping = aes(label = sprintf('~italic(p)~"="~%.2g',
      #                                          stat(p.value))),
      #           parse = TRUE)
      
      
      
      png(paste0("Figures/clr_vs_delta_noaxes/",subsamples[i],"/",len[j],".png"), 
          width = 3, height = 3, units = "in", res = 300)
      print(p1)
      dev.off()    
      
      
      
      #clr abundance vs delta for each taxon
      formula = y ~ x
      
      p1=ggplot(data=temp, aes(x=clr,y=`delta`))+ 
        geom_point(position=pd, size=3, shape=21, fill="black") + 
        xlab("CLR Abundance (t1)") +
        ylab("CLR Delta (t2-t1)") +
        theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
              legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
              legend.text=element_text(size=0),axis.text.x = element_text(size=10.5),
              axis.text.y=element_text(size=11),legend.title=element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.title = element_text(size = 14, face = "bold"))+
        labs(title = len[j])+
        geom_smooth(method = "lm", formula = formula, se = F) #+
      #xlim(ranges$xmin[ranges$taxa == len[j]], ranges$xmax[ranges$taxa == len[j]])+
      #ylim(ranges$ymin[ranges$taxa == len[j]], ranges$ymax[ranges$taxa == len[j]])
      #stat_poly_eq(formula = formula, 
      #            aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
      #           parse = TRUE)+
      #stat_fit_glance(method = "lm",
      #               label.y = "bottom",
      #              method.args = list(formula = y ~ x),
      #             mapping = aes(label = sprintf('~italic(p)~"="~%.2g',
      #                                          stat(p.value))),
      #           parse = TRUE)
      
      
      
      png(paste0("Figures/clr_vs_delta_with_axes/",subsamples[i],"/",len[j],".png"), 
          width = 3, height = 3, units = "in", res = 300)
      print(p1)
      dev.off()    
      
      #estimate average
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
      #capacity_sub = capacity_sub[1:j]##
    }
    
    else{
      len_new[j] <- NA
      next
    }
    
    
  }
  
  k = data.frame(k = unlist(meta_del))
  c = data.frame(mean_clr = unlist(mean_clr), mean_ptr = unlist(mean_ptr))
  
  k = bind_cols(k,c)
  rownames(k) = len_new[complete.cases(len_new)]
  meta_k[[i]] = k
  
  d[[i]] = unlist(deltas)
  
  
}


#Fig3A
load("cv_list.RData")
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
      fit = lm(delta~clr, data = temp)
      
      coeff = fit$coefficients
      
      k = as.numeric(-coeff[1]/coeff[2])
      
      
      
      meta_del[[j]] = k
      
      mean_clr[[j]] = mean(temp$clr)
      mean_ptr[[j]] = mean(temp$log2_ptr)
      deltas[[j]] = temp$delta
      
      
      temp$t = 1:nrow(temp)
      
      p1=ggplot(data=temp, aes(x=clr,y=log2_ptr))+ 
        geom_point(position=pd, size=3, shape=21, fill="black") +
        theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
              legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
              legend.text=element_text(size=0),axis.text.x = element_text(size=10.5),
              axis.text.y=element_text(size=11),legend.title=element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.title = element_text(size = 14, face = "bold"))+
        xlab("CLR Abundance") +
        ylab("log2(PTR)") 
      
      #if (len[j] %in% ranges$taxa) {
      #  p1 = p1+xlim(ranges$xmin[ranges$taxa == len[j]], ranges$xmax[ranges$taxa == len[j]])+
      #    ylim(ranges$ptrmin[ranges$taxa == len[j]], ranges$ptrmax[ranges$taxa == len[j]])
      #}
        
      if (len[j] %in% coef$Taxon) {
        if (coef$adj_p[j] < 0.05 & coef$Slope[j] > 0) {
          p1 = p1 + geom_smooth(method = "lm", formula = formula, color = "orange", se = F)
        }else if (coef$adj_p[j] < 0.05 & coef$Slope[j] < 0) {
          p1 = p1 + geom_smooth(method = "lm", formula = formula, color = "royalblue", se = F)
        }else{
          p1 = p1 + geom_smooth(method = "lm", formula = formula, color = "grey70", se = F)
        }
      }
      
      
      png(paste0("Figures/clr_vs_ptr_noaxes/",subsamples[i],"/",taxon_name[j],"_colored.png"), 
          width = 3, height = 3, units = "in", res = 300)
      print(p1)
      dev.off()    
      
      
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
