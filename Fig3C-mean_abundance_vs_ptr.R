setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggpmisc)
library(broom)
library(ggplot2)


load("shotgun.RData")


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
      #calculate x intercept
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
      #capacity_sub = capacity_sub[1:j]##
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

#mean abundance vs mean ptr
for (i in 1:length(meta_k)) {
  temp = meta_k[[i]]
  meta_k[[i]]$donor = subsamples[i]
  
}

temp = bind_rows(meta_k)
summary(lm(mean_clr~0+mean_ptr+donor, temp))

summary(lm(mean_clr~0+mean_ptr+donor, temp))

p1=ggplot(data=temp, aes(x=`k`,y=`mean_ptr`,color = donor))+ 
  geom_point(size=3, alpha = 0.2) +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.key = element_blank(),
        legend.text=element_text(size=20),axis.text.x = element_text(size=10),
        axis.text.y=element_text(size=10),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 20, face = "bold"),plot.title = element_text(size = 15, face = "bold")) + 
  ylab("Mean log2(PTR)") + xlab("Mean CLR Abundance")+ 
  geom_smooth(method = "lm", se = F, show.legend = F)+
  geom_smooth(method = "lm", se = F)+
  #scale_color_manual(values = c("orangered", "saddlebrown", "peru", "red3"))
  #scale_color_manual(values = c("burlywood", "chocolate4", "gold2", "orangered"))
  scale_color_manual(values = c("dodgerblue3", "#00AFBB", "gold2", "orangered"))
print(p1)

png(paste0("mean_ptr_clr_4_donors.png"), width = 4, height = 4, units = "in", res = 300)
print(p1)
dev.off()






for (i in 1:length(meta_k)) {
  temp = meta_k[[i]]
  p1=ggplot(data=temp, aes(x=`k`,y=`mean_ptr`))+ 
    geom_point(position=pd, size=3, shape=21, fill="black") +
    theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.text=element_text(size=10),axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=10),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 12, face = "bold"),plot.title = element_text(size = 15, face = "bold")) + 
    ylab("log2(PTR)") + xlab("Normalized Abundance")+
    labs(title = subsamples[i])
  
  p1 = p1 + geom_smooth(method = "lm", formula = formula, color = "orange", se = F)
  
  png(paste0("Figures/k_vs_ptr/", subsamples[i], "_v2.png"))
  print(p1)
  dev.off()
}

