setwd("~/Desktop/2021_jlim_rotation/")

source("For_Github_upload/scripts/scripts_for_simulation/plot_multi_histogram.R")
source("For_Github_upload/scripts/scripts_for_simulation/summarize_results.R")

plot_corr_coef = function(df, sliced = F){
  require(dplyr)
  require(ggplot2)
  
  if (sliced == T){
    temp = df
    
    bind = lapply(temp, function(x) do.call(bind_rows, x))
    bind = bind_rows(bind)
  }
  if (sliced == F) {
    bind = df
  }
  
  
  
  #separate phases - correlation coefficient
  e = as.data.frame(bind[seq(3, nrow(bind), 5),4])
  s = as.data.frame(bind[seq(5, nrow(bind), 5),4])
  a = as.data.frame(bind[seq(2, nrow(bind), 5),4])
  d = as.data.frame(bind[seq(4, nrow(bind), 5),4])
  l = as.data.frame(bind[seq(1, nrow(bind), 5),4])
  
  #name columns occurance and phase
  colnames(e) = 'occ'
  colnames(s) = 'occ'
  colnames(a) = 'occ'
  colnames(d) = 'occ'
  colnames(l) = 'occ'
  
  e$phase = 'Mid-log'
  s$phase = 'Stationary'
  a$phase = 'Acceleration'
  d$phase = 'Deceleration'
  l$phase = 'Acceleration'
  
  all = rbind(l, a, e, d, s)
  all = all[complete.cases(all$occ),]
  
  
  all$phase = factor(all$phase, levels = c("Acceleration", "Mid-log", "Deceleration", "Stationary"))
  
  p1 <- ggplot(all, aes(x = factor(`phase`), y = `occ`, fill = `phase`)) +
    ggdist::stat_halfeye(
      adjust = .5,
      width = .6, 
      ## set slab interval to show IQR and 95% data range
      .width = c(.05, .95)
    ) + 
    
    coord_cartesian(xlim = c(1.2, NA))+
    #geom_boxplot(alpha = 0.80,  size = 2, color = "grey") +
    #geom_dotplot(binaxis='y', stackdir='center')+
    xlab("") + ylab(expression(paste("Pearson's r"))) + 
    theme(axis.text.y=element_text(size=8), 
          legend.text = element_text(size = 0),
          axis.text.x=element_text(size=12, face = "bold", angle = 45, vjust = 1, hjust = 1),
          axis.title.y=element_text(size=15, face = "bold"),
          legend.position=c(10,1),legend.justification=c(1,1),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")
    )+
    scale_fill_manual(values = c("orange", "grey70", "royalblue", "navy"))
  
  
  
  return(p1)
}

#load simulation results
#outputs from completed simulations saved as different RData files. 
load("coef_var_r_noise_extlow.RData") #noise from 0.001 to 0.01
load("coef_var_r_noise_low.RData") #noise from 0.01 to 0.1
load("coef_var_r_noise_high.RData") #noise from 0.1 to 1

combined = c(coef_var_r_noise_extlow, coef_var_r_noise_low, coef_var_r_noise_high)

#Fig. 5B

# r = 1.2, k = 100, noise = 0.1
df = coef_var_r_noise_high[[1]][[2]][[6]]
p1 = plot_corr_coef(df, sliced = T)

png(filename = "Pearson's_r_all_phases_noise=0.1_r=1.2_k=100.png", width = 4.2, height = 4.2, units = "in", res = 300)
print(p1)
dev.off()



#Fig. S1B
#noise from 0.001 to 1
temp = list()
ind = c(1,2,3,5,6,7,9,10,11,12)
for (i in 1:length(ind)) {
  bind = lapply(combined[[i]], function(x) do.call(bind_rows, x))
  bind = bind_rows(bind)
  temp[[i]] = bind
}

df = bind_rows(temp)
p1 = plot_corr_coef(df, sliced = F)

png(filename = "Pearson's_r_all_phases_all_noise.png", width = 4.2, height = 4.2, units = "in", res = 300)
print(p1)
dev.off()







