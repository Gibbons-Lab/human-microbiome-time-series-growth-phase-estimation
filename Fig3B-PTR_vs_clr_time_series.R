#growth rates changes 
setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggpmisc)
library(broom)
library(ggplot2)
library(dplyr)

load("shotgun.RData")

clr_ptr = list()
for(x in 1:length(shotgun)){
  
  mf = data.frame(clr = shotgun[[x]]$clr, log2_ptr = shotgun[[x]]$log2_ptr, 
                  species = shotgun[[x]]$species_name, time = shotgun[[x]]$time)
  
  #remove taxa that appears less than 5 times for all time
  sdf = mf[mf$species %in% names(which(table(mf$species) > 5)), ]
  
  lev = levels(factor(mf$species))

  betas = list()
  cors = list()
  for (i in 1:length(lev)) {
    if (nrow(sdf[sdf$species == lev[i],])<2) {
      next
    }
    tdf = sdf[sdf$species == lev[i],]
    std = sd(tdf$clr)
    mu = mean(tdf$clr)
    
    std.r = sd(tdf$log2_ptr)
    mu.r = mean(tdf$log2_ptr)
    
    fit = lm(log2_ptr~clr, data = tdf)
    
    coeff = fit$coefficients
    
    beta = as.numeric(coeff[2])

    betas[[i]] = data.frame(norm_var_a = (std^2)/mu, var_a = std^2, mu_a = mu, beta = beta, 
                            cor = cor(x = tdf$clr, y = tdf$log2_ptr, method = "pearson"))
    
    
    
  }
  cdf = bind_rows(betas)
  
  clr_ptr[[x]] = cdf
  
  
}

names(clr_ptr) = c("ae", "an", "ao", "am")




#plot results as a boxplot
don1 = data.frame(group = "ae", value = clr_ptr[[1]]$beta)
don2 = data.frame(group = "an", value = clr_ptr[[2]]$beta)
don3 = data.frame(group = "ao", value = clr_ptr[[3]]$beta)
don4 = data.frame(group = "am", value = clr_ptr[[4]]$beta)

plot_data = bind_rows(don1, don2, don3, don4)
require(latex2exp)

p1 = ggplot(plot_data, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(alpha = 0.80,  size = 2, color = "black", fill = "white") +
  #geom_point(size = 3, shape = 21,alpha = 0) +
  #geom_jitter(size = 2, alpha = 0.3, width = 0.2)+
  xlab("") + ylab("Linear regression coefficient")+ #ylab(TeX("Cor. coef $\\beta_1$")) + 
  theme(axis.text.y=element_text(size=15), 
        legend.text = element_text(size = 0),
        axis.text.x=element_text(size=25, face = "bold"),
        axis.title.y=element_text(size=20, face = "bold"),
        legend.position="None",legend.justification=c(1,1),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
  )+
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", size = 2)

png(filename = "Cor_coef_PTR_clr_4_timeseries.png", width = 4, height = 4.5, res = 300, units = "in")
print(p1)
dev.off()


#plot results as a boxplot - Pearson coefficient
don1 = data.frame(group = "ae", value = clr_ptr[[1]]$cor)
don2 = data.frame(group = "an", value = clr_ptr[[2]]$cor)
don3 = data.frame(group = "ao", value = clr_ptr[[3]]$cor)
don4 = data.frame(group = "am", value = clr_ptr[[4]]$cor)

plot_data = bind_rows(don1, don2, don3, don4)
require(latex2exp)

p1 = ggplot(plot_data, aes(x = group, y = value, fill = group)) + 
  geom_boxplot(alpha = 0.80,  size = 2, color = "black", fill = "white") +
  #geom_point(size = 3, shape = 21,alpha = 0) +
  #geom_jitter(size = 2, alpha = 0.3, width = 0.2)+
  xlab("") + ylab("Pearson r")+ #ylab(TeX("Cor. coef $\\beta_1$")) + 
  theme(axis.text.y=element_text(size=15), 
        legend.text = element_text(size = 0),
        axis.text.x=element_text(size=25, face = "bold"),
        axis.title.y=element_text(size=20, face = "bold"),
        legend.position="None",legend.justification=c(1,1),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
  )+
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", size = 2)

png(filename = "Pearson_r_PTR_clr_4_timeseries.png", width = 4, height = 4.5, res = 300, units = "in")
print(p1)
dev.off()

#plot results as a boxplot - Pearson coefficient with all four donors combined - Fig. SX.
don1 = data.frame(group = "ae", value = clr_ptr[[1]]$cor)
don2 = data.frame(group = "an", value = clr_ptr[[2]]$cor)
don3 = data.frame(group = "ao", value = clr_ptr[[3]]$cor)
don4 = data.frame(group = "am", value = clr_ptr[[4]]$cor)

combined = do.call("bind_rows", list(don1, don2, don3, don4))
combined$group = rep("Donors", nrow(combined))

p1 = ggplot(combined, aes(x = group, y = value)) + 
  geom_boxplot(alpha = 0.80,  size = 2, color = "black", fill = "white") +
  #geom_point(size = 3, shape = 21,alpha = 0) +
  #geom_jitter(size = 2, alpha = 0.3, width = 0.2)+
  xlab("") + ylab("Pearson r")+ #ylab(TeX("Cor. coef $\\beta_1$")) + 
  theme(axis.text.y=element_text(size=15), 
        legend.text = element_text(size = 0),
        axis.text.x=element_text(size=25, face = "bold"),
        axis.title.y=element_text(size=20, face = "bold"),
        legend.position="None",legend.justification=c(1,1),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
  )+
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", size = 2)

png(filename = "Pearson_r_PTR_clr_all_donors_combined.png", width = 3, height = 4, res = 300, units = "in")
print(p1)
dev.off()
