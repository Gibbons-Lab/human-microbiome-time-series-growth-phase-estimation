#Simulate and show range of Pearson r values for each growth phase
setwd("~/Desktop/2021_jlim_rotation/")

library(sde)
library(ggplot2)
library(vegan)
library(reshape)
library(tidyr)
library(ggpmisc)
library(broom)
library(dplyr)
library(matrixStats)
library(ggrepel)

source("For_Github_upload/scripts/scripts_for_simulation/append_seventh.R")
source("For_Github_upload/scripts/scripts_for_simulation/coef_var_cor.R")
source("For_Github_upload/scripts/scripts_for_simulation/CV_cor.R")
source("For_Github_upload/scripts/scripts_for_simulation/second_der.R")
source("For_Github_upload/scripts/scripts_for_simulation/growth_phases.R")
source("For_Github_upload/scripts/scripts_for_simulation/plot_sim.R")
source("For_Github_upload/scripts/scripts_for_simulation/summarize_results.R")

set.seed(123)
r = 1.2
k = 100
w = 0.2
plot_figs = F

res = list()
temp_all = list()
for (i in 1:50) {
  df = coef_var_cor(r = r, k = k, w = w, plot_figs = plot_figs)
  temp = bind_rows(df)
  temp_all[[i]] = temp
  res[[i]] = summarize_results(temp)[[2]]
}

temp = temp_all
bind = lapply(temp, function(x) do.call(bind_rows, x))
bind = bind_rows(bind)

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
all$all = rep("Simulation", nrow(all))

#Plot entire range
ggplot(all, aes(x = all, y = occ)) + 
  geom_boxplot(alpha = 0.80,  size = 2, color = "black", fill = "white") +
  xlab("") + ylab("Pearson r")+ 
  theme(axis.text.y=element_text(size=15), 
        legend.text = element_text(size = 0),
        axis.text.x=element_text(size=25, face = "bold"),
        axis.title.y=element_text(size=20, face = "bold"),
        legend.position="None",legend.justification=c(1,1),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
  )+
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed", size = 2)

dev.copy(png, "Pearson_r_entire_growth_phase_1.2_100_0.2.png", width = 3, height = 4, res = 300, unit = "in")
dev.off()


#Plot distribution from simulation
all$phase = factor(all$phase, 
                   levels = c("Acceleration", "Mid-log", "Deceleration", "Stationary"),
                   ordered = T)
ggplot(all, aes(x = factor(`phase`), y = `occ`, fill = `phase`)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6, 
    .width = c(.10, .90)
  ) + 
  coord_cartesian(xlim = c(1.2, NA))+
  xlab("") + ylab(expression(paste("Pearson r"))) + 
  theme(axis.text.y=element_text(size=8), 
        legend.text = element_text(size = 0),
        axis.text.x=element_text(size=12, face = "bold", angle = 45, vjust = 1, hjust = 1),
        axis.title.y=element_text(size=15, face = "bold"),
        legend.position=c(10,1),legend.justification=c(1,1),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
  )+
  scale_fill_manual(values = c("orange", "grey70", "royalblue", "navy"))

dev.copy(png, "Pearson_r_separated_by_growth_phase_1.2_100_0.2.png", 
    width = 4.5, height = 3, units = "in", res = 300)
dev.off()

