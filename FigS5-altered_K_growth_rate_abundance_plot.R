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

setwd("~/Desktop/2021_jlim_rotation/")
source("For_Github_upload/scripts/scripts_for_simulation/append_seventh.R")
source("For_Github_upload/scripts/scripts_for_simulation/coef_var_cor.R")
source("For_Github_upload/scripts/scripts_for_simulation/CV_cor.R")
source("For_Github_upload/scripts/scripts_for_simulation/second_der.R")
source("For_Github_upload/scripts/scripts_for_simulation/growth_phases.R")
source("For_Github_upload/scripts/scripts_for_simulation/plot_sim.R")
source("For_Github_upload/scripts/scripts_for_simulation/summarize_results.R")

##Iterate each 100 times
library(dplyr)
library(ggplot2)
source("~/Desktop/2021_jlim_rotation/For_Github_upload/scripts/scripts_for_simulation/stochastic_K.R")
source("~/Desktop/2021_jlim_rotation/For_Github_upload/scripts/scripts_for_simulation/stochastic_K_X.R")

#medium r, medium k
k = 100
r = 1.2
sim = list()
for (i in 1:100) {
  temp = stochastic_K(r = r, k = k, N = 100)
  temp$group = i
  sim[[i]] = temp
  
}
sim = dplyr::bind_rows(sim)

#simulate non-stochastic growth 
k = 100
r = 1.2

lg = stochastic_K(r = r, k = k, mean = 0, sd = 0, N = 100)
lg$group = rep(i, 101)


#convert to wide table
wide_sim = reshape(sim, timevar = "group", idvar = c("t"), direction = "wide")
wide_sim = wide_sim[,-1]

#calculate growth phases of averaged growth from acceleration curve
calc = lg[-nrow(lg),]

#deterministic dxdt
gr = c()
for (i in 1:(nrow(lg)-1)) {
  gr[i] = lg$sol[i+1] - lg$sol[i]
}
gr = data.frame(t = 1:100, sol = gr, group = lg$group[1:(nrow(lg)-1)])

acl = second_der(calc, k)
acl = data.frame(t = 1:length(acl$sol), acl = acl$sol)

phase = growth_phases(acl)

lag = phase[[1]]
acceleration = phase[[2]]
exponential = phase[[3]]
deceleration = phase[[4]]
stationary = phase[[5]]

#calcuate dA/dt
da = list()
for (i in 1:nrow(wide_sim)) {
  da[[i]] = as.data.frame(wide_sim[i+1,] - wide_sim[i,])
}
da = bind_rows(da)
da = da[-nrow(da),]


#plot growth and acceleration showing each growth phase
mean_da = rowMeans(da)
mean_da = data.frame(mean_da = mean_da, t = 1:100)

l_a = c(lag, acceleration)
lp = exponential
d = deceleration
s = stationary

sim$col = ifelse(sim$t %in% l_a, "Acceleration", 
                 ifelse(sim$t %in% lp, "Mid-log", 
                        ifelse(sim$t %in% d, "Deceleration", 
                               ifelse(sim$t %in% s, "Stationary", "Acceleration"))))

sim$col = factor(sim$col, levels = c("Acceleration", "Mid-log", "Deceleration", "Stationary"), ordered = T)

scaleFactor <- max(sim$sol) / max(gr$sol)
p1=ggplot(data=sim, aes(x=`t`))+ 
  geom_line(aes(y=`sol`, group = group, color = `col`))+
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="none",
        legend.text=element_text(size=15),axis.text.x = element_text(size=10),
        axis.text.y=element_text(size=10),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 14, face = "bold")) + 
  xlab("Time")+
  labs(title = paste0("Simulated growth varying k,x"))+
  scale_color_manual(values = c("orange", "grey70", "royalblue", "navy"))+
  geom_line(data = lg, aes(y = `sol`), color = "orangered", size = 1, linetype = 2)+
  geom_line(data = gr, aes(y = `sol`*scaleFactor), color = "black", size = 1, linetype = 2)+
  scale_y_continuous(
    name = "Abundance", 
    sec.axis = sec_axis(~./scaleFactor, name="Growth rate")
  )
p1



df = data.frame(dx = unlist(da[l_a,]), x = unlist(wide_sim[l_a,]))
df = log1p(df)
df = df[complete.cases(df$dx),]
p=ggplot(data=df, aes(x=`x`,y=`dx`))+ 
  geom_point(size=1, shape=21, fill="black") +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.text=element_text(size=15),axis.text.x = element_text(size=0),
        axis.text.y=element_text(size=0),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=20, face = "bold")) + 
  ylab("dx/dt") + xlab("x")+
  geom_smooth(method = "lm", formula = y~x, colour = "orange", size = 2, se = F)

png(filename = paste0("altered_K_lag_acceleration_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()


df = data.frame(dx = unlist(da[exponential,]), x = unlist(wide_sim[exponential,]))
df = log1p(df)
df = df[complete.cases(df$dx),]
p=ggplot(data=df, aes(x=`x`,y=`dx`))+ 
  geom_point(size=1, shape=21, fill="black") +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.text=element_text(size=15),axis.text.x = element_text(size=0),
        axis.text.y=element_text(size=0),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=20, face = "bold")) + 
  ylab("dx/dt") + xlab("x")+
  geom_smooth(method = "lm", formula = y~x, colour = "grey70", size = 2, se = F)

png(filename = paste0("altered_K_log_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()


df = data.frame(dx = unlist(da[deceleration,]), x = unlist(wide_sim[deceleration,]))
df = log1p(df)
df = df[complete.cases(df$dx),]
p=ggplot(data=df, aes(x=`x`,y=`dx`))+ 
  geom_point(size=1, shape=21, fill="black") +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.text=element_text(size=15),axis.text.x = element_text(size=0),
        axis.text.y=element_text(size=0),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=20, face = "bold")) + 
  ylab("dx/dt") + xlab("x")+
  geom_smooth(method = "lm", formula = y~x, colour = "royalblue", size = 2, se = F)

png(filename = paste0("altered_K_deceleration_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()



df = data.frame(dx = unlist(da[stationary,]), x = unlist(wide_sim[stationary,]))
df = log1p(df)
df = df[complete.cases(df$dx),]
p=ggplot(data=df, aes(x=`x`,y=`dx`))+ 
  geom_point(size=1, shape=21, fill="black") +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.text=element_text(size=15),axis.text.x = element_text(size=0),
        axis.text.y=element_text(size=0),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=20, face = "bold")) + 
  ylab("dx/dt") + xlab("x")+
  geom_smooth(method = "lm", formula = y~x, colour = "navy", size = 2, se = F)

png(filename = paste0("altered_K_stationary_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()


df = data.frame(dx = unlist(gr$sol[deceleration]), x = unlist(lg$sol[deceleration]))
ggplot(data=df, aes(x=`x`,y=`dx`))+ 
  geom_point(size=1, shape=21, fill="black") +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.text=element_text(size=15),axis.text.x = element_text(size=10),
        axis.text.y=element_text(size=10),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=20, face = "bold")) + 
  ylab("dx/dt") + xlab("x")+
  geom_smooth(method = "lm", formula = y~x, colour = "royalblue", size = 2, se = F)


df = data.frame(dx = unlist(gr$sol[stationary]), x = unlist(lg$sol[stationary]))
ggplot(data=df, aes(x=`x`,y=`dx`))+ 
  geom_point(size=1, shape=21, fill="black") +
  theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
        legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
        legend.text=element_text(size=15),axis.text.x = element_text(size=10),
        axis.text.y=element_text(size=10),legend.title=element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        title = element_text(size = 16, face = "bold"),
        axis.title=element_text(size=20, face = "bold")) + 
  ylab("dx/dt") + xlab("x")+
  geom_smooth(method = "lm", formula = y~x, colour = "navy", size = 2, se = F)

