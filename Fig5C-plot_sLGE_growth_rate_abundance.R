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

r = 1.2
k = 100
w = 0.1


sim = list()
for (i in 1:100) {
  f = expression(r*((k-x)/k)*x)
  s = expression(w*x)
  sol = sde.sim(t0=0, T=10, N=100, X0=1, drift=f, sigma=s)
  t = seq(0, 100, length.out=101)  # 100 intervals, 101 boundary values
  
  df = data.frame(sol = sol, t = t, group = rep(i, 101))
  
  sim[[i]] = df
}

sim = dplyr::bind_rows(sim)

plot(sim$t, sim$sol)

#simulate non-stochastic growth 
f = expression(r*((k-x)/k)*x)
s = expression(0*x)
sol = sde.sim(t0=0, T=10, N=100, X0=1, drift=f, sigma=s)
t = seq(0, 100, length.out=101)  # 100 intervals, 101 boundary values

lg = data.frame(sol = sol, t = t, group = rep(i, 101))

#convert to wide table
wide_sim = reshape(sim, timevar = "group", idvar = c("t"), direction = "wide")
wide_sim = wide_sim[,-1]

#calculate growth phases of averaged growth from acceleration curve
calc = lg[-nrow(lg),]


acl = second_der(calc, k)
acl = data.frame(t = 1:length(acl$sol), acl = acl$sol)


phase = growth_phases(acl)

lag = phase[[1]]
acceleration = phase[[2]]
exponential = phase[[3]]
deceleration = phase[[4]]
stationary = phase[[5]]

l_a = c(lag, acceleration)
#d_s = c(deceleration, stationary)

#calcuate dA/dt
da = list()
for (i in 1:nrow(wide_sim)) {
  da[[i]] = as.data.frame(wide_sim[i+1,] - wide_sim[i,])
}
da = bind_rows(da)
da = da[-nrow(da),]

#Identify whole range of dx
min = as.numeric(summary(unlist(da))[1])
max = as.numeric(summary(unlist(da))[6])

df = data.frame(dx = unlist(da[l_a,]), x = unlist(wide_sim[l_a,]))
df = log1p(df)
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
  #ylim(min = min*1.05, max = max*1.05)+
  geom_smooth(method = "lm", formula = y~x, colour = "orange", size = 2, se = F)

png(filename = paste0("sde_lag_acceleration_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()


df = data.frame(dx = unlist(da[exponential,]), x = unlist(wide_sim[exponential,]))
df = log1p(df)
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

png(filename = paste0("sde_log_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()


df = data.frame(dx = unlist(da[deceleration,]), x = unlist(wide_sim[deceleration,]))
df = log1p(df)
min = as.numeric(summary(df$dx)[2])
max = as.numeric(summary(df$dx)[5])

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
  #ylim(min*1.1, max*1.1)+
  geom_smooth(method = "lm", formula = y~x, colour = "royalblue", size = 2, se = F)

png(filename = paste0("sde_deceleration_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()



df = data.frame(dx = unlist(da[stationary,]), x = unlist(wide_sim[stationary,]))
df = log1p(df)
min = as.numeric(summary(df$dx)[2])
max = as.numeric(summary(df$dx)[5])

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
  #ylim(min*1.1, max*1.1)+
  geom_smooth(method = "lm", formula = y~x, colour = "navy", size = 2, se = F)

png(filename = paste0("sde_stationary_0.1noise.png"), width = 3, height = 3, units = "in", res = 400)
p
dev.off()
