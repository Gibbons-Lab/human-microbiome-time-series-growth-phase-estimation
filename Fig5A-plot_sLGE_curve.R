setwd("~/Desktop/2021_jlim_rotation/")
source("scripts_for_simulation/append_seventh.R")
source("scripts_for_simulation/coef_var_cor.R")
source("scripts_for_simulation/CV_cor.R")
source("scripts_for_simulation/second_der.R")
source("scripts_for_simulation/growth_phases.R")
source("scripts_for_simulation/plot_sim.R")
source("scripts_for_simulation/summarize_results.R")

library(sde)
library(ggplot2)
library(vegan)
library(reshape)
library(tidyr)
library(ggpmisc)
library(broom)
library(dplyr)

#sLGE simulation
r = 1.2
k = 100
w = 0.1
set.seed(123)
plot_figs = T
df = coef_var_cor(r = r, k = k, w = w, plot_figs = plot_figs)

