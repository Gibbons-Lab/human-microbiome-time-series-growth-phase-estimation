#plot number of donations/sample collections
setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggpmisc)
library(broom)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

load("shotgun.RData")

don_ae = data.frame(clr = shotgun[[1]]$clr, log2_ptr = shotgun[[1]]$log2_ptr, species = shotgun[[1]]$species_name, time = shotgun[[1]]$time)
don_an = data.frame(clr = shotgun[[2]]$clr, log2_ptr = shotgun[[2]]$log2_ptr, species = shotgun[[2]]$species_name, time = shotgun[[2]]$time)
don_ao = data.frame(clr = shotgun[[3]]$clr, log2_ptr = shotgun[[3]]$log2_ptr, species = shotgun[[3]]$species_name, time = shotgun[[3]]$time)
don_am = data.frame(clr = shotgun[[4]]$clr, log2_ptr = shotgun[[4]]$log2_ptr, species = shotgun[[4]]$species_name, time = shotgun[[4]]$time)


time_ae = ifelse(1:max(don_ae$time) %in% unique(don_ae$time), "1", "0")
time_an = ifelse(1:max(don_an$time) %in% unique(don_an$time), "1", "0")
time_ao = ifelse(1:max(don_ao$time) %in% unique(don_ao$time), "1", "0")
time_am = ifelse(1:max(don_am$time) %in% unique(don_am$time), "1", "0")

Heatmap(time_ae, show_row_names = F, show_column_names = F, width = unit(2,"mm"), name = " ", 
        col = c("grey90", "red"),
        rect_gp = gpar(col = "white", lwd = 0.8), show_heatmap_legend = F)

dev.copy(png, "Time_heatmap_ae.png", width = 5, height = 5,units = "in", res = 300)
dev.off()

Heatmap(time_an, show_row_names = F, show_column_names = F, width = unit(2,"mm"), name = " ", 
        col = c("grey90", "red"),
        rect_gp = gpar(col = "white", lwd = 0.8), show_heatmap_legend = F)

dev.copy(png, "Time_heatmap_an.png", width = 5, height = 5,units = "in", res = 300)
dev.off()

Heatmap(time_ao, show_row_names = F, show_column_names = F, width = unit(2,"mm"), name = " ", 
        col = c("grey90", "red"),
        rect_gp = gpar(col = "white", lwd = 0.8), show_heatmap_legend = F)

dev.copy(png, "Time_heatmap_ao.png", width = 5, height = 5,units = "in", res = 300)
dev.off()

Heatmap(time_am, show_row_names = F, show_column_names = F, width = unit(2,"mm"), name = " ", 
        col = c("grey90", "red"),
        rect_gp = gpar(col = "white", lwd = 0.8), show_heatmap_legend = F)

dev.copy(png, "Time_heatmap_am.png", width = 5, height = 5,units = "in", res = 300)
dev.off()

dev.off()




lgd = Legend(labels = c("Absence", "Presence"),title = "Collection", nrow = 1, title_position = "topcenter",
             graphics = list(
               function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "grey90", col = "grey90")),
               function(x, y, w, h) grid.rect(x, y, w, h, gp = gpar(fill = "red", col = "red"))
             ))

draw(lgd)

dev.copy(png, "Time_heatmap_legend.png", width = 5, height = 5,units = "in", res = 300)
dev.off()
