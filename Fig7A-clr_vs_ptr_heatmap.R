setwd("~/Desktop/2021_jlim_rotation")
library(ggpmisc)
library(broom)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
#combine as one data frame
load("for_viewing.rds")

ae = for_viewing[[1]]
an = for_viewing[[2]]
ao = for_viewing[[3]]
am = for_viewing[[4]]

#Rename B. vulgatus
ao$Taxon[ao$Taxon == "Bacteroides vulgatus"] <- "Phocaeicola vulgatus"

taxa = unique(bind_rows(for_viewing)$Taxon)
taxa = taxa[order(taxa)]

df_for_hm = function(df, taxa){
  
  df = df[,c(2,7)]
  diff = setdiff(taxa, df$Taxon)
  diff = data.frame(Taxon = diff, assignment = rep("Not detected", length(diff)))
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

colors = c("Not detected" = "floralwhite", "Acceleration" = "orange", "Non-stationary" = "grey80", 
           "Deceleration" = "royalblue", "Stationary" = "navy")

column_ha = HeatmapAnnotation('Avg defecation/day' = anno_barplot(c(2.5, 2.5, 1.5, 1)))
(p1 = Heatmap(for_heatmap, show_row_names = T, show_column_names = T, name = "Phase", col = colors,
              rect_gp = gpar(col = "white", lwd = 0.5), row_names_gp = gpar(fontsize = 12, fontface = "bold"),
              column_names_gp = gpar(fontsize = 18, fontface = "bold"), column_names_rot = 0, column_names_centered = T,
              top_annotation = column_ha))

png(filename = "Fig7A_significant_correlations_PTR_CLR.png", width = 6, height = 6, units = "in", res = 300)
print(p1)
dev.off()
