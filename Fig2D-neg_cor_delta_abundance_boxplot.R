#negative correlation between delta vs abundance as boxplots for each donor/time series
setwd("~/Desktop/2021_jlim_rotation/11-12-21_copy/metagenomics/")
library(ggpmisc)
library(broom)
library(ggplot2)


load("all_clr.RData")
#name B.ovatus subspecies
all_clr$species_name[all_clr$species_id == "OTU-04707"] <- "Bacteroides ovatus_1"
all_clr$species_name[all_clr$species_id == "OTU-04709"] <- "Bacteroides ovatus_2"

head(all_clr)

meta = read.csv("metadata.csv")


sub = data.frame(sample_id = meta$sample_id, Donor = meta$Donor,  time = meta$time)

ind = match(all_clr$sample_id, sub$sample_id)
out = sub[ind,]


ordered = subset(all_clr, select = c("sample_id","species_name", "reads", "clr"))
ordered = bind_cols(ordered, out)

subsamples = c("ae","an","ao","am") #these samples have more than 50 collections

all = list()
for (i in 1:length(subsamples)) {
    temp = ordered[ordered$Donor == subsamples[i],]
    temp = temp[order(temp$time),]
    all[[i]] = temp
    
    all[[i]] = all[[i]][!all[[i]]$species_name == "",]
    
}
names(all) = subsamples

#calculate delta (abundance at t+delta t -abundance at t)
meta_b_all = list()
b_all = list()
taxon_save = list()
for (i in 1:length(subsamples)) {
    
    mf = data.frame(clr = all[[i]]$clr, species = all[[i]]$species_name, time = all[[i]]$time)
    
    #remove taxa that appears less than 5 times for all time
    mf = mf[mf$species %in% names(which(table(mf$species) > 5)), ]
    
    len = levels(factor(mf$species))
    len_new = len
    
    #calculate delta for each taxon
    meta_del = list()
    mean_clr = list()
    betas = list()
    saved = list()
    
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
        
        
        #only calculate correlation coefficients between abundance and delta if taxon is frequently detected
        if (nrow(temp)>20) {
            #estimate carrying capacity
            fit = lm(delta~clr, data = temp)
            
            coeff = fit$coefficients
            
            beta = as.numeric(coeff[2])
            
            
            
            betas[[j]] = beta
            
            #stop and check if correlation coefficient is positive - this indicates that distribution of delta is not normal
            if (beta >= 0) {
                break
            }
            mean_clr[[j]] = mean(temp$clr)
            deltas[[j]] = temp$delta
            
            #calculate x intercept
            fit = lm(delta~clr, data = temp)
            
            coeff = fit$coefficients
            
            K = as.numeric(-coeff[1]/coeff[2])
            
            
            
            meta_del[[j]] = K
            
            mean_clr[[j]] = mean(temp$clr)

            saved[[j]] = temp
            
        }
        
        else{
            len_new[j] <- NA
            next
        }
        
        
    }
    
    
    
    taxon_save[[i]] = saved
    cor_coef = data.frame(beta = unlist(betas))
    c = data.frame(mean_clr = unlist(mean_clr))
    
    b = bind_cols(cor_coef,c)
    rownames(b) = len_new[complete.cases(len_new)]
    meta_b_all[[i]] = b
    
    b_all[[i]] = unlist(betas)
    
    
}




#check results - all good

summary(meta_b_all[[1]]$beta)
summary(meta_b_all[[2]]$beta)
summary(meta_b_all[[3]]$beta)
summary(meta_b_all[[4]]$beta)



#plot results as a boxplot
don1 = data.frame(group = "ae", value = meta_b_all[[1]]$beta)
don2 = data.frame(group = "an", value = meta_b_all[[2]]$beta)
don3 = data.frame(group = "ao", value = meta_b_all[[3]]$beta)
don4 = data.frame(group = "am", value = meta_b_all[[4]]$beta)

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

png(filename = "Cor_coef_delta_clr_4_timeseries.png", width = 4, height = 4.5, res = 300, units = "in")
print(p1)
dev.off()

