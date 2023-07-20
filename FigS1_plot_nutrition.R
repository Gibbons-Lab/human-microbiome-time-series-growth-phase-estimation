#plot whether diet is a stable factor over time
setwd("~/Desktop/2021_jlim_rotation/autoregressive_diet/")

library(openxlsx)
library(ggplot2)
library(tidyr)
library(R.utils)
library(tseries)

df1 = read.xlsx("13059_2013_3286_MOESM18_ESM.xlsx", sheet = 1)
df2 = read.xlsx("13059_2013_3286_MOESM18_ESM.xlsx", sheet = 2)

head(df1)
head(df2)


#Everything together
df2 = df2[,c(1,grep("food_Root", colnames(df2)),
            grep("nutrition", colnames(df2)),
            grep("diet", colnames(df2)))]
df2 = df2[3:nrow(df2),]

days = df2$day
df2 = gather(df2, condition, measurement, names(df2)[2:ncol(df2)], factor_key=TRUE)
df2$measurement = as.numeric(df2$measurement)
ggplot(df2, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`)) + theme(legend.position = "none")

#Food only
food = df2[grep("food", df2$condition),]
food$measurement = as.numeric(food$measurement)
ggplot(food, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`)) + theme(legend.position = "none")

lv = levels(factor(food$condition))

all = food[food$condition == lv[1],]
alcohol = food[food$condition == lv[3],]
sugary_drinks = food[food$condition == lv[16],]
starch = food[food$condition == lv[47],]
dairy = food[food$condition == lv[82],]
eggsnmeat = food[food$condition == lv[103],]
eggs = food[food$condition == lv[104],]
meat = food[food$condition == lv[107],]
nuts = food[food$condition == lv[125],]
seafood = food[food$condition == lv[135],]
fruits = food[food$condition == lv[152],]
sweets = food[food$condition == lv[193],]
veggie = food[food$condition == lv[226],]

#ADF and KPSS test on food items
lst = list(all, alcohol, sugary_drinks, starch, dairy, eggs, meat, nuts, seafood, fruits, sweets, veggie)
adf = list()
kpss = list()
for (i in 1:length(lst)) {
  temp = lst[[i]]
  temp = temp[complete.cases(temp$measurement),]
  
  adf[[i]] = adf.test(temp$measurement)
  kpss[[i]] = kpss.test(temp$measurement, "Trend")
}

labs = c("All Food", "Alcohol", "Sugary drinks", "Starch", "Dairy", "Eggs", 
         "Meat", "Nuts","Seafood", "Fuits","Sweets", "Vegetable")
for (i in 1:length(lst)) {
  temp = lst[[i]]
  temp1 = temp[complete.cases(temp$measurement),]
  
  p1=ggplot(data=temp, aes(x=`day`,y=`measurement`))+ 
    geom_line(size=2, shape=21, fill="black") + 
    theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.key = element_blank(),
          legend.text=element_text(size=0),axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=10),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 15, face = "bold"),plot.title = element_text(size = 20, face = "bold")) + 
    ylab("Measurement") + xlab("Days")+
    labs(title = labs[i])+
    geom_hline(yintercept = mean(temp1$measurement), colour = "firebrick1", linetype = 2)
  
  png(filename = paste0("Lineplot_",labs[i],"_throughout.png"), res = 300, width = 3, height = 3, units = "in")
  print(p1)
  dev.off()
}
#Plot as heatmap
df = data.frame(alcohol = alcohol$measurement, sugary_drinks = sugary_drinks$measurement, starch = starch$measurement,
                dairy = dairy$measurement, eggs = eggs$measurement, meat = meat$measurement, nuts = nuts$measurement, 
                seafood = seafood$measurement, fruits = fruits$measurement, sweets = sweets$measurement, vegetables = veggie$measurement)

library(ComplexHeatmap)
df = df[complete.cases(df),]
Heatmap(df, col = c("white", "red"), cluster_rows = F, cluster_columns = F, show_row_names = F)

#Nutrition
nutrition = df2[grep("nutrition", df2$condition),]
nutrition$measurement = as.numeric(nutrition$measurement)
nutrition$condition = gsub("nutrition_", "", nutrition$condition)
nutrition$condition = capitalize(tolower(nutrition$condition))
ggplot(nutrition, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`))# + theme(legend.position = "none")

##subset sodium
na = nutrition[grep("Sodium", nutrition$condition),]
ggplot(na, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`))# + theme(legend.position = "none")

delta = c()
for (i in 1:(nrow(na)-1)) {
  delta[i] = na$measurement[i+1] - na$measurement[i]
}
temp = data.frame(delta = delta, ref = na$measurement[1:nrow(na)-1])
plot(temp$ref, temp$delta)

##subset calories
cal = nutrition[grep("Calorie", nutrition$condition),]
ggplot(cal, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`))# + theme(legend.position = "none")

##subset calcium
ca = nutrition[grep("Calcium", nutrition$condition),]
ggplot(ca, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`))# + theme(legend.position = "none")

##Nutrition without calories, calcium, and sodium
ind = nutrition$condition == "Sodium" | nutrition$condition == "Calorie" | nutrition$condition == "Calcium"
rest = nutrition[!ind,]
ggplot(rest, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`))# + theme(legend.position = "none")

#KPSS and ADF test on nutrition
nutrition$condition[nutrition$condition == "Carb"] <- "Carbohydrate"
nutrition$condition[nutrition$condition == "Satfat"] <- "Saturated fat"
lv = levels(factor(nutrition$condition))
res1 = list()
res2 = list()
adf1 = list()
adf2 = list()
for (i in 1:length(lv)) {
  name = lv[i]
  temp = nutrition[grep(paste0(name), nutrition$condition),]
  
  temp1 = temp[temp$day < 71,]
  res1[[i]] = kpss.test(temp1$measurement, null = "Trend")
  temp1 = temp1[complete.cases(temp1),]
  adf1[[i]] = adf.test(temp1$measurement)
  
  
  p1=ggplot(data=temp1, aes(x=`day`,y=`measurement`))+ 
    geom_line(size=2, shape=21, fill="black") + 
    theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.key = element_blank(),
          legend.text=element_text(size=0),axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=10),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 15, face = "bold"),plot.title = element_text(size = 20, face = "bold")) + 
    ylab("Measurement") + xlab("Days")+
    labs(title = lv[i])+
    geom_hline(yintercept = mean(temp1$measurement), colour = "firebrick1", linetype = 2)
  
  
  png(filename = paste0("Lineplot_",lv[i],"_before.png"), res = 300, width = 3, height = 3, units = "in")
  print(p1)
  dev.off()
  
  temp2 = temp[temp$day > 122,]
  res2[[i]] = kpss.test(temp2$measurement, null = "Trend")
  temp2 = temp2[complete.cases(temp2),]
  adf2[[i]] = adf.test(temp2$measurement, )
  
  
  p2=ggplot(data=temp2, aes(x=`day`,y=`measurement`))+ 
    geom_line(size=2, shape=21, fill="black") + 
    theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.key = element_blank(),
          legend.text=element_text(size=0),axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=10),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 15, face = "bold"),plot.title = element_text(size = 20, face = "bold")) + 
    ylab("Measurement") + xlab("Days")+
    labs(title = lv[i])+
    geom_hline(yintercept = mean(temp2$measurement), colour = "firebrick1", linetype = 2, size = 2)
  
  png(filename = paste0("Lineplot_",lv[i],"_after.png"), res = 300, width = 3, height = 3, units = "in")
  print(p2)
  dev.off()
}
names(res1) = lv
names(adf1) = lv
names(res2) = lv
names(adf2) = lv

#KPSS and ADF test on nutrition -- before and after combined
nutrition$condition[nutrition$condition == "Carb"] <- "Carbohydrate"
nutrition$condition[nutrition$condition == "Satfat"] <- "Saturated fat"

lv = levels(factor(nutrition$condition))
res3 = list()
adf3 = list()
for (i in 1:length(lv)) {
  name = lv[i]
  temp = nutrition[grep(paste0(name), nutrition$condition),]
  
  temp1 = temp[temp$day < 123,]
  temp2 = temp[temp$day > 122,]
  temp2 = temp2[complete.cases(temp2$measurement),]
  temp3 = rbind(temp1, temp2)
  temp4 = temp3[complete.cases(temp3),]
  
  p1 = ggplot(data=temp3, aes(x=`day`,y=`measurement`))+ 
    geom_line(size=2, shape=21, fill="black") + 
    theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.key = element_blank(),
          legend.text=element_text(size=0),axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=10),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 15, face = "bold"),plot.title = element_text(size = 20, face = "bold")) + 
    ylab("Measurement") + xlab("Days")+
    labs(title = lv[i]) +
    geom_hline(yintercept = mean(temp4$measurement), colour = "firebrick1", linetype = 2)
  
  
  png(filename = paste0("Lineplot_",lv[i],"_before_after_combined.png"), res = 300, width = 3, height = 3, units = "in")
  print(p1)
  dev.off()
  
  res3[[i]] = kpss.test(temp3$measurement, null = "Trend")
  adf3[[i]] = adf.test(temp4$measurement)
  
}
  


#Diet
diet = df2[grep("diet", df2$condition),]
diet$measurement = as.numeric(diet$measurement)
ggplot(diet, aes(x = `day`, y = `measurement`)) + geom_line(aes(group = `condition`, color = `condition`)) + theme(legend.position = "none")
