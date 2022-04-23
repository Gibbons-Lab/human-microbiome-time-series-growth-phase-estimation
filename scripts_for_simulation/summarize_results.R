summarize_results = function(df){
  temp1 = df[df$phase == "Lag",]
  temp2 = df[df$phase == "Acceleration",]
  temp3 = df[df$phase == "Exponential",]
  temp4 = df[df$phase == "Deceleration",]
  temp5 = df[df$phase == "Stationary",]
  
  #summarize coefficient of variation and correlation coefficient
  #got a bit tired
  lag.cv = summary(temp1$coef_var)
  lag.cc = summary(temp1$cor_coef)
  
  acc.cv = summary(temp2$coef_var)
  acc.cc = summary(temp2$cor_coef)
  
  exp.cv = summary(temp3$coef_var)
  exp.cc = summary(temp3$cor_coef)
  
  dec.cv = summary(temp4$coef_var)
  dec.cc = summary(temp4$cor_coef)
  
  sta.cv = summary(temp5$coef_var)
  sta.cc = summary(temp5$cor_coef)
  
  lag.cv = t(as.data.frame(as.matrix(lag.cv)))
  lag.cc = t(as.data.frame(as.matrix(lag.cc)))
  
  acc.cv = t(as.data.frame(as.matrix(acc.cv)))
  acc.cc = t(as.data.frame(as.matrix(acc.cc)))
  
  exp.cv = t(as.data.frame(as.matrix(exp.cv)))
  exp.cc = t(as.data.frame(as.matrix(exp.cc)))
  
  dec.cv = t(as.data.frame(as.matrix(dec.cv)))
  dec.cc = t(as.data.frame(as.matrix(dec.cc)))
  
  sta.cv = t(as.data.frame(as.matrix(sta.cv)))
  sta.cc = t(as.data.frame(as.matrix(sta.cc)))
  
  
  
  
  coef_var_summ = as.data.frame(bind_rows(append_seventh(lag.cv), append_seventh(acc.cv), 
                                          append_seventh(exp.cv), append_seventh(dec.cv), 
                                          append_seventh(sta.cv)))
  
  cor_coef_summ = as.data.frame(bind_rows(append_seventh(lag.cc), append_seventh(acc.cc), 
                                          append_seventh(exp.cc), append_seventh(dec.cc), 
                                          append_seventh(sta.cc)))
  
  rownames(coef_var_summ) = NULL
  coef_var_summ$phase = c("Lag", "Acceleration", "Exponential", "Deceleration", "Stationary")
  
  rownames(cor_coef_summ) = NULL
  cor_coef_summ$phase = c("Lag", "Acceleration", "Exponential", "Deceleration", "Stationary")
  
  return(list(coef_var_summ, cor_coef_summ))
}
