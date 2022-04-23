#coefficient of variation and correlation coefficient
CV_cor = function(df, da, kurt){
  
  if (kurt == F) {
    std = sd(df)
    mu = mean(df)
    
    coef_var = std/mu
    
    cor_coef = cor(da, df)
    
    coef_var_cor = data.frame("var" = std^2, "mu" = mu, "coef_var" = coef_var, "cor_coef" = cor_coef)#, "sample" = 1:100)
  }
  if (kurt == T) {
    std = sd(df)
    mu = mean(df)
    
    coef_var = std/mu
    
    cor_coef = cor(da, df)
    
    kurt = kurtosis(df, na.rm = T)
    
    coef_var_cor = data.frame("var" = std^2, "mu" = mu, "coef_var" = coef_var, "cor_coef" = cor_coef, "kurt" = kurt)#, "sample" = 1:100)
  }
  
  
  return(coef_var_cor)
}
