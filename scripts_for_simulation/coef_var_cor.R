coef_var_cor = function(r, k, w, plot_figs = F){
  #simulate stochastic growth
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
  
  #calcuate dA/dt
  da = list()
  for (i in 1:nrow(wide_sim)) {
    da[[i]] = as.data.frame(wide_sim[i+1,] - wide_sim[i,])
  }
  da = bind_rows(da)
  da = da[-nrow(da),]
  
  if( plot_figs == T ){
    #plot growth and acceleration showing each growth phase
    mean_da = rowMeans(da)
    mean_da = data.frame(mean_da = mean_da, t = 1:100)
    
    l_a = c(lag, acceleration)
    lp = exponential
    d_s = c(deceleration, stationary)
    
    sim$col = ifelse(sim$t %in% l_a, "L_A", 
                     ifelse(sim$t %in% lp, "Log", 
                            ifelse(sim$t %in% d_s, "D_S", "L_A")))
    
    
    p1=ggplot(data=sim, aes(x=`t`,y=`sol`))+ 
      geom_line(aes(group = group, color = sim$col))+
      theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
            legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="none",
            legend.text=element_text(size=15),axis.text.x = element_text(size=10),
            axis.text.y=element_text(size=0),legend.title=element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            title = element_text(size = 14, face = "bold"), ) + 
      ylab("Abundance") + xlab("Time")+
      labs(title = paste0("Simulated growth"))+
      scale_color_manual(values = c("royalblue", "orange", "grey70"))
      #geom_vline(xintercept = lag[length(lag)], linetype = "dashed", size = 2)+
      #geom_vline(xintercept = acceleration[length(acceleration)], linetype = "dashed", size = 2)+
      #geom_vline(xintercept = exponential[length(exponential)], linetype = "dashed", size = 2)+
      #geom_vline(xintercept = deceleration[length(deceleration)], linetype = "dashed", size = 2)
    
    
    
    p2=ggplot(data=lg, aes(x=`t`,y=`sol`))+ 
      geom_line(size = 3)+
      theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
            legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
            legend.text=element_text(size=15),axis.text.x = element_text(size=10),
            axis.text.y=element_text(size=0),legend.title=element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            title = element_text(size = 14, face = "bold"), ) + 
      ylab("Abundance") + xlab("Time")+
      labs(title = paste0("Nonstochastic growth r=", r, ", k=", k, ", n=",w))+
      geom_vline(xintercept = lag[length(lag)], linetype = "dashed", size = 2)+
      geom_vline(xintercept = acceleration[length(acceleration)], linetype = "dashed", size = 2)+
      geom_vline(xintercept = exponential[length(exponential)], linetype = "dashed", size = 2)+
      geom_vline(xintercept = deceleration[length(deceleration)], linetype = "dashed", size = 2)
    
    
    
    
    p3=ggplot(data=mean_da, aes(x=`t`,y=`mean_da`))+ 
      geom_line()+
      theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
            legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
            legend.text=element_text(size=15),axis.text.x = element_text(size=10),
            axis.text.y=element_text(size=0),legend.title=element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            title = element_text(size = 15, face = "bold"), ) + 
      ylab("dA/dt") + xlab("Time")+
      labs(title = "Mean growth rate")+
      geom_smooth(method = "loess")
    
    tmp = abs(acl$acl-0.5*max(acl$acl))
    tmp2 = abs(acl$acl-0.5*min(acl$acl))
    
    p4 = ggplot(acl, aes(x = t, y = `acl`)) +
      geom_line(color = "steelblue", size = 3)+
      #geom_point(aes(y = acl), colour = "navy", size = 3)+
      xlab("Time") +
      ylab("Acceleration") +
      theme_bw()+
      theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
            legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
            legend.text=element_text(size=15),axis.text.x = element_text(size=10),
            axis.text.y=element_text(size=0),legend.title=element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"), 
            title = element_text(size = 14, face = "bold"))+
      geom_hline(yintercept = acl$acl[which(tmp == min(abs(acl$acl-0.5*max(acl$acl))))], linetype = "dashed", size = 2)+
      geom_hline(yintercept = acl$acl[which(tmp2 == min(abs(acl$acl-0.5*min(acl$acl))))], linetype = "dashed", size = 2)+
      labs(title = paste0("Acc. curve r=", r, ", k=", k, ", n=",w))
    
    
    plots = list(p1, p2, p3, p4)
    
    plot_names = c("Simulated_growth", "Nonstochastic_growth", "Growth_rate_curve", "Acc_curve")
    for (a in 1:length(plots)) {
      png(paste0(plot_names[a]," r = ", r, ", k = ", k, ", n = ", w, ".png"), width = 5, height = 5, res = 300, units = "in")
      print(plots[[a]])
      dev.off()
    }
  }
  
  
  #calculate coefficient of variation and correlation coefficient
  calc = wide_sim[-nrow(wide_sim),]
  
  #geo = list()
  #for (z in 1:ncol(calc)) {
  #  geo[[z]] = log(as.data.frame(calc[,z]/geometricmean(calc[,z])))
  #}
  
  #calc = bind_cols(geo)
  
  
  prop = list()
  for (z in 1:ncol(calc)) {
    prop[[z]] = as.data.frame(log1p((calc[,z]/geometricmean(calc[,z]))))
  }
  
  calc = bind_cols(prop)
  
  #calc = log(calc)
  
  result = list()
  
  for (i in 1:ncol(wide_sim)) {
    coef_var1 = data.frame(CV_cor(calc[lag,i], da[lag,i], kurt = T))
    
    coef_var2 = data.frame(CV_cor(calc[acceleration,i], da[acceleration,i], kurt = T))
    
    coef_var3 = data.frame(CV_cor(calc[exponential,i], da[exponential,i], kurt = T))
    
    coef_var4 = data.frame(CV_cor(calc[deceleration,i], da[deceleration,i], kurt = T))
    
    coef_var5 = data.frame(CV_cor(calc[stationary,i], da[stationary,i], kurt = T))
    
    
    coef_var = dplyr::bind_rows(coef_var1, coef_var2,coef_var3, coef_var4, coef_var5)
    coef_var$phase = c("Lag", "Acceleration", "Exponential", "Deceleration", "Stationary")
    coef_var$phase = factor(coef_var$phase, levels = coef_var$phase, ordered = T)
    
    result[[i]] = coef_var
  }
  
  
  
  
  
  return(result)
}
