plot_sim = function(df, df1, df2, df3, lag = lag, acceleration = acceleration, exponential = exponential, deceleration = deceleration, acl = acl, tmp = tmp, tmp2 = tmp2){
  p1=ggplot(data=df, aes(x=`t`,y=`sol`))+ 
    geom_line(aes(group = group))+
    theme(legend.key.width=unit(0.15,"cm"),legend.key.height=unit(0.3,"cm"),
          legend.justification=c(0,1),panel.spacing=unit(-0.07,"cm"),legend.position="top",
          legend.text=element_text(size=15),axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=0),legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          title = element_text(size = 14, face = "bold"), ) + 
    ylab("Abundance") + xlab("Time")+
    labs(title = paste0("Simulated growth r=", r, ", k=", k, ", n=",w))+
    geom_vline(xintercept = lag[length(lag)], linetype = "dashed", size = 2)+
    geom_vline(xintercept = acceleration[length(acceleration)], linetype = "dashed", size = 2)+
    geom_vline(xintercept = exponential[length(exponential)], linetype = "dashed", size = 2)+
    geom_vline(xintercept = deceleration[length(deceleration)], linetype = "dashed", size = 2)
  

  
  p2=ggplot(data=df1, aes(x=`t`,y=`sol`))+ 
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
  
 
  

  p3=ggplot(data=df2, aes(x=`t`,y=`mean_da`))+ 
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
  
  tmp = abs(df3$acl-0.5*max(df3$acl))
  tmp2 = abs(df3$acl-0.5*min(df3$acl))
  
  p4 = ggplot(df3, aes(x = t, y = `acl`)) +
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
 
  
  return(list(p1, p2, p3, p4))
}
