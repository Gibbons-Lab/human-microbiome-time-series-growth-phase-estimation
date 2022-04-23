plot_multi_histogram <- function(df, feature, label_column, bins, xlab, ylab, density = T, color) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.4, position="identity", aes(y = ..density..), color="black", bins = bins) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")+
    theme(axis.text.y=element_text(size=10), 
          legend.text = element_text(size = 8),
          axis.text.x=element_text(size=10, face = "bold"),#, angle = 45, vjust = 1, hjust = 1),
          axis.title.y=element_text(size=15, face = "bold"),
          axis.title.x = element_text(size = 15, face = "bold"),
          legend.position="top",legend.justification="top",legend.title=element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          #plot.title = element_text(size = 30, face = "bold")
          legend.key.size = unit(0.3, "cm")
          )+
    scale_fill_manual(values = alpha(color, 0.2), breaks = c(levels(factor(df[,2]))))+
    xlab(paste(xlab))+
    ylab(paste(ylab))
  plt + guides(fill=guide_legend(title=label_column))
  
  if (density == T) {
    plt+geom_density(alpha=0.5)
  }
  
  
}
