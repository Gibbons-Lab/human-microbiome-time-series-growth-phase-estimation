growth_phases = function(df){
  
  if (sum(df$acl == Inf) > 0 | sum(df$acl == -Inf) > 0) {
    #remove |inf| and borrow information from adjacent time points
    ind <- abs(df$acl) == Inf
    df = df[!ind,]
    
    tmp = abs(df$acl-0.5*max(df$acl))
    tmp2 = abs(df$acl-0.5*min(df$acl))
    
  }else{
    tmp = abs(df$acl-0.5*max(df$acl))
    tmp2 = abs(df$acl-0.5*min(df$acl))
  }
  
  if (length(tmp) == 0 | length(tmp2) == 0) {
    lag = NULL
    acceleration = NULL
    exponential = NULL
    deceleration = NULL
    stationary = NULL
    
  }else{
    half.max = df$acl[which(tmp == min(abs(df$acl-0.5*max(df$acl))))]
    half.min = df$acl[which(tmp2 == min(abs(df$acl-0.5*min(df$acl))))]
    
    half.max.line = df$acl - half.max
    half.min.line = df$acl - half.min
    
    half.max.points = order(abs(half.max.line))
    half.min.points = order(abs(half.min.line))
    
    ordered.max.acl = df$acl[half.max.points]
    ordered.min.acl = df$acl[half.min.points]
    
    if (length(ordered.min.acl)==0 & length(ordered.max.acl)!=0) {
      ##should be no first or second half min
      first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                         ordered.max.acl >= half.max][1]
      second.half.max = half.max.points[half.max.points > which(df$acl == max(df$acl)) & 
                                          ordered.max.acl >= half.max][1]
      
      lag = 1:first.half.max
      acceleration = first.half.max+1:second.half.max
      exponential = ((second.half.max+1)):nrow(df)
      deceleration = NULL
      stationary = NULL
    }
    
    
    
    #from wolfram alpha for solving for x for intersection point - didn't work
    #y = half.max
    #x=0.50000*(0.48075*(1.7321*sqrt((108*y^2 - k^6)) + 18*y)^(1/3) + (0.69336*k^2)/(1.7321*sqrt((108*y^2 - k^6)) + 18*y)^(1/3) + k)
    
    ######have half.max.points[2] be later than the peak after half.max.points[1]
    #same for half.min.points[2]
    if ( df$acl[1] < round(max(df$acl)) ) {
      
      if( r >= k ){
        first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                           ordered.max.acl >= half.max][1]
        first.half.min = half.min.points[half.min.points <= which(df$acl == min(df$acl)) & 
                                           ordered.min.acl <= half.min][1]
        
        second.half.max = half.max.points[half.max.points > which(df$acl == max(df$acl)) & 
                                            ordered.max.acl >= half.max][1]
        second.half.min = half.min.points[half.min.points > which(df$acl == min(df$acl)) & 
                                            ordered.min.acl <= half.min][1]
      }else if ( which(df$acl == max(df$acl)) > half.max.points[1] & 
                 (half.max.points > which(df$acl == max(df$acl)) & ordered.max.acl >= half.max)[1] == F & 
                 df$acl[half.max.points[half.max.points > which(df$acl == max(df$acl))][1]] <= half.max ) {
        #growth rate is too high -> time between half max/min points not separated
        
        first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                           ordered.max.acl >= half.max][1]
        second.half.max =  which(df$acl == max(df$acl))+1
        
      }else{
        first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                           ordered.max.acl >= half.max][1]
        first.half.min = half.min.points[half.min.points <= which(df$acl == min(df$acl)) & 
                                           ordered.min.acl <= half.min][1]
        
        second.half.max = half.max.points[half.max.points > which(df$acl == max(df$acl)) & 
                                            ordered.max.acl >= half.max][1]
        second.half.min = half.min.points[half.min.points > which(df$acl == min(df$acl)) & 
                                            ordered.min.acl <= half.min][1]
      }
      
      if ( which(df$acl == min(df$acl)) < half.min.points[1] & 
           (half.min.points < which(df$acl == min(df$acl)) & ordered.min.acl <= half.min)[1] == F & 
           df$acl[half.min.points[half.min.points < which(df$acl == min(df$acl))][1]] >=half.min & 
           (which(df$acl <= half.min))[1] == which(df$acl == min(df$acl)) ) {
        #growth rate is too high -> time between half max/min points not separated
        
        first.half.min =  (which(df$acl >= min(df$acl) & df$acl <= half.min)[1]-1)
        second.half.min =   half.min.points[half.min.points >= which(df$acl == min(df$acl)) & 
                                              ordered.min.acl <= half.min & half.min.points > second.half.max][1]
        
      }else{
        first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                           ordered.max.acl >= half.max][1]
        first.half.min = half.min.points[half.min.points <= which(df$acl == min(df$acl)) & 
                                           ordered.min.acl <= half.min][1]
        
        half.max.points[2] = half.max.points[half.max.points > which(df$acl == max(df$acl)) & 
                                               ordered.max.acl >= half.max][1]
        second.half.min = half.min.points[half.min.points > which(df$acl == min(df$acl)) & 
                                            ordered.min.acl <= half.min][1]
      }
    }else if( df$acl[1] >= round(max(df$acl)) ){
      if( r >= k ){
        first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                           ordered.max.acl >= half.max][1]
        first.half.min = half.min.points[half.min.points <= which(df$acl == min(df$acl)) & 
                                           ordered.min.acl <= half.min][1]
        
        half.max.points[2] = half.max.points[half.max.points > which(df$acl == max(df$acl)) & 
                                               ordered.max.acl >= half.max][1]
        second.half.min = half.min.points[half.min.points > which(df$acl == min(df$acl)) & 
                                            ordered.min.acl <= half.min][1]
      }else if ( which(df$acl == max(df$acl)) > half.max.points[1] & (half.max.points > which(df$acl == max(df$acl)) & ordered.max.acl >= half.max)[1] == F & 
                 df$acl[half.max.points[half.max.points > which(df$acl == max(df$acl))][1]] <= half.max ) {
        #growth rate is too high -> time between half max/min points not separated
        
        first.half.max = half.max.points[half.max.points <= which(df$acl == max(df$acl)) & 
                                           ordered.max.acl >= half.max][1]
        second.half.max = which(df$acl == max(df$acl))+1
        
      }else{
        first.half.max = 0
        first.half.min = half.min.points[half.min.points < which(df$acl == min(df$acl))][1]
        
        second.half.max = half.max.points[half.max.points > which(df$acl == max(df$acl)) & df$acl >= half.max][1]
        second.half.min = half.min.points[half.min.points > which(df$acl == min(df$acl)) & df$acl <= half.min][1]
      }
      if ( which(df$acl == min(df$acl)) < half.min.points[1] & 
           (half.min.points < which(df$acl == min(df$acl)) & ordered.min.acl <= half.min)[1] == F & 
           df$acl[half.min.points[half.min.points < which(df$acl == min(df$acl))][1]] >=half.min & 
           (which(df$acl <= half.min))[1] == which(df$acl == min(df$acl))  ) {
        #growth rate is too high -> time between half max/min points not separated
        
        first.half.min = (which(df$acl >= min(df$acl) & df$acl <= half.min)[1]-1)
        
        second.half.min = half.min.points[half.min.points >= which(df$acl == min(df$acl)) & ordered.min.acl <= half.min & half.min.points > second.half.max][1]
        
        
      }else{
        first.half.max = 0
        first.half.min = half.min.points[half.min.points < which(df$acl == min(df$acl))][1]
        
        second.half.max = half.max.points[half.max.points > which(df$acl == max(df$acl)) & df$acl >= half.max][1]
        second.half.min = half.min.points[half.min.points > which(df$acl == min(df$acl)) & df$acl <= half.min][1]
      }
    }
    
    #define growth phases
    if ( r >= k ) {
      #if growth passes carrying capacity, taxon is always at carrying capacity
      lag = NULL
      acceleration = NULL
      exponential = NULL
      deceleration = NULL
      stationary = 1:nrow(df)
    }else if ( (df$acl[1] >= half.max & df$acl[1] != 0 & first.half.min < which(df$acl == min(df$acl))) & 
               second.half.min > which(df$acl == min(df$acl)) | 
               (first.half.max == 0 & first.half.min < which(df$acl == min(df$acl)) &  
                second.half.min > which(df$acl == min(df$acl)) & 
                first.half.min >= second.half.max & is.na(second.half.max) == F) ) {
      #if there's no lag phase
      lag = NULL
      acceleration = 1:second.half.max
      exponential = ((second.half.max+1)):first.half.min
      deceleration = (first.half.min+1):second.half.min
      stationary = (second.half.min+1):nrow(df)
    }else if( df$acl[1] == half.max & df$acl[first.half.min] < which(df$acl == min(df$acl)) & 
              df$acl[second.half.min] > which(df$acl == min(df$acl))| 
              which(df$acl == max(df$acl)) == 1 & df$acl[1] < max(df$acl) & 
              df$acl[first.half.min] < which(df$acl == min(df$acl)) &
              df$acl[second.half.min] > which(df$acl == min(df$acl)) |
              (first.half.max == 0 & is.na(second.half.max) == T & is.na(first.half.min) == F & 
               is.na(second.half.min) == F) ){
      #if there's no lag nor acceleration phase
      lag = NULL
      acceleration = NULL
      exponential = 1:first.half.min
      deceleration = (first.half.min+1):second.half.min
      stationary = (second.half.min+1):nrow(df)
    }else if( df$acl[1] > half.max & df$acl[1] <= round(half.min) | df$acl[1] < half.min ){
      #if there's no lag, acceleration phase
      lag = NULL
      acceleration = NULL
      exponential = NULL
      deceleration = 1:first.half.min
      stationary = (second.half.min+1):nrow(df)
    }else if( is.na(second.half.min) == F & (second.half.min+1) >= nrow(df) | 
              (is.na(second.half.min) == T & df$acl[length(df$acl)] > min(df$acl) & 
               first.half.max < which(df$acl == max(df$acl)) & 
               second.half.max > which(df$acl == max(df$acl)) & is.na(second.half.max) == F & 
               length(df$acl) > second.half.max & is.na(first.half.min) == F)|
              (is.na(second.half.min) == T & is.na(first.half.min) == F & 
               first.half.min < which(df$acl == min(df$acl)) & 
               first.half.max < which(df$acl == max(df$acl)) & 
               second.half.max > which(df$acl == max(df$acl)) & is.na(second.half.max) == F & 
               length(df$acl) > second.half.max) ){
      #if there's no stationary phase
      lag = 1:first.half.max
      acceleration = (first.half.max+1):second.half.max
      exponential = (second.half.max):first.half.min
      deceleration = (first.half.min+1):nrow(df)
      stationary = NULL
    }else if( (which(df$acl == max(df$acl)) > which(df$acl == min(df$acl)) & is.na(first.half.min) == F) | 
              (is.na(first.half.min) == T & is.na(second.half.max) == F & length(df$acl) >  second.half.max) ){
      #if there's only lag, acceleration, exponential phases
      lag = 1:first.half.max
      acceleration = (first.half.max+1):second.half.max
      exponential = (second.half.max):nrow(df)
      deceleration = NULL
      stationary = NULL
    }else if( (max(df$acl) == df$acl[length(df$acl)] & is.na(second.half.max) != T) |
              (first.half.max < which(df$acl == max(df$acl)) & is.na(first.half.min) == T & 
               is.na(second.half.max) == T) |
              (is.na(first.half.min) == T & is.na(second.half.min) == T & df$acl[length(df$acl)] >= half.max) ){
      #if there's only lag and acceleration phases 
      lag = 1:first.half.max
      acceleration = (first.half.max+1):nrow(df)
      exponential = NULL
      deceleration = NULL
      stationary = NULL
    }else if( is.na(first.half.max) == T ){
      #if there's only lag, acceleration, exponential phases
      lag = 1:nrow(df)
      acceleration = NULL
      exponential = NULL
      deceleration = NULL
      stationary = NULL
    }else{
      lag = 1:first.half.max
      acceleration = (first.half.max+1):second.half.max
      exponential = (second.half.max):first.half.min
      deceleration = (first.half.min+1):second.half.min
      stationary = (second.half.min+1):nrow(df)
    }
    
    #if there's only one time point for a particular phase, calculating the coefficient of variation or correlation coefficient is not feasible.
    #in this circumstance, we borrow information of the growth rate and abundance from neighboring values to proceed.
    if (length(lag) == 1) {
      lag = c(lag, lag[length(lag)]+1)
    }
    if(length(acceleration) == 1){
      acceleration = c(acceleration[length(acceleration)]-1, acceleration, acceleration[length(acceleration)]+1)
    }
    if(length(exponential) == 1){
      exponential = c(exponential[length(exponential)]-1, exponential, exponential[length(exponential)]+1)
    }
    if(length(deceleration) == 1){
      deceleration = c(deceleration[length(deceleration)]-1, deceleration, deceleration[length(deceleration)]+1)
    }
    if(length(stationary) == 1){
      stationary = c(stationary[length(stationary)]-1, stationary)
    }
    
    #order growth phases
    if ( is.null(lag) == F ) {
      lag = lag[order(lag)]
    }
    if( is.null(acceleration) == F){
      acceleration = acceleration[order(acceleration)]
    }
    if( is.null(exponential) == F ){
      exponential = exponential[order(exponential)]
    }
    if( is.null(deceleration) == F ){
      deceleration = deceleration[order(deceleration)]
    }
    if( is.null(stationary) == F ){
      stationary = stationary[order(stationary)]
    }
    
    
  }
  
  return(list(lag, acceleration, exponential, deceleration, stationary))
}
