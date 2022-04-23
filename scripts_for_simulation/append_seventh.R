append_seventh = function(df){
  if (ncol(df) != 7) {
    df[7] = 0
    df = as.data.frame(t(df))
    colnames(df) = c("Min", "1st.Qu", "Median", "Mean", "3rd.Qu", "Max", "NA's")
  }else if(ncol(df) == 7){
    df = unclass(df)
    #df = df[1]$V1
    df = as.data.frame(df)
    colnames(df) = c("Min", "1st.Qu", "Median", "Mean", "3rd.Qu", "Max", "NA's")
  }
  return(df)
}