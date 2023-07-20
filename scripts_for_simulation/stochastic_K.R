stochastic_K = function(r, k, N = 100, d = 0.1, mean = 0, sd = 1){
  r = r
  k = k
  t = seq(0, N, d)
  w = rnorm(n = length(t), mean = mean, sd = sd)
  x = ((k+w)*exp(r*t))/((k+w) + exp(r*t)-1)
  x = x[1:(N+1)]
  
  df = data.frame(t = 0:N, sol = x)
  return(df)
}
