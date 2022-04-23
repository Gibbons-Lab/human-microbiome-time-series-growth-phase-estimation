
#fit acceleration equation
second_der = function(x,k){
  k^2*x*(1-(x/k))*(1-(2*x/k))
}