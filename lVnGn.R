# Implementation of Vn using gnedin prior

l.Vn.Gn <- function(n,t, gamma.gn){
  return(lfactorial(t-1)+lgamma(2-gamma.gn)-lgamma(3-gamma.gn-t)+lgamma(gamma.gn+1)-lgamma(gamma.gn-n+t+1)-(lfactorial(n-1)+lgamma(gamma.gn+2)-lgamma(gamma.gn-n+3)))
}