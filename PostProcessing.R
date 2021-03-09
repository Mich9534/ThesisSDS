# Post-processing graph/probability p(k|t) (3.7) formula
source('ApproximationVn.R')

pk.t <- function(k, t, n, gamma_){
  lgk.n <- log(gamma_ * k)
  for(i in 1:(n-1)){
    lgk.n <- lgk.n + log(gamma_* k + i)
  }
  lres <- -flv_n.app(n, t, gamma_) + log(fk_t(t, k))-lgk.n+log(dtpois(k, 1, 0))
  return(exp(lres))
}

pk.t(2,2, 50,1)

